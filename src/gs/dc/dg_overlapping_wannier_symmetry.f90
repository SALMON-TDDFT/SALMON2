module dg_overlapping_wannier_symmetry
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use iso_fortran_env,only:int64
  implicit none
  private
  integer,parameter::maximum_crystallographic_point_group_order=48
  integer,parameter::maximum_enumerated_subgroups=65536
  public::select_dg_exact_fragment_subgroup
  public::build_dg_fragment_group_representation
  public::project_dg_fragment_covariant_operators
  public::promote_dg_exact_global_subgroup
  public::build_dg_fragment_site_stabilizer
  public::evaluate_dg_covariance_residuals_by_operation
  public::fingerprint_dg_exact_fragment_symmetry
  public::build_dg_fragment_permuted_representation
  public::build_dg_fragment_symmetry_orbits
  public::build_dg_symmetry_constrained_pair_generator
  public::factor_dg_affine_translation_cocycle
  public::measure_dg_hamiltonian_density_commutators
  public::symmetrize_dg_distributed_pencil_rows
contains
  subroutine symmetrize_dg_distributed_pencil_rows(comm,row_ids,h_rows,s_rows,rho_rows,artifact_rows,&
      generator_representation,generator_operations,product_table,translation_operations,&
      coset_representatives,tolerance,sym_h_rows,sym_s_rows,&
      sym_rho_rows,before_residual,after_residual,artifact_change,artifact_magnitude,&
      workspace_peak_elements,ok,message,component_rows,component_residual)
    use mpi
    integer,intent(in)::comm,generator_operations(:),product_table(:,:),translation_operations(:),&
      coset_representatives(:)
    integer(int64),intent(in)::row_ids(:)
    complex(8),intent(in)::h_rows(:,:),s_rows(:,:),rho_rows(:,:),artifact_rows(:,:),&
      generator_representation(:,:,:)
    real(8),intent(in)::tolerance
    complex(8),allocatable,intent(out)::sym_h_rows(:,:),sym_s_rows(:,:),sym_rho_rows(:,:)
    real(8),intent(out)::before_residual(3),after_residual(3),artifact_change,artifact_magnitude
    integer(int64),intent(out)::workspace_peak_elements
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),intent(in),optional::component_rows(:,:,:)
    real(8),intent(out),optional::component_residual(:)
    integer::rank,nproc,ierr,n,noperation,ngenerator,identity,total_rows,r,i,j,operation,component,&
      generator,parent_operation,path_length,local_bad,global_bad
    integer,allocatable::row_counts(:),row_displs(:),parent(:),parent_generator(:),queue(:),path(:),&
      seen(:),group_seen(:)
    integer(int64),allocatable::all_row_ids(:)
    complex(8),allocatable::representation(:,:),identity_matrix(:,:),transformed(:,:),&
      sym_artifact_rows(:,:),check_rows(:,:),translation_h_rows(:,:),translation_s_rows(:,:),&
      translation_rho_rows(:,:),translation_artifact_rows(:,:)
    real(8)::local_change,local_magnitude,scale(3)
    ok=.false.;message='';before_residual=huge(1d0);after_residual=huge(1d0)
    artifact_change=huge(1d0);artifact_magnitude=huge(1d0);workspace_peak_elements=0_int64
    if(present(component_residual))component_residual=huge(1d0)
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    n=size(h_rows,2);noperation=size(product_table,1);ngenerator=size(generator_operations)
    local_bad=0
    if(ierr/=MPI_SUCCESS.or.n<1.or.noperation<1.or.ngenerator<1.or.&
        size(s_rows,2)/=n.or.size(rho_rows,2)/=n.or.size(artifact_rows,2)/=n.or.&
        size(h_rows,1)/=size(row_ids).or.size(s_rows,1)/=size(row_ids).or.&
        size(rho_rows,1)/=size(row_ids).or.size(artifact_rows,1)/=size(row_ids).or.&
        any(shape(product_table)/=[noperation,noperation]).or.&
        any(shape(generator_representation)/=[n,n,ngenerator]).or.&
        any(generator_operations<1).or.any(generator_operations>noperation).or.&
        size(translation_operations)<1.or.size(coset_representatives)<1.or.&
        any(translation_operations<1).or.any(translation_operations>noperation).or.&
        any(coset_representatives<1).or.any(coset_representatives>noperation).or.&
        any(product_table<1).or.any(product_table>noperation).or.tolerance<=0d0.or.&
        any(row_ids<1_int64).or.any(row_ids>int(n,int64)))local_bad=1
    if(present(component_rows).neqv.present(component_residual))local_bad=1
    if(present(component_rows))then
      if(size(component_rows,1)/=size(row_ids).or.size(component_rows,2)/=n.or.&
          size(component_rows,3)/=size(component_residual))local_bad=1
      if(.not.all(ieee_is_finite(real(component_rows))).or.&
          .not.all(ieee_is_finite(aimag(component_rows))))local_bad=1
    endif
    if(.not.all(ieee_is_finite(real(h_rows))).or..not.all(ieee_is_finite(aimag(h_rows))).or.&
        .not.all(ieee_is_finite(real(s_rows))).or..not.all(ieee_is_finite(aimag(s_rows))).or.&
        .not.all(ieee_is_finite(real(rho_rows))).or..not.all(ieee_is_finite(aimag(rho_rows))).or.&
        .not.all(ieee_is_finite(real(artifact_rows))).or.&
        .not.all(ieee_is_finite(aimag(artifact_rows))).or.&
        .not.all(ieee_is_finite(real(generator_representation))).or.&
        .not.all(ieee_is_finite(aimag(generator_representation))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='invalid distributed pencil symmetry contract';return;endif
    identity=0
    do operation=1,noperation
      if(all(product_table(operation,:)==[(i,i=1,noperation)]).and.&
          all(product_table(:,operation)==[(i,i=1,noperation)]))then
        identity=operation;exit
      endif
    enddo
    if(identity==0)then;message='pencil symmetry product table has no identity';return;endif
    allocate(group_seen(noperation));group_seen=0
    do i=1,size(translation_operations);do j=1,size(coset_representatives)
      operation=product_table(translation_operations(i),coset_representatives(j))
      group_seen(operation)=group_seen(operation)+1
    enddo;enddo
    if(any(group_seen/=1))then
      message='translation-coset factorization does not cover the affine group exactly once';return
    endif
    allocate(row_counts(nproc),row_displs(nproc))
    call MPI_Allgather(size(row_ids),1,MPI_INTEGER,row_counts,1,MPI_INTEGER,comm,ierr)
    total_rows=0
    do r=1,nproc;row_displs(r)=total_rows;total_rows=total_rows+row_counts(r);enddo
    if(total_rows/=n)local_bad=1
    allocate(all_row_ids(total_rows),seen(n));seen=0
    call MPI_Allgatherv(row_ids,size(row_ids),MPI_INTEGER8,all_row_ids,row_counts,row_displs,&
      MPI_INTEGER8,comm,ierr)
    do i=1,total_rows
      if(all_row_ids(i)<1_int64.or.all_row_ids(i)>int(n,int64))then
        local_bad=1
      else
        seen(int(all_row_ids(i)))=seen(int(all_row_ids(i)))+1
      endif
    enddo
    if(any(seen/=1))local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,local_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(local_bad/=0.or.ierr/=MPI_SUCCESS)then;message='duplicate or missing pencil row owner';return;endif
    allocate(parent(noperation),parent_generator(noperation),queue(noperation),path(noperation))
    parent=0;parent_generator=0;parent(identity)=identity;queue(1)=identity;i=1;j=1
    do while(i<=j)
      parent_operation=queue(i);i=i+1
      do generator=1,ngenerator
        operation=product_table(parent_operation,generator_operations(generator))
        if(parent(operation)/=0)cycle
        parent(operation)=parent_operation;parent_generator(operation)=generator
        j=j+1;queue(j)=operation
      enddo
    enddo
    if(any(parent==0))then;message='pencil symmetry generators do not close the full affine group';return;endif
    allocate(identity_matrix(n,n),representation(n,n),sym_h_rows(size(row_ids),n),&
      sym_s_rows(size(row_ids),n),sym_rho_rows(size(row_ids),n),sym_artifact_rows(size(row_ids),n))
    identity_matrix=0d0;do i=1,n;identity_matrix(i,i)=1d0;enddo
    do generator=1,ngenerator
      transformed=matmul(conjg(transpose(generator_representation(:,:,generator))),&
        generator_representation(:,:,generator))-identity_matrix
      if(maxval(abs(transformed))>tolerance)then
        message='pencil symmetry generator representation is not unitary';return
      endif
    enddo
    allocate(translation_h_rows(size(row_ids),n),translation_s_rows(size(row_ids),n),&
      translation_rho_rows(size(row_ids),n),translation_artifact_rows(size(row_ids),n))
    translation_h_rows=0d0;translation_s_rows=0d0;translation_rho_rows=0d0
    translation_artifact_rows=0d0
    do j=1,size(translation_operations)
      call build_representation(translation_operations(j),representation)
      call transform_rows(h_rows,representation,transformed);translation_h_rows=translation_h_rows+transformed
      call transform_rows(s_rows,representation,transformed);translation_s_rows=translation_s_rows+transformed
      call transform_rows(rho_rows,representation,transformed)
      translation_rho_rows=translation_rho_rows+transformed
      call transform_rows(artifact_rows,representation,transformed)
      translation_artifact_rows=translation_artifact_rows+transformed
    enddo
    translation_h_rows=translation_h_rows/real(size(translation_operations),8)
    translation_s_rows=translation_s_rows/real(size(translation_operations),8)
    translation_rho_rows=translation_rho_rows/real(size(translation_operations),8)
    translation_artifact_rows=translation_artifact_rows/real(size(translation_operations),8)
    sym_h_rows=0d0;sym_s_rows=0d0;sym_rho_rows=0d0;sym_artifact_rows=0d0
    do j=1,size(coset_representatives)
      call build_representation(coset_representatives(j),representation)
      call transform_rows(translation_h_rows,representation,transformed);sym_h_rows=sym_h_rows+transformed
      call transform_rows(translation_s_rows,representation,transformed);sym_s_rows=sym_s_rows+transformed
      call transform_rows(translation_rho_rows,representation,transformed)
      sym_rho_rows=sym_rho_rows+transformed
      call transform_rows(translation_artifact_rows,representation,transformed)
      sym_artifact_rows=sym_artifact_rows+transformed
    enddo
    sym_h_rows=sym_h_rows/real(size(coset_representatives),8)
    sym_s_rows=sym_s_rows/real(size(coset_representatives),8)
    sym_rho_rows=sym_rho_rows/real(size(coset_representatives),8)
    sym_artifact_rows=sym_artifact_rows/real(size(coset_representatives),8)
    scale=[global_scale(h_rows),global_scale(s_rows),global_scale(rho_rows)]
    before_residual=0d0;after_residual=0d0
    if(present(component_residual))component_residual=0d0
    do generator=1,ngenerator
      call transform_rows(h_rows,generator_representation(:,:,generator),check_rows)
      before_residual(1)=max(before_residual(1),global_difference(check_rows,h_rows)/scale(1))
      call transform_rows(s_rows,generator_representation(:,:,generator),check_rows)
      before_residual(2)=max(before_residual(2),global_difference(check_rows,s_rows)/scale(2))
      call transform_rows(rho_rows,generator_representation(:,:,generator),check_rows)
      before_residual(3)=max(before_residual(3),global_difference(check_rows,rho_rows)/scale(3))
      if(present(component_rows))then
        do component=1,size(component_residual)
          call transform_rows(component_rows(:,:,component),&
            generator_representation(:,:,generator),check_rows)
          component_residual(component)=max(component_residual(component),&
            global_difference(check_rows,component_rows(:,:,component))/&
            global_scale(component_rows(:,:,component)))
        enddo
      endif
      call transform_rows(sym_h_rows,generator_representation(:,:,generator),check_rows)
      after_residual(1)=max(after_residual(1),global_difference(check_rows,sym_h_rows)/scale(1))
      call transform_rows(sym_s_rows,generator_representation(:,:,generator),check_rows)
      after_residual(2)=max(after_residual(2),global_difference(check_rows,sym_s_rows)/scale(2))
      call transform_rows(sym_rho_rows,generator_representation(:,:,generator),check_rows)
      after_residual(3)=max(after_residual(3),global_difference(check_rows,sym_rho_rows)/scale(3))
    enddo
    local_change=0d0;local_magnitude=0d0
    if(size(artifact_rows)>0)then
      local_change=maxval(abs(sym_artifact_rows-artifact_rows))
      local_magnitude=maxval(abs(artifact_rows))
    endif
    call MPI_Allreduce(local_change,artifact_change,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_magnitude,artifact_magnitude,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,workspace_peak_elements,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.maxval(after_residual)>tolerance)then
      message='full-group pencil average is not generator invariant';return
    endif
    ok=.true.
  contains
    subroutine build_representation(target,d)
      integer,intent(in)::target
      complex(8),intent(out)::d(:,:)
      integer::step
      d=identity_matrix;path_length=0;parent_operation=target
      do while(parent_operation/=identity)
        path_length=path_length+1;path(path_length)=parent_generator(parent_operation)
        parent_operation=parent(parent_operation)
      enddo
      do step=path_length,1,-1
        d=matmul(d,generator_representation(:,:,path(step)))
      enddo
    end subroutine build_representation
    subroutine transform_rows(input_rows,d,output_rows)
      complex(8),intent(in)::input_rows(:,:),d(:,:)
      complex(8),allocatable,intent(out)::output_rows(:,:)
      complex(8),allocatable::right_rows(:,:),partial(:,:),reduced(:,:)
      integer,parameter::batch_size=32
      integer::owner_rank,owner_rows,first,count_rows,a,b
      allocate(output_rows(size(row_ids),n),right_rows(size(row_ids),n));output_rows=0d0
      ! Spatial basis rows obey G phi = d phi.  Operator matrix elements
      ! therefore transform as conjg(d) A transpose(d), not d^H A d.
      right_rows=matmul(input_rows,transpose(d))
      do owner_rank=0,nproc-1
        owner_rows=row_counts(owner_rank+1)
        do first=1,owner_rows,batch_size
          count_rows=min(batch_size,owner_rows-first+1)
          allocate(partial(count_rows,n),reduced(count_rows,n));partial=0d0
          do b=1,n;do a=1,count_rows
            partial(a,b)=sum(conjg(d(int(all_row_ids(&
              row_displs(owner_rank+1)+first+a-1)),int(row_ids)))*right_rows(:,b))
          enddo;enddo
          call MPI_Reduce(partial,reduced,count_rows*n,MPI_DOUBLE_COMPLEX,MPI_SUM,owner_rank,comm,ierr)
          if(rank==owner_rank)output_rows(first:first+count_rows-1,:)=reduced
          workspace_peak_elements=max(workspace_peak_elements,int(n*n+size(right_rows)+&
            size(output_rows)+size(partial)+size(reduced)+size(sym_h_rows)+size(sym_s_rows)+&
            size(sym_rho_rows)+size(sym_artifact_rows)+size(translation_h_rows)+&
            size(translation_s_rows)+size(translation_rho_rows)+size(translation_artifact_rows),int64))
          deallocate(partial,reduced)
        enddo
      enddo
    end subroutine transform_rows
    real(8) function global_scale(rows)
      complex(8),intent(in)::rows(:,:)
      real(8)::local
      local=0d0;if(size(rows)>0)local=maxval(abs(rows))
      call MPI_Allreduce(local,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      global_scale=max(1d0,global_scale)
    end function global_scale
    real(8) function global_difference(left,right)
      complex(8),intent(in)::left(:,:),right(:,:)
      real(8)::local
      local=0d0;if(size(left)>0)local=maxval(abs(left-right))
      call MPI_Allreduce(local,global_difference,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    end function global_difference
  end subroutine symmetrize_dg_distributed_pencil_rows

  subroutine measure_dg_hamiltonian_density_commutators(representation,hamiltonian,density_projector,&
      tolerance,hamiltonian_residual,density_residual,ok,message)
    complex(8),intent(in)::representation(:,:,:),hamiltonian(:,:),density_projector(:,:)
    real(8),intent(in)::tolerance
    real(8),allocatable,intent(out)::hamiltonian_residual(:),density_residual(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::work(:,:),identity(:,:)
    real(8)::hamiltonian_scale,density_scale,unitarity_defect
    integer::n,noperation,operation,i

    ok=.false.;message='';n=size(representation,1);noperation=size(representation,3)
    if(n<1.or.size(representation,2)/=n.or.noperation<1.or.&
        any(shape(hamiltonian)/=[n,n]).or.any(shape(density_projector)/=[n,n]).or.&
        tolerance<=0d0.or..not.ieee_is_finite(tolerance).or.&
        .not.all(ieee_is_finite(real(representation))).or.&
        .not.all(ieee_is_finite(aimag(representation))).or.&
        .not.all(ieee_is_finite(real(hamiltonian))).or.&
        .not.all(ieee_is_finite(aimag(hamiltonian))).or.&
        .not.all(ieee_is_finite(real(density_projector))).or.&
        .not.all(ieee_is_finite(aimag(density_projector))))then
      message='invalid Hamiltonian-density commutator contract';return
    end if
    allocate(work(n,n),identity(n,n),hamiltonian_residual(noperation),density_residual(noperation))
    identity=(0d0,0d0);do i=1,n;identity(i,i)=1d0;end do
    unitarity_defect=0d0
    do operation=1,noperation
      work=matmul(conjg(transpose(representation(:,:,operation))),representation(:,:,operation))-identity
      unitarity_defect=max(unitarity_defect,maxval(abs(work)))
    end do
    if(unitarity_defect>tolerance)then
      message='Hamiltonian-density symmetry representation is not unitary';return
    end if
    hamiltonian_scale=max(1d0,maxval(abs(hamiltonian)))
    density_scale=max(1d0,maxval(abs(density_projector)))
    do operation=1,noperation
      work=matmul(hamiltonian,representation(:,:,operation))-&
        matmul(representation(:,:,operation),hamiltonian)
      hamiltonian_residual(operation)=maxval(abs(work))/hamiltonian_scale
      work=matmul(density_projector,representation(:,:,operation))-&
        matmul(representation(:,:,operation),density_projector)
      density_residual(operation)=maxval(abs(work))/density_scale
    end do
    ok=.true.
  end subroutine measure_dg_hamiltonian_density_commutators

  subroutine factor_dg_affine_translation_cocycle(rotations,translations,product_table,tolerance,&
      translation_subgroup,point_representatives,point_product,translation_cocycle,ok,message)
    integer,intent(in)::rotations(:,:,:),product_table(:,:)
    real(8),intent(in)::translations(:,:),tolerance
    integer,allocatable,intent(out)::translation_subgroup(:),point_representatives(:),&
      point_product(:,:),translation_cocycle(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::operation_class(:),inverse(:),work_representatives(:),work_translations(:)
    integer::nop,operation,left,right,product,identity,npoint,ntranslation,point_left,point_right,&
      point_target,cocycle_operation,cocycle_position,i,j,k
    integer::identity_rotation(3,3)
    real(8)::composed_translation(3),difference(3)

    ok=.false.;message='';nop=size(rotations,3)
    identity_rotation=reshape([1,0,0,0,1,0,0,0,1],[3,3])
    if(nop<1.or.size(rotations,1)/=3.or.size(rotations,2)/=3.or.&
        any(shape(translations)/=[3,nop]).or.any(shape(product_table)/=[nop,nop]).or.&
        tolerance<=0d0.or..not.ieee_is_finite(tolerance).or.&
        .not.all(ieee_is_finite(translations)).or.any(product_table<1).or.&
        any(product_table>nop))then
      message='invalid affine translation-cocycle contract';return
    end if
    identity=0
    do operation=1,nop
      if(all(product_table(operation,:)==[(i,i=1,nop)]).and.&
          all(product_table(:,operation)==[(i,i=1,nop)]))then
        identity=operation;exit
      end if
    end do
    if(identity==0.or.any(rotations(:,:,identity)/=identity_rotation).or.&
        maxval(abs(translations(:,identity)-anint(translations(:,identity))))>tolerance)then
      message='affine group has no valid identity';return
    end if
    allocate(inverse(nop));inverse=0
    do left=1,nop
      do right=1,nop
        if(product_table(left,right)==identity.and.product_table(right,left)==identity)then
          inverse(left)=right;exit
        end if
      end do
      if(inverse(left)==0)then;message='affine product table has no inverse';return;end if
    end do
    do i=1,nop;do j=1,nop;do k=1,nop
      if(product_table(product_table(i,j),k)/=product_table(i,product_table(j,k)))then
        message='affine product table is not associative';return
      end if
    end do;end do;end do
    do left=1,nop;do right=1,nop
      product=product_table(left,right)
      if(any(matmul(rotations(:,:,left),rotations(:,:,right))/=rotations(:,:,product)))then
        message='affine rotations disagree with the product table';return
      end if
      composed_translation=translations(:,left)+matmul(real(rotations(:,:,left),8),translations(:,right))
      difference=composed_translation-translations(:,product)
      if(maxval(abs(difference-anint(difference)))>tolerance)then
        message='affine translations disagree with the product table';return
      end if
    end do;end do
    allocate(work_translations(nop),work_representatives(nop),operation_class(nop))
    ntranslation=0
    do operation=1,nop
      if(all(rotations(:,:,operation)==identity_rotation))then
        ntranslation=ntranslation+1;work_translations(ntranslation)=operation
      end if
    end do
    if(ntranslation<1)then;message='affine group has no pure-translation subgroup';return;end if
    npoint=0;operation_class=0
    do operation=1,nop
      do point_target=1,npoint
        if(all(rotations(:,:,operation)==rotations(:,:,work_representatives(point_target))))then
          operation_class(operation)=point_target;exit
        end if
      end do
      if(operation_class(operation)==0)then
        npoint=npoint+1;work_representatives(npoint)=operation;operation_class(operation)=npoint
      end if
    end do
    do point_target=1,npoint
      if(count(operation_class==point_target)/=ntranslation)then
        message='affine rotation classes are not complete translation cosets';return
      end if
    end do
    allocate(translation_subgroup(ntranslation),point_representatives(npoint),&
      point_product(npoint,npoint),translation_cocycle(npoint,npoint))
    translation_subgroup=work_translations(1:ntranslation)
    point_representatives=work_representatives(1:npoint)
    do point_left=1,npoint;do point_right=1,npoint
      product=product_table(point_representatives(point_left),point_representatives(point_right))
      point_target=operation_class(product);point_product(point_left,point_right)=point_target
      cocycle_operation=product_table(product,inverse(point_representatives(point_target)))
      cocycle_position=findloc(translation_subgroup,cocycle_operation,dim=1)
      if(cocycle_position==0)then
        message='affine representative product does not yield a pure-translation cocycle';return
      end if
      translation_cocycle(point_left,point_right)=cocycle_position
    end do;end do
    ok=.true.
  end subroutine factor_dg_affine_translation_cocycle

  subroutine build_dg_symmetry_constrained_pair_generator(first,second,amplitude,representation,&
      product_table,tolerance,generator,antihermiticity_defect,commutator_defect,&
      active_indices,ok,message)
    integer,intent(in)::first,second
    complex(8),intent(in)::amplitude,representation(:,:,:)
    integer,intent(in)::product_table(:,:)
    real(8),intent(in)::tolerance
    complex(8),allocatable,intent(out)::generator(:,:)
    real(8),intent(out)::antihermiticity_defect,commutator_defect
    integer,allocatable,intent(out)::active_indices(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::identity(:,:),difference(:,:),work(:,:),block_representation(:,:,:),seed(:,:)
    logical,allocatable::active(:)
    real(8)::unitarity_defect,closure_defect,leakage
    integer::n,noperation,operation,left,right,product,i,source,target,nactive,first_local,second_local
    logical::changed

    ok=.false.;message='';antihermiticity_defect=huge(1d0);commutator_defect=huge(1d0)
    n=size(representation,1);noperation=size(representation,3)
    if(n<2.or.size(representation,2)/=n.or.noperation<1.or.first<1.or.second<1.or.&
        first>n.or.second>n.or.first==second.or.&
        .not.ieee_is_finite(real(amplitude)).or..not.ieee_is_finite(aimag(amplitude)).or.&
        abs(amplitude)<=tiny(1d0).or.any(shape(product_table)/=[noperation,noperation]).or.&
        tolerance<=0d0.or..not.ieee_is_finite(tolerance).or.&
        .not.all(ieee_is_finite(real(representation))).or.&
        .not.all(ieee_is_finite(aimag(representation))))then
      message='invalid symmetry-constrained pair-generator contract';return
    end if
    if(any(product_table<1).or.any(product_table>noperation))then
      message='symmetry-constrained generator product table is invalid';return
    end if
    allocate(active(n));active=.false.;active([first,second])=.true.
    do
      changed=.false.
      do operation=1,noperation;do source=1,n
        if(.not.active(source))cycle
        do target=1,n
          if(abs(representation(target,source,operation))<=tolerance.or.active(target))cycle
          active(target)=.true.;changed=.true.
        end do
      end do;end do
      if(.not.changed)exit
    end do
    nactive=count(active);allocate(active_indices(nactive));active_indices=pack([(i,i=1,n)],active)
    allocate(identity(nactive,nactive),difference(nactive,nactive),work(nactive,nactive),&
      generator(nactive,nactive),block_representation(nactive,nactive,noperation),seed(nactive,nactive))
    identity=(0d0,0d0);do i=1,nactive;identity(i,i)=1d0;end do
    do operation=1,noperation
      block_representation(:,:,operation)=representation(active_indices,active_indices,operation)
    end do
    leakage=0d0
    do operation=1,noperation;do source=1,nactive;do target=1,n
      if(active(target))cycle
      leakage=max(leakage,abs(representation(target,active_indices(source),operation)))
    end do;end do;end do
    if(leakage>tolerance)then
      message='symmetry-constrained pair support block is not invariant';return
    end if
    unitarity_defect=0d0
    do operation=1,noperation
      difference=matmul(conjg(transpose(block_representation(:,:,operation))),&
        block_representation(:,:,operation))-identity
      unitarity_defect=max(unitarity_defect,maxval(abs(difference)))
    end do
    if(unitarity_defect>tolerance)then
      message='symmetry-constrained representation is not unitary';return
    end if
    closure_defect=0d0
    do left=1,noperation;do right=1,noperation
      product=product_table(left,right)
      difference=matmul(block_representation(:,:,left),block_representation(:,:,right))-&
        block_representation(:,:,product)
      closure_defect=max(closure_defect,maxval(abs(difference)))
    end do;end do
    if(closure_defect>tolerance)then
      message='symmetry-constrained representation is not group closed';return
    end if
    seed=(0d0,0d0);first_local=findloc(active_indices,first,dim=1)
    second_local=findloc(active_indices,second,dim=1)
    seed(first_local,second_local)=amplitude;seed(second_local,first_local)=-conjg(amplitude)
    generator=(0d0,0d0)
    do operation=1,noperation
      generator=generator+matmul(block_representation(:,:,operation),&
        matmul(seed,conjg(transpose(block_representation(:,:,operation)))))
    end do
    generator=generator/real(noperation,8)
    antihermiticity_defect=maxval(abs(generator+conjg(transpose(generator))))/&
      max(1d0,maxval(abs(generator)))
    commutator_defect=0d0
    do operation=1,noperation
      work=matmul(generator,block_representation(:,:,operation))-&
        matmul(block_representation(:,:,operation),generator)
      commutator_defect=max(commutator_defect,maxval(abs(work))/&
        max(1d0,maxval(abs(generator))))
    end do
    if(antihermiticity_defect>tolerance)then
      message='symmetry-constrained generator lost anti-Hermiticity';return
    end if
    if(commutator_defect>tolerance)then
      message='symmetry-constrained generator does not commute with the group';return
    end if
    ok=.true.
  end subroutine build_dg_symmetry_constrained_pair_generator

  subroutine build_dg_fragment_symmetry_orbits(fragment_permutation,fragment_orbit,&
      orbit_representative,ok,message)
    integer,intent(in)::fragment_permutation(:,:)
    integer,allocatable,intent(out)::fragment_orbit(:),orbit_representative(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::representative(:)
    integer::nfragment,noperation,operation,source,target,root
    logical::changed

    ok=.false.;message='';nfragment=size(fragment_permutation,1)
    noperation=size(fragment_permutation,2)
    if(nfragment<1.or.noperation<1)then
      message='invalid fragment symmetry-orbit contract';return
    end if
    do operation=1,noperation
      if(any(fragment_permutation(:,operation)<1).or.&
          any(fragment_permutation(:,operation)>nfragment))then
        message='fragment symmetry operation is not a permutation';return
      end if
      do target=1,nfragment
        if(count(fragment_permutation(:,operation)==target)/=1)then
          message='fragment symmetry operation is not a permutation';return
        end if
      end do
    end do
    allocate(representative(nfragment));representative=[(source,source=1,nfragment)]
    do
      changed=.false.
      do operation=1,noperation;do source=1,nfragment
        target=fragment_permutation(source,operation)
        root=min(representative(source),representative(target))
        if(representative(source)/=root.or.representative(target)/=root)changed=.true.
        representative(source)=root;representative(target)=root
      end do;end do
      do source=1,nfragment
        root=representative(source)
        do while(representative(root)/=root);root=representative(root);end do
        if(representative(source)/=root)changed=.true.
        representative(source)=root
      end do
      if(.not.changed)exit
    end do
    allocate(fragment_orbit(nfragment),orbit_representative(nfragment))
    do source=1,nfragment
      root=representative(source);orbit_representative(source)=root
      fragment_orbit(source)=count([(representative(target)==target,target=1,root)])
    end do
    ok=.true.
  end subroutine build_dg_fragment_symmetry_orbits

  subroutine build_dg_fragment_permuted_representation(local_representation,rotations,&
      fragment_centers,tolerance,global_representation,fragment_permutation,ok,message)
    complex(8),intent(in)::local_representation(:,:,:)
    real(8),intent(in)::rotations(:,:,:),fragment_centers(:,:),tolerance
    complex(8),allocatable,intent(out)::global_representation(:,:,:)
    integer,allocatable,intent(out)::fragment_permutation(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8)::mapped(3),difference(3),residual,best
    integer::nlocal,nfragment,nop,operation,source,target,best_target
    ok=.false.;message='';nlocal=size(local_representation,1)
    nfragment=size(fragment_centers,2);nop=size(local_representation,3)
    if(nlocal<1.or.size(local_representation,2)/=nlocal.or.nfragment<1.or.nop<1.or. &
        size(fragment_centers,1)/=3.or.size(rotations,1)/=3.or.size(rotations,2)/=3.or. &
        size(rotations,3)/=nop.or.tolerance<=0d0.or..not.ieee_is_finite(tolerance).or. &
        .not.all(ieee_is_finite(fragment_centers)).or..not.all(ieee_is_finite(rotations)).or. &
        .not.all(ieee_is_finite(real(local_representation))).or. &
        .not.all(ieee_is_finite(aimag(local_representation))))then
      message='invalid fragment-permuted representation contract';return
    end if
    allocate(global_representation(nlocal*nfragment,nlocal*nfragment,nop),&
      fragment_permutation(nfragment,nop));global_representation=(0d0,0d0)
    do operation=1,nop;do source=1,nfragment
      mapped=modulo(matmul(rotations(:,:,operation),fragment_centers(:,source)-0.5d0)+0.5d0,1d0)
      best=huge(1d0);best_target=0
      do target=1,nfragment
        difference=mapped-fragment_centers(:,target);difference=difference-anint(difference)
        residual=maxval(abs(difference))
        if(residual<best)then;best=residual;best_target=target;end if
      end do
      if(best_target==0.or.best>tolerance)then
        message='point operation does not permute fragment centers';return
      end if
      fragment_permutation(source,operation)=best_target
      global_representation((best_target-1)*nlocal+1:best_target*nlocal,&
        (source-1)*nlocal+1:source*nlocal,operation)=local_representation(:,:,operation)
    end do;end do
    do operation=1,nop
      if(any([(count(fragment_permutation(:,operation)==target),target=1,nfragment)]/=1))then
        message='point operation fragment map is not a permutation';return
      end if
    end do
    ok=.true.
  end subroutine build_dg_fragment_permuted_representation

  integer(int64) function fingerprint_dg_exact_fragment_symmetry(rotations,product_table,tolerance,&
      fractional_translations)
    integer,intent(in)::rotations(:,:,:),product_table(:,:)
    real(8),intent(in)::tolerance
    real(8),intent(in),optional::fractional_translations(:,:)
    integer(int64)::word
    integer::i,j,k,shift,nop
    nop=size(rotations,3)
    if(nop<1.or.size(rotations,1)/=3.or.size(rotations,2)/=3.or. &
        size(product_table,1)/=nop.or.size(product_table,2)/=nop.or. &
        .not.ieee_is_finite(tolerance).or.tolerance<=0d0)then
      fingerprint_dg_exact_fragment_symmetry=0_int64;return
    end if
    if(present(fractional_translations))then
      if(any(shape(fractional_translations)/=[3,nop]))then
        fingerprint_dg_exact_fragment_symmetry=0_int64;return
      end if
    end if
    fingerprint_dg_exact_fragment_symmetry=ieor(int(nop,int64),int(z'243F6A8885A308D3',int64))
    do k=1,nop;do j=1,3;do i=1,3
      shift=modulo(11*i+17*j+23*k,63)
      word=int(rotations(i,j,k),int64)
      fingerprint_dg_exact_fragment_symmetry=ieor(fingerprint_dg_exact_fragment_symmetry, &
        ishftc(ieor(word,int(97*i+193*j+389*k,int64)),shift))
    end do;end do;end do
    do j=1,nop;do i=1,nop
      if(product_table(i,j)<1.or.product_table(i,j)>nop)then
        fingerprint_dg_exact_fragment_symmetry=0_int64;return
      end if
      shift=modulo(7*i+29*j,63)
      word=int(product_table(i,j),int64)
      fingerprint_dg_exact_fragment_symmetry=ieor(fingerprint_dg_exact_fragment_symmetry, &
        ishftc(ieor(word,int(521*i+1031*j,int64)),shift))
    end do;end do
    if(present(fractional_translations))then
      if(.not.all(ieee_is_finite(fractional_translations)))then
        fingerprint_dg_exact_fragment_symmetry=0_int64;return
      end if
      do k=1,nop;do i=1,3
        shift=modulo(19*i+31*k,63)
        word=transfer(modulo(fractional_translations(i,k),1d0),word)
        fingerprint_dg_exact_fragment_symmetry=ieor(fingerprint_dg_exact_fragment_symmetry,&
          ishftc(ieor(word,int(2053*i+4099*k,int64)),shift))
      end do;end do
    end if
    word=transfer(tolerance,word)
    fingerprint_dg_exact_fragment_symmetry=ieor(fingerprint_dg_exact_fragment_symmetry,ishftc(word,37))
    if(fingerprint_dg_exact_fragment_symmetry==0_int64)&
      fingerprint_dg_exact_fragment_symmetry=int(z'13198A2E03707344',int64)
  end function fingerprint_dg_exact_fragment_symmetry

  subroutine evaluate_dg_covariance_residuals_by_operation(representation,rotations,scalars,vectors, &
      scalar_residual,vector_residual,ok,message)
    complex(8),intent(in)::representation(:,:,:),scalars(:,:,:),vectors(:,:,:,:)
    real(8),intent(in)::rotations(:,:,:)
    real(8),allocatable,intent(out)::scalar_residual(:),vector_residual(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::transformed(:,:),target(:,:)
    real(8)::scalar_scale,vector_scale
    integer::n,nop,iop,i,a,b
    ok=.false.;message='';n=size(representation,1);nop=size(representation,3)
    if(n<1.or.size(representation,2)/=n.or.nop<1.or.size(rotations,1)/=3.or. &
        size(rotations,2)/=3.or.size(rotations,3)/=nop.or.size(scalars,1)/=n.or. &
        size(scalars,2)/=n.or.size(scalars,3)<1.or.size(vectors,1)/=n.or. &
        size(vectors,2)/=n.or.size(vectors,3)/=3.or.size(vectors,4)<1)then
      message='operation covariance residual arrays have inconsistent dimensions';return
    end if
    if(.not.all(ieee_is_finite(real(representation))).or. &
        .not.all(ieee_is_finite(aimag(representation))).or..not.all(ieee_is_finite(rotations)).or. &
        .not.all(ieee_is_finite(real(scalars))).or..not.all(ieee_is_finite(aimag(scalars))).or. &
        .not.all(ieee_is_finite(real(vectors))).or..not.all(ieee_is_finite(aimag(vectors))))then
      message='operation covariance residual payload is not finite';return
    end if
    allocate(scalar_residual(nop),vector_residual(nop),transformed(n,n),target(n,n))
    scalar_residual=0d0;vector_residual=0d0
    do iop=1,nop
      do i=1,size(scalars,3)
        scalar_scale=max(1d0,maxval(abs(scalars(:,:,i))))
        transformed=matmul(conjg(transpose(representation(:,:,iop))), &
          matmul(scalars(:,:,i),representation(:,:,iop)))
        scalar_residual(iop)=max(scalar_residual(iop),maxval(abs(transformed-scalars(:,:,i)))/scalar_scale)
      end do
      do i=1,size(vectors,4)
        vector_scale=max(1d0,maxval(abs(vectors(:,:,:,i))))
        do a=1,3
        transformed=matmul(conjg(transpose(representation(:,:,iop))), &
          matmul(vectors(:,:,a,i),representation(:,:,iop)))
        target=(0d0,0d0)
        do b=1,3;target=target+rotations(a,b,iop)*vectors(:,:,b,i);end do
        vector_residual(iop)=max(vector_residual(iop),maxval(abs(transformed-target))/vector_scale)
        end do
      end do
    end do
    ok=.true.
  end subroutine evaluate_dg_covariance_residuals_by_operation

  subroutine build_dg_fragment_site_stabilizer(rotations,translations,fragment_center,allowed, &
      tolerance,selected,product_table,maximum_site_residual,ok,message)
    integer,intent(in)::rotations(:,:,:)
    real(8),intent(in)::translations(:,:),fragment_center(3),tolerance
    logical,intent(in)::allowed(:)
    integer,allocatable,intent(out)::selected(:),product_table(:,:)
    real(8),intent(out)::maximum_site_residual
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::work_selected(:)
    integer(8)::product_rotation(3,3)
    real(8)::mapped_center(3),residual,product_translation(3)
    integer::nop,iop,jop,kop,nselected,identity_position,i,j,k
    logical::duplicate,matched
    ok=.false.;message='';maximum_site_residual=0d0;nop=size(rotations,3)
    if(nop<1.or.size(rotations,1)/=3.or.size(rotations,2)/=3.or. &
        size(translations,1)/=3.or.size(translations,2)/=nop.or.size(allowed)/=nop)then
      message='fragment site-stabilizer arrays have inconsistent dimensions';return
    end if
    if(.not.ieee_is_finite(tolerance).or.tolerance<=0d0.or. &
        .not.all(ieee_is_finite(translations)).or..not.all(ieee_is_finite(fragment_center)))then
      message='fragment site-stabilizer coordinates or tolerance are not finite';return
    end if
    allocate(work_selected(min(nop,maximum_crystallographic_point_group_order)))
    nselected=0
    do iop=1,nop
      if(.not.allowed(iop))cycle
      mapped_center=matmul(real(rotations(:,:,iop),8),fragment_center)+translations(:,iop)-fragment_center
      residual=maxval(abs(mapped_center-anint(mapped_center)))
      if(residual>tolerance)cycle
      duplicate=.false.
      do i=1,nselected
        if(all(rotations(:,:,work_selected(i))==rotations(:,:,iop)))then
          duplicate=.true.;exit
        end if
      end do
      if(duplicate)cycle
      if(nselected>=maximum_crystallographic_point_group_order)then
        message='fragment site stabilizer exceeds crystallographic point-group order 48';return
      end if
      nselected=nselected+1;work_selected(nselected)=iop
      maximum_site_residual=max(maximum_site_residual,residual)
    end do
    if(nselected<1)then;message='fragment site stabilizer has no admissible operation';return;end if
    identity_position=0
    do i=1,nselected
      iop=work_selected(i)
      if(all(rotations(:,:,iop)==reshape([1,0,0,0,1,0,0,0,1],[3,3])).and. &
          maxval(abs(translations(:,iop)-anint(translations(:,iop))))<=tolerance)then
        identity_position=i;exit
      end if
    end do
    if(identity_position==0)then;message='fragment site stabilizer is missing identity';return;end if
    if(identity_position/=1)then
      iop=work_selected(1);work_selected(1)=work_selected(identity_position)
      work_selected(identity_position)=iop
    end if
    allocate(selected(nselected),product_table(nselected,nselected));selected=work_selected(1:nselected)
    do i=1,nselected;do j=1,nselected
      iop=selected(i);jop=selected(j)
      product_rotation=matmul(int(rotations(:,:,iop),8),int(rotations(:,:,jop),8))
      product_translation=translations(:,iop)+matmul(real(rotations(:,:,iop),8),translations(:,jop))
      matched=.false.;kop=0
      do k=1,nselected
        if(any(product_rotation/=int(rotations(:,:,selected(k)),8)))cycle
        residual=maxval(abs(product_translation-translations(:,selected(k))- &
          anint(product_translation-translations(:,selected(k)))))
        if(residual<=tolerance)then;matched=.true.;kop=k;exit;end if
      end do
      if(.not.matched)then
        message='fragment site stabilizer affine operations are not closed';return
      end if
      product_table(i,j)=kop
    end do;end do
    ok=.true.
  end subroutine build_dg_fragment_site_stabilizer

  subroutine promote_dg_exact_global_subgroup(product_table,fragment_exact,scalar_block_residual, &
      vector_block_residual,tolerance,subgroup,ok,message)
    integer,intent(in)::product_table(:,:)
    logical,intent(in)::fragment_exact(:,:)
    real(8),intent(in)::scalar_block_residual(:,:),vector_block_residual(:,:),tolerance
    integer,allocatable,intent(out)::subgroup(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),allocatable::fragment_residual(:),scalar_residual(:),vector_residual(:),zero_residual(:)
    integer::nop,iop
    nop=size(product_table,1);ok=.false.;message=''
    if(nop<1.or.size(product_table,2)/=nop.or.size(fragment_exact,1)<1.or. &
        size(fragment_exact,2)/=nop.or.size(scalar_block_residual,1)<1.or. &
        size(scalar_block_residual,2)/=nop.or.size(vector_block_residual,1)<1.or. &
        size(vector_block_residual,2)/=nop)then
      message='global symmetry promotion arrays have inconsistent dimensions';return
    end if
    allocate(fragment_residual(nop),scalar_residual(nop),vector_residual(nop),zero_residual(nop))
    zero_residual=0d0
    do iop=1,nop
      fragment_residual(iop)=merge(0d0,2d0*tolerance,all(fragment_exact(:,iop)))
      scalar_residual(iop)=maxval(scalar_block_residual(:,iop))
      vector_residual(iop)=maxval(vector_block_residual(:,iop))
    end do
    call select_dg_exact_fragment_subgroup(product_table,fragment_residual,scalar_residual, &
      vector_residual,zero_residual,tolerance,tolerance,tolerance,tolerance,subgroup,ok,message)
    if(.not.ok)message='global symmetry promotion failed: '//trim(message)
  end subroutine promote_dg_exact_global_subgroup

  subroutine project_dg_fragment_covariant_operators(representation,rotations,scalars,vectors, &
      tolerance,projected_scalars,projected_vectors,pre_projection_defect,post_projection_defect,ok,message,&
      maximum_pre_projection_defect)
    complex(8),intent(in)::representation(:,:,:),scalars(:,:,:),vectors(:,:,:,:)
    real(8),intent(in)::rotations(:,:,:),tolerance
    complex(8),allocatable,intent(out)::projected_scalars(:,:,:),projected_vectors(:,:,:,:)
    real(8),intent(out)::pre_projection_defect,post_projection_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),intent(in),optional::maximum_pre_projection_defect
    complex(8),allocatable::transformed(:,:),target(:,:)
    real(8)::rotation_defect,determinant,correction_limit
    integer::n,nop,nscalar,nvector,iop,i,j,a,b
    ok=.false.;message='';pre_projection_defect=huge(1d0);post_projection_defect=huge(1d0)
    n=size(representation,1);nop=size(representation,3);nscalar=size(scalars,3);nvector=size(vectors,4)
    if(n<1.or.size(representation,2)/=n.or.nop<1.or.size(rotations,1)/=3.or. &
        size(rotations,2)/=3.or.size(rotations,3)/=nop.or.size(scalars,1)/=n.or. &
        size(scalars,2)/=n.or.nscalar<1.or.size(vectors,1)/=n.or.size(vectors,2)/=n.or. &
        size(vectors,3)/=3.or.nvector<1)then
      message='fragment covariant operator arrays have inconsistent dimensions';return
    end if
    if(.not.ieee_is_finite(tolerance).or.tolerance<=0d0.or. &
        .not.all(ieee_is_finite(real(representation))).or. &
        .not.all(ieee_is_finite(aimag(representation))).or. &
        .not.all(ieee_is_finite(rotations)).or..not.all(ieee_is_finite(real(scalars))).or. &
        .not.all(ieee_is_finite(aimag(scalars))).or..not.all(ieee_is_finite(real(vectors))).or. &
        .not.all(ieee_is_finite(aimag(vectors))))then
      message='fragment covariant operators or tolerance are not finite';return
    end if
    do iop=1,nop
      rotation_defect=maxval(abs(matmul(transpose(rotations(:,:,iop)),rotations(:,:,iop))-identity3()))
      determinant=determinant3(rotations(:,:,iop))
      if(rotation_defect>tolerance.or.abs(abs(determinant)-1d0)>tolerance)then
        message='fragment covariant operator Cartesian rotation is not orthogonal';return
      end if
    end do
    if(maxval(abs(representation(:,:,1)-complex_identity(n)))>tolerance.or. &
        maxval(abs(rotations(:,:,1)-identity3()))>tolerance)then
      message='fragment covariant operator group is not identity-first';return
    end if
    allocate(projected_scalars(n,n,nscalar),projected_vectors(n,n,3,nvector), &
      transformed(n,n),target(n,n))
    if(nop==1)then
      projected_scalars=scalars;projected_vectors=vectors
      pre_projection_defect=0d0;post_projection_defect=0d0;ok=.true.;return
    end if
    call covariance_defect(representation,rotations,scalars,vectors,pre_projection_defect)
    correction_limit=sqrt(tolerance)
    if(present(maximum_pre_projection_defect))correction_limit=maximum_pre_projection_defect
    if(.not.ieee_is_finite(correction_limit).or.correction_limit<=0d0)then
      message='fragment operator correction limit must be finite and positive';return
    end if
    if(pre_projection_defect>correction_limit)then
      message='fragment operator pre-projection covariance defect exceeds correction limit';return
    end if
    projected_scalars=(0d0,0d0);projected_vectors=(0d0,0d0)
    do iop=1,nop
      do i=1,nscalar
        projected_scalars(:,:,i)=projected_scalars(:,:,i)+matmul(conjg(transpose( &
          representation(:,:,iop))),matmul(scalars(:,:,i),representation(:,:,iop)))
      end do
      do i=1,nvector
        do a=1,3
          target=(0d0,0d0)
          do b=1,3
            transformed=matmul(conjg(transpose(representation(:,:,iop))), &
              matmul(vectors(:,:,b,i),representation(:,:,iop)))
            target=target+rotations(b,a,iop)*transformed
          end do
          projected_vectors(:,:,a,i)=projected_vectors(:,:,a,i)+target
        end do
      end do
    end do
    projected_scalars=projected_scalars/real(nop,8);projected_vectors=projected_vectors/real(nop,8)
    call covariance_defect(representation,rotations,projected_scalars,projected_vectors, &
      post_projection_defect)
    if(post_projection_defect>tolerance)then
      message='fragment operator post-projection covariance exceeds tolerance';return
    end if
    ok=.true.
  end subroutine project_dg_fragment_covariant_operators

  subroutine covariance_defect(representation,rotations,scalars,vectors,defect)
    complex(8),intent(in)::representation(:,:,:),scalars(:,:,:),vectors(:,:,:,:)
    real(8),intent(in)::rotations(:,:,:)
    real(8),intent(out)::defect
    complex(8),allocatable::transformed(:,:),target(:,:)
    integer::n,iop,i,a,b
    real(8)::scale
    n=size(representation,1);allocate(transformed(n,n),target(n,n));defect=0d0
    do iop=1,size(representation,3)
      do i=1,size(scalars,3)
        scale=max(1d0,maxval(abs(scalars(:,:,i))))
        transformed=matmul(conjg(transpose(representation(:,:,iop))), &
          matmul(scalars(:,:,i),representation(:,:,iop)))
        defect=max(defect,maxval(abs(transformed-scalars(:,:,i)))/scale)
      end do
      do i=1,size(vectors,4)
        scale=max(1d0,maxval(abs(vectors(:,:,:,i))))
        do a=1,3
        transformed=matmul(conjg(transpose(representation(:,:,iop))), &
          matmul(vectors(:,:,a,i),representation(:,:,iop)))
        target=(0d0,0d0)
        do b=1,3;target=target+rotations(a,b,iop)*vectors(:,:,b,i);end do
        defect=max(defect,maxval(abs(transformed-target))/scale)
        end do
      end do
    end do
  end subroutine covariance_defect

  function complex_identity(n)result(identity)
    integer,intent(in)::n
    complex(8)::identity(n,n)
    integer::i
    identity=(0d0,0d0);do i=1,n;identity(i,i)=1d0;end do
  end function complex_identity

  function identity3()result(identity)
    real(8)::identity(3,3)
    identity=0d0;identity(1,1)=1d0;identity(2,2)=1d0;identity(3,3)=1d0
  end function identity3

  real(8) function determinant3(matrix)
    real(8),intent(in)::matrix(3,3)
    determinant3=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))- &
      matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))+ &
      matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
  end function determinant3

  subroutine build_dg_fragment_group_representation(metric,raw,product_table,tolerance, &
      representation,raw_unitarity_defect,unitarity_defect,closure_defect,ok,message,&
      maximum_raw_unitarity_defect)
    complex(8),intent(in)::metric(:,:),raw(:,:,:)
    integer,intent(in)::product_table(:,:)
    real(8),intent(in)::tolerance
    complex(8),allocatable,intent(out)::representation(:,:,:)
    real(8),intent(out)::raw_unitarity_defect,unitarity_defect,closure_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),intent(in),optional::maximum_raw_unitarity_defect
    complex(8),allocatable::metric_sqrt(:,:),metric_inverse_sqrt(:,:),gram_inverse_sqrt(:,:), &
      transformed(:,:),unitary(:,:),difference(:,:),orthogonal_representation(:,:,:),&
      synchronized_representation(:,:,:),initial_orthogonal_representation(:,:,:)
    real(8)::metric_scale,representation_scale,correction_limit,synchronization_change,&
      synchronization_correction
    integer::n,nop,iop,jop,kop,iteration
    logical::power_ok,synchronization_converged
    ok=.false.;message='';raw_unitarity_defect=huge(1d0);unitarity_defect=huge(1d0)
    closure_defect=huge(1d0);n=size(metric,1);nop=size(raw,3)
    if(n<1.or.size(metric,2)/=n.or.size(raw,1)/=n.or.size(raw,2)/=n.or.nop<1.or. &
        size(product_table,1)/=nop.or.size(product_table,2)/=nop)then
      message='fragment group representation arrays have inconsistent dimensions';return
    end if
    if(.not.ieee_is_finite(tolerance).or.tolerance<=0d0)then
      message='fragment group representation tolerance must be finite and positive';return
    end if
    if(any(product_table<1).or.any(product_table>nop))then
      message='fragment group representation product table has an invalid index';return
    end if
    metric_scale=max(1d0,maxval(abs(metric)))
    if(maxval(abs(metric-conjg(transpose(metric))))>tolerance*metric_scale)then
      message='fragment group representation metric is not Hermitian';return
    end if
    call hermitian_matrix_power(metric,0.5d0,tolerance,metric_sqrt,power_ok)
    if(.not.power_ok)then;message='fragment group representation metric is not positive definite';return;end if
    call hermitian_matrix_power(metric,-0.5d0,tolerance,metric_inverse_sqrt,power_ok)
    if(.not.power_ok)then;message='fragment group representation metric inverse failed';return;end if
    allocate(representation(n,n,nop),transformed(n,n),unitary(n,n),difference(n,n),&
      orthogonal_representation(n,n,nop),synchronized_representation(n,n,nop),&
      initial_orthogonal_representation(n,n,nop))
    raw_unitarity_defect=0d0
    do iop=1,nop
      difference=matmul(conjg(transpose(raw(:,:,iop))),matmul(metric,raw(:,:,iop)))-metric
      raw_unitarity_defect=max(raw_unitarity_defect,maxval(abs(difference))/metric_scale)
      transformed=matmul(metric_sqrt,matmul(raw(:,:,iop),metric_inverse_sqrt))
      call hermitian_matrix_power(matmul(conjg(transpose(transformed)),transformed),-0.5d0, &
        tolerance,gram_inverse_sqrt,power_ok)
      if(.not.power_ok)then;message='fragment group representation polar correction is singular';return;end if
      unitary=matmul(transformed,gram_inverse_sqrt)
      orthogonal_representation(:,:,iop)=unitary
    end do
    correction_limit=sqrt(tolerance)
    if(present(maximum_raw_unitarity_defect))correction_limit=maximum_raw_unitarity_defect
    if(.not.ieee_is_finite(correction_limit).or.correction_limit<=0d0)then
      message='fragment group representation correction limit must be finite and positive';return
    end if
    if(raw_unitarity_defect>correction_limit)then
      message='fragment group representation raw unitarity defect exceeds correction limit';return
    end if
    initial_orthogonal_representation=orthogonal_representation
    synchronization_converged=.false.
    do iteration=1,64
      synchronized_representation=(0d0,0d0)
      do iop=1,nop;do jop=1,nop
        kop=product_table(iop,jop)
        synchronized_representation(:,:,iop)=synchronized_representation(:,:,iop)+&
          matmul(orthogonal_representation(:,:,kop),&
          conjg(transpose(orthogonal_representation(:,:,jop))))/real(nop,8)
      end do;end do
      do iop=1,nop
        call hermitian_matrix_power(matmul(conjg(transpose(synchronized_representation(:,:,iop))),&
          synchronized_representation(:,:,iop)),-0.5d0,tolerance,gram_inverse_sqrt,power_ok)
        if(.not.power_ok)then;message='fragment group synchronization polar factor is singular';return;end if
        synchronized_representation(:,:,iop)=matmul(synchronized_representation(:,:,iop),gram_inverse_sqrt)
        if(.not.all(ieee_is_finite(real(synchronized_representation(:,:,iop)))).or.&
            .not.all(ieee_is_finite(aimag(synchronized_representation(:,:,iop)))))then
          message='fragment group synchronization produced a nonfinite representation';return
        end if
      end do
      synchronization_change=maxval(abs(synchronized_representation-orthogonal_representation))
      orthogonal_representation=synchronized_representation
      if(synchronization_change<=0.1d0*tolerance)then
        synchronization_converged=.true.;exit
      end if
    end do
    if(.not.synchronization_converged)then
      message='fragment group synchronization did not converge';return
    end if
    synchronization_correction=maxval(abs(orthogonal_representation-initial_orthogonal_representation))
    if(synchronization_correction>correction_limit)then
      message='fragment group synchronization correction exceeds correction limit';return
    end if
    do iop=1,nop
      representation(:,:,iop)=matmul(metric_inverse_sqrt,&
        matmul(orthogonal_representation(:,:,iop),metric_sqrt))
    end do
    unitarity_defect=0d0
    do iop=1,nop
      difference=matmul(conjg(transpose(representation(:,:,iop))), &
        matmul(metric,representation(:,:,iop)))-metric
      unitarity_defect=max(unitarity_defect,maxval(abs(difference))/metric_scale)
    end do
    representation_scale=max(1d0,maxval(abs(representation)))
    closure_defect=0d0
    do iop=1,nop;do jop=1,nop
      kop=product_table(iop,jop)
      difference=matmul(representation(:,:,iop),representation(:,:,jop))-representation(:,:,kop)
      closure_defect=max(closure_defect,maxval(abs(difference))/representation_scale)
    end do;end do
    if(unitarity_defect>tolerance)then
      message='fragment group representation metric unitarity exceeds tolerance';return
    end if
    if(closure_defect>tolerance)then
      message='fragment group representation closure exceeds tolerance';return
    end if
    ok=.true.
  end subroutine build_dg_fragment_group_representation

  subroutine hermitian_matrix_power(matrix,power,tolerance,result,ok)
    complex(8),intent(in)::matrix(:,:)
    real(8),intent(in)::power,tolerance
    complex(8),allocatable,intent(out)::result(:,:)
    logical,intent(out)::ok
    complex(8),allocatable::vectors(:,:),work(:)
    real(8),allocatable::eigenvalues(:),rwork(:)
    complex(8)::work_query(1)
    integer::n,info,lwork,i
    n=size(matrix,1);ok=.false.
    allocate(vectors(n,n),eigenvalues(n),rwork(max(1,3*n-2)))
    vectors=matrix
    call zheev('V','U',n,vectors,n,eigenvalues,work_query,-1,rwork,info)
    if(info/=0)return
    lwork=max(1,int(real(work_query(1))));allocate(work(lwork));vectors=matrix
    call zheev('V','U',n,vectors,n,eigenvalues,work,lwork,rwork,info)
    if(info/=0.or.minval(eigenvalues)<=tolerance*max(1d0,maxval(abs(eigenvalues))))return
    allocate(result(n,n));result=(0d0,0d0)
    do i=1,n
      result=result+eigenvalues(i)**power*spread(vectors(:,i),2,n)* &
        spread(conjg(vectors(:,i)),1,n)
    end do
    ok=.true.
  end subroutine hermitian_matrix_power

  subroutine select_dg_exact_fragment_subgroup(product_table,atom_residual,boundary_residual, &
      grid_residual,center_residual,atom_tolerance,boundary_tolerance,grid_tolerance, &
      center_tolerance,subgroup,ok,message)
    integer,intent(in)::product_table(:,:)
    real(8),intent(in)::atom_residual(:),boundary_residual(:),grid_residual(:),center_residual(:)
    real(8),intent(in)::atom_tolerance,boundary_tolerance,grid_tolerance,center_tolerance
    integer,allocatable,intent(out)::subgroup(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,allocatable::accepted(:),groups(:,:),candidate(:)
    integer::n,igroup,iop,jop,ngroups,best_group,best_size,candidate_size
    logical::has_inverse

    ok=.false.;message='';n=size(product_table,1)
    if(n<1.or.size(product_table,2)/=n.or.size(atom_residual)/=n.or. &
        size(boundary_residual)/=n.or.size(grid_residual)/=n.or.size(center_residual)/=n)then
      message='exact fragment subgroup arrays have inconsistent dimensions';return
    end if
    if(n>maximum_crystallographic_point_group_order)then
      message='exact fragment subgroup exceeds crystallographic point-group order 48';return
    end if
    if(any(product_table<1).or.any(product_table>n))then
      message='exact fragment subgroup product table has an invalid index';return
    end if
    do iop=1,n
      if(product_table(1,iop)/=iop.or.product_table(iop,1)/=iop)then
        message='exact fragment subgroup operation 1 is not identity';return
      end if
      has_inverse=.false.
      do jop=1,n
        if(product_table(iop,jop)==1.and.product_table(jop,iop)==1)then
          has_inverse=.true.;exit
        end if
      end do
      if(.not.has_inverse)then
        message='exact fragment subgroup product table is not a group';return
      end if
    end do
    if(.not.valid_tolerance(atom_tolerance).or..not.valid_tolerance(boundary_tolerance).or. &
        .not.valid_tolerance(grid_tolerance).or..not.valid_tolerance(center_tolerance))then
      message='exact fragment subgroup tolerances must be finite and positive';return
    end if
    if(.not.valid_residuals(atom_residual).or..not.valid_residuals(boundary_residual).or. &
        .not.valid_residuals(grid_residual).or..not.valid_residuals(center_residual))then
      message='exact fragment subgroup residuals must be finite and nonnegative';return
    end if
    allocate(accepted(n))
    accepted=atom_residual<=atom_tolerance.and.boundary_residual<=boundary_tolerance.and. &
      grid_residual<=grid_tolerance.and.center_residual<=center_tolerance
    if(.not.accepted(1))then
      message='exact fragment subgroup identity operation failed a mapping tolerance';return
    end if

    allocate(groups(n,1));groups=.false.;groups(1,1)=.true.;ngroups=1
    igroup=1
    do while(igroup<=ngroups)
      do iop=2,n
        if(.not.accepted(iop).or.groups(iop,igroup))cycle
        candidate=groups(:,igroup);candidate(iop)=.true.
        call close_generated_subgroup(product_table,candidate)
        if(any(candidate.and..not.accepted))cycle
        if(.not.subgroup_is_known(groups,ngroups,candidate))then
          if(ngroups>=maximum_enumerated_subgroups)then
            message='exact fragment subgroup enumeration limit exceeded';return
          end if
          call append_subgroup(groups,ngroups,candidate)
        end if
      end do
      igroup=igroup+1
    end do

    best_group=1;best_size=1
    do igroup=2,ngroups
      candidate_size=count(groups(:,igroup))
      if(candidate_size>best_size.or.(candidate_size==best_size.and. &
          subgroup_lexically_precedes(groups(:,igroup),groups(:,best_group))))then
        best_group=igroup;best_size=candidate_size
      end if
    end do
    allocate(subgroup(best_size));subgroup=pack([(iop,iop=1,n)],groups(:,best_group))
    ok=.true.
  end subroutine select_dg_exact_fragment_subgroup

  logical function valid_tolerance(value)
    real(8),intent(in)::value
    valid_tolerance=ieee_is_finite(value).and.value>0d0
  end function valid_tolerance

  logical function valid_residuals(values)
    real(8),intent(in)::values(:)
    valid_residuals=all(ieee_is_finite(values)).and.all(values>=0d0)
  end function valid_residuals

  subroutine close_generated_subgroup(product_table,mask)
    integer,intent(in)::product_table(:,:)
    logical,intent(inout)::mask(:)
    integer::i,j
    logical::changed
    changed=.true.
    do while(changed)
      changed=.false.
      do i=1,size(mask)
        if(.not.mask(i))cycle
        do j=1,size(mask)
          if(.not.mask(j))cycle
          if(.not.mask(product_table(i,j)))then
            mask(product_table(i,j))=.true.;changed=.true.
          end if
        end do
      end do
    end do
  end subroutine close_generated_subgroup

  logical function subgroup_is_known(groups,ngroups,candidate)
    logical,intent(in)::groups(:,:),candidate(:)
    integer,intent(in)::ngroups
    integer::i
    subgroup_is_known=.false.
    do i=1,ngroups
      if(all(groups(:,i).eqv.candidate))then
        subgroup_is_known=.true.;return
      end if
    end do
  end function subgroup_is_known

  subroutine append_subgroup(groups,ngroups,candidate)
    logical,allocatable,intent(inout)::groups(:,:)
    integer,intent(inout)::ngroups
    logical,intent(in)::candidate(:)
    logical,allocatable::expanded(:,:)
    allocate(expanded(size(groups,1),ngroups+1))
    expanded(:,1:ngroups)=groups(:,1:ngroups);expanded(:,ngroups+1)=candidate
    call move_alloc(expanded,groups);ngroups=ngroups+1
  end subroutine append_subgroup

  logical function subgroup_lexically_precedes(left,right)
    logical,intent(in)::left(:),right(:)
    integer::i
    subgroup_lexically_precedes=.false.
    do i=1,size(left)
      if(left(i).eqv.right(i))cycle
      subgroup_lexically_precedes=left(i)
      return
    end do
  end function subgroup_lexically_precedes
end module dg_overlapping_wannier_symmetry
