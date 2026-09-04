module dg_hybrid_projected_fragment_pipeline
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_get_halting_mode,ieee_set_halting_mode,&
    ieee_set_flag,ieee_invalid,ieee_divide_by_zero,ieee_overflow
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog
  use dg_hybrid_windowed_pw_basis,only:materialize_dg_hybrid_windowed_pw_columns
  use dg_hybrid_wannier_complement,only:s_dg_hybrid_generalized_metric_factor,&
    prepare_dg_hybrid_generalized_wannier_metric,apply_dg_hybrid_generalized_wannier_projection_tile,&
    materialize_dg_hybrid_projected_pw_tile,compute_dg_hybrid_union_to_complete_binding
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_fragment_basis_stream,only:s_dg_hybrid_fragment_basis_stream,&
    initialize_dg_hybrid_fragment_basis_stream,append_dg_hybrid_projected_pw_tile,&
    finalize_dg_hybrid_fragment_basis_stream
  implicit none
  private
  type,public::s_dg_hybrid_support_samples
    ! Rows are required boundary samples, derivative samples or projector
    ! overlaps; columns are active basis functions / original DC states.
    complex(real64),allocatable::basis(:,:),reference(:,:)
    real(real64),allocatable::weights(:)
  end type
  type,public::s_dg_hybrid_core_projection_report
    logical::measured=.false.
    integer::metric_rank=0
    real(real64)::orbital_residual=huge(1d0),density_defect=huge(1d0),electron_defect=huge(1d0)
  end type
  type,public::s_dg_hybrid_projection_factorization_receipt
    logical::valid=.false.
    integer::basis_generation=0,metric_rank=0,metric_factorization_count=0,projected_tile_count=0
    integer(int64)::wannier_fingerprint=0_int64,metric_fingerprint=0_int64
  end type s_dg_hybrid_projection_factorization_receipt
  type,public::s_dg_hybrid_dual_basis_catalog
    logical::valid=.false.
    type(s_dg_hybrid_fragment_basis),allocatable::fragment_bases(:)
    complex(real64),allocatable::union_to_complete(:,:)
    integer::uncompressed_rank=0,complete_rank=0
    integer(int64)::fragment_catalog_fingerprint=0_int64,complete_map_fingerprint=0_int64,&
      complete_transform_binding_fingerprint=0_int64
    integer(int64),allocatable::uncompressed_global_basis_ids(:)
    integer,allocatable::uncompressed_owner_ranks(:),uncompressed_fragment_ids(:),&
      uncompressed_local_slots(:),uncompressed_sectors(:),uncompressed_generations(:)
  end type s_dg_hybrid_dual_basis_catalog
  public::build_dg_hybrid_projected_fragment_basis,finalize_dg_hybrid_dual_basis_catalog
  public::build_dg_hybrid_projected_local_fragment_basis
  public::project_dg_hybrid_core_seeds
  public::check_dg_hybrid_seed_support
  interface
    subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
      import::real64
      character(1),intent(in)::jobz,uplo
      integer,intent(in)::n,lda,lwork
      complex(real64),intent(inout)::a(lda,*),work(*)
      real(real64),intent(out)::w(*),rwork(*)
      integer,intent(out)::info
    end subroutine zheev
  end interface
contains
  ! Numerical admission of supplied support evidence, not a proof of inventory
  ! completeness. The operator adapter must supply required_counts independently
  ! of these arrays and bind sample identities/order to the same raw reference.
  ! Channel order: boundary, derivative, nonlocal projector. Each tolerance is
  ! an absolute weighted L2 orbital error in that channel's physical units,
  ! matching the raw seed span norm convention. Every seed, not only occupied
  ! seeds, is checked. This routine neither changes C nor publishes solver state.
  subroutine check_dg_hybrid_seed_support(comm,coefficients,samples,required_counts,tolerances,&
      selected_count,pw_cutoff,defects,measured,ok,message)
    integer,intent(in)::comm,required_counts(3),selected_count
    complex(real64),intent(in)::coefficients(:,:)
    type(s_dg_hybrid_support_samples),intent(in)::samples(3)
    real(real64),intent(in)::tolerances(3),pw_cutoff
    real(real64),intent(out)::defects(3)
    logical,intent(out)::measured,ok
    character(*),intent(out)::message
    logical::halting(3)
    call ieee_get_halting_mode(ieee_invalid,halting(1))
    call ieee_get_halting_mode(ieee_divide_by_zero,halting(2))
    call ieee_get_halting_mode(ieee_overflow,halting(3))
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.)
    call execute()
    call ieee_set_flag(ieee_invalid,.false.);call ieee_set_flag(ieee_divide_by_zero,.false.)
    call ieee_set_flag(ieee_overflow,.false.)
    call ieee_set_halting_mode(ieee_invalid,halting(1))
    call ieee_set_halting_mode(ieee_divide_by_zero,halting(2))
    call ieee_set_halting_mode(ieee_overflow,halting(3))
  contains
    subroutine execute()
      complex(real64),allocatable::reconstructed(:,:)
      real(real64)::metadata(4),minimum(4),maximum(4),error2
      integer::channel,n,m,j,status,ierr
      logical::valid
      character(256)::why
      ok=.false.;measured=.false.;message='';defects=huge(1d0)
      n=size(coefficients,1);m=size(coefficients,2)
      valid=n>0.and.m>0.and.selected_count>0.and.selected_count<=n.and.&
        finite_complex_matrix(coefficients).and.all(required_counts>=0).and.&
        all(ieee_is_finite(tolerances)).and.all(tolerances>0d0).and.&
        ieee_is_finite(pw_cutoff).and.pw_cutoff>=0d0
      do channel=1,3
        if(.not.allocated(samples(channel)%basis).or..not.allocated(samples(channel)%reference).or.&
            .not.allocated(samples(channel)%weights))then
          valid=.false.;cycle
        endif
        valid=valid.and.all(shape(samples(channel)%basis)==[required_counts(channel),n]).and.&
          all(shape(samples(channel)%reference)==[required_counts(channel),m]).and.&
          size(samples(channel)%weights)==required_counts(channel).and.&
          finite_complex_matrix(samples(channel)%basis).and.finite_complex_matrix(samples(channel)%reference).and.&
          all(ieee_is_finite(samples(channel)%weights)).and.all(samples(channel)%weights>0d0)
      enddo
      call synchronize_status(comm,valid,'invalid or incomplete required support evidence',ok,message)
      if(.not.ok)return
      metadata=[tolerances,pw_cutoff]
      call MPI_Allreduce(metadata,minimum,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      call synchronize_status(comm,ierr==MPI_SUCCESS,'support controls exchange failed',ok,message)
      if(.not.ok)return
      call MPI_Allreduce(metadata,maximum,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      call synchronize_status(comm,ierr==MPI_SUCCESS.and.all(minimum==maximum),&
        'support controls differ between ranks',ok,message);if(.not.ok)return
      defects=0d0
      do channel=1,3
        allocate(reconstructed(required_counts(channel),m),stat=status)
        call synchronize_status(comm,status==0,'support reconstruction allocation failed',ok,message)
        if(.not.ok)return
        reconstructed=matmul(samples(channel)%basis,coefficients)
        valid=finite_complex_matrix(reconstructed)
        do j=1,m
          error2=sum(samples(channel)%weights*abs(reconstructed(:,j)-samples(channel)%reference(:,j))**2)
          valid=valid.and.ieee_is_finite(error2)
          defects(channel)=max(defects(channel),sqrt(error2))
        enddo
        call synchronize_status(comm,valid,'nonfinite required support reconstruction',ok,message)
        if(.not.ok)return
        deallocate(reconstructed)
      enddo
      measured=.true.
      write(why,'(a,i0,a,es12.4,a,3es12.4)')'required support mismatch: selected=',selected_count,&
        ' cutoff=',pw_cutoff,' boundary/derivative/projector=',defects
      call synchronize_status(comm,all(defects<=tolerances),why,ok,message)
    end subroutine
  end subroutine check_dg_hybrid_seed_support

  ! Local core quadrature only; no complete-system Hamiltonian diagonalization.
  ! limits = metric rank, max relative orbital norm, relative density L1,
  ! absolute electron-number defect. A rank loss is never compressed away.
  subroutine project_dg_hybrid_core_seeds(comm,basis,weights,seeds,occupations,limits,&
      selected_count,pw_cutoff,coefficients,report,ok,message)
    integer,intent(in)::comm,selected_count
    complex(real64),intent(in)::basis(:,:),seeds(:,:)
    real(real64),intent(in)::weights(:),occupations(:),limits(4),pw_cutoff
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    type(s_dg_hybrid_core_projection_report),intent(out)::report
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical::halting(3)
    call ieee_get_halting_mode(ieee_invalid,halting(1))
    call ieee_get_halting_mode(ieee_divide_by_zero,halting(2))
    call ieee_get_halting_mode(ieee_overflow,halting(3))
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.)
    call execute()
    call ieee_set_flag(ieee_invalid,.false.);call ieee_set_flag(ieee_divide_by_zero,.false.)
    call ieee_set_flag(ieee_overflow,.false.)
    call ieee_set_halting_mode(ieee_invalid,halting(1))
    call ieee_set_halting_mode(ieee_divide_by_zero,halting(2))
    call ieee_set_halting_mode(ieee_overflow,halting(3))
  contains
    subroutine execute()
      complex(real64),allocatable::weighted(:,:),target(:,:),gram(:,:),inverse(:,:),work(:,:),reconstructed(:,:)
      real(real64),allocatable::spectrum(:),rho(:),reference(:)
      real(real64)::minimum(5),maximum(5),metadata(5),norm2,error2,electrons
      integer::n,m,p,j,status,ierr
      logical::valid,stage_ok
      character(256)::why
      ok=.false.;message='';n=size(basis,2);m=size(seeds,2);p=size(basis,1)
      valid=n>0.and.n<=huge(0)/3.and.m>0.and.p>0.and.size(seeds,1)==p.and.size(weights)==p.and.&
        size(occupations)==m.and.selected_count>0.and.selected_count<=n.and.&
        finite_complex_matrix(basis).and.finite_complex_matrix(seeds).and.&
        all(ieee_is_finite(weights)).and.all(weights>0d0).and.&
        all(ieee_is_finite(occupations)).and.all(occupations>=0d0).and.&
        all(ieee_is_finite(limits)).and.all(limits>0d0).and.ieee_is_finite(pw_cutoff).and.pw_cutoff>=0d0
      call synchronize_status(comm,valid,'invalid core seed projection input',ok,message);if(.not.ok)return
      metadata=[limits,pw_cutoff]
      call MPI_Allreduce(metadata,minimum,5,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      call synchronize_status(comm,ierr==MPI_SUCCESS,'core projection controls exchange failed',ok,message)
      if(.not.ok)return
      call MPI_Allreduce(metadata,maximum,5,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      call synchronize_status(comm,ierr==MPI_SUCCESS.and.all(minimum==maximum),&
        'core projection controls differ between ranks',ok,message);if(.not.ok)return
      allocate(weighted(p,n),target(p,m),gram(n,n),work(n,m),reconstructed(p,m),rho(p),reference(p),stat=status)
      call synchronize_status(comm,status==0,'core projection allocation failed',ok,message);if(.not.ok)return
      weighted=basis*spread(sqrt(weights),2,n);target=seeds*spread(sqrt(weights),2,m)
      gram=matmul(conjg(transpose(weighted)),weighted)
      call hermitian_pseudoinverse(gram,limits(1),inverse,report%metric_rank,spectrum,stage_ok,why)
      valid=stage_ok.and.report%metric_rank==n
      call synchronize_status(comm,valid,'unresolved or dependent core metric: '//trim(why),ok,message)
      if(.not.ok)return
      work=matmul(inverse,matmul(conjg(transpose(weighted)),target))
      reconstructed=matmul(basis,work)
      call synchronize_status(comm,finite_complex_matrix(work).and.finite_complex_matrix(reconstructed),&
        'nonfinite core seed projection result',ok,message);if(.not.ok)return
      report%orbital_residual=0d0;rho=0d0;reference=0d0;valid=.true.
      do j=1,m
        norm2=sum(weights*abs(seeds(:,j))**2)
        error2=sum(weights*abs(reconstructed(:,j)-seeds(:,j))**2)
        valid=valid.and.norm2>tiny(1d0).and.ieee_is_finite(norm2).and.ieee_is_finite(error2)
        report%orbital_residual=max(report%orbital_residual,sqrt(error2/max(norm2,tiny(1d0))))
        rho=rho+occupations(j)*abs(reconstructed(:,j))**2
        reference=reference+occupations(j)*abs(seeds(:,j))**2
      enddo
      electrons=sum(weights*reference)
      report%density_defect=sum(weights*abs(rho-reference))/max(1d0,electrons)
      report%electron_defect=abs(sum(weights*(rho-reference)))
      valid=valid.and.all(ieee_is_finite([report%orbital_residual,report%density_defect,&
        report%electron_defect,electrons])).and.all(ieee_is_finite(rho)).and.all(ieee_is_finite(reference))
      call synchronize_status(comm,valid,'unresolved core projection diagnostics',ok,message);if(.not.ok)return
      report%measured=.true.
      valid=report%orbital_residual<=limits(2).and.report%density_defect<=limits(3).and.&
        report%electron_defect<=limits(4)
      write(why,'(a,i0,a,es12.4,a,3es12.4)')'insufficient core span: selected=',selected_count,&
        ' cutoff=',pw_cutoff,' orbital/density/electron=',report%orbital_residual,&
        report%density_defect,report%electron_defect
      call synchronize_status(comm,valid,why,ok,message);if(.not.ok)return
      call move_alloc(work,coefficients)
    end subroutine
  end subroutine project_dg_hybrid_core_seeds

  subroutine build_dg_hybrid_projected_local_fragment_basis(comm,global_point_count,fragment_count,&
      fragment_id,core_ids,weights,core_coordinates,core_windows,buffer_ids,local_wannier,&
      buffer_coordinates,buffer_windows,catalog,g_vectors,wannier_owner,tile_width,tolerance,&
      wannier_fingerprint,basis,workspace_peak_bytes,fingerprint,ok,message,basis_generation,projection_receipt)
    integer,intent(in)::comm,global_point_count,fragment_count,fragment_id,tile_width,wannier_owner(:)
    integer(int64),intent(in)::core_ids(:),buffer_ids(:),wannier_fingerprint
    real(real64),intent(in)::weights(:),core_coordinates(:,:),core_windows(:,:),buffer_coordinates(:,:),&
      buffer_windows(:,:),g_vectors(:,:),tolerance
    complex(real64),intent(in)::local_wannier(:,:)
    type(s_dg_hybrid_basis_catalog),intent(in)::catalog
    type(s_dg_hybrid_fragment_basis),intent(out)::basis
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,optional,intent(in)::basis_generation
    type(s_dg_hybrid_projection_factorization_receipt),optional,intent(out)::projection_receipt
    integer::rank,nrank,ierr,metadata(4),minimum(4),maximum(4),status,root,npoint,nwf,&
      first,width,i,j,k,slot,nowner
    integer,allocatable::fragments(:),point_counts(:),owner_min(:),owner_max(:),&
      columns(:),core_map(:),buffer_map(:),point_slot(:)
    integer(int64),allocatable::source_ids(:)
    complex(real64),allocatable::core_union(:,:),buffer_union(:,:),tile(:,:)
    integer(int64)::assembly_elements,assembly_bytes,stage_bytes,working_bytes,index_elements,index_bytes
    logical::valid,stage_ok
    character(256)::stage_message
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64
    call clear_fragment_basis(basis)
    if(present(projection_receipt))projection_receipt=s_dg_hybrid_projection_factorization_receipt()
    call MPI_Comm_rank(comm,rank,ierr)
    call MPI_Comm_size(comm,nrank,ierr)
    metadata=[global_point_count,fragment_count,size(wannier_owner),tile_width]
    call MPI_Allreduce(metadata,minimum,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='local WF metadata agreement failed';return;endif
    call MPI_Allreduce(metadata,maximum,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum/=maximum).or.any(minimum<1))then
      message='invalid or inconsistent local WF metadata';return
    endif
    valid=fragment_count==nrank.and.fragment_id>=1.and.fragment_id<=fragment_count.and.&
      size(local_wannier,1)==count(wannier_owner==fragment_id).and.&
      size(local_wannier,2)==size(buffer_ids).and.all(wannier_owner>=1).and.&
      all(wannier_owner<=fragment_count).and.all(buffer_ids>=1_int64).and.&
      all(buffer_ids<=int(global_point_count,int64)).and.all(core_ids>=1_int64).and.&
      all(core_ids<=int(global_point_count,int64)).and.&
      all(ieee_is_finite(real(local_wannier))).and.all(ieee_is_finite(aimag(local_wannier)))
    call synchronize_status(comm,valid,'invalid single-owner local WF support',stage_ok,stage_message)
    if(.not.stage_ok)then;message=stage_message;return;endif
    allocate(point_slot(global_point_count),stat=status)
    call synchronize_status(comm,status==0,'cannot allocate physical WF point lookup',stage_ok,stage_message)
    if(.not.stage_ok)then;message=stage_message;return;endif
    point_slot=0
    do i=1,size(buffer_ids)
      if(point_slot(buffer_ids(i))/=0)valid=.false.
      point_slot(buffer_ids(i))=i
    enddo
    call synchronize_status(comm,valid,'duplicate local WF physical support ID',stage_ok,stage_message)
    if(.not.stage_ok)then;message=stage_message;return;endif
    nowner=size(wannier_owner)
    allocate(fragments(nrank),point_counts(nrank),owner_min(nowner),owner_max(nowner),stat=status)
    call synchronize_status(comm,status==0,'cannot allocate local WF directory',stage_ok,stage_message)
    if(.not.stage_ok)then;message=stage_message;return;endif
    call MPI_Allreduce(wannier_owner,owner_min,nowner,MPI_INTEGER,MPI_MIN,comm,ierr)
    valid=ierr==MPI_SUCCESS
    call MPI_Allreduce(wannier_owner,owner_max,nowner,MPI_INTEGER,MPI_MAX,comm,ierr)
    valid=valid.and.ierr==MPI_SUCCESS.and.all(owner_min==owner_max)
    call MPI_Allgather(fragment_id,1,MPI_INTEGER,fragments,1,MPI_INTEGER,comm,ierr)
    valid=valid.and.ierr==MPI_SUCCESS
    call MPI_Allgather(size(buffer_ids),1,MPI_INTEGER,point_counts,1,MPI_INTEGER,comm,ierr)
    valid=valid.and.ierr==MPI_SUCCESS
    do i=1,fragment_count
      if(count(fragments==i)/=1)valid=.false.
    enddo
    call synchronize_status(comm,valid,'inconsistent local WF ownership directory',stage_ok,stage_message)
    if(.not.stage_ok)then;message=stage_message;return;endif
    ! Only the accepted physical buffer support is exported; elsewhere each WF is zero.
    ! This conversion does not certify omitted tails or alter the retained WF columns.
    allocate(core_union(nowner,size(core_ids)),buffer_union(nowner,size(buffer_ids)),&
      core_map(size(core_ids)),buffer_map(size(buffer_ids)),columns(nowner),stat=status)
    call synchronize_status(comm,status==0,'cannot allocate local WF union views',stage_ok,stage_message)
    if(.not.stage_ok)then;message=stage_message;return;endif
    core_union=(0d0,0d0);buffer_union=(0d0,0d0);working_bytes=0_int64
    do root=0,nrank-1
      npoint=point_counts(root+1);nwf=0
      do i=1,nowner
        if(wannier_owner(i)/=fragments(root+1))cycle
        nwf=nwf+1;columns(nwf)=i
      enddo
      allocate(source_ids(npoint),stat=status)
      call synchronize_status(comm,status==0,'cannot allocate WF support IDs',stage_ok,stage_message)
      if(.not.stage_ok)then;message=stage_message;return;endif
      if(rank==root)source_ids=buffer_ids
      call MPI_Bcast(source_ids,npoint,MPI_INTEGER8,root,comm,ierr)
      call synchronize_status(comm,ierr==MPI_SUCCESS,'WF support ID exchange failed',stage_ok,stage_message)
      if(.not.stage_ok)then;message=stage_message;return;endif
      point_slot=0
      do j=1,npoint
        point_slot(source_ids(j))=j
      enddo
      core_map=point_slot(core_ids);buffer_map=point_slot(buffer_ids)
      working_bytes=max(working_bytes,8_int64*npoint)
      do first=1,nwf,tile_width
        width=min(tile_width,nwf-first+1)
        valid=int(width,int64)*int(npoint,int64)<=int(huge(0),int64)
        call synchronize_status(comm,valid,'WF tile MPI count overflows',stage_ok,stage_message)
        if(.not.stage_ok)then;message=stage_message;return;endif
        allocate(tile(width,npoint),stat=status)
        call synchronize_status(comm,status==0,'cannot allocate WF support tile',stage_ok,stage_message)
        if(.not.stage_ok)then;message=stage_message;return;endif
        if(rank==root)tile=local_wannier(first:first+width-1,:)
        call MPI_Bcast(tile,size(tile),MPI_DOUBLE_COMPLEX,root,comm,ierr)
        call synchronize_status(comm,ierr==MPI_SUCCESS,'WF support tile exchange failed',stage_ok,stage_message)
        if(.not.stage_ok)then;message=stage_message;return;endif
        do k=1,width
          slot=columns(first+k-1)
          do i=1,size(core_ids)
            if(core_map(i)>0)core_union(slot,i)=tile(k,core_map(i))
          enddo
          do i=1,size(buffer_ids)
            if(buffer_map(i)>0)buffer_union(slot,i)=tile(k,buffer_map(i))
          enddo
        enddo
        working_bytes=max(working_bytes,16_int64*size(tile,kind=int64)+8_int64*npoint)
        deallocate(tile)
      enddo
      deallocate(source_ids)
    enddo
    valid=.true.;assembly_elements=size(core_union,kind=int64)
    call checked_add_nonnegative_int64(assembly_elements,size(buffer_union,kind=int64),valid)
    call checked_multiply_nonnegative_int64(assembly_elements,16_int64,assembly_bytes,valid)
    index_elements=2_int64*nrank+3_int64*nowner+size(core_map,kind=int64)+&
      size(buffer_map,kind=int64)+size(point_slot,kind=int64)
    call checked_multiply_nonnegative_int64(index_elements,int(storage_size(0)/8,int64),index_bytes,valid)
    call checked_add_nonnegative_int64(assembly_bytes,index_bytes,valid)
    call checked_add_nonnegative_int64(working_bytes,assembly_bytes,valid)
    call synchronize_status(comm,valid,'WF union workspace receipt overflows',stage_ok,stage_message)
    if(.not.stage_ok)then;message=stage_message;return;endif
    call build_dg_hybrid_projected_fragment_basis(comm,global_point_count,fragment_count,fragment_id,&
      core_ids,weights,core_union,core_coordinates,core_windows,buffer_ids,buffer_union,&
      buffer_coordinates,buffer_windows,catalog,g_vectors,wannier_owner,tile_width,tolerance,&
      wannier_fingerprint,basis,stage_bytes,fingerprint,ok,message,basis_generation,projection_receipt)
    if(.not.ok)return
    call checked_add_nonnegative_int64(assembly_bytes,stage_bytes,valid)
    call synchronize_status(comm,valid,'local WF pipeline workspace receipt overflows',stage_ok,stage_message)
    if(.not.stage_ok)then
      call clear_fragment_basis(basis)
      if(present(projection_receipt))projection_receipt=s_dg_hybrid_projection_factorization_receipt()
      ok=.false.;fingerprint=0_int64;message=stage_message;return
    endif
    workspace_peak_bytes=max(working_bytes,assembly_bytes)
  end subroutine build_dg_hybrid_projected_local_fragment_basis

  subroutine build_dg_hybrid_projected_fragment_basis(comm,global_point_count,fragment_count,&
      fragment_id,core_ids,weights,core_wannier,core_coordinates,core_windows,buffer_ids,&
      buffer_wannier,buffer_coordinates,buffer_windows,catalog,g_vectors,wannier_owner,tile_width,&
      tolerance,wannier_fingerprint,basis,workspace_peak_bytes,fingerprint,ok,message,&
      basis_generation,projection_receipt)
    integer,intent(in)::comm,global_point_count,fragment_count,fragment_id,tile_width,wannier_owner(:)
    integer(int64),intent(in)::core_ids(:),buffer_ids(:),wannier_fingerprint
    real(real64),intent(in)::weights(:),core_coordinates(:,:),core_windows(:,:),buffer_coordinates(:,:),&
      buffer_windows(:,:),g_vectors(:,:),tolerance
    complex(real64),intent(in)::core_wannier(:,:),buffer_wannier(:,:)
    type(s_dg_hybrid_basis_catalog),intent(in)::catalog
    type(s_dg_hybrid_fragment_basis),intent(out)::basis
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,intent(in),optional::basis_generation
    type(s_dg_hybrid_projection_factorization_receipt),intent(out),optional::projection_receipt
    type(s_dg_hybrid_fragment_basis_stream)::stream
    type(s_dg_hybrid_fragment_basis)::working_basis
    type(s_dg_hybrid_generalized_metric_factor)::metric_factor
    type(s_dg_hybrid_projection_factorization_receipt)::working_receipt
    integer::packet,g,npw,column,width,i,requested_generation,allocation_status,ierr,&
      minimum_generation,maximum_generation,nrank
    integer,allocatable::pw_owner(:)
    real(real64),allocatable::normalized_buffer_windows(:,:)
    complex(real64),allocatable::core_raw(:,:),buffer_raw(:,:),coefficients(:,:),projected_core(:,:),&
      projected_buffer(:,:)
    integer(int64)::stage_workspace,stage_fingerprint,stream_workspace,stream_fingerprint,&
      metric_workspace,metric_fingerprint,working_workspace,working_fingerprint,npw_count,&
      tile_elements,tile_bytes
    real(real64)::scale,orthogonality_defect
    logical::stage_ok,local_ok
    character(256)::stage_message,local_message
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64;npw=0
    working_workspace=0_int64;working_fingerprint=0_int64;npw_count=0_int64
    requested_generation=1;if(present(basis_generation))requested_generation=basis_generation
    call clear_fragment_basis(basis)
    working_receipt=s_dg_hybrid_projection_factorization_receipt()
    if(present(projection_receipt))projection_receipt=s_dg_hybrid_projection_factorization_receipt()
    call MPI_Allreduce(requested_generation,minimum_generation,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='projected basis-generation agreement failed';return;endif
    call MPI_Allreduce(requested_generation,maximum_generation,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_generation/=maximum_generation)then
      message='projected basis generation differs between ranks';return
    endif
    local_ok=.true.;local_message=''
    if(.not.catalog%valid.or.tile_width<1.or.requested_generation<1.or.&
        size(core_wannier,1)/=size(wannier_owner).or.&
        size(core_wannier,2)/=size(core_ids).or.size(buffer_wannier,1)/=size(wannier_owner).or.&
        size(buffer_wannier,2)/=size(buffer_ids).or.size(weights)/=size(core_ids).or.&
        any(shape(core_coordinates)/=[3,size(core_ids)]).or.any(shape(buffer_coordinates)/=[3,size(buffer_ids)]).or.&
        size(core_windows,2)/=size(core_ids).or.size(buffer_windows,2)/=size(buffer_ids).or.&
        size(core_windows,1)/=fragment_count.or.size(buffer_windows,1)/=fragment_count)then
      local_ok=.false.;local_message='invalid projected fragment pipeline shape or generation'
    endif
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    do packet=1,size(catalog%packets)
      call checked_add_nonnegative_int64(npw_count,&
        size(catalog%packets(packet)%g_indices,kind=int64),local_ok)
    enddo
    if(npw_count>int(huge(npw),int64))local_ok=.false.
    if(local_ok)npw=int(npw_count)
    local_ok=local_ok.and.npw>=1;local_message='empty or overflowing projected fragment PW catalog'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    allocate(pw_owner(npw),stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate projected PW ownership workspace'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    allocate(normalized_buffer_windows,source=buffer_windows,stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate normalized buffer-window workspace'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    column=0
    do packet=1,size(catalog%packets)
      if(catalog%packets(packet)%fragment_id<1.or.catalog%packets(packet)%fragment_id>fragment_count)then
        local_ok=.false.;local_message='invalid projected fragment packet owner';exit
      endif
      do g=1,size(catalog%packets(packet)%g_indices)
        column=column+1;pw_owner(column)=catalog%packets(packet)%fragment_id
      enddo
    enddo
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    do i=1,size(buffer_ids)
      scale=sqrt(sum(normalized_buffer_windows(:,i)**2))
      if(.not.ieee_is_finite(scale).or.scale<=0d0)then
        local_ok=.false.;local_message='fragment buffer lies outside all PW windows';exit
      endif
      normalized_buffer_windows(:,i)=normalized_buffer_windows(:,i)/scale
    enddo
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif

    call prepare_dg_hybrid_generalized_wannier_metric(comm,global_point_count,core_ids,weights,&
      core_wannier,wannier_fingerprint,tolerance,metric_factor,metric_workspace,metric_fingerprint,&
      stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    working_workspace=max(working_workspace,metric_workspace)
    call MPI_Comm_size(comm,nrank,ierr)
    if(nrank==fragment_count)then
      call initialize_dg_hybrid_fragment_basis_stream(comm,fragment_count,fragment_id,buffer_ids,&
        buffer_wannier(pack([(i,i=1,size(wannier_owner))],wannier_owner==fragment_id),:),&
        wannier_owner,pw_owner,stream,working_basis,stream_workspace,stream_fingerprint,&
        stage_ok,stage_message,local_wannier_only=.true.,basis_generation=requested_generation)
    else
      call initialize_dg_hybrid_fragment_basis_stream(comm,fragment_count,fragment_id,buffer_ids,&
        buffer_wannier,wannier_owner,pw_owner,stream,working_basis,stream_workspace,stream_fingerprint,&
        stage_ok,stage_message,basis_generation=requested_generation)
    endif
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    do column=1,npw,tile_width
      width=min(tile_width,npw-column+1)
      call materialize_dg_hybrid_windowed_pw_columns(catalog,g_vectors,core_coordinates,core_windows,&
        column,width,core_raw,stage_ok,stage_message)
      call synchronize_status(comm,stage_ok,stage_message,local_ok,local_message)
      if(.not.local_ok)then;message=trim(local_message);return;endif
      call apply_dg_hybrid_generalized_wannier_projection_tile(comm,core_ids,weights,core_wannier,&
        core_raw,catalog%packet_fingerprint,column,metric_factor,coefficients,projected_core,&
        orthogonality_defect,stage_workspace,stage_fingerprint,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      working_workspace=max(working_workspace,stage_workspace)
      call materialize_dg_hybrid_windowed_pw_columns(catalog,g_vectors,buffer_coordinates,&
        normalized_buffer_windows,column,width,buffer_raw,stage_ok,stage_message)
      call synchronize_status(comm,stage_ok,stage_message,local_ok,local_message)
      if(.not.local_ok)then;message=trim(local_message);return;endif
      if(size(buffer_ids)==0)then
        allocate(projected_buffer(width,0),stat=allocation_status)
        stage_ok=allocation_status==0;stage_message='cannot allocate idle projected PW buffer'
      else
        call materialize_dg_hybrid_projected_pw_tile(buffer_wannier,buffer_raw,coefficients,&
          projected_buffer,stage_ok,stage_message)
      endif
      call synchronize_status(comm,stage_ok,stage_message,local_ok,local_message)
      if(.not.local_ok)then;message=trim(local_message);return;endif
      call append_dg_hybrid_projected_pw_tile(stream,column,projected_buffer,working_basis,stage_ok,stage_message)
      call synchronize_status(comm,stage_ok,stage_message,local_ok,local_message)
      if(.not.local_ok)then;message=trim(local_message);return;endif
      tile_elements=0_int64;local_ok=.true.
      call checked_add_nonnegative_int64(tile_elements,size(core_raw,kind=int64),local_ok)
      call checked_add_nonnegative_int64(tile_elements,size(buffer_raw,kind=int64),local_ok)
      call checked_add_nonnegative_int64(tile_elements,size(coefficients,kind=int64),local_ok)
      call checked_add_nonnegative_int64(tile_elements,size(projected_core,kind=int64),local_ok)
      call checked_add_nonnegative_int64(tile_elements,size(projected_buffer,kind=int64),local_ok)
      call checked_multiply_nonnegative_int64(tile_elements,16_int64,tile_bytes,local_ok)
      call synchronize_status(comm,local_ok,'projected fragment workspace receipt overflow',stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      working_workspace=max(working_workspace,tile_bytes)
      working_fingerprint=ieor(working_fingerprint,ishftc(stage_fingerprint,modulo(column,63)))
      deallocate(core_raw,buffer_raw,coefficients,projected_core,projected_buffer)
      working_receipt%projected_tile_count=working_receipt%projected_tile_count+1
    enddo
    call finalize_dg_hybrid_fragment_basis_stream(comm,stream,working_basis,stream_workspace,&
      stream_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    working_workspace=max(working_workspace,stream_workspace)
    working_fingerprint=ieor(working_fingerprint,stream_fingerprint)
    if(working_fingerprint==0_int64)working_fingerprint=1_int64
    if(fragment_id==0)then
      working_basis%generation=0
      working_basis%provenance_fingerprint=0_int64
    else
      working_basis%generation=requested_generation
    endif
    call move_fragment_basis(working_basis,basis)
    working_receipt%valid=.true.;working_receipt%basis_generation=requested_generation
    working_receipt%metric_rank=metric_factor%metric_rank
    working_receipt%metric_factorization_count=1
    working_receipt%wannier_fingerprint=wannier_fingerprint
    working_receipt%metric_fingerprint=metric_fingerprint
    if(present(projection_receipt))projection_receipt=working_receipt
    workspace_peak_bytes=working_workspace;fingerprint=working_fingerprint
    deallocate(pw_owner,normalized_buffer_windows);ok=.true.
  end subroutine build_dg_hybrid_projected_fragment_basis

  subroutine finalize_dg_hybrid_dual_basis_catalog(comm,global_row_count,row_ids,weights,fragment_bases,&
      uncompressed_global_basis_ids,uncompressed_basis_values,union_to_complete,expected_fragment_ranks,&
      maximum_terminal_rank_loss,seed_fragment_owner,seed_coefficients_in_uncompressed,metric_tolerance,&
      required_interface_fragment_ids,required_interface_row_ids,required_periodic_wrap_fragment_ids,&
      required_periodic_wrap_row_ids,required_projector_fragment_ids,required_projector_row_ids,packet_ids,&
      packet_neighbor_offsets,packet_neighbor_basis_ids,required_neighbor_packet_ids,&
      required_neighbor_basis_ids,tail_tolerance,catalog,preserved_seed_fragment_owner,&
      preserved_seed_coefficients_in_uncompressed,complete_seed_coefficients,seed_reconstruction_defect,ok,message,&
      expected_fragment_wannier_ranks)
    integer,intent(in)::comm,global_row_count,expected_fragment_ranks(:),maximum_terminal_rank_loss,&
      seed_fragment_owner(:),required_interface_fragment_ids(:),required_periodic_wrap_fragment_ids(:),&
      required_projector_fragment_ids(:),packet_ids(:),packet_neighbor_offsets(:),required_neighbor_packet_ids(:)
    integer(int64),intent(in)::row_ids(:),uncompressed_global_basis_ids(:),required_interface_row_ids(:),&
      required_periodic_wrap_row_ids(:),required_projector_row_ids(:),packet_neighbor_basis_ids(:),&
      required_neighbor_basis_ids(:)
    real(real64),intent(in)::weights(:),metric_tolerance,tail_tolerance
    type(s_dg_hybrid_fragment_basis),intent(in)::fragment_bases(:)
    complex(real64),intent(in)::uncompressed_basis_values(:,:),union_to_complete(:,:),&
      seed_coefficients_in_uncompressed(:,:)
    type(s_dg_hybrid_dual_basis_catalog),intent(out)::catalog
    integer,allocatable,intent(out)::preserved_seed_fragment_owner(:)
    complex(real64),allocatable,intent(out)::preserved_seed_coefficients_in_uncompressed(:,:),&
      complete_seed_coefficients(:,:)
    real(real64),intent(out)::seed_reconstruction_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,intent(in),optional::expected_fragment_wannier_ranks(:)

    type(s_dg_hybrid_dual_basis_catalog)::working_catalog
    integer::rank,nproc,ierr,nfragment,nuncompressed,ncomplete,nseed,nlocal,b,j,p,f,slot,entry,&
      root,npoint,ncolumn,allocation_status,global_bad,gram_rank,complete_gram_rank,seed_rank,&
      expected_wannier_size,has_expected_wannier,publisher_value_count,gram_mpi_count
    integer::metadata(23),metadata_min(23),metadata_max(23),publisher_metadata(2)
    integer,allocatable::row_count(:),local_publishers(:),global_publishers(:),&
      local_id_count(:),global_id_count(:),local_owner_code(:),local_fragment(:),local_slot(:),&
      local_sector(:),local_generation(:),global_owner_code(:),global_fragment(:),global_slot(:),&
      global_sector(:),global_generation(:),working_seed_owner(:),fragment_slots(:)
    integer(int64),allocatable::publisher_point_ids(:),local_row_hash(:),global_row_hash(:)
    complex(real64),allocatable::publisher_values(:,:),local_gram(:,:),gram(:,:),gram_inverse(:,:),&
      complete_gram(:,:),complete_gram_inverse(:,:),seed_gram(:,:),seed_gram_inverse(:,:),rhs(:,:),&
      working_complete_seed(:,:),working_preserved_seed(:,:),residual(:),retained_projector(:,:),&
      metric_projector(:,:)
    real(real64),allocatable::gram_eigenvalues(:),complete_eigenvalues(:),seed_eigenvalues(:),&
      local_omitted_by_column(:),global_omitted_by_column(:)
    real(real64)::total_omitted,local_norm,seed_norm,residual_norm,working_seed_defect,column_defect,&
      orthogonality_tolerance,retained_span_defect,amplitude,weight_root,omitted_term
    integer(int64)::fragment_hash,map_hash,transform_binding_hash,provenance,bits,expected_rank_total,gram_count64,&
      publisher_value_count64
    logical::local_ok,stage_ok,identity_map,arithmetic_ok
    character(256)::local_message,stage_message

    ok=.false.;message='';seed_reconstruction_defect=0d0
    call clear_dual_catalog(catalog)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then;message='dual catalog communicator failure';return;endif
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)then;message='dual catalog communicator failure';return;endif
    nfragment=size(expected_fragment_ranks);nuncompressed=size(uncompressed_global_basis_ids)
    ncomplete=size(union_to_complete,2);nseed=size(seed_fragment_owner);nlocal=size(row_ids)
    has_expected_wannier=merge(1,0,present(expected_fragment_wannier_ranks));expected_wannier_size=0
    if(present(expected_fragment_wannier_ranks))expected_wannier_size=size(expected_fragment_wannier_ranks)
    metadata=[global_row_count,nfragment,nuncompressed,ncomplete,nseed,maximum_terminal_rank_loss,&
      size(fragment_bases),size(union_to_complete,1),size(seed_coefficients_in_uncompressed,1),&
      size(seed_coefficients_in_uncompressed,2),size(required_interface_fragment_ids),&
      size(required_interface_row_ids),size(required_periodic_wrap_fragment_ids),&
      size(required_periodic_wrap_row_ids),size(required_projector_fragment_ids),&
      size(required_projector_row_ids),size(required_neighbor_packet_ids),size(packet_ids),&
      size(packet_neighbor_offsets),size(packet_neighbor_basis_ids),size(required_neighbor_basis_ids),&
      has_expected_wannier,expected_wannier_size]
    call MPI_Allreduce(metadata,metadata_min,size(metadata),MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(metadata,metadata_max,size(metadata),MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(metadata_min/=metadata_max))then
      message='inconsistent dual catalog replicated dimensions';return
    endif
    expected_rank_total=0_int64;arithmetic_ok=.true.
    do f=1,nfragment
      call checked_add_nonnegative_int64(expected_rank_total,int(expected_fragment_ranks(f),int64),arithmetic_ok)
    enddo
    local_ok=arithmetic_ok.and.expected_rank_total==int(nuncompressed,int64).and.&
      global_row_count>0.and.nfragment>0.and.nuncompressed>0.and.ncomplete>0.and.nseed>0.and.&
      ncomplete<=nuncompressed.and.nseed<=ncomplete.and.&
      maximum_terminal_rank_loss>=0.and.size(fragment_bases)==nfragment.and.&
      size(union_to_complete,1)==nuncompressed.and.&
      all(shape(uncompressed_basis_values)==[nuncompressed,nlocal]).and.size(weights)==nlocal.and.&
      all(shape(seed_coefficients_in_uncompressed)==[nuncompressed,nseed]).and.&
      size(required_interface_fragment_ids)==size(required_interface_row_ids).and.&
      size(required_periodic_wrap_fragment_ids)==size(required_periodic_wrap_row_ids).and.&
      size(required_projector_fragment_ids)==size(required_projector_row_ids).and.&
      size(required_neighbor_packet_ids)==size(required_neighbor_basis_ids).and.&
      size(packet_neighbor_offsets)==size(packet_ids)+1
    local_ok=local_ok.and.ieee_is_finite(metric_tolerance).and.metric_tolerance>0d0.and.&
      ieee_is_finite(tail_tolerance).and.tail_tolerance>=0d0
    if(present(expected_fragment_wannier_ranks))then
      local_ok=local_ok.and.size(expected_fragment_wannier_ranks)==nfragment
      if(local_ok)local_ok=all(expected_fragment_wannier_ranks>=0).and.&
        all(expected_fragment_wannier_ranks<=expected_fragment_ranks)
    endif
    if(.not.local_ok)local_message='invalid dual catalog shape, rank, or tolerance contract'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif

    call validate_replicated_dual_inputs(comm,expected_fragment_ranks,uncompressed_global_basis_ids,&
      union_to_complete,seed_fragment_owner,seed_coefficients_in_uncompressed,metric_tolerance,&
      required_interface_fragment_ids,required_interface_row_ids,required_periodic_wrap_fragment_ids,&
      required_periodic_wrap_row_ids,required_projector_fragment_ids,required_projector_row_ids,packet_ids,&
      packet_neighbor_offsets,packet_neighbor_basis_ids,required_neighbor_packet_ids,required_neighbor_basis_ids,&
      tail_tolerance,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    if(present(expected_fragment_wannier_ranks))then
      stage_ok=replicated_integer_array(comm,expected_fragment_wannier_ranks)
      if(.not.stage_ok)then;message='expected fragment Wannier ranks differ between ranks';return;endif
    endif

    local_ok=.true.;local_message=''
    if(any(uncompressed_global_basis_ids<=0_int64))then
      local_ok=.false.;local_message='canonical Hybrid basis IDs must be positive'
    else
      do j=2,nuncompressed
        if(uncompressed_global_basis_ids(j)==uncompressed_global_basis_ids(j-1))then
          local_ok=.false.;local_message='duplicate canonical Hybrid basis ID';exit
        elseif(uncompressed_global_basis_ids(j)<uncompressed_global_basis_ids(j-1))then
          local_ok=.false.;local_message='canonical Hybrid basis ID order must be strictly increasing';exit
        endif
      enddo
    endif
    if(local_ok.and.(any(expected_fragment_ranks<1).or.any(seed_fragment_owner<1).or.&
        any(seed_fragment_owner>nfragment)))then
      local_ok=.false.;local_message='invalid fragment rank or seed owner metadata'
    endif
    if(local_ok.and.(.not.finite_complex_matrix(uncompressed_basis_values).or.&
        .not.finite_complex_matrix(union_to_complete).or.&
        .not.finite_complex_matrix(seed_coefficients_in_uncompressed)))then
      local_ok=.false.;local_message='nonfinite dual catalog numerical payload'
    endif
    if(local_ok.and.(.not.all(ieee_is_finite(weights)).or.any(weights<=0d0)))then
      local_ok=.false.;local_message='dual catalog weights must be positive and finite'
    endif
    if(local_ok.and.nuncompressed-ncomplete>maximum_terminal_rank_loss)then
      local_ok=.false.;local_message='terminal rank loss exceeds declared maximum'
    endif
    if(local_ok.and.(any(row_ids<1_int64).or.any(row_ids>int(global_row_count,int64))))then
      local_ok=.false.;local_message='dual catalog spatial row ID is out of range'
    endif
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif

    allocate(row_count(global_row_count),stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate dual catalog row ownership workspace'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    row_count=0
    do p=1,nlocal
      row_count(int(row_ids(p)))=row_count(int(row_ids(p)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,row_count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    local_ok=ierr==MPI_SUCCESS.and.all(row_count==1)
    local_message='duplicate or missing dual catalog spatial row'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif

    allocate(local_publishers(nfragment),global_publishers(nfragment),local_id_count(nuncompressed),&
      global_id_count(nuncompressed),local_owner_code(nuncompressed),local_fragment(nuncompressed),&
      local_slot(nuncompressed),local_sector(nuncompressed),local_generation(nuncompressed),&
      global_owner_code(nuncompressed),global_fragment(nuncompressed),global_slot(nuncompressed),&
      global_sector(nuncompressed),global_generation(nuncompressed),stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate dual catalog ownership tuple workspace'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    local_publishers=0;local_id_count=0;local_owner_code=0;local_fragment=0;local_slot=0
    local_sector=0;local_generation=0;local_ok=.true.;local_message=''
    do b=1,size(fragment_bases)
      f=fragment_bases(b)%fragment_id
      if(f==0)then
        if(fragment_bases(b)%generation/=0.or.fragment_bases(b)%provenance_fingerprint/=0_int64.or.&
            allocated_nonempty_fragment_payload(fragment_bases(b)))then
          local_ok=.false.;local_message='idle rank published a fragment basis payload';exit
        endif
        cycle
      endif
      if(f<1.or.f>nfragment)then
        local_ok=.false.;local_message='invalid fragment publisher ID';exit
      endif
      local_publishers(f)=local_publishers(f)+1
      if(.not.allocated(fragment_bases(b)%global_ids).or..not.allocated(fragment_bases(b)%sector).or.&
          .not.allocated(fragment_bases(b)%buffer_point_ids).or.&
          .not.allocated(fragment_bases(b)%buffer_values))then
        local_ok=.false.;local_message='fragment publisher basis payload is missing';exit
      endif
      if(size(fragment_bases(b)%global_ids)/=expected_fragment_ranks(f).or.&
          size(fragment_bases(b)%sector)/=size(fragment_bases(b)%global_ids))then
        local_ok=.false.;local_message='fragment publisher rank does not match expected rank';exit
      endif
      if(any(shape(fragment_bases(b)%buffer_values)/=&
          [size(fragment_bases(b)%buffer_point_ids),size(fragment_bases(b)%global_ids)]))then
        local_ok=.false.;local_message='fragment publisher buffer shape mismatch';exit
      endif
      if(fragment_bases(b)%generation<1.or.fragment_bases(b)%provenance_fingerprint==0_int64)then
        local_ok=.false.;local_message='invalid fragment publisher generation or provenance';exit
      endif
      if(any(fragment_bases(b)%sector<1).or.any(fragment_bases(b)%sector>2).or.&
          any(fragment_bases(b)%buffer_point_ids<1_int64).or.&
          any(fragment_bases(b)%buffer_point_ids>int(global_row_count,int64)).or.&
          .not.finite_complex_matrix(fragment_bases(b)%buffer_values))then
        local_ok=.false.;local_message='invalid fragment publisher sector, row, or value payload';exit
      endif
      if(present(expected_fragment_wannier_ranks))then
        if(count(fragment_bases(b)%sector==1)/=expected_fragment_wannier_ranks(f).or.&
            count(fragment_bases(b)%sector==2)/=&
              expected_fragment_ranks(f)-expected_fragment_wannier_ranks(f))then
          local_ok=.false.;local_message='fragment publisher WF/PW sector cardinality mismatch';exit
        endif
      endif
      do p=1,size(fragment_bases(b)%buffer_point_ids)
        if(any(fragment_bases(b)%buffer_point_ids(p+1:)==fragment_bases(b)%buffer_point_ids(p)))then
          local_ok=.false.;local_message='duplicate fragment publisher buffer row';exit
        endif
      enddo
      if(.not.local_ok)exit
      do j=1,size(fragment_bases(b)%global_ids)
        slot=find_sorted_int64(uncompressed_global_basis_ids,fragment_bases(b)%global_ids(j))
        if(slot==0)then
          local_ok=.false.;local_message='fragment publisher basis ID is absent from canonical union';exit
        endif
        local_id_count(slot)=local_id_count(slot)+1;local_owner_code(slot)=rank+1
        local_fragment(slot)=f;local_slot(slot)=j;local_sector(slot)=fragment_bases(b)%sector(j)
        local_generation(slot)=fragment_bases(b)%generation
      enddo
      if(.not.local_ok)exit
    enddo
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call MPI_Allreduce(local_publishers,global_publishers,nfragment,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment publisher count reduction failed';return;endif
    call MPI_Allreduce(local_id_count,global_id_count,nuncompressed,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='canonical basis ID count reduction failed';return;endif
    call MPI_Allreduce(local_owner_code,global_owner_code,nuncompressed,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='canonical basis owner reduction failed';return;endif
    call MPI_Allreduce(local_fragment,global_fragment,nuncompressed,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='canonical basis fragment reduction failed';return;endif
    call MPI_Allreduce(local_slot,global_slot,nuncompressed,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='canonical basis slot reduction failed';return;endif
    call MPI_Allreduce(local_sector,global_sector,nuncompressed,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='canonical basis sector reduction failed';return;endif
    call MPI_Allreduce(local_generation,global_generation,nuncompressed,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='canonical basis generation reduction failed';return;endif
    local_ok=.true.;local_message='dual catalog publisher tuple reduction failed'
    if(local_ok.and.any(global_publishers/=1))then
      local_ok=.false.;local_message='missing or duplicate fragment publisher'
    endif
    if(local_ok.and.any(global_id_count>1))then
      local_ok=.false.;local_message='duplicate publisher basis ID in canonical union'
    endif
    if(local_ok.and.any(global_id_count==0))then
      local_ok=.false.;local_message='missing publisher basis ID in canonical union'
    endif
    if(local_ok.and.minval(global_generation)/=maxval(global_generation))then
      local_ok=.false.;local_message='fragment publishers mix basis generations'
    endif
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif

    call validate_coverage_pairs(comm,fragment_bases,required_interface_fragment_ids,&
      required_interface_row_ids,'interface',stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call validate_coverage_pairs(comm,fragment_bases,required_periodic_wrap_fragment_ids,&
      required_periodic_wrap_row_ids,'wrap',stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call validate_coverage_pairs(comm,fragment_bases,required_projector_fragment_ids,&
      required_projector_row_ids,'projector',stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call validate_packet_neighbor_evidence(uncompressed_global_basis_ids,packet_ids,packet_neighbor_offsets,&
      packet_neighbor_basis_ids,required_neighbor_packet_ids,required_neighbor_basis_ids,stage_ok,stage_message)
    call synchronize_status(comm,stage_ok,stage_message,local_ok,local_message)
    if(.not.local_ok)then;message=trim(local_message);return;endif

    fragment_hash=int(z'6A09E667F3BCC909',int64);total_omitted=0d0
    call hash_accumulate(fragment_hash,int(nproc,int64))
    call hash_accumulate(fragment_hash,int(nfragment,int64))
    do f=1,nfragment
      ncolumn=expected_fragment_ranks(f);root=-1
      allocate(fragment_slots(ncolumn),stat=allocation_status)
      local_ok=allocation_status==0;local_message='cannot allocate fragment publisher slot map'
      call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      fragment_slots=0;local_ok=.true.;local_message=''
      do slot=1,nuncompressed
        if(global_fragment(slot)/=f)cycle
        j=global_slot(slot)
        if(j<1.or.j>ncolumn)then
          local_ok=.false.;local_message='fragment publisher local slot is out of range';exit
        endif
        if(fragment_slots(j)/=0)then
          local_ok=.false.;local_message='duplicate fragment publisher local slot';exit
        endif
        fragment_slots(j)=slot
        if(root<0)root=global_owner_code(slot)-1
        if(root/=global_owner_code(slot)-1)then
          local_ok=.false.;local_message='fragment basis columns have different publisher ranks';exit
        endif
      enddo
      if(local_ok.and.(root<0.or.any(fragment_slots==0)))then
        local_ok=.false.;local_message='fragment publisher slot map is incomplete'
      endif
      call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif

      npoint=0;provenance=0_int64;entry=0;publisher_metadata=[0,ncolumn]
      if(rank==root)then
        do b=1,size(fragment_bases)
          if(fragment_bases(b)%fragment_id==f)then;entry=b;exit;endif
        enddo
        if(entry>0)then
          npoint=size(fragment_bases(entry)%buffer_point_ids)
          provenance=fragment_bases(entry)%provenance_fingerprint
          publisher_metadata=[npoint,size(fragment_bases(entry)%global_ids)]
        endif
      endif
      call MPI_Bcast(publisher_metadata,size(publisher_metadata),MPI_INTEGER,root,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='fragment publisher metadata broadcast failed';return;endif
      call MPI_Bcast(provenance,1,MPI_INTEGER8,root,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='fragment publisher provenance broadcast failed';return;endif
      npoint=publisher_metadata(1)
      local_ok=npoint>=0.and.publisher_metadata(2)==ncolumn.and.provenance/=0_int64
      local_message='fragment publisher metadata is invalid'
      call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif

      publisher_value_count64=0_int64;arithmetic_ok=.true.
      call checked_multiply_nonnegative_int64(int(npoint,int64),int(ncolumn,int64),&
        publisher_value_count64,arithmetic_ok)
      local_ok=arithmetic_ok.and.publisher_value_count64<=int(huge(publisher_value_count),int64)
      if(local_ok)publisher_value_count=int(publisher_value_count64)
      local_message='fragment publisher MPI count overflow'
      call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif

      allocate(publisher_point_ids(npoint),publisher_values(npoint,ncolumn),&
        local_omitted_by_column(ncolumn),global_omitted_by_column(ncolumn),stat=allocation_status)
      local_ok=allocation_status==0;local_message='cannot allocate batched fragment publisher workspace'
      call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      if(rank==root)then
        publisher_point_ids=fragment_bases(entry)%buffer_point_ids
        publisher_values=fragment_bases(entry)%buffer_values
      endif
      call MPI_Bcast(publisher_point_ids,npoint,MPI_INTEGER8,root,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='fragment publisher row-ID broadcast failed';return;endif
      call MPI_Bcast(publisher_values,publisher_value_count,MPI_DOUBLE_COMPLEX,root,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='fragment publisher value broadcast failed';return;endif

      local_omitted_by_column=0d0;local_ok=.true.;local_message=''
      do j=1,ncolumn
        slot=fragment_slots(j)
        do p=1,nlocal
          entry=find_int64(publisher_point_ids,row_ids(p))
          if(entry==0)then
            amplitude=abs(uncompressed_basis_values(slot,p));weight_root=sqrt(weights(p))
            if(.not.ieee_is_finite(amplitude).or.amplitude>huge(1d0)/weight_root)then
              local_ok=.false.;local_message='fragment publisher omitted tail norm overflow';exit
            endif
            local_norm=amplitude*weight_root
            if(local_norm>sqrt(huge(1d0)))then
              local_ok=.false.;local_message='fragment publisher omitted tail norm overflow';exit
            endif
            omitted_term=local_norm*local_norm
            if(local_omitted_by_column(j)>huge(1d0)-omitted_term)then
              local_ok=.false.;local_message='fragment publisher omitted tail sum overflow';exit
            endif
            local_omitted_by_column(j)=local_omitted_by_column(j)+omitted_term
          elseif(.not.bitwise_complex_scalar_equal(publisher_values(entry,j),&
              uncompressed_basis_values(slot,p)))then
            local_ok=.false.;local_message='fragment publisher present value differs from canonical basis';exit
          endif
        enddo
        if(.not.local_ok)exit
      enddo
      call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      call MPI_Allreduce(local_omitted_by_column,global_omitted_by_column,ncolumn,&
        MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='fragment publisher tail reduction failed';return;endif
      if(.not.all(ieee_is_finite(global_omitted_by_column)).or.any(global_omitted_by_column<0d0).or.&
          any(sqrt(global_omitted_by_column)>tail_tolerance))then
        message='fragment publisher omitted tail exceeds tolerance';return
      endif
      do j=1,ncolumn
        if(total_omitted>huge(1d0)-global_omitted_by_column(j))then
          message='aggregate fragment publisher omitted tail overflow';return
        endif
        total_omitted=total_omitted+global_omitted_by_column(j)
      enddo

      call hash_accumulate(fragment_hash,int(f,int64))
      call hash_accumulate(fragment_hash,int(root,int64))
      call hash_accumulate(fragment_hash,int(ncolumn,int64))
      call hash_accumulate(fragment_hash,int(npoint,int64))
      call hash_accumulate(fragment_hash,provenance)
      do p=1,npoint;call hash_accumulate(fragment_hash,publisher_point_ids(p));enddo
      do j=1,ncolumn
        slot=fragment_slots(j)
        call hash_accumulate(fragment_hash,uncompressed_global_basis_ids(slot))
        call hash_accumulate(fragment_hash,int(global_slot(slot),int64))
        call hash_accumulate(fragment_hash,int(global_sector(slot),int64))
        call hash_accumulate(fragment_hash,int(global_generation(slot),int64))
        do p=1,npoint
          bits=transfer(real(publisher_values(p,j),real64),bits);call hash_accumulate(fragment_hash,bits)
          bits=transfer(aimag(publisher_values(p,j)),bits);call hash_accumulate(fragment_hash,bits)
        enddo
      enddo
      deallocate(fragment_slots,publisher_point_ids,publisher_values,local_omitted_by_column,&
        global_omitted_by_column)
    enddo
    if(.not.ieee_is_finite(total_omitted).or.sqrt(total_omitted)>tail_tolerance)then
      message='aggregate fragment publisher omitted tail exceeds tolerance';return
    endif
    do j=1,size(packet_ids);call hash_accumulate(fragment_hash,int(packet_ids(j),int64));enddo
    do j=1,size(packet_neighbor_offsets)
      call hash_accumulate(fragment_hash,int(packet_neighbor_offsets(j),int64))
    enddo
    do j=1,size(packet_neighbor_basis_ids);call hash_accumulate(fragment_hash,packet_neighbor_basis_ids(j));enddo
    if(fragment_hash==0_int64)fragment_hash=1_int64

    gram_count64=0_int64;arithmetic_ok=.true.
    call checked_multiply_nonnegative_int64(int(nuncompressed,int64),int(nuncompressed,int64),&
      gram_count64,arithmetic_ok)
    local_ok=arithmetic_ok.and.gram_count64<=int(huge(gram_mpi_count),int64)
    if(local_ok)gram_mpi_count=int(gram_count64)
    local_message='dual catalog Gram MPI count overflow'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    allocate(local_gram(nuncompressed,nuncompressed),gram(nuncompressed,nuncompressed),&
      stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate dual catalog Gram workspace'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    do j=1,nuncompressed;do b=1,nuncompressed
      local_gram(b,j)=sum(weights*conjg(uncompressed_basis_values(b,:))*uncompressed_basis_values(j,:))
    enddo;enddo
    call MPI_Allreduce(local_gram,gram,gram_mpi_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='dual catalog Gram reduction failed';return;endif
    call hermitian_pseudoinverse(gram,metric_tolerance,gram_inverse,gram_rank,gram_eigenvalues,&
      stage_ok,stage_message)
    call synchronize_status(comm,stage_ok,stage_message,local_ok,local_message)
    if(.not.local_ok)then;message=trim(local_message);return;endif
    if(gram_rank/=ncomplete)then;message='terminal map rank disagrees with physical union Gram rank';return;endif

    orthogonality_tolerance=max(100d0*metric_tolerance,&
      100d0*epsilon(1d0)*real(max(nuncompressed,ncomplete),real64))
    column_defect=0d0
    do j=1,ncomplete
      do b=1,ncomplete
        column_defect=max(column_defect,abs(sum(conjg(union_to_complete(:,b))*&
          union_to_complete(:,j))-merge(1d0,0d0,b==j)))
      enddo
    enddo
    if(column_defect>orthogonality_tolerance)then
      message='terminal map columns are not orthonormal';return
    endif
    allocate(retained_projector(nuncompressed,nuncompressed),&
      metric_projector(nuncompressed,nuncompressed),stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate terminal span projector workspace'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    retained_projector=matmul(union_to_complete,conjg(transpose(union_to_complete)))
    metric_projector=matmul(gram,gram_inverse)
    retained_span_defect=maxval(abs(retained_projector-metric_projector))
    if(.not.ieee_is_finite(retained_span_defect).or.retained_span_defect>orthogonality_tolerance)then
      message='terminal map retained span disagrees with the physical union Gram range';return
    endif
    identity_map=ncomplete==nuncompressed.and.bitwise_identity_matrix(union_to_complete)
    if(ncomplete==nuncompressed.and..not.identity_map)then
      message='full-rank terminal map must be bitwise identity';return
    endif
    allocate(complete_gram(ncomplete,ncomplete),stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate retained terminal Gram'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    complete_gram=matmul(conjg(transpose(union_to_complete)),matmul(gram,union_to_complete))
    call hermitian_pseudoinverse(complete_gram,metric_tolerance,complete_gram_inverse,complete_gram_rank,&
      complete_eigenvalues,stage_ok,stage_message)
    call synchronize_status(comm,stage_ok,stage_message,local_ok,local_message)
    if(.not.local_ok)then;message=trim(local_message);return;endif
    if(complete_gram_rank/=ncomplete)then;message='terminal map retained metric is rank deficient';return;endif
    allocate(seed_gram(nseed,nseed),stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate physical seed Gram'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    seed_gram=matmul(conjg(transpose(seed_coefficients_in_uncompressed)),&
      matmul(gram,seed_coefficients_in_uncompressed))
    call hermitian_pseudoinverse(seed_gram,metric_tolerance,seed_gram_inverse,seed_rank,seed_eigenvalues,&
      stage_ok,stage_message)
    call synchronize_status(comm,stage_ok,'physical seed metric failure: '//trim(stage_message),&
      local_ok,local_message)
    if(.not.local_ok)then;message=trim(local_message);return;endif
    if(seed_rank/=nseed)then;message='physical seed rank is deficient before terminal composition';return;endif
    allocate(working_complete_seed(ncomplete,nseed),stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate complete seed coefficient workspace'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    if(identity_map)then
      working_complete_seed=seed_coefficients_in_uncompressed
    else
      allocate(rhs(ncomplete,nseed),stat=allocation_status)
      local_ok=allocation_status==0;local_message='cannot allocate terminal seed right-hand side'
      call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      rhs=matmul(conjg(transpose(union_to_complete)),matmul(gram,seed_coefficients_in_uncompressed))
      working_complete_seed=matmul(complete_gram_inverse,rhs)
    endif
    allocate(residual(nuncompressed),stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate seed reconstruction workspace'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    working_seed_defect=0d0
    do j=1,nseed
      residual=seed_coefficients_in_uncompressed(:,j)-matmul(union_to_complete,working_complete_seed(:,j))
      seed_norm=max(0d0,real(dot_product(seed_coefficients_in_uncompressed(:,j),&
        matmul(gram,seed_coefficients_in_uncompressed(:,j))),real64))
      residual_norm=max(0d0,real(dot_product(residual,matmul(gram,residual)),real64))
      if(seed_norm<=0d0.or..not.ieee_is_finite(seed_norm).or..not.ieee_is_finite(residual_norm))then
        message='physical seed reconstruction norm is invalid';return
      endif
      working_seed_defect=max(working_seed_defect,sqrt(residual_norm/seed_norm))
    enddo
    if(working_seed_defect>orthogonality_tolerance)then
      message='terminal compression loses a normalized physical seed';return
    endif

    map_hash=int(z'BB67AE8584CAA73B',int64)
    do j=1,nuncompressed;call hash_accumulate(map_hash,uncompressed_global_basis_ids(j));enddo
    call hash_accumulate(map_hash,int(ncomplete,int64))
    do j=1,nuncompressed;do b=1,nuncompressed
      bits=transfer(real(retained_projector(b,j),real64),bits)
      call hash_accumulate(map_hash,bits)
      bits=transfer(aimag(retained_projector(b,j)),bits)
      call hash_accumulate(map_hash,bits)
    enddo;enddo
    ! The raw terminal map deliberately remains a separate, gauge-sensitive
    ! part of the receipt in addition to the physical retained projector.
    do j=1,ncomplete;do b=1,nuncompressed
      bits=transfer(real(union_to_complete(b,j),real64),bits)
      call hash_accumulate(map_hash,bits)
      bits=transfer(aimag(union_to_complete(b,j)),bits)
      call hash_accumulate(map_hash,bits)
    enddo;enddo
    allocate(local_row_hash(global_row_count),global_row_hash(global_row_count),stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate terminal map row fingerprint workspace'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    local_row_hash=0_int64
    do p=1,nlocal
      local_row_hash(int(row_ids(p)))=raw_basis_row_hash(row_ids(p),weights(p),&
        uncompressed_basis_values(:,p))
    enddo
    call MPI_Allreduce(local_row_hash,global_row_hash,global_row_count,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='terminal map row fingerprint reduction failed';return;endif
    do p=1,global_row_count
      call hash_accumulate(map_hash,global_row_hash(p))
    enddo
    if(map_hash==0_int64)map_hash=1_int64
    call compute_dg_hybrid_union_to_complete_binding(comm,union_to_complete,map_hash,&
      metric_tolerance,transform_binding_hash,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif

    call copy_fragment_basis_array_collective(comm,fragment_bases,working_catalog%fragment_bases,&
      stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    allocate(working_catalog%union_to_complete,source=union_to_complete,stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate immutable terminal map'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    allocate(working_catalog%uncompressed_global_basis_ids,source=uncompressed_global_basis_ids,&
      stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate canonical dual basis IDs'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    allocate(working_catalog%uncompressed_owner_ranks(nuncompressed),&
      working_catalog%uncompressed_fragment_ids(nuncompressed),&
      working_catalog%uncompressed_local_slots(nuncompressed),&
      working_catalog%uncompressed_sectors(nuncompressed),&
      working_catalog%uncompressed_generations(nuncompressed),stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate canonical ownership tuples'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    working_catalog%uncompressed_owner_ranks=global_owner_code-1
    working_catalog%uncompressed_fragment_ids=global_fragment
    working_catalog%uncompressed_local_slots=global_slot
    working_catalog%uncompressed_sectors=global_sector
    working_catalog%uncompressed_generations=global_generation
    working_catalog%uncompressed_rank=nuncompressed;working_catalog%complete_rank=ncomplete
    working_catalog%fragment_catalog_fingerprint=fragment_hash
    working_catalog%complete_map_fingerprint=map_hash
    working_catalog%complete_transform_binding_fingerprint=transform_binding_hash
    working_catalog%valid=.true.
    allocate(working_seed_owner,source=seed_fragment_owner,stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate preserved seed owners'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    allocate(working_preserved_seed,source=seed_coefficients_in_uncompressed,stat=allocation_status)
    local_ok=allocation_status==0;local_message='cannot allocate preserved seed coordinates'
    call synchronize_status(comm,local_ok,local_message,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif

    call move_dual_catalog(working_catalog,catalog)
    call move_alloc(working_seed_owner,preserved_seed_fragment_owner)
    call move_alloc(working_preserved_seed,preserved_seed_coefficients_in_uncompressed)
    call move_alloc(working_complete_seed,complete_seed_coefficients)
    seed_reconstruction_defect=working_seed_defect;ok=.true.;message=''
  end subroutine finalize_dg_hybrid_dual_basis_catalog

  subroutine validate_replicated_dual_inputs(comm,expected_fragment_ranks,basis_ids,terminal_map,&
      seed_owner,seed_coefficients,metric_tolerance,interface_fragments,interface_rows,wrap_fragments,&
      wrap_rows,projector_fragments,projector_rows,packet_ids,packet_offsets,neighbor_ids,&
      required_packets,required_neighbors,tail_tolerance,ok,message)
    integer,intent(in)::comm,expected_fragment_ranks(:),seed_owner(:),interface_fragments(:),&
      wrap_fragments(:),projector_fragments(:),packet_ids(:),packet_offsets(:),required_packets(:)
    integer(int64),intent(in)::basis_ids(:),interface_rows(:),wrap_rows(:),projector_rows(:),&
      neighbor_ids(:),required_neighbors(:)
    complex(real64),intent(in)::terminal_map(:,:),seed_coefficients(:,:)
    real(real64),intent(in)::metric_tolerance,tail_tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64)::scalar_bits(2)
    ok=.false.;message=''
    if(.not.replicated_integer_array(comm,expected_fragment_ranks))then
      message='expected fragment ranks differ between ranks';return
    endif
    if(.not.replicated_int64_array(comm,basis_ids))then
      message='canonical Hybrid basis IDs differ between ranks';return
    endif
    if(.not.replicated_complex_matrix(comm,terminal_map))then
      message='terminal map differs between ranks';return
    endif
    if(.not.replicated_integer_array(comm,seed_owner))then
      message='seed fragment owners differ between ranks';return
    endif
    if(.not.replicated_complex_matrix(comm,seed_coefficients))then
      message='seed coordinates differ between ranks';return
    endif
    scalar_bits=[transfer(metric_tolerance,scalar_bits(1)),transfer(tail_tolerance,scalar_bits(2))]
    if(.not.replicated_int64_array(comm,scalar_bits))then
      message='dual catalog tolerances differ between ranks';return
    endif
    if(.not.replicated_integer_array(comm,interface_fragments))then
      message='interface coverage evidence differs between ranks';return
    endif
    if(.not.replicated_int64_array(comm,interface_rows))then
      message='interface coverage evidence differs between ranks';return
    endif
    if(.not.replicated_integer_array(comm,wrap_fragments))then
      message='periodic wrap coverage evidence differs between ranks';return
    endif
    if(.not.replicated_int64_array(comm,wrap_rows))then
      message='periodic wrap coverage evidence differs between ranks';return
    endif
    if(.not.replicated_integer_array(comm,projector_fragments))then
      message='projector coverage evidence differs between ranks';return
    endif
    if(.not.replicated_int64_array(comm,projector_rows))then
      message='projector coverage evidence differs between ranks';return
    endif
    if(.not.replicated_integer_array(comm,packet_ids))then
      message='packet neighbor provenance differs between ranks';return
    endif
    if(.not.replicated_integer_array(comm,packet_offsets))then
      message='packet neighbor provenance differs between ranks';return
    endif
    if(.not.replicated_int64_array(comm,neighbor_ids))then
      message='packet neighbor provenance differs between ranks';return
    endif
    if(.not.replicated_integer_array(comm,required_packets))then
      message='packet neighbor provenance differs between ranks';return
    endif
    if(.not.replicated_int64_array(comm,required_neighbors))then
      message='packet neighbor provenance differs between ranks';return
    endif
    ok=.true.
  end subroutine validate_replicated_dual_inputs

  logical function replicated_integer_array(comm,values)result(equal)
    integer,intent(in)::comm,values(:)
    integer::rank,ierr,allocation_status,local_bad,global_bad,count
    integer(int64)::count64
    integer,allocatable::reference(:)
    equal=.false.;count64=size(values,kind=int64)
    local_bad=merge(0,1,count64<=int(huge(count),int64))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    count=int(count64)
    allocate(reference(count),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    if(rank==0)reference=values
    call MPI_Bcast(reference,count,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    local_bad=merge(0,1,all(values==reference))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    equal=ierr==MPI_SUCCESS.and.global_bad==0
  end function replicated_integer_array

  logical function replicated_int64_array(comm,values)result(equal)
    integer,intent(in)::comm
    integer(int64),intent(in)::values(:)
    integer::rank,ierr,allocation_status,local_bad,global_bad,count
    integer(int64)::count64
    integer(int64),allocatable::reference(:)
    equal=.false.;count64=size(values,kind=int64)
    local_bad=merge(0,1,count64<=int(huge(count),int64))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    count=int(count64)
    allocate(reference(count),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    if(rank==0)reference=values
    call MPI_Bcast(reference,count,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    local_bad=merge(0,1,all(values==reference))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    equal=ierr==MPI_SUCCESS.and.global_bad==0
  end function replicated_int64_array

  logical function replicated_complex_matrix(comm,values)result(equal)
    integer,intent(in)::comm
    complex(real64),intent(in)::values(:,:)
    integer::rank,i,j,ierr,allocation_status,local_bad,global_bad,count
    integer(int64)::count64
    complex(real64),allocatable::reference(:,:)
    equal=.false.;count64=size(values,kind=int64)
    local_bad=merge(0,1,count64<=int(huge(count),int64))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    count=int(count64)
    allocate(reference(size(values,1),size(values,2)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    if(rank==0)reference=values
    call MPI_Bcast(reference,count,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    local_bad=0
    do j=1,size(values,2);do i=1,size(values,1)
      if(.not.bitwise_complex_scalar_equal(values(i,j),reference(i,j)))local_bad=1
    enddo;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    equal=ierr==MPI_SUCCESS.and.global_bad==0
  end function replicated_complex_matrix

  subroutine validate_coverage_pairs(comm,bases,fragment_ids,row_ids,label,ok,message)
    integer,intent(in)::comm,fragment_ids(:)
    integer(int64),intent(in)::row_ids(:)
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    character(*),intent(in)::label
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,b,p,ierr,allocation_status
    integer,allocatable::local_hit(:),global_hit(:)
    logical::local_ok,stage_ok
    character(256)::stage_message
    ok=.false.;message=''
    allocate(local_hit(size(fragment_ids)),global_hit(size(fragment_ids)),stat=allocation_status)
    local_ok=allocation_status==0
    call synchronize_status(comm,local_ok,'cannot allocate coverage evidence workspace',stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    local_hit=0
    do i=1,size(fragment_ids)
      do b=1,size(bases)
        if(bases(b)%fragment_id/=fragment_ids(i))cycle
        if(.not.allocated(bases(b)%buffer_point_ids))cycle
        p=find_int64(bases(b)%buffer_point_ids,row_ids(i))
        if(p>0)local_hit(i)=1
      enddo
    enddo
    call MPI_Allreduce(local_hit,global_hit,size(fragment_ids),MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(global_hit/=1))then
      message=trim(label)//' coverage is incomplete';return
    endif
    ok=.true.
  end subroutine validate_coverage_pairs

  subroutine validate_packet_neighbor_evidence(basis_ids,packet_ids,offsets,neighbor_ids,&
      required_packets,required_neighbors,ok,message)
    integer(int64),intent(in)::basis_ids(:),neighbor_ids(:),required_neighbors(:)
    integer,intent(in)::packet_ids(:),offsets(:),required_packets(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,p
    logical::found
    ok=.false.;message=''
    if(size(packet_ids)<1.or.size(offsets)/=size(packet_ids)+1.or.offsets(1)/=1.or.&
        offsets(size(offsets))/=size(neighbor_ids)+1.or.any(packet_ids<1).or.&
        any(offsets(2:)<offsets(:size(offsets)-1)))then
      message='invalid packet neighbor CSR provenance';return
    endif
    do i=1,size(packet_ids)
      if(any(packet_ids(i+1:)==packet_ids(i)))then;message='duplicate packet provenance ID';return;endif
      do j=offsets(i),offsets(i+1)-1
        if(find_sorted_int64(basis_ids,neighbor_ids(j))==0)then
          message='packet neighbor basis ID is absent from canonical union';return
        endif
        if(j>offsets(i))then
          if(any(neighbor_ids(offsets(i):j-1)==neighbor_ids(j)))then
            message='duplicate packet neighbor basis ID';return
          endif
        endif
      enddo
    enddo
    do i=1,size(required_packets)
      p=0
      do j=1,size(packet_ids)
        if(packet_ids(j)==required_packets(i))then;p=j;exit;endif
      enddo
      found=.false.
      if(p>0)found=any(neighbor_ids(offsets(p):offsets(p+1)-1)==required_neighbors(i))
      if(.not.found)then;message='required packet neighbor provenance is incomplete';return;endif
    enddo
    ok=.true.
  end subroutine validate_packet_neighbor_evidence

  subroutine hermitian_pseudoinverse(matrix,tolerance,inverse,retained_rank,eigenvalues,ok,message)
    complex(real64),intent(in)::matrix(:,:)
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(out)::inverse(:,:)
    integer,intent(out)::retained_rank
    real(real64),allocatable,intent(out)::eigenvalues(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::vectors(:,:),work(:)
    real(real64),allocatable::rwork(:)
    real(real64)::scale,cutoff,negative_limit,roundoff_floor,hermitian_defect
    integer::n,i,j,k,info,allocation_status
    ok=.false.;message='';retained_rank=0;n=size(matrix,1)
    if(n<1.or.size(matrix,2)/=n.or..not.finite_complex_matrix(matrix).or.&
        .not.ieee_is_finite(tolerance).or.tolerance<=0d0)then
      message='invalid Hermitian metric matrix';return
    endif
    scale=max(1d0,maxval(abs(matrix)));hermitian_defect=maxval(abs(matrix-conjg(transpose(matrix))))
    if(hermitian_defect>100d0*epsilon(1d0)*real(n,real64)*scale)then
      message='metric matrix is not Hermitian';return
    endif
    allocate(vectors(n,n),eigenvalues(n),work(max(1,2*n-1)),rwork(max(1,3*n-2)),inverse(n,n),&
      stat=allocation_status)
    if(allocation_status/=0)then;message='cannot allocate Hermitian eigensolver workspace';return;endif
    vectors=0.5d0*(matrix+conjg(transpose(matrix)))
    call zheev('V','U',n,vectors,n,eigenvalues,work,size(work),rwork,info)
    if(info/=0.or..not.all(ieee_is_finite(eigenvalues)))then
      message='Hermitian metric eigensolver failed';return
    endif
    scale=max(1d0,maxval(abs(eigenvalues)))
    roundoff_floor=64d0*epsilon(1d0)*scale*real(max(1,n),real64)
    cutoff=max(tolerance*scale,roundoff_floor);negative_limit=roundoff_floor
    if(minval(eigenvalues)<-negative_limit)then;message='indefinite Hermitian metric';return;endif
    ! Match the generalized projection policy: roundoff-sized modes are null.
    if(any(eigenvalues>roundoff_floor.and.abs(eigenvalues-cutoff)<=16d0*roundoff_floor))then
      message='ambiguous Hermitian metric rank at cutoff';return
    endif
    retained_rank=count(eigenvalues>cutoff);inverse=(0d0,0d0)
    do k=1,n
      if(eigenvalues(k)<=cutoff)cycle
      do j=1,n;do i=1,n
        inverse(i,j)=inverse(i,j)+vectors(i,k)*conjg(vectors(j,k))/eigenvalues(k)
      enddo;enddo
    enddo
    ok=.true.
  end subroutine hermitian_pseudoinverse

  subroutine synchronize_status(comm,input_ok,input_message,output_ok,output_message)
    integer,intent(in)::comm
    logical,intent(in)::input_ok
    character(*),intent(in)::input_message
    logical,intent(out)::output_ok
    character(*),intent(out)::output_message
    integer::rank,ierr,candidate,failure_rank
    character(512)::buffer
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;output_ok=.false.;output_message='MPI communicator status failure';return;endif
    candidate=merge(huge(0),rank,input_ok)
    call MPI_Allreduce(candidate,failure_rank,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;output_ok=.false.;output_message='MPI collective status failure';return;endif
    if(failure_rank==huge(0))then;output_ok=.true.;output_message='';return;endif
    buffer='';if(rank==failure_rank)buffer=trim(input_message)
    call MPI_Bcast(buffer,len(buffer),MPI_CHARACTER,failure_rank,comm,ierr)
    output_ok=.false.
    if(ierr/=MPI_SUCCESS)then
      output_message='MPI collective diagnostic failure'
    elseif(len_trim(buffer)==0)then
      output_message='collective pipeline failure'
    else
      output_message=trim(buffer)
    endif
  end subroutine synchronize_status

  subroutine clear_fragment_basis(basis)
    type(s_dg_hybrid_fragment_basis),intent(inout)::basis
    if(allocated(basis%global_ids))deallocate(basis%global_ids)
    if(allocated(basis%buffer_point_ids))deallocate(basis%buffer_point_ids)
    if(allocated(basis%sector))deallocate(basis%sector)
    if(allocated(basis%buffer_values))deallocate(basis%buffer_values)
    basis%fragment_id=0;basis%generation=0;basis%provenance_fingerprint=0_int64
  end subroutine clear_fragment_basis

  subroutine move_fragment_basis(source,destination)
    type(s_dg_hybrid_fragment_basis),intent(inout)::source,destination
    call clear_fragment_basis(destination)
    destination%fragment_id=source%fragment_id;destination%generation=source%generation
    destination%provenance_fingerprint=source%provenance_fingerprint
    call move_alloc(source%global_ids,destination%global_ids)
    call move_alloc(source%buffer_point_ids,destination%buffer_point_ids)
    call move_alloc(source%sector,destination%sector)
    call move_alloc(source%buffer_values,destination%buffer_values)
    source%fragment_id=0;source%generation=0;source%provenance_fingerprint=0_int64
  end subroutine move_fragment_basis

  subroutine copy_fragment_basis_array_collective(comm,source,destination,ok,message)
    integer,intent(in)::comm
    type(s_dg_hybrid_fragment_basis),intent(in)::source(:)
    type(s_dg_hybrid_fragment_basis),allocatable,intent(out)::destination(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::b,allocation_status
    logical::local_ok,stage_ok
    character(256)::stage_message
    ok=.false.;message=''
    allocate(destination(size(source)),stat=allocation_status)
    local_ok=allocation_status==0
    call synchronize_status(comm,local_ok,'cannot allocate final fragment catalog',stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    do b=1,size(source)
      destination(b)%fragment_id=source(b)%fragment_id;destination(b)%generation=source(b)%generation
      destination(b)%provenance_fingerprint=source(b)%provenance_fingerprint
      allocation_status=0
      if(allocated(source(b)%global_ids))allocate(destination(b)%global_ids(size(source(b)%global_ids)),&
        stat=allocation_status)
      local_ok=allocation_status==0
      call synchronize_status(comm,local_ok,'cannot copy fragment global IDs',stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      if(allocated(source(b)%global_ids))destination(b)%global_ids=source(b)%global_ids
      allocation_status=0
      if(allocated(source(b)%buffer_point_ids))allocate(&
        destination(b)%buffer_point_ids(size(source(b)%buffer_point_ids)),stat=allocation_status)
      local_ok=allocation_status==0
      call synchronize_status(comm,local_ok,'cannot copy fragment buffer point IDs',stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      if(allocated(source(b)%buffer_point_ids))destination(b)%buffer_point_ids=source(b)%buffer_point_ids
      allocation_status=0
      if(allocated(source(b)%sector))allocate(destination(b)%sector(size(source(b)%sector)),&
        stat=allocation_status)
      local_ok=allocation_status==0
      call synchronize_status(comm,local_ok,'cannot copy fragment sectors',stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      if(allocated(source(b)%sector))destination(b)%sector=source(b)%sector
      allocation_status=0
      if(allocated(source(b)%buffer_values))allocate(destination(b)%buffer_values(&
        size(source(b)%buffer_values,1),size(source(b)%buffer_values,2)),stat=allocation_status)
      local_ok=allocation_status==0
      call synchronize_status(comm,local_ok,'cannot copy fragment buffer values',stage_ok,stage_message)
      if(.not.stage_ok)then;message=trim(stage_message);return;endif
      if(allocated(source(b)%buffer_values))destination(b)%buffer_values=source(b)%buffer_values
    enddo
    ok=.true.
  end subroutine copy_fragment_basis_array_collective

  subroutine clear_dual_catalog(catalog)
    type(s_dg_hybrid_dual_basis_catalog),intent(inout)::catalog
    if(allocated(catalog%fragment_bases))deallocate(catalog%fragment_bases)
    if(allocated(catalog%union_to_complete))deallocate(catalog%union_to_complete)
    if(allocated(catalog%uncompressed_global_basis_ids))deallocate(catalog%uncompressed_global_basis_ids)
    if(allocated(catalog%uncompressed_owner_ranks))deallocate(catalog%uncompressed_owner_ranks)
    if(allocated(catalog%uncompressed_fragment_ids))deallocate(catalog%uncompressed_fragment_ids)
    if(allocated(catalog%uncompressed_local_slots))deallocate(catalog%uncompressed_local_slots)
    if(allocated(catalog%uncompressed_sectors))deallocate(catalog%uncompressed_sectors)
    if(allocated(catalog%uncompressed_generations))deallocate(catalog%uncompressed_generations)
    catalog%valid=.false.;catalog%uncompressed_rank=0;catalog%complete_rank=0
    catalog%fragment_catalog_fingerprint=0_int64;catalog%complete_map_fingerprint=0_int64
    catalog%complete_transform_binding_fingerprint=0_int64
  end subroutine clear_dual_catalog

  subroutine move_dual_catalog(source,destination)
    type(s_dg_hybrid_dual_basis_catalog),intent(inout)::source
    type(s_dg_hybrid_dual_basis_catalog),intent(inout)::destination
    call clear_dual_catalog(destination)
    destination%valid=source%valid;destination%uncompressed_rank=source%uncompressed_rank
    destination%complete_rank=source%complete_rank
    destination%fragment_catalog_fingerprint=source%fragment_catalog_fingerprint
    destination%complete_map_fingerprint=source%complete_map_fingerprint
    destination%complete_transform_binding_fingerprint=source%complete_transform_binding_fingerprint
    call move_alloc(source%fragment_bases,destination%fragment_bases)
    call move_alloc(source%union_to_complete,destination%union_to_complete)
    call move_alloc(source%uncompressed_global_basis_ids,destination%uncompressed_global_basis_ids)
    call move_alloc(source%uncompressed_owner_ranks,destination%uncompressed_owner_ranks)
    call move_alloc(source%uncompressed_fragment_ids,destination%uncompressed_fragment_ids)
    call move_alloc(source%uncompressed_local_slots,destination%uncompressed_local_slots)
    call move_alloc(source%uncompressed_sectors,destination%uncompressed_sectors)
    call move_alloc(source%uncompressed_generations,destination%uncompressed_generations)
    source%valid=.false.;source%uncompressed_rank=0;source%complete_rank=0
    source%fragment_catalog_fingerprint=0_int64;source%complete_map_fingerprint=0_int64
    source%complete_transform_binding_fingerprint=0_int64
  end subroutine move_dual_catalog

  logical function allocated_nonempty_fragment_payload(basis)result(nonempty)
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    nonempty=.false.
    if(allocated(basis%global_ids))nonempty=nonempty.or.size(basis%global_ids)>0
    if(allocated(basis%buffer_point_ids))nonempty=nonempty.or.size(basis%buffer_point_ids)>0
    if(allocated(basis%sector))nonempty=nonempty.or.size(basis%sector)>0
    if(allocated(basis%buffer_values))nonempty=nonempty.or.size(basis%buffer_values)>0
  end function allocated_nonempty_fragment_payload

  logical function finite_complex_matrix(values)result(finite)
    complex(real64),intent(in)::values(:,:)
    finite=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex_matrix

  integer function find_sorted_int64(values,target)result(position)
    integer(int64),intent(in)::values(:),target
    integer::low,high,middle
    position=0;low=1;high=size(values)
    do while(low<=high)
      middle=low+(high-low)/2
      if(values(middle)==target)then;position=middle;return
      elseif(values(middle)<target)then;low=middle+1
      else;high=middle-1
      endif
    enddo
  end function find_sorted_int64

  integer function find_int64(values,target)result(position)
    integer(int64),intent(in)::values(:),target
    integer::i
    position=0
    do i=1,size(values)
      if(values(i)==target)then;position=i;return;endif
    enddo
  end function find_int64

  logical function bitwise_complex_scalar_equal(left,right)result(equal)
    complex(real64),intent(in)::left,right
    integer(int64)::left_bits,right_bits
    left_bits=transfer(real(left,real64),left_bits);right_bits=transfer(real(right,real64),right_bits)
    equal=left_bits==right_bits;if(.not.equal)return
    left_bits=transfer(aimag(left),left_bits);right_bits=transfer(aimag(right),right_bits)
    equal=left_bits==right_bits
  end function bitwise_complex_scalar_equal

  logical function bitwise_identity_matrix(matrix)result(identity)
    complex(real64),intent(in)::matrix(:,:)
    complex(real64)::expected
    integer::i,j
    identity=size(matrix,1)==size(matrix,2);if(.not.identity)return
    do j=1,size(matrix,2);do i=1,size(matrix,1)
      expected=(0d0,0d0);if(i==j)expected=(1d0,0d0)
      if(.not.bitwise_complex_scalar_equal(matrix(i,j),expected))then;identity=.false.;return;endif
    enddo;enddo
  end function bitwise_identity_matrix

  subroutine hash_accumulate(hash,word)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::word
    hash=ieor(ishftc(hash,13),word)
    hash=ieor(hash,ishftc(word,29))
    hash=ieor(ishftc(hash,7),int(z'9E3779B97F4A7C15',int64))
  end subroutine hash_accumulate

  integer(int64) function raw_basis_row_hash(row_id,weight,values)result(hash)
    integer(int64),intent(in)::row_id
    real(real64),intent(in)::weight
    complex(real64),intent(in)::values(:)
    integer::j
    integer(int64)::bits
    hash=int(z'3C6EF372FE94F82B',int64)
    call hash_accumulate(hash,row_id)
    bits=transfer(weight,bits);call hash_accumulate(hash,bits)
    do j=1,size(values)
      bits=transfer(real(values(j),real64),bits);call hash_accumulate(hash,bits)
      bits=transfer(aimag(values(j)),bits);call hash_accumulate(hash,bits)
    enddo
  end function raw_basis_row_hash

  subroutine checked_add_nonnegative_int64(accumulator,addend,valid)
    integer(int64),intent(inout)::accumulator
    integer(int64),intent(in)::addend
    logical,intent(inout)::valid
    if(.not.valid)return
    if(accumulator<0_int64.or.addend<0_int64.or.accumulator>huge(accumulator)-addend)then
      valid=.false.;return
    endif
    accumulator=accumulator+addend
  end subroutine checked_add_nonnegative_int64

  subroutine checked_multiply_nonnegative_int64(left,right,product,valid)
    integer(int64),intent(in)::left,right
    integer(int64),intent(out)::product
    logical,intent(inout)::valid
    product=0_int64
    if(.not.valid)return
    if(left<0_int64.or.right<0_int64)then
      valid=.false.;return
    endif
    if(left/=0_int64)then
      if(right>huge(product)/left)then;valid=.false.;return;endif
    endif
    product=left*right
  end subroutine checked_multiply_nonnegative_int64
end module dg_hybrid_projected_fragment_pipeline
