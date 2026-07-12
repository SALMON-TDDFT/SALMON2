module lcfo_wannier_sawf_templates
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter :: sawf_template_schema_version=1
  character(24),parameter :: sawf_template_magic='SALMON_SAWF_TEMPLATE_V1'

  type, public :: t_sawf_template_fingerprint
    character(256) :: geometry='', pseudopotential='', grid='', band_window=''
    character(256) :: complete_projection_shell='', symmetry='', buffer='', generator=''
    integer :: schema_version=sawf_template_schema_version
  end type

  type, public :: t_sawf_template_checkpoint
    type(t_sawf_template_fingerprint) :: fingerprint
    real(8), allocatable :: centers(:,:),spreads(:)
    real(8), allocatable :: basis(:,:),buffer_basis(:,:)
    complex(8), allocatable :: orbitals(:,:),d_band(:,:,:),d_wann(:,:,:)
    complex(8), allocatable :: gauge_unitary(:,:)
    real(8) :: gauge_residual=huge(0d0)
  end type

  public :: build_sawf_environment_orbits, validate_sawf_template_fingerprint
  public :: validate_sawf_actual_group_operation, replicate_sawf_operator_block
  public :: stitch_sawf_neighbor_gauge, validate_sawf_buffer_convergence
  public :: validate_sawf_global_local_equivalence
  public :: write_sawf_template_checkpoint, read_sawf_template_checkpoint
  public :: materialize_sawf_local_bases
  public :: apply_sawf_gauge_connection
  public :: stitch_and_apply_sawf_neighbor_pair
  public :: stitch_and_apply_sawf_neighbor_pair_real
  public :: build_sawf_neighbor_trace_overlap
  public :: whiten_sawf_buffered_overlap
  public :: stitch_sawf_buffered_neighbor_gauge

contains
  subroutine stitch_sawf_buffered_neighbor_gauge(cross_overlap,left_gram,right_gram, &
      relative_cutoff,alignment_tolerance,gauge_unitary,residual,ok,message)
    complex(8),intent(in)::cross_overlap(:,:),left_gram(:,:),right_gram(:,:)
    real(8),intent(in)::relative_cutoff,alignment_tolerance
    complex(8),intent(out)::gauge_unitary(size(cross_overlap,2),size(cross_overlap,1))
    real(8),intent(out)::residual;logical,intent(out)::ok;character(*),intent(out)::message
    complex(8),allocatable::whitened(:,:)
    allocate(whitened(size(cross_overlap,1),size(cross_overlap,2)))
    call whiten_sawf_buffered_overlap(cross_overlap,left_gram,right_gram,relative_cutoff, &
      whitened,ok,message)
    if(.not.ok)return
    call stitch_sawf_neighbor_gauge(whitened,alignment_tolerance,gauge_unitary,residual,ok,message)
  end subroutine

  subroutine whiten_sawf_buffered_overlap(cross_overlap,left_gram,right_gram,relative_cutoff, &
      whitened_overlap,ok,message)
    complex(8),intent(in)::cross_overlap(:,:),left_gram(:,:),right_gram(:,:)
    real(8),intent(in)::relative_cutoff
    complex(8),intent(out)::whitened_overlap(size(cross_overlap,1),size(cross_overlap,2))
    logical,intent(out)::ok;character(*),intent(out)::message
    complex(8),allocatable::left_inverse_sqrt(:,:),right_inverse_sqrt(:,:)
    logical::left_ok,right_ok;character(256)::detail
    ok=.false.;message='';whitened_overlap=(0d0,0d0)
    if(size(cross_overlap,1)/=size(cross_overlap,2).or. &
       any(shape(left_gram)/=shape(cross_overlap)).or.any(shape(right_gram)/=shape(cross_overlap)))then
      message='SAWF buffered overlap Gram dimensions are inconsistent';return
    end if
    call sawf_hermitian_inverse_sqrt(left_gram,relative_cutoff,left_inverse_sqrt,left_ok,detail)
    if(.not.left_ok)then;message='left '//trim(detail);return;end if
    call sawf_hermitian_inverse_sqrt(right_gram,relative_cutoff,right_inverse_sqrt,right_ok,detail)
    if(.not.right_ok)then;message='right '//trim(detail);return;end if
    whitened_overlap=matmul(left_inverse_sqrt,matmul(cross_overlap,right_inverse_sqrt))
    ok=.true.
  end subroutine

  subroutine sawf_hermitian_inverse_sqrt(gram,relative_cutoff,inverse_sqrt,ok,message)
    complex(8),intent(in)::gram(:,:);real(8),intent(in)::relative_cutoff
    complex(8),allocatable,intent(out)::inverse_sqrt(:,:)
    logical,intent(out)::ok;character(*),intent(out)::message
    complex(8),allocatable::eigenvectors(:,:),work(:),scaled(:,:)
    complex(8)::query(1);real(8),allocatable::eigenvalues(:),rwork(:)
    integer::n,info,lwork,i;external::zheev
    ok=.false.;message=''
    if(size(gram,1)<=0.or.size(gram,1)/=size(gram,2).or.relative_cutoff<=0d0)then
      message='SAWF buffered Gram dimensions or cutoff are invalid';return
    end if
    n=size(gram,1);allocate(eigenvectors(n,n),eigenvalues(n),rwork(max(1,3*n-2)))
    eigenvectors=gram;lwork=-1
    call zheev('V','U',n,eigenvectors,n,eigenvalues,query,lwork,rwork,info)
    if(info/=0)then;message='SAWF buffered Gram eigensolver query failed';return;end if
    lwork=max(1,int(real(query(1),8)));allocate(work(lwork));eigenvectors=gram
    call zheev('V','U',n,eigenvectors,n,eigenvalues,work,lwork,rwork,info)
    if(info/=0.or.maxval(eigenvalues)<=0d0)then
      message='SAWF buffered Gram eigensolver failed';return
    end if
    if(minval(eigenvalues)<=relative_cutoff*maxval(eigenvalues))then
      message='SAWF buffered Gram is rank deficient';return
    end if
    allocate(scaled(n,n),inverse_sqrt(n,n));scaled=eigenvectors
    do i=1,n;scaled(:,i)=scaled(:,i)/sqrt(eigenvalues(i));end do
    inverse_sqrt=matmul(scaled,conjg(transpose(eigenvectors)))
    ok=.true.
  end subroutine

  subroutine build_sawf_neighbor_trace_overlap(left_trace,right_trace,area_weight,overlap,ok,message)
    real(8),intent(in)::left_trace(:,:),right_trace(:,:),area_weight
    real(8),intent(out)::overlap(size(left_trace,2),size(right_trace,2))
    logical,intent(out)::ok; character(*),intent(out)::message
    ok=.false.;message='';overlap=0d0
    if(size(left_trace,1)<=0.or.size(left_trace,1)/=size(right_trace,1).or. &
       size(left_trace,2)/=size(right_trace,2).or.area_weight<=0d0)then
      message='SAWF neighbor face traces have inconsistent dimensions or weight';return
    end if
    if(.not.all(ieee_is_finite(left_trace)).or..not.all(ieee_is_finite(right_trace)).or. &
       .not.ieee_is_finite(area_weight))then
      message='SAWF neighbor face trace overlap input is non-finite';return
    end if
    overlap=area_weight*matmul(transpose(left_trace),right_trace)
    ok=.true.
  end subroutine

  subroutine stitch_and_apply_sawf_neighbor_pair_real(overlap,tolerance,neighbor_gauge, &
      basis,ww_block,wp_block,face_self_block,face_neighbor_block,gauge,residual,ok,message)
    real(8),intent(in)::overlap(:,:),neighbor_gauge(:,:),tolerance
    real(8),intent(inout)::basis(:,:),ww_block(:,:),wp_block(:,:),face_self_block(:,:), &
      face_neighbor_block(:,:)
    real(8),intent(out)::gauge(size(overlap,2),size(overlap,1)),residual
    logical,intent(out)::ok; character(*),intent(out)::message
    complex(8),allocatable::overlap_c(:,:),neighbor_c(:,:),basis_c(:,:),ww_c(:,:),wp_c(:,:), &
      self_c(:,:),face_c(:,:),gauge_c(:,:)
    real(8)::imaginary_residual
    allocate(overlap_c(size(overlap,1),size(overlap,2)),neighbor_c(size(neighbor_gauge,1), &
      size(neighbor_gauge,2)),basis_c(size(basis,1),size(basis,2)), &
      ww_c(size(ww_block,1),size(ww_block,2)),wp_c(size(wp_block,1),size(wp_block,2)), &
      self_c(size(face_self_block,1),size(face_self_block,2)), &
      face_c(size(face_neighbor_block,1),size(face_neighbor_block,2)), &
      gauge_c(size(overlap,2),size(overlap,1)))
    overlap_c=cmplx(overlap,0d0,8);neighbor_c=cmplx(neighbor_gauge,0d0,8)
    basis_c=cmplx(basis,0d0,8);ww_c=cmplx(ww_block,0d0,8);wp_c=cmplx(wp_block,0d0,8)
    self_c=cmplx(face_self_block,0d0,8);face_c=cmplx(face_neighbor_block,0d0,8)
    call stitch_and_apply_sawf_neighbor_pair(overlap_c,tolerance,neighbor_c,basis_c,ww_c,wp_c, &
      self_c,face_c,gauge_c,residual,ok,message)
    if(.not.ok)return
    imaginary_residual=maxval(abs(aimag(gauge_c)))
    if(imaginary_residual>tolerance)then
      ok=.false.; write(message,'(a,es13.5)') &
        'SAWF real GS gauge has non-negligible imaginary residual=',imaginary_residual
      return
    end if
    gauge=real(gauge_c,8);basis=real(basis_c,8);ww_block=real(ww_c,8);wp_block=real(wp_c,8)
    face_self_block=real(self_c,8);face_neighbor_block=real(face_c,8)
  end subroutine

  subroutine stitch_and_apply_sawf_neighbor_pair(overlap,tolerance,neighbor_gauge_unitary, &
      basis,ww_block,wp_block,face_self_block,face_neighbor_block,gauge_unitary,residual,ok,message)
    complex(8),intent(in)::overlap(:,:),neighbor_gauge_unitary(:,:)
    real(8),intent(in)::tolerance
    complex(8),intent(inout)::basis(:,:),ww_block(:,:),wp_block(:,:), &
      face_self_block(:,:),face_neighbor_block(:,:)
    complex(8),intent(out)::gauge_unitary(size(overlap,2),size(overlap,1))
    real(8),intent(out)::residual; logical,intent(out)::ok; character(*),intent(out)::message
    call stitch_sawf_neighbor_gauge(overlap,tolerance,gauge_unitary,residual,ok,message)
    if(.not.ok)return
    call apply_sawf_gauge_connection(gauge_unitary,neighbor_gauge_unitary,basis,ww_block, &
      wp_block,face_self_block,face_neighbor_block)
  end subroutine

  subroutine apply_sawf_gauge_connection(gauge_unitary,neighbor_gauge_unitary,basis,ww_block,wp_block, &
      face_self_block,face_neighbor_block)
    complex(8),intent(in)::gauge_unitary(:,:),neighbor_gauge_unitary(:,:)
    complex(8),intent(inout)::basis(:,:),ww_block(:,:),wp_block(:,:), &
      face_self_block(:,:),face_neighbor_block(:,:)
    complex(8),allocatable::left(:,:)
    basis=matmul(basis,gauge_unitary)
    left=conjg(transpose(gauge_unitary))
    ww_block=matmul(left,matmul(ww_block,gauge_unitary))
    wp_block=matmul(left,wp_block)
    face_self_block=matmul(left,matmul(face_self_block,gauge_unitary))
    face_neighbor_block=matmul(left,matmul(face_neighbor_block,neighbor_gauge_unitary))
  end subroutine

  subroutine materialize_sawf_local_bases(representative_basis,representative_index, &
      operation_index,independent_local,point_map,d_wann,local_basis,ok,message)
    complex(8),intent(in)::representative_basis(:,:,:),d_wann(:,:,:)
    integer,intent(in)::representative_index(:),operation_index(:),point_map(:,:)
    logical,intent(in)::independent_local(:)
    complex(8),intent(out)::local_basis(:,:,:)
    logical,intent(out)::ok; character(*),intent(out)::message
    integer::environment,source_point,target_point,representative,operation,npoint,nwann
    logical,allocatable::seen(:)
    ok=.false.; message=''; local_basis=(0d0,0d0)
    npoint=size(representative_basis,1); nwann=size(representative_basis,2)
    if(size(local_basis,1)/=npoint.or.size(local_basis,2)/=nwann.or. &
       size(local_basis,3)/=size(representative_index).or. &
       size(operation_index)/=size(representative_index).or. &
       size(independent_local)/=size(representative_index).or. &
       size(point_map,1)/=npoint.or.size(d_wann,1)/=nwann.or.size(d_wann,2)/=nwann)then
      message='SAWF local materialization dimensions are inconsistent'; return
    end if
    allocate(seen(npoint))
    do environment=1,size(representative_index)
      representative=representative_index(environment)
      if(representative<1.or.representative>size(representative_basis,3))then
        message='SAWF representative environment index is out of range'; return
      end if
      if(independent_local(environment))then
        local_basis(:,:,environment)=representative_basis(:,:,representative)
        cycle
      end if
      operation=operation_index(environment)
      if(operation<1.or.operation>size(point_map,2).or.operation>size(d_wann,3))then
        message='SAWF replicated environment operation index is invalid'; return
      end if
      seen=.false.
      do source_point=1,npoint
        target_point=point_map(source_point,operation)
        if(target_point<1.or.target_point>npoint.or.seen(target_point))then
          message='SAWF symmetry point map is not a permutation'; return
        end if
        seen(target_point)=.true.
        local_basis(target_point,:,environment)=matmul( &
          representative_basis(source_point,:,representative),d_wann(:,:,operation))
      end do
      if(.not.all(seen))then; message='SAWF symmetry point map is incomplete'; return; end if
    end do
    ok=.true.
  end subroutine

  subroutine clear_sawf_template_checkpoint(template)
    type(t_sawf_template_checkpoint),intent(inout)::template
    if(allocated(template%centers))deallocate(template%centers)
    if(allocated(template%spreads))deallocate(template%spreads)
    if(allocated(template%basis))deallocate(template%basis)
    if(allocated(template%buffer_basis))deallocate(template%buffer_basis)
    if(allocated(template%orbitals))deallocate(template%orbitals)
    if(allocated(template%d_band))deallocate(template%d_band)
    if(allocated(template%d_wann))deallocate(template%d_wann)
    if(allocated(template%gauge_unitary))deallocate(template%gauge_unitary)
    template%gauge_residual=huge(0d0)
  end subroutine

  subroutine write_sawf_template_checkpoint(filename,template,ok,message)
    character(*),intent(in)::filename
    type(t_sawf_template_checkpoint),intent(in)::template
    logical,intent(out)::ok; character(*),intent(out)::message
    integer::unit,ios,dims(15); character(512)::iomsg
    ok=.false.; message=''; dims=0
    if(.not.allocated(template%centers).or..not.allocated(template%spreads).or. &
       .not.allocated(template%basis).or..not.allocated(template%buffer_basis).or. &
       .not.allocated(template%d_band).or..not.allocated(template%d_wann).or. &
       .not.allocated(template%gauge_unitary))then
      message='SAWF template checkpoint is missing required basis metadata'; return
    end if
    dims=[size(template%centers,1),size(template%centers,2),size(template%spreads), &
      size(template%basis,1),size(template%basis,2),size(template%buffer_basis,1), &
      size(template%buffer_basis,2),0,0, &
      size(template%d_band,1),size(template%d_band,2),size(template%d_band,3), &
      size(template%d_wann,1),size(template%d_wann,2),size(template%d_wann,3)]
    if(allocated(template%orbitals))dims(8:9)=shape(template%orbitals)
    open(newunit=unit,file=filename,access='stream',form='unformatted',status='replace', &
      action='write',iostat=ios,iomsg=iomsg)
    if(ios/=0)then; message='SAWF template checkpoint open failed: '//trim(iomsg); return; end if
    write(unit,iostat=ios,iomsg=iomsg)sawf_template_magic,template%fingerprint%schema_version, &
      template%fingerprint%geometry,template%fingerprint%pseudopotential, &
      template%fingerprint%grid,template%fingerprint%band_window, &
      template%fingerprint%complete_projection_shell,template%fingerprint%symmetry, &
      template%fingerprint%buffer,template%fingerprint%generator,dims, &
      size(template%gauge_unitary,1),size(template%gauge_unitary,2),template%gauge_residual, &
      template%centers,template%spreads,template%basis,template%buffer_basis
    if(ios==0.and.allocated(template%orbitals))write(unit,iostat=ios,iomsg=iomsg)template%orbitals
    if(ios==0)write(unit,iostat=ios,iomsg=iomsg)template%d_band,template%d_wann,template%gauge_unitary
    close(unit)
    if(ios/=0)then; message='SAWF template checkpoint write failed: '//trim(iomsg); return; end if
    ok=.true.
  end subroutine

  subroutine read_sawf_template_checkpoint(filename,expected,template,reuse,ok,message)
    character(*),intent(in)::filename
    type(t_sawf_template_fingerprint),intent(in)::expected
    type(t_sawf_template_checkpoint),intent(inout)::template
    logical,intent(out)::reuse,ok; character(*),intent(out)::message
    type(t_sawf_template_fingerprint)::stored
    integer::unit,ios,dims(15),gauge_dims(2); character(512)::iomsg
    character(24)::magic; logical::exists
    call clear_sawf_template_checkpoint(template)
    ok=.false.; reuse=.false.; message=''; inquire(file=filename,exist=exists)
    if(.not.exists)then; ok=.true.; message='SAWF cache absent; regeneration required'; return; end if
    open(newunit=unit,file=filename,access='stream',form='unformatted',status='old', &
      action='read',iostat=ios,iomsg=iomsg)
    if(ios/=0)then; message='SAWF template checkpoint open failed: '//trim(iomsg); return; end if
    read(unit,iostat=ios,iomsg=iomsg)magic,stored%schema_version,stored%geometry, &
      stored%pseudopotential,stored%grid,stored%band_window,stored%complete_projection_shell, &
      stored%symmetry,stored%buffer,stored%generator,dims,gauge_dims,template%gauge_residual
    if(ios/=0.or.magic/=sawf_template_magic)then
      close(unit); message='SAWF template checkpoint header is invalid'; return
    end if
    if(any(dims<0).or.any(gauge_dims<=0).or.dims(1)/=3.or.dims(3)/=dims(2).or. &
       any(dims(4:7)<=0))then
      close(unit); message='SAWF template checkpoint dimensions are invalid'; return
    end if
    template%fingerprint=stored
    call validate_sawf_template_fingerprint(stored,expected,reuse)
    if(.not.reuse)then
      close(unit); ok=.true.; message='SAWF fingerprint mismatch; regeneration required'; return
    end if
    allocate(template%centers(dims(1),dims(2)),template%spreads(dims(3)), &
      template%basis(dims(4),dims(5)),template%buffer_basis(dims(6),dims(7)), &
      template%d_band(dims(10),dims(11),dims(12)),template%d_wann(dims(13),dims(14),dims(15)), &
      template%gauge_unitary(gauge_dims(1),gauge_dims(2)),stat=ios)
    if(ios/=0)then; close(unit); message='SAWF template checkpoint allocation failed'; return; end if
    if(dims(8)>0.and.dims(9)>0)allocate(template%orbitals(dims(8),dims(9)),stat=ios)
    if(ios==0)read(unit,iostat=ios,iomsg=iomsg)template%centers,template%spreads, &
      template%basis,template%buffer_basis
    if(ios==0.and.allocated(template%orbitals))read(unit,iostat=ios,iomsg=iomsg)template%orbitals
    if(ios==0)read(unit,iostat=ios,iomsg=iomsg)template%d_band,template%d_wann,template%gauge_unitary
    close(unit)
    if(ios/=0)then
      call clear_sawf_template_checkpoint(template)
      reuse=.false.; message='SAWF template checkpoint payload read failed: '//trim(iomsg); return
    end if
    ok=.true.
  end subroutine

  subroutine build_sawf_environment_orbits(equivalent,defect_intersects,orbit,regenerate,ok,message)
    logical,intent(in)::equivalent(:,:),defect_intersects(:)
    integer,intent(out)::orbit(size(defect_intersects))
    logical,intent(out)::regenerate(size(defect_intersects)),ok
    character(*),intent(out)::message
    integer::i,j,norb,head,tail,current
    integer,allocatable::queue(:)
    ok=.false.; message=''; orbit=0; regenerate=.false.
    if(size(equivalent,1)/=size(defect_intersects).or.size(equivalent,2)/=size(defect_intersects))then
      message='SAWF environment equivalence matrix has invalid dimensions'; return
    end if
    if(any(equivalent.neqv.transpose(equivalent)).or.any(.not.[(equivalent(i,i),i=1,size(orbit))]))then
      message='SAWF environment equivalence relation is not symmetric/reflexive'; return
    end if
    allocate(queue(size(orbit)));norb=0
    do i=1,size(orbit)
      if(orbit(i)==0)then
        norb=norb+1; orbit(i)=norb; regenerate(i)=.true.;head=1;tail=1;queue(1)=i
        do while(head<=tail)
          current=queue(head);head=head+1
          do j=1,size(orbit)
            if(equivalent(current,j).and.orbit(j)==0)then
              orbit(j)=norb;tail=tail+1;queue(tail)=j
            end if
          end do
        end do
      end if
    end do
    regenerate=regenerate.or.defect_intersects
    ok=.true.
  end subroutine

  subroutine validate_sawf_template_fingerprint(cached,current,reuse)
    type(t_sawf_template_fingerprint),intent(in)::cached,current
    logical,intent(out)::reuse
    reuse=cached%schema_version==current%schema_version .and. cached%geometry==current%geometry .and. &
      cached%pseudopotential==current%pseudopotential .and. cached%grid==current%grid .and. &
      cached%band_window==current%band_window .and. &
      cached%complete_projection_shell==current%complete_projection_shell .and. &
      cached%symmetry==current%symmetry .and. cached%buffer==current%buffer .and. &
      cached%generator==current%generator
  end subroutine

  subroutine validate_sawf_actual_group_operation(in_actual_group,is_parent_only,ok,message)
    logical,intent(in)::in_actual_group,is_parent_only
    logical,intent(out)::ok; character(*),intent(out)::message
    ok=in_actual_group.and..not.is_parent_only; message=''
    if(.not.ok)message='parent-crystal operation is not in the actual defective-supercell group'
  end subroutine

  subroutine replicate_sawf_operator_block(block,d_left,d_right,replica)
    complex(8),intent(in)::block(:,:),d_left(:,:),d_right(:,:)
    complex(8),intent(out)::replica(size(d_left,2),size(d_right,2))
    replica=matmul(conjg(transpose(d_left)),matmul(block,d_right))
  end subroutine

  subroutine stitch_sawf_neighbor_gauge(overlap,tolerance,gauge_unitary,residual,ok,message)
    complex(8),intent(in)::overlap(:,:); real(8),intent(in)::tolerance
    complex(8),intent(out)::gauge_unitary(size(overlap,2),size(overlap,1))
    real(8),intent(out)::residual; logical,intent(out)::ok; character(*),intent(out)::message
    complex(8),allocatable::a(:,:),u(:,:),vt(:,:),work(:),aligned(:,:)
    real(8),allocatable::singular(:),rwork(:); integer::n,info,lwork,i
    ok=.false.; message=''; residual=huge(0d0); gauge_unitary=(0d0,0d0)
    if(size(overlap,1)/=size(overlap,2).or.tolerance<=0d0)then
      message='SAWF neighbor overlap has ambiguous orbital counts'; return
    end if
    n=size(overlap,1); lwork=max(1,8*n)
    allocate(a(n,n),u(n,n),vt(n,n),work(lwork),singular(n),rwork(max(1,5*n)),aligned(n,n))
    a=overlap
    call zgesvd('A','A',n,n,a,n,singular,u,n,vt,n,work,lwork,rwork,info)
    if(info/=0)then; message='SAWF neighbor overlap SVD failed'; return; end if
    if(minval(singular)<=tolerance)then; message='SAWF neighbor overlap is rank deficient'; return; end if
    gauge_unitary=matmul(conjg(transpose(vt)),conjg(transpose(u)))
    aligned=matmul(overlap,gauge_unitary)
    do i=1,n; aligned(i,i)=aligned(i,i)-(1d0,0d0); end do
    residual=sqrt(sum(abs(aligned)**2))/sqrt(dble(n))
    ok=residual<=tolerance
    if(.not.ok)message='SAWF post-alignment gauge residual exceeds tolerance'
  end subroutine

  subroutine validate_sawf_buffer_convergence(center_residual,projector_residual,overlap_residual, &
      ww_residual,wp_residual,face_residual,tolerance,ok,message)
    real(8),intent(in)::center_residual,projector_residual,overlap_residual,ww_residual,wp_residual
    real(8),intent(in)::face_residual(:),tolerance; logical,intent(out)::ok; character(*),intent(out)::message
    ok=max(center_residual,projector_residual,overlap_residual,ww_residual,wp_residual, &
      maxval(face_residual))<=tolerance; message=''
    if(.not.ok)message='SAWF buffer convergence failed for orbital or operator/DG face data'
  end subroutine

  subroutine validate_sawf_global_local_equivalence(occupied_projector_residual,overlap_residual, &
      fixed_hamiltonian_residual,face_residual,tolerance,ok,message)
    real(8),intent(in)::occupied_projector_residual,overlap_residual,fixed_hamiltonian_residual
    real(8),intent(in)::face_residual(:),tolerance; logical,intent(out)::ok; character(*),intent(out)::message
    ok=max(occupied_projector_residual,overlap_residual,fixed_hamiltonian_residual, &
      maxval(face_residual))<=tolerance; message=''
    if(.not.ok)message='hierarchical SAWF is not operator-equivalent to monolithic global SAWF'
  end subroutine
end module lcfo_wannier_sawf_templates
