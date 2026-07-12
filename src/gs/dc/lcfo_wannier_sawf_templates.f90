module lcfo_wannier_sawf_templates
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
    complex(8), allocatable :: orbitals(:,:),d_band(:,:,:),d_wann(:,:,:)
    complex(8), allocatable :: gauge_unitary(:,:)
    real(8) :: gauge_residual=huge(0d0)
  end type

  public :: build_sawf_environment_orbits, validate_sawf_template_fingerprint
  public :: validate_sawf_actual_group_operation, replicate_sawf_operator_block
  public :: stitch_sawf_neighbor_gauge, validate_sawf_buffer_convergence
  public :: validate_sawf_global_local_equivalence
  public :: write_sawf_template_checkpoint, read_sawf_template_checkpoint

contains
  subroutine clear_sawf_template_checkpoint(template)
    type(t_sawf_template_checkpoint),intent(inout)::template
    if(allocated(template%centers))deallocate(template%centers)
    if(allocated(template%spreads))deallocate(template%spreads)
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
    integer::unit,ios,dims(11); character(512)::iomsg
    ok=.false.; message=''; dims=0
    if(.not.allocated(template%centers).or..not.allocated(template%spreads).or. &
       .not.allocated(template%d_band).or..not.allocated(template%d_wann).or. &
       .not.allocated(template%gauge_unitary))then
      message='SAWF template checkpoint is missing required basis metadata'; return
    end if
    dims=[size(template%centers,1),size(template%centers,2),size(template%spreads),0,0, &
      size(template%d_band,1),size(template%d_band,2),size(template%d_band,3), &
      size(template%d_wann,1),size(template%d_wann,2),size(template%d_wann,3)]
    if(allocated(template%orbitals))dims(4:5)=shape(template%orbitals)
    open(newunit=unit,file=filename,access='stream',form='unformatted',status='replace', &
      action='write',iostat=ios,iomsg=iomsg)
    if(ios/=0)then; message='SAWF template checkpoint open failed: '//trim(iomsg); return; end if
    write(unit,iostat=ios,iomsg=iomsg)sawf_template_magic,template%fingerprint%schema_version, &
      template%fingerprint%geometry,template%fingerprint%pseudopotential, &
      template%fingerprint%grid,template%fingerprint%band_window, &
      template%fingerprint%complete_projection_shell,template%fingerprint%symmetry, &
      template%fingerprint%buffer,template%fingerprint%generator,dims, &
      size(template%gauge_unitary,1),size(template%gauge_unitary,2),template%gauge_residual, &
      template%centers,template%spreads
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
    integer::unit,ios,dims(11),gauge_dims(2); character(512)::iomsg
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
    if(any(dims<0).or.any(gauge_dims<=0).or.dims(1)/=3.or.dims(3)/=dims(2))then
      close(unit); message='SAWF template checkpoint dimensions are invalid'; return
    end if
    template%fingerprint=stored
    call validate_sawf_template_fingerprint(stored,expected,reuse)
    if(.not.reuse)then
      close(unit); ok=.true.; message='SAWF fingerprint mismatch; regeneration required'; return
    end if
    allocate(template%centers(dims(1),dims(2)),template%spreads(dims(3)), &
      template%d_band(dims(6),dims(7),dims(8)),template%d_wann(dims(9),dims(10),dims(11)), &
      template%gauge_unitary(gauge_dims(1),gauge_dims(2)),stat=ios)
    if(ios/=0)then; close(unit); message='SAWF template checkpoint allocation failed'; return; end if
    if(dims(4)>0.and.dims(5)>0)allocate(template%orbitals(dims(4),dims(5)),stat=ios)
    if(ios==0)read(unit,iostat=ios,iomsg=iomsg)template%centers,template%spreads
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
    integer::i,j,norb
    ok=.false.; message=''; orbit=0; regenerate=.false.
    if(size(equivalent,1)/=size(defect_intersects).or.size(equivalent,2)/=size(defect_intersects))then
      message='SAWF environment equivalence matrix has invalid dimensions'; return
    end if
    if(any(equivalent.neqv.transpose(equivalent)).or.any(.not.[(equivalent(i,i),i=1,size(orbit))]))then
      message='SAWF environment equivalence relation is not symmetric/reflexive'; return
    end if
    norb=0
    do i=1,size(orbit)
      if(orbit(i)==0)then
        norb=norb+1; orbit(i)=norb; regenerate(i)=.true.
        do j=i+1,size(orbit)
          if(equivalent(i,j))orbit(j)=norb
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
