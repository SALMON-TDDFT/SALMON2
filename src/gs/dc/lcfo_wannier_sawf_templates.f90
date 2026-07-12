module lcfo_wannier_sawf_templates
  implicit none
  private
  integer, parameter :: sawf_template_schema_version=1

  type, public :: t_sawf_template_fingerprint
    character(256) :: geometry='', pseudopotential='', grid='', band_window=''
    character(256) :: complete_projection_shell='', symmetry='', buffer='', generator=''
    integer :: schema_version=sawf_template_schema_version
  end type

  public :: build_sawf_environment_orbits, validate_sawf_template_fingerprint
  public :: validate_sawf_actual_group_operation, replicate_sawf_operator_block
  public :: stitch_sawf_neighbor_gauge, validate_sawf_buffer_convergence
  public :: validate_sawf_global_local_equivalence

contains
  subroutine build_sawf_environment_orbits(equivalent,defect_intersects,orbit,regenerate,ok,message)
    logical,intent(in)::equivalent(:,:),defect_intersects(:)
    integer,intent(out)::orbit(size(defect_intersects))
    logical,intent(out)::regenerate(size(defect_intersects)),ok
    character(*),intent(out)::message
    integer::i,j,norb
    ok=.false.; message=''; orbit=0; regenerate=defect_intersects
    if(size(equivalent,1)/=size(defect_intersects).or.size(equivalent,2)/=size(defect_intersects))then
      message='SAWF environment equivalence matrix has invalid dimensions'; return
    end if
    if(any(equivalent.neqv.transpose(equivalent)).or.any(.not.[(equivalent(i,i),i=1,size(orbit))]))then
      message='SAWF environment equivalence relation is not symmetric/reflexive'; return
    end if
    norb=0
    do i=1,size(orbit)
      if(orbit(i)==0)then
        norb=norb+1; orbit(i)=norb
        do j=i+1,size(orbit)
          if(equivalent(i,j))orbit(j)=norb
        end do
      end if
    end do
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
