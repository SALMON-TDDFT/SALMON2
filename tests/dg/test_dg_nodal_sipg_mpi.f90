program test_dg_nodal_sipg_mpi
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use dg_nodal_sipg
  implicit none
  type(s_dg_nodal_sipg_face), allocatable :: faces(:), relabeled(:)
  type(s_dg_nodal_sipg_action), allocatable :: action(:), saved(:)
  type(s_dg_nodal_sipg_diagnostics) :: diagnostics
  complex(8) :: x(4), y(4), hx(4), hy(4), dense(4,4), lhs, rhs
  real(8) :: eta, h
  integer :: ierr, rank, nproc
  logical :: ok
  character(len=256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc /= 2) error stop 'test requires two ranks'
  eta=8.0d0; h=0.5d0
  dense=cmplx(0.0d0,0.0d0,8)
  dense(1,:)=[eta/h,-eta/h,-0.5d0,-0.5d0]
  dense(2,:)=[-eta/h,eta/h,0.5d0,0.5d0]
  dense(3,:)=[-0.5d0,0.5d0,0.0d0,0.0d0]
  dense(4,:)=[-0.5d0,0.5d0,0.0d0,0.0d0]
  x=[cmplx(1.0d0,0.2d0,8),cmplx(-0.3d0,0.4d0,8),cmplx(0.7d0,-0.1d0,8),cmplx(-0.2d0,0.5d0,8)]
  y=[cmplx(0.2d0,-0.6d0,8),cmplx(0.8d0,0.1d0,8),cmplx(-0.4d0,0.3d0,8),cmplx(0.5d0,0.7d0,8)]
  call evaluate(x,h,eta,hx)
  call evaluate(y,h,eta,hy)
  call require(maxval(abs(hx-matmul(dense,x)))<1.0d-13,'independent dense SIPG action')
  lhs=sum(conjg(x)*hy); rhs=conjg(sum(conjg(y)*hx))
  call require(abs(lhs-rhs)<1.0d-13,'Hermiticity')
  call require(abs(hx(1)+hx(2))<1.0d-13,'left-right consistency cancellation')

  allocate(faces(1))
  call fill_face(faces(1),rank)
  call apply_dg_nodal_sipg_faces_mpi(faces,7_8,eta,MPI_COMM_WORLD,action,diagnostics,ok,message)
  call require(ok,'valid physical face batch')
  call require(diagnostics%ownership_count==1,'single canonical owner per physical face')
  call require(diagnostics%hermiticity_defect<1.0d-13,'face Hermiticity diagnostic')
  call require(diagnostics%internal_cancellation_defect<1.0d-13,'face cancellation diagnostic')
  call require(diagnostics%jump_norm>0.0d0 .and. diagnostics%penalty_energy>0.0d0,'jump diagnostics')
  call require(diagnostics%trace_epoch==7_8,'trace epoch diagnostic')
  call require(faces(1)%periodic_shift(1)/=0,'physical-supercell periodic face')
  call require(actions_agree(action(1)),'owner action delivered to nonowner rank')

  saved=action
  relabeled=faces
  if(rank==1)then
    relabeled(1)%fragment_minus=faces(1)%fragment_plus
    relabeled(1)%fragment_plus=faces(1)%fragment_minus
    relabeled(1)%normal=-faces(1)%normal
    relabeled(1)%periodic_shift=-faces(1)%periodic_shift
    relabeled(1)%u_minus=faces(1)%u_plus
    relabeled(1)%u_plus=faces(1)%u_minus
    relabeled(1)%dn_minus=-faces(1)%dn_plus
    relabeled(1)%dn_plus=-faces(1)%dn_minus
  endif
  call apply_dg_nodal_sipg_faces_mpi(relabeled,7_8,eta,MPI_COMM_WORLD,action,diagnostics,ok,message)
  call require(.not.ok.and.actions_equal(action,saved),'rank-disagreeing canonical orientation rejection')
  relabeled=faces
  if(rank==1)relabeled(1)%h_normal=2.0d0*faces(1)%h_normal
  call apply_dg_nodal_sipg_faces_mpi(relabeled,7_8,eta,MPI_COMM_WORLD,action,diagnostics,ok,message)
  call require(.not.ok.and.actions_equal(action,saved),'rank-disagreeing face metadata rejection')

  relabeled=faces
  relabeled(1)%fragment_minus=faces(1)%fragment_plus
  relabeled(1)%fragment_plus=faces(1)%fragment_minus
  relabeled(1)%u_minus=faces(1)%u_plus
  relabeled(1)%u_plus=faces(1)%u_minus
  relabeled(1)%dn_minus=-faces(1)%dn_plus
  relabeled(1)%dn_plus=-faces(1)%dn_minus
  relabeled(1)%normal=-faces(1)%normal
  call apply_dg_nodal_sipg_faces_mpi(relabeled,7_8,eta,MPI_COMM_WORLD,action,diagnostics,ok,message)
  call require(ok .and. abs(action(1)%total_value(1)-saved(1)%total_value(2))<1.0d-13, &
               'fragment relabeling invariance')

  relabeled=faces
  relabeled(1)%physical_face=.false.
  relabeled(1)%physical_supercell_periodic=.false.
  relabeled(1)%auxiliary_periodic_wrap=.true.
  call apply_dg_nodal_sipg_faces_mpi(relabeled,7_8,eta,MPI_COMM_WORLD,action,diagnostics,ok,message)
  call require(ok .and. maxval(abs(action(1)%total_value))==0.0d0 .and. &
               diagnostics%ownership_count==0,'auxiliary local-box wrap excluded')

  action=saved
  relabeled=faces
  relabeled(1)%normal=[2.0d0,0.0d0,0.0d0]
  call apply_dg_nodal_sipg_faces_mpi(relabeled,7_8,eta,MPI_COMM_WORLD,action,diagnostics,ok,message)
  call require(.not.ok .and. actions_equal(action,saved),'invalid normal transactional failure')
  relabeled=faces
  relabeled(1)%normal=[0.0d0,1.0d0,0.0d0]
  call apply_dg_nodal_sipg_faces_mpi(relabeled,7_8,eta,MPI_COMM_WORLD,action,diagnostics,ok,message)
  call require(.not.ok .and. actions_equal(action,saved),'axis-normal mismatch rejection')
  relabeled=faces
  relabeled(1)%periodic_shift=[0,0,0]
  call apply_dg_nodal_sipg_faces_mpi(relabeled,7_8,eta,MPI_COMM_WORLD,action,diagnostics,ok,message)
  call require(.not.ok .and. actions_equal(action,saved),'physical periodic shift rejection')
  relabeled=faces
  relabeled(1)%trace_epoch=6_8
  call apply_dg_nodal_sipg_faces_mpi(relabeled,7_8,eta,MPI_COMM_WORLD,action,diagnostics,ok,message)
  call require(.not.ok .and. actions_equal(action,saved),'stale trace transactional failure')
  relabeled=faces
  relabeled(1)%u_minus=cmplx(ieee_value(0.0d0,ieee_quiet_nan),0.0d0,8)
  call apply_dg_nodal_sipg_faces_mpi(relabeled,7_8,eta,MPI_COMM_WORLD,action,diagnostics,ok,message)
  call require(.not.ok .and. actions_equal(action,saved),'nonfinite transactional failure')
  relabeled=faces
  relabeled(1)%canonical_owner=.false.
  call apply_dg_nodal_sipg_faces_mpi(relabeled,7_8,eta,MPI_COMM_WORLD,action,diagnostics,ok,message)
  call require(.not.ok .and. actions_equal(action,saved),'missing canonical owner rejection')
  relabeled=faces
  relabeled(1)%canonical_owner=.true.
  call apply_dg_nodal_sipg_faces_mpi(relabeled,7_8,eta,MPI_COMM_WORLD,action,diagnostics,ok,message)
  call require(.not.ok .and. actions_equal(action,saved),'duplicate canonical ownership rejection')

  if(rank==0) print '(a)','PASS analytic complete Hermitian nodal SIPG face operator'
  call MPI_Finalize(ierr)
contains
  subroutine evaluate(values,h_local,eta_local,result)
    complex(8),intent(in)::values(4)
    real(8),intent(in)::h_local,eta_local
    complex(8),intent(out)::result(4)
    type(s_dg_nodal_sipg_action)::piece
    integer::info
    call evaluate_dg_nodal_sipg_face(values(1),values(2),values(3),values(4),h_local,1.0d0,eta_local,piece,info)
    if(info/=0)error stop 'analytic evaluation failed'
    result=[piece%total_value,piece%total_normal]
    call require(maxval(abs(piece%consistency_value-[-0.5d0*(values(3)+values(4)), &
      0.5d0*(values(3)+values(4))]))<1.0d-13,'consistency term')
    call require(maxval(abs(piece%symmetry_normal-[-0.5d0*(values(1)-values(2)), &
      -0.5d0*(values(1)-values(2))]))<1.0d-13,'symmetric term')
    call require(maxval(abs(piece%penalty_value-[eta_local/h_local*(values(1)-values(2)), &
      -eta_local/h_local*(values(1)-values(2))]))<1.0d-13,'penalty term')
  end subroutine evaluate
  subroutine fill_face(face,irank)
    type(s_dg_nodal_sipg_face),intent(out)::face
    integer,intent(in)::irank
    face%global_face_id=99_8
    face%fragment_minus=1;face%fragment_plus=2
    face%axis=1;face%normal=[1.0d0,0.0d0,0.0d0]
    face%periodic_shift=[1,0,0]
    face%h_normal=h;face%face_weight=1.0d0;face%canonical_owner=irank==0
    face%physical_supercell_periodic=.true.
    face%physical_face=.true.;face%auxiliary_periodic_wrap=.false.;face%trace_epoch=7_8
    face%u_minus=x(1);face%u_plus=x(2)
    face%dn_minus=x(3);face%dn_plus=x(4)
  end subroutine fill_face
  logical function actions_agree(value)
    type(s_dg_nodal_sipg_action),intent(in)::value
    complex(8)::local(4),global_sum(4)
    integer::mpi_error
    local=[value%total_value,value%total_normal]
    call MPI_Allreduce(local,global_sum,4,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpi_error)
    actions_agree=maxval(abs(local-global_sum/dble(nproc)))<1.0d-14
  end function actions_agree
  logical function actions_equal(left,right)
    type(s_dg_nodal_sipg_action),intent(in)::left(:),right(:)
    integer::i
    actions_equal=size(left)==size(right)
    if(.not.actions_equal)return
    do i=1,size(left)
      actions_equal=actions_equal.and.all(left(i)%total_value==right(i)%total_value).and. &
        all(left(i)%total_normal==right(i)%total_normal)
    enddo
  end function actions_equal
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    logical::agreed
    integer::mpi_error
    call MPI_Allreduce(condition,agreed,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,mpi_error)
    if(.not.agreed)then
      if(rank==0)write(*,'(a,1x,a)')'FAIL',trim(label)
      call MPI_Abort(MPI_COMM_WORLD,1,mpi_error)
    endif
  end subroutine require
end program test_dg_nodal_sipg_mpi
