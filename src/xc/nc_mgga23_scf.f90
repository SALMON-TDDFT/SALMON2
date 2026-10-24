! SCF-only damping; RT always uses the instantaneous variational operator.
module nc_mgga23_scf
 use structures,only:s_nc_mgga23_payload
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 use communication,only:comm_get_max,comm_summation
 implicit none
 private
 public::nc23_scf_prepare,nc23_scf_residual,nc23_spin_density_residual
 interface finite_array
  module procedure finite4,finite5,finite_complex5
 end interface
contains
 subroutine nc23_scf_prepare(p,history,alpha)
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  type(s_nc_mgga23_payload),intent(inout)::p,history
  real(8),intent(in)::alpha
  if(.not.p%active)return
  if(.not.ieee_is_finite(alpha))error stop 'nc_mgga23 SCF: nonfinite mixing rate'
  if(alpha<=0d0.or.alpha>1d0)error stop 'nc_mgga23 SCF: mixing rate must be in (0,1]'
  if(history%active)then
   p%local=(1d0-alpha)*history%local+alpha*p%local
   p%tau=(1d0-alpha)*history%tau+alpha*p%tau
   p%current=(1d0-alpha)*history%current+alpha*p%current
  endif
  ! Energy band contraction and explicit core field derivative remain physical.
  history=p
 end subroutine
 subroutine nc23_scf_residual(p,history,comm,residual)
  implicit none
  type(s_nc_mgga23_payload),intent(in)::p,history
  integer,intent(in)::comm
  real(8),intent(out)::residual
  real(8)::v(6),g(6)
  integer::bad,global_bad,i
  residual=0d0
  if(.not.p%active.or..not.history%active)return
  bad=0
  if(.not.finite_array(p%local))bad=1
  if(.not.finite_array(p%tau))bad=1
  if(.not.finite_array(p%current))bad=1
  if(.not.finite_array(history%local))bad=1
  if(.not.finite_array(history%tau))bad=1
  if(.not.finite_array(history%current))bad=1
  call comm_summation(bad,global_bad,comm)
  if(global_bad/=0)error stop 'nc_mgga23 SCF: nonfinite coefficient'
  v(1)=maxval(abs(p%local-history%local))
  v(2)=maxval(abs(p%tau-history%tau))
  v(3)=maxval(abs(p%current-history%current))
  v(4)=max(1d0,maxval(abs(p%local)),maxval(abs(history%local)))
  v(5)=max(1d0,maxval(abs(p%tau)),maxval(abs(history%tau)))
  v(6)=max(1d0,maxval(abs(p%current)),maxval(abs(history%current)))
  do i=1,6
   if(.not.ieee_is_finite(v(i)))bad=1
  enddo
  call comm_summation(bad,global_bad,comm)
  if(global_bad/=0)error stop 'nc_mgga23 SCF: nonfinite residual'
  call comm_get_max(v,g,6,comm)
  residual=maxval(g(1:3)/g(4:6))
 end subroutine
 subroutine nc23_spin_density_residual(now,before,hvol,alpha,comm,residual)
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  complex(8),intent(in)::now(:,:,:,:,:),before(:,:,:,:,:)
  real(8),intent(in)::hvol,alpha
  integer,intent(in)::comm
  real(8),intent(out)::residual
  real(8)::v(2),g(2)
  integer::bad,global_bad,i
  if(.not.ieee_is_finite(alpha))error stop 'SOC SCF: nonfinite mixing rate'
  if(alpha<=0d0.or.alpha>1d0)error stop 'SOC SCF: mixing rate must be in (0,1]'
  ! For simple_dm only: undo its exact linear alpha suppression. The full
  ! density matrix includes charge-invisible transverse spin changes.
  bad=0
  if(.not.finite_array(now))bad=1
  if(.not.finite_array(before))bad=1
  if(.not.ieee_is_finite(hvol))bad=1
  if(hvol<=0d0)bad=1
  call comm_summation(bad,global_bad,comm)
  if(global_bad/=0)error stop 'SOC SCF: nonfinite density or invalid volume'
  v(1)=sum(abs(now-before))*hvol/alpha
  v(2)=sum(real(now(:,:,:,1,1)+now(:,:,:,2,2),8))*hvol
  do i=1,2
   if(.not.ieee_is_finite(v(i)))bad=1
  enddo
  call comm_summation(bad,global_bad,comm)
  if(global_bad/=0)error stop 'SOC SCF: nonfinite density residual'
  call comm_summation(v,g,2,comm)
  do i=1,2
   if(.not.ieee_is_finite(g(i)))error stop 'SOC SCF: nonfinite global density residual'
  enddo
  if(g(2)<=0d0)error stop 'SOC SCF: nonpositive electron count'
  residual=g(1)/g(2)
 end subroutine
 logical function finite4(a) result(ok)
  implicit none
  real(8),intent(in)::a(:,:,:,:)
  integer::i,j,k,l
  ok=.true.
  do l=1,size(a,4)
  do k=1,size(a,3)
  do j=1,size(a,2)
  do i=1,size(a,1)
   if(.not.ieee_is_finite(a(i,j,k,l)))ok=.false.
  enddo
  enddo
  enddo
  enddo
 end function
 logical function finite5(a) result(ok)
  implicit none
  real(8),intent(in)::a(:,:,:,:,:)
  integer::i,j,k,l,m
  ok=.true.
  do m=1,size(a,5)
  do l=1,size(a,4)
  do k=1,size(a,3)
  do j=1,size(a,2)
  do i=1,size(a,1)
   if(.not.ieee_is_finite(a(i,j,k,l,m)))ok=.false.
  enddo
  enddo
  enddo
  enddo
  enddo
 end function
 logical function finite_complex5(a) result(ok)
  implicit none
  complex(8),intent(in)::a(:,:,:,:,:)
  integer::i,j,k,l,m
  ok=.true.
  do m=1,size(a,5)
  do l=1,size(a,4)
  do k=1,size(a,3)
  do j=1,size(a,2)
  do i=1,size(a,1)
   if(.not.ieee_is_finite(real(a(i,j,k,l,m),8)))ok=.false.
   if(.not.ieee_is_finite(aimag(a(i,j,k,l,m))))ok=.false.
  enddo
  enddo
  enddo
  enddo
  enddo
 end function
end module
