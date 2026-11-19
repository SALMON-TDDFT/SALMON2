!
!  Copyright 2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
module automatic_distribution
  implicit none
  private
  public :: select_numcpu_dft
contains
! Select a decomposition independently of MPI so the policy can be tested directly.
subroutine select_numcpu_dft(nproc,nk,no,ng,nd,allow_rspace,use_ffte,pk,po,pr,found)
  implicit none
  integer,intent(in) :: nproc,nk,no,ng(3),nd
  logical,intent(in) :: allow_rspace,use_ffte
  integer,intent(out) :: pk,po,pr(3)
  logical,intent(out) :: found
  integer :: k,r,left,o,trial(3),best_r,best_grid(3)
  real(8) :: volume,local_volume,target
  logical :: grid_ok,small

  found=.false.
  pk=1; po=1; pr=1
  if(nproc<1.or.nk<1.or.no<1.or.any(ng<1)) return
  volume=product(real(ng,8))
  small=volume<=4096d0
  target=4096d0
  if(small) target=huge(target)
  ! Maximize k parallelism, but backtrack if the remaining ranks cannot fit.
  do k=min(nproc,nk),1,-1
    if(mod(nproc,k)/=0) cycle
    left=nproc/k
    best_r=0
    best_grid=1
    do r=1,left
      if(mod(left,r)/=0) cycle
      o=left/r
      if(o>no) cycle
      if(.not.allow_rspace.and.r/=1) cycle
      trial=1
      call distribute_grid_factors(r,ng,nd,use_ffte,target,trial,grid_ok)
      if(.not.grid_ok.and..not.small) then
        trial=1
        call distribute_grid_factors(r,ng,nd,use_ffte,huge(target),trial,grid_ok)
      end if
      if(.not.grid_ok) cycle
      best_r=r
      best_grid=trial
      local_volume=product(real((ng-1)/trial+1,8))
      ! Small cells maximize orbital parallelism. Large cells use the fewest
      ! spatial ranks reaching the target, or the largest feasible split.
      if(small.or.local_volume<=4096d0.or..not.allow_rspace) exit
    end do
    if(best_r==0) cycle
    pk=k
    po=left/best_r
    pr=best_grid
    found=.true.
    return
  end do
end subroutine select_numcpu_dft

recursive subroutine distribute_grid_factors(remaining,ng,nd,use_ffte,target,pr,found)
  implicit none
  integer,intent(in) :: remaining,ng(3),nd
  logical,intent(in) :: use_ffte
  real(8),intent(in) :: target
  integer,intent(inout) :: pr(3)
  logical,intent(out) :: found
  integer :: factor,j,axis,local_size(3),trial(3)

  found=remaining==1
  if(found) then
    found=product(real((ng-1)/pr+1,8))<=target
    return
  end if
  factor=2
  do while(factor<=remaining/factor)
    if(mod(remaining,factor)==0) exit
    factor=factor+1
  end do
  if(mod(remaining,factor)/=0) factor=remaining
  local_size=(ng-1)/pr+1
  ! Try the longest current local direction first (x, y, z for ties).
  ! Backtracking preserves valid layouts when a greedy choice hits a constraint.
  do j=1,3
    axis=maxloc(local_size,dim=1)
    local_size(axis)=-1
    trial=pr
    trial(axis)=trial(axis)*factor
    if(.not.valid_grid_split(ng,nd,use_ffte,trial)) cycle
    call distribute_grid_factors(remaining/factor,ng,nd,use_ffte,target,trial,found)
    if(found) then
      pr=trial
      return
    end if
  end do
end subroutine distribute_grid_factors

logical function valid_grid_split(ng,nd,use_ffte,pr) result(ok)
  implicit none
  integer,intent(in) :: ng(3),nd,pr(3)
  logical,intent(in) :: use_ffte
  integer :: j,width,last

  ok=.false.
  do j=1,3
    if(pr(j)==1) cycle
    width=(ng(j)-1)/pr(j)+1
    last=ng(j)-width*(pr(j)-1)
    ! init_grid_parallel uses ceiling-sized blocks, including a shorter last
    ! block. Every split block must supply the finite-difference halo.
    if(last<nd) return
  end do
  if(use_ffte) then
    if(mod(ng(1),pr(2))/=0.or.mod(ng(2),pr(2))/=0) return
    if(mod(ng(2),pr(3))/=0.or.mod(ng(3),pr(3))/=0) return
  end if
  ok=.true.
end function valid_grid_split

end module automatic_distribution
