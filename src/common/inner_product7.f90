!
!  Copyright 2018 SALMON developers
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
!=======================================================================
subroutine inner_product7(mg,iparaway_ob,itotmst,mst,iobnum,rmatbox1,rmatbox2,rbox2,elp3,hvol)
  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_group_korbital
  use salmon_communication, only: comm_summation
  use misc_routines, only: get_wtime
  implicit none
  type(s_rgrid),intent(in) :: mg
  integer,intent(in)  :: iparaway_ob
  integer,intent(in)  :: itotmst
  integer,intent(in)  :: mst(2)
  integer,intent(in)  :: iobnum
  real(8),intent(in)  :: rmatbox1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:iobnum)
  real(8),intent(in)  :: rmatbox2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:iobnum)
  real(8),intent(out) :: rbox2(itotmst)
  real(8),intent(out) :: elp3(3000)
  real(8),intent(in)  :: hvol
  integer :: ix,iy,iz,iob,iob_allob
  real(8) :: rbox,rbox1(itotmst)
  
  rbox1(:)=0.d0
 
  do iob=1,iobnum
    call calc_allob(iob,iob_allob,iparaway_ob,itotmst,mst,iobnum)
    rbox=0.d0
    !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : rbox)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rbox=rbox+rmatbox1(ix,iy,iz,iob)*rmatbox2(ix,iy,iz,iob)
    end do
    end do
    end do
    rbox1(iob_allob)=rbox*hvol
  end do
  
  elp3(186)=get_wtime()
  call comm_summation(rbox1,rbox2,itotmst,nproc_group_korbital)
  elp3(187)=get_wtime()
  elp3(190)=elp3(190)+elp3(187)-elp3(186)
  
end subroutine inner_product7
