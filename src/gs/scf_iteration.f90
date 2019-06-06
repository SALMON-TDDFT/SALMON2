!
!  Copyright 2019 SALMON developers
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
module scf_iteration_sub
  implicit none

contains

subroutine scf_iteration(mg,system,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,ilsda,nproc_ob,iparaway_ob, &
               rxk_ob,rhxk_ob,rgk_ob,rpk_ob,   &
               info_ob,bnmat,cnmat,ppg,vlocal)
  use inputoutput, only: ispin
  use structures
  use gscg_sub
  implicit none

  type(s_rgrid),         intent(in)    :: mg
  type(s_system),        intent(in)    :: system
  type(s_wf_info),       intent(in)    :: info
  type(s_wavefunction),  intent(inout) :: spsi
  type(s_stencil),       intent(in)    :: stencil
  type(s_sendrecv_grid), intent(inout) :: srg_ob_1
  type(s_pp_grid),       intent(in)    :: ppg
  integer,               intent(inout) :: iflag
  integer,               intent(in)    :: itotmst
  integer,               intent(in)    :: mst(2)
  integer,               intent(in)    :: ilsda
  integer,               intent(in)    :: nproc_ob
  integer,               intent(in)    :: iparaway_ob
  real(8),               intent(inout) :: rxk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  real(8),               intent(inout) :: rhxk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  real(8),               intent(inout) :: rgk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  real(8),               intent(inout) :: rpk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:system%nspin*info%numo)
  type(s_wf_info),       intent(in)    :: info_ob
  real(8),               intent(in)    :: cnmat(0:12,0:12),bnmat(0:12,0:12)
  real(8),               intent(in)    :: vlocal(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),ispin+1)

  call sgscg(mg,system,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,ilsda,nproc_ob,iparaway_ob, &
             rxk_ob,rhxk_ob,rgk_ob,rpk_ob,   &
             info_ob,bnmat,cnmat,ppg,vlocal)

end subroutine scf_iteration

end module scf_iteration_sub
