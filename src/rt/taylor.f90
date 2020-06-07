!
!  Copyright 2019-2020 SALMON developers
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
module taylor_sub
  implicit none

contains

subroutine taylor(mg,system,info,stencil,srg,tspsi_in,tspsi_out,sshtpsi,   &
                  ppg,V_local,rt)
  use inputoutput, only: n_hamil
  use structures
  use hamiltonian, only: hpsi
  use sendrecv_grid, only: s_sendrecv_grid
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_stencil),intent(in) :: stencil
  type(s_sendrecv_grid),intent(inout) :: srg
  type(s_orbital),intent(inout) :: tspsi_in
  type(s_orbital),intent(inout) :: tspsi_out
  type(s_orbital),intent(inout) :: sshtpsi
  type(s_pp_grid),intent(in) :: ppg
  type(s_scalar) ,intent(in) :: V_local(system%nspin)
  type(s_rt),     intent(in) :: rt
  integer :: nn,ix,iy,iz
  integer :: ik,io,is,nspin
  nspin = system%nspin

  do nn=1,n_hamil
    if(mod(nn,2)==1)then
      call hpsi(tspsi_in,sshtpsi,info,mg,V_local,system,stencil,srg,ppg)
      if (nn==1) then
!$OMP parallel do collapse(5) private(ik,io,is,iz,iy,ix)
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          do is=1,nspin
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              tspsi_out%zwf(ix,iy,iz,is,io,ik,1)=tspsi_in%zwf(ix,iy,iz,is,io,ik,1)+ &
                                                   rt%zc(nn)*sshtpsi%zwf(ix,iy,iz,is,io,ik,1)
            end do
            end do
            end do
          end do
        end do
        end do
      else
!$OMP parallel do collapse(5) private(ik,io,is,iz,iy,ix)
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          do is=1,nspin
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              tspsi_out%zwf(ix,iy,iz,is,io,ik,1)=tspsi_out%zwf(ix,iy,iz,is,io,ik,1)+ &
                                                   rt%zc(nn)*sshtpsi%zwf(ix,iy,iz,is,io,ik,1)
            end do
            end do
            end do
          end do
        end do
        end do
      end if
    else
      call hpsi(sshtpsi,tspsi_in,info,mg,V_local,system,stencil,srg,ppg)
!$OMP parallel do collapse(5) private(ik,io,is,iz,iy,ix)
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
        do is=1,nspin
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            tspsi_out%zwf(ix,iy,iz,is,io,ik,1)=tspsi_out%zwf(ix,iy,iz,is,io,ik,1)+ &
                                               rt%zc(nn)*tspsi_in%zwf(ix,iy,iz,is,io,ik,1)
          end do
          end do
          end do
        end do
      end do
      end do
    end if
  end do

end subroutine taylor

end module taylor_sub

