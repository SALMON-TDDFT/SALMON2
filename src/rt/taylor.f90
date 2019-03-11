!
!  Copyright 2017 SALMON developers
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

subroutine taylor(mg,nspin,info,itotmst,mst,lg_sta,lg_end,ilsda,stencil,tspsi_in,tspsi_out,sshtpsi,   &
                  ppg,vlocal,vbox,num_kpoints_rd,k_rd,zc,ihpsieff,rocc,wtk,iparaway_ob)
  use inputoutput, only: iperiodic,ispin,natom,n_hamil
  use structures, only: s_rgrid,s_wf_info,s_wavefunction,s_stencil,s_scalar,s_pp_grid
  use hpsi_sub
  use calc_allob_sub
  use sendrecv_grid, only: s_sendrecv_grid
  implicit none
  integer,parameter     :: nd=4 
  type(s_rgrid),intent(in) :: mg
  integer,intent(in)    :: nspin
  type(s_wf_info),intent(in) :: info
  integer,intent(in) :: itotmst
  integer,intent(in) :: mst(2)
  integer,intent(in) :: lg_sta(3)
  integer,intent(in) :: lg_end(3)
  integer,intent(in)    :: ilsda
  type(s_stencil),intent(inout) :: stencil
  type(s_sendrecv_grid),intent(in) :: srg
  type(s_wavefunction),intent(inout) :: tspsi_in
  type(s_wavefunction),intent(inout) :: tspsi_out
  type(s_wavefunction),intent(inout) :: sshtpsi
  type(s_pp_grid),intent(inout) :: ppg
  real(8),intent(in)    :: vlocal(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),ispin+1)
  real(8),intent(in)    :: vbox(lg_sta(1)-nd:lg_end(1)+nd,  &
                                lg_sta(2)-nd:lg_end(2)+nd,  &
                                lg_sta(3)-nd:lg_end(3)+nd)
  integer,intent(in)    :: num_kpoints_rd
  real(8),intent(in)    :: k_rd(3,num_kpoints_rd)
  complex(8),intent(in) :: zc(n_hamil)
  integer,intent(in)    :: ihpsieff
  real(8),intent(in)    :: rocc(itotmst,num_kpoints_rd)
  real(8),intent(in)    :: wtk(num_kpoints_rd)
  integer,intent(in)    :: iparaway_ob
  type(s_scalar),allocatable :: v(:)
  integer :: nn,ix,iy,iz
  integer :: ik,io,io_allob
  complex(8) :: ekr(ppg%nps,natom)
  integer :: a,iatom
  integer :: ilma,j
  real(8) :: x,y,z
  complex(8),parameter :: zi=(0.d0,1.d0)
  integer :: is
  
  allocate(v(nspin))
  do is=1,nspin
    allocate(v(is)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  end do

  if(iperiodic==0.and.ihpsieff==1)then
!$OMP parallel do collapse(3) private(is,iz,iy,ix)
    do is=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        v(is)%f(ix,iy,iz)=vlocal(ix,iy,iz,is)+vbox(ix,iy,iz)
      end do
      end do
      end do
    end do
  else
!$OMP parallel do collapse(3) private(is,iz,iy,ix)
    do is=1,nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        v(is)%f(ix,iy,iz)=vlocal(ix,iy,iz,is)
      end do
      end do
      end do
    end do
  end if

  do ik=info%ik_s,info%ik_e
    if(iperiodic==3)then
      do j=1,3
        stencil%kAc(ik,j) = k_rd(j,ik)
      end do
    end if
  end do
  if(iperiodic==3) then
    call update_kvector_nonlocalpt(ppg,stencil%kAc,info%ik_s,info%ik_e)
  end if

!$OMP parallel do collapse(5) private(ik,io,is,iz,iy,ix)
  do ik=info%ik_s,info%ik_e
  do io=info%io_s,info%io_e
    do is=1,nspin
      do iz=mg%is_array(3),mg%ie_array(3)
      do iy=mg%is_array(2),mg%ie_array(2)
      do ix=mg%is_array(1),mg%ie_array(1)
        tspsi_out%zwf(ix,iy,iz,is,io,ik,1)=tspsi_in%zwf(ix,iy,iz,is,io,ik,1)
      end do
      end do
      end do
    end do
  end do
  end do

  do nn=1,n_hamil
    if(mod(nn,2)==1)then
      call hpsi(tspsi_in,sshtpsi,info,mg,v,nspin,stencil,srg,ppg)
!$OMP parallel do collapse(5) private(ik,io,is,iz,iy,ix)
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
        do is=1,nspin
          do iz=mg%is_array(3),mg%ie_array(3)
          do iy=mg%is_array(2),mg%ie_array(2)
          do ix=mg%is_array(1),mg%ie_array(1)
            tspsi_out%zwf(ix,iy,iz,is,io,ik,1)=tspsi_out%zwf(ix,iy,iz,is,io,ik,1)+ &
                                                 zc(nn)*sshtpsi%zwf(ix,iy,iz,is,io,ik,1)
          end do
          end do
          end do
        end do
      end do
      end do
    else
      call hpsi(sshtpsi,tspsi_in,info,mg,v,nspin,stencil,srg,ppg)
!$OMP parallel do collapse(5) private(ik,io,is,iz,iy,ix)
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
        do is=1,nspin
          do iz=mg%is_array(3),mg%ie_array(3)
          do iy=mg%is_array(2),mg%ie_array(2)
          do ix=mg%is_array(1),mg%ie_array(1)
            tspsi_out%zwf(ix,iy,iz,is,io,ik,1)=tspsi_out%zwf(ix,iy,iz,is,io,ik,1)+ &
                                               zc(nn)*tspsi_in%zwf(ix,iy,iz,is,io,ik,1)
          end do
          end do
          end do
        end do
      end do
      end do
    end if
  end do

  do is=1,nspin
    deallocate(v(is)%f)
  end do
  deallocate(v)
  if(allocated(ppg%zproj)) deallocate(ppg%zproj)

end subroutine taylor

end module taylor_sub

