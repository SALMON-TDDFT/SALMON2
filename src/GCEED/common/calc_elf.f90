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
subroutine calc_elf(lg,mg,ng,srg,info,srho,ttmp)
use structures
use salmon_parallel, only: nproc_group_global
use salmon_communication, only: comm_summation
use misc_routines, only: get_wtime
use calc_allob_sub
use scf_data
use gradient_sub
use allocate_mat_sub
use new_world_sub
use sendrecv_grid, only: init_sendrecv_grid
implicit none
type(s_rgrid),intent(in)            :: lg
type(s_rgrid),intent(in)            :: mg
type(s_rgrid),intent(in)            :: ng
type(s_sendrecv_grid),intent(in)    :: srg
type(s_orbital_parallel),intent(in) :: info
type(s_scalar),intent(in) :: srho
integer :: iob,ix,iy,iz
integer :: p_allob
integer :: ttmp

type(s_sendrecv_grid) :: srg_ob_1 ! experimental implementation

real(8) :: elftau(mg%is(1):mg%ie(1),   &
                  mg%is(2):mg%ie(2),   &
                  mg%is(3):mg%ie(3))
real(8) :: mrelftau(mg%is(1):mg%ie(1),   &
                    mg%is(2):mg%ie(2),   &
                    mg%is(3):mg%ie(3))
real(8) :: curden(mg%is(1):mg%ie(1),   &
                  mg%is(2):mg%ie(2),   &
                  mg%is(3):mg%ie(3))
real(8) :: mrcurden(mg%is(1):mg%ie(1),   &
                    mg%is(2):mg%ie(2),   &
                    mg%is(3):mg%ie(3))
real(8) :: gradpsi(3,mg%is(1):mg%ie(1),   &
                     mg%is(2):mg%ie(2),   &
                     mg%is(3):mg%ie(3))
complex(8) :: gradzpsi(3,mg%is(1):mg%ie(1),   &
                         mg%is(2):mg%ie(2),   &
                         mg%is(3):mg%ie(3))
real(8) :: gradrho(3,mg%is(1):mg%ie(1),   &
                     mg%is(2):mg%ie(2),   &
                     mg%is(3):mg%ie(3))
real(8) :: gradrho2(mg%is(1):mg%ie(1),   &
                    mg%is(2):mg%ie(2),   &
                    mg%is(3):mg%ie(3))
real(8) :: elfc(mg%is(1):mg%ie(1),   &
                mg%is(2):mg%ie(2),   &
                mg%is(3):mg%ie(3))
real(8) :: elfcuni(mg%is(1):mg%ie(1),   &
                   mg%is(2):mg%ie(2),   &
                   mg%is(3):mg%ie(3))
real(8) :: rho_half(mg%is(1):mg%ie(1),   &
                    mg%is(2):mg%ie(2),   &
                    mg%is(3):mg%ie(3))

! Create communicator for single orbital (experimental implementation):
call init_sendrecv_grid(srg_ob_1, mg, 1, info%icomm_r, srg%neig)

!$OMP parallel do private(iz,iy,ix)
do iz=mg%is(3),mg%ie(3)
do iy=mg%is(2),mg%ie(2)
do ix=mg%is(1),mg%ie(1)
  rho_half(ix,iy,iz)=srho%f(ix,iy,iz)/2.d0
end do
end do
end do
mrelftau=0.d0
mrcurden=0.d0


if(iSCFRT==1)then
  if(iperiodic==0)then

    do iob=1,iobnum
      call calc_allob(iob,info,p_allob,itotmst,mst,iobnum)
      if((ilsda==0.and.p_allob<=ifMST(1)).or.  &
         (ilsda==1.and.(p_allob<=ifMST(1).or.(p_allob>=MST(1)+1.and.p_allob<=MST(1)+ifMST(2)))))then
  
        call calc_gradient(mg, srg_ob_1, mg%is, mg%ie, psi(:,:,:,iob,1),gradpsi(:,:,:,:))
  
  !$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          mrelftau(ix,iy,iz)=mrelftau(ix,iy,iz)+abs(gradpsi(1,ix,iy,iz))**2      &
                             +abs(gradpsi(2,ix,iy,iz))**2      &
                             +abs(gradpsi(3,ix,iy,iz))**2
        end do
        end do
        end do
      end if
    end do

  else if(iperiodic==3)then

    do iob=1,iobnum
      call calc_allob(iob,info,p_allob,itotmst,mst,iobnum)
      if((ilsda==0.and.p_allob<=ifMST(1)).or.  &
         (ilsda==1.and.(p_allob<=ifMST(1).or.(p_allob>=MST(1)+1.and.p_allob<=MST(1)+ifMST(2)))))then
  
        call calc_gradient(mg, srg_ob_1, mg%is, mg%ie, zpsi(:,:,:,iob,1),gradzpsi(:,:,:,:))
  
  !$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          mrelftau(ix,iy,iz)=mrelftau(ix,iy,iz)+abs(gradzpsi(1,ix,iy,iz))**2      &
                             +abs(gradzpsi(2,ix,iy,iz))**2      &
                             +abs(gradzpsi(3,ix,iy,iz))**2
        end do
        end do
        end do
      end if
    end do

  end if

  call comm_summation(mrelftau,elftau,mg_num(1)*mg_num(2)*mg_num(3),info%icomm_o)


  call calc_gradient(mg, srg_ob_1, mg%is, mg%ie, rho_half(:,:,:),gradrho(:,:,:,:))
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    gradrho2(ix,iy,iz)=gradrho(1,ix,iy,iz)**2      &
          +gradrho(2,ix,iy,iz)**2      &
          +gradrho(3,ix,iy,iz)**2
    elfc(ix,iy,iz)=elftau(ix,iy,iz)-gradrho2(ix,iy,iz)/rho_half(ix,iy,iz)/4.d0
  end do
  end do
  end do

else 

  do iob=1,iobnum

    call calc_allob(iob,info,p_allob,itotmst,mst,iobnum)
    if((ilsda==0.and.p_allob<=ifMST(1)).or.   &
       (ilsda==1.and.(p_allob<=ifMST(1).or.(p_allob>=MST(1)+1.and.p_allob<=MST(1)+ifMST(2)))))then

      cmatbox_m=0.d0
      if(iSCFRT==1)then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          cmatbox_m(ix,iy,iz)=zpsi(ix,iy,iz,iob,1)
        end do
        end do
        end do
      else if(iSCFRT==2)then
        if(mod(ttmp,2)==1)then
!$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            cmatbox_m(ix,iy,iz)=zpsi_out(ix,iy,iz,iob,1)
          end do
          end do
          end do
        else if(mod(ttmp,2)==0)then
!$OMP parallel do private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            cmatbox_m(ix,iy,iz)=zpsi_in(ix,iy,iz,iob,1)
          end do
          end do
          end do
        end if
      end if
 
      call calc_gradient(mg, srg_ob_1, mg%is, mg%ie, cmatbox_m(:,:,:),gradzpsi(:,:,:,:))



!$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
  
        mrelftau(ix,iy,iz)=mrelftau(ix,iy,iz)+abs(gradzpsi(1,ix,iy,iz))**2      &
                           +abs(gradzpsi(2,ix,iy,iz))**2      &
                           +abs(gradzpsi(3,ix,iy,iz))**2
  
        mrcurden(ix,iy,iz)=mrcurden(ix,iy,iz)      &
             +( abs(conjg(cmatbox_m(ix,iy,iz))*gradzpsi(1,ix,iy,iz)      &
                  -cmatbox_m(ix,iy,iz)*conjg(gradzpsi(1,ix,iy,iz)))**2      &
               +abs(conjg(cmatbox_m(ix,iy,iz))*gradzpsi(2,ix,iy,iz)      &
                  -cmatbox_m(ix,iy,iz)*conjg(gradzpsi(2,ix,iy,iz)))**2      &
               +abs(conjg(cmatbox_m(ix,iy,iz))*gradzpsi(3,ix,iy,iz)      &
                  -cmatbox_m(ix,iy,iz)*conjg(gradzpsi(3,ix,iy,iz)))**2 )/2.d0
  
      end do
      end do
      end do
      



    end if
  end do



  
  call comm_summation(mrelftau,elftau,mg_num(1)*mg_num(2)*mg_num(3),info%icomm_o)
  call comm_summation(mrcurden,curden,mg_num(1)*mg_num(2)*mg_num(3),info%icomm_o)

  


  
  call calc_gradient(mg, srg_ob_1, mg%is, mg%ie, rho_half(:,:,:),gradrho(:,:,:,:))



  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    gradrho2(ix,iy,iz)=gradrho(1,ix,iy,iz)**2      &
          +gradrho(2,ix,iy,iz)**2      &
          +gradrho(3,ix,iy,iz)**2
    elfc(ix,iy,iz)=elftau(ix,iy,iz)-gradrho2(ix,iy,iz)/rho_half(ix,iy,iz)/4.d0  &
                                   -curden(ix,iy,iz)/rho_half(ix,iy,iz)
  end do
  end do
  end do




end if

! matbox_l stores ELF
matbox_l=0.d0
do iz=ng%is(3),ng%ie(3)
do iy=ng%is(2),ng%ie(2)
do ix=ng%is(1),ng%ie(1)
  elfcuni(ix,iy,iz)=3.d0/5.d0*(6.d0*Pi**2)**(2.d0/3.d0)      &
            *rho_half(ix,iy,iz)**(5.d0/3.d0)
  matbox_l(ix,iy,iz)=1.d0/(1.d0+elfc(ix,iy,iz)**2/elfcuni(ix,iy,iz)**2)
end do
end do
end do

call comm_summation(matbox_l,elf,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)


call dealloc_cache(srg_ob_1)

end subroutine calc_elf
