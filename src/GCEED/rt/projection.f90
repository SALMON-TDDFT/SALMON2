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
subroutine projection(itt,mg,system,info,tpsi,tpsi0)
use structures
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_global, only: dt,iwrite_projnum,iwrite_projection_k,iwrite_projection_ob,num_projection
use pack_unpack, only: copy_data
use scf_data, only: iSCFRT
implicit none
integer                 ,intent(in) :: itt
type(s_rgrid)           ,intent(in) :: mg
type(s_dft_system)      ,intent(in) :: system
type(s_orbital_parallel),intent(in) :: info
type(s_orbital)         ,intent(in) :: tpsi,tpsi0
!
integer :: nspin,no,nk,ik_s,ik_e,io_s,io_e,is(3),ie(3)
integer :: ix,iy,iz,io1,io2,io,ik,ispin
complex(8),dimension(system%no,system%no,system%nspin,system%nk) :: mat1,mat2
complex(8) :: wf_io1(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3))
real(8) :: coef(system%no,system%nk,system%nspin)
complex(8) :: cbox
character(100) :: projOutFile
character(20) :: fileNumber

if(iSCFRT==2)then
  if(iwrite_projnum==1)then
    write(fileNumber, '(i8)') itt
    projOutFile = trim("proj.")//adjustl(fileNumber)
    open(61,file=projOutFile)
  end if
end if

if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ projection"

nspin = system%nspin
no = system%no
nk = system%nk
is = mg%is
ie = mg%ie
ik_s = info%ik_s
ik_e = info%ik_e
io_s = info%io_s
io_e = info%io_e

! copied from subspace_diagonalization.f90
mat1 = 0d0
if(info%if_divide_orbit) then
  do ik=ik_s,ik_e
  do ispin = 1, nspin
    do io1 = 1, no ! future work: no --> no0 (# of GS orbitals)
      if (io_s<= io1 .and. io1 <= io_e) then
        call copy_data(tpsi0%zwf(:, :, :, ispin, io1, ik, 1),wf_io1)
      end if
      call comm_bcast(wf_io1, info%icomm_o, info%irank_io(io1))
      do io2 = 1, no
        if (io_s<= io2 .and. io2 <= io_e) then
          cbox = 0d0
          !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+:cbox)
          do iz=is(3),ie(3)
          do iy=is(2),ie(2)
          do ix=is(1),ie(1)
            cbox = cbox + conjg(wf_io1(ix,iy,iz)) * tpsi%zwf(ix,iy,iz,ispin,io2,ik,1)
          end do
          end do
          end do
          mat1(io1,io2,ispin,ik) = cbox * system%hvol
        end if
      end do
    end do !io1
  end do !ispin
  end do
else
  !$omp parallel do private(ik,io1,io2,ispin,cbox,iz,iy,ix) collapse(4)
  do ik=ik_s,ik_e
  do ispin=1,nspin
  do io1=io_s,io_e
  do io2=io_s,io_e
    cbox = 0d0
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      cbox = cbox + conjg(tpsi0%zwf(ix,iy,iz,ispin,io1,ik,1)) * tpsi%zwf(ix,iy,iz,ispin,io2,ik,1)
    end do
    end do
    end do
    mat1(io1,io2,ispin,ik) = cbox * system%hvol
  end do
  end do
  end do
  end do
end if

call comm_summation(mat1,mat2,no**2*nspin*nk,info%icomm_rko)

coef=0.d0
do ispin=1,nspin
do ik=1,nk
do io1=1,no
  do io2=1,no
    coef(io1,ik,ispin) = coef(io1,ik,ispin)+abs(mat2(io2,io1,ispin,ik)*system%Hvol)**2
  end do
end do
end do
end do

if(comm_is_root(nproc_id_global))then
  write(41,'(200f14.8)') dble(itt)*dt*2.41888d-2, &
  & ((coef(iwrite_projection_ob(io),iwrite_projection_k(io),ispin),io=1,num_projection),ispin=1,nspin),  &
    sum(coef(1:no,:,1)),sum(coef(1:no,:,1)) ! no ---> no0
end if
if(mod(itt,100)==0)then
  if(comm_is_root(nproc_id_global))then
    do ik=1,nk
    do io=1,no ! no --> no0
      write(*,'(a12,3i6,f16.8)') "projection",io,ik,(coef(io,ik,ispin),ispin=1,nspin)
    end do
    end do
  end if
end if

if(iSCFRT==2)then
  if(iwrite_projnum==1)then
    close(61)
  end if
end if

end subroutine projection
