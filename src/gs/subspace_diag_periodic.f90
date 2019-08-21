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
module subspace_diag_periodic_sub
  implicit none

contains

subroutine subspace_diag_periodic(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
  use structures
  use salmon_communication, only: comm_bcast, comm_summation
  use timer
  use hpsi_sub
  use eigen_subdiag_periodic_sub
  use sendrecv_grid, only: s_sendrecv_grid
  implicit none
  type(s_rgrid)           ,intent(in) :: mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_stencil),intent(in) :: stencil
  type(s_pp_grid),intent(in) :: ppg
  type(s_scalar) ,intent(in) :: vlocal(system%nspin)
  type(s_orbital)            :: spsi,shpsi
  type(s_sendrecv_grid)      :: srg
  !
  integer :: nspin,no,nk,ik,io,io1,io2,ispin,ik_s,ik_e,io_s,io_e,is(3),ie(3),ix,iy,iz,ierr
  real(8)   ,dimension(system%nspin,system%no,system%nk) :: rbox1,rbox2
  complex(8),dimension(system%no,system%no,system%nspin,system%nk) :: mat1,mat2,evec
  complex(8) :: cbox
  complex(8),parameter :: zi=(0.d0,1.d0)
  
  call timer_begin(LOG_DIAG_TOTAL)

  if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ subspace_diag"
  if(info%if_divide_orbit) stop "error: nproc_ob/=1 @ subspace_diag"

  nspin = system%nspin
  no = system%no
  nk = system%nk
  is = mg%is
  ie = mg%ie
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e

  call hpsi(spsi,shpsi,info,mg,vlocal,system,stencil,srg,ppg)

  mat1 = 0d0

  !$omp parallel do private(ik,io1,io2,ispin,cbox,iz,iy,ix) collapse(4)
  do ik=ik_s,ik_e
  do ispin=1,nspin
  do io1=io_s,io_e
  do io2=io_s,io_e
    cbox = 0d0
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      cbox = cbox + conjg(spsi%zwf(ix,iy,iz,ispin,io1,ik,1)) * shpsi%zwf(ix,iy,iz,ispin,io2,ik,1)
    end do
    end do
    end do
    mat1(io1,io2,ispin,ik) = cbox * system%hvol
  end do
  end do
  end do
  end do

  call timer_begin(LOG_DIAG_ALLREDUCE)
  call comm_summation(mat1,mat2,no**2*nspin*nk,info%icomm_rko)
  call timer_end(LOG_DIAG_ALLREDUCE)
    
  do ik=ik_s,ik_e
  do ispin=1,nspin
    call timer_begin(LOG_DIAG_EIGEN)
    call eigen_subdiag_periodic(mat2(:,:,ispin,ik),evec(:,:,ispin,ik),no,ierr)
    call timer_end(LOG_DIAG_EIGEN)
  end do
  end do

!$omp workshare
  shpsi%zwf = 0d0
!$omp end workshare

  call timer_begin(LOG_DIAG_UPDATE)
  !$omp parallel do private(ik,io1,io2,ispin) collapse(3)
  do ik=ik_s,ik_e
  do ispin=1,nspin
  do io1=io_s,io_e
  do io2=io_s,io_e
    shpsi%zwf(:,:,:,ispin,io1,ik,1) = shpsi%zwf(:,:,:,ispin,io1,ik,1) + evec(io2,io1,ispin,ik) * spsi%zwf(:,:,:,ispin,io2,ik,1)
  end do
  end do
  end do
  end do

! normalization
  rbox1 = 0d0
  !$omp parallel do private(ik,io,ispin,iz,iy,ix) collapse(2)
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      rbox1(ispin,io,ik) = rbox1(ispin,io,ik) + abs(shpsi%zwf(ix,iy,iz,ispin,io,ik,1))**2
    end do
    end do
    end do
  end do
  end do
  end do
  call comm_summation(rbox1,rbox2,nspin*no*nk,info%icomm_rko)

  !$omp parallel do private(ik,io,ispin,iz,iy,ix) collapse(2)
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      spsi%zwf(ix,iy,iz,ispin,io,ik,1) = shpsi%zwf(ix,iy,iz,ispin,io,ik,1) / sqrt(rbox2(ispin,io,ik)*system%hvol)
    end do
    end do
    end do
  end do
  end do
  end do
  call timer_end(LOG_DIAG_UPDATE)

  call timer_end(LOG_DIAG_TOTAL)

end subroutine subspace_diag_periodic

end module subspace_diag_periodic_sub
