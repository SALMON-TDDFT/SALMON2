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

#include "config.h"

module subspace_diagonalization

  use subspace_diagonalization_so, only: ssdg_periodic_so, SPIN_ORBIT_ON

  implicit none

contains

subroutine ssdg_rwf(mg,system,info,pinfo,stencil,spsi,shpsi,ppg,vlocal,srg)
  use structures
  use communication, only: comm_summation,comm_bcast
  use timer
  use hamiltonian, only: hpsi
  use eigen_subdiag_sub
  use sendrecv_grid, only: s_sendrecv_grid
  use pack_unpack, only: copy_data
  implicit none
  type(s_rgrid)           ,intent(in) :: mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_process_info)    ,intent(in) :: pinfo
  type(s_stencil),intent(in) :: stencil
  type(s_pp_grid),intent(in) :: ppg
  type(s_scalar) ,intent(in) :: vlocal(system%nspin)
  type(s_orbital)            :: spsi,shpsi
  type(s_sendrecv_grid)      :: srg
  !
  integer :: nspin,no,io,io1,io2,ispin,io_s,io_e,is(3),ie(3),ix,iy,iz,ierr
  real(8)   ,dimension(system%nspin,system%no) :: rbox1,rbox2
  real(8),dimension(system%no,system%no,system%nspin) :: mat1,mat2,evec
  real(8) :: cbox
  real(8) :: wf_io1(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3))
  real(8) :: wf_io2(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3))

  if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ subspace_diag"

  call timer_begin(LOG_SSDG_ISOLATED_CALC)
  nspin = system%nspin
  no = system%no
  is = mg%is
  ie = mg%ie
  io_s = info%io_s
  io_e = info%io_e
  call timer_end(LOG_SSDG_ISOLATED_CALC)

  call timer_begin(LOG_SSDG_ISOLATED_HPSI)
  call hpsi(spsi,shpsi,info,mg,vlocal,system,stencil,srg,ppg)
  call timer_end(LOG_SSDG_ISOLATED_HPSI)

  call timer_begin(LOG_SSDG_ISOLATED_CALC)
  mat1 = 0d0
  if(info%if_divide_orbit) then
    do ispin = 1, nspin
      do io1 = 1, no
        if (io_s<= io1 .and. io1 <= io_e) then
          call copy_data(spsi%rwf(:, :, :, ispin, io1, 1, 1),wf_io1)
        end if
        call timer_end(LOG_SSDG_ISOLATED_CALC)

        call timer_begin(LOG_SSDG_ISOLATED_COMM_COLL)
        call comm_bcast(wf_io1, info%icomm_o, info%irank_io(io1))
        call timer_end(LOG_SSDG_ISOLATED_COMM_COLL)

        call timer_begin(LOG_SSDG_ISOLATED_CALC)
        do io2 = 1, no
          if (io_s<= io2 .and. io2 <= io_e) then
            cbox = 0d0
            !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+:cbox)
            do iz=is(3),ie(3)
            do iy=is(2),ie(2)
            do ix=is(1),ie(1)
              cbox = cbox + wf_io1(ix,iy,iz) * shpsi%rwf(ix,iy,iz,ispin,io2,1,1)
            end do
            end do
            end do
            mat1(io1,io2,ispin) = cbox * system%hvol
          end if
        end do
      end do !io1
    end do !ispin
  else
    !$omp parallel do private(io1,io2,ispin,cbox,iz,iy,ix) collapse(3)
    do ispin=1,nspin
    do io1=io_s,io_e
    do io2=io_s,io_e
      cbox = 0d0
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        cbox = cbox + spsi%rwf(ix,iy,iz,ispin,io1,1,1) * shpsi%rwf(ix,iy,iz,ispin,io2,1,1)
      end do
      end do
      end do
      mat1(io1,io2,ispin) = cbox * system%hvol
    end do
    end do
    end do
  end if
  call timer_end(LOG_SSDG_ISOLATED_CALC)

  call timer_begin(LOG_SSDG_ISOLATED_COMM_COLL)
  call comm_summation(mat1,mat2,no**2*nspin,info%icomm_rko)
  call timer_end(LOG_SSDG_ISOLATED_COMM_COLL)

  call timer_begin(LOG_SSDG_ISOLATED_CALC)
  do ispin=1,nspin
    call eigen_subdiag(mat2(:,:,ispin),evec(:,:,ispin),no,ierr)
  end do

!$omp workshare
  shpsi%rwf = 0d0
!$omp end workshare

  if(info%if_divide_orbit) then
    do ispin=1,nspin
    do io2 = 1, no
      if (io_s<= io2 .and. io2 <= io_e) then
        call copy_data(spsi%rwf(:, :, :, ispin, io2, 1, 1),wf_io2)
      end if
      call timer_end(LOG_SSDG_ISOLATED_CALC)

      call timer_begin(LOG_SSDG_ISOLATED_COMM_COLL)
      call comm_bcast(wf_io2, info%icomm_o, info%irank_io(io2))
      call timer_end(LOG_SSDG_ISOLATED_COMM_COLL)

      call timer_begin(LOG_SSDG_ISOLATED_CALC)
      !$omp parallel do private(io1,iz,iy,ix) collapse(3)
      do io1=io_s,io_e
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
          shpsi%rwf(ix,iy,iz,ispin,io1,1,1) = shpsi%rwf(ix,iy,iz,ispin,io1,1,1) + evec(io2,io1,ispin) * wf_io2(ix,iy,iz)
        end do
        end do
        end do
      end do
    end do
    end do
  else
  !$omp parallel do private(io1,io2,ispin) collapse(2)
    do ispin=1,nspin
    do io1=io_s,io_e
    do io2=io_s,io_e
      shpsi%rwf(:,:,:,ispin,io1,1,1) = shpsi%rwf(:,:,:,ispin,io1,1,1) + evec(io2,io1,ispin) * spsi%rwf(:,:,:,ispin,io2,1,1)
    end do
    end do
    end do
  end if

! normalization
  rbox1 = 0d0
  !$omp parallel do private(io,ispin,iz,iy,ix)
  do io=io_s,io_e
  do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      rbox1(ispin,io) = rbox1(ispin,io) + abs(shpsi%rwf(ix,iy,iz,ispin,io,1,1))**2
    end do
    end do
    end do
  end do
  end do
  call timer_end(LOG_SSDG_ISOLATED_CALC)

  call timer_begin(LOG_SSDG_ISOLATED_COMM_COLL)
  call comm_summation(rbox1,rbox2,nspin*no,info%icomm_rko)
  call timer_end(LOG_SSDG_ISOLATED_COMM_COLL)

  call timer_begin(LOG_SSDG_ISOLATED_CALC)
  !$omp parallel do private(io,ispin,iz,iy,ix)
  do io=io_s,io_e
  do ispin=1,nspin
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      spsi%rwf(ix,iy,iz,ispin,io,1,1) = shpsi%rwf(ix,iy,iz,ispin,io,1,1) / sqrt(rbox2(ispin,io)*system%hvol)
    end do
    end do
    end do
  end do
  end do
  call timer_end(LOG_SSDG_ISOLATED_CALC)

end subroutine ssdg_rwf

!===================================================================================================================================

subroutine ssdg_zwf(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg,pinfo)
  use salmon_global, only: yn_gbp
  use structures
  use communication, only: comm_summation,comm_bcast
  use timer
  use hamiltonian, only: hpsi
  use eigen_subdiag_sub
  use sendrecv_grid, only: s_sendrecv_grid
  use pack_unpack, only: copy_data
  implicit none
  type(s_rgrid)           ,intent(in) :: mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_stencil),intent(in) :: stencil
  type(s_pp_grid),intent(in) :: ppg
  type(s_scalar) ,intent(in) :: vlocal(system%nspin)
  type(s_orbital)            :: spsi,shpsi
  type(s_sendrecv_grid)      :: srg
  type(s_process_info),intent(in) :: pinfo

  if(yn_gbp=='c') then
    call ssdg_periodic_cblas(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg,pinfo)
  elseif (yn_gbp=='r') then
    call ssdg_periodic_rblas(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg,pinfo)
  else
    call ssdg_periodic_org(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
  end if

  return
end subroutine

  !
  !===================================================================================================================================

  subroutine ssdg_periodic_cblas(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg,pinfo)
    use structures
    use communication, only: comm_summation,comm_bcast
    use timer
    use hamiltonian, only: hpsi
    use eigen_subdiag_sub
    use sendrecv_grid, only: s_sendrecv_grid
    use pack_unpack, only: copy_data
    implicit none
    type(s_rgrid)           ,intent(in) :: mg
    type(s_dft_system)      ,intent(in) :: system
    type(s_orbital_parallel),intent(in) :: info
    type(s_stencil),intent(in) :: stencil
    type(s_pp_grid),intent(in) :: ppg
    type(s_scalar) ,intent(in) :: vlocal(system%nspin)
    type(s_orbital)            :: spsi,shpsi
    type(s_sendrecv_grid)      :: srg
    type(s_process_info),intent(in) :: pinfo

  complex(8),parameter :: zero = 0d0, one = 1d0
  integer :: im,ispin,ik,io,jo,io1,io2,nsize_rg,ierr,m
  complex(8),dimension(system%no,system%no) :: hmat,hmat_tmp,evec
  complex(8) :: wf1_block(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), info%numo)
  complex(8) :: wf2_block(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), info%numo)
  complex(8) :: wf_block_send(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), info%numo_max)
  complex(8) :: hmat_block(info%numo_max, info%numo)
  complex(8) :: hmat_block_tmp(info%numo_max, info%numo)


  if ( SPIN_ORBIT_ON ) then
    call ssdg_periodic_so(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
    return
  end if

  call timer_begin(LOG_SSDG_PERIODIC_HPSI)
  call hpsi(spsi,shpsi,info,mg,vlocal,system,stencil,srg,ppg)
  call timer_end(LOG_SSDG_PERIODIC_HPSI)

  nsize_rg = (mg%ie(1)-mg%is(1)+1)*(mg%ie(2)-mg%is(2)+1)*(mg%ie(3)-mg%is(3)+1)

  do im = info%im_s, info%im_e
  do ik = info%ik_s, info%ik_e
  do ispin = 1, system%nspin

    call timer_begin(LOG_SSDG_PERIODIC_CALC)
    ! Copy wave function
    do io = info%io_s, info%io_e
      jo = io - info%io_s + 1
      call copy_data( &
        & spsi%zwf(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), ispin, io, ik, im), &
        & wf1_block(:, :, :, jo))
      call copy_data( &
        & shpsi%zwf(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), ispin, io, ik, im), &
        & wf2_block(:, :, :, jo))
    end do

    wf_block_send = 0d0
    hmat_tmp = 0d0
    do m = 0, pinfo%nporbital - 1

      if(m == info%id_o ) then
        call copy_data( wf1_block(:,:,:,:), wf_block_send(:,:,:,1:info%numo) )
      end if
      call timer_end(LOG_SSDG_PERIODIC_CALC)

      call timer_begin(LOG_SSDG_PERIODIC_COMM_COLL)
      if (info%if_divide_orbit) then
        call comm_bcast( wf_block_send(:,:,:,1:info%numo_all(m)), info%icomm_o, info%irank_io(info%io_s_all(m)))
      end if
      call timer_end(LOG_SSDG_PERIODIC_COMM_COLL)

      call timer_begin(LOG_SSDG_PERIODIC_CALC)
      hmat_block_tmp = 0d0
      call zgemm('C', 'N', info%numo_all(m), info%numo, nsize_rg,  &
        &   dcmplx(system%hvol), wf_block_send(:,:,:,1:info%numo_all(m)), nsize_rg,  &
        &                        wf2_block(:,:,:,1), nsize_rg,  &
        &                  zero, hmat_block_tmp(1:info%numo_all(m),1:info%numo), info%numo_all(m) )
      call timer_end(LOG_SSDG_PERIODIC_CALC)

      call timer_begin(LOG_SSDG_PERIODIC_COMM_COLL)
      if (info%if_divide_rspace) then
        call comm_summation(hmat_block_tmp, hmat_block, info%numo_max*info%numo, info%icomm_r)
      else
         hmat_block = hmat_block_tmp
      endif
      call timer_end(LOG_SSDG_PERIODIC_COMM_COLL)

      call timer_begin(LOG_SSDG_PERIODIC_CALC)
      do io1 = info%io_s_all(m), info%io_e_all(m)
      do io2 = info%io_s, info%io_e
        hmat_tmp(io1,io2) = hmat_block(io1-info%io_s_all(m)+1, io2-info%io_s+1)
      end do
      end do

    end do ! m
    call timer_end(LOG_SSDG_PERIODIC_CALC)

    call timer_begin(LOG_SSDG_PERIODIC_COMM_COLL)
    if (info%if_divide_orbit) then
      call comm_summation(hmat_tmp, hmat, system%no*system%no, info%icomm_o)
    else
      hmat = hmat_tmp
    end if
    call timer_end(LOG_SSDG_PERIODIC_COMM_COLL)


    call timer_begin(LOG_SSDG_PERIODIC_EIGEN)
    call eigen_subdiag_periodic(hmat, evec, system%no, ierr)
    call timer_end(LOG_SSDG_PERIODIC_EIGEN)


    call timer_begin(LOG_SSDG_PERIODIC_CALC)
    !$omp workshare
    wf2_block = 0d0
    !$omp end workshare

    wf_block_send = 0d0

    do m = 0, pinfo%nporbital - 1

      if(m == info%id_o ) then
        call copy_data( wf1_block(:,:,:,:), wf_block_send(:,:,:,1:info%numo) )
      end if
      call timer_end(LOG_SSDG_PERIODIC_CALC)

      call timer_begin(LOG_SSDG_PERIODIC_COMM_COLL)
      if (info%if_divide_orbit) then
        call comm_bcast( wf_block_send(:,:,:,1:info%numo_all(m)), info%icomm_o, info%irank_io(info%io_s_all(m)))
      end if
      call timer_end(LOG_SSDG_PERIODIC_COMM_COLL)

      call timer_begin(LOG_SSDG_PERIODIC_CALC)
      call zgemm('N', 'N', nsize_rg, info%numo, info%numo_all(m),  &
        &                one, wf_block_send(:,:,:,1:info%numo_all(m)), nsize_rg,  &
        &                     evec(info%io_s_all(m):info%io_e_all(m), info%io_s:info%io_e), info%numo_all(m),  &
        &                one, wf2_block(:,:,:,1), nsize_rg )

    end do ! m

    ! Copy wave function
    do io = info%io_s, info%io_e
      jo = io - info%io_s + 1
      call copy_data( &
        & wf2_block(:, :, :, jo), &
        & spsi%zwf(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), ispin, io, ik, im) )
    end do
    call timer_end(LOG_SSDG_PERIODIC_CALC)

  enddo ! ispin
  enddo ! ik
  enddo ! im

return

end subroutine ssdg_periodic_cblas
!===================================================================================================================================

!===================================================================================================================================

subroutine ssdg_periodic_rblas(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg,pinfo)
  use structures
  use communication, only: comm_summation,comm_bcast
  use timer
  use hamiltonian, only: hpsi
  use eigen_subdiag_sub
  use sendrecv_grid, only: s_sendrecv_grid
  use pack_unpack, only: copy_data
  use salmon_global, only: yn_pdsyev
  implicit none
  type(s_rgrid)           ,intent(in) :: mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_stencil),intent(in) :: stencil
  type(s_pp_grid),intent(in) :: ppg
  type(s_scalar) ,intent(in) :: vlocal(system%nspin)
  type(s_orbital)            :: spsi,shpsi
  type(s_sendrecv_grid)      :: srg
  type(s_process_info),intent(in) :: pinfo

real(8),parameter :: zero = 0d0, one = 1d0
integer :: im,ispin,ik,io,jo,io1,io2,nsize_rg,m !,ierr
real(8),dimension(system%no,system%no) :: hmat,hmat_tmp,evec
real(8) :: wf1_block(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), info%numo)
real(8) :: wf2_block(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), info%numo)
real(8) :: wf_block_send(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), info%numo_max)
real(8) :: hmat_block(info%numo_max, info%numo)
real(8) :: hmat_block_tmp(info%numo_max, info%numo)
!complex(8),dimension(system%no,system%no) :: zhmat, zevec
real(8) :: eval(system%no)

if ( SPIN_ORBIT_ON ) then
  call ssdg_periodic_so(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
  return
end if

call timer_begin(LOG_SSDG_PERIODIC_HPSI)
call hpsi(spsi,shpsi,info,mg,vlocal,system,stencil,srg,ppg)
call timer_end(LOG_SSDG_PERIODIC_HPSI)

nsize_rg = (mg%ie(1)-mg%is(1)+1)*(mg%ie(2)-mg%is(2)+1)*(mg%ie(3)-mg%is(3)+1)

do im = info%im_s, info%im_e
do ik = info%ik_s, info%ik_e
do ispin = 1, system%nspin

  call timer_begin(LOG_SSDG_PERIODIC_CALC)
  ! Copy wave function
  do io = info%io_s, info%io_e
    jo = io - info%io_s + 1
    call copy_data( &
      &  spsi%rwf(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), ispin, io, ik, im) , &
      & wf1_block(:, :, :, jo))
    call copy_data( &
      & shpsi%rwf(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), ispin, io, ik, im) , &
      & wf2_block(:, :, :, jo))
  end do

  wf_block_send = 0d0 !to fix NaN for large calc, but this may not real solution
  hmat_tmp = 0d0
  do m = 0, pinfo%nporbital - 1

    if(m == info%id_o ) then
      call copy_data( wf1_block(:,:,:,:), wf_block_send(:,:,:,1:info%numo) )
    end if
    call timer_end(LOG_SSDG_PERIODIC_CALC)

    call timer_begin(LOG_SSDG_PERIODIC_COMM_COLL)
    if (info%if_divide_orbit) then
      call comm_bcast( wf_block_send(:,:,:,1:info%numo_all(m)), info%icomm_o, info%irank_io(info%io_s_all(m)))
    end if
    call timer_end(LOG_SSDG_PERIODIC_COMM_COLL)

    call timer_begin(LOG_SSDG_PERIODIC_CALC)
    hmat_block_tmp = 0d0
    call dgemm('T', 'N', info%numo_all(m), info%numo, nsize_rg,  &
      &           system%hvol, wf_block_send(:,:,:,1:info%numo_all(m)), nsize_rg,  &
      &                        wf2_block(:,:,:,1), nsize_rg,  &
      &                  zero, hmat_block_tmp(1:info%numo_all(m),1:info%numo), info%numo_all(m) )
    call timer_end(LOG_SSDG_PERIODIC_CALC)

    call timer_begin(LOG_SSDG_PERIODIC_COMM_COLL)
    if (info%if_divide_rspace) then
      call comm_summation(hmat_block_tmp, hmat_block, info%numo_max*info%numo, info%icomm_r)
    else
      hmat_block = hmat_block_tmp
    endif
    call timer_end(LOG_SSDG_PERIODIC_COMM_COLL)

    call timer_begin(LOG_SSDG_PERIODIC_CALC)
    do io1 = info%io_s_all(m), info%io_e_all(m)
    do io2 = info%io_s, info%io_e
      hmat_tmp(io1,io2) = hmat_block(io1-info%io_s_all(m)+1, io2-info%io_s+1)
    end do
    end do

  end do ! m
  call timer_end(LOG_SSDG_PERIODIC_CALC)

  call timer_begin(LOG_SSDG_PERIODIC_COMM_COLL)
  if (info%if_divide_orbit) then
    call comm_summation(hmat_tmp, hmat, system%no*system%no, info%icomm_o)
  else
    hmat = hmat_tmp
  end if
  call timer_end(LOG_SSDG_PERIODIC_COMM_COLL)


  call timer_begin(LOG_SSDG_PERIODIC_EIGEN)
!  zhmat = hmat

!  call eigen_subdiag_periodic(zhmat, zevec, system%no, ierr)
  if(yn_pdsyev=='y') then
#ifdef USE_SCALAPACK
     call eigen_real8_pdsyev(pinfo, info, hmat, eval, evec)
#else
     stop "ScaLAPACK does not enabled, please check your build configuration."
#endif
  else
     call eigen_real8_dsyev(hmat, eval, evec)
  endif

!  evec = real(zevec)
  call timer_end(LOG_SSDG_PERIODIC_EIGEN)

  call timer_begin(LOG_SSDG_PERIODIC_CALC)
  !$omp workshare
  wf2_block = 0d0
  !$omp end workshare

  wf_block_send = 0d0 !to fix NaN for large calc, but this may not real solution

  do m = 0, pinfo%nporbital - 1

    if(m == info%id_o ) then
      call copy_data( wf1_block(:,:,:,:), wf_block_send(:,:,:,1:info%numo) )
    end if
    call timer_end(LOG_SSDG_PERIODIC_CALC)

    call timer_begin(LOG_SSDG_PERIODIC_COMM_COLL)
    if (info%if_divide_orbit) then
      call comm_bcast( wf_block_send(:,:,:,1:info%numo_all(m)), info%icomm_o, info%irank_io(info%io_s_all(m)))
    end if
    call timer_end(LOG_SSDG_PERIODIC_COMM_COLL)

    call timer_begin(LOG_SSDG_PERIODIC_CALC)
    call dgemm('N', 'N', nsize_rg, info%numo, info%numo_all(m),  &
      &                one, wf_block_send(:,:,:,1:info%numo_all(m)), nsize_rg,  &
      &                     evec(info%io_s_all(m):info%io_e_all(m), info%io_s:info%io_e), info%numo_all(m),  &
      &                one, wf2_block(:,:,:,1), nsize_rg )
  end do ! m

  ! Copy wave function
  do io = info%io_s, info%io_e
    jo = io - info%io_s + 1
    call copy_data( &
      &  wf2_block(:, :, :, jo) , &
      & spsi%rwf(mg%is(1):mg%ie(1), mg%is(2):mg%ie(2), mg%is(3):mg%ie(3), ispin, io, ik, im) )
  end do
  call timer_end(LOG_SSDG_PERIODIC_CALC)

enddo ! ispin
enddo ! ik
enddo ! im

return

contains
  subroutine eigen_real8_dsyev(h,e,v)
    implicit none
    real(8), intent(in) :: h(:,:)
    real(8), intent(out) :: e(:)
    real(8), intent(out) :: v(:,:)
    real(8), allocatable :: work(:)
    integer :: n, lwork, info

    n = ubound(h,1)
    lwork = 3*n-1
    allocate(work(lwork))
    v=h
    call dsyev('V', 'U', n, v, n, e, work, lwork, info)
    deallocate(work)
    return
  end subroutine eigen_real8_dsyev

end subroutine ssdg_periodic_rblas
!===================================================================================================================================



subroutine ssdg_periodic_org(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
  use structures
  use communication, only: comm_summation,comm_bcast
  use timer
  use hamiltonian, only: hpsi
  use eigen_subdiag_sub
  use sendrecv_grid, only: s_sendrecv_grid
  use pack_unpack, only: copy_data
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
  complex(8) :: wf_io1(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3))
  complex(8) :: wf_io2(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3))

  if ( SPIN_ORBIT_ON ) then
    call ssdg_periodic_so(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
    return
  end if

  if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ subspace_diag"

  call timer_begin(LOG_SSDG_PERIODIC_CALC)
  nspin = system%nspin
  no = system%no
  nk = system%nk
  is = mg%is
  ie = mg%ie
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  call timer_end(LOG_SSDG_PERIODIC_CALC)

  call timer_begin(LOG_SSDG_PERIODIC_HPSI)
  call hpsi(spsi,shpsi,info,mg,vlocal,system,stencil,srg,ppg)
  call timer_end(LOG_SSDG_PERIODIC_HPSI)

  call timer_begin(LOG_SSDG_PERIODIC_CALC)
  mat1 = 0d0
  if(info%if_divide_orbit) then
    do ik=ik_s,ik_e
    do ispin = 1, nspin
      do io1 = 1, no
        if (io_s<= io1 .and. io1 <= io_e) then
          call copy_data(spsi%zwf(:, :, :, ispin, io1, ik, 1),wf_io1)
        end if
        call timer_end(LOG_SSDG_PERIODIC_CALC)

        call timer_begin(LOG_SSDG_PERIODIC_COMM_COLL)
        call comm_bcast(wf_io1, info%icomm_o, info%irank_io(io1))
        call timer_end(LOG_SSDG_PERIODIC_COMM_COLL)

        call timer_begin(LOG_SSDG_PERIODIC_CALC)
        do io2 = 1, no
          if (io_s<= io2 .and. io2 <= io_e) then
            cbox = 0d0
            !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+:cbox)
            do iz=is(3),ie(3)
            do iy=is(2),ie(2)
            do ix=is(1),ie(1)
              cbox = cbox + conjg(wf_io1(ix,iy,iz)) * shpsi%zwf(ix,iy,iz,ispin,io2,ik,1)
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
        cbox = cbox + conjg(spsi%zwf(ix,iy,iz,ispin,io1,ik,1)) * shpsi%zwf(ix,iy,iz,ispin,io2,ik,1)
      end do
      end do
      end do
      mat1(io1,io2,ispin,ik) = cbox * system%hvol
    end do
    end do
    end do
    end do
  end if
  call timer_end(LOG_SSDG_PERIODIC_CALC)

  call timer_begin(LOG_SSDG_PERIODIC_COMM_COLL)
  call comm_summation(mat1,mat2,no**2*nspin*nk,info%icomm_rko)
  call timer_end(LOG_SSDG_PERIODIC_COMM_COLL)


  call timer_begin(LOG_SSDG_PERIODIC_EIGEN)
  do ik=ik_s,ik_e
  do ispin=1,nspin
    call eigen_subdiag_periodic(mat2(:,:,ispin,ik),evec(:,:,ispin,ik),no,ierr)
  end do
  end do
  call timer_end(LOG_SSDG_PERIODIC_EIGEN)


  call timer_begin(LOG_SSDG_PERIODIC_CALC)
!$omp workshare
  shpsi%zwf = 0d0
!$omp end workshare

  if(info%if_divide_orbit) then
    do ik=ik_s,ik_e
    do ispin=1,nspin
    do io2 = 1, no
      if (io_s<= io2 .and. io2 <= io_e) then
        call copy_data(spsi%zwf(:, :, :, ispin, io2, ik, 1),wf_io2)
      end if
      call timer_end(LOG_SSDG_PERIODIC_CALC)

      call timer_begin(LOG_SSDG_PERIODIC_COMM_COLL)
      call comm_bcast(wf_io2, info%icomm_o, info%irank_io(io2))
      call timer_end(LOG_SSDG_PERIODIC_COMM_COLL)

      call timer_begin(LOG_SSDG_PERIODIC_CALC)
      !$omp parallel do private(io1,iz,iy,ix) collapse(3)
      do io1=io_s,io_e
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
          shpsi%zwf(ix,iy,iz,ispin,io1,ik,1) = shpsi%zwf(ix,iy,iz,ispin,io1,ik,1) + evec(io2,io1,ispin,ik) * wf_io2(ix,iy,iz)
        end do
        end do
        end do
      end do
    end do
    end do
    end do
  else
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
  end if

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
  call timer_end(LOG_SSDG_PERIODIC_CALC)

  call timer_begin(LOG_SSDG_PERIODIC_COMM_COLL)
  call comm_summation(rbox1,rbox2,nspin*no*nk,info%icomm_rko)
  call timer_end(LOG_SSDG_PERIODIC_COMM_COLL)

  call timer_begin(LOG_SSDG_PERIODIC_CALC)
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
  call timer_end(LOG_SSDG_PERIODIC_CALC)

end subroutine ssdg_periodic_org

end module subspace_diagonalization
