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
module subspace_diagonalization_so

  use spin_orbit_global, only: SPIN_ORBIT_ON

  implicit none

  private
  public :: ssdg_periodic_so
  public :: SPIN_ORBIT_ON

contains

  subroutine ssdg_isolated(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
    use structures
    use salmon_communication, only: comm_summation,comm_bcast
    use timer
    use hamiltonian
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
    integer :: nspin,no,io,io1,io2,ispin,io_s,io_e,is(3),ie(3),ix,iy,iz,ierr
    real(8)   ,dimension(system%nspin,system%no) :: rbox1,rbox2
    real(8),dimension(system%no,system%no,system%nspin) :: mat1,mat2,evec
    real(8) :: cbox
    real(8) :: wf_io1(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3))
    real(8) :: wf_io2(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2),mg%is_array(3):mg%ie_array(3))

    call timer_begin(LOG_DIAG_TOTAL)

    if ( info%im_s /= 1 .or. info%im_e /= 1 ) stop "error: im/=1 @ subspace_diag"

    nspin = system%nspin
    no = system%no
    is = mg%is
    ie = mg%ie
    io_s = info%io_s
    io_e = info%io_e

  call hpsi(spsi,shpsi,info,mg,vlocal,system,stencil,srg,ppg)

  mat1 = 0d0

  if(info%if_divide_orbit) then
    do ispin = 1, nspin
      do io1 = 1, no
        if (io_s<= io1 .and. io1 <= io_e) then
          call copy_data(spsi%rwf(:, :, :, ispin, io1, 1, 1),wf_io1)
        end if
        call comm_bcast(wf_io1, info%icomm_o, info%irank_io(io1))
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

  call timer_begin(LOG_DIAG_ALLREDUCE)
  call comm_summation(mat1,mat2,no**2*nspin,info%icomm_rko)
  call timer_end(LOG_DIAG_ALLREDUCE)

  do ispin=1,nspin
    call timer_begin(LOG_DIAG_EIGEN)
    call eigen_subdiag(mat2(:,:,ispin),evec(:,:,ispin),no,ierr)
    call timer_end(LOG_DIAG_EIGEN)
  end do

!$omp workshare
  shpsi%rwf = 0d0
!$omp end workshare

  call timer_begin(LOG_DIAG_UPDATE)
  if(info%if_divide_orbit) then
    do ispin=1,nspin
    do io2 = 1, no
      if (io_s<= io2 .and. io2 <= io_e) then
        call copy_data(spsi%rwf(:, :, :, ispin, io2, 1, 1),wf_io2)
      end if
      call comm_bcast(wf_io2, info%icomm_o, info%irank_io(io2))
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
  call comm_summation(rbox1,rbox2,nspin*no,info%icomm_rko)

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
  call timer_end(LOG_DIAG_UPDATE)

  call timer_end(LOG_DIAG_TOTAL)

  end subroutine ssdg_isolated

!===================================================================================================================================

  subroutine ssdg_periodic_so(mg,system,info,stencil,spsi,shpsi,ppg,vlocal,srg)
    use structures
    use salmon_communication, only: comm_summation,comm_bcast
    use timer
    use hamiltonian
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
    integer :: nspin,no,nk,ik,io1,io2,ispin,ik_s,ik_e,io_s,io_e,is(3),ie(3),ix,iy,iz,ierr
    complex(8) :: cbox
    complex(8),allocatable :: mat1(:,:,:), mat2(:,:,:), evec(:,:,:)
    complex(8),allocatable :: wf_io1(:,:,:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)

    call timer_begin(LOG_DIAG_TOTAL)

    if ( info%im_s /= 1 .or. info%im_e /= 1 ) stop "error: im/=1 @ subspace_diag"

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

    allocate( wf_io1(size(spsi%zwf,1),size(spsi%zwf,2),size(spsi%zwf,3),size(spsi%zwf,4)) )
    wf_io1=zero
    allocate( mat1(no,no,ik_s:ik_e) ); mat1=zero
    allocate( mat2(no,no,ik_s:ik_e) ); mat2=zero
    allocate( evec(no,no,ik_s:ik_e) ); evec=zero

    if ( info%if_divide_orbit ) then

      do ik = ik_s, ik_e

        do io1 = 1, no

          if ( io_s<= io1 .and. io1 <= io_e ) then
            call copy_data( spsi%zwf(:,:,:,:,io1,ik,1), wf_io1 )
          end if
          do ispin=1,nspin
            call comm_bcast( wf_io1(:,:,:,ispin), info%icomm_o, info%irank_io(io1) )
          end do

          do io2 = 1, no

            if ( io_s <= io2 .and. io2 <= io_e) then
              cbox=zero
              !$omp parallel do private(ispin,iz,iy,ix) collapse(2) reduction(+:cbox)
              do ispin=1,nspin
              do iz=is(3),ie(3)
              do iy=is(2),ie(2)
              do ix=is(1),ie(1)
                cbox = cbox + conjg( wf_io1(ix,iy,iz,ispin) ) * shpsi%zwf(ix,iy,iz,ispin,io2,ik,1)
              end do
              end do
              end do
              end do
              mat1(io1,io2,ik) = cbox * system%hvol
            end if

          end do !io2

        end do !io1

      end do !ik

    else ![ .not. info%if_divide_orbit ]

      !!$omp parallel do private(ik,io1,io2,ispin,cbox,iz,iy,ix) collapse(4)
      do ik=ik_s,ik_e
      do io1=io_s,io_e
      do io2=io_s,io_e
        cbox=zero
        do ispin=1,nspin
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
          cbox = cbox + conjg( spsi%zwf(ix,iy,iz,ispin,io1,ik,1) ) * shpsi%zwf(ix,iy,iz,ispin,io2,ik,1)
        end do
        end do
        end do
        end do
        mat1(io1,io2,ik) = cbox * system%hvol
      end do
      end do
      end do

    end if

    call timer_begin(LOG_DIAG_ALLREDUCE)
    call comm_summation( mat1, mat2, size(mat1) ,info%icomm_rko )
    call timer_end(LOG_DIAG_ALLREDUCE)
 
    do ik=ik_s,ik_e
      call timer_begin(LOG_DIAG_EIGEN)
      call eigen_subdiag_periodic( mat2(:,:,ik), evec(:,:,ik), no, ierr )
      call timer_end(LOG_DIAG_EIGEN)
    end do

!$omp workshare
    shpsi%zwf = zero
!$omp end workshare

    call timer_begin(LOG_DIAG_UPDATE)
  
    if ( info%if_divide_orbit ) then

      do ik=ik_s,ik_e

        do io1 = 1, no

          if ( io_s <= io1 .and. io1 <= io_e ) then
            call copy_data( spsi%zwf(:,:,:,:,io1,ik,1), wf_io1 )
          end if
          do ispin=1,nspin
            call comm_bcast( wf_io1(:,:,:,ispin), info%icomm_o, info%irank_io(io1) )
          end do

          !$omp parallel do private(io1,ispin,iz,iy,ix) collapse(3)
          do io2=io_s,io_e
            do ispin=1,nspin
            do iz=is(3),ie(3)
            do iy=is(2),ie(2)
            do ix=is(1),ie(1)
              shpsi%zwf(ix,iy,iz,ispin,io2,ik,1) = shpsi%zwf(ix,iy,iz,ispin,io2,ik,1) &
                   + evec(io1,io2,ik) * wf_io1(ix,iy,iz,ispin)
            end do
            end do
            end do
            end do
          end do !io2

        end do !io1

      end do !ik

    else

      do ik=ik_s,ik_e 
        !$omp parallel do private(io1,io2,ispin,iz,iy,ix) collapse(3)
        do io1=io_s,io_e
          do io2=io_s,io_e
            do ispin=1,nspin
            do iz=is(3),ie(3)
            do iy=is(2),ie(2)
            do ix=is(1),ie(1)
            shpsi%zwf(ix,iy,iz,ispin,io1,ik,1) = shpsi%zwf(ix,iy,iz,ispin,io1,ik,1) &
                 + evec(io2,io1,ik) * spsi%zwf(ix,iy,iz,ispin,io2,ik,1)
            end do
            end do
            end do
            end do
          end do !io2
        end do !io1
      end do !ik

    end if

    call timer_end(LOG_DIAG_UPDATE)

    !$omp parallel workshare
    spsi%zwf(:,:,:,:,:,:,1) = shpsi%zwf(:,:,:,:,:,:,1)
    !$omp end parallel workshare

    deallocate( evec )
    deallocate( mat2 )
    deallocate( mat1 )
    deallocate( wf_io1 )

    call timer_end(LOG_DIAG_TOTAL)

  end subroutine ssdg_periodic_so

end module subspace_diagonalization_so
