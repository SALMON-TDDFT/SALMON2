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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

#include "config.h"

module density_matrix
  implicit none
  integer,private,parameter :: Nd = 4

contains

  subroutine calc_density(system,rho,psi,info,mg)
    use structures
    use communication, only: comm_summation
    use parallelization, only: get_thread_id,get_nthreads
    use misc_routines, only: ceiling_pow2
    use timer
    use sym_rho_sub, only: sym_rho
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_orbital),intent(in) :: psi
    type(s_scalar) :: rho(system%nspin,info%im_s:info%im_e)
    !
    integer :: im,ispin,ik,io,is(3),ie(3),nsize,nspin,tid,ix,iy,iz,nthreads
    real(8) :: wrk2
    real(8),allocatable :: wrk(:,:,:,:)

    call timer_begin(LOG_DENSITY_CALC)
    nspin = system%nspin

#ifdef FORTRAN_COMPILER_HAS_2MB_ALIGNED_ALLOCATION
!dir$ attributes align : 2097152 :: wrk
#endif

    is = mg%is
    ie = mg%ie
    nsize = mg%num(1) * mg%num(2) * mg%num(3)
    nthreads = get_nthreads()

    allocate(wrk(is(1):ie(1),is(2):ie(2),is(3):ie(3),0:ceiling_pow2(nthreads)-1))
    call timer_end(LOG_DENSITY_CALC)

    if(allocated(psi%rwf)) then

      do im=info%im_s,info%im_e
      do ispin=1,nspin
        call timer_begin(LOG_DENSITY_CALC)
        tid = 0
!$omp parallel private(ik,io,iz,iy,ix,wrk2) firstprivate(tid)
!$      tid = get_thread_id()
        wrk(:,:,:,tid) = 0.d0

!$omp do collapse(4)
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
          wrk2 = abs( psi%rwf(ix,iy,iz,ispin,io,ik,im) )**2
          wrk(ix,iy,iz,tid) = wrk(ix,iy,iz,tid) + wrk2 * system%rocc(io,ik,ispin)*system%wtk(ik)
        end do
        end do
        end do
        end do
        end do
!$omp end do

        ix = size(wrk,4)/2
        do while(ix > 0)
          if(tid < ix .and. tid + ix < nthreads) then
            wrk(:,:,:,tid) = wrk(:,:,:,tid) + wrk(:,:,:,tid + ix)
          end if
          ix = ix/2
!$omp barrier
        end do

!$omp end parallel
        call timer_end(LOG_DENSITY_CALC)

        call timer_begin(LOG_DENSITY_COMM_COLL)
        call comm_summation(wrk(:,:,:,0),rho(ispin,im)%f(:,:,:),nsize,info%icomm_ko)
        call timer_end(LOG_DENSITY_COMM_COLL)
      end do
      end do

    else

      do im=info%im_s,info%im_e
      do ispin=1,nspin
        call timer_begin(LOG_DENSITY_CALC)
        tid = 0
!$omp parallel private(ik,io,iz,iy,ix,wrk2) firstprivate(tid)
!$      tid = get_thread_id()
        wrk(:,:,:,tid) = 0.d0

!$omp do collapse(4)
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
          wrk2 = abs( psi%zwf(ix,iy,iz,ispin,io,ik,im) )**2
          wrk(ix,iy,iz,tid) = wrk(ix,iy,iz,tid) + wrk2 * system%rocc(io,ik,ispin)*system%wtk(ik)
        end do
        end do
        end do
        end do
        end do
!$omp end do

        ix = size(wrk,4)/2
        do while(ix > 0)
          if(tid < ix .and. tid + ix < nthreads) then
            wrk(:,:,:,tid) = wrk(:,:,:,tid) + wrk(:,:,:,tid + ix)
          end if
          ix = ix/2
!$omp barrier
        end do

!$omp end parallel
        call timer_end(LOG_DENSITY_CALC)

        call timer_begin(LOG_DENSITY_COMM_COLL)
        call comm_summation(wrk(:,:,:,0),rho(ispin,im)%f(:,:,:),nsize,info%icomm_ko)
        call timer_end(LOG_DENSITY_COMM_COLL)

        call timer_begin(LOG_DENSITY_CALC)
        call sym_rho( rho(ispin,im)%f(:,:,:) )
        call timer_end(LOG_DENSITY_CALC)
      end do
      end do

    end if

    deallocate(wrk)
    return
  end subroutine calc_density

!===================================================================================================================================

  subroutine calc_current(system,mg,stencil,info,srg,psi,ppg,curr)
    use structures
    use salmon_global, only: yn_jm
    use sendrecv_grid, only: update_overlap_complex8
    use communication, only: comm_summation
    use nonlocal_potential, only: calc_uVpsi_rdivided
    use pseudo_pt_current_so, only: SPIN_ORBIT_ON, calc_current_nonlocal_so &
                                  , calc_current_nonlocal_rdivided_so
    use sym_vector_sub, only: sym_vector_xyz
    use code_optimization, only: current_omp_mode
    use timer
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid)  ,intent(in) :: mg
    type(s_stencil),intent(in) :: stencil
    type(s_parallel_info),intent(in) :: info
    type(s_sendrecv_grid)      :: srg
    type(s_orbital)            :: psi
    type(s_pp_grid),intent(in) :: ppg
    real(8) :: curr(3,system%nspin,info%im_s:info%im_e)
    !
    integer :: ispin,im,ik,io,nspin,ngrid
    real(8),dimension(3) :: wrk1,wrk2,wrk3,wrk4
    real(8) :: BT(3,3),kAc(3)
    complex(8),allocatable :: uVpsibox (:,:,:,:,:)
    complex(8),allocatable :: uVpsibox2(:,:,:,:,:)
    complex(8),allocatable :: uVpsi(:)

    call timer_begin(LOG_CURRENT_CALC)
#ifdef FORTRAN_COMPILER_HAS_2MB_ALIGNED_ALLOCATION
!dir$ attributes align : 2097152 :: uVpsibox, uVpsibox2
#endif

    nspin = system%nspin
    ngrid = system%ngrid

    BT = transpose(system%rmatrix_B)
    call timer_end(LOG_CURRENT_CALC)

    call timer_begin(LOG_CURRENT_CALC_UVPSI_RDIVIDED)
    if (info%if_divide_rspace .and. yn_jm=='n' .and. .not.SPIN_ORBIT_ON) then
      call calc_uVpsi_rdivided(nspin,info,ppg,psi,uVpsibox,uVpsibox2)
      allocate(uVpsi(ppg%Nlma))
    end if
    call timer_end(LOG_CURRENT_CALC_UVPSI_RDIVIDED)

  ! overlap region communication
    call timer_begin(LOG_CURRENT_COMM_HALO)
    if(info%if_divide_rspace) then
      call update_overlap_complex8(srg, mg, psi%zwf)
    end if
    call timer_end(LOG_CURRENT_COMM_HALO)

    do im=info%im_s,info%im_e
    do ispin=1,nspin
      call timer_begin(LOG_CURRENT_CALC)
      wrk4 = 0d0
!$omp parallel do collapse(2) default(none) &
!$omp             private(ik,io,kAc,wrk1,wrk2,wrk3,uVpsi) &
!$omp             shared(info,system,mg,stencil,ppg,psi,uVpsibox2,BT,im,ispin,yn_jm) &
!$omp             shared(SPIN_ORBIT_ON) &
!$omp             reduction(+:wrk4) if(current_omp_mode)
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e

        kAc(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
        call stencil_current(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz,stencil%coef_nab &
                            ,kAc,psi%zwf(:,:,:,ispin,io,ik,im),wrk1,wrk2)
        wrk2 = matmul(BT,wrk2)

        if ( yn_jm == 'n' ) then
          if ( SPIN_ORBIT_ON ) then
            if ( info%if_divide_rspace ) then
              call calc_current_nonlocal_rdivided_so &
                   ( wrk3,psi%zwf(:,:,:,:,io,ik,im),ppg,mg%is_array,mg%ie_array,ik,info%icomm_r )
            else
              call calc_current_nonlocal_so &
                   ( wrk3,psi%zwf(:,:,:,:,io,ik,im),ppg,mg%is_array,mg%ie_array,ik )
            end if
          else
            if ( info%if_divide_rspace)then
              uVpsi(:) = uVpsibox2(ispin,io,ik,im,:)
              call calc_current_nonlocal_rdivided(wrk3,psi%zwf(:,:,:,ispin,io,ik,im),ppg,mg%is_array,mg%ie_array,ik,uVpsi)
            else
              call calc_current_nonlocal         (wrk3,psi%zwf(:,:,:,ispin,io,ik,im),ppg,mg%is_array,mg%ie_array,ik)
            end if
          end if
        else
          wrk3=0d0
        end if

        wrk4 = wrk4 + (wrk1 + wrk2 + wrk3) * system%rocc(io,ik,ispin)*system%wtk(ik)

      end do
      end do
!$omp end parallel do
      call timer_end(LOG_CURRENT_CALC)

      call timer_begin(LOG_CURRENT_COMM_COLL)
      call comm_summation(wrk4,wrk1,3,info%icomm_rko)

      curr(:,ispin,im) = wrk1 / dble(ngrid) ! ngrid = aLxyz/Hxyz
      call timer_end(LOG_CURRENT_COMM_COLL)

      call timer_begin(LOG_CURRENT_CALC)
      call sym_vector_xyz( curr(:,ispin,im) )
      call timer_end(LOG_CURRENT_CALC)
    end do
    end do

    if (info%if_divide_rspace .and. yn_jm=='n' .and. .not.SPIN_ORBIT_ON) deallocate(uVpsibox,uVpsibox2,uVpsi)

    return

  contains

    subroutine stencil_current(is_array,ie_array,is,ie,idx,idy,idz,nabt,kAc,psi,j1,j2)
      integer   ,intent(in) :: is_array(3),ie_array(3),is(3),ie(3) &
                              ,idx(is(1)-Nd:ie(1)+Nd),idy(is(2)-Nd:ie(2)+Nd),idz(is(3)-Nd:ie(3)+Nd)
      real(8)   ,intent(in) :: nabt(Nd,3),kAc(3)
      complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
      real(8)               :: j1(3),j2(3)
      !
      integer :: ix,iy,iz
      real(8) :: rtmp
      complex(8) :: cpsi,xtmp,ytmp,ztmp
      rtmp = 0d0
      xtmp = 0d0
      ytmp = 0d0
      ztmp = 0d0
!$omp parallel do collapse(2) private(iz,iy,ix,cpsi) reduction(+:rtmp,xtmp,ytmp,ztmp)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)

!OCL swp
      do ix=is(1),ie(1)
        rtmp = rtmp + abs(psi(ix,iy,iz))**2
      end do

!OCL swp
      do ix=is(1),ie(1)
        cpsi = conjg(psi(ix,iy,iz))
        xtmp = xtmp + nabt(1,1) * cpsi * psi(idx(ix+1),iy,iz) &
                    + nabt(2,1) * cpsi * psi(idx(ix+2),iy,iz) &
                    + nabt(3,1) * cpsi * psi(idx(ix+3),iy,iz) &
                    + nabt(4,1) * cpsi * psi(idx(ix+4),iy,iz)
      end do

!OCL swp
      do ix=is(1),ie(1)
        cpsi = conjg(psi(ix,iy,iz))
        ytmp = ytmp + nabt(1,2) * cpsi * psi(ix,idy(iy+1),iz) &
                    + nabt(2,2) * cpsi * psi(ix,idy(iy+2),iz) &
                    + nabt(3,2) * cpsi * psi(ix,idy(iy+3),iz) &
                    + nabt(4,2) * cpsi * psi(ix,idy(iy+4),iz)
      end do

!OCL swp
      do ix=is(1),ie(1)
        cpsi = conjg(psi(ix,iy,iz))
        ztmp = ztmp + nabt(1,3) * cpsi * psi(ix,iy,idz(iz+1)) &
                    + nabt(2,3) * cpsi * psi(ix,iy,idz(iz+2)) &
                    + nabt(3,3) * cpsi * psi(ix,iy,idz(iz+3)) &
                    + nabt(4,3) * cpsi * psi(ix,iy,idz(iz+4))
      end do

      end do
      end do
!$omp end parallel do
      j1 = kAc(:) * rtmp
      j2(1) = aimag(xtmp * 2d0)
      j2(2) = aimag(ytmp * 2d0)
      j2(3) = aimag(ztmp * 2d0)
      return
    end subroutine stencil_current

  end subroutine calc_current

  subroutine calc_current_nonlocal(jw,psi,ppg,is_array,ie_array,ik)
    use structures
    implicit none
    integer   ,intent(in) :: is_array(3),ie_array(3),ik
    complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
    type(s_pp_grid),intent(in) :: ppg
    real(8)               :: jw(3)
    !
    integer    :: ilma,ia,j,ix,iy,iz
    real(8)    :: x,y,z
    complex(8) :: uVpsi,uVpsi_r(3)
    jw = 0d0
!$omp parallel do private(ilma,ia,uVpsi,uVpsi_r,j,x,y,z,ix,iy,iz) reduction(+:jw)
    do ilma=1,ppg%Nlma
      ia=ppg%ia_tbl(ilma)
      uVpsi = 0d0
      uVpsi_r = 0d0
      do j=1,ppg%Mps(ia)
        x = ppg%Rxyz(1,j,ia)
        y = ppg%Rxyz(2,j,ia)
        z = ppg%Rxyz(3,j,ia)
        ix = ppg%Jxyz(1,j,ia)
        iy = ppg%Jxyz(2,j,ia)
        iz = ppg%Jxyz(3,j,ia)
        uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * psi(ix,iy,iz)
        uVpsi_r(1) = uVpsi_r(1) + conjg(ppg%zekr_uV(j,ilma,ik)) * x * psi(ix,iy,iz)
        uVpsi_r(2) = uVpsi_r(2) + conjg(ppg%zekr_uV(j,ilma,ik)) * y * psi(ix,iy,iz)
        uVpsi_r(3) = uVpsi_r(3) + conjg(ppg%zekr_uV(j,ilma,ik)) * z * psi(ix,iy,iz)
      end do
      uVpsi = uVpsi * ppg%rinv_uvu(ilma)
      jw = jw + aimag(conjg(uVpsi_r)*uVpsi)
    end do
!$omp end parallel do
    jw = jw * 2d0
    return
  end subroutine calc_current_nonlocal

  subroutine calc_current_nonlocal_rdivided(jw,psi,ppg,is_array,ie_array,ik,uVpsibox)
    use structures
    implicit none
    integer   ,intent(in) :: is_array(3),ie_array(3),ik
    complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
    type(s_pp_grid),intent(in) :: ppg
    complex(8),intent(in) :: uVpsibox(ppg%Nlma)
    real(8)               :: jw(3)
    !
    integer    :: ilocal,ilma,ia,j,ix,iy,iz
    real(8)    :: x,y,z
    complex(8) :: uVpsi,uVpsi_r(3)
    jw = 0d0
!$omp parallel do private(ilocal,ilma,ia,uVpsi,uVpsi_r,j,x,y,z,ix,iy,iz) reduction(+:jw)
    do ilocal=1,ppg%ilocal_nlma
      ilma=ppg%ilocal_nlma2ilma(ilocal)
      ia  =ppg%ilocal_nlma2ia  (ilocal)
      uVpsi_r = 0d0
      do j=1,ppg%Mps(ia)
        x = ppg%Rxyz(1,j,ia)
        y = ppg%Rxyz(2,j,ia)
        z = ppg%Rxyz(3,j,ia)
        ix = ppg%Jxyz(1,j,ia)
        iy = ppg%Jxyz(2,j,ia)
        iz = ppg%Jxyz(3,j,ia)
        uVpsi_r(1) = uVpsi_r(1) + conjg(ppg%zekr_uV(j,ilma,ik)) * x * psi(ix,iy,iz)
        uVpsi_r(2) = uVpsi_r(2) + conjg(ppg%zekr_uV(j,ilma,ik)) * y * psi(ix,iy,iz)
        uVpsi_r(3) = uVpsi_r(3) + conjg(ppg%zekr_uV(j,ilma,ik)) * z * psi(ix,iy,iz)
      end do
      uVpsi = uVpsibox(ilma)
      jw = jw + aimag(conjg(uVpsi_r)*uVpsi)
    end do
!$omp end parallel do
    jw = jw * 2d0
    return
  end subroutine calc_current_nonlocal_rdivided

!===================================================================================================================================

  subroutine calc_microscopic_current(system,mg,stencil,info,psi,curr)
    use structures
    use communication, only: comm_summation
    use timer
    implicit none
    type(s_dft_system)      ,intent(in) :: system
    type(s_rgrid)           ,intent(in) :: mg
    type(s_stencil)         ,intent(in) :: stencil
    type(s_parallel_info),intent(in) :: info
    type(s_orbital)         ,intent(in) :: psi
    type(s_vector)                      :: curr ! electron number current density (without rho*A/c)
    !
    integer :: ispin,im,ik,io,is(3),ie(3),nsize,nspin,ix,iy,iz
    real(8) :: kAc(3)
    real(8),allocatable :: wrk(:,:,:,:),wrk2(:,:,:,:)

    call timer_begin(LOG_MCURRENT_CALC)
    nspin = system%nspin

    if(.not. psi%update_zwf_overlap) stop "halo of wavefunction is not updated"
    if(info%im_s/=1 .or. info%im_e/=1) stop "error: im_s, im_e @ calc_microscopic_current"
    im = 1

    is = mg%is
    ie = mg%ie
    allocate(wrk(3,is(1):ie(1),is(2):ie(2),is(3):ie(3)),wrk2(3,is(1):ie(1),is(2):ie(2),is(3):ie(3)))
    nsize = 3* mg%num(1) * mg%num(2) * mg%num(3)

    wrk2 = 0d0
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do ispin=1,nspin

      kAc(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
      call micro_current(mg%is_array,mg%ie_array,is,ie,mg%idx,mg%idy,mg%idz, &
      & stencil%coef_nab,kAc,psi%zwf(:,:,:,ispin,io,ik,im),wrk)

!$omp parallel do collapse(2) private(ix,iy,iz)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        wrk2(1,ix,iy,iz) = wrk2(1,ix,iy,iz) + wrk(1,ix,iy,iz) * system%rocc(io,ik,ispin)*system%wtk(ik)
        wrk2(2,ix,iy,iz) = wrk2(2,ix,iy,iz) + wrk(2,ix,iy,iz) * system%rocc(io,ik,ispin)*system%wtk(ik)
        wrk2(3,ix,iy,iz) = wrk2(3,ix,iy,iz) + wrk(3,ix,iy,iz) * system%rocc(io,ik,ispin)*system%wtk(ik)
      end do
      end do
      end do

!     call nonlocal_part

    end do
    end do
    end do

    call timer_end(LOG_MCURRENT_CALC)

    call timer_begin(LOG_MCURRENT_COMM_COLL)
    call comm_summation(wrk2,curr%v,nsize,info%icomm_ko)
    call timer_end(LOG_MCURRENT_COMM_COLL)

    deallocate(wrk,wrk2)
    return

  contains
  
# define DX(dt) idx(ix+(dt)),iy,iz
# define DY(dt) ix,idy(iy+(dt)),iz
# define DZ(dt) ix,iy,idz(iz+(dt))

    subroutine micro_current(is_array,ie_array,is,ie,idx,idy,idz,nabt,kAc,tpsi,jw)
      implicit none
      integer   ,intent(in) :: is_array(3),ie_array(3),is(3),ie(3), &
                             & idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
      real(8)   ,intent(in) :: nabt(Nd,3),kAc(3)
      complex(8),intent(in) :: tpsi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
      real(8)               :: jw(3,is(1):ie(1),is(2):ie(2),is(3):ie(3))
      !
      integer :: ix,iy,iz
      complex(8) :: xtmp,ytmp,ztmp,px,py,pz
!$omp parallel do collapse(2) private(iz,iy,ix,xtmp,ytmp,ztmp,px,py,pz)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)

!OCL swp
        do ix=is(1),ie(1)
          px = conjg(tpsi(ix,iy,iz))
          xtmp = nabt(1,1) * ( tpsi(DX(1)) - tpsi(DX(-1)) ) &
               + nabt(2,1) * ( tpsi(DX(2)) - tpsi(DX(-2)) ) &
               + nabt(3,1) * ( tpsi(DX(3)) - tpsi(DX(-3)) ) &
               + nabt(4,1) * ( tpsi(DX(4)) - tpsi(DX(-4)) )
          jw(1,ix,iy,iz) = aimag(px*xtmp) + kAc(1) * abs(px)**2
        end do

!OCL swp
        do ix=is(1),ie(1)
          py = conjg(tpsi(ix,iy,iz))
          ytmp = nabt(1,2) * ( tpsi(DY(1)) - tpsi(DY(-1)) ) &
               + nabt(2,2) * ( tpsi(DY(2)) - tpsi(DY(-2)) ) &
               + nabt(3,2) * ( tpsi(DY(3)) - tpsi(DY(-3)) ) &
               + nabt(4,2) * ( tpsi(DY(4)) - tpsi(DY(-4)) )
          jw(2,ix,iy,iz) = aimag(py*ytmp) + kAc(2) * abs(py)**2
        end do

!OCL swp
        do ix=is(1),ie(1)
          pz = conjg(tpsi(ix,iy,iz))
          ztmp = nabt(1,3) * ( tpsi(DZ(1)) - tpsi(DZ(-1)) ) &
               + nabt(2,3) * ( tpsi(DZ(2)) - tpsi(DZ(-2)) ) &
               + nabt(3,3) * ( tpsi(DZ(3)) - tpsi(DZ(-3)) ) &
               + nabt(4,3) * ( tpsi(DZ(4)) - tpsi(DZ(-4)) )
          jw(3,ix,iy,iz) = aimag(pz*ztmp) + kAc(3) * abs(pz)**2
        end do

      end do
      end do
!$omp end parallel do
      return
    end subroutine micro_current

    subroutine kvec_part(is_array,ie_array,is,ie,k,psi,jw)
      implicit none
      integer   ,intent(in) :: is_array(3),ie_array(3),is(3),ie(3)
      real(8)   ,intent(in) :: k(3)
      complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
      real(8)               :: jw(3,is(1):ie(1),is(2):ie(2),is(3):ie(3))
      !
      integer :: ik,io,ix,iy,iz
!$omp parallel do collapse(2) private(iz,iy,ix)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        jw(:,ix,iy,iz) = k(:) * abs(psi(ix,iy,iz))**2
      end do
      end do
      end do
!$omp end parallel do
      return
    end subroutine kvec_part

    subroutine stencil_current(jw,zdm,nabt,is,ie,ndir)
      implicit none
      integer   ,intent(in) :: is(3),ie(3),ndir
      real(8)   ,intent(in) :: nabt(Nd,3)
      complex(8),intent(in) :: zdm(Nd,ndir,is(1)-Nd:ie(1),is(2)-Nd:ie(2),is(3)-Nd:ie(3))
      real(8)               :: jw(3,is(1):ie(1),is(2):ie(2),is(3):ie(3))
      !
      integer :: ix,iy,iz
      complex(8) :: xtmp,ytmp,ztmp
!$omp parallel do collapse(2) private(iz,iy,ix,xtmp,ytmp,ztmp)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)

!OCL swp
        do ix=is(1),ie(1)
          xtmp = nabt(1,1) * ( zdm(1,1,ix,iy,iz) - conjg(zdm(1,1,ix-1,iy,iz)) ) & 
               + nabt(2,1) * ( zdm(2,1,ix,iy,iz) - conjg(zdm(2,1,ix-2,iy,iz)) ) &
               + nabt(3,1) * ( zdm(3,1,ix,iy,iz) - conjg(zdm(3,1,ix-3,iy,iz)) ) &
               + nabt(4,1) * ( zdm(4,1,ix,iy,iz) - conjg(zdm(4,1,ix-4,iy,iz)) )
          jw(1,ix,iy,iz) = aimag(xtmp)
        end do

!OCL swp
        do ix=is(1),ie(1)
          ytmp = nabt(1,2) * ( zdm(1,2,ix,iy,iz) - conjg(zdm(1,2,ix,iy-1,iz)) ) &
               + nabt(2,2) * ( zdm(2,2,ix,iy,iz) - conjg(zdm(2,2,ix,iy-2,iz)) ) &
               + nabt(3,2) * ( zdm(3,2,ix,iy,iz) - conjg(zdm(3,2,ix,iy-3,iz)) ) &
               + nabt(4,2) * ( zdm(4,2,ix,iy,iz) - conjg(zdm(4,2,ix,iy-4,iz)) )
          jw(2,ix,iy,iz) = aimag(ytmp)
        end do

!OCL swp
        do ix=is(1),ie(1)
          ztmp = nabt(1,3) * ( zdm(1,3,ix,iy,iz) - conjg(zdm(1,3,ix,iy,iz-1)) ) &
               + nabt(2,3) * ( zdm(2,3,ix,iy,iz) - conjg(zdm(2,3,ix,iy,iz-2)) ) &
               + nabt(3,3) * ( zdm(3,3,ix,iy,iz) - conjg(zdm(3,3,ix,iy,iz-3)) ) &
               + nabt(4,3) * ( zdm(4,3,ix,iy,iz) - conjg(zdm(4,3,ix,iy,iz-4)) )
          jw(3,ix,iy,iz) = aimag(ztmp)
        end do

      end do
      end do
!$omp end parallel do
      return
    end subroutine stencil_current

  end subroutine calc_microscopic_current

end module
