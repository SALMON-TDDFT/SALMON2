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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

module density_matrix
  implicit none
  integer,private,parameter :: Nd = 4

contains

! density matrix: dmat(r,dr) = sum occ psi(r+dr) * conjg(psi(r))
! dmat(r,-dr) = conjg(dmat(r-dr,dr))
! j(r) = sum occ aimag( conjg(psi(r))* (sum nabt*psi)(r) ) = aimag( sum_dr nabt(dr)* dmat(r,dr) )

  subroutine calc_density_matrix(nspin,info,mg,srg,psi,dmat)
    use structures
    use sendrecv_grid, only: s_sendrecv_grid, update_overlap_real8, update_overlap_complex8
    use salmon_communication, only: comm_summation
    implicit none
    integer        ,intent(in) :: nspin
    type(s_orbital_parallel),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_sendrecv_grid)      :: srg
    type(s_orbital)            :: psi
    type(s_dmatrix)            :: dmat
    !
    integer :: im,ispin,ik,io,is(3),ie(3),nsize
    complex(8),allocatable :: wrk(:,:,:,:,:),wrk2(:,:,:,:,:)

    is = mg%is
    ie = mg%ie
    nsize = Nd* mg%ndir * (mg%num(1)+Nd) * (mg%num(2)+Nd) * (mg%num(3)+Nd)

    if(allocated(psi%rwf)) then
      allocate(psi%zwf(mg%is_array(1):mg%ie_array(1) &
                      ,mg%is_array(2):mg%ie_array(2) &
                      ,mg%is_array(3):mg%ie_array(3) &
                      ,nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e))
      psi%zwf = cmplx(psi%rwf)
    end if

    allocate( wrk(Nd,mg%ndir,is(1)-Nd:ie(1),is(2)-Nd:ie(2),is(3)-Nd:ie(3)) &
            ,wrk2(Nd,mg%ndir,is(1)-Nd:ie(1),is(2)-Nd:ie(2),is(3)-Nd:ie(3)) )

  ! overlap region communication
    if(info%if_divide_rspace) then
      call update_overlap_complex8(srg, mg, psi%zwf)
    end if

    do im=info%im_s,info%im_e
    do ispin=1,nspin
      wrk = 0d0
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
        call calc_dm(wrk2,psi%zwf(:,:,:,ispin,io,ik,im),mg%is_array,mg%ie_array,is,ie,mg%idx,mg%idy,mg%idz,mg%ndir)
        wrk = wrk + wrk2 * info%occ(io,ik,ispin,im)
      end do
      end do
      call comm_summation(wrk,wrk2,nsize,info%icomm_ko)
      dmat%zrho_mat(:,:,:,:,:,ispin,im) = wrk2(:,:,:,:,:)
    end do
    end do

    if(allocated(psi%rwf)) deallocate(psi%zwf)
    deallocate(wrk,wrk2)
    return

  contains
    subroutine calc_dm(zdm,psi,is_array,ie_array,is,ie,idx,idy,idz,ndir)
      implicit none
      integer   ,intent(in) :: is_array(3),ie_array(3),is(3),ie(3) &
                              ,idx(is(1)-Nd:ie(1)+Nd),idy(is(2)-Nd:ie(2)+Nd),idz(is(3)-Nd:ie(3)+Nd),ndir
      complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
      complex(8)            :: zdm(Nd,ndir,is(1)-Nd:ie(1),is(2)-Nd:ie(2),is(3)-Nd:ie(3))
      !
      integer :: ix,iy,iz,ii
      zdm = 0d0
      do iz=is(3)-Nd,ie(3)
      do iy=is(2)-Nd,ie(2)
      do ix=is(1)-Nd,ie(1)
        do ii=1,Nd
        ! dir = 1,2,3 = xx,yy,zz (yz,zx,xy)
          zdm(ii,1,ix,iy,iz) = psi(idx(ix+ii),iy,iz) * conjg( psi(ix,iy,iz) )
          zdm(ii,2,ix,iy,iz) = psi(ix,idy(iy+ii),iz) * conjg( psi(ix,iy,iz) )
          zdm(ii,3,ix,iy,iz) = psi(ix,iy,idz(iz+ii)) * conjg( psi(ix,iy,iz) )
        end do
      end do
      end do
      end do
      return
    end subroutine calc_dm
  end subroutine calc_density_matrix

!===================================================================================================================================

  subroutine calc_density(rho,psi,info,mg,nspin)
    use structures
    use salmon_communication, only: comm_summation
    use salmon_parallel, only: get_thread_id,get_nthreads
    use misc_routines, only: ceiling_pow2
    use timer
    implicit none
    integer        ,intent(in) :: nspin
    type(s_orbital_parallel),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_orbital),intent(in) :: psi
    type(s_scalar) :: rho(nspin,info%im_s:info%im_e)
    !
    integer :: im,ispin,ik,io,is(3),ie(3),nsize,tid,ix,iy,iz,nthreads
    real(8) :: wrk2
    real(8),allocatable :: wrk(:,:,:,:)
    is = mg%is
    ie = mg%ie
    nsize = mg%num(1) * mg%num(2) * mg%num(3)
    nthreads = get_nthreads()

    allocate(wrk(is(1):ie(1),is(2):ie(2),is(3):ie(3),0:ceiling_pow2(nthreads)-1))

    if(allocated(psi%rwf)) then

      do im=info%im_s,info%im_e
      do ispin=1,nspin
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
          wrk(ix,iy,iz,tid) = wrk(ix,iy,iz,tid) + wrk2 * info%occ(io,ik,ispin,im)
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
        call timer_begin(LOG_ALLREDUCE_RHO)
        call comm_summation(wrk(:,:,:,0),rho(ispin,im)%f(:,:,:),nsize,info%icomm_ko)
        call timer_end(LOG_ALLREDUCE_RHO)
      end do
      end do

    else

      do im=info%im_s,info%im_e
      do ispin=1,nspin
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
          wrk(ix,iy,iz,tid) = wrk(ix,iy,iz,tid) + wrk2 * info%occ(io,ik,ispin,im)
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
        call timer_begin(LOG_ALLREDUCE_RHO)
        call comm_summation(wrk(:,:,:,0),rho(ispin,im)%f(:,:,:),nsize,info%icomm_ko)
        call timer_end(LOG_ALLREDUCE_RHO)
      end do
      end do

    end if

    deallocate(wrk)
    return
  end subroutine calc_density

!===================================================================================================================================

  subroutine calc_current(nspin,ngrid,mg,stencil,info,psi,ppg,dmat,curr)
    use structures
    use salmon_communication, only: comm_summation
    implicit none
    integer        ,intent(in) :: nspin,ngrid
    type(s_rgrid)  ,intent(in) :: mg
    type(s_stencil),intent(in) :: stencil
    type(s_orbital_parallel),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    type(s_pp_grid),intent(in) :: ppg
    type(s_dmatrix),intent(in) :: dmat
    real(8) :: curr(3,nspin,info%im_s:info%im_e)
    !
    integer :: ispin,im,ik,io
    real(8) :: wrk(3),wrk2(3)
    complex(8),allocatable :: uVpsibox (:,:,:,:,:)
    complex(8),allocatable :: uVpsibox2(:,:,:,:,:)

    if(info%if_divide_rspace) call nonlocal_part_rdivided1(nspin,info,ppg,psi,uVpsibox,uVpsibox2)

    do im=info%im_s,info%im_e
    do ispin=1,nspin

      wrk2 = 0d0
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e

        call kvec_part(wrk,psi%zwf(:,:,:,ispin,io,ik,im),stencil%vec_kAc(:,ik),mg%is_array,mg%ie_array,mg%is,mg%ie)
        wrk2 = wrk2 + wrk * info%occ(io,ik,ispin,im)

        if(info%if_divide_rspace) then
          call nonlocal_part_rdivided2(wrk,psi%zwf(:,:,:,ispin,io,ik,im) &
          & ,ppg,mg%is_array,mg%ie_array,ik,uVpsibox2(:,ispin,io,ik,im))
        else
          call nonlocal_part(wrk,psi%zwf(:,:,:,ispin,io,ik,im),ppg,mg%is_array,mg%ie_array,ik)
        end if
        wrk2 = wrk2 + wrk * info%occ(io,ik,ispin,im)

      end do
      end do
      call comm_summation(wrk2,wrk,3,info%icomm_ko)

      call stencil_current(wrk2,dmat%zrho_mat(:,:,:,:,:,ispin,im),stencil%coef_nab,mg%is,mg%ie,mg%ndir)
      wrk2 = wrk + wrk2

      call comm_summation(wrk2,wrk,3,info%icomm_r)

      curr(:,ispin,im) = wrk / dble(ngrid) ! ngrid = aLxyz/Hxyz
    end do
    end do

    if(info%if_divide_rspace) deallocate(uVpsibox,uVpsibox2)

    return

  contains

    subroutine kvec_part(jw,psi,kAc,is_array,ie_array,is,ie)
      implicit none
      integer   ,intent(in) :: is_array(3),ie_array(3),is(3),ie(3)
      real(8)   ,intent(in) :: kAc(3)
      complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
      real(8)               :: jw(3)
      !
      integer :: ik,io,ix,iy,iz
      real(8) :: tmp
      tmp = 0d0
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        tmp = tmp + abs(psi(ix,iy,iz))**2
      end do
      end do
      end do
      jw = kAc(:) * tmp
      return
    end subroutine kvec_part

    subroutine stencil_current(jw,zdm,nabt,is,ie,ndir)
      implicit none
      integer   ,intent(in) :: is(3),ie(3),ndir
      real(8)   ,intent(in) :: nabt(Nd,3)
      complex(8),intent(in) :: zdm(Nd,ndir,is(1)-Nd:ie(1),is(2)-Nd:ie(2),is(3)-Nd:ie(3))
      real(8)               :: jw(3)
      !
      integer :: ix,iy,iz
      complex(8) :: tmp(3)
      tmp = 0d0
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        tmp(1) = tmp(1) + nabt(1,1) * zdm(1,1,ix,iy,iz) * 2d0 &
                        + nabt(2,1) * zdm(2,1,ix,iy,iz) * 2d0 &
                        + nabt(3,1) * zdm(3,1,ix,iy,iz) * 2d0 &
                        + nabt(4,1) * zdm(4,1,ix,iy,iz) * 2d0

        tmp(2) = tmp(2) + nabt(1,2) * zdm(1,2,ix,iy,iz) * 2d0 &
                        + nabt(2,2) * zdm(2,2,ix,iy,iz) * 2d0 &
                        + nabt(3,2) * zdm(3,2,ix,iy,iz) * 2d0 &
                        + nabt(4,2) * zdm(4,2,ix,iy,iz) * 2d0

        tmp(3) = tmp(3) + nabt(1,3) * zdm(1,3,ix,iy,iz) * 2d0 &
                        + nabt(2,3) * zdm(2,3,ix,iy,iz) * 2d0 &
                        + nabt(3,3) * zdm(3,3,ix,iy,iz) * 2d0 &
                        + nabt(4,3) * zdm(4,3,ix,iy,iz) * 2d0
      end do
      end do
      end do
      jw = aimag(tmp)
      return
    end subroutine stencil_current

    subroutine nonlocal_part(jw,psi,ppg,is_array,ie_array,ik)
      implicit none
      integer   ,intent(in) :: is_array(3),ie_array(3),ik
      complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
      type(s_pp_grid),intent(in) :: ppg
      real(8)               :: jw(3)
      !
      integer    :: ilma,ia,j,i,ix,iy,iz
      real(8)    :: x,y,z
      complex(8) :: uVpsi,uVpsi_r(3)
      jw = 0d0
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
        jw = jw + 2d0* aimag(conjg(uVpsi_r)*uVpsi)
      end do
      return
    end subroutine nonlocal_part

    subroutine nonlocal_part_rdivided1(nspin,info,ppg,tpsi,uVpsibox,uVpsibox2)
      implicit none
      integer        ,intent(in) :: nspin
      type(s_orbital_parallel),intent(in) :: info
      type(s_pp_grid),intent(in) :: ppg
      type(s_orbital),intent(in) :: tpsi
      complex(8)    ,allocatable :: uVpsibox (:,:,:,:,:)
      complex(8)    ,allocatable :: uVpsibox2(:,:,:,:,:)
      !
      integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
      integer :: ilma,ia,j,ix,iy,iz,Nlma
      complex(8) :: uVpsi,wrk

      im_s = info%im_s
      im_e = info%im_e
      ik_s = info%ik_s
      ik_e = info%ik_e
      io_s = info%io_s
      io_e = info%io_e
      norb = Nspin* info%numo * info%numk * info%numm

      Nlma = ppg%Nlma

      allocate(uVpsibox (Nlma,Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e))
      allocate(uVpsibox2(Nlma,Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e))

!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz)
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin

        do ilma=1,Nlma
          ia = ppg%ia_tbl(ilma)
          uVpsi = 0.d0
          do j=1,ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
          end do
          uVpsi = uVpsi * ppg%rinv_uvu(ilma)
          uVpsibox(ilma,ispin,io,ik,im) = uVpsi
        end do

      end do
      end do
      end do
      end do
!$omp end parallel do

      call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,info%icomm_r)

      return
    end subroutine nonlocal_part_rdivided1

    subroutine nonlocal_part_rdivided2(jw,psi,ppg,is_array,ie_array,ik,uVpsibox)
      implicit none
      integer   ,intent(in) :: is_array(3),ie_array(3),ik
      complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
      type(s_pp_grid),intent(in) :: ppg
      complex(8),intent(in) :: uVpsibox(ppg%Nlma)
      real(8)               :: jw(3)
      !
      integer    :: ilma,ia,j,i,ix,iy,iz
      real(8)    :: x,y,z
      complex(8) :: uVpsi,uVpsi_r(3)
      jw = 0d0
      do ilma=1,ppg%Nlma
        ia=ppg%ia_tbl(ilma)
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
        jw = jw + 2d0* aimag(conjg(uVpsi_r)*uVpsi)
      end do
      return
    end subroutine nonlocal_part_rdivided2

  end subroutine calc_current

!===================================================================================================================================

  subroutine calc_microscopic_current(nspin,mg,stencil,info,psi,dmat,curr)
    use structures
    use salmon_communication, only: comm_summation
    implicit none
    integer,intent(in) :: nspin
    type(s_rgrid)  ,intent(in) :: mg
    type(s_stencil),intent(in) :: stencil
    type(s_orbital_parallel),intent(in) :: info
    type(s_orbital),intent(in) :: psi
    type(s_dmatrix),intent(in) :: dmat
    type(s_vector)             :: curr(nspin,info%im_s:info%im_e)
    !
    integer :: ispin,im,ik,io,is(3),ie(3),nsize
    real(8),allocatable :: wrk(:,:,:,:),wrk2(:,:,:,:)
    is = mg%is
    ie = mg%ie
    allocate(wrk(3,is(1):ie(1),is(2):ie(2),is(3):ie(3)),wrk2(3,is(1):ie(1),is(2):ie(2),is(3):ie(3)))
    nsize = 3* mg%num(1) * mg%num(2) * mg%num(3)

    do im=info%im_s,info%im_e
    do ispin=1,nspin

      wrk2 = 0d0
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e

        call kvec_part(wrk,psi%zwf(:,:,:,ispin,io,ik,im),stencil%vec_kAc(:,ik),mg%is_array,mg%ie_array,is,ie)
        wrk2 = wrk2 + wrk * info%occ(io,ik,ispin,im)

!       call nonlocal_part

      end do
      end do
      call comm_summation(wrk2,wrk,nsize,info%icomm_ko)

      call stencil_current(wrk2,dmat%zrho_mat(:,:,:,:,:,ispin,im),stencil%coef_nab,is,ie,mg%idx,mg%idy,mg%idz,mg%ndir)

      curr(ispin,im)%v = wrk + wrk2
    end do
    end do

    deallocate(wrk,wrk2)
    return

  contains

    subroutine kvec_part(jw,psi,kAc,is_array,ie_array,is,ie)
      implicit none
      integer   ,intent(in) :: is_array(3),ie_array(3),is(3),ie(3)
      real(8)   ,intent(in) :: kAc(3)
      complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
      real(8)               :: jw(3,is(1):ie(1),is(2):ie(2),is(3):ie(3))
      !
      integer :: ik,io,ix,iy,iz
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        jw(:,ix,iy,iz) = kAc(:) * abs(psi(ix,iy,iz))**2
      end do
      end do
      end do
      return
    end subroutine kvec_part

    subroutine stencil_current(jw,zdm,nabt,is,ie,idx,idy,idz,ndir)
      implicit none
      integer   ,intent(in) :: is(3),ie(3) &
                              ,idx(is(1)-Nd:ie(1)+Nd),idy(is(2)-Nd:ie(2)+Nd),idz(is(3)-Nd:ie(3)+Nd),ndir
      real(8)   ,intent(in) :: nabt(Nd,3)
      complex(8),intent(in) :: zdm(Nd,ndir,is(1)-Nd:ie(1),is(2)-Nd:ie(2),is(3)-Nd:ie(3))
      real(8)               :: jw(3,is(1):ie(1),is(2):ie(2),is(3):ie(3))
      !
      integer :: ix,iy,iz
      complex(8) :: tmp(3)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        tmp(1) = nabt(1,1) * ( zdm(1,1,ix,iy,iz) - conjg(zdm(1,1,idx(ix-1),iy,iz)) ) & ! dmat(x,-dx)==conjg(dmat(x-dx,dx))
               + nabt(2,1) * ( zdm(2,1,ix,iy,iz) - conjg(zdm(2,1,idx(ix-2),iy,iz)) ) &
               + nabt(3,1) * ( zdm(3,1,ix,iy,iz) - conjg(zdm(3,1,idx(ix-3),iy,iz)) ) &
               + nabt(4,1) * ( zdm(4,1,ix,iy,iz) - conjg(zdm(4,1,idx(ix-4),iy,iz)) )

        tmp(2) = nabt(1,2) * ( zdm(1,2,ix,iy,iz) - conjg(zdm(1,2,ix,idy(iy-1),iz)) ) &
               + nabt(2,2) * ( zdm(2,2,ix,iy,iz) - conjg(zdm(2,2,ix,idy(iy-2),iz)) ) &
               + nabt(3,2) * ( zdm(3,2,ix,iy,iz) - conjg(zdm(3,2,ix,idy(iy-3),iz)) ) &
               + nabt(4,2) * ( zdm(4,2,ix,iy,iz) - conjg(zdm(4,2,ix,idy(iy-4),iz)) )

        tmp(3) = nabt(1,3) * ( zdm(1,3,ix,iy,iz) - conjg(zdm(1,3,ix,iy,idz(iz-1))) ) &
               + nabt(2,3) * ( zdm(2,3,ix,iy,iz) - conjg(zdm(2,3,ix,iy,idz(iz-2))) ) &
               + nabt(3,3) * ( zdm(3,3,ix,iy,iz) - conjg(zdm(3,3,ix,iy,idz(iz-3))) ) &
               + nabt(4,3) * ( zdm(4,3,ix,iy,iz) - conjg(zdm(4,3,ix,iy,idz(iz-4))) )

        jw(:,ix,iy,iz) = aimag(tmp)
      end do
      end do
      end do
      return
    end subroutine stencil_current

  end subroutine calc_microscopic_current

end module
