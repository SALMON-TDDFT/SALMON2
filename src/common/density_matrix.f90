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
    use timer
    implicit none
    integer        ,intent(in) :: nspin
    type(s_orbital_parallel),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_sendrecv_grid)      :: srg
    type(s_orbital)            :: psi
    type(s_dmatrix)            :: dmat
    !
    integer :: im,ispin,ik,io,is(3),ie(3),nsize
    integer :: iz,iy,ix,ii
    complex(8) :: pocc
    complex(8),allocatable :: wrk(:,:,:,:,:)

    call timer_begin(LOG_CALC_DENSITY_MATRIX)

    is = mg%is
    ie = mg%ie
    nsize = Nd* mg%ndir * (mg%num(1)+Nd) * (mg%num(2)+Nd) * (mg%num(3)+Nd)

    allocate(wrk(Nd,mg%ndir,is(1)-Nd:ie(1),is(2)-Nd:ie(2),is(3)-Nd:ie(3)))

    if(allocated(psi%rwf)) then
      allocate(psi%zwf(mg%is_array(1):mg%ie_array(1) &
                      ,mg%is_array(2):mg%ie_array(2) &
                      ,mg%is_array(3):mg%ie_array(3) &
                      ,nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e))
      psi%zwf = cmplx(psi%rwf)
    end if

  ! overlap region communication
    if(info%if_divide_rspace) then
      call update_overlap_complex8(srg, mg, psi%zwf)
    end if

    do im=info%im_s,info%im_e
    do ispin=1,nspin
!$omp parallel private(ik,io,iz,iy,ix,ii,pocc)
!$omp do collapse(2)
      do iz=is(3)-Nd,ie(3)
      do iy=is(2)-Nd,ie(2)
      do ix=is(1)-Nd,ie(1)
        wrk(:,:,ix,iy,iz) = 0d0
      end do
      end do
      end do
!$omp end do

      ! ik and io have a data dependency.
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
!$omp do collapse(2)
      do iz=is(3)-Nd,ie(3)
      do iy=is(2)-Nd,ie(2)
      do ix=is(1)-Nd,ie(1)

        ! dir = 1,2,3 = xx,yy,zz (yz,zx,xy)
        pocc = conjg( psi%zwf(ix,iy,iz,ispin,io,ik,im) ) * info%occ(io,ik,ispin,im)
!dir$ unroll
        do ii=1,Nd
          wrk(ii,1,ix,iy,iz) = wrk(ii,1,ix,iy,iz) + psi%zwf(mg%idx(ix+ii),iy,iz,ispin,io,ik,im) * pocc
          wrk(ii,2,ix,iy,iz) = wrk(ii,2,ix,iy,iz) + psi%zwf(ix,mg%idy(iy+ii),iz,ispin,io,ik,im) * pocc
          wrk(ii,3,ix,iy,iz) = wrk(ii,3,ix,iy,iz) + psi%zwf(ix,iy,mg%idz(iz+ii),ispin,io,ik,im) * pocc
        end do

      end do
      end do
      end do
!$omp end do
      end do
      end do
!$omp end parallel
      call comm_summation(wrk(:,:,:,:,:),dmat%zrho_mat(:,:,:,:,:,ispin,im),nsize,info%icomm_ko)
    end do
    end do

    if(allocated(psi%rwf)) deallocate(psi%zwf)

    deallocate(wrk)

    call timer_end(LOG_CALC_DENSITY_MATRIX)
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

  subroutine calc_current(nspin,ngrid,mg,stencil,info,srg,psi,ppg,curr)
    use structures
    use sendrecv_grid, only: update_overlap_complex8
    use salmon_communication, only: comm_summation
    use pseudo_pt_sub, only: calc_uVpsi_rdivided
    implicit none
    integer        ,intent(in) :: nspin,ngrid
    type(s_rgrid)  ,intent(in) :: mg
    type(s_stencil),intent(in) :: stencil
    type(s_orbital_parallel),intent(in) :: info
    type(s_sendrecv_grid)      :: srg
    type(s_orbital)            :: psi
    type(s_pp_grid),intent(in) :: ppg
    real(8) :: curr(3,nspin,info%im_s:info%im_e)
    !
    integer :: ispin,im,ik,io
    real(8),dimension(3) :: wrk1,wrk2,wrk3,wrk4
    real(8) :: BT(3,3)
    complex(8),allocatable :: uVpsibox (:,:,:,:,:)
    complex(8),allocatable :: uVpsibox2(:,:,:,:,:)

    BT = transpose(stencil%rmatrix_B)
    if(info%if_divide_rspace) call calc_uVpsi_rdivided(nspin,info,ppg,psi,uVpsibox,uVpsibox2)

  ! overlap region communication
    if(info%if_divide_rspace) then
      call update_overlap_complex8(srg, mg, psi%zwf)
    end if

    do im=info%im_s,info%im_e
    do ispin=1,nspin

      wrk4 = 0d0
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e

        call stencil_current(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz,stencil%coef_nab &
                            ,stencil%vec_kAc(:,ik),psi%zwf(:,:,:,ispin,io,ik,im),wrk1,wrk2)
        wrk2 = matmul(BT,wrk2)

        if(info%if_divide_rspace) then
          call calc_current_nonlocal_rdivided(wrk3,psi%zwf(:,:,:,ispin,io,ik,im),ppg,mg%is_array,mg%ie_array,ik &
                                             ,uVpsibox2(:,ispin,io,ik,im))
        else
          call calc_current_nonlocal         (wrk3,psi%zwf(:,:,:,ispin,io,ik,im),ppg,mg%is_array,mg%ie_array,ik)
        end if

        wrk4 = wrk4 + (wrk1 + wrk2 + wrk3) * info%occ(io,ik,ispin,im)

      end do
      end do

      call comm_summation(wrk4,wrk1,3,info%icomm_rko)

      curr(:,ispin,im) = wrk1 / dble(ngrid) ! ngrid = aLxyz/Hxyz
    end do
    end do

    if(info%if_divide_rspace) deallocate(uVpsibox,uVpsibox2)

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
      complex(8) :: cpsi,tmp(3)
      rtmp = 0d0
      tmp = 0d0
!$omp parallel do collapse(2) private(iz,iy,ix,cpsi) reduction(+:rtmp,tmp)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        cpsi = conjg(psi(ix,iy,iz))
        rtmp = rtmp + abs(psi(ix,iy,iz))**2

        tmp(1) = tmp(1) + nabt(1,1) * cpsi * psi(idx(ix+1),iy,iz) &
                        + nabt(2,1) * cpsi * psi(idx(ix+2),iy,iz) &
                        + nabt(3,1) * cpsi * psi(idx(ix+3),iy,iz) &
                        + nabt(4,1) * cpsi * psi(idx(ix+4),iy,iz)

        tmp(2) = tmp(2) + nabt(1,2) * cpsi * psi(ix,idy(iy+1),iz) &
                        + nabt(2,2) * cpsi * psi(ix,idy(iy+2),iz) &
                        + nabt(3,2) * cpsi * psi(ix,idy(iy+3),iz) &
                        + nabt(4,2) * cpsi * psi(ix,idy(iy+4),iz)

        tmp(3) = tmp(3) + nabt(1,3) * cpsi * psi(ix,iy,idz(iz+1)) &
                        + nabt(2,3) * cpsi * psi(ix,iy,idz(iz+2)) &
                        + nabt(3,3) * cpsi * psi(ix,iy,idz(iz+3)) &
                        + nabt(4,3) * cpsi * psi(ix,iy,idz(iz+4))
      end do
      end do
      end do
!$omp end parallel do
      j1 = kAc(:) * rtmp
      j2 = aimag(tmp * 2d0)
      return
    end subroutine stencil_current

  end subroutine calc_current

  subroutine calc_current_use_dmat(nspin,ngrid,mg,stencil,info,psi,ppg,dmat,curr)
    use structures
    use salmon_communication, only: comm_summation
    use pseudo_pt_sub, only: calc_uVpsi_rdivided
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
    real(8) :: wrk1(3),wrk2(3),wrk3(3),BT(3,3)
    complex(8),allocatable :: uVpsibox (:,:,:,:,:)
    complex(8),allocatable :: uVpsibox2(:,:,:,:,:)

    BT = transpose(stencil%rmatrix_B)
    if(info%if_divide_rspace) call calc_uVpsi_rdivided(nspin,info,ppg,psi,uVpsibox,uVpsibox2)

    do im=info%im_s,info%im_e
    do ispin=1,nspin

      wrk3 = 0d0
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e

        call kvec_part(wrk1,psi%zwf(:,:,:,ispin,io,ik,im),stencil%vec_kAc(:,ik),mg%is_array,mg%ie_array,mg%is,mg%ie)

        if(info%if_divide_rspace) then
          call calc_current_nonlocal_rdivided(wrk2,psi%zwf(:,:,:,ispin,io,ik,im),ppg,mg%is_array,mg%ie_array,ik &
                                             ,uVpsibox2(:,ispin,io,ik,im))
        else
          call calc_current_nonlocal         (wrk2,psi%zwf(:,:,:,ispin,io,ik,im),ppg,mg%is_array,mg%ie_array,ik)
        end if

        wrk3 = wrk3 + (wrk1 + wrk2) * info%occ(io,ik,ispin,im)

      end do
      end do

      call stencil_current(wrk2,dmat%zrho_mat(:,:,:,:,:,ispin,im),stencil%coef_nab,mg%is,mg%ie,mg%ndir)

      call comm_summation(wrk3,wrk1,3,info%icomm_ko)

      wrk2 = wrk1 + matmul(BT,wrk2)
      call comm_summation(wrk2,wrk1,3,info%icomm_r)

      curr(:,ispin,im) = wrk1 / dble(ngrid) ! ngrid = aLxyz/Hxyz
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
!$omp parallel do collapse(2) private(iz,iy,ix) reduction(+:tmp)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        tmp = tmp + abs(psi(ix,iy,iz))**2
      end do
      end do
      end do
!$omp end parallel do
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
!$omp parallel do collapse(2) private(iz,iy,ix) reduction(+:tmp)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        tmp(1) = tmp(1) + nabt(1,1) * zdm(1,1,ix,iy,iz) &
                        + nabt(2,1) * zdm(2,1,ix,iy,iz) &
                        + nabt(3,1) * zdm(3,1,ix,iy,iz) &
                        + nabt(4,1) * zdm(4,1,ix,iy,iz)

        tmp(2) = tmp(2) + nabt(1,2) * zdm(1,2,ix,iy,iz) &
                        + nabt(2,2) * zdm(2,2,ix,iy,iz) &
                        + nabt(3,2) * zdm(3,2,ix,iy,iz) &
                        + nabt(4,2) * zdm(4,2,ix,iy,iz)

        tmp(3) = tmp(3) + nabt(1,3) * zdm(1,3,ix,iy,iz) &
                        + nabt(2,3) * zdm(2,3,ix,iy,iz) &
                        + nabt(3,3) * zdm(3,3,ix,iy,iz) &
                        + nabt(4,3) * zdm(4,3,ix,iy,iz)
      end do
      end do
      end do
!$omp end parallel do
      jw = aimag(tmp * 2d0)
      return
    end subroutine stencil_current

  end subroutine calc_current_use_dmat

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
    integer    :: ilma,ia,j,ix,iy,iz
    real(8)    :: x,y,z
    complex(8) :: uVpsi,uVpsi_r(3)
    jw = 0d0
!$omp parallel do private(ilma,ia,uVpsi,uVpsi_r,j,x,y,z,ix,iy,iz) reduction(+:jw)
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
      jw = jw + aimag(conjg(uVpsi_r)*uVpsi)
    end do
!$omp end parallel do
    jw = jw * 2d0
    return
  end subroutine calc_current_nonlocal_rdivided

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
