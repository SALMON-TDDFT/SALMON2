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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

module density_matrix
  implicit none

contains

! dm(r,dr) = sum occ psi(r+dr) * conjg(psi(r))
! j(r) = sum occ aimag( conjg(psi(r))* (sum nabt*psi)(r) ) = aimag( sum_dr nabt(dr)* dm(r,dr) )

  subroutine calc_density_matrix(dmat,psi,info,rg,nspin,occ)
    use structures
    use salmon_communication, only: comm_summation
    implicit none
    integer        ,intent(in) :: nspin
    type(s_wf_info),intent(in) :: info
    real(8)        ,intent(in) :: occ(info%io_s:info%io_e,info%ik_s:info%ik_e)
    type(s_rgrid)  ,intent(in) :: rg
    type(s_wavefunction),intent(in) :: psi
    type(s_dmatrix)            :: dmat
    !
    integer :: im,ispin,ik,io,is(3),ie(3),nsize
    complex(8),allocatable :: wrk(:,:,:,:,:),wrk2(:,:,:,:,:)

    is = rg%is
    ie = rg%ie
    nsize = 9* rg%ndir * rg%num(1) * rg%num(2) * rg%num(3)

! real(rwf) & complex(zwf) ?

    allocate( wrk(-4:4,rg%ndir,is(1):ie(1),is(2):ie(2),is(3):ie(3)) &
            ,wrk2(-4:4,rg%ndir,is(1):ie(1),is(2):ie(2),is(3):ie(3)) )

    do im=info%im_s,info%im_e
    do ispin=1,nspin
      wrk = 0d0
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
        call calc_dm(wrk2,psi%zwf(:,:,:,ispin,io,ik,im),rg%is_array,rg%ie_array,is,ie,rg%idx,rg%idy,rg%idz,rg%ndir)
        wrk = wrk + wrk2 * occ(io,ik)
      end do
      end do
      call comm_summation(wrk,wrk2,nsize,info%icomm_pseudo) !??????? info%icomm_pseudo ?
      dmat%rho(:,:,:,:,:,ispin,im) = wrk2(:,:,:,:,:)
    end do
    end do

    deallocate(wrk,wrk2)
    return

  contains
    subroutine calc_dm(zdm,psi,is_array,ie_array,is,ie,idx,idy,idz,ndir)
      implicit none
      integer   ,intent(in) :: is_array(3),ie_array(3),is(3),ie(3) &
                              ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4),ndir
      complex(8),intent(in) :: psi(is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3))
      complex(8)            :: zdm(-4:4,ndir,is(1):ie(1),is(2):ie(2),is(3):ie(3))
      !
      integer :: ix,iy,iz,ii
      zdm = 0d0
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
        do ii=-4,4 ! 4=Nd
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

  subroutine density(rho,dmat,nspin,im_s,im_e)
    use structures
    implicit none
    integer,intent(in) :: nspin,im_s,im_e
    type(s_dmatrix),intent(in) :: dmat
    type(s_scalar) :: rho(nspin,im_s:im_e)
    !
    integer :: ispin,im
    do im=im_s,im_e
      do ispin=1,nspin
        rho(ispin,im)%f(:,:,:) = (1d0/3d0)* ( dmat%rho(0,1,:,:,:,ispin,im) &
                                            + dmat%rho(0,2,:,:,:,ispin,im) &
                                            + dmat%rho(0,3,:,:,:,ispin,im) )
      end do
    end do
    return
  end subroutine density

!===================================================================================================================================

  subroutine calc_current(curr,nspin,ngrid,rg,stencil,info,psi,ppg,occ,dmat)
    use structures
    use salmon_communication, only: comm_summation
    implicit none
    integer,intent(in) :: nspin,ngrid
    type(s_rgrid)  ,intent(in) :: rg
    type(s_stencil),intent(in) :: stencil
    type(s_wf_info),intent(in) :: info
    type(s_wavefunction),intent(in) :: psi
    type(s_pp_grid),intent(in) :: ppg
    real(8)        ,intent(in) :: occ(info%io_s:info%io_e,info%ik_s:info%ik_e)
    type(s_dmatrix),intent(in) :: dmat
    real(8) :: curr(3,nspin,info%im_s:info%im_e)
    !
    integer :: ispin,im,ik,io
    real(8) :: wrk(3),wrk_sum(3)

    do im=info%im_s,info%im_e
    do ispin=1,nspin
      wrk_sum = 0d0

      call stencil_current(wrk,dmat%rho(:,:,:,:,:,ispin,im),stencil%nabt,rg%is,rg%ie,rg%idx,rg%idy,rg%idz,rg%ndir)
      wrk_sum = wrk_sum + wrk

      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e

        call kvec_part(wrk,psi%zwf(:,:,:,ispin,io,ik,im),stencil%kAc(ik,:),rg%is_array,rg%ie_array,rg%is,rg%ie)
        wrk_sum = wrk_sum + wrk * occ(io,ik)

        call nonlocal_part(wrk,psi%zwf(:,:,:,ispin,io,ik,im),ppg,rg%is_array,rg%ie_array,ik)
        wrk_sum = wrk_sum + wrk * occ(io,ik)

      end do
      end do
!      call comm_summation(,,,info%icomm_pseudo) !??????? info%icomm_pseudo ?

      curr(:,ispin,im) = wrk_sum / ngrid ! ngrid = aLxyz/Hxyz
    end do
    end do
!    call comm_summation(,,,) !??????? icomm ?

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

    subroutine stencil_current(jw,zdm,nabt,is,ie,idx,idy,idz,ndir)
      implicit none
      integer   ,intent(in) :: is(3),ie(3) &
                              ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4),ndir
      real(8)   ,intent(in) :: nabt(4,3)
      complex(8),intent(in) :: zdm(-4:4,ndir,is(1):ie(1),is(2):ie(2),is(3):ie(3))
      real(8)               :: jw(3)
      !
      integer :: ix,iy,iz
      complex(8) :: tmp(3)
      tmp = 0d0
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
      !?????? yz,zx,xy ?
        tmp(1) = tmp(1) + nabt(1,1) * ( zdm(1,1,ix,iy,iz) - zdm(-1,1,ix,iy,iz) ) &
                        + nabt(2,1) * ( zdm(2,1,ix,iy,iz) - zdm(-2,1,ix,iy,iz) ) &
                        + nabt(3,1) * ( zdm(3,1,ix,iy,iz) - zdm(-3,1,ix,iy,iz) ) &
                        + nabt(4,1) * ( zdm(4,1,ix,iy,iz) - zdm(-4,1,ix,iy,iz) )

        tmp(2) = tmp(2) + nabt(1,2) * ( zdm(1,2,ix,iy,iz) - zdm(-1,2,ix,iy,iz) ) &
                        + nabt(2,2) * ( zdm(2,2,ix,iy,iz) - zdm(-2,2,ix,iy,iz) ) &
                        + nabt(3,2) * ( zdm(3,2,ix,iy,iz) - zdm(-3,2,ix,iy,iz) ) &
                        + nabt(4,2) * ( zdm(4,2,ix,iy,iz) - zdm(-4,2,ix,iy,iz) )

        tmp(3) = tmp(3) + nabt(1,3) * ( zdm(1,3,ix,iy,iz) - zdm(-1,3,ix,iy,iz) ) &
                        + nabt(2,3) * ( zdm(2,3,ix,iy,iz) - zdm(-2,3,ix,iy,iz) ) &
                        + nabt(3,3) * ( zdm(3,3,ix,iy,iz) - zdm(-3,3,ix,iy,iz) ) &
                        + nabt(4,4) * ( zdm(4,3,ix,iy,iz) - zdm(-4,3,ix,iy,iz) )
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
!          ix=ppg%Jxx(j,ia); x=Lx(i)*Hx-ix*aLx
!          iy=ppg%Jyy(j,ia); y=Ly(i)*Hy-iy*aLy
!          iz=ppg%Jzz(j,ia); z=Lz(i)*Hz-iz*aLz

          x = ppg%Rxyz(1,j,ia)
          y = ppg%Rxyz(2,j,ia)
          z = ppg%Rxyz(3,j,ia)

          ix = ppg%Jxyz(1,j,ia)
          iy = ppg%Jxyz(2,j,ia)
          iz = ppg%Jxyz(3,j,ia)

          uVpsi = uVpsi +conjg(ppg%zproj(j,ilma,ik))  *psi(ix,iy,iz)
          uVpsi_r(1) = uVpsi_r(1) + conjg(ppg%zproj(j,ilma,ik))*x*psi(ix,iy,iz)
          uVpsi_r(2) = uVpsi_r(2) + conjg(ppg%zproj(j,ilma,ik))*y*psi(ix,iy,iz)
          uVpsi_r(3) = uVpsi_r(3) + conjg(ppg%zproj(j,ilma,ik))*z*psi(ix,iy,iz)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
        jw = jw + 2d0* aimag(conjg(uVpsi_r)*uVpsi)
!        uVpsix=uVpsix*Hxyz
!        uVpsiy=uVpsiy*Hxyz
!        uVpsiz=uVpsiz*Hxyz
!        jxt=jxt+occ(ib,ik)*IaLxyz*2*aimag(conjg(uVpsix)*uVpsi)
!        jyt=jyt+occ(ib,ik)*IaLxyz*2*aimag(conjg(uVpsiy)*uVpsi)
!        jzt=jzt+occ(ib,ik)*IaLxyz*2*aimag(conjg(uVpsiz)*uVpsi)

      end do

!      jx=jx*Hxyz*IaLxyz+jxt
!      jy=jy*Hxyz*IaLxyz+jyt
!      jz=jz*Hxyz*IaLxyz+jzt

      return
    end subroutine nonlocal_part

  end subroutine calc_current

!===================================================================================================================================

  subroutine calc_microscopic_current(curr,nspin,ngrid,rg,stencil,info,psi,ppg,occ,dmat)
    use structures
    use salmon_communication, only: comm_summation
    implicit none
    integer,intent(in) :: nspin,ngrid
    type(s_rgrid)  ,intent(in) :: rg
    type(s_stencil),intent(in) :: stencil
    type(s_wf_info),intent(in) :: info
    type(s_wavefunction),intent(in) :: psi
    type(s_pp_grid),intent(in) :: ppg
    real(8)        ,intent(in) :: occ(info%io_s:info%io_e,info%ik_s:info%ik_e)
    type(s_dmatrix),intent(in) :: dmat
    type(s_vector)             :: curr(nspin,info%im_s:info%im_e)
    !
    integer :: ispin,im,ik,io,is(3),ie(3)
    real(8),allocatable :: wrk(:,:,:,:),wrk_sum(:,:,:,:)
    is = rg%is
    ie = rg%ie
    allocate(wrk(3,is(1):ie(1),is(2):ie(2),is(3):ie(3)),wrk_sum(3,is(1):ie(1),is(2):ie(2),is(3):ie(3)))

    do im=info%im_s,info%im_e
    do ispin=1,nspin
      wrk_sum = 0d0

      call stencil_current(wrk,dmat%rho(:,:,:,:,:,ispin,im),stencil%nabt,is,ie,rg%idx,rg%idy,rg%idz,rg%ndir)
      wrk_sum = wrk_sum + wrk

      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e

        call kvec_part(wrk,psi%zwf(:,:,:,ispin,io,ik,im),stencil%kAc(ik,:),rg%is_array,rg%ie_array,is,ie)
        wrk_sum = wrk_sum + wrk * occ(io,ik)

!       call nonlocal_part

      end do
      end do
!      call comm_summation(,,,info%icomm_pseudo) !??????? info%icomm_pseudo ?

      curr(ispin,im)%v = wrk_sum
    end do
    end do

    deallocate(wrk,wrk_sum)
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
                              ,idx(is(1)-4:ie(1)+4),idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4),ndir
      real(8)   ,intent(in) :: nabt(4,3)
      complex(8),intent(in) :: zdm(-4:4,ndir,is(1):ie(1),is(2):ie(2),is(3):ie(3))
      real(8)               :: jw(3,is(1):ie(1),is(2):ie(2),is(3):ie(3))
      !
      integer :: ix,iy,iz
      complex(8) :: tmp(3)
      do iz=is(3),ie(3)
      do iy=is(2),ie(2)
      do ix=is(1),ie(1)
      !?????? yz,zx,xy ?
        tmp(1) = nabt(1,1) * ( zdm(1,1,ix,iy,iz) - zdm(-1,1,ix,iy,iz) ) &
               + nabt(2,1) * ( zdm(2,1,ix,iy,iz) - zdm(-2,1,ix,iy,iz) ) &
               + nabt(3,1) * ( zdm(3,1,ix,iy,iz) - zdm(-3,1,ix,iy,iz) ) &
               + nabt(4,1) * ( zdm(4,1,ix,iy,iz) - zdm(-4,1,ix,iy,iz) )

        tmp(2) = nabt(1,2) * ( zdm(1,2,ix,iy,iz) - zdm(-1,2,ix,iy,iz) ) &
               + nabt(2,2) * ( zdm(2,2,ix,iy,iz) - zdm(-2,2,ix,iy,iz) ) &
               + nabt(3,2) * ( zdm(3,2,ix,iy,iz) - zdm(-3,2,ix,iy,iz) ) &
               + nabt(4,2) * ( zdm(4,2,ix,iy,iz) - zdm(-4,2,ix,iy,iz) )

        tmp(3) = nabt(1,3) * ( zdm(1,3,ix,iy,iz) - zdm(-1,3,ix,iy,iz) ) &
               + nabt(2,3) * ( zdm(2,3,ix,iy,iz) - zdm(-2,3,ix,iy,iz) ) &
               + nabt(3,3) * ( zdm(3,3,ix,iy,iz) - zdm(-3,3,ix,iy,iz) ) &
               + nabt(4,4) * ( zdm(4,3,ix,iy,iz) - zdm(-4,3,ix,iy,iz) )

        jw(:,ix,iy,iz) = aimag(tmp)
      end do
      end do
      end do
      return
    end subroutine stencil_current

  end subroutine calc_microscopic_current

end module
