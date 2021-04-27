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
module write_sub
  use math_constants,only : zi,pi
  implicit none

contains

!===================================================================================================================================

  !! export SYSNAME_k.data file
  subroutine write_k_data(system,stencil)
    use structures
    use salmon_global, only: sysname
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root,comm_sync_all
    use filesystem, only: open_filehandle
    implicit none
    type(s_dft_system) ,intent(in) :: system
    type(s_stencil),intent(in) :: stencil
    !
    integer :: fh_k
    integer :: ik,NK
    character(100) :: file_k_data

    NK = system%nk
    file_k_data = trim(sysname)//'_k.data'

    if (comm_is_root(nproc_id_global)) then
      fh_k = open_filehandle(file_k_data, status="replace")
      if(stencil%if_orthogonal) then
        write(fh_k, '("#",1X,A)') "k-point distribution"
        write(fh_k, '("#",1X,A,":",1X,A)') "ik", "k-point index"
        write(fh_k, '("#",1X,A,":",1X,A)') "kx,ky,kz", "Reduced coordinate of k-points"
        write(fh_k, '("#",1X,A,":",1X,A)') "wk", "Weight of k-point"
        write(fh_k, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "ik", "none", &
          & 2, "kx", "none", &
          & 3, "ky", "none", &
          & 4, "kz", "none", &
          & 5, "wk", "none"
        do ik = 1, NK
          write(fh_k, '(I6,99(1X,E23.15E3))') &
            & ik, &
            & system%vec_k(1,ik) / system%primitive_b(1,1), &
            & system%vec_k(2,ik) / system%primitive_b(2,2), &
            & system%vec_k(3,ik) / system%primitive_b(3,3), &
            & system%wtk(ik)
        end do !ik
        write(fh_k, '("#",1X,A)') "coefficients (2*pi/a [a.u.]) in kx, ky, kz"
        write(fh_k, '(3E23.15E3)') system%primitive_b(1,1), system%primitive_b(2,2), system%primitive_b(3,3) 

      else
        write(fh_k, '("#",1X,A)') "k-point distribution (nonorthogonal coordinate)"
        write(fh_k, '("#",1X,A,":",1X,A)') "brl", "reciprocal primitive vectors"
        write(fh_k, '("#",1X,A,":",1X,A)') "ik", "k-point index"
        write(fh_k, '("#",1X,A,":",1X,A)') "kx,ky,kz", "k-vectors"
        write(fh_k, '("#",1X,A,":",1X,A)') "wk", "Weight of k-point"
        write(fh_k, '("#",A,"[",A,"]")') "brl(1:3,1:3)", "a.u."
        write(fh_k, '(9(1X,E23.15E3))') system%primitive_b(1:3,1:3)
        write(fh_k, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "ik", "none", &
          & 2, "kx", "a.u.", &
          & 3, "ky", "a.u.", &
          & 4, "kz", "a.u.", &
          & 5, "wk", "none"
        do ik = 1, NK
          write(fh_k, '(I6,99(1X,E23.15E3))') &
            & ik, &
            & system%vec_k(1,ik), &
            & system%vec_k(2,ik), &
            & system%vec_k(3,ik), &
            & system%wtk(ik)
        end do !ik
      end if
      close(fh_k)
    end if
    call comm_sync_all
  end subroutine write_k_data

!===================================================================================================================================

  !! export SYSNAME_tm.data file
  subroutine write_tm_data(tpsi,system,info,mg,stencil,srg,ppg)
    use structures
    use stencil_sub
    use sendrecv_grid
    use salmon_global, only: sysname
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root,comm_summation,comm_sync_all
    use filesystem, only: open_filehandle
    implicit none
    type(s_dft_system) ,intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_pp_grid),intent(in) :: ppg
    type(s_orbital)       :: tpsi
    !
    integer :: fh_tm, narray
    integer :: i,j,ik,ib,ib1,ib2,ilma,nlma,ia,ix,iy,iz,NB,NK,im,ispin,ik_s,ik_e,is(3),ie(3)
    real(8) :: x,y,z
    complex(8),allocatable :: upu(:,:,:,:),upu_l(:,:,:,:)
    complex(8),allocatable :: gtpsi(:,:,:,:),gtpsi_l(:,:,:,:)
    complex(8),allocatable :: uVpsi(:,:),uVpsix(:,:),uVpsiy(:,:),uVpsiz(:,:)
    complex(8),allocatable :: uVpsixx(:,:),uVpsixy(:,:),uVpsixz(:,:)
    complex(8),allocatable :: uVpsiyy(:,:),uVpsiyz(:,:),uVpsizz(:,:)
    complex(8),allocatable :: u_rVnl_Vnlr_u(:,:,:,:),u_rVnl_Vnlr_u_l(:,:,:,:)
    complex(8) :: u_rVnl_u(3),u_Vnlr_u(3),veik
    complex(8),allocatable ::  u_rVnlr_Vnlrr_u(:,:,:,:),u_rVnlr_Vnlrr_u_l(:,:,:,:)
!    complex(8) :: ctmp1,ctmp2
    complex(8) :: wrk(3)
    character(100) :: file_tm_data

    if(info%im_s/=1 .or. info%im_e/=1) then!??????
      write(*,*) "error @ write_tm_data: im/=1"
      return
    endif
    if(system%Nspin/=1) then!??????
      write(*,*) "error @ write_tm_data: nspin/=1"
      return
    endif
    if(info%io_s/=1 .or. info%io_e/=system%no) then!??????
      if (comm_is_root(nproc_id_global)) then
        write(*,*) "error @ write_tm_data: do not use orbital parallelization"
        write(*,*) "Only <u|p|u> is printed (pseudopotential terms are not printed)"
      endif
     !return
    endif
    if(.not. allocated(tpsi%zwf)) then!??????
      write(*,*) "error @ write_tm_data: do not use real wavefunction (iperiodic=0)"
      return
    endif

    im = 1
    ispin = 1

    NB = system%no
    NK = system%nk

    ik_s = info%ik_s
    ik_e = info%ik_e

    is = mg%is
    ie = mg%ie

    Nlma = ppg%Nlma

    allocate(upu_l(3,NB,NB,NK),upu(3,NB,NB,NK),uVpsi(NB,NK),uVpsix(NB,NK),uVpsiy(NB,NK),uVpsiz(NB,NK), &
    uVpsixx(NB,NK),uVpsixy(NB,NK),uVpsixz(NB,NK),uVpsiyy(NB,NK),uVpsiyz(NB,NK),uVpsizz(NB,NK), &
    u_rVnl_Vnlr_u(3,NB,NB,NK),u_rVnl_Vnlr_u_l(3,NB,NB,NK),u_rVnlr_Vnlrr_u(3,3,NB,NK),u_rVnlr_Vnlrr_u_l(3,3,NB,NK))

    !calculate <u_nk|p_j|u_mk>  (j=x,y,z)

    upu_l(:,:,:,:) = 0d0

    allocate(gtpsi_l(3,mg%is_array(1):mg%ie_array(1) &
                      ,mg%is_array(2):mg%ie_array(2) &
                      ,mg%is_array(3):mg%ie_array(3)))
    allocate(gtpsi(3,mg%is_array(1):mg%ie_array(1) &
                    ,mg%is_array(2):mg%ie_array(2) &
                    ,mg%is_array(3):mg%ie_array(3)))
    narray= 3 * ( mg%ie_array(1) - mg%is_array(1) + 1 ) &
              * ( mg%ie_array(2) - mg%is_array(2) + 1 ) &
              * ( mg%ie_array(3) - mg%is_array(3) + 1 )

  ! overlap region communication
    if(info%if_divide_rspace) then
      call update_overlap_complex8(srg, mg, tpsi%zwf)
    end if

    do ik=ik_s,ik_e
    do ib2=1,NB

      ! gtpsi = (nabla) psi
      gtpsi_l(:,:,:,:) = 0d0
      if(ib2.ge.info%io_s .and. ib2.le.info%io_e) &
        call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,ib2,ik,im),gtpsi_l,mg%is_array,mg%ie_array,mg%is,mg%ie &
            ,mg%idx,mg%idy,mg%idz,stencil%coef_nab,system%rmatrix_B)
      call comm_summation(gtpsi_l,gtpsi,narray,info%icomm_o)

!     do ib1=1,NB
      do ib1=info%io_s,info%io_e
        do i=1,3
          wrk(i) = sum(conjg(tpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,ib1,ik,im)) &
                            * gtpsi(i,is(1):ie(1),is(2):ie(2),is(3):ie(3)) )
        end do
        upu_l(:,ib1,ib2,ik) = - zI * wrk * system%Hvol
      end do
    end do
    end do
    call comm_summation(upu_l,upu,3*NB*NB*NK,info%icomm_rko)
    deallocate(gtpsi)

    !--(Print only for orbital parallelization: just temporal)
    if(info%io_s/=1 .or. info%io_e/=system%no) then

       file_tm_data = trim(sysname)//'_tm.data'
       if (comm_is_root(nproc_id_global)) then
          fh_tm = open_filehandle(file_tm_data, status="replace")
          write(fh_tm, '("#",1X,A)') "#Transition Moment between occupied and unocupied orbitals in GS"
          write(fh_tm, '("#",1X,A)') "# (Separated analysis tool is available)"

          !<u_nk|p_j|u_mk>  (j=x,y,z)
          write(fh_tm,*) "#<u_nk|p_j|u_mk>  (j=x,y,z)"
          do ik =1,NK
          do ib1=1,NB
          do ib2=1,NB
             write(fh_tm,9000) ik,ib1,ib2,(upu(j,ib1,ib2,ik),j=1,3)
          enddo
          enddo
          enddo
          close(fh_tm)
       endif
       return
    endif
    !-------------------------------------------------

    !calculate <u_mk|[r_j,dVnl^(0)]|u_nk>  (j=x,y,z)
    u_rVnl_Vnlr_u_l = 0d0
    do ik=ik_s,ik_e
    do ilma=1,Nlma
       ia=ppg%ia_tbl(ilma)
       uVpsi=0d0;  uVpsix=0d0;  uVpsiy=0d0;  uVpsiz=0d0
       do j=1,ppg%Mps(ia)
          x = ppg%Rxyz(1,j,ia)
          y = ppg%Rxyz(2,j,ia)
          z = ppg%Rxyz(3,j,ia)
          ix = ppg%Jxyz(1,j,ia)
          iy = ppg%Jxyz(2,j,ia)
          iz = ppg%Jxyz(3,j,ia)
          veik = conjg(ppg%zekr_uV(j,ilma,ik))
          do ib=1,NB
          uVpsi( ib,ik) =uVpsi( ib,ik)+ veik*    tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik|u>
          uVpsix(ib,ik) =uVpsix(ib,ik)+ veik* x *tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*x|u>
          uVpsiy(ib,ik) =uVpsiy(ib,ik)+ veik* y *tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*y|u>
          uVpsiz(ib,ik) =uVpsiz(ib,ik)+ veik* z *tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*z|u>
          enddo
       end do
       uVpsi  = uVpsi * ppg%rinv_uvu(ilma)

!$omp parallel
!$omp do private(ib1,ib2,u_rVnl_u,u_Vnlr_u) collapse(2)
       do ib1=1,NB
       do ib2=1,NB
          !<u|e^-ik*r|v><v|e^ik|u>
          u_rVnl_u(1)= conjg(uVpsix(ib1,ik))*uVpsi(ib2,ik)
          u_rVnl_u(2)= conjg(uVpsiy(ib1,ik))*uVpsi(ib2,ik)
          u_rVnl_u(3)= conjg(uVpsiz(ib1,ik))*uVpsi(ib2,ik)
          !<u|e^-ik|v><v|e^ik*r|u>
          u_Vnlr_u(1)= conjg(uVpsi(ib1,ik))*uVpsix(ib2,ik)
          u_Vnlr_u(2)= conjg(uVpsi(ib1,ik))*uVpsiy(ib2,ik)
          u_Vnlr_u(3)= conjg(uVpsi(ib1,ik))*uVpsiz(ib2,ik)

          u_rVnl_Vnlr_u_l(:,ib1,ib2,ik) = u_rVnl_Vnlr_u_l(:,ib1,ib2,ik)  &
          &                               + u_rVnl_u(:) - u_Vnlr_u(:)
       enddo
       enddo
!$omp end do
!$omp end parallel
    enddo  !ilma
    enddo  !ik
    call comm_summation(u_rVnl_Vnlr_u_l,u_rVnl_Vnlr_u,3*NB*NB*NK,info%icomm_rko)


!    !calculate <u_nk|[r_j,dVnl^(0)]r|u_nk>  (j=x,y,z)
!    u_rVnlr_Vnlrr_u_l(:,:,:,:) = 0d0
!    do ik=ik_s,ik_e
!    do ilma=1,Nlma
!       ia=ppg%ia_tbl(ilma)
!       uVpsi=0d0;  uVpsix=0d0;  uVpsiy=0d0;  uVpsiz=0d0
!       uVpsixx=0d0;  uVpsixy=0d0;  uVpsixz=0d0
!                     uVpsiyy=0d0;  uVpsiyz=0d0
!                                   uVpsizz=0d0
!       do j=1,ppg%Mps(ia)
!          x = ppg%Rxyz(1,j,ia)
!          y = ppg%Rxyz(2,j,ia)
!          z = ppg%Rxyz(3,j,ia)
!          ix = ppg%Jxyz(1,j,ia)
!          iy = ppg%Jxyz(2,j,ia)
!          iz = ppg%Jxyz(3,j,ia)
!          veik = conjg(ppg%zekr_uV(j,ilma,ik))
!          do ib=1,NB
!          uVpsi(  ib,ik)=uVpsi(  ib,ik)+veik*    tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik|u>
!          uVpsix( ib,ik)=uVpsix( ib,ik)+veik* x *tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*x|u>
!          uVpsiy( ib,ik)=uVpsiy( ib,ik)+veik* y *tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*y|u>
!          uVpsiz( ib,ik)=uVpsiz( ib,ik)+veik* z *tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*z|u>
!          uVpsixx(ib,ik)=uVpsixx(ib,ik)+veik*x*x*tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*xx|u>
!          uVpsixy(ib,ik)=uVpsixy(ib,ik)+veik*x*y*tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*xy|u>
!          uVpsixz(ib,ik)=uVpsixz(ib,ik)+veik*x*z*tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*xz|u>
!          uVpsiyy(ib,ik)=uVpsiyy(ib,ik)+veik*y*y*tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*yy|u>
!          uVpsiyz(ib,ik)=uVpsiyz(ib,ik)+veik*y*z*tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*yz|u>
!          uVpsizz(ib,ik)=uVpsizz(ib,ik)+veik*z*z*tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*zz|u>
!          enddo
!
!       end do
!
!       do ib=1,NB
!          !xx
!          ctmp1 = conjg(uVpsix(ib,ik))*uVpsix( ib,ik)
!          ctmp2 = conjg(uVpsi( ib,ik))*uVpsixx(ib,ik)
!          u_rVnlr_Vnlrr_u_l(1,1,ib,ik) = &
!          u_rVnlr_Vnlrr_u_l(1,1,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
!          !xy
!          ctmp1 = conjg(uVpsix(ib,ik))*uVpsiy( ib,ik)
!          ctmp2 = conjg(uVpsi( ib,ik))*uVpsixy(ib,ik)
!          u_rVnlr_Vnlrr_u_l(1,2,ib,ik) = &
!          u_rVnlr_Vnlrr_u_l(1,2,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
!          !xz
!          ctmp1 = conjg(uVpsix(ib,ik))*uVpsiz( ib,ik)
!          ctmp2 = conjg(uVpsi( ib,ik))*uVpsixz(ib,ik)
!          u_rVnlr_Vnlrr_u_l(1,3,ib,ik) = &
!          u_rVnlr_Vnlrr_u_l(1,3,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
!          !yx
!          ctmp1 = conjg(uVpsiy(ib,ik))*uVpsix( ib,ik)
!          ctmp2 = conjg(uVpsi( ib,ik))*uVpsixy(ib,ik)
!          u_rVnlr_Vnlrr_u_l(2,1,ib,ik) = &
!          u_rVnlr_Vnlrr_u_l(2,1,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
!          !yy
!          ctmp1 = conjg(uVpsiy(ib,ik))*uVpsiy( ib,ik)
!          ctmp2 = conjg(uVpsi( ib,ik))*uVpsiyy(ib,ik)
!          u_rVnlr_Vnlrr_u_l(2,2,ib,ik) = &
!          u_rVnlr_Vnlrr_u_l(2,2,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
!          !yz
!          ctmp1 = conjg(uVpsiy(ib,ik))*uVpsiz( ib,ik)
!          ctmp2 = conjg(uVpsi( ib,ik))*uVpsiyz(ib,ik)
!          u_rVnlr_Vnlrr_u_l(2,3,ib,ik) = &
!          u_rVnlr_Vnlrr_u_l(2,3,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
!          !zx
!          ctmp1 = conjg(uVpsiz(ib,ik))*uVpsix( ib,ik)
!          ctmp2 = conjg(uVpsi( ib,ik))*uVpsixz(ib,ik)
!          u_rVnlr_Vnlrr_u_l(3,1,ib,ik) = &
!          u_rVnlr_Vnlrr_u_l(3,1,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
!          !zy
!          ctmp1 = conjg(uVpsiz(ib,ik))*uVpsiy( ib,ik)
!          ctmp2 = conjg(uVpsi( ib,ik))*uVpsiyz(ib,ik)
!          u_rVnlr_Vnlrr_u_l(3,2,ib,ik) = &
!          u_rVnlr_Vnlrr_u_l(3,2,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
!          !zz
!          ctmp1 = conjg(uVpsiz(ib,ik))*uVpsiz( ib,ik)
!          ctmp2 = conjg(uVpsi( ib,ik))*uVpsizz(ib,ik)
!          u_rVnlr_Vnlrr_u_l(3,3,ib,ik) = &
!          u_rVnlr_Vnlrr_u_l(3,3,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
!
!       enddo
!
!    enddo  !ilma
!    enddo  !ik
!    call comm_summation(u_rVnlr_Vnlrr_u_l,u_rVnlr_Vnlrr_u,3*3*NB*NK,info%icomm_rko)

    file_tm_data = trim(sysname)//'_tm.data'!??????

    if (comm_is_root(nproc_id_global)) then
      fh_tm = open_filehandle(file_tm_data, status="replace")
      write(fh_tm, '("#",1X,A)') "#Transition Moment between occupied and unocupied orbitals in GS"
      write(fh_tm, '("#",1X,A)') "# (Separated analysis tool is available)"


      !Currently, TEST: print format is not decided

       !<u_nk|p_j|u_mk>  (j=x,y,z)
       write(fh_tm,*) "#<u_nk|p_j|u_mk>  (j=x,y,z)"
       do ik =1,NK
       do ib1=1,NB
       do ib2=1,NB
          write(fh_tm,9000) ik,ib1,ib2,(upu(j,ib1,ib2,ik),j=1,3)
       enddo
       enddo
       enddo
!9000     format(3i8,6e18.10)
9000     format(3i8,6e18.5)

!       !<u_mk|[r_j,dVnl^(0)]r_i|u_nk>  (j,i=x,y,z)
!       write(fh_tm,*) "#<u_mk|[r_j,dVnl^(0)]r_i|u_nk>  (j,i=x,y,z)"
!       do ik=1,NK
!       do ib=1,NB
!          do i=1,3
!             write(fh_tm,9000) ik,ib,i,(u_rVnlr_Vnlrr_u(i,j,ib,ik),j=1,3)
!          enddo
!       enddo
!       enddo

       !<u_mk|[r_j,dVnl^(0)]|u_nk>  (j=x,y,z)
       write(fh_tm,*) "#<u_mk|[r_j,dVnl^(0)]|u_nk>  (j=x,y,z)"
       do ik =1,NK
       do ib1=1,NB
       do ib2=1,NB
          write(fh_tm,9000) ik,ib1,ib2,(u_rVnl_Vnlr_u(j,ib1,ib2,ik),j=1,3)
       enddo
       enddo
       enddo

    end if

    if (comm_is_root(nproc_id_global)) then
      close(fh_tm)
    end if
    call comm_sync_all
    return
  end subroutine write_tm_data

!===================================================================================================================================

  subroutine write_xyz(comment,action,rvf,system)
  ! Write xyz in xyz format but also velocity and force are printed if necessary
  ! (these can be used for restart of opt and md)
    use structures, only: s_dft_system
    use inputoutput, only: au_length_aa
    use salmon_global, only: SYSname,atom_name
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    implicit none

    type(s_dft_system),intent(in) :: system

    integer :: ia,unit_xyz=200
    character(3) :: action,rvf
    character(1024) :: file_trj
    character(*) :: comment

    if(.not. comm_is_root(nproc_id_global)) return

    if(action=='new') then

       file_trj=trim(SYSname)//'_trj.xyz'
       open(unit_xyz,file=trim(file_trj),status="unknown")

    else if(action=='add') then

       write(unit_xyz,*) system%nion
       write(unit_xyz,*) trim(comment)
       do ia=1,system%nion
          if(      rvf=="r  " ) then
             write(unit_xyz,100) trim(atom_name(ia)),system%Rion(1:3,ia)*au_length_aa
          else if( rvf=="rv " ) then
             write(unit_xyz,110) trim(atom_name(ia)),system%Rion(1:3,ia)*au_length_aa,system%Velocity(1:3,ia)
          else if( rvf=="rvf" ) then
             write(unit_xyz,120) trim(atom_name(ia)),system%Rion(1:3,ia)*au_length_aa,system%Velocity(1:3,ia),system%Force(1:3,ia)
          endif
       enddo

    else if(action=='end') then
       close(unit_xyz)
    endif

100 format(a2,3f18.10)
110 format(a2,3f18.10, "  #v=",3f18.10)
120 format(a2,3f18.10, "  #v=",3f18.10, "  #f=",3f18.10)

  end subroutine write_xyz
  
!===================================================================================================================================

  subroutine write_rt_data_0d(it,ofl,dt,system,rt)
    use structures, only: s_ofile, s_dft_system, s_rt
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use filesystem, only: open_filehandle
    use inputoutput, only: t_unit_length,t_unit_time,t_unit_ac,t_unit_elec
    implicit none
    integer, intent(in) :: it
    type(s_ofile) :: ofl
    type(s_dft_system), intent(in) :: system
    type(s_rt),intent(in) :: rt
    integer :: uid, jt
    real(8) :: dm_e(3,system%nspin),ddm_e(3,system%nspin), dm(3), dt

    jt = it
    if(it.lt.0) jt=0
    ddm_e(:,1)= rt%dDp_e(:,jt)
    dm_e(:,1) = rt%Dp_e(:,jt)
    dm(:)     = rt%Dp_e(:,jt) + rt%Dp_i(:,jt) 


    if (comm_is_root(nproc_id_global)) then

    if(it.lt.0) then  !print header
       ofl%fh_rt        = open_filehandle(ofl%file_rt_data)
       uid = ofl%fh_rt

10     format("#",1X,A,":",1X,A)
       write(uid,10) "Real time calculation",""
       write(uid,10) "Ac_ext", "External vector potential field"
       write(uid,10) "E_ext", "External electric field"
       write(uid,10) "Ac_tot", "Total vector potential field"
       write(uid,10) "E_tot", "Total electric field"
       write(uid,10) "ddm_e", "Change of dipole moment (electrons/plus definition)"
       write(uid,10) "dm",    "Total dipole moment (electrons/minus + ions/plus)"

       write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
         & 1, "Time", trim(t_unit_time%name), &
         & 2, "Ac_ext_x", trim(t_unit_ac%name), &
         & 3, "Ac_ext_y", trim(t_unit_ac%name), &
         & 4, "Ac_ext_z", trim(t_unit_ac%name), &
         & 5, "E_ext_x", trim(t_unit_elec%name), &
         & 6, "E_ext_y", trim(t_unit_elec%name), &
         & 7, "E_ext_z", trim(t_unit_elec%name), &
         & 8, "Ac_tot_x", trim(t_unit_ac%name), &
         & 9, "Ac_tot_y", trim(t_unit_ac%name), &
         & 10, "Ac_tot_z", trim(t_unit_ac%name), &
         & 11, "E_tot_x", trim(t_unit_elec%name), &
         & 12, "E_tot_y", trim(t_unit_elec%name), &
         & 13, "E_tot_z", trim(t_unit_elec%name), &
         & 14, "ddm_e_x", trim(t_unit_length%name), &
         & 15, "ddm_e_y", trim(t_unit_length%name), &
         & 16, "ddm_e_z", trim(t_unit_length%name), &
         & 17, "dm_x", trim(t_unit_length%name), &
         & 18, "dm_y", trim(t_unit_length%name), &
         & 19, "dm_z", trim(t_unit_length%name)
       write(uid,*)
       flush(uid)

    else  !it>=0
       uid = ofl%fh_rt
       write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
          & it * dt * t_unit_time%conv,    &
          & system%vec_Ac_ext(1:3) * t_unit_ac%conv, &
          & system%vec_E_ext(1:3) * t_unit_elec%conv, &
          & system%vec_Ac(1:3) * t_unit_ac%conv, &
          & system%vec_E(1:3) * t_unit_elec%conv, &
          & ddm_e(1:3,1) * t_unit_length%conv, &
          & dm(1:3)      * t_unit_length%conv
       write(uid,*)
    endif
    endif

  end subroutine
!===================================================================================================================================
  subroutine write_rt_data_3d(it,ofl,dt,system,curr_e,curr_i)
    use structures, only: s_ofile, s_dft_system
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use filesystem, only: open_filehandle
    use inputoutput, only: yn_md,t_unit_time,t_unit_current,t_unit_ac,t_unit_elec
    implicit none
    type(s_ofile) :: ofl
    integer, intent(in) :: it
    type(s_dft_system), intent(in) :: system
    real(8),intent(in) :: curr_e(3,2), curr_i(3)
    integer :: uid
    real(8) :: dt

    if (comm_is_root(nproc_id_global)) then

    if(it.lt.0) then  !print header
       ofl%fh_rt        = open_filehandle(ofl%file_rt_data)
       uid = ofl%fh_rt

10     format("#",1X,A,":",1X,A)
       write(uid,10) "Real time calculation",""
       write(uid,10) "Ac_ext", "External vector potential field"
       write(uid,10) "E_ext", "External electric field"
       write(uid,10) "Ac_tot", "Total vector potential field"
       write(uid,10) "E_tot", "Total electric field"
       if(system%nspin==1) then
         write(uid,10) "Jm", "Matter current density (electrons)"
       else if(system%nspin==2) then
         write(uid,10) "Jm_u", "Matter current density for spin-up electrons"
         write(uid,10) "Jm_d", "Matter current density for spin-down electrons"
       end if
       if(yn_md=='y') then
         write(uid,10) "Jmi","Matter current density (ions)"
       endif
       write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
         & 1, "Time", trim(t_unit_time%name), &
         & 2, "Ac_ext_x", trim(t_unit_ac%name), &
         & 3, "Ac_ext_y", trim(t_unit_ac%name), &
         & 4, "Ac_ext_z", trim(t_unit_ac%name), &
         & 5, "E_ext_x", trim(t_unit_elec%name), &
         & 6, "E_ext_y", trim(t_unit_elec%name), &
         & 7, "E_ext_z", trim(t_unit_elec%name), &
         & 8, "Ac_tot_x", trim(t_unit_ac%name), &
         & 9, "Ac_tot_y", trim(t_unit_ac%name), &
         & 10, "Ac_tot_z", trim(t_unit_ac%name), &
         & 11, "E_tot_x", trim(t_unit_elec%name), &
         & 12, "E_tot_y", trim(t_unit_elec%name), &
         & 13, "E_tot_z", trim(t_unit_elec%name)
       if(system%nspin==1) then
         write(uid, '(99(1X,I0,":",A,"[",A,"]"))',advance='no') &
         & 14, "Jm_x", trim(t_unit_current%name), &
         & 15, "Jm_y", trim(t_unit_current%name), &
         & 16, "Jm_z", trim(t_unit_current%name)
         if(yn_md=='y') then
            write(uid, '(99(1X,I0,":",A,"[",A,"]"))',advance='no') &
                 & 17, "Jmi_x", trim(t_unit_current%name), &
                 & 18, "Jmi_y", trim(t_unit_current%name), &
                 & 19, "Jmi_z", trim(t_unit_current%name)
         endif
       else if(system%nspin==2) then
         write(uid, '(99(1X,I0,":",A,"[",A,"]"))',advance='no') &
         & 14, "Jm_u_x", trim(t_unit_current%name), &
         & 15, "Jm_u_y", trim(t_unit_current%name), &
         & 16, "Jm_u_z", trim(t_unit_current%name), &
         & 17, "Jm_d_x", trim(t_unit_current%name), &
         & 18, "Jm_d_y", trim(t_unit_current%name), &
         & 19, "Jm_d_z", trim(t_unit_current%name)
         if(yn_md=='y') then
            write(uid, '(99(1X,I0,":",A,"[",A,"]"))',advance='no') &
                 & 20, "Jmi_x", trim(t_unit_current%name), &
                 & 21, "Jmi_y", trim(t_unit_current%name), &
                 & 22, "Jmi_z", trim(t_unit_current%name)
         endif
       end if
       write(uid,*)
       flush(uid)

    else  !it>=0
       uid = ofl%fh_rt
       write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
          & it * dt * t_unit_time%conv,    &
          & system%vec_Ac_ext(1:3) * t_unit_ac%conv, &
          & system%vec_E_ext(1:3) * t_unit_elec%conv, &
          & system%vec_Ac(1:3) * t_unit_ac%conv, &
          & system%vec_E(1:3) * t_unit_elec%conv
       if(system%nspin==1) then
          write(uid, "(99(1X,E23.15E3))",advance='no') &
          & curr_e(1:3,1) * t_unit_current%conv
       else if(system%nspin==2) then
          write(uid, "(99(1X,E23.15E3))",advance='no') &
          & curr_e(1:3,1) * t_unit_current%conv, &
          & curr_e(1:3,2) * t_unit_current%conv
       end if
       if(yn_md=='y') then
          write(uid, "(99(1X,E23.15E3))",advance='no') &
          & curr_i(1:3) * t_unit_current%conv
       endif
       write(uid,*)
    endif
    endif

  end subroutine
  
!===================================================================================================================================

  subroutine write_rt_energy_data(it,ofl,dt,energy,md)
    use structures, only: s_ofile,s_dft_energy,s_md
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use salmon_global, only: ensemble, thermostat, out_rt_energy_step, yn_periodic, yn_jm
    use filesystem, only: open_filehandle
    use inputoutput, only: yn_md,t_unit_time,t_unit_energy
    implicit none
    type(s_dft_energy) :: energy
    type(s_md) :: md
    type(s_ofile) :: ofl
    integer, intent(in) :: it
    integer :: uid
    real(8) :: dt

    if(comm_is_root(nproc_id_global)) then

    if(it.lt.0) then  !print header
       ofl%fh_rt_energy = open_filehandle(ofl%file_rt_energy_data)
       uid = ofl%fh_rt_energy

10     format("#",1X,A,":",1X,A)
       write(uid,10) "Real time calculation",""
       write(uid,10) "Eall", "Total energy"
       write(uid,10) "Eall0", "Initial energy"
       if(yn_md=='y') then
       write(uid,10) "Tion", "Kinetic energy of ions"
       write(uid,10) "Temperature_ion", "Temperature of ions"
       write(uid,10) "E_work", "Work energy of ions(sum f*dr)"
       if(ensemble=="NVT".and.thermostat=="nose-hoover")then
       write(uid,10) "Enh", "NH thermostat energy (MD)"
       write(uid,10) "Hnvt", "Hamiltonian with NH thermostat(MD)"
       write(uid,10) "Hnvt'","Hnvt using E_work"
       endif
       endif

       if(yn_periodic=='y' .and. yn_jm=='y') then
         write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
          & 1, "Time", trim(t_unit_time%name), &
          & 2, "Eall-Eall0", trim(t_unit_energy%name)
       else
         write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
          & 1, "Time", trim(t_unit_time%name), &
          & 2, "Eall", trim(t_unit_energy%name), &
          & 3, "Eall-Eall0", trim(t_unit_energy%name)
       end if

       if(yn_md=='y') then
       write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
        & 4, "Tion", trim(t_unit_energy%name), &
        & 5, "Temperature_ion", "K", &
        & 6, "E_work", trim(t_unit_energy%name)
       if(ensemble=="NVT".and.thermostat=="nose-hoover")then
       write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
        & 7, "Enh",  trim(t_unit_energy%name), &
        & 8, "Hnvt", trim(t_unit_energy%name), &
        & 9, "Hnvt'",trim(t_unit_energy%name)
       endif
       endif

       write(uid,*)
       flush(uid)
       
       if(yn_periodic=='y' .and. yn_jm=='y') then
         write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
             & 0d0,        &
             & 0d0
       else
         write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
             & 0d0,        &
             & energy%E_tot0 * t_unit_energy%conv, &
             & 0d0
       end if
       if(yn_md=='y') then
         write(uid, "(99(1X,E23.15E3))",advance='no') &
             & md%Tene * t_unit_energy%conv, &
             & md%Temperature,               &
             & md%E_work * t_unit_energy%conv
       endif
       
       write(uid,*)
       flush(uid)

    else  !it>=0
       if(mod(it,out_rt_energy_step)==0)then
          uid = ofl%fh_rt_energy
   
          if(yn_periodic=='y' .and. yn_jm=='y') then
            write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
               & it * dt * t_unit_time%conv,        &
               & (energy%E_tot-energy%E_tot0) * t_unit_energy%conv
          else
            write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
               & it * dt * t_unit_time%conv,        &
               & energy%E_tot * t_unit_energy%conv, &
               & (energy%E_tot-energy%E_tot0) * t_unit_energy%conv
          end if
          if(yn_md=='y') then
          write(uid, "(99(1X,E23.15E3))",advance='no') &
             & md%Tene * t_unit_energy%conv, &
             & md%Temperature,               &
             & md%E_work * t_unit_energy%conv
          endif
   
   
   !        write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
   !          & it * dt * t_unit_time%conv, &
   !          & Eall_t(it) * t_unit_energy%conv, &
   !          & (Eall_t(it) - Eall0) * t_unit_energy%conv
   !        if(yn_md=='y') then
   !        write(uid, "(99(1X,E23.15E3))",advance='no') &
   !          & Tion_t(it) * t_unit_energy%conv, &
   !          & Temperature_ion_t(it), &
   !          & Ework_integ_fdR(it) * t_unit_energy%conv
   !        if(ensemble=="NVT".and.thermostat=="nose-hoover")then
   !        write(uid, "(99(1X,E23.15E3))",advance='no') &
   !          & Enh_t(it) * t_unit_energy%conv, &
   !          & Hnvt_t(it) * t_unit_energy%conv, &
   !          & (Tion_t(it)+Ework_integ_fdR(it)+Enh_t(it)) * t_unit_energy%conv
   !        endif
   !        endif
          write(uid,*)
          flush(uid)  !for debug
       end if
    endif
    endif

  end subroutine
  
!===================================================================================================================================
  subroutine write_response_0d(ofl,rt)
    use salmon_global, only: e_impulse, nt, dt, nenergy, de
    use inputoutput, only: t_unit_energy,t_unit_polarizability
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use structures, only: s_ofile, s_rt
    use filesystem, only: open_filehandle
    implicit none
    type(s_ofile) :: ofl
    type(s_rt),intent(in) :: rt
    integer :: uid
    integer :: ihw,n,ixyz
    real(8) :: tt,hw,t2
    complex(8) :: zalpha(3)
    real(8) :: sf(3)

    if (comm_is_root(nproc_id_global)) then
      ofl%fh_response        = open_filehandle(ofl%file_response_data)
      uid = ofl%fh_response

10    format("#",1X,A,":",1X,A)
      write(uid,10) "Fourier-transform spectra",""
      write(uid,10) "alpha", "Polarizability"
      write(uid,10) "df/dE", "Strength function"

      write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "Energy", trim(t_unit_energy%name), &
        & 2, "Re(alpha_x)", trim(t_unit_polarizability%name), &
        & 3, "Re(alpha_y)", trim(t_unit_polarizability%name), &
        & 4, "Re(alpha_z)", trim(t_unit_polarizability%name), &
        & 5, "Im(alpha_x)", trim(t_unit_polarizability%name), &
        & 6, "Im(alpha_y)", trim(t_unit_polarizability%name), &
        & 7, "Im(alpha_z)", trim(t_unit_polarizability%name), &
        & 8, "df_x/dE", "none", &
        & 9, "df_y/dE", "none", &
        & 10,"df_z/dE", "none"

      tt = dt*dble(nt)

      do ihw=1,nenergy
        hw=dble(ihw)*de ; zalpha(:)=(0.d0,0.d0)  
        do n=1,nt
          t2=dble(n)*dt ; zalpha(:)=zalpha(:)+exp(zi*hw*t2)*rt%dDp_e(:,n) & 
                                             *(1-3*(t2/tt)**2+2*(t2/tt)**3)
        end do

        zalpha(:)=zalpha(:)/e_impulse*dt
        sf(:)=-2*hw/pi*aimag(zalpha(:))

        write(uid,'(F16.8,99(1X,E23.15E3))') hw * t_unit_energy%conv &
             &,(real(zalpha(ixyz))*t_unit_polarizability%conv,ixyz=1,3)&
             &,(aimag(zalpha(ixyz))*t_unit_polarizability%conv,ixyz=1,3)&
             &,(sf(ixyz),ixyz=1,3)

      end do

      flush(uid)  !for debug

    end if

  end subroutine

!===================================================================================================================================
  subroutine write_response_3d(ofl,rt)
    use salmon_global, only: e_impulse, trans_longi, nt, dt, nenergy, de, temperature, yn_lr_w0_correction
    use inputoutput, only: t_unit_energy,t_unit_conductivity
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use structures, only: s_ofile, s_rt
    use filesystem, only: open_filehandle
    implicit none
    type(s_ofile) :: ofl
    type(s_rt),intent(in) :: rt
    integer :: uid
    integer :: ihw,n,ixyz, it
    real(8) :: tt,hw,t2, smoothing, jav_d(3), jav_s(3), smt_s
    complex(8) :: zsigma(3),zeps(3)
    real(8) :: rtdata(1:3,1:nt)

    if (comm_is_root(nproc_id_global)) then
      ofl%fh_response  = open_filehandle(ofl%file_response_data)
      uid = ofl%fh_response

10    format("#",1X,A,":",1X,A)
      write(uid,10) "Fourier-transform spectra",""
      write(uid,10) "sigma", "Conductivity"
      write(uid,10) "eps", "Dielectric constant"

      write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1,  "Energy", trim(t_unit_energy%name), &
        & 2,  "Re(sigma_x)", trim(t_unit_conductivity%name), &
        & 3,  "Re(sigma_y)", trim(t_unit_conductivity%name), &
        & 4,  "Re(sigma_z)", trim(t_unit_conductivity%name), &
        & 5,  "Im(sigma_x)", trim(t_unit_conductivity%name), &
        & 6,  "Im(sigma_y)", trim(t_unit_conductivity%name), &
        & 7,  "Im(sigma_z)", trim(t_unit_conductivity%name), &
        & 8,  "Re(eps_x)", "none", &
        & 9,  "Re(eps_y)", "none", &
        & 10, "Re(eps_z)", "none", &
        & 11, "Im(eps_x)", "none", &
        & 12, "Im(eps_y)", "none", &
        & 13, "Im(eps_z)", "none"

     
      tt = dt*dble(nt)
      if(trans_longi=="tr") then
        rtdata(1:3,1:nt) = rt%curr(1:3,1:nt)
      else if(trans_longi=="lo")then
        rtdata(1:3,1:nt) = rt%E_tot(1:3,1:nt)
      end if

      ! sigma(omega=0) correcton
      jav_d(:)=0d0
      if( yn_lr_w0_correction=='y' .and. temperature < 0d0) then
         jav_s = 0d0
         smt_s = 0d0;
         do it = 1,nt
            t2 = dble(it)*dt
            smoothing = 1d0 - 3d0*(t2/tt)**2 + 2d0*(t2/tt)**3
            jav_s(:)  = jav_s(:) + rtdata(:,it) * smoothing
            smt_s     = smt_s + smoothing
         end do
         jav_d(:) = jav_s(:)/smt_s
      end if

      do ihw=1,nenergy
        hw=dble(ihw)*de
        zsigma(:)=(0.d0,0.d0)
        do n=1,nt
          t2=dble(n)*dt
          smoothing = 1d0 - 3d0*(t2/tt)**2 + 2d0*(t2/tt)**3
         !zsigma(:)=zsigma(:) + exp(zi*hw*t2)* rtdata(:,n) * smoothing
          zsigma(:)=zsigma(:) + exp(zi*hw*t2)* (rtdata(:,n)-jav_d(:)) * smoothing
        end do

        zsigma(:) = (zsigma(:)/e_impulse)*dt
        if(trans_longi=="tr")then
          zeps(:)=1.d0+4.d0*pi*zi*zsigma(:)/hw
        else if(trans_longi=="lo")then
          zeps(:)=1.d0/(1.d0-zsigma(:))
          zsigma(:)=(zeps(:)-1.d0)/(4.d0*pi*zi/hw)
        end if

        write(uid,'(F16.8,99(1X,E23.15E3))') hw * t_unit_energy%conv &
             &,(real(zsigma(ixyz))*t_unit_conductivity%conv,ixyz=1,3)&
             &,(aimag(zsigma(ixyz))*t_unit_conductivity%conv,ixyz=1,3)&
             &,(real(zeps(ixyz)),ixyz=1,3)&
             &,(aimag(zeps(ixyz)),ixyz=1,3)

      end do

      flush(uid)  !for debug

    end if

  end subroutine

!===================================================================================================================================
  subroutine write_dft_md_data(it,ofl,md)
    use structures, only: s_md, s_ofile
    use inputoutput, only: t_unit_time,t_unit_energy
    use salmon_global, only: dt,nt,sysname
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root,comm_sync_all
    use filesystem, only: open_filehandle
    implicit none
    type(s_ofile) :: ofl
    type(s_md) :: md
    integer :: uid, it

    if(it==0 .and. comm_is_root(nproc_id_global)) then
       ofl%file_dft_md = trim(sysname)//'_dft_md.data'
       ofl%fh_dft_md   = open_filehandle(ofl%file_dft_md)
       uid = ofl%fh_dft_md
       open(uid,file=trim(ofl%file_dft_md),status="unknown")
    endif
    uid = ofl%fh_dft_md

    if(comm_is_root(nproc_id_global)) then
      if(it==0) then
         write(uid,'("#",1X,A)') "DFT-MD: adiabatic (ground state) molecular dynamics"
         write(uid,'("#",1X,A,":",1X,A)') "Tene", "Kinetic energy of atoms(ions)"
         write(uid,'("#",1X,A,":",1X,A)') "Uene", "Potential energy"
         write(uid,'("#",1X,A,":",1X,A)') "Uene0", "Initial potential energy"
         write(uid,'("#",1X,A,":",1X,A)') "E_work", "Work energy (sum F*dR)"
         write(uid,'("#",1X,A,":",1X,A)') "E_nh", "Energy of NH thermostat"
         write(uid,'("#",1X,A,":",1X,A)') "E_tot", "Total energy(=Tene+Uene)"
         write(uid,'("#",1X,A,":",1X,A)') "E_tot0", "Initial total energy"
         write(uid,'("#",1X,A,":",1X,A)') "Hnvt", "Hamiltonian with NH thermostat"
         write(uid,'("#",1X,A,":",1X,A)') "Hnvt'", "Hnvt using E_work"
         write(uid,'("#",1X,A,":",1X,A)') "Temperature", "Temperature of atoms"
         write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))') &
              &  1, "Time",   trim(t_unit_time%name), &
              &  2, "Tene",   trim(t_unit_energy%name), &
              &  3, "Uene",   trim(t_unit_energy%name), &
              &  4, "Uene-Uene0", trim(t_unit_energy%name), &
              &  5, "E_work", trim(t_unit_energy%name), &
              &  6, "E_nh",   trim(t_unit_energy%name), &
              &  7, "E_tot",  trim(t_unit_energy%name), &
              &  8, "E_tot-E_tot0", trim(t_unit_energy%name), &
              &  9, "Hnvt(2+3+6)",  trim(t_unit_energy%name), &
              & 10, "Hnvt(2+5+6)'", trim(t_unit_energy%name), &
              & 11, "Temperature_ion", "K"
      endif

      write(uid, "(F16.8,99(1X,E23.15E3))") &
          & it * dt    * t_unit_time%conv, &
          & md%Tene    * t_unit_energy%conv, &
          & md%Uene    * t_unit_energy%conv, &
          & (md%Uene-md%Uene0) * t_unit_energy%conv, &
          & md%E_work  * t_unit_energy%conv, &
          & md%E_nh    * t_unit_energy%conv, &
          & md%E_tot   * t_unit_energy%conv, &
          & (md%E_tot-md%E_tot0) * t_unit_energy%conv, &
          & md%Htot    * t_unit_energy%conv, &
          & (md%Tene+md%E_work+md%E_nh) * t_unit_energy%conv, &
          & md%Temperature
      flush(uid)

      if(it==Nt) close(uid)
    end if

    call comm_sync_all
    return
  end subroutine write_dft_md_data
!===================================================================================================================================

  subroutine write_pulse_0d(ofl,rt)
    use inputoutput, only: nt, dt, nenergy, de,  &
                           t_unit_energy,  &
                           t_unit_spectrum_dipole,  &
                           t_unit_spectrum_dipole_square
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use structures, only: s_ofile, s_rt
    use filesystem, only: open_filehandle
    implicit none
    type(s_ofile) :: ofl
    type(s_rt),intent(in) :: rt
    integer :: uid
    integer :: ihw,n,ixyz
    real(8) :: tt,hw,t2
    complex(8) :: zdDp_e(3)

    if (comm_is_root(nproc_id_global)) then
      ofl%fh_pulse           = open_filehandle(ofl%file_pulse_data)
      uid = ofl%fh_pulse

10    format("#",1X,A,":",1X,A)
      write(uid,10) "Fourier-transform spectra",""
      write(uid,10) "energy", "Frequency"
      write(uid,10) "dm", "Dopile moment"

      write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1,  "energy", trim(t_unit_energy%name), &
        & 2,  "Re(dm_x)", trim(t_unit_spectrum_dipole%name), &
        & 3,  "Re(dm_y)", trim(t_unit_spectrum_dipole%name), &
        & 4,  "Re(dm_z)", trim(t_unit_spectrum_dipole%name), &
        & 5,  "Im(dm_x)", trim(t_unit_spectrum_dipole%name), &
        & 6,  "Im(dm_y)", trim(t_unit_spectrum_dipole%name), &
        & 7,  "Im(dm_z)", trim(t_unit_spectrum_dipole%name), &
        & 8,  "|dm_x|^2", trim(t_unit_spectrum_dipole_square%name), &
        & 9,  "|dm_y|^2", trim(t_unit_spectrum_dipole_square%name), &
        & 10, "|dm_z|^2", trim(t_unit_spectrum_dipole_square%name)

      tt = dt*dble(nt)

      do ihw=1,nenergy
        hw=dble(ihw)*de 
        zdDp_e(:)=(0.d0,0.d0) 
        do n=1,nt
          t2=dble(n)*dt ; zdDp_e(:)=zdDp_e(:)+exp(zi*hw*t2)*rt%dDp_e(:,n) & 
                                             *(1-3*(t2/tt)**2+2*(t2/tt)**3)
        end do
        zdDp_e(:)=zdDp_e(:)*dt

        write(uid,'(F16.8,99(1X,E23.15E3))') hw * t_unit_energy%conv &
             &,(real(zdDp_e(ixyz))*t_unit_spectrum_dipole%conv,ixyz=1,3)&
             &,(aimag(zdDp_e(ixyz))*t_unit_spectrum_dipole%conv,ixyz=1,3)&
             &,(abs(zdDp_e(ixyz))**2*t_unit_spectrum_dipole_square%conv,ixyz=1,3)
      end do

    end if

  end subroutine

!===================================================================================================================================
  subroutine write_pulse_3d(ofl,rt)
    use inputoutput, only: nt, dt, nenergy, de,  &
                           t_unit_energy,  &
                           t_unit_spectrum_current,  &
                           t_unit_spectrum_current_square,  &
                           t_unit_spectrum_elec,  &
                           t_unit_spectrum_elec_square
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use structures, only: s_ofile, s_rt
    use filesystem, only: open_filehandle
    implicit none
    type(s_ofile) :: ofl
    type(s_rt),intent(in) :: rt
    integer :: uid
    integer :: ihw,n,ixyz
    real(8) :: tt,hw,t2
    complex(8) :: zcurr(3),zE_ext(3),zE_tot(3)

    if (comm_is_root(nproc_id_global)) then
      ofl%fh_pulse           = open_filehandle(ofl%file_pulse_data)
      uid = ofl%fh_pulse

10    format("#",1X,A,":",1X,A)
      write(uid,10) "Fourier-transform spectra",""
      write(uid,10) "energy", "Frequency"
      write(uid,10) "Jm", "Matter current"
      write(uid,10) "E_ext", "External electric field"
      write(uid,10) "E_tot", "Total electric field"

      write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1,  "energy", trim(t_unit_energy%name), &
        & 2,  "Re(Jm_x)", trim(t_unit_spectrum_current%name), &
        & 3,  "Re(Jm_y)", trim(t_unit_spectrum_current%name), &
        & 4,  "Re(Jm_z)", trim(t_unit_spectrum_current%name), &
        & 5,  "Im(Jm_x)", trim(t_unit_spectrum_current%name), &
        & 6,  "Im(Jm_y)", trim(t_unit_spectrum_current%name), &
        & 7,  "Im(Jm_z)", trim(t_unit_spectrum_current%name), &
        & 8,  "|Jm_x|^2", trim(t_unit_spectrum_current_square%name), &
        & 9,  "|Jm_y|^2", trim(t_unit_spectrum_current_square%name), &
        & 10, "|Jm_z|^2", trim(t_unit_spectrum_current_square%name), &
        & 11, "Re(E_ext_x)", trim(t_unit_spectrum_elec%name), &
        & 12, "Re(E_ext_y)", trim(t_unit_spectrum_elec%name), &
        & 13, "Re(E_ext_z)", trim(t_unit_spectrum_elec%name), &
        & 14, "Im(E_ext_x)", trim(t_unit_spectrum_elec%name), &
        & 15, "Im(E_ext_y)", trim(t_unit_spectrum_elec%name), &
        & 16, "Im(E_ext_z)", trim(t_unit_spectrum_elec%name), &
        & 17, "|E_ext_x|^2", trim(t_unit_spectrum_elec_square%name), &
        & 18, "|E_ext_y|^2", trim(t_unit_spectrum_elec_square%name), &
        & 19, "|E_ext_z|^2", trim(t_unit_spectrum_elec_square%name), &
        & 20, "Re(E_tot_x)", trim(t_unit_spectrum_elec%name), &
        & 21, "Re(E_tot_y)", trim(t_unit_spectrum_elec%name), &
        & 22, "Re(E_tot_z)", trim(t_unit_spectrum_elec%name), &
        & 23, "Im(E_tot_x)", trim(t_unit_spectrum_elec%name), &
        & 24, "Im(E_tot_y)", trim(t_unit_spectrum_elec%name), &
        & 25, "Im(E_tot_z)", trim(t_unit_spectrum_elec%name), &
        & 26, "|E_tot_x|^2", trim(t_unit_spectrum_elec_square%name), &
        & 27, "|E_tot_y|^2", trim(t_unit_spectrum_elec_square%name), &
        & 28, "|E_tot_z|^2", trim(t_unit_spectrum_elec_square%name)

      tt = dt*dble(nt)

      do ihw=1,nenergy
        hw=dble(ihw)*de 
        zcurr(:)=(0.d0,0.d0) 
        zE_ext(:)=(0.d0,0.d0)  
        zE_tot(:)=(0.d0,0.d0)  
        do n=1,nt
          t2=dble(n)*dt 
          zcurr(:)=zcurr(:)+exp(zi*hw*t2)*rt%curr(:,n) & 
                                      *(1-3*(t2/tt)**2+2*(t2/tt)**3)
          zE_ext(:)=zE_ext(:)+exp(zi*hw*t2)*rt%E_ext(:,n) & 
                                      *(1-3*(t2/tt)**2+2*(t2/tt)**3)
          zE_tot(:)=zE_tot(:)+exp(zi*hw*t2)*rt%E_tot(:,n) & 
                                      *(1-3*(t2/tt)**2+2*(t2/tt)**3)
        end do
        zcurr(:)=zcurr(:)*dt
        zE_ext(:)=zE_ext(:)*dt
        zE_tot(:)=zE_tot(:)*dt

        write(uid,'(F16.8,99(1X,E23.15E3))') hw * t_unit_energy%conv &
             &,(real(zcurr(ixyz))*t_unit_spectrum_current%conv,ixyz=1,3)&
             &,(aimag(zcurr(ixyz))*t_unit_spectrum_current%conv,ixyz=1,3)&
             &,(abs(zcurr(ixyz))**2*t_unit_spectrum_current_square%conv,ixyz=1,3)&
             &,(real(zE_ext(ixyz))*t_unit_spectrum_elec%conv,ixyz=1,3)&
             &,(aimag(zE_ext(ixyz))*t_unit_spectrum_elec%conv,ixyz=1,3)&
             &,(abs(zE_ext(ixyz))**2*t_unit_spectrum_elec_square%conv,ixyz=1,3)&
             &,(real(zE_tot(ixyz))*t_unit_spectrum_elec%conv,ixyz=1,3)&
             &,(aimag(zE_tot(ixyz))*t_unit_spectrum_elec%conv,ixyz=1,3)&
             &,(abs(zE_tot(ixyz))**2*t_unit_spectrum_elec_square%conv,ixyz=1,3)

      end do

      flush(uid)  !for debug

    end if

  end subroutine

!===================================================================================================================================

  subroutine write_prod_dk_data(rgrid_lg, rgrid_mg, system, wf_info, wavefunction)
    use structures,           only: s_rgrid, s_dft_system, s_parallel_info, s_orbital
    use parallelization,      only: nproc_id_global
    use communication, only: comm_is_root
    use filesystem,          only: open_filehandle
    use inputoutput,          only: sysname, base_directory, num_kgrid
    use band,                 only: calc_kgrid_prod
    implicit none
    type(s_rgrid),        intent(in) :: rgrid_lg, rgrid_mg
    type(s_dft_system),       intent(in) :: system
    type(s_parallel_info),      intent(in) :: wf_info
    type(s_orbital), intent(in) :: wavefunction

    ! Specify the neighboring k-grid region to consider:
    integer, parameter :: ndk = 1 
    ! (ndk=1 corresponds to first nearlest neighbors)

    integer :: ik, ik1, ik2, ik3 
    integer :: jdk1, jdk2, jdk3, io, jo
    integer :: fh
    character(256) :: file_prod_dk_data
    integer :: ik3d_tbl(1:3, 1:system%nk)
    complex(8) :: prod_dk( &
      & 1:system%no, 1:system%no, -ndk:ndk, -ndk:ndk, -ndk:ndk, &
      & 1:system%nk)

    ! Export filename: project_directory/sysname_kprod_dk.data
    file_prod_dk_data = trim(base_directory) // trim(sysname) // "_prod_dk.data"

    ! If k-point is distributed as uniform rectangular grid:
    if (0 < minval(num_kgrid)) then
      ! Calculate inner-product table: prod_dk
      call calc_kgrid_prod( &
        & system, rgrid_lg, rgrid_mg, wf_info, wavefunction, &
        & num_kgrid(1), num_kgrid(2), num_kgrid(3), ndk, &
        & ik3d_tbl, prod_dk)
      
      if(comm_is_root(nproc_id_global)) then
        fh = open_filehandle(trim(file_prod_dk_data))
        write(fh, '(a)') "# 1:ik 2:ik1 3:ik2 4:ik3 5:jk1-ik1 6:jk2-ik2 7:jk3-ik3 8:io 9:jo 10:re 11:im"
        do ik = 1, system%nk
          ik1 = ik3d_tbl(1, ik)
          ik2 = ik3d_tbl(2, ik)
          ik3 = ik3d_tbl(3, ik)
          do jdk3 = -ndk, ndk
            do jdk2 = -ndk, ndk
              do jdk1 = -ndk, ndk
                do jo = 1, system%no
                  do io = 1, system%no
                    write(fh, '(9(i10),2(e25.16e3))') &
                      & ik, ik1, ik2, ik3, &
                      & jdk1, jdk2, jdk3, io, jo, &
                      & real(prod_dk(io, jo, jdk1, jdk2, jdk3, ik)), &
                      & aimag(prod_dk(io, jo, jdk1, jdk2, jdk3, ik))
                  end do
                end do
              end do
            end do
          end do
        end do
        close(fh)
      end if
    end if
    return
  end subroutine write_prod_dk_data
  
!===================================================================================================================================

  !! export SYSNAME_info.data file (GS info)
  subroutine write_info_data(Miter,system,energy,pp)
    use structures
    use salmon_global,       only: natom,nelem,iZatom,nelec,sysname,nstate,nelec_spin,unit_system, &
                                   yn_jm, yn_periodic
    use parallelization,     only: nproc_id_global
    use communication,only: comm_is_root
    use filesystem,         only: open_filehandle
    use inputoutput,         only: au_length_aa,au_energy_ev
    implicit none
    integer           ,intent(in) :: Miter
    type(s_dft_energy),intent(in) :: energy
    type(s_dft_system),intent(in) :: system
    type(s_pp_info)   ,intent(in) :: pp
    !
    integer :: fh,is,p1,p2,p5,iob,jj,ik,ikoa,iatom,ix
    character(100) :: file_gs_info

    file_gs_info = trim(sysname)//"_info.data"
    fh = open_filehandle(trim(file_gs_info))

    if(comm_is_root(nproc_id_global)) then

       write(fh,*) "Total number of iteration = ", Miter
       write(fh,*)
       write(fh,*) "Number of states = ", nstate
       if(sum(nelec_spin(:))==0) then
          write(fh,*) "Number of electrons = ", nelec
       else
          write(fh,*) "Number of electrons = ", (nelec_spin(is),is=1,2)
       end if
       write(fh,*)
       if(yn_jm=='y'.and.yn_periodic=='y') then
         write(fh,*) "For yn_jm = y and yn_periodic=y, this version still cannot output Total Energy."
       else
         write(fh,*) "Total energy (eV) = ", energy%E_tot*au_energy_ev
       end if
       write(fh,*) "1-particle energies (eV)"
       select case (system%nspin)
       case(1)
          do p5=1,(nstate+3)/4
             p1=4*(p5-1)+1
             p2=4*p5 ; if ( p2 > nstate ) p2=nstate
             write(fh,100) (iob,energy%esp(iob,1,1)*au_energy_ev,iob=p1,p2)
          end do
       case(2)
          write(fh,*) "for up-spin"
          do p5=1,(nstate+3)/4
            p1=4*(p5-1)+1
            p2=4*p5 ; if ( p2 > nstate ) p2=nstate
            write(fh,100) (iob,energy%esp(iob,1,1)*au_energy_ev,iob=p1,p2)
          end do
          write(fh,*) "for down-spin"
          do p5=1,(nstate+3)/4
            p1=4*(p5-1)+1
            p2=4*p5 ; if ( p2 > nstate ) p2=nstate
            write(fh,100) (iob,energy%esp(iob,1,2)*au_energy_ev,iob=p1,p2)
          end do
       end select
       write(fh,*)       
100    format(1x,4(i5,f15.4,2x))

       write(fh,200) "Size of the box (A) = ", system%primitive_a*au_length_aa
       write(fh,200) "Grid spacing (A)    = ", (system%Hgs(jj)*au_length_aa,jj=1,3)
       write(fh,*)
200    format(1x,a,30f14.8)

       if(yn_jm=='n')then
         write(fh,'(1x,"Number of atoms = ",i8)') natom
         do ik=1,nelem
            write(fh,'(1x,"iZatom(",i3,")     = ",i8)') ik, iZatom(ik)
         end do
         write(fh,*)
         write(fh,*) "Ref. and max angular momentum",  &
                     " and pseudo-core radius of PP (A)"
         do ikoa=1,nelem
            write(fh,'(1x,"(",i3,")  "," Ref, Max, Rps =",2i4,f8.3)') &
                  ikoa,pp%Lref(ikoa),pp%Mlps(ikoa),pp%Rps(ikoa)*au_length_aa
         end do
         
         write(fh,*)
         select case(unit_system)
         case('au','a.u.')
            write(fh,*) "Force [au] "
            do iatom=1,natom
               write(fh,300) iatom,(system%Force(ix,iatom),ix=1,3)
            end do
         case('A_eV_fs')
            write(fh,*) "Force [eV/A] "
            do iatom=1,natom
               write(fh,300) iatom,(system%Force(ix,iatom)*au_energy_ev/au_length_aa,ix=1,3)
            end do
         end select
300    format(i6,3e16.8)
       end if

       close(fh)
    endif

  end subroutine write_info_data
  
!===================================================================================================================================

  subroutine write_eigen(ofl,system,energy)
    use structures, only: s_ofile, s_dft_system, s_dft_energy
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use inputoutput, only: uenergy_from_au,iperiodic,unit_energy,sysname
    use filesystem, only: open_filehandle
    implicit none
    type(s_ofile),intent(inout) :: ofl
    type(s_dft_system),intent(in) :: system
    type(s_dft_energy),intent(in) :: energy
    integer :: iob,iik,is, uid

    if(comm_is_root(nproc_id_global))then
       ofl%file_eigen_data=trim(sysname)//"_eigen.data"
       ofl%fh_eigen = open_filehandle(trim(ofl%file_eigen_data))
       uid = ofl%fh_eigen
       open(uid,file=ofl%file_eigen_data)
       write(uid,'("#esp: single-particle energies (eigen energies)")')
       write(uid,'("#occ: occupation numbers, io: orbital index")')
       select case(unit_energy)
       case('au','a.u.')
          write(uid,'("# 1:io, 2:esp[a.u.], 3:occ")')
       case('ev','eV')
          write(uid,'("# 1:io, 2:esp[eV], 3:occ")')
       end select
       do is=1,system%nspin
       do iik=1,system%nk
          if(iperiodic==3)then
             write(uid,'("k=",1x,i5,",  spin=",1x,i5)') iik,is
          end if
          do iob=1,system%no
             write(uid,'(1x,i5,2(e26.16e3))') iob, energy%esp(iob,iik,is)*uenergy_from_au, system%rocc(iob,iik,is)
          end do
       end do
       end do
       close(uid)
    end if

  end subroutine write_eigen
  
!===================================================================================================================================

  subroutine write_dos(system,energy)
    use structures
    use math_constants, only: pi
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use inputoutput, only: uenergy_from_au
    use salmon_global, only: out_dos_start, out_dos_end, out_dos_function, &
                           out_dos_width, out_dos_nenergy, yn_out_dos_set_fe_origin, unit_energy, &
                           nelec,nstate,temperature
    use spin_orbit_global, only: SPIN_ORBIT_ON
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_dft_energy),intent(in) :: energy
    !
    integer :: iob,iik,is
    real(8) :: dos_l(1:out_dos_nenergy,system%nspin)
    real(8) :: fk,ww,dw
    integer :: iw,index_vbm
    real(8) :: ene_min,ene_max,eshift

    ene_min = minval(energy%esp)
    ene_max = maxval(energy%esp)
    if(yn_out_dos_set_fe_origin=='y') then
      if(temperature>=0.d0) then
        eshift = system%mu
      else
        if( SPIN_ORBIT_ON )then
          index_vbm=nelec
        else
          index_vbm=nelec/2
        end if
        eshift = maxval(energy%esp(1:index_vbm,:,:)) ! valence band maximum (HOMO)
      end if
    else
      eshift = 0.d0
    endif
    out_dos_start = max(out_dos_start,ene_min-0.25d0*(ene_max-ene_min)-eshift)
    out_dos_end = min(out_dos_end,ene_max+0.25d0*(ene_max-ene_min)-eshift)
    dw=(out_dos_end-out_dos_start)/dble(out_dos_nenergy-1)

    dos_l = 0.d0

    do is=1,system%nspin
    do iik=1,system%nk
    do iob=1,system%no
      select case (out_dos_function)
      case('lorentzian')
        fk=2.d0*out_dos_width/pi
        do iw=1,out_dos_nenergy
          ww=out_dos_start+dble(iw-1)*dw+eshift-energy%esp(iob,iik,is)
          dos_l(iw,is)=dos_l(iw,is)+system%wtk(iik)*fk/(ww**2+out_dos_width**2)
        end do
      case('gaussian')
        fk=2.d0/(sqrt(2.d0*pi)*out_dos_width)
        do iw=1,out_dos_nenergy
          ww=out_dos_start+dble(iw-1)*dw+eshift-energy%esp(iob,iik,is)
          dos_l(iw,is)=dos_l(iw,is)+system%wtk(iik)*fk*exp(-(0.5d0/out_dos_width**2)*ww**2)
        end do
      end select
    end do
    end do
    end do

    if(comm_is_root(nproc_id_global))then
      open(101,file="dos.data")
      write(101,'("# Density of States")')
      select case(unit_energy)
      case('au','a.u.')
        write(101,'("# Energy[a.u.] DOS[a.u.]")')
      case('ev','eV')
        write(101,'("# Energy[eV]  DOS[1/eV]")')
      end select
      write(101,'("#-----------------------")')
      do iw=1,out_dos_nenergy
        ww=out_dos_start+dble(iw-1)*dw
        write(101,'(F16.8,99(1X,E23.15E3))') ww*uenergy_from_au, ( dos_l(iw,is)/uenergy_from_au, is=1,system%nspin )
      end do
      close(101)
    end if

  end subroutine write_dos
  
!===================================================================================================================================
  
  subroutine write_pdos(lg,mg,system,info,pp,energy,tpsi)
    use structures
    use math_constants      ,only: pi,zi
    use salmon_math         ,only: ylm
    use parallelization     ,only: nproc_id_global
    use communication       ,only: comm_is_root, comm_summation
    use salmon_global       ,only: out_dos_start, out_dos_end, out_dos_function, &
                                   out_dos_width, out_dos_nenergy, yn_out_dos_set_fe_origin, &
                                   nelec, kion, natom, nstate, unit_energy, temperature
    use inputoutput         ,only: uenergy_from_au
    use prep_pp_sub         ,only: bisection
    use spin_orbit_global, only: SPIN_ORBIT_ON
    implicit none
    type(s_rgrid)           ,intent(in) :: lg,mg
    type(s_dft_system)      ,intent(in) :: system
    type(s_parallel_info)   ,intent(in) :: info
    type(s_pp_info)         ,intent(in) :: pp
    type(s_dft_energy)      ,intent(in) :: energy
    type(s_orbital)         ,intent(in) :: tpsi
    !
    integer :: iob,iatom,L,m,ix,iy,iz,iik,ispin
    integer :: ikoa
    integer :: intr
    real(8) :: phi_r
    real(8) :: rr
    real(8) :: ratio1,ratio2
    real(8) :: xx,yy,zz
    real(8) :: xxxx,yyyy,zzzz,rinv
    integer :: lm
    real(8) :: rbox_pdos(25,natom)
    real(8) :: rbox_pdos2(25,natom)
    real(8) :: pdos_l_tmp(out_dos_nenergy,0:4,natom)
    real(8) :: pdos_l(out_dos_nenergy,0:4,natom)
    character(100) :: Outfile
    real(8) :: fk,ww,dw
    integer :: iw,index_vbm
    real(8) :: ene_min,ene_max,eshift
    character(20) :: fileNumber

    if( all(pp%upp_f==0.0d0) )then
      write(*,*) "@calc_pdos: Pseudoatom wave function is not available"
      return
    end if

    ene_min = minval(energy%esp(:,:,:))
    ene_max = maxval(energy%esp(:,:,:))
    if(yn_out_dos_set_fe_origin=='y') then
      if(temperature>=0.d0) then
        eshift = system%mu
      else
        if( SPIN_ORBIT_ON )then
          index_vbm=nelec
        else
          index_vbm=nelec/2
        end if
        eshift = maxval(energy%esp(1:index_vbm,:,:)) ! valence band maximum (HOMO)
      end if
    else
      eshift = 0.d0
    endif
    out_dos_start = max(out_dos_start,ene_min-0.25d0*(ene_max-ene_min)-eshift)
    out_dos_end = min(out_dos_end,ene_max+0.25d0*(ene_max-ene_min)-eshift)
    dw=(out_dos_end-out_dos_start)/dble(out_dos_nenergy-1)

    pdos_l_tmp=0.d0

    do ispin=1,system%nspin
    do iik=info%ik_s,info%ik_e
    do iob=info%io_s,info%io_e
      rbox_pdos=0.d0
      do iatom=1,natom
        ikoa=Kion(iatom)
        do L=0,pp%mlps(ikoa)
          do m=-L,L
            lm=L*L+L+1+m
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              xx=lg%coordinate(ix,1)-system%Rion(1,iatom)
              yy=lg%coordinate(iy,2)-system%Rion(2,iatom)
              zz=lg%coordinate(iz,3)-system%Rion(3,iatom)
              rr=sqrt(xx**2+yy**2+zz**2)+1.d-50
              rinv=1.0d0/rr
              xxxx=xx*rinv
              yyyy=yy*rinv
              zzzz=zz*rinv
              call bisection(rr,intr,ikoa,pp%nrmax,pp%rad)
              if(intr==1) intr=2
              ratio1=(rr-pp%rad(intr,ikoa))/(pp%rad(intr+1,ikoa)-pp%rad(intr,ikoa)) ; ratio2=1.d0-ratio1
              phi_r= ratio1*pp%upp_f(intr,pp%lref(ikoa),ikoa)/rr**(pp%lref(ikoa)+1)*sqrt((2*pp%lref(ikoa)+1)/(4*Pi)) +  &
                     ratio2*pp%upp_f(intr-1,pp%lref(ikoa),ikoa)/rr**(pp%lref(ikoa)+1)*sqrt((2*pp%lref(ikoa)+1)/(4*Pi))
                                            !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
              if(allocated(tpsi%rwf)) then
                rbox_pdos(lm,iatom)=rbox_pdos(lm,iatom)+tpsi%rwf(ix,iy,iz,ispin,iob,iik,1)*phi_r*Ylm(xxxx,yyyy,zzzz,L,m)*system%Hvol
              else
                rbox_pdos(lm,iatom)=rbox_pdos(lm,iatom)+tpsi%zwf(ix,iy,iz,ispin,iob,iik,1)*phi_r*Ylm(xxxx,yyyy,zzzz,L,m)*system%Hvol
              end if
            end do
            end do
            end do
          end do
        end do
      end do
      call comm_summation(rbox_pdos,rbox_pdos2,25*natom,info%icomm_r)
      do iatom=1,natom
        ikoa=Kion(iatom)
        do L=0,pp%mlps(ikoa)
          do lm=L**2+1,(L+1)**2
            select case (out_dos_function)
            case('lorentzian')
              fk=2.d0*out_dos_width/pi
              do iw=1,out_dos_nenergy
                ww=out_dos_start+dble(iw-1)*dw+eshift-energy%esp(iob,iik,ispin)
                pdos_l_tmp(iw,L,iatom)=pdos_l_tmp(iw,L,iatom)  &
                  +abs(rbox_pdos2(lm,iatom))**2*fk/(ww**2+out_dos_width**2)
              end do
            case('gaussian')
              fk=2.d0/(sqrt(2.d0*pi)*out_dos_width)
              do iw=1,out_dos_nenergy
                ww=out_dos_start+dble(iw-1)*dw+eshift-energy%esp(iob,iik,ispin)
                pdos_l_tmp(iw,L,iatom)=pdos_l_tmp(iw,L,iatom)  &
                  +abs(rbox_pdos2(lm,iatom))**2*fk*exp(-(0.5d0/out_dos_width**2)*ww**2)
              end do
            end select
          end do
        end do
      end do
    end do
    end do
    end do
    call comm_summation(pdos_l_tmp,pdos_l,out_dos_nenergy*5*natom,info%icomm_ko)

    if(comm_is_root(nproc_id_global))then
      do iatom=1,natom
        ikoa=Kion(iatom)
        write(fileNumber, '(i8)') iatom
        OutFile = "pdos"//trim(adjustl(fileNumber))//".data"
        open(101,file=OutFile)
        write(101,'("# Projected Density of States")')
        select case(unit_energy)
        case('au','a.u.')
          if(pp%mlps(ikoa)==0)then
            write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.]")')
          else if(pp%mlps(ikoa)==1)then
            write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.]")')
          else if(pp%mlps(ikoa)==2)then
            write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.] PDOS(l=2)[a.u.]")')
          else if(pp%mlps(ikoa)==3)then
            write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.] PDOS(l=2)[a.u.] PDOS(l=3)[a.u.]")')
          end if
        case('ev','eV')
          if(pp%mlps(ikoa)==0)then
            write(101,'("# Energy[eV]  PDOS(l=0)[1/eV]")')
          else if(pp%mlps(ikoa)==1)then
            write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV]")')
          else if(pp%mlps(ikoa)==2)then
            write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV] PDOS(l=2)[1/eV]")')
          else if(pp%mlps(ikoa)==3)then
            write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV] PDOS(l=2)[1/eV] PDOS(l=3)[1/eV]")')
          end if
        end select
        write(101,'("#-----------------------")')
        if(pp%mlps(ikoa)==0)then
          do iw=1,out_dos_nenergy
            ww=out_dos_start+dble(iw-1)*dw
            write(101,'(f10.5,f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,pp%mlps(ikoa))
          end do
        else if(pp%mlps(ikoa)==1)then
          do iw=1,out_dos_nenergy
            ww=out_dos_start+dble(iw-1)*dw
            write(101,'(f10.5,2f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,pp%mlps(ikoa))
          end do
        else if(pp%mlps(ikoa)==2)then
          do iw=1,out_dos_nenergy
            ww=out_dos_start+dble(iw-1)*dw
            write(101,'(f10.5,3f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,pp%mlps(ikoa))
          end do
        else if(pp%mlps(ikoa)==3)then
          do iw=1,out_dos_nenergy
            ww=out_dos_start+dble(iw-1)*dw
            write(101,'(f10.5,4f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,pp%mlps(ikoa))
          end do
        end if
        close(101)
      end do
    end if
    return
  end subroutine write_pdos
  
!===================================================================================================================================

  subroutine write_band_information(system,energy)
    use structures
    use salmon_global, only: nelec
    use inputoutput, only: au_energy_ev
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use spin_orbit_global, only: SPIN_ORBIT_ON
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_dft_energy),intent(in) :: energy
    !
    integer :: ik,index_vbm
    real(8),dimension(system%nk) :: esp_vb_min,esp_vb_max,esp_cb_min,esp_cb_max
    !
    if( SPIN_ORBIT_ON )then
      index_vbm=nelec
    else
      index_vbm=nelec/2
    end if
    if(comm_is_root(nproc_id_global) .and. index_vbm<system%no) then
      do ik=1,system%nk
        esp_vb_min(ik)=minval(energy%esp(1:index_vbm,ik,:))
        esp_vb_max(ik)=maxval(energy%esp(1:index_vbm,ik,:))
        esp_cb_min(ik)=minval(energy%esp(index_vbm+1:system%no,ik,:))
        esp_cb_max(ik)=maxval(energy%esp(index_vbm+1:system%no,ik,:))
      end do
      write(*,*) 'band information-----------------------------------------'
      write(*,*) 'Bottom of VB',minval(esp_vb_min(:))
      write(*,*) 'Top of VB',maxval(esp_vb_max(:))
      write(*,*) 'Bottom of CB',minval(esp_cb_min(:))
      write(*,*) 'Top of CB',maxval(esp_cb_max(:))
      write(*,*) 'Fundamental gap',minval(esp_cb_min(:))-maxval(esp_vb_max(:))
      write(*,*) 'BG between same k-point',minval(esp_cb_min(:)-esp_vb_max(:))
      write(*,*) 'Physicaly upper bound of CB for DOS',minval(esp_cb_max(:))
      write(*,*) 'Physicaly upper bound of CB for eps(omega)',minval(esp_cb_max(:)-esp_vb_min(:))
      write(*,*) '---------------------------------------------------------'
      write(*,*) 'Bottom of VB[eV]',minval(esp_vb_min(:))*au_energy_ev
      write(*,*) 'Top of VB[eV]',maxval(esp_vb_max(:))*au_energy_ev
      write(*,*) 'Bottom of CB[eV]',minval(esp_cb_min(:))*au_energy_ev
      write(*,*) 'Top of CB[eV]',maxval(esp_cb_max(:))*au_energy_ev
      write(*,*) 'Fundamental gap[eV]',(minval(esp_cb_min(:))-maxval(esp_vb_max(:)))*au_energy_ev
      write(*,*) 'BG between same k-point[eV]',(minval(esp_cb_min(:)-esp_vb_max(:)))*au_energy_ev
      write(*,*) '---------------------------------------------------------'
    end if
    return
  end subroutine write_band_information
  
!===================================================================================================================================
  subroutine init_projection(system,lg,mg,info,ttpsi,V_local,rt)
    use structures
    use communication, only: comm_is_root, comm_bcast, comm_summation, comm_sync_all
    use parallelization, only: nproc_id_global,nproc_group_global
    use salmon_global, only: projection_option,directory_read_data
    use checkpoint_restart_sub, only: read_wavefunction
    implicit none
    type(s_rgrid)           ,intent(in) :: lg,mg
    type(s_dft_system)                  :: system
    type(s_parallel_info)               :: info
    type(s_orbital)                     :: ttpsi, spsi_tmp !ttpsi= spsi_in
    type(s_scalar)          ,intent(in) :: V_local(system%nspin)
    type(s_rt)                          :: rt
    integer :: jspin

    character(256) :: dir_gs
    character(1024) :: dir_file_in
    integer :: m, iu1, iu2, mk, mo,io, itt_tmp, nprocs_tmp, comm, im, nproc_ob
    logical :: if_real_orbital, iself


    if(projection_option=='gs') then

       call allocate_orbital_complex(system%nspin,mg,info,rt%tpsi0)

       !$omp workshare
       rt%tpsi0%zwf = ttpsi%zwf
       !$omp end workshare

    else if(projection_option=='gx') then

       ! read gs data 
       dir_gs = directory_read_data   ! default= ./restart/
       comm = nproc_group_global
       iu1 = 96
       if(comm_is_root(nproc_id_global))then
          dir_file_in = trim(dir_gs)//"info.bin"
          open(iu1,file=dir_file_in,form='unformatted')
          read(iu1) mk
          read(iu1) mo
          read(iu1) itt_tmp
          read(iu1) nprocs_tmp
          read(iu1) if_real_orbital
          close(iu1)
          rt%no0 = mo
       end if
       call comm_bcast(mk,comm)
       call comm_bcast(mo,comm)
       call comm_bcast(rt%no0,comm)
       call comm_bcast(if_real_orbital,comm)

       if(comm_is_root(nproc_id_global)) then
          write(*,*) "  initialize projection option(gs)"
          write(*,*) "    num. of rt orbitals=", system%no
          write(*,*) "    num. of gs orbitals=", rt%no0
       endif
    
       !occupation
       allocate(rt%rocc0(rt%no0, system%nk, system%nspin))
       if(comm_is_root(nproc_id_global))then
          dir_file_in = trim(dir_gs)//"occupation.bin"
          open(iu1,file=dir_file_in,form='unformatted')
          read(iu1) rt%rocc0(1:mo,1:mk,1:system%nspin)
          close(iu1)
       end if
       call comm_bcast(rt%rocc0,comm)

       ! make orbital info for gs
       nproc_ob = info%nporbital
       rt%io_s0 = (info%id_o * rt%no0) / nproc_ob + 1
       rt%io_e0 = ((info%id_o+1) * rt%no0) / nproc_ob
       rt%numo0 = rt%io_e0 - rt%io_s0 + 1

       allocate(rt%io_s_all0(0:nproc_ob-1))
       allocate(rt%io_e_all0(0:nproc_ob-1))
       allocate(rt%numo_all0(0:nproc_ob-1))
       do m=0,nproc_ob-1
          rt%io_s_all0(m) = m * rt%no0 / nproc_ob + 1
          rt%io_e_all0(m) = (m+1) * rt%no0 / nproc_ob
          rt%numo_all0(m) = rt%io_e_all0(m) - rt%io_s_all0(m) + 1
       end do
       rt%numo_max0 = maxval(rt%numo_all0(0:nproc_ob-1))

       if ( allocated(rt%irank_io0) ) deallocate(rt%irank_io0)
       allocate(rt%irank_io0(1:rt%no0))
       do io=1, rt%no0
          if(mod(io*nproc_ob,rt%no0)==0)then
             rt%irank_io0(io) = io*nproc_ob/rt%no0 - 1
          else
             rt%irank_io0(io) = io*nproc_ob/rt%no0
          end if
       end do

       ! store original orbital info values of rt
       rt%no_org   = system%no
       rt%io_s_org = info%io_s
       rt%io_e_org = info%io_e
       rt%numo_org     = info%numo
       rt%numo_max_org = info%numo_max

       allocate(rt%io_s_all_org(0:nproc_ob-1))
       allocate(rt%io_e_all_org(0:nproc_ob-1))
       allocate(rt%numo_all_org(0:nproc_ob-1))
       allocate(rt%irank_io_org(1:system%no))
       rt%io_s_all_org(:) = info%io_s_all(:)
       rt%io_e_all_org(:) = info%io_e_all(:)
       rt%numo_all_org(:) = info%numo_all(:)
       rt%irank_io_org(:) = info%irank_io(:)

       ! switch no,io_s,io_e,numo,numo_max : rt => gs
       call norb_rt2gs_for_proj(system,info,rt)

       !wavefunction
       iself=.false. ! assuming yn_restart=='n' .and. yn_self_checkpoint=='n'
       call allocate_orbital_complex(system%nspin,mg,info,rt%tpsi0)
       call allocate_orbital_complex(system%nspin,mg,info,spsi_tmp)
       call read_wavefunction(dir_gs,lg,mg,system,info,spsi_tmp,mk,mo,if_real_orbital,iself)

       !$omp workshare
       rt%tpsi0%zwf(:,:,:,:,:,:,:) = spsi_tmp%zwf(:,:,:,:,:,:,:)
       !$omp end workshare

       ! switch back no,io_s,io_e,numo,numo_max : gs => rt
       call norb_gs2rt_for_proj(system,info,rt)

    endif  ! projection_option


    allocate(rt%vloc0(system%nspin))
    do jspin=1,system%nspin
      call allocate_scalar(mg,rt%vloc0(jspin))
      !$omp workshare
      rt%vloc0(jspin)%f = V_local(jspin)%f
      !$omp end workshare
    end do

  end subroutine init_projection

  subroutine norb_rt2gs_for_proj(system,info,rt)
    use structures, only: s_dft_system, s_parallel_info, s_rt
    implicit none
    type(s_dft_system)   ,intent(inout) :: system
    type(s_parallel_info),intent(inout) :: info
    type(s_rt)           ,intent(in) :: rt
    system%no = rt%no0  
    info%io_s = rt%io_s0
    info%io_e = rt%io_e0
    info%numo     = rt%numo0
    info%numo_max = rt%numo_max0
    info%io_s_all(:) = rt%io_s_all0(:)
    info%io_e_all(:) = rt%io_e_all0(:)
    info%numo_all(:) = rt%numo_all0(:)
    deallocate(info%irank_io)
    allocate(info%irank_io(1:system%no))
    info%irank_io(:) = rt%irank_io0(:)
  end subroutine norb_rt2gs_for_proj

  subroutine norb_gs2rt_for_proj(system,info,rt)
    use structures, only: s_dft_system, s_parallel_info, s_rt
    implicit none
    type(s_dft_system)   ,intent(inout) :: system
    type(s_parallel_info),intent(inout) :: info
    type(s_rt)           ,intent(in)  :: rt
    system%no = rt%no_org
    info%io_s = rt%io_s_org
    info%io_e = rt%io_e_org
    info%numo     = rt%numo_org
    info%numo_max = rt%numo_max_org
    info%io_s_all(:) = rt%io_s_all_org(:)
    info%io_e_all(:) = rt%io_e_all_org(:)
    info%numo_all(:) = rt%numo_all_org(:)
    deallocate(info%irank_io)
    allocate(info%irank_io(1:system%no))
    info%irank_io(:) = rt%irank_io_org(:)
  end subroutine norb_gs2rt_for_proj

  subroutine projection(itt,ofl,dt,mg,system,info,stencil,ppg,psi_t,tpsi,ttpsi,srg,energy,rt)
    use structures
    use salmon_global, only: projection_option
    implicit none
    integer                 ,intent(in) :: itt
    type(s_ofile)           ,intent(in) :: ofl
    real(8)                 ,intent(in) :: dt
    type(s_rgrid)           ,intent(in) :: mg
    type(s_dft_system)                  :: system
    type(s_parallel_info)               :: info
    type(s_stencil)         ,intent(in) :: stencil
    type(s_pp_grid)         ,intent(in) :: ppg
    type(s_orbital)         ,intent(in) :: psi_t ! | u_{n,k}(t) >
    type(s_orbital)                     :: tpsi,ttpsi ! temporary arrays
    type(s_sendrecv_grid)               :: srg
    type(s_dft_energy)                  :: energy
    type(s_rt)                          :: rt

    select case(projection_option)
    case('gs')
       call  projection_gs(itt,ofl,dt,mg,system,info,stencil,ppg,psi_t,tpsi,ttpsi,srg,energy,rt)
    case('gx')
       call  projection_gx(itt,ofl,dt,mg,system,info,stencil,ppg,psi_t,tpsi,ttpsi,srg,energy,rt)
    end select

    return
  end subroutine projection

  subroutine projection_gs(itt,ofl,dt,mg,system,info,stencil,ppg,psi_t,tpsi,ttpsi,srg,energy,rt)
    use structures
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global
    use salmon_global, only: ncg,nelec,yn_spinorbit,nscf
    use inputoutput, only: t_unit_time
    use subspace_diagonalization, only: ssdg
    use gram_schmidt_orth, only: gram_schmidt
    use Conjugate_Gradient, only: gscg_zwf
    use Total_Energy, only: calc_eigen_energy
    implicit none
    integer                 ,intent(in) :: itt
    type(s_ofile)           ,intent(in) :: ofl
    real(8)                 ,intent(in) :: dt
    type(s_rgrid)           ,intent(in) :: mg
    type(s_dft_system)      ,intent(in) :: system
    type(s_parallel_info)   ,intent(in) :: info
    type(s_stencil)         ,intent(in) :: stencil
    type(s_pp_grid)         ,intent(in) :: ppg
    type(s_orbital)         ,intent(in) :: psi_t ! | u_{n,k}(t) >
    type(s_orbital)                     :: tpsi,ttpsi ! temporary arrays
    type(s_sendrecv_grid)               :: srg
    type(s_dft_energy)                  :: energy
    type(s_rt)                          :: rt
    !
    integer :: nspin,nspin_tmp,no,nk,ik_s,ik_e,io_s,io_e,is(3),ie(3)
    integer :: ix,iy,iz,io1,io2,io,ik,ispin,iter_GS,niter
    complex(8),dimension(system%no,system%no,system%nspin,system%nk) :: mat
    real(8) :: coef(system%no,system%nk,system%nspin)
    real(8) :: nee, neh, wspin, dE
    complex(8) :: cbox
      
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
    
    if(nspin==1) then
      wspin = 2d0
    else if(nspin==2) then
      wspin = 1d0
    end if
    
    nspin_tmp = nspin
    if(yn_spinorbit=='y') then
      nspin_tmp = 1
    end if
    
    if(nscf==0) then
      niter = 10
    else
      niter = nscf
    end if
    
  ! rt%tpsi0 = | u_{n,k+A(t)/c} >, the ground-state wavefunction whose k-point is shifted by A(t).
    call calc_eigen_energy(energy,rt%tpsi0,tpsi,ttpsi,system,info,mg,rt%vloc0,stencil,srg,ppg)
    dE = energy%E_kin - rt%E_old
    if(abs(dE) < 1e-12) then
      if(comm_is_root(nproc_id_global)) write(*,*) "projection: already converged, E_kin(new)-E_kin(old)=",dE
    else
      do iter_GS=1,niter
        call ssdg(mg,system,info,stencil,rt%tpsi0,tpsi,ppg,rt%vloc0,srg)
        call gscg_zwf(ncg,mg,system,info,stencil,ppg,rt%vloc0,srg,rt%tpsi0,rt%cg)
        call gram_schmidt(system, mg, info, rt%tpsi0)
        call calc_eigen_energy(energy,rt%tpsi0,tpsi,ttpsi,system,info,mg,rt%vloc0,stencil,srg,ppg)
        dE = energy%E_kin - rt%E_old
        if(comm_is_root(nproc_id_global)) write(*,'(a,i6,e20.10)') "projection: ",iter_GS,dE
        if(abs(dE) < 1e-6) exit
        rt%E_old = energy%E_kin
      end do
    end if
    
    call inner_product(rt%tpsi0,psi_t,mat) ! mat(n,m) = < u_{n,k+A(t)/c} | u_{m,k}(t) >

    if(yn_spinorbit=='y') then
      mat(1:no,1:no,1,1:nk) = mat(1:no,1:no,1,1:nk) + mat(1:no,1:no,2,1:nk)
      mat(1:no,1:no,2,1:nk) = mat(1:no,1:no,1,1:nk)
    end if

    coef=0.d0
    do ispin=1,nspin
    do ik=1,nk
    do io1=1,no
      do io2=1,no
        coef(io1,ik,ispin) = coef(io1,ik,ispin) &
        & + system%rocc(io2,ik,ispin)*system%wtk(ik)* abs(mat(io1,io2,ispin,ik))**2
      end do
    end do
    end do
    end do

    nee = 0d0
    neh = dble(nelec)
    do ispin=1,nspin_tmp
    do ik=1,nk
    do io=1,no
    ! /wspin: for canceling double counting of 2 in rocc.
      nee = nee + ((wspin-system%rocc(io,ik,ispin))/wspin) * coef(io,ik,ispin)
      neh = neh - (system%rocc(io,ik,ispin)/wspin) * coef(io,ik,ispin)
    end do
    end do
    end do
   !nee  = sum(ovlp_occ(NBoccmax+1:NB,:))
   !neh  = sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
    if(comm_is_root(nproc_id_global))then
      write(ofl%fh_nex,'(99(1X,E23.15E3))') dble(itt)*dt*t_unit_time%conv, nee, neh
      write(ofl%fh_ovlp,'(i11)') itt
      do ispin=1,nspin_tmp
      do ik=1,nk
        write(ofl%fh_ovlp,'(i6,1000(1X,E23.15E3))') ik,(coef(io,ik,ispin)*nk,io=1,no)
      end do
      end do
    end if
    
!        if(action=="proj_last ") then
!       tconv = t_unit_energy%conv
!       do ik=1,NK
!          write(409,10)ik,(esp_all(ib,ik)*tconv,ovlp_occ(ib,ik)*NKxyz,ib=1,NB)
!       end do

    return
    
  contains
  
    subroutine inner_product(psi1,psi2,mat)
      use communication, only: comm_summation, comm_bcast
      use pack_unpack, only: copy_data
      implicit none
      type(s_orbital), intent(in) :: psi1,psi2
      complex(8),dimension(system%no,system%no,system%nspin,system%nk) :: mat
      !
      complex(8),dimension(system%no,system%no,system%nspin,system%nk) :: mat1
      complex(8) :: wf_io1(mg%is_array(1):mg%ie_array(1) &
                        & ,mg%is_array(2):mg%ie_array(2) &
                        & ,mg%is_array(3):mg%ie_array(3))
      ! copied from subspace_diagonalization.f90
      mat1 = 0d0
      if(info%if_divide_orbit) then
        do ik=ik_s,ik_e
        do ispin = 1, nspin
          do io1 = 1, no
            if (io_s<= io1 .and. io1 <= io_e) then
              call copy_data(psi1%zwf(:, :, :, ispin, io1, ik, 1),wf_io1)
            end if
            call comm_bcast(wf_io1, info%icomm_o, info%irank_io(io1))
            do io2 = 1, no
              if (io_s<= io2 .and. io2 <= io_e) then
                cbox = 0d0
                !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+:cbox)
                do iz=is(3),ie(3)
                do iy=is(2),ie(2)
                do ix=is(1),ie(1)
                  cbox = cbox + conjg(wf_io1(ix,iy,iz)) * psi2%zwf(ix,iy,iz,ispin,io2,ik,1)
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
            cbox = cbox + conjg(psi1%zwf(ix,iy,iz,ispin,io1,ik,1)) * psi2%zwf(ix,iy,iz,ispin,io2,ik,1)
          end do
          end do
          end do
          mat1(io1,io2,ispin,ik) = cbox * system%hvol
        end do
        end do
        end do
        end do
      end if

      call comm_summation(mat1,mat,no**2*nspin*nk,info%icomm_rko)

    end subroutine inner_product
  end subroutine projection_gs


  ! Trial: num. of gs orbital is taken from restart file (full gsorbitals) 
  !        while that of rt orbital is from nstate in input (or Nelec/2 for temperature<0)
  subroutine projection_gx(itt,ofl,dt,mg,system,info,stencil,ppg,psi_t,tpsi,ttpsi,srg,energy,rt)
    use structures
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global
    use salmon_global, only: ncg,nelec,nscf
    use inputoutput, only: t_unit_time
    use subspace_diagonalization, only: ssdg
    use gram_schmidt_orth, only: gram_schmidt
    use Conjugate_Gradient, only: gscg_zwf
    use Total_Energy, only: calc_eigen_energy
    use sendrecv_grid
    implicit none
    integer                 ,intent(in) :: itt
    type(s_ofile)           ,intent(in) :: ofl
    real(8)                 ,intent(in) :: dt
    type(s_rgrid)           ,intent(in) :: mg
    type(s_dft_system)                  :: system
    type(s_parallel_info)               :: info
    type(s_stencil)         ,intent(in) :: stencil
    type(s_pp_grid)         ,intent(in) :: ppg
    type(s_orbital)         ,intent(in) :: psi_t ! | u_{n,k}(t) >
    type(s_orbital)                     :: tpsi,ttpsi ! temporary arrays
    type(s_sendrecv_grid)               :: srg
    type(s_dft_energy)                  :: energy
    type(s_rt)                          :: rt
    !
    integer :: nspin,no,nk,ik_s,ik_e,io_s,io_e,is(3),ie(3)
    integer :: ix,iy,iz,io1,io2,io,ik,ispin,iter_GS,niter
    complex(8),dimension(rt%no0,system%no,system%nspin,system%nk) :: mat
    real(8) :: coef(rt%no0,system%nk,system%nspin)
    real(8) :: nee, neh, wspin, dE
    complex(8) :: cbox
    !
    integer :: mo
    real(8) :: rocc_rt_bk(1:system%no,1:system%nk,1:system%nspin)
    integer,dimension(2,3) :: neig
     
    if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ projection"

    if(nscf==0) then
      niter = 10
    else
      niter = nscf
    end if

    ! switch no,io_s,io_e,numo : rt => gs
    call norb_rt2gs_for_proj(system,info,rt)

    ! change size of rocc : rt => gs
    call change_size_rt2gs_for_proj


    ! GS iteration
    ! rt%tpsi0 = | u_{n,k+A(t)/c} >, the ground-state wavefunction whose k-point is shifted by A(t).
    call calc_eigen_energy(energy,rt%tpsi0,tpsi,ttpsi,system,info,mg,rt%vloc0,stencil,srg,ppg)
    dE = energy%E_kin - rt%E_old
    if(abs(dE) < 1e-12) then
      if(comm_is_root(nproc_id_global)) write(*,*) "projection: already converged, E_kin(new)-E_kin(old)=",dE
    else
      do iter_GS=1,niter
        call ssdg(mg,system,info,stencil,rt%tpsi0,tpsi,ppg,rt%vloc0,srg)
        call gscg_zwf(ncg,mg,system,info,stencil,ppg,rt%vloc0,srg,rt%tpsi0,rt%cg)
        call gram_schmidt(system, mg, info, rt%tpsi0)
        call calc_eigen_energy(energy,rt%tpsi0,tpsi,ttpsi,system,info,mg,rt%vloc0,stencil,srg,ppg)
        dE = energy%E_kin - rt%E_old
        if(comm_is_root(nproc_id_global)) write(*,'(a,i6,e20.10)') "projection: ",iter_GS,dE
        if(abs(dE) < 1e-6) exit
        rt%E_old = energy%E_kin
      end do
    end if


    ! switch back no,io_s,io_e,numo : gs => rt
    call norb_gs2rt_for_proj(system,info,rt)

    ! change back size of rocc : gs => rt
    call change_size_gs2rt_for_proj

    nspin = system%nspin
    no = system%no
    nk = system%nk
    is = mg%is
    ie = mg%ie
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e
    if(nspin==1) then
      wspin = 2d0
    else if(nspin==2) then
      wspin = 1d0
    end if
    mo = rt%no0

    call inner_product(rt%tpsi0,psi_t,mat) ! mat(n,m) = < u_{n,k+A(t)/c} | u_{m,k}(t) >
    
    coef=0.d0
    do ispin=1,nspin
    do ik=1,nk
    do io1=1,mo
      do io2=1,no
        coef(io1,ik,ispin) = coef(io1,ik,ispin) &
        & + system%rocc(io2,ik,ispin)*system%wtk(ik)* abs(mat(io1,io2,ispin,ik))**2
      end do
    end do
    end do
    end do

    nee = 0d0
    neh = dble(nelec)
    do ispin=1,nspin
    do ik=1,nk
    do io=1,mo
      nee = nee + ((wspin-rt%rocc0(io,ik,ispin))/wspin) * coef(io,ik,ispin)
      neh = neh - (rt%rocc0(io,ik,ispin)/wspin) * coef(io,ik,ispin)
    end do
    end do
    end do
   !nee  = sum(ovlp_occ(NBoccmax+1:NB,:))
   !neh  = sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
    if(comm_is_root(nproc_id_global))then
      write(ofl%fh_nex,'(99(1X,E23.15E3))') dble(itt)*dt*t_unit_time%conv, nee, neh
      write(ofl%fh_ovlp,'(i11)') itt
      do ispin=1,nspin
      do ik=1,nk
        write(ofl%fh_ovlp,'(i6,5000(1X,E23.15E3))') ik,(coef(io,ik,ispin)*nk,io=1,mo)
      end do
      end do
    end if

    return
    
  contains


    subroutine change_size_rt2gs_for_proj
      implicit none

      !$omp workshare
      rocc_rt_bk(:,:,:) = system%rocc(:,:,:)
      !$omp end workshare
      deallocate(system%rocc)
      allocate(system%rocc(1:system%no,1:system%nk,1:system%nspin)) !size of gs
      !$omp workshare
      system%rocc = rt%rocc0
      !$omp end workshare
      deallocate(energy%esp,tpsi%zwf,ttpsi%zwf)
      deallocate(tpsi%zwf_real,tpsi%zwf_imag,ttpsi%zwf_real,ttpsi%zwf_imag)
      allocate(energy%esp(system%no,system%nk,system%nspin))
      call allocate_orbital_complex(system%nspin,mg,info, tpsi)
      call allocate_orbital_complex(system%nspin,mg,info,ttpsi)
      ! sendrecv_grid object for wavefunction updates
      call create_sendrecv_neig(neig, info) ! neighboring node array
      call init_sendrecv_grid(srg, mg, info%numo*info%numk*system%nspin, info%icomm_rko, neig)
      call dealloc_cache(srg)

    end subroutine change_size_rt2gs_for_proj

    subroutine change_size_gs2rt_for_proj
      implicit none

      !$omp workshare
      rt%rocc0 = system%rocc 
      !$omp end workshare
      deallocate(system%rocc)
      allocate(system%rocc(1:system%no,1:system%nk,1:system%nspin)) !size of rt
      !$omp workshare
      system%rocc(:,:,:) = rocc_rt_bk(:,:,:)
      !$omp end workshare
      deallocate(energy%esp,tpsi%zwf,ttpsi%zwf)
      deallocate(tpsi%zwf_real,tpsi%zwf_imag,ttpsi%zwf_real,ttpsi%zwf_imag)
      allocate(energy%esp(system%no,system%nk,system%nspin))
      call allocate_orbital_complex(system%nspin,mg,info, tpsi)
      call allocate_orbital_complex(system%nspin,mg,info,ttpsi)
      ! sendrecv_grid object for wavefunction updates
      call create_sendrecv_neig(neig, info) ! neighboring node array
      call init_sendrecv_grid(srg, mg, info%numo*info%numk*system%nspin, info%icomm_rko, neig)
      call dealloc_cache(srg)

    end subroutine change_size_gs2rt_for_proj

  
    subroutine inner_product(psi1,psi2,mat)
      use communication, only: comm_summation, comm_bcast
      use pack_unpack, only: copy_data
      implicit none
      type(s_orbital), intent(in) :: psi1,psi2
      complex(8),dimension(rt%no0,system%no,system%nspin,system%nk) :: mat
      complex(8),dimension(rt%no0,system%no,system%nspin,system%nk) :: mat1
      complex(8) :: wf_io1(mg%is_array(1):mg%ie_array(1) &
                        & ,mg%is_array(2):mg%ie_array(2) &
                        & ,mg%is_array(3):mg%ie_array(3))
      ! copied from subspace_diagonalization.f90
      mat1 = 0d0
      if(info%if_divide_orbit) then
        do ik=ik_s,ik_e
        do ispin = 1, nspin
          do io1 = 1, mo
            if (rt%io_s0<= io1 .and. io1 <= rt%io_e0) then
              call copy_data(psi1%zwf(:, :, :, ispin, io1, ik, 1),wf_io1)
            end if
            call comm_bcast(wf_io1, info%icomm_o, rt%irank_io0(io1))
            do io2 = 1, no
              if (io_s<= io2 .and. io2 <= io_e) then
                cbox = 0d0
                !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+:cbox)
                do iz=is(3),ie(3)
                do iy=is(2),ie(2)
                do ix=is(1),ie(1)
                  cbox = cbox + conjg(wf_io1(ix,iy,iz)) * psi2%zwf(ix,iy,iz,ispin,io2,ik,1)
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
        do io1=rt%io_s0,rt%io_e0
        do io2=io_s,io_e
          cbox = 0d0
          do iz=is(3),ie(3)
          do iy=is(2),ie(2)
          do ix=is(1),ie(1)
            cbox = cbox + conjg(psi1%zwf(ix,iy,iz,ispin,io1,ik,1)) * psi2%zwf(ix,iy,iz,ispin,io2,ik,1)
          end do
          end do
          end do
          mat1(io1,io2,ispin,ik) = cbox * system%hvol
        end do
        end do
        end do
        end do
      end if

      call comm_summation(mat1,mat,mo*no*nspin*nk,info%icomm_rko)

    end subroutine inner_product

  end subroutine projection_gx

!===================================================================================================================================

end module write_sub
