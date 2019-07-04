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
module print_sub
  use math_constants,only : zi
  implicit none

contains

!--------------------------------------------------------------------------------
!! export SYSNAME_k.data file
  subroutine write_k_data(k_rd,system,stencil)
    use structures
    use salmon_global, only: sysname
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root,comm_sync_all
    use salmon_file, only: open_filehandle
    implicit none
    type(s_system) ,intent(in) :: system
    type(s_stencil),intent(in) :: stencil
    real(8),intent(in) :: k_rd(3,system%nk)
    !
    integer :: fh_k
    integer :: ik,NK
    real(8) :: wk
    character(100) :: file_k_data

    NK = system%nk
    file_k_data = trim(sysname)//'_k.data'!??????
    wk = 1.0 !??????? wk=1 only (symmetry weight)

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
            & k_rd(1,ik) / system%brl(1,1), &
            & k_rd(2,ik) / system%brl(2,2), &
            & k_rd(3,ik) / system%brl(3,3), &
            & wk !??????? wk=1 only (symmetry weight)
        end do !ik
      else
        write(fh_k, '("#",1X,A)') "k-point distribution (nonorthogonal coordinate)"
        write(fh_k, '("#",1X,A,":",1X,A)') "brl", "reciprocal primitive vectors"
        write(fh_k, '("#",1X,A,":",1X,A)') "ik", "k-point index"
        write(fh_k, '("#",1X,A,":",1X,A)') "kx,ky,kz", "k-vectors"
        write(fh_k, '("#",1X,A,":",1X,A)') "wk", "Weight of k-point"
        write(fh_k, '("#",A,"[",A,"]")') "brl(1:3,1:3)", "a.u."
        write(fh_k, '(9(1X,E23.15E3))') system%brl(1:3,1:3)
        write(fh_k, '("#",99(1X,I0,":",A,"[",A,"]"))') &
          & 1, "ik", "none", &
          & 2, "kx", "a.u.", &
          & 3, "ky", "a.u.", &
          & 4, "kz", "a.u.", &
          & 5, "wk", "none"
        do ik = 1, NK
          write(fh_k, '(I6,99(1X,E23.15E3))') &
            & ik, &
            & k_rd(1,ik), &
            & k_rd(2,ik), &
            & k_rd(3,ik), &
            & wk !??????? wk=1 only (symmetry weight)
        end do !ik
      end if
      close(fh_k)
    end if
    call comm_sync_all
  end subroutine write_k_data

!--------------------------------------------------------------------------------
!! export SYSNAME_tm.data file
  subroutine write_tm_data(tpsi,system,info,mg,stencil,srg,ppg)
    use structures
    use stencil_sub
    use sendrecv_grid
    use salmon_global, only: sysname
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root,comm_summation,comm_sync_all
    use salmon_file, only: open_filehandle
    implicit none
    type(s_system) ,intent(in) :: system
    type(s_wf_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_pp_grid),intent(in) :: ppg
    type(s_wavefunction)       :: tpsi
    !
    integer :: fh_tm
    integer :: i,j,ik,ib,ib1,ib2,ilma,nlma,ia,ix,iy,iz,NB,NK,im,ispin,ik_s,ik_e,is(3),ie(3)
    real(8) :: x,y,z
    complex(8),allocatable :: upu(:,:,:,:),upu_l(:,:,:,:),gtpsi(:,:,:,:)
    complex(8),allocatable :: uVpsi(:,:),uVpsix(:,:),uVpsiy(:,:),uVpsiz(:,:)
    complex(8),allocatable :: uVpsixx(:,:),uVpsixy(:,:),uVpsixz(:,:)
    complex(8),allocatable :: uVpsiyy(:,:),uVpsiyz(:,:),uVpsizz(:,:)
    complex(8),allocatable :: u_rVnl_Vnlr_u(:,:,:,:),u_rVnl_Vnlr_u_l(:,:,:,:)
    complex(8) :: u_rVnl_u(3),u_Vnlr_u(3),veik
    complex(8),allocatable ::  u_rVnlr_Vnlrr_u(:,:,:,:),u_rVnlr_Vnlrr_u_l(:,:,:,:)
    complex(8) :: ctmp1,ctmp2,wrk(3)
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
      write(*,*) "error @ write_tm_data: do not use orbital parallelization"
      return
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

    allocate(gtpsi(3,mg%is_array(1):mg%ie_array(1) &
                    ,mg%is_array(2):mg%ie_array(2) &
                    ,mg%is_array(3):mg%ie_array(3)))
  ! overlap region communication
    if(info%if_divide_rspace) then
      !call update_overlap_C(tpsi%zwf,mg%is_array,mg%ie_array,norb,Nd & !????????
      !                     ,mg%is,mg%ie,info%irank_r,info%icomm_r)
      call update_overlap_complex8(srg, mg, tpsi%zwf)
    end if

    do ik=ik_s,ik_e
    do ib2=1,NB

    ! gtpsi = (nabla) psi
      call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,ib2,ik,im),gtpsi,mg%is_array,mg%ie_array,mg%is,mg%ie &
          ,mg%idx,mg%idy,mg%idz,stencil%nabt,stencil%matrix_B)
      do ib1=1,NB
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
          veik = conjg(ppg%ekr_uV(j,ilma,ik))
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


    !calculate <u_nk|[r_j,dVnl^(0)]r|u_nk>  (j=x,y,z)
    u_rVnlr_Vnlrr_u_l(:,:,:,:) = 0d0
    do ik=ik_s,ik_e
    do ilma=1,Nlma
       ia=ppg%ia_tbl(ilma)
       uVpsi=0d0;  uVpsix=0d0;  uVpsiy=0d0;  uVpsiz=0d0
       uVpsixx=0d0;  uVpsixy=0d0;  uVpsixz=0d0
                     uVpsiyy=0d0;  uVpsiyz=0d0
                                   uVpsizz=0d0
       do j=1,ppg%Mps(ia)
          x = ppg%Rxyz(1,j,ia)
          y = ppg%Rxyz(2,j,ia)
          z = ppg%Rxyz(3,j,ia)
          ix = ppg%Jxyz(1,j,ia)
          iy = ppg%Jxyz(2,j,ia)
          iz = ppg%Jxyz(3,j,ia)
          veik = conjg(ppg%ekr_uV(j,ilma,ik))
          do ib=1,NB
          uVpsi(  ib,ik)=uVpsi(  ib,ik)+veik*    tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik|u>
          uVpsix( ib,ik)=uVpsix( ib,ik)+veik* x *tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*x|u>
          uVpsiy( ib,ik)=uVpsiy( ib,ik)+veik* y *tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*y|u>
          uVpsiz( ib,ik)=uVpsiz( ib,ik)+veik* z *tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*z|u>
          uVpsixx(ib,ik)=uVpsixx(ib,ik)+veik*x*x*tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*xx|u>
          uVpsixy(ib,ik)=uVpsixy(ib,ik)+veik*x*y*tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*xy|u>
          uVpsixz(ib,ik)=uVpsixz(ib,ik)+veik*x*z*tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*xz|u>
          uVpsiyy(ib,ik)=uVpsiyy(ib,ik)+veik*y*y*tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*yy|u>
          uVpsiyz(ib,ik)=uVpsiyz(ib,ik)+veik*y*z*tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*yz|u>
          uVpsizz(ib,ik)=uVpsizz(ib,ik)+veik*z*z*tpsi%zwf(ix,iy,iz,ispin,ib,ik,im) !=<v|e^ik*zz|u>
          enddo

       end do

       do ib=1,NB
          !xx
          ctmp1 = conjg(uVpsix(ib,ik))*uVpsix( ib,ik)
          ctmp2 = conjg(uVpsi( ib,ik))*uVpsixx(ib,ik)
          u_rVnlr_Vnlrr_u_l(1,1,ib,ik) = &
          u_rVnlr_Vnlrr_u_l(1,1,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
          !xy
          ctmp1 = conjg(uVpsix(ib,ik))*uVpsiy( ib,ik)
          ctmp2 = conjg(uVpsi( ib,ik))*uVpsixy(ib,ik)
          u_rVnlr_Vnlrr_u_l(1,2,ib,ik) = &
          u_rVnlr_Vnlrr_u_l(1,2,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
          !xz
          ctmp1 = conjg(uVpsix(ib,ik))*uVpsiz( ib,ik)
          ctmp2 = conjg(uVpsi( ib,ik))*uVpsixz(ib,ik)
          u_rVnlr_Vnlrr_u_l(1,3,ib,ik) = &
          u_rVnlr_Vnlrr_u_l(1,3,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
          !yx
          ctmp1 = conjg(uVpsiy(ib,ik))*uVpsix( ib,ik)
          ctmp2 = conjg(uVpsi( ib,ik))*uVpsixy(ib,ik)
          u_rVnlr_Vnlrr_u_l(2,1,ib,ik) = &
          u_rVnlr_Vnlrr_u_l(2,1,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
          !yy
          ctmp1 = conjg(uVpsiy(ib,ik))*uVpsiy( ib,ik)
          ctmp2 = conjg(uVpsi( ib,ik))*uVpsiyy(ib,ik)
          u_rVnlr_Vnlrr_u_l(2,2,ib,ik) = &
          u_rVnlr_Vnlrr_u_l(2,2,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
          !yz
          ctmp1 = conjg(uVpsiy(ib,ik))*uVpsiz( ib,ik)
          ctmp2 = conjg(uVpsi( ib,ik))*uVpsiyz(ib,ik)
          u_rVnlr_Vnlrr_u_l(2,3,ib,ik) = &
          u_rVnlr_Vnlrr_u_l(2,3,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
          !zx
          ctmp1 = conjg(uVpsiz(ib,ik))*uVpsix( ib,ik)
          ctmp2 = conjg(uVpsi( ib,ik))*uVpsixz(ib,ik)
          u_rVnlr_Vnlrr_u_l(3,1,ib,ik) = &
          u_rVnlr_Vnlrr_u_l(3,1,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
          !zy
          ctmp1 = conjg(uVpsiz(ib,ik))*uVpsiy( ib,ik)
          ctmp2 = conjg(uVpsi( ib,ik))*uVpsiyz(ib,ik)
          u_rVnlr_Vnlrr_u_l(3,2,ib,ik) = &
          u_rVnlr_Vnlrr_u_l(3,2,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)
          !zz
          ctmp1 = conjg(uVpsiz(ib,ik))*uVpsiz( ib,ik)
          ctmp2 = conjg(uVpsi( ib,ik))*uVpsizz(ib,ik)
          u_rVnlr_Vnlrr_u_l(3,3,ib,ik) = &
          u_rVnlr_Vnlrr_u_l(3,3,ib,ik) + (ctmp1 - ctmp2)*ppg%rinv_uvu(ilma)

       enddo

    enddo  !ilma
    enddo  !ik
    call comm_summation(u_rVnlr_Vnlrr_u_l,u_rVnlr_Vnlrr_u,3*3*NB*NK,info%icomm_rko)

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

       !<u_mk|[r_j,dVnl^(0)]r_i|u_nk>  (j,i=x,y,z)
       write(fh_tm,*) "#<u_mk|[r_j,dVnl^(0)]r_i|u_nk>  (j,i=x,y,z)"
       do ik=1,NK
       do ib=1,NB
          do i=1,3
             write(fh_tm,9000) ik,ib,i,(u_rVnlr_Vnlrr_u(i,j,ib,ik),j=1,3)
          enddo
       enddo
       enddo

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

!--------------------------------------------------------------------------------
  subroutine write_xyz(comment,action,rvf,system,force)
  ! Write xyz in xyz format but also velocity and force are printed if necessary
  ! (these can be used for restart of opt and md)
    use structures, only: s_system, s_force
    use inputoutput, only: au_length_aa
    use salmon_global, only: SYSname,atom_name
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    implicit none

    type(s_system),intent(in) :: system
    type(s_force),intent(in) :: force

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
             write(unit_xyz,120) trim(atom_name(ia)),system%Rion(1:3,ia)*au_length_aa,system%Velocity(1:3,ia),force%F(1:3,ia)
          endif
       enddo

    else if(action=='end') then
       close(unit_xyz)
    endif

100 format(a2,3f18.10)
110 format(a2,3f18.10, "  #v=",3f18.10)
120 format(a2,3f18.10, "  #v=",3f18.10, "  #f=",3f18.10)

  end subroutine write_xyz

  subroutine write_rt_data_3d(it,ofl,iflag_md,dt)
    use structures, only: s_ofile
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    use salmon_file, only: open_filehandle
    use inputoutput, only: t_unit_time,t_unit_current,t_unit_ac,t_unit_elec
    implicit none
    type(s_ofile) :: ofl
    integer, intent(in) :: it,iflag_md
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
       if(iflag_md==1) then
          write(uid,10) "Jm", "Matter current density(electrons)"
          write(uid,10) "Jmi","Matter current density(ions)"
       else
          write(uid,10) "Jm", "Matter current density"
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
         & 13, "E_tot_z", trim(t_unit_elec%name), &
         & 14, "Jm_x", trim(t_unit_current%name), &
         & 15, "Jm_y", trim(t_unit_current%name), &
         & 16, "Jm_z", trim(t_unit_current%name)
       if(iflag_md==1) then
          write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
               & 17, "Jmi_x", trim(t_unit_current%name), &
               & 18, "Jmi_y", trim(t_unit_current%name), &
               & 19, "Jmi_z", trim(t_unit_current%name)
       endif
       write(uid,*)
       flush(uid)

    else  !it>=0
       uid = ofl%fh_rt
       write(uid,"(F16.8,a)",advance='no') &   !just temporal
          & it * dt * t_unit_time%conv,    &   !just temporal
          & "   under construction"            !just temporal
!        write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
!          & it * dt * t_unit_time%conv, &
!          & Ac_ext(it, 1:3) * t_unit_ac%conv, &
!          & E_ext(it, 1:3) * t_unit_elec%conv, &
!          & Ac_tot(it, 1:3) * t_unit_ac%conv, &
!          & E_tot(it, 1:3) * t_unit_elec%conv, &
!          & javt(it, 1:3) * t_unit_current%conv
!        if(iflag_md==1) then
!           write(uid, "(99(1X,E23.15E3))",advance='no') &
!                & javt_ion(it, 1:3) * t_unit_current%conv
!        endif
       write(uid,*)
    endif
    endif

  end subroutine

  subroutine write_rt_energy_data(it,ofl,iflag_md,dt,energy,md)
    use structures, only: s_ofile,s_energy,s_md
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    use salmon_global, only: ensemble, thermostat
    use salmon_file, only: open_filehandle
    use inputoutput, only: t_unit_time,t_unit_energy
    implicit none
    type(s_energy) :: energy
    type(s_md) :: md
    type(s_ofile) :: ofl
    integer, intent(in) :: it,iflag_md
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
       if(iflag_md==1) then
       write(uid,10) "Tion", "Kinetic energy of ions"
       write(uid,10) "Temperature_ion", "Temperature of ions"
       write(uid,10) "E_work", "Work energy of ions(sum f*dr)"
       if(ensemble=="NVT".and.thermostat=="nose-hoover")then
       write(uid,10) "Enh", "NH thermostat energy (MD)"
       write(uid,10) "Hnvt", "Hamiltonian with NH thermostat(MD)"
       write(uid,10) "Hnvt'","Hnvt using E_work"
       endif
       endif

       write(uid, '("#",99(1X,I0,":",A,"[",A,"]"))',advance='no') &
        & 1, "Time", trim(t_unit_time%name), &
        & 2, "Eall", trim(t_unit_energy%name), &
        & 3, "Eall-Eall0", trim(t_unit_energy%name)

       if(iflag_md==1) then
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

    else  !it>=0
       uid = ofl%fh_rt_energy

       write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
          & it * dt * t_unit_time%conv,        &
          & energy%E_tot * t_unit_energy%conv, &
          & energy%E_tot * t_unit_energy%conv     !just temporal
!          & (Eall_t(it) - Eall0) * t_unit_energy%conv
       if(iflag_md==1) then
       write(uid, "(99(1X,E23.15E3))",advance='no') &
          & md%Tene * t_unit_energy%conv, &
          & md%Temperature,               &
          & md%E_work * t_unit_energy%conv
       endif


!        write(uid, "(F16.8,99(1X,E23.15E3))",advance='no') &
!          & it * dt * t_unit_time%conv, &
!          & Eall_t(it) * t_unit_energy%conv, &
!          & (Eall_t(it) - Eall0) * t_unit_energy%conv
!        if(iflag_md==1) then
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

    endif
    endif

  end subroutine

end module print_sub
