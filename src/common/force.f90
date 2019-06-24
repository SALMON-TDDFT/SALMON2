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
module force_sub
  use math_constants,only : pi,zi
  implicit none

contains

!===================================================================================================================================

  subroutine calc_force_salmon(force,system,pp,fg,info,mg,stencil,srg,ppg,tpsi)
    use structures
    use hpsi_sub
    use stencil_sub
    use update_overlap_sub
    use salmon_communication, only: comm_summation
    implicit none
    type(s_system) ,intent(in) :: system
    type(s_pp_info),intent(in) :: pp
    type(s_fourier_grid),intent(in) :: fg
    type(s_wf_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(in) :: srg
    type(s_pp_grid),intent(in) :: ppg
    type(s_wavefunction)       :: tpsi
    type(s_force)              :: force
    !
    integer :: ix,iy,iz,ia,nion,im,Nspin,ik_s,ik_e,io_s,io_e,norb,iorb,nlma,ik,io,ispin,ilma,j
    real(8) :: kAc(3)
    real(8),allocatable :: F_tmp(:,:),F_sum(:,:)
    complex(8) :: w(3),uVpsi,duVpsi(3)
    complex(8),allocatable :: gtpsi(:,:,:,:),uVpsibox(:,:),uVpsibox2(:,:)

    nion = system%nion
    if(.not.allocated(force%F)) allocate(force%F(3,nion))
    allocate(F_tmp(3,nion),F_sum(3,nion))

    if(info%im_s/=1 .or. info%im_e/=1) stop "error: calc_force_periodic" !??????
    im = 1

    Nspin = system%Nspin
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e
    norb = Nspin* info%numo * info%numk * info%numm

    Nlma = ppg%Nlma

  ! Ewald sum of ion-ion interaction
    call force_ion_ion(F_sum,F_tmp,system,pp,fg,nion)
    force%F = F_sum

  ! electron-ion interaction

    if(allocated(tpsi%rwf)) then
      allocate(tpsi%zwf(mg%is_array(1):mg%ie_array(1) &
                       ,mg%is_array(2):mg%ie_array(2) &
                       ,mg%is_array(3):mg%ie_array(3) &
                       ,nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e))
      tpsi%zwf = cmplx(tpsi%rwf)
    end if

    allocate(gtpsi(3,mg%is_array(1):mg%ie_array(1) &
                    ,mg%is_array(2):mg%ie_array(2) &
                    ,mg%is_array(3):mg%ie_array(3)))

  ! uVpsibox2 = < uV | exp(ikr) | psi >
    allocate(uVpsibox(Nlma,Norb),uVpsibox2(Nlma,Norb))
    uVpsibox = 0d0
    iorb = 0
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      iorb = iorb + 1
      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + conjg(ppg%ekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
        uVpsibox(ilma,iorb) = uVpsi
      end do
    end do
    end do
    end do
    call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,info%icomm_r)

    if(info%if_divide_rspace) then
      call update_overlap_C(tpsi%zwf,mg%is_array,mg%ie_array,norb,4 & !????????
                           ,mg%is,mg%ie,info%irank_r,info%icomm_r)
      !call update_overlap_complex8(srg, mg, tpsi%rwf)
    end if

    kAc = 0d0
    F_tmp = 0d0
    iorb = 0
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      iorb = iorb + 1

    ! gtpsi = (nabla) psi
      call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,io,ik,im),gtpsi,mg%is_array,mg%ie_array,mg%is,mg%ie &
          ,mg%idx,mg%idy,mg%idz,stencil%nabt,stencil%matrix_B)

    ! local part
      do ia=1,nion
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          w = conjg(gtpsi(:,ix,iy,iz)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
          F_tmp(:,ia) = F_tmp(:,ia) - 2d0*info%occ(io,ik,ispin) * dble(w(:))* ppg%Vpsl_atom(ix,iy,iz,ia) * system%Hvol
        end do
        end do
        end do
      end do

    ! nonlocal part
      if(allocated(stencil%kAc)) kAc(1:3) = stencil%kAc(ik,1:3)
      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        duVpsi = 0d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          w = gtpsi(:,ix,iy,iz) + zI* kAc(:) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
          duVpsi = duVpsi + conjg(ppg%ekr_uV(j,ilma,ik)) * w ! < uV | exp(ikr) (nabla) | psi >
        end do
        F_tmp(:,ia) = F_tmp(:,ia) - 2d0*info%occ(io,ik,ispin) * dble( conjg(duVpsi(:)) * uVpsibox2(ilma,iorb) ) * system%Hvol
      end do

    end do
    end do
    end do
    call comm_summation(F_tmp,F_sum,3*nion,info%icomm_rko)
    force%F = force%F + F_sum

    if(allocated(tpsi%rwf)) deallocate(tpsi%zwf)
    deallocate(F_tmp,F_sum,gtpsi,uVpsibox,uVpsibox2)
    return
  end subroutine calc_force_salmon

  subroutine force_ion_ion(F_sum,F_tmp,system,pp,fg,nion)
    use structures
    use salmon_math
    use salmon_global, only: kion,NEwald,aEwald
    use salmon_communication, only: comm_summation
    implicit none
    type(s_system) ,intent(in) :: system
    type(s_pp_info),intent(in) :: pp
    type(s_fourier_grid),intent(in) :: fg
    integer,intent(in) :: nion
    real(8) :: F_tmp(3,nion),F_sum(3,nion)
    !
    integer :: ix,iy,iz,ia,ib,ig
    real(8) :: rr,rab(3),r(3),g(3),G2,Gd
    complex(8) :: rho_i

    select case(system%iperiodic)
    case(0)

      F_sum = 0d0
      do ia=1,nion
        do ib=1,nion
          if(ia==ib) cycle
          rab(:) = system%Rion(:,ia) - system%Rion(:,ib)
          rr = sum(rab(:)**2)
          F_sum(:,ia) = F_sum(:,ia) + pp%Zps(Kion(ia))*pp%Zps(Kion(ib)) * rab(:)/sqrt(rr)**3
        end do
      end do

    case(3)
      F_tmp = 0d0
      do ia=1,nion
        r = system%Rion(1:3,ia)
        do ig=fg%ig_s,fg%ig_e
          if(ig == fg%iGzero ) cycle
          g(1) = fg%Gx(ig)
          g(2) = fg%Gy(ig)
          g(3) = fg%Gz(ig)
          G2 = g(1)**2 + g(2)**2 + g(3)**2
          Gd = g(1)*r(1) + g(2)*r(2) + g(3)*r(3)
          rho_i = fg%rhoG_ion(ig)
          F_tmp(:,ia) = F_tmp(:,ia) + g(:)*(4*Pi/G2)*exp(-G2/(4*aEwald))*pp%Zps(Kion(ia)) &
                                      *zI*0.5d0*(conjg(rho_i)*exp(-zI*Gd)-rho_i*exp(zI*Gd))
        end do
      end do
      call comm_summation(F_tmp,F_sum,3*nion,fg%icomm_fourier)

      F_tmp = 0d0
      do ia=1,nion
        do ix=-NEwald,NEwald
        do iy=-NEwald,NEwald
        do iz=-NEwald,NEwald
          do ib=1,nion
            if (ix**2+iy**2+iz**2 == 0 .and. ia == ib) cycle
            r(1) = ix*system%al(1,1) + iy*system%al(1,2) + iz*system%al(1,3)
            r(2) = ix*system%al(2,1) + iy*system%al(2,2) + iz*system%al(2,3)
            r(3) = ix*system%al(3,1) + iy*system%al(3,2) + iz*system%al(3,3)
            rab(1) = system%Rion(1,ia)-r(1) - system%Rion(1,ib)
            rab(2) = system%Rion(2,ia)-r(2) - system%Rion(2,ib)
            rab(3) = system%Rion(3,ia)-r(3) - system%Rion(3,ib)
            rr = sum(rab(:)**2)
            F_tmp(:,ia) = F_tmp(:,ia) - pp%Zps(Kion(ia))*pp%Zps(Kion(ib))*rab(:)/sqrt(rr)*(-erfc_salmon(sqrt(aEwald*rr))/rr &
                                        -2*sqrt(aEwald/(rr*Pi))*exp(-aEwald*rr))
          end do
        end do
        end do
        end do
      end do
      F_sum = F_sum + F_tmp
    end select
  end subroutine force_ion_ion
end module force_sub
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
