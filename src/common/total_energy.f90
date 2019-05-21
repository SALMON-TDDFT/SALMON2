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
MODULE salmon_Total_Energy
implicit none
real(8),parameter :: pi=3.141592653589793d0 ! copied from salmon_math !?????????
complex(8),parameter :: zI=(0d0,1d0)

CONTAINS

!===================================================================================================================================

  SUBROUTINE calc_Total_Energy_isolated(energy,system,info,ng,pp,rho,Vh,Vxc)
    use structures
    use salmon_global, only: kion
    use salmon_communication, only: comm_summation
    implicit none
    type(s_system) ,intent(in) :: system
    type(s_wf_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: ng
    type(s_pp_info),intent(in) :: pp
    type(s_scalar) ,intent(in) :: rho(system%Nspin),Vh,Vxc(system%Nspin)
    type(s_energy)             :: energy
    !
    integer :: io,ik,ispin,Nspin
    integer :: ix,iy,iz,is,ia,ib
    real(8) :: sum1,sum2,Eion,Etot,r
    Nspin = system%Nspin

!    if (Rion_update) then
      Eion = 0d0
      do ia=1,system%nion
        do ib=1,ia-1
          r = sqrt((system%Rion(1,ia)-system%Rion(1,ib))**2      &
                  +(system%Rion(2,ia)-system%Rion(2,ib))**2      &
                  +(system%Rion(3,ia)-system%Rion(3,ib))**2)
          Eion = Eion + pp%Zps(Kion(ia)) * pp%Zps(Kion(ib)) /r
        end do
      end do
!    end if

    Etot = 0d0
    do ispin=1,Nspin
    do ik=1,system%nk
    do io=1,system%no
      Etot = Etot + system%rocc(io,ik,ispin) * energy%esp(io,ik,ispin)
    end do
    end do
    end do

    sum1 = 0d0
    do is=1,Nspin
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        sum1 = sum1 - 0.5d0* Vh%f(ix,iy,iz) * rho(is)%f(ix,iy,iz)    &
                    - ( Vxc(is)%f(ix,iy,iz) * rho(is)%f(ix,iy,iz) )
      end do
      end do
      end do
    end do
    
    call comm_summation(sum1,sum2,info%icomm_r)

    Etot = Etot + sum2*system%Hvol + energy%E_xc + Eion

    energy%E_ion_ion = Eion
    energy%E_tot = Etot

    return
  end SUBROUTINE calc_Total_Energy_isolated

!===================================================================================================================================

  SUBROUTINE calc_Total_Energy_periodic(energy,system,pp,fg)
    use structures
    use salmon_math
    use salmon_global, only: kion,NEwald,aEwald
    use salmon_communication, only: comm_summation
    implicit none
    type(s_system) ,intent(in) :: system
    type(s_pp_info),intent(in) :: pp
    type(s_fourier_grid),intent(in) :: fg
    type(s_energy)             :: energy
    !
    integer :: ix,iy,iz,ia,ib,ig
    real(8) :: rr,rab(3),r(3),E_tmp,g(3),G2,Gd,sysvol,E_wrk(4),E_sum(4)
    complex(8) :: rho_e,rho_i

    sysvol = system%det_al

!    if (Rion_update) then ! Ewald sum
    E_tmp = 0d0
    do ia=1,system%nion
      do ix=-NEwald,NEwald
      do iy=-NEwald,NEwald
      do iz=-NEwald,NEwald
        do ib=1,system%nion
          if (ix**2+iy**2+iz**2 == 0 .and. ia == ib) cycle
          r(1) = ix*system%al(1,1) + iy*system%al(1,2) + iz*system%al(1,3)
          r(2) = ix*system%al(2,1) + iy*system%al(2,2) + iz*system%al(2,3)
          r(3) = ix*system%al(3,1) + iy*system%al(3,2) + iz*system%al(3,3)
          rab(1) = system%Rion(1,ia)-r(1) - system%Rion(1,ib)
          rab(2) = system%Rion(2,ia)-r(2) - system%Rion(2,ib)
          rab(3) = system%Rion(3,ia)-r(3) - system%Rion(3,ib)
          rr = sum(rab(:)**2)
          E_tmp = E_tmp + 0.5d0*pp%Zps(Kion(ia))*pp%Zps(Kion(ib))*erfc_salmon(sqrt(aEwald*rr))/sqrt(rr)
        end do
      end do
      end do
      end do
    end do
    E_tmp = E_tmp - Pi*sum(pp%Zps(Kion(:)))**2/(2*aEwald*sysvol)-sqrt(aEwald/Pi)*sum(pp%Zps(Kion(:))**2)
!    end if

    E_wrk = 0d0
    do ig=fg%ig_s,fg%ig_e
      if(ig == fg%iGzero ) cycle
      g(1) = fg%Gx(ig)
      g(2) = fg%Gy(ig)
      g(3) = fg%Gz(ig)
      G2 = g(1)**2 + g(2)**2 + g(3)**2
      rho_i = fg%rhoG_ion(ig)
      rho_e = fg%rhoG_elec(ig)
      E_wrk(1) = E_wrk(1) + sysvol*(4*Pi/G2)*(abs(rho_i)**2*exp(-G2/(4*aEwald))*0.5d0) ! ewald (--> Rion_update)
      E_wrk(2) = E_wrk(2) + sysvol*(4*Pi/G2)*(abs(rho_e)**2*0.5d0)                     ! Hartree
      E_wrk(3) = E_wrk(3) + sysvol*(4*Pi/G2)*(-rho_e*conjg(rho_i))                     ! electron-ion (valence)
      do ia=1,system%nion
        r = system%Rion(1:3,ia)
        Gd = g(1)*r(1) + g(2)*r(2) + g(3)*r(3)
        E_wrk(4) = E_wrk(4) + conjg(rho_e)*fg%dVG_ion(ig,Kion(ia))*exp(-zI*Gd)         ! electron-ion (core)
      end do
    enddo
    call comm_summation(E_wrk,E_sum,4,fg%icomm_fourier)

  ! ion-ion energy
    energy%E_ion_ion = E_tmp + E_sum(1)

  ! Hartree energy
    energy%E_h = E_sum(2)

  ! electron-ion energy (local part)
    energy%E_ion_loc = E_sum(3) + E_sum(4)

  ! total energy
    energy%E_tot = energy%E_kin + energy%E_h + energy%E_ion_loc + energy%E_ion_nloc + energy%E_xc + energy%E_ion_ion

    return
  end SUBROUTINE calc_Total_Energy_periodic

!===================================================================================================================================

! eigen energies (esp), kinetic energy (E_kin), & nonlocal part of electron-ion energy (E_ion_nloc)
  Subroutine calc_eigen_energy(energy,tpsi,htpsi,ttpsi,system,info,mg,V_local,stencil,srg,ppg)
    use structures
    use salmon_communication, only: comm_summation
    use hpsi_sub
    use pseudo_pt_sub, only: pseudo_C
    implicit none
    type(s_energy)             :: energy
    type(s_wavefunction)       :: tpsi,htpsi,ttpsi
    type(s_system) ,intent(in) :: system
    type(s_wf_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_scalar) ,intent(in) :: V_local(system%Nspin)
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_pp_grid),intent(in) :: ppg
    !
    integer :: ik,io,jo,ispin,im,nk,no,is(3),ie(3),Nspin
    real(8) :: E_kin,E_ion_nloc,E_tmp
    real(8),allocatable :: wrk1(:,:),wrk2(:,:)

    if(info%im_s/=1 .or. info%im_e/=1) stop "error: calc_eigen_energy"
    im = 1

 write(*,*) "gomi1", real(tpsi%rwf(-31,-31,34,1,1,1,1))
    Nspin = system%Nspin
    is = mg%is
    ie = mg%ie
    no = system%no
    nk = system%nk
    allocate(wrk1(no,nk),wrk2(no,nk))
    wrk1 = 0d0
    call hpsi(tpsi,htpsi,info,mg,V_local,Nspin,stencil,srg,ppg,ttpsi)

 write(*,*) "gomi2", real(tpsi%rwf(-31,-31,34,1,1,1,1))
    if(allocated(tpsi%rwf)) then

      do ispin=1,Nspin
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          jo = info%io_tbl(io)
          wrk1(jo,ik) = sum( tpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
                        * htpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) * system%Hvol
        end do
        end do
        call comm_summation(wrk1,wrk2,no*nk,info%icomm_rko)
        energy%esp(:,:,ispin) = wrk2
      end do
    else
    ! eigen energies (esp)
      do ispin=1,Nspin
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          jo = info%io_tbl(io)
          wrk1(jo,ik) = sum( conjg( tpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) &
                                 * htpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) * system%Hvol
        end do
        end do
        call comm_summation(wrk1,wrk2,no*nk,info%icomm_rko)
        energy%esp(:,:,ispin) = wrk2
      end do

    ! kinetic energy (E_kin)
      E_tmp = 0d0
      do ispin=1,Nspin
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          E_tmp = E_tmp + info%occ(io,ik,ispin) &
                      * sum( conjg( tpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) &
                                 * ttpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) * system%Hvol
        end do
        end do
      end do
      call comm_summation(E_tmp,E_kin,info%icomm_rko)

    ! nonlocal part (E_ion_nloc)
      ttpsi%zwf = 0d0
      call pseudo_C(tpsi,ttpsi,info,nspin,ppg)
      E_tmp = 0d0
      do ispin=1,Nspin
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          E_tmp = E_tmp + info%occ(io,ik,ispin) &
                      * sum( conjg( tpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) &
                                 * ttpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) * system%Hvol
        end do
        end do
      end do
      call comm_summation(E_tmp,E_ion_nloc,info%icomm_rko)

      energy%E_kin = E_kin
      energy%E_ion_nloc = E_ion_nloc

    end if

    deallocate(wrk1,wrk2)
    return
  End Subroutine calc_eigen_energy

END MODULE salmon_Total_Energy
