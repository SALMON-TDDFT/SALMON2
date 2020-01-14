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

CONTAINS

!===================================================================================================================================

  SUBROUTINE calc_Total_Energy_isolated(energy,system,info,ng,pp,rho,Vh,Vxc)
    use structures
    use salmon_global, only: kion
    use communication, only: comm_summation
    use timer
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_orbital_parallel),intent(in) :: info
    type(s_rgrid)     ,intent(in) :: ng
    type(s_pp_info)   ,intent(in) :: pp
    type(s_scalar)    ,intent(in) :: rho(system%Nspin),Vh,Vxc(system%Nspin)
    type(s_dft_energy)            :: energy
    !
    integer :: io,ik,ispin,Nspin
    integer :: ix,iy,iz,ia,ib
    real(8) :: sum1,sum2,Eion,Etot,r

    call timer_begin(LOG_TE_ISOLATED_CALC)

    Nspin = system%Nspin

!    if (Rion_update) then
      Eion = 0d0
!$omp parallel do default(none) &
!$omp          reduction(+:Eion) &
!$omp          private(ia,ib,r) &
!$omp          shared(system,pp,Kion)
      do ia=1,system%nion
        do ib=1,ia-1
          r = sqrt((system%Rion(1,ia)-system%Rion(1,ib))**2      &
                  +(system%Rion(2,ia)-system%Rion(2,ib))**2      &
                  +(system%Rion(3,ia)-system%Rion(3,ib))**2)
          Eion = Eion + pp%Zps(Kion(ia)) * pp%Zps(Kion(ib)) /r
        end do
      end do
!$omp end parallel do
!    end if

    Etot = 0d0
!$omp parallel do collapse(3) default(none) &
!$omp          reduction(+:Etot) &
!$omp          private(ispin,ik,io) &
!$omp          shared(Nspin,system,energy)
    do ispin=1,Nspin
    do ik=1,system%nk
    do io=1,system%no
      Etot = Etot + system%rocc(io,ik,ispin) * energy%esp(io,ik,ispin)
    end do
    end do
    end do
!$omp end parallel do

    sum1 = 0d0
!$omp parallel do collapse(4) default(none) &
!$omp          reduction(+:sum1) &
!$omp          private(ispin,ix,iy,iz) &
!$omp          shared(Nspin,ng,Vh,rho,Vxc)
    do ispin=1,Nspin
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        sum1 = sum1 - 0.5d0* Vh%f(ix,iy,iz) * rho(ispin)%f(ix,iy,iz)    &
                    - ( Vxc(ispin)%f(ix,iy,iz) * rho(ispin)%f(ix,iy,iz) )
      end do
      end do
      end do
    end do
!$omp end parallel do
    call timer_end(LOG_TE_ISOLATED_CALC)

    call timer_begin(LOG_TE_ISOLATED_COMM_COLL)

    call comm_summation(sum1,sum2,info%icomm_rko)

    Etot = Etot + sum2*system%Hvol + energy%E_xc + Eion

    energy%E_ion_ion = Eion
    energy%E_tot = Etot

    call timer_end(LOG_TE_ISOLATED_COMM_COLL)

    return
  end SUBROUTINE calc_Total_Energy_isolated

!===================================================================================================================================

  SUBROUTINE calc_Total_Energy_periodic(energy,system,pp,fg,rion_update)
    use structures
    use salmon_math
    use math_constants,only : pi,zi
    use salmon_global, only: kion,NEwald,aEwald
    use communication, only: comm_summation,comm_get_groupinfo
    use timer
    implicit none
    type(s_dft_system) ,intent(in) :: system
    type(s_pp_info)    ,intent(in) :: pp
    type(s_reciprocal_grid),intent(in) :: fg
    type(s_dft_energy)             :: energy
    logical,intent(in)             :: rion_update
    !
    integer :: ix,iy,iz,ia,ib,ig,zps1,zps2
    real(8) :: rr,rab(3),r(3),E_tmp,E_tmp_l,g(3),G2,Gd,sysvol,E_wrk(5),E_sum(5)
    real(8) :: etmp
    complex(8) :: rho_e,rho_i
    integer :: irank,nproc,nion_r,nion_s,nion_e

    call timer_begin(LOG_TE_PERIODIC_CALC)

    sysvol = system%det_a

    E_tmp = 0d0
    E_tmp_l = 0d0
    if (rion_update) then ! Ewald sum
      call comm_get_groupinfo(fg%icomm_G,irank,nproc)
      nion_r = (system%nion + 1) / nproc
      nion_s = nion_r * irank + 1
      nion_e = nion_s + nion_r - 1
      if (irank == nproc-1) nion_e = system%nion
!$omp parallel do collapse(4) default(none) &
!$omp          reduction(+:E_tmp) &
!$omp          private(ia,ib,ix,iy,iz,r,rab,rr) &
!$omp          shared(NEwald,system,pp,Kion,aEwald,nion_s,nion_e)
      do ia=1,system%nion
      do ix=-NEwald,NEwald
      do iy=-NEwald,NEwald
      do iz=-NEwald,NEwald
        do ib=nion_s,nion_e
          !if (ix**2+iy**2+iz**2 == 0 .and. ia == ib) cycle ! iwata
          r(1) = ix*system%primitive_a(1,1) + iy*system%primitive_a(1,2) + iz*system%primitive_a(1,3)
          r(2) = ix*system%primitive_a(2,1) + iy*system%primitive_a(2,2) + iz*system%primitive_a(2,3)
          r(3) = ix*system%primitive_a(3,1) + iy*system%primitive_a(3,2) + iz*system%primitive_a(3,3)
          rab(1) = system%Rion(1,ia)-r(1) - system%Rion(1,ib)
          rab(2) = system%Rion(2,ia)-r(2) - system%Rion(2,ib)
          rab(3) = system%Rion(3,ia)-r(3) - system%Rion(3,ib)
          rr = sum(rab(:)**2) ; if ( rr == 0.0d0 ) cycle ! iwata 
          E_tmp = E_tmp + 0.5d0*pp%Zps(Kion(ia))*pp%Zps(Kion(ib))*erfc_salmon(sqrt(aEwald*rr))/sqrt(rr)
        end do
      end do
      end do
      end do
      end do
!$omp end parallel do

!$omp parallel do default(none) &
!$omp          reduction(+:E_tmp_l) &
!$omp          private(ig,g,G2,rho_i) &
!$omp          shared(fg,aEwald,sysvol)
      do ig=fg%ig_s,fg%ig_e
        if(ig == fg%iGzero ) cycle
        g(1) = fg%Gx(ig)
        g(2) = fg%Gy(ig)
        g(3) = fg%Gz(ig)
        G2 = g(1)**2 + g(2)**2 + g(3)**2
        rho_i = fg%zrhoG_ion(ig)
        E_tmp_l = E_tmp_l + sysvol*(4*Pi/G2)*(abs(rho_i)**2*exp(-G2/(4*aEwald))*0.5d0) ! ewald (--> Rion_update)
      end do
!$omp end parallel do
    end if

    E_wrk = 0d0
!$omp parallel do default(none) &
!$omp          reduction(+:E_wrk) &
!$omp          private(ig,g,G2,rho_i,rho_e,ia,r,Gd,etmp) &
!$omp          shared(fg,aEwald,system,sysvol,kion)
    do ig=fg%ig_s,fg%ig_e
      g(1) = fg%Gx(ig)
      g(2) = fg%Gy(ig)
      g(3) = fg%Gz(ig)
      rho_e = fg%zrhoG_ele(ig)

      if (ig /= fg%iGzero) then
        G2 = g(1)**2 + g(2)**2 + g(3)**2
        rho_i = fg%zrhoG_ion(ig)
        E_wrk(1) = E_wrk(1) + sysvol*(4*Pi/G2)*(abs(rho_e)**2*0.5d0)     ! Hartree
        E_wrk(2) = E_wrk(2) + sysvol*(4*Pi/G2)*(-rho_e*conjg(rho_i))     ! electron-ion (valence)
      end if

      etmp = 0d0

#if _OPENMP >= 201307
!$omp simd reduction(+:etmp)
#endif
      do ia=1,system%nion
        r = system%Rion(1:3,ia)
        Gd = g(1)*r(1) + g(2)*r(2) + g(3)*r(3)
        etmp = etmp + conjg(rho_e)*fg%zdVG_ion(ig,Kion(ia))*exp(-zI*Gd)  ! electron-ion (core)
      end do
      E_wrk(3) = E_wrk(3) + etmp
    enddo
!$omp end parallel do
    call timer_end(LOG_TE_PERIODIC_CALC)

    call timer_begin(LOG_TE_PERIODIC_COMM_COLL)

    if (rion_update) then
      E_wrk(4) = E_tmp_l
      E_wrk(5) = E_tmp
      call comm_summation(E_wrk,E_sum,5,fg%icomm_G)

  ! ion-ion energy
      zps1 = 0
      zps2 = 0
!$omp parallel do default(none) private(ia) shared(system,pp,Kion) reduction(+:zps1,zps2)
      do ia=1,system%nion
        zps1 = zps1 + pp%Zps(Kion(ia))
        zps2 = zps2 + pp%Zps(Kion(ia))**2
      end do
!$omp end parallel do

      E_sum(5) = E_sum(5) - Pi*zps1**2/(2*aEwald*sysvol) - sqrt(aEwald/Pi)*zps2
      energy%E_ion_ion = E_sum(5) + E_sum(4)
    else
      call comm_summation(E_wrk,E_sum,3,fg%icomm_G)
    end if

  ! Hartree energy
    energy%E_h = E_sum(1)

  ! electron-ion energy (local part)
    energy%E_ion_loc = E_sum(2) + E_sum(3)

  ! total energy
    energy%E_tot = energy%E_kin + energy%E_h + energy%E_ion_loc + energy%E_ion_nloc + energy%E_xc + energy%E_ion_ion

    call timer_end(LOG_TE_PERIODIC_COMM_COLL)

    return
  end SUBROUTINE calc_Total_Energy_periodic

!===================================================================================================================================

! eigen energies (esp), kinetic energy (E_kin), & nonlocal part of electron-ion energy (E_ion_nloc)
  Subroutine calc_eigen_energy(energy,tpsi,htpsi,ttpsi,system,info,mg,V_local,stencil,srg,ppg)
    use structures
    use communication, only: comm_summation
    use hamiltonian, only: hpsi
    use spin_orbit_global, only: SPIN_ORBIT_ON
    use timer
    implicit none
    type(s_dft_energy)         :: energy
    type(s_orbital)            :: tpsi,htpsi,ttpsi
    type(s_dft_system),intent(in) :: system
    type(s_orbital_parallel),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_scalar) ,intent(in) :: V_local(system%Nspin)
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_pp_grid),intent(in) :: ppg
    !
    integer :: ik,io,ispin,im,nk,no,is(3),ie(3),Nspin
    real(8) :: E_tmp,E_local(2),E_sum(2)
    real(8),allocatable :: wrk1(:,:),wrk2(:,:)

    call timer_begin(LOG_EIGEN_ENERGY_CALC)
    if(info%im_s/=1 .or. info%im_e/=1) stop "error: calc_eigen_energy"
    im = 1

    Nspin = system%Nspin
    is = mg%is
    ie = mg%ie
    no = system%no
    nk = system%nk
    allocate(wrk1(no,nk),wrk2(no,nk))
    wrk1 = 0d0
    call timer_end(LOG_EIGEN_ENERGY_CALC)

    call timer_begin(LOG_EIGEN_ENERGY_HPSI)
    call hpsi(tpsi,htpsi,info,mg,V_local,system,stencil,srg,ppg,ttpsi)
    call timer_end(LOG_EIGEN_ENERGY_HPSI)

    if(allocated(tpsi%rwf)) then
      do ispin=1,Nspin
        call timer_begin(LOG_EIGEN_ENERGY_CALC)
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          wrk1(io,ik) = sum( tpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
                        * htpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) * system%Hvol
        end do
        end do
        call timer_end(LOG_EIGEN_ENERGY_CALC)

        call timer_begin(LOG_EIGEN_ENERGY_COMM_COLL)
        call comm_summation(wrk1,wrk2,no*nk,info%icomm_rko)
        energy%esp(:,:,ispin) = wrk2
        call timer_end(LOG_EIGEN_ENERGY_COMM_COLL)
      end do
    else
    ! eigen energies (esp)
      do ispin=1,Nspin
        call timer_begin(LOG_EIGEN_ENERGY_CALC)
!$omp parallel do collapse(2) default(none) &
!$omp          private(ik,io) &
!$omp          shared(info,wrk1,tpsi,htpsi,system,is,ie,ispin,im)
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          wrk1(io,ik) = sum( conjg( tpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) &
                                 * htpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) * system%Hvol
        end do
        end do
!$omp end parallel do
        call timer_end(LOG_EIGEN_ENERGY_CALC)

        call timer_begin(LOG_EIGEN_ENERGY_COMM_COLL)
        call comm_summation(wrk1,wrk2,no*nk,info%icomm_rko)
        energy%esp(:,:,ispin) = wrk2
        call timer_end(LOG_EIGEN_ENERGY_COMM_COLL)
      end do

      call timer_begin(LOG_EIGEN_ENERGY_CALC)
      if ( SPIN_ORBIT_ON ) then
        energy%esp(:,:,1) = energy%esp(:,:,1) + energy%esp(:,:,2)
        energy%esp(:,:,2) = energy%esp(:,:,1)
      end if

    ! kinetic energy (E_kin)
      E_tmp = 0d0
!$omp parallel do collapse(3) default(none) &
!$omp          reduction(+:E_tmp) &
!$omp          private(ispin,ik,io) &
!$omp          shared(Nspin,info,tpsi,ttpsi,system,is,ie,im)
      do ispin=1,Nspin
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          E_tmp = E_tmp + system%rocc(io,ik,ispin)*system%wtk(ik) &
                      * sum( conjg( tpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) &
                                 * ttpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) * system%Hvol
        end do
        end do
      end do
!$omp end parallel do
      E_local(1) = E_tmp

    ! nonlocal part (E_ion_nloc)
      E_tmp = 0d0
!$omp parallel do collapse(3) default(none) &
!$omp          reduction(+:E_tmp) &
!$omp          private(ispin,ik,io) &
!$omp          shared(Nspin,info,tpsi,htpsi,ttpsi,system,is,ie,im,V_local)
      do ispin=1,Nspin
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e

          E_tmp = E_tmp + system%rocc(io,ik,ispin)*system%wtk(ik) * system%hvol &
            * sum( conjg(tpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im)) &
              * (htpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
                - (ttpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
                   + V_local(ispin)%f(is(1):ie(1),is(2):ie(2),is(3):ie(3)) &
                   * tpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
                ) &
              ) &
            )

        end do
        end do
      end do
!$omp end parallel do
      E_local(2) = E_tmp
      call timer_end(LOG_EIGEN_ENERGY_CALC)

      call timer_begin(LOG_EIGEN_ENERGY_COMM_COLL)
      call comm_summation(E_local,E_sum,2,info%icomm_rko)

      energy%E_kin      = E_sum(1)
      energy%E_ion_nloc = E_sum(2)
      call timer_end(LOG_EIGEN_ENERGY_COMM_COLL)

    end if

    deallocate(wrk1,wrk2)
    return
  End Subroutine calc_eigen_energy

!===================================================================================================================================

  function check_rion_update() result(rion_update)
    use salmon_global, only: theory,yn_opt,yn_md
    implicit none
    logical :: rion_update

    select case(theory)
    case('dft','dft_band','dft_md')
      rion_update = (yn_opt == 'y' .or. theory == 'dft_md')
    case('tddft_response','tddft_pulse','single_scale_maxwell_tddft','multiscale_experiment')
      rion_update = (yn_md == 'y')
    case default
      rion_update = .false.
    end select
  end function

END MODULE salmon_Total_Energy
