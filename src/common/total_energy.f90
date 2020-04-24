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

  SUBROUTINE calc_Total_Energy_isolated(system,info,mg,pp,rho,Vh,Vxc,rion_update,energy)
    use structures
    use salmon_global, only: kion
    use communication, only: comm_summation
    use timer
    implicit none
    type(s_dft_system)      ,intent(in) :: system
    type(s_parallel_info)   ,intent(in) :: info
    type(s_rgrid)           ,intent(in) :: mg
    type(s_pp_info)         ,intent(in) :: pp
    type(s_scalar)          ,intent(in) :: rho(system%Nspin),Vh,Vxc(system%Nspin)
    logical                 ,intent(in) :: rion_update
    type(s_dft_energy)                  :: energy
    !
    integer :: io,ik,ispin,Nspin
    integer :: ix,iy,iz,ia,ib
    real(8) :: sum1,sum2,Eion,Etot,r

    call timer_begin(LOG_TE_ISOLATED_CALC)

    Nspin = system%Nspin

    if (Rion_update) then
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
      energy%E_ion_ion = Eion
    end if

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
!$omp          shared(Nspin,mg,Vh,rho,Vxc)
    do ispin=1,Nspin
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        sum1 = sum1 - 0.5d0* Vh%f(ix,iy,iz) * rho(ispin)%f(ix,iy,iz)    &
                    - ( Vxc(ispin)%f(ix,iy,iz) * rho(ispin)%f(ix,iy,iz) )
      end do
      end do
      end do
    end do
!$omp end parallel do
    call timer_end(LOG_TE_ISOLATED_CALC)

    call timer_begin(LOG_TE_ISOLATED_COMM_COLL)

    call comm_summation(sum1,sum2,info%icomm_r)

    Etot = Etot + sum2*system%Hvol + energy%E_xc + energy%E_ion_ion

    energy%E_tot = Etot

    call timer_end(LOG_TE_ISOLATED_COMM_COLL)

    return
  end SUBROUTINE calc_Total_Energy_isolated

!===================================================================================================================================

  SUBROUTINE calc_Total_Energy_periodic(mg,ewald,system,info,pp,ppg,fg,poisson,rion_update,energy)
    use structures
    use salmon_math
    use math_constants,only : pi,zi
    use salmon_global, only: kion,aEwald, cutoff_r
    use communication, only: comm_summation,comm_is_root
    use timer
    implicit none
    type(s_rgrid)           ,intent(in) :: mg
    type(s_ewald_ion_ion)   ,intent(in) :: ewald
    type(s_dft_system)      ,intent(in) :: system
    type(s_parallel_info)   ,intent(in) :: info
    type(s_pp_info)         ,intent(in) :: pp
    type(s_pp_grid)         ,intent(in) :: ppg
    type(s_reciprocal_grid) ,intent(in) :: fg
    type(s_poisson)         ,intent(in) :: poisson
    logical                 ,intent(in) :: rion_update
    type(s_dft_energy)                  :: energy
    !
    integer :: ix,iy,iz,iia,ia,ib,zps1,zps2,ipair
    real(8) :: rr,rab(3),r(3),E_tmp,E_tmp_l,g(3),Gd,sysvol,E_wrk(5),E_sum(5)
    real(8) :: etmp
    complex(8) :: rho_e,rho_i

    call timer_begin(LOG_TE_PERIODIC_CALC)

    sysvol = system%det_a

    E_tmp = 0d0
    E_tmp_l = 0d0
    if (rion_update) then ! Ewald sum

      if(ewald%yn_bookkeep=='y') then

!$omp parallel do private(iia,ia,ipair,ix,iy,iz,ib,r,rab,rr) reduction(+:E_tmp)
         do iia=1,info%nion_mg
        !do ia=1,system%nion
            ia = info%ia_mg(iia)
            do ipair = 1,ewald%npair_bk(iia)
               ix = ewald%bk(1,ipair,iia)
               iy = ewald%bk(2,ipair,iia)
               iz = ewald%bk(3,ipair,iia)
               ib = ewald%bk(4,ipair,iia)
              !if (ix**2+iy**2+iz**2 == 0 .and. ia == ib) cycle
               r(1) = ix*system%primitive_a(1,1) &
                    + iy*system%primitive_a(1,2) &
                    + iz*system%primitive_a(1,3)
               r(2) = ix*system%primitive_a(2,1) &
                    + iy*system%primitive_a(2,2) &
                    + iz*system%primitive_a(2,3)
               r(3) = ix*system%primitive_a(3,1) &
                    + iy*system%primitive_a(3,2) &
                    + iz*system%primitive_a(3,3)
               rab(1) = system%Rion(1,ia)-r(1) - system%Rion(1,ib)
               rab(2) = system%Rion(2,ia)-r(2) - system%Rion(2,ib)
               rab(3) = system%Rion(3,ia)-r(3) - system%Rion(3,ib)
               rr = sum(rab(:)**2)
               if(rr .gt. cutoff_r**2) cycle
               E_tmp = E_tmp + 0.5d0*pp%Zps(Kion(ia))*pp%Zps(Kion(ib))*erfc_salmon(sqrt(aEwald*rr))/sqrt(rr)

            end do  !ipair
         end do     !ia
!$omp end parallel do


      else
         stop "must use book-keeping method for periodic condition"

      endif

!$omp parallel do collapse(2) default(none) &
!$omp          reduction(+:E_tmp_l) &
!$omp          private(ix,iy,iz,rho_i) &
!$omp          shared(fg,aEwald,sysvol,mg,ppg)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        rho_i = ppg%zrhoG_ion(ix,iy,iz)
        E_tmp_l = E_tmp_l + sysvol* fg%coef(ix,iy,iz) * (abs(rho_i)**2 * fg%exp_ewald(ix,iy,iz) *0.5d0) ! ewald (--> Rion_update)
      end do
      end do
      end do
!$omp end parallel do
    end if

    etmp = 0d0
    E_wrk = 0d0
!$omp parallel do collapse(2) default(none) &
!$omp          reduction(+:E_wrk,etmp) &
!$omp          private(ix,iy,iz,g,rho_i,rho_e,ia,r,Gd) &
!$omp          shared(mg,fg,aEwald,system,sysvol,kion,poisson,ppg,info)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      g(1) = fg%vec_G(1,ix,iy,iz)
      g(2) = fg%vec_G(2,ix,iy,iz)
      g(3) = fg%vec_G(3,ix,iy,iz)
      
      rho_e = poisson%zrhoG_ele(ix,iy,iz)
      rho_i = ppg%zrhoG_ion(ix,iy,iz)
      
      E_wrk(1) = E_wrk(1) + sysvol* fg%coef(ix,iy,iz) * (abs(rho_e)**2*0.5d0)     ! Hartree
      E_wrk(2) = E_wrk(2) + sysvol* fg%coef(ix,iy,iz) * (-rho_e*conjg(rho_i))     ! electron-ion (valence)

      do ia=info%ia_s,info%ia_e
        r = system%Rion(1:3,ia)
        Gd = g(1)*r(1) + g(2)*r(2) + g(3)*r(3)
        etmp = etmp + conjg(rho_e)*ppg%zVG_ion(ix,iy,iz,Kion(ia))*exp(-zI*Gd)  ! electron-ion (core)
      end do
    end do
    end do
    end do
!$omp end parallel do
    call timer_end(LOG_TE_PERIODIC_CALC)

    call timer_begin(LOG_TE_PERIODIC_COMM_COLL)
    
    call comm_summation(etmp,E_wrk(3),info%icomm_ko) ! for atom index #ia

    if (rion_update) then
      E_wrk(4) = E_tmp_l
      E_wrk(5) = E_tmp
      call comm_summation(E_wrk,E_sum,5,info%icomm_r)

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
      call comm_summation(E_wrk,E_sum,3,info%icomm_r)
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
    type(s_dft_energy)                     :: energy
    type(s_orbital)                        :: tpsi,htpsi,ttpsi
    type(s_dft_system)      ,intent(in)    :: system
    type(s_parallel_info)   ,intent(in)    :: info
    type(s_rgrid)           ,intent(in)    :: mg
    type(s_scalar)          ,intent(in)    :: V_local(system%Nspin)
    type(s_stencil)         ,intent(in)    :: stencil
    type(s_sendrecv_grid)   ,intent(inout) :: srg
    type(s_pp_grid)         ,intent(in)    :: ppg
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
!$omp parallel do collapse(2) default(none) &
!$omp          private(ik,io) &
!$omp          shared(info,wrk1,tpsi,htpsi,system,is,ie,ispin,im)
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          wrk1(io,ik) = sum( tpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
                        * htpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) * system%Hvol
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
                      * sum(  tpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
                           * ttpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) * system%Hvol
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
            * sum( tpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
              * (htpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
                - (ttpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
                   + V_local(ispin)%f(is(1):ie(1),is(2):ie(2),is(3):ie(3)) &
                   * tpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
                ) &
              ) &
            )

        end do
        end do
      end do
!$omp end parallel do
      E_local(2) = E_tmp
      call timer_end(LOG_EIGEN_ENERGY_CALC)
      
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

    end if
    
    call timer_begin(LOG_EIGEN_ENERGY_COMM_COLL)
    call comm_summation(E_local,E_sum,2,info%icomm_rko)

    energy%E_kin      = E_sum(1)
    energy%E_ion_nloc = E_sum(2)
    call timer_end(LOG_EIGEN_ENERGY_COMM_COLL)

    deallocate(wrk1,wrk2)
    return
  End Subroutine calc_eigen_energy
  
!===================================================================================================================================

  subroutine init_ewald(system,info,ewald)
    use structures
    use salmon_math
!    use math_constants,only : pi,zi
    use salmon_global, only: NEwald,aEwald, cutoff_r,cutoff_r_buff, cutoff_g
    use communication, only: comm_is_root,comm_summation,comm_get_groupinfo
    use parallelization, only: nproc_id_global
    use inputoutput, only: au_length_aa
    use timer
    implicit none
    type(s_dft_system) ,intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_ewald_ion_ion) :: ewald
    !
    integer :: ix,iy,iz,iia,ia,ib,ir,ipair   !,ig
    integer :: npair_bk_max, npair_bk_loc
   !integer :: k,irank, nproc, ig_tmp,ig_sum
    real(8) :: rr,rab(3),r(3) !,g(3),G2
    real(8) :: r1, cutoff_erfc_r, tmp

    !(find cut off length)
    cutoff_erfc_r = 1d-10*au_length_aa  !cut-off threshold of erfc(ar)/r [1/bohr]
    if(cutoff_r .lt. 0d0) then
       do ir=1,100
          r1=dble(ir)/au_length_aa  ![bohr]
          tmp = erfc_salmon(sqrt(aEwald)*r1)/r1
          if(tmp .le. cutoff_erfc_r) then
             cutoff_r = r1
             exit
          endif
       enddo
    endif
    if(cutoff_g .lt. 0d0) cutoff_g = 99d99 ![1/Bohr]   !cutoff in G space

    if(comm_is_root(nproc_id_global)) then
       write(*,900) " == Ewald =="
       write(*,800) " cutoff length in real-space in ewald =", cutoff_r*au_length_aa, " [A]"
       write(*,800) " (buffer length in bookkeeping =", cutoff_r_buff*au_length_aa, " [A])"
       write(*,810) " cutoff length in G-space in ewald =", cutoff_g/au_length_aa, " [1/A]"
    endif

800 format(a,f6.2,a)
810 format(a,e18.10,a)
900 format(a)

    !Book-keeping in ewald(ion-ion)

    !(check maximum number of pairs and allocate)
    npair_bk_max = 0
!$omp parallel do private(iia,ia,ix,iy,iz,ib,r,rab,rr,npair_bk_loc) &
!$omp             reduction(max:npair_bk_max)
    do iia=1,info%nion_mg
   !do ia=1,system%nion
       ia = info%ia_mg(iia)
       npair_bk_loc = 0
       do ix=-NEwald,NEwald
       do iy=-NEwald,NEwald
       do iz=-NEwald,NEwald
          do ib=1,system%nion
             if (ix**2+iy**2+iz**2 == 0 .and. ia == ib) cycle
             r(1) = ix*system%primitive_a(1,1) + &
                    iy*system%primitive_a(1,2) + &
                    iz*system%primitive_a(1,3)
             r(2) = ix*system%primitive_a(2,1) + &
                    iy*system%primitive_a(2,2) + &
                    iz*system%primitive_a(2,3)
             r(3) = ix*system%primitive_a(3,1) + &
                    iy*system%primitive_a(3,2) + &
                    iz*system%primitive_a(3,3)
             rab(1) = system%Rion(1,ia)-r(1) - system%Rion(1,ib)
             rab(2) = system%Rion(2,ia)-r(2) - system%Rion(2,ib)
             rab(3) = system%Rion(3,ia)-r(3) - system%Rion(3,ib)
             rr = sum(rab(:)**2)
             if(rr .le. (cutoff_r+cutoff_r_buff)**2) then
                npair_bk_loc = npair_bk_loc + 1
             endif
          end do
        end do
        end do
        end do
        npair_bk_max = max(npair_bk_max,npair_bk_loc)
      end do
!$omp end parallel do

      ewald%nmax_pair_bk = npair_bk_max
      ewald%nmax_pair_bk = nint(ewald%nmax_pair_bk * 1.5d0)
      allocate( ewald%bk(4,ewald%nmax_pair_bk,info%nion_mg) )
      allocate( ewald%npair_bk(info%nion_mg) )

      if(comm_is_root(nproc_id_global)) then
         write(*,820) " number of ion-ion pair(/atom) used for allocation of bookkeeping=", ewald%nmax_pair_bk
         write(*,*)"==========="
820      format(a,i6)
      endif

!$omp parallel do private(iia,ia,ipair,ix,iy,iz,ib,r,rab,rr)
    do iia=1,info%nion_mg
   !do ia=1,system%nion
       ia = info%ia_mg(iia)
       ipair = 0
       do ix=-NEwald,NEwald
       do iy=-NEwald,NEwald
       do iz=-NEwald,NEwald
          do ib=1,system%nion
             if (ix**2+iy**2+iz**2 == 0 .and. ia == ib) cycle
             r(1) = ix*system%primitive_a(1,1) + &
                    iy*system%primitive_a(1,2) + &
                    iz*system%primitive_a(1,3)
             r(2) = ix*system%primitive_a(2,1) + &
                    iy*system%primitive_a(2,2) + &
                    iz*system%primitive_a(2,3)
             r(3) = ix*system%primitive_a(3,1) + &
                    iy*system%primitive_a(3,2) + &
                    iz*system%primitive_a(3,3)
             rab(1) = system%Rion(1,ia)-r(1) - system%Rion(1,ib)
             rab(2) = system%Rion(2,ia)-r(2) - system%Rion(2,ib)
             rab(3) = system%Rion(3,ia)-r(3) - system%Rion(3,ib)
             rr = sum(rab(:)**2)
             if(rr .le. cutoff_r**2) then
                ipair = ipair + 1
                ewald%bk(1,ipair,iia) = ix
                ewald%bk(2,ipair,iia) = iy
                ewald%bk(3,ipair,iia) = iz
                ewald%bk(4,ipair,iia) = ib
             endif
          end do
        end do
        end do
        end do
        ewald%npair_bk(iia) = ipair
      end do
!$omp end parallel do

      return
      !xxxxxxxxxxxxxx
      !currently, following part is under construction

!      !for G-space
!      ig_sum = 0
!      ig_tmp = 0
!      do ig=fg%ig_s,fg%ig_e
!         if(ig == fg%iGzero ) cycle
!         g(1) = fg%Gx(ig)
!         g(2) = fg%Gy(ig)
!         g(3) = fg%Gz(ig)
!         G2   = sum(g(:)**2)
!         if(G2 .gt. cutoff_g**2) cycle
!         ig_tmp = ig_tmp + 1
!      enddo
!      !this cause MPI communitation error ---- why??
!      call comm_summation(ig_tmp,ig_sum,1,nproc_group_global)
!      ewald%ng_bk = ig_sum
!
!      call comm_get_groupinfo(fg%icomm_G,irank,nproc)
!
!      if(nproc .le. ewald%ng_bk) then
!         k = mod(ewald%ng_bk,nproc)
!         if(k==0) then
!            ewald%ng_r = ewald%ng_bk / nproc
!         else
!            ewald%ng_r = ewald%ng_bk / nproc + 1
!         endif
!         ewald%ng_s = ewald%ng_r * irank + 1
!         ewald%ng_e = ewald%ng_s + ewald%ng_r - 1
!         if (irank == nproc-1) ewald%ng_e = ewald%ng_bk
!         if (ewald%ng_e .gt. ewald%ng_bk) ewald%ng_e = -1
!         if (ewald%ng_s .gt. ewald%ng_bk) then
!            ewald%ng_s =  0
!            ewald%ng_e = -1
!         endif
!
!      else
!         if(irank+1.le.ewald%ng_bk) then
!            ewald%ng_s = irank + 1
!            ewald%ng_e = ewald%ng_s
!         else
!            ewald%ng_s = 0
!            ewald%ng_e = -1
!         endif
!      endif
!
!      if(comm_is_root(nproc_id_global)) &
!
!      write(*,'(a,i8)') " number of G-points in ewald", ewald%ng_bk
!
!      if(comm_is_root(nproc_id_global)) then
!         write(*,*) "  #irank, ng_s, ng_e"
!      endif
!      write(*,'(3i6)')  irank, ewald%ng_s, ewald%ng_e
!

  end subroutine init_ewald

!===================================================================================================================================

  function check_rion_update() result(rion_update)
    use salmon_global, only: theory,yn_opt,yn_md
    implicit none
    logical :: rion_update

    select case(theory)
    case('dft','dft_band','dft_md')
      rion_update = (yn_opt == 'y' .or. theory == 'dft_md')
    case('tddft_response','tddft_pulse','single_scale_maxwell_tddft','multi_scale_maxwell_tddft')
      rion_update = (yn_md == 'y')
    case default
      rion_update = .false.
    end select
  end function

END MODULE salmon_Total_Energy
