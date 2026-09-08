!
!  Copyright 2019-2024 SALMON developers
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

module dcdft_soi
  implicit none

contains

  subroutine init_dcdft_soi(dc,pp,mixing,ewald)
    use structures
    use dcdft, only: init_dcdft
    implicit none
    type(s_dcdft)        ,intent(inout) :: dc
    type(s_pp_info)      ,intent(inout) :: pp
    type(s_mixing)       ,intent(inout) :: mixing
    type(s_ewald_ion_ion),intent(inout) :: ewald

    call init_dcdft(dc,pp,mixing,ewald)
  end subroutine init_dcdft_soi

  subroutine ne2mu_dcdft_soi(mg,info,energy,spsi,dc,system)
    use structures
    use communication, only: comm_summation
    use salmon_global, only: temperature
    implicit none
    type(s_rgrid),        intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_dft_energy),   intent(in) :: energy
    type(s_orbital),      intent(in) :: spsi
    type(s_dcdft),        intent(in) :: dc
    type(s_dft_system),   intent(inout) :: system
    !
    real(8) :: wspin,emax,emin
    real(8) :: ne_each(system%no,system%nspin)
    real(8),dimension(system%no,system%nspin,dc%n_frag) :: rocc,esp,ne_frag_orb,wrk1,wrk2

    if(system%nspin==1) then
      wspin = 2d0
    else if(system%nspin==2) then
      wspin = 1d0
    else
      stop "ne2mu_dcdft_soi: unsupported nspin."
    end if

    call calc_ne_each
    wrk1 = 0d0
    wrk2 = 0d0
    if(info%id_rko==0) then
      wrk1(1:system%no,1:system%nspin,dc%i_frag) = energy%esp(1:system%no,1,1:system%nspin)
      wrk2(1:system%no,1:system%nspin,dc%i_frag) = ne_each(1:system%no,1:system%nspin)
    end if
    call comm_summation(wrk1,esp,        system%no*system%nspin*dc%n_frag,dc%icomm_tot)
    call comm_summation(wrk2,ne_frag_orb,system%no*system%nspin*dc%n_frag,dc%icomm_tot)

    emin = minval(esp)
    emax = maxval(esp)
    call ne2mu_core(dc%elec_num_tot,emax,emin,system%mu)

    system%rocc(1:system%no,1,1:system%nspin) = rocc(1:system%no,1:system%nspin,dc%i_frag)
    return

  contains
    subroutine calc_ne_each
      implicit none
      integer :: io,ispin,ix,iy,iz
      real(8) :: wrk(system%no,system%nspin)

      wrk = 0d0
      if(allocated(spsi%zwf)) then
        ! A SOI orbital is one two-component spinor.  The duplicated spin
        ! slots below keep the existing rocc/esp storage convention intact.
        do io=info%io_s,info%io_e
          do iz=mg%is(3),min(mg%ie(3),dc%nxyz_domain(3))
          do iy=mg%is(2),min(mg%ie(2),dc%nxyz_domain(2))
          do ix=mg%is(1),min(mg%ie(1),dc%nxyz_domain(1))
            do ispin=1,system%nspin
              wrk(io,1) = wrk(io,1) + (abs(spsi%zwf(ix,iy,iz,ispin,io,1,1))**2) * system%hvol
            end do
          end do
          end do
          end do
        end do
        do ispin=2,system%nspin
          wrk(:,ispin) = wrk(:,1)
        end do
      else if(allocated(spsi%rwf)) then
        stop "ne2mu_dcdft_soi: real wavefunction is invalid for SOI occupation."
      else
        stop "ne2mu_dcdft_soi: neither rwf nor zwf is allocated."
      end if
      call comm_summation(wrk,ne_each,system%no*system%nspin,info%icomm_rko)
    end subroutine calc_ne_each

    subroutine mu2ne(muin,neout)
      implicit none
      real(8), intent(in)  :: muin
      real(8), intent(out) :: neout
      integer :: i_frag,ispin,io
      real(8) :: fact
      neout=0d0

      if(temperature==0d0) then
        do i_frag=1,dc%n_frag
        do ispin=1,system%nspin
        do io=1,system%no
          fact = esp(io,ispin,i_frag) - muin
          if(fact > 0d0) then
            rocc(io,ispin,i_frag) = 0d0
          else
            rocc(io,ispin,i_frag) = wspin
          endif
          neout = neout + rocc(io,ispin,i_frag) * ne_frag_orb(io,ispin,i_frag)
        end do
        end do
        end do
      else
        do i_frag=1,dc%n_frag
        do ispin=1,system%nspin
        do io=1,system%no
          fact = (esp(io,ispin,i_frag)-muin)/temperature
          if(fact.ge.40.d0) then
            rocc(io,ispin,i_frag) = 0d0
          else
            rocc(io,ispin,i_frag) = wspin/(1d0 + exp(fact))
          endif
          neout = neout + rocc(io,ispin,i_frag) * ne_frag_orb(io,ispin,i_frag)
        end do
        end do
        end do
      end if
      ! SO mode: two spin components are duplicate copy in current DC occupation model.
      neout = neout*0.5d0
    end subroutine mu2ne

    subroutine ne2mu_core(nein,emax,emin,muout)
      implicit none
      real(8),intent(in)  :: nein,emax,emin
      real(8),intent(out) :: muout
      ! Parameters
      integer, parameter :: MAX_BRACKET = 50
      integer, parameter :: MAX_BISECT  = 1000
      real(8), parameter :: CONV_NE_REL = 1d-10
      real(8), parameter :: CONV_MU     = 1d-9
      real(8), parameter :: CONV_NE_ABS = 1d-9
      ! Locals
      integer :: nc,ii
      real(8) :: mu1,mu2,mu3,ne1,ne3,ne3o,diff_ne,diff_mu,muo,diff_ne2

      ! --- Step 1: find a valid bracket [mu1,mu2] such that ne(mu1) >= nein >= ne(mu2) ---
      mu1 = emin
      mu2 = emax
      call mu2ne(mu1,ne1)
      nc  = 0
      do while ((ne1 - nein) <= 0d0 .and. nc <= MAX_BRACKET)
        nc  = nc + 1
        mu1 = emin - 0.2d0*dble(nc)
        mu2 = emax + 0.2d0*dble(nc)
        call mu2ne(mu1,ne1)
      end do
      if ((ne1 - nein) <= 0d0) then
        write(*,'(1x,a,1x,es23.15,1x,a,1x,es23.15,1x,a,1x,i0)') &
          'ne2mu_dcdft_soi: failed to bracket chemical potential (ne1,nein,nc)=', ne1, ',', nein, ',', nc
        stop "ne2mu_dcdft_soi: bracket search failed in ne2mu_core."
      end if

      ! --- Step 2: bisection to converge mu ---
      mu3  = 0.5d0*(mu1+mu2)
      ne3  = 0d0
      ne3o = 0d0
      muo  = 0d0
      do ii=1,MAX_BISECT
        mu3 = mu1 + (mu2-mu1)*0.5d0
        call mu2ne(mu3,ne3)

        if (abs(ne3) > tiny(1d0)) then
          diff_ne = abs((ne3-ne3o)/ne3)
        else
          diff_ne = abs(ne3-ne3o)
        end if
        diff_ne2 = abs(nein-ne3)
        diff_mu  = abs(mu3 - muo)

        if (diff_ne < CONV_NE_REL .and. diff_mu < CONV_MU .and. diff_ne2 < CONV_NE_ABS) exit

        if ((ne1-nein)*(ne3-nein) > 0d0) then
          mu1 = mu3
          ne1 = ne3
        else
          mu2 = mu3
        end if
        ne3o = ne3
        muo  = mu3
      end do

      muout = mu3
      call mu2ne(muout,ne3)
    end subroutine ne2mu_core
  end subroutine ne2mu_dcdft_soi

  subroutine calc_rho_total_dcdft_soi(nspin,lg,mg,info,rho_s,dc)
    use structures
    use communication, only: comm_summation
    use noncollinear_module, only: get_den_mat_noncollinear, set_den_mat_noncollinear, &
      & rot_dm_noncollinear
    implicit none
    integer,              intent(in) :: nspin
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_parallel_info),intent(in) :: info
    type(s_scalar),       intent(in) :: rho_s(nspin)
    type(s_dcdft),        intent(inout) :: dc
    !
    integer :: ix,iy,iz,is,js,ix_tot,iy_tot,iz_tot
    real(8) :: rho_int
    complex(8),allocatable :: frg_tmp(:,:,:,:,:),frg(:,:,:,:,:)
    complex(8),allocatable :: tot_tmp(:,:,:,:,:),tot(:,:,:,:,:)

    if(nspin /= 2) stop "calc_rho_total_dcdft_soi: nspin must be 2."

    allocate(frg_tmp(lg%num(1),lg%num(2),lg%num(3),2,2))
    allocate(frg    (lg%num(1),lg%num(2),lg%num(3),2,2))
    allocate(tot_tmp(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),2,2))
    allocate(tot    (dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),2,2))

    frg_tmp = (0.0d0,0.0d0)
    call get_den_mat_noncollinear(mg, &
      & frg_tmp(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:2,1:2))
    call comm_summation(frg_tmp,frg,lg%num(1)*lg%num(2)*lg%num(3)*4,info%icomm_r)

    tot_tmp = (0.0d0,0.0d0)
    if(info%id_rko==0) then
      do js=1,2
      do is=1,2
      do iz=1,dc%nxyz_domain(3); iz_tot = dc%jxyz_tot(iz,3)
      do iy=1,dc%nxyz_domain(2); iy_tot = dc%jxyz_tot(iy,2)
      do ix=1,dc%nxyz_domain(1); ix_tot = dc%jxyz_tot(ix,1)
        tot_tmp(ix_tot,iy_tot,iz_tot,is,js) = frg(ix,iy,iz,is,js)
      end do
      end do
      end do
      end do
      end do
    end if
    call comm_summation(tot_tmp,tot,dc%lg_tot%num(1)*dc%lg_tot%num(2)*dc%lg_tot%num(3)*4,dc%icomm_tot)

    call set_den_mat_noncollinear(dc%mg_tot, &
      & tot(dc%mg_tot%is(1):dc%mg_tot%ie(1),dc%mg_tot%is(2):dc%mg_tot%ie(2), &
      &     dc%mg_tot%is(3):dc%mg_tot%ie(3),1:2,1:2))
    call rot_dm_noncollinear(dc%rho_tot_s,dc%system_tot,dc%mg_tot)

    if(dc%id_tot==0) then
      rho_int = sum(dble(tot(:,:,:,1,1)+tot(:,:,:,2,2)))*dc%system_tot%hvol
      write(*,*) "integral(rho_tot)=",rho_int," Ne=",dc%elec_num_tot
    end if

    deallocate(frg_tmp,frg,tot_tmp,tot)
  end subroutine calc_rho_total_dcdft_soi

  subroutine calc_vlocal_fragment_dcdft_soi(nspin,mg,vloc,dc)
    use structures
    use dcdft, only: calc_vlocal_fragment_dcdft
    use communication, only: comm_summation
    use noncollinear_module, only: get_vxc_mat_noncollinear, set_vxc_mat_noncollinear
    implicit none
    integer,      intent(in) :: nspin
    type(s_rgrid),intent(in) :: mg
    type(s_scalar),intent(inout) :: vloc(nspin)
    type(s_dcdft), intent(inout) :: dc
    !
    integer :: ix,iy,iz,is,js,ix_tot,iy_tot,iz_tot
    complex(8),allocatable :: tot_tmp(:,:,:,:,:),tot(:,:,:,:,:),frg(:,:,:,:,:)

    if(nspin /= 2) stop "calc_vlocal_fragment_dcdft_soi: nspin must be 2."
    call calc_vlocal_fragment_dcdft(nspin,mg,vloc,dc)

    allocate(tot_tmp(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),2,2))
    allocate(tot    (dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),2,2))
    allocate(frg(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),2,2))

    tot_tmp = (0.0d0,0.0d0)
    call get_vxc_mat_noncollinear(dc%mg_tot, &
      & tot_tmp(dc%mg_tot%is(1):dc%mg_tot%ie(1),dc%mg_tot%is(2):dc%mg_tot%ie(2), &
      &         dc%mg_tot%is(3):dc%mg_tot%ie(3),1:2,1:2))
    call comm_summation(tot_tmp,tot,dc%lg_tot%num(1)*dc%lg_tot%num(2)*dc%lg_tot%num(3)*4,dc%icomm_tot)

    do js=1,2
    do is=1,2
    do iz=mg%is(3),mg%ie(3); iz_tot = dc%jxyz_tot(iz,3)
    do iy=mg%is(2),mg%ie(2); iy_tot = dc%jxyz_tot(iy,2)
    do ix=mg%is(1),mg%ie(1); ix_tot = dc%jxyz_tot(ix,1)
      frg(ix,iy,iz,is,js) = tot(ix_tot,iy_tot,iz_tot,is,js)
    end do
    end do
    end do
    end do
    end do

    call set_vxc_mat_noncollinear(mg,frg)

    deallocate(tot_tmp,tot,frg)
  end subroutine calc_vlocal_fragment_dcdft_soi

  subroutine calc_total_energy_dcdft_soi(mg,system,info,v_local,spsi,shpsi,sttpsi,ewald,pp,ppg,rion_update,dc,energy)
    use structures
    use communication, only: comm_summation
    use Total_Energy, only: calc_Total_Energy_periodic
    use pseudo_pt_so_sub, only: pseudo_so
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none
    type(s_rgrid),        intent(in) :: mg
    type(s_dft_system),   intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_scalar)       ,intent(in) :: v_local(system%nspin)
    type(s_orbital),      intent(in) :: spsi,shpsi,sttpsi
    type(s_ewald_ion_ion),intent(in) :: ewald
    type(s_pp_info)      ,intent(in) :: pp
    type(s_pp_grid)      ,intent(in) :: ppg
    logical              ,intent(in) :: rion_update
    type(s_dcdft),        intent(in) :: dc
    type(s_dft_energy)               :: energy
    type(s_orbital) :: nloc_psi
    integer :: ispin,io
    integer,dimension(3) :: is,ie
    real(8) :: E_tmp,E_local(2),E_sum(2)
    complex(8) :: ztmp

    is(1:3) = mg%is(1:3)
    ie(1:3) = min(mg%ie(1:3),dc%nxyz_domain(1:3))

    E_tmp = 0d0
    if(allocated(spsi%zwf)) then
      if(dc%id_tot==0) then
        do ispin=1,system%Nspin
          if(any(spsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,info%io_s:info%io_e,1,1) /= &
            & spsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,info%io_s:info%io_e,1,1))) then
            write(*,'(1x,a,i0)') '[DC-SOI][WARN] NaN in spsi%zwf before E_ion_nloc at spin=',ispin
          end if
          if(any(shpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,info%io_s:info%io_e,1,1) /= &
            & shpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,info%io_s:info%io_e,1,1))) then
            write(*,'(1x,a,i0)') '[DC-SOI][WARN] NaN in shpsi%zwf before E_ion_nloc at spin=',ispin
          end if
          if(any(sttpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,info%io_s:info%io_e,1,1) /= &
            & sttpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,info%io_s:info%io_e,1,1))) then
            write(*,'(1x,a,i0)') '[DC-SOI][WARN] NaN in sttpsi%zwf before E_ion_nloc at spin=',ispin
          end if
          if(any(v_local(ispin)%f(is(1):ie(1),is(2):ie(2),is(3):ie(3)) /= &
            & v_local(ispin)%f(is(1):ie(1),is(2):ie(2),is(3):ie(3)))) then
            write(*,'(1x,a,i0)') '[DC-SOI][WARN] NaN in v_local before E_ion_nloc at spin=',ispin
          end if
        end do
      end if
      do ispin=1,system%Nspin
      do io=info%io_s,info%io_e
        ztmp = sum( conjg(spsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1)) &
            &     * sttpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) )
        E_tmp = E_tmp + system%rocc(io,1,ispin) * dble(ztmp) * system%Hvol
      end do
      end do
    else if(allocated(spsi%rwf)) then
      do ispin=1,system%Nspin
      do io=info%io_s,info%io_e
        E_tmp = E_tmp + system%rocc(io,1,ispin) &
            & * sum(spsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) &
            &     * sttpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1)) * system%Hvol
      end do
      end do
    else
      stop "calc_total_energy_dcdft_soi: neither rwf nor zwf is allocated (E_kin)."
    end if
    E_local(1) = E_tmp

    E_tmp = 0d0
    if(allocated(spsi%zwf)) then
      call allocate_orbital_complex(system%nspin,mg,info,nloc_psi)
      ! shpsi also contains the noncollinear XC operator; evaluate SO nonlocal
      ! energy from the spin-orbit pseudopotential alone, as in calc_eigen_energy.
      call pseudo_so(spsi,nloc_psi,info,system%Nspin,ppg,mg)
      do ispin=1,system%Nspin
      do io=info%io_s,info%io_e
        ztmp = sum( conjg(spsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1)) &
            &     * nloc_psi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) )
        E_tmp = E_tmp + system%rocc(io,1,ispin) * dble(ztmp) * system%hvol
      end do
      end do
      call deallocate_orbital(nloc_psi)
    else if(allocated(spsi%rwf)) then
      stop "calc_total_energy_dcdft_soi: real wavefunction is invalid for SOI E_ion_nloc."
    else
      stop "calc_total_energy_dcdft_soi: neither rwf nor zwf is allocated (E_ion_nloc)."
    end if
    E_local(2) = E_tmp

    call comm_summation(E_local,E_sum,2,info%icomm_rko)
    E_local = 0d0
    if(info%id_rko == 0) E_local = E_sum
    call comm_summation(E_local,E_sum,2,dc%icomm_tot)

    energy%E_kin = E_sum(1)
    energy%E_ion_nloc = E_sum(2)
    if(dc%id_tot==0) then
      if(ieee_is_nan(energy%E_kin) .or. ieee_is_nan(energy%E_ion_nloc)) then
        write(*,'(1x,a,2(1x,es23.15))') '[DC-SOI][WARN] pre-periodic NaN: E_kin,E_ion_nloc=', &
          energy%E_kin, energy%E_ion_nloc
      end if
    end if

    call calc_Total_Energy_periodic(dc%mg_tot,ewald,dc%system_tot,dc%info_tot,pp,dc%ppg_tot, &
      & dc%fg_tot,dc%poisson_tot,rion_update,energy)
    if(dc%id_tot==0) then
      if(ieee_is_nan(energy%E_tot)) then
        write(*,'(1x,a,5(1x,es23.15))') '[DC-SOI][WARN] post-periodic NaN: Etot,Ekin,Enloc,Eh,Exc=', &
          energy%E_tot, energy%E_kin, energy%E_ion_nloc, energy%E_h, energy%E_xc
      end if
    end if
  end subroutine calc_total_energy_dcdft_soi

  subroutine write_total_dcdft_soi(system,dc)
    use structures
    use dcdft, only: write_total_dcdft
    implicit none
    type(s_dcdft) :: dc
    type(s_dft_system),intent(in) :: system

    call write_total_dcdft(system,dc)
  end subroutine write_total_dcdft_soi

end module dcdft_soi
