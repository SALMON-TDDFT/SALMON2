!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!
! Scaffold driver for theory='gw'.
! Loads the ground-state restart, then emits quasiparticle (QP) energies
! as a plain text file.  In this scaffold the QP energies equal the KS
! eigenvalues (passthrough); no self-energy physics is implemented yet.
!
! Restart-load sequence (mirrors main_dft / initialization_rt paths):
!
!   1. init_dft  -- sets up real-space grids, k-mesh, parallel distribution,
!                   and Poisson solver (same call as main_dft line 107).
!   2. initialization1_dft -- allocates energy%esp, wavefunctions, scalar
!                   fields; reads the pseudopotential file.
!   3. initialization2_dft -- calls restart_gs (when yn_restart='y') which
!                   reads occupation (system%rocc) and wavefunctions from
!                   the binary restart directory.  Then builds V_local from
!                   the restarted density and calls calc_eigen_energy once
!                   to populate energy%esp from the loaded wavefunctions.
!
! Simplification / concern:
!   When yn_restart='n' no prior GS exists; initialization2_dft will
!   call init_wf (random wavefunction) and the resulting energy%esp will
!   be meaningless.  A real GW run must always supply a converged restart.
!   This is not guarded against here; the controller should enforce it.
!
! No proper nouns appear in this file per the project constraint.
!
subroutine main_gw
  use structures
  use inputoutput
  use parallelization,    only: nproc_id_global, nproc_group_global
  use communication,      only: comm_is_root
  use salmon_xc
  use initialization_sub
  use initialization_dft
  use gw_qp_output_sub,   only: write_qp_energies
  use sendrecv_grid
  implicit none

  ! ---- local derived-type objects (same set as main_dft) ----
  type(s_rgrid)          :: lg, mg
  type(s_parallel_info)  :: info
  type(s_sendrecv_grid)  :: srg, srg_scalar
  type(s_orbital)        :: spsi, shpsi, sttpsi
  type(s_dft_system)     :: system
  type(s_poisson)        :: poisson
  type(s_stencil)        :: stencil
  type(s_xc_functional)  :: xc_func
  type(s_scalar)         :: rho, rho_jm, Vh, Vpsl
  type(s_scalar),allocatable :: V_local(:), rho_s(:), Vxc(:)
  type(s_reciprocal_grid):: fg
  type(s_pp_info)        :: pp
  type(s_pp_grid)        :: ppg
  type(s_pp_nlcc)        :: ppn
  type(s_dft_energy)     :: energy
  type(s_ewald_ion_ion)  :: ewald
  type(s_mixing)         :: mixing
  type(s_ofile)          :: ofl

  ! QP arrays (band, k, spin) — passthrough scaffold values
  real(8),allocatable :: eqp (:,:,:)
  real(8),allocatable :: zfac(:,:,:)
  real(8),allocatable :: sigx(:,:,:)
  real(8),allocatable :: sigc(:,:,:)
  real(8),allocatable :: vxc_arr(:,:,:)

  integer :: Miter, nspin
  logical :: rion_update
  integer :: ib_min, ib_max

  ! ----------------------------------------------------------------
  ! Step 0: report active GW input parameters (root only)
  ! ----------------------------------------------------------------
  if (comm_is_root(nproc_id_global)) then
    write(*,*)
    write(*,*) "  [gw] theory=gw  (scaffold: QP = KS passthrough)"
    write(*,*) "  [gw] epsilon_cutoff =", epsilon_cutoff, " Ry"
    write(*,*) "  [gw] n_empty_gw     =", n_empty_gw
    write(*,*) "  [gw] sigma_type     = ", trim(sigma_type)
    write(*,*) "  [gw] nband_qp_min   =", nband_qp_min
    write(*,*) "  [gw] nband_qp_max   =", nband_qp_max
  end if

  ! ----------------------------------------------------------------
  ! Step 1: initialise XC functional (identical to main_dft line 100)
  ! ----------------------------------------------------------------
  call init_xc(xc_func, spin, cval, xcname=xc, xname=xname, cname=cname)

  ! ----------------------------------------------------------------
  ! Step 2: init_dft — grids, k-mesh, parallel info, Poisson setup
  !         (identical to main_dft line 107)
  ! ----------------------------------------------------------------
  call init_dft(nproc_group_global, info, lg, mg, system, stencil, fg, &
                poisson, srg, srg_scalar, ofl)

  allocate(rho_s  (system%nspin))
  allocate(V_local(system%nspin))
  allocate(Vxc    (system%nspin))

  ! ----------------------------------------------------------------
  ! Step 3: initialization1_dft — allocates energy%esp, wavefunctions,
  !         reads pseudopotential, builds ppg
  !         (identical to main_dft lines 110-117)
  ! ----------------------------------------------------------------
  call initialization1_dft(system, energy, stencil, fg, poisson, &
                            lg, mg, info, srg, srg_scalar,        &
                            rho, rho_jm, rho_s, Vh, V_local, Vpsl, Vxc, &
                            spsi, shpsi, sttpsi, pp, ppg, ppn, ofl)

  ! ----------------------------------------------------------------
  ! Step 4: initialization2_dft — restart_gs (loads system%rocc and
  !         wavefunctions from directory_read_data), builds V_local,
  !         then calls calc_eigen_energy to populate energy%esp.
  !         (identical to main_dft lines 119-126)
  ! ----------------------------------------------------------------
  nspin = system%nspin
  rion_update = .false.
  call initialization2_dft(Miter, nspin, rion_update,            &
                            system, energy, ewald, stencil, fg, poisson, &
                            lg, mg, info, srg, srg_scalar,               &
                            rho, rho_jm, rho_s, Vh, V_local, Vpsl, Vxc, &
                            spsi, shpsi, sttpsi, pp, ppg, ppn,           &
                            xc_func, mixing)

  ! After the above call:
  !   energy%esp   (no, nk, nspin) -- KS eigenvalues in a.u.
  !   system%rocc  (no, nk, nspin) -- occupations
  !   system%nk, system%no, system%nspin -- dimensions

  ! ----------------------------------------------------------------
  ! Step 5: QP passthrough — allocate and fill
  ! ----------------------------------------------------------------
  allocate(eqp    (system%no, system%nk, system%nspin))
  allocate(zfac   (system%no, system%nk, system%nspin))
  allocate(sigx   (system%no, system%nk, system%nspin))
  allocate(sigc   (system%no, system%nk, system%nspin))
  allocate(vxc_arr(system%no, system%nk, system%nspin))

  eqp     = energy%esp  ! QP = KS passthrough
  zfac    = 1.0d0
  sigx    = 0.0d0
  sigc    = 0.0d0
  vxc_arr = 0.0d0

  ! Determine band output window; nband_qp_max=0 means all bands
  ib_min = max(1, nband_qp_min)
  if (nband_qp_max < 1) then
    ib_max = system%no
  else
    ib_max = min(nband_qp_max, system%no)
  end if

  ! ----------------------------------------------------------------
  ! Step 6: write output
  ! ----------------------------------------------------------------
  call write_qp_energies(system, energy%esp, eqp, zfac, sigx, sigc, vxc_arr)

  if (comm_is_root(nproc_id_global)) then
    write(*,*) "  [gw] wrote QP energies (scaffold passthrough)"
    write(*,*) "  [gw] band window: ", ib_min, " to ", ib_max
  end if

  deallocate(eqp, zfac, sigx, sigc, vxc_arr)

end subroutine main_gw
