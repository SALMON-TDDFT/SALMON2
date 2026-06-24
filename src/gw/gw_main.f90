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
  use gw_qp_output_sub,   only: write_qp_energies, write_qp_dos, read_qp_energies
  use gw_coulomb_sub,     only: build_gvectors, build_vcoul
  use gw_mtxel_sub,       only: calc_mtxel, replicate_orbitals_k
  use gw_epsilon_sub,     only: calc_epsinv, calc_chi0_freq, calc_w_freq
  use gw_absorption_sub,  only: calc_absorption, calc_absorption_velocity
  use gw_sigma_x_sub,     only: calc_sigma_x
  use gw_sigma_cohsex_sub,only: calc_sigma_cohsex
  use gw_sigma_gpp_sub,   only: calc_sigma_gpp, calc_sigma_gpp_qcache, calc_sigma_gpp_sym
  use gw_band_caps_sub,   only: resolve_band_caps
  use gw_sigma_c_real_sub,only: calc_sigma_c_real, calc_sigma_c_real_qcache, calc_sigma_c_real_sym
  use gw_qp_sub,          only: calc_vxc_expect, solve_qp
  use gw_symmetry_sub,    only: gw_sym_selftest, gw_grid_perm_selftest, gw_symmetrize_orbitals
  use sendrecv_grid
  implicit none

  ! ---- local derived-type objects (same set as main_dft) ----
  type(s_rgrid)          :: lg, mg
  type(s_parallel_info)  :: info
  type(s_sendrecv_grid)  :: srg, srg_scalar
  type(s_orbital)        :: spsi, shpsi, sttpsi
  type(s_orbital)        :: spsi_full   ! k-replicated orbitals for node-parallel sigma
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
  ! t2/t3/t4 sanity blocks are validated; compile them out of production runs
  ! (set .true. to re-enable the per-stage diagnostics).
  logical, parameter :: run_sanity_t234 = .false.
  ! [gw][sym] self-test: validate the eps^{-1} G-rotation against a direct solve
  ! (set .true. + provide sym.dat + run on one rank).
  logical, parameter :: run_sanity_sym  = .false.
  ! [gw][chi0w] self-test (spec-b1): omega=0,eta=0 limit of the real-frequency
  ! chi0 reproduces calc_epsinv's static weight (set .true., run -n1).
  logical, parameter :: run_sanity_chi0w = .false.
  ! --- variables for the [gw][chi0w] self-test ---
  integer               :: ng_cw, iq_cw, ismall_cw
  real(8),  allocatable :: gvec_cw(:,:), gg_cw(:)
  real(8)               :: qvec_cw(3), qabs_cw, qabs_min_cw, eta_cw, omg_cw(1)
  complex(8),allocatable:: chi0_cw(:,:,:)
  logical               :: ok_cw
  ! --- W-level gate (b1-2): epsinv_w(w=0) vs calc_epsinv static ---
  complex(8),allocatable:: epsinv_w_cw(:,:,:), epsinv_ref_cw(:,:)
  real(8),  allocatable :: vcoul_cw(:)
  real(8)               :: wmax_cw
  integer               :: igw_cw, jgw_cw
  ! --- absorption mode (b1-3): full-frequency eps_M(w) at q->0 proxy ---
  integer               :: ng_ab, ig0_ab, ismall_ab, iq_ab, iw_ab
  integer               :: n_occ_ab, io_ab, is_ab          ! gw_response scissors
  real(8),  allocatable :: esp_ab(:,:,:)                   ! KS or QP-scissored energies
  real(8),  allocatable :: eqp_read_ab(:,:,:)              ! per-state QP (b1-3b inject)
  logical               :: ok_qp_ab
  real(8),  allocatable :: gvec_ab(:,:), gg_ab(:), omega_grid_ab(:), vcoul_ab(:)
  real(8)               :: qvec_ab(3), qabs_ab, qabsmin_ab, eta_ab
  complex(8),allocatable:: chi0_ab(:,:,:), epsinv_ab(:,:,:), eps_macro_ab(:)
  logical               :: ok_ab
  integer               :: fh_ab

  ! --- variables for the task-2 sanity block ---
  integer,  parameter   :: ngmax_t2 = 100000
  integer               :: ng_t2
  real(8),  allocatable :: gvec_t2(:,:), gg_t2(:), vcoul_t2(:)
  real(8)               :: qzero(3)
  integer               :: ig_t2, nprint_t2, nsub1, nsub2
  real(8)               :: fourpi_t2, analytic_t2

  ! --- variables for the task-3 sanity block (plane-wave matrix elements) ---
  integer               :: ng_t3, ik_t3, nb_t3, ig0_t3, ib, jb
  real(8),  allocatable :: gvec_t3(:,:), gg_t3(:)
  complex(8),allocatable:: mtxel_t3(:,:,:)
  real(8)               :: errmax_t3, offmax_t3

  ! --- variables for the task-4 sanity block (static RPA dielectric) ---
  integer               :: ng_t4, iq_t4, nq_t4, ismallest
  real(8),  allocatable :: gvec_t4(:,:), gg_t4(:)
  complex(8),allocatable:: epsinv_t4(:,:)
  real(8),  allocatable :: epsd_t4(:)
  real(8)               :: qvec_t4(3), resid_t4, qabs, qabs_min
  real(8)               :: eps_M, eps_inf
  logical               :: ok_t4
  integer               :: nskip_t4

  ! --- variables for task 5 (bare exchange Sigma_x + QP solve) ---
  integer               :: ng_t5
  real(8),  allocatable :: gvec_t5(:,:), gg_t5(:)
  real(8),  allocatable :: sigx_w(:,:), vxc_w(:,:), sigc_w(:,:), eqp_w(:,:)
  integer               :: is_t5, ik_t5, ivtop, icbot, io
  real(8)               :: gap_ks, gap_qp
  real(8),  parameter   :: hartree2ev = 27.21138505d0  ! eV per Hartree (display only)
  real(8),  allocatable :: sigx_w0(:,:)

  ! --- variables for task 6 (static COHSEX self-energy + QP solve) ---
  integer               :: ng_t6
  real(8),  allocatable :: gvec_t6(:,:), gg_t6(:)
  real(8),  allocatable :: sigc_w6(:,:), sex_w6(:,:), coh_w6(:,:)
  real(8),  allocatable :: vxc_w6(:,:), eqp_w6(:,:), sigx_w6(:,:)
  integer               :: is_t6, ik_t6, ivtop6, icbot6
  real(8)               :: gap_ks6, gap_cohsex6, skipfrac6

  ! --- variables for task 7 (dynamic G0W0 via the GPP self-energy) ---
  integer               :: ng_t7
  real(8),  allocatable :: gvec_t7(:,:), gg_t7(:)
  real(8),  allocatable :: sigc_w7(:,:), zfac_w7(:,:)
  real(8),  allocatable :: vxc_w7(:,:), eqp_w7(:,:), sigx_w7(:,:)
  integer               :: is_t7, ik_t7, ivtop7, icbot7
  integer               :: nocc_t7, io_t7, neps_t7, nsig_t7   ! eps/sigma band caps
  ! b1-4 real-axis Sigma_c (sigma_type=='gpp_real'): own omega' grid
  real(8),  allocatable :: omega_grid_t7(:)
  real(8)               :: eta_t7
  real(8)               :: ebar_max_t7, ev_min_t7, om_need_t7   ! extrapolar grid-cover
  integer               :: iw_t7
  ! sp3 spectral function A(k,w) / Im Sigma (yn_out_gw_spectral)
  real(8),    allocatable :: w_scan_t7(:)
  complex(8), allocatable :: sigc_scan_t7(:,:)
  integer                 :: nw_scan_t7, iws_t7, io_spec, fh_spec
  real(8)                 :: resc7, imsc7, ek7, sx7, vx7, aval7
  ! sp2 optical constants n, k, R, alpha derived from eps_M(w)
  real(8)                 :: ree_ab, ime_ab, mag_ab, nopt_ab, kopt_ab, refl_ab, alph_ab
  real(8)               :: gap_ks7, gap_g0w0_7, skipfrac7

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
  ! Guard: GW reads a converged GS through the restart path
  ! (yn_restart='y' + directory_read_data). Otherwise initialization2_dft
  ! falls back to random orbitals (init_wf) and the eigenvalues are
  ! meaningless. Enforce the precondition with a clear message.
  ! ----------------------------------------------------------------
  if (yn_restart /= 'y') then
    if (comm_is_root(nproc_id_global)) then
      write(*,*) "  [gw] ERROR: theory=gw requires yn_restart='y' with"
      write(*,*) "  [gw]        directory_read_data pointing at a converged GS."
    end if
    stop
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

  if (run_sanity_t234) then
  ! ----------------------------------------------------------------
  ! Task-2 sanity block: verify build_gvectors and build_vcoul.
  ! All output prefixed with [gw][t2] for grep.
  ! Runs on the root process only; allocations are root-local.
  ! ----------------------------------------------------------------
  if (comm_is_root(nproc_id_global)) then
    allocate(gvec_t2(3, ngmax_t2))
    allocate(gg_t2(ngmax_t2))
    allocate(vcoul_t2(ngmax_t2))

    call build_gvectors(system%primitive_b, epsilon_cutoff, ngmax_t2, &
                        ng_t2, gvec_t2, gg_t2)
    write(*,*)
    write(*,*) "[gw][t2] G-vector count  ng =", ng_t2
    write(*,*) "[gw][t2] ecut (Ry)          =", epsilon_cutoff
    nprint_t2 = min(5, ng_t2)
    write(*,*) "[gw][t2] smallest |G|^2 values (bohr^-2):"
    do ig_t2 = 1, nprint_t2
      write(*,'(A,I4,A,ES14.6)') "  [gw][t2]   ig=", ig_t2, "  |G|^2=", gg_t2(ig_t2)
    end do

    ! v(q=0, G) should equal 4*pi/|G|^2 for G /= 0.
    fourpi_t2 = 4.0d0 * acos(-1.0d0)
    qzero     = 0.0d0
    nsub1     = 10
    nsub2     = 30
    call build_vcoul(ng_t2, gvec_t2, gg_t2, qzero, system%det_a, system%nk, nsub1, vcoul_t2)
    write(*,*)
    write(*,*) "[gw][t2] v(q=0,G) vs 4pi/|G|^2 for 3 non-zero G (nsub=", nsub1, "):"
    nprint_t2 = min(4, ng_t2)  ! print ig=2..4 (skip ig=1 which is G=0 head)
    do ig_t2 = 2, nprint_t2
      if (gg_t2(ig_t2) > 0.0d0) then
        analytic_t2 = fourpi_t2 / gg_t2(ig_t2)
        write(*,'(A,I4,A,ES14.6,A,ES14.6)') &
          "  [gw][t2]   ig=", ig_t2, &
          "  vcoul=", vcoul_t2(ig_t2), "  4pi/|G|^2=", analytic_t2
      end if
    end do

    ! q->0 head at two sampling densities to show convergence
    write(*,*)
    write(*,*) "[gw][t2] q->0 head (mini-BZ average) convergence:"
    write(*,'(A,I3,A,ES14.6)') "  [gw][t2]   nsub=", nsub1, &
      "  head=", vcoul_t2(1)
    call build_vcoul(ng_t2, gvec_t2, gg_t2, qzero, system%det_a, system%nk, nsub2, vcoul_t2)
    write(*,'(A,I3,A,ES14.6)') "  [gw][t2]   nsub=", nsub2, &
      "  head=", vcoul_t2(1)
    write(*,*)

    deallocate(gvec_t2, gg_t2, vcoul_t2)
  end if

  ! ----------------------------------------------------------------
  ! Task-3 sanity block: plane-wave matrix elements M_{nn'}(k,q=0,G).
  ! All output prefixed with [gw][t3].  calc_mtxel performs collective
  ! MPI (comm_summation over icomm_rko), so the build/compute runs on ALL
  ! ranks; only the printing is root-only.
  !
  !   M_{nn}(k,0,G=0)  must equal 1   (pins the normalisation)
  !   M_{nn'}(k,0,G=0) (n/=n') must be ~0 (orthonormality)
  !
  ! The G-vector list is rebuilt here on every rank (cheap, deterministic).
  ! q=0 means ikq = ik (no umklapp).
  ! ----------------------------------------------------------------
  allocate(gvec_t3(3, ngmax_t2))
  allocate(gg_t3(ngmax_t2))
  call build_gvectors(system%primitive_b, epsilon_cutoff, ngmax_t2, &
                      ng_t3, gvec_t3, gg_t3)

  ! index of the G=0 entry in the gvec list (build_gvectors sorts |G|^2
  ! ascending, so ig=1 is G=0).
  ig0_t3 = 1

  ! Use the first k-point and the first few bands; spin 1 (non-magnetic test).
  ik_t3 = 1
  nb_t3 = min(5, system%no)

  allocate(mtxel_t3(ng_t3, nb_t3, nb_t3))
  call calc_mtxel(system, info, mg, lg, spsi, &
                  gvec_t3, ng_t3, ik_t3, ik_t3, 1, nb_t3, nb_t3, mtxel_t3)

  if (comm_is_root(nproc_id_global)) then
    write(*,*)
    write(*,*) "[gw][t3] plane-wave matrix elements M_{nn'}(k,q=0,G)"
    write(*,*) "[gw][t3] ng =", ng_t3, "  k-index =", ik_t3, "  bands =", nb_t3
    write(*,*) "[gw][t3] diagonal M_{nn}(k,0,G=0) (target 1.0):"
    errmax_t3 = 0.0d0
    do ib = 1, nb_t3
      write(*,'(A,I3,A,ES16.8,A,ES12.4)') "  [gw][t3]   n=", ib, &
        "  Re=", dble(mtxel_t3(ig0_t3,ib,ib)), &
        "  Im=", aimag(mtxel_t3(ig0_t3,ib,ib))
      errmax_t3 = max(errmax_t3, abs(mtxel_t3(ig0_t3,ib,ib) - (1.0d0,0.0d0)))
    end do
    write(*,'(A,ES12.4)') "  [gw][t3]   max |M_nn - 1| =", errmax_t3

    write(*,*) "[gw][t3] off-diagonal M_{nn'}(k,0,G=0) (target 0.0):"
    offmax_t3 = 0.0d0
    do ib = 1, nb_t3
    do jb = 1, nb_t3
      if (ib == jb) cycle
      offmax_t3 = max(offmax_t3, abs(mtxel_t3(ig0_t3,ib,jb)))
    end do
    end do
    ! print a couple of representative off-diagonal entries
    if (nb_t3 >= 2) then
      write(*,'(A,ES16.8,A,ES16.8)') "  [gw][t3]   M_{12}(G=0): Re=", &
        dble(mtxel_t3(ig0_t3,1,2)), "  Im=", aimag(mtxel_t3(ig0_t3,1,2))
    end if
    if (nb_t3 >= 3) then
      write(*,'(A,ES16.8,A,ES16.8)') "  [gw][t3]   M_{13}(G=0): Re=", &
        dble(mtxel_t3(ig0_t3,1,3)), "  Im=", aimag(mtxel_t3(ig0_t3,1,3))
    end if
    write(*,'(A,ES12.4)') "  [gw][t3]   max |M_nn'(n/=n')| =", offmax_t3
    write(*,*)
  end if

  deallocate(gvec_t3, gg_t3, mtxel_t3)

  ! ----------------------------------------------------------------
  ! Task-4 sanity block: static (omega=0) RPA dielectric matrix.
  ! All output prefixed with [gw][t4].  calc_epsinv calls calc_mtxel
  ! (collective), so it is entered on ALL ranks; only printing is root-only.
  !
  ! q-set: q_iq = vec_k(:,iq) - vec_k(:,1), iq = 1..nk (the mesh-point
  ! differences relative to the first k-point).  iq=1 is q=0 (head).  For each
  ! q we report the macroscopic dielectric eps_M(q) = 1/epsinv(1,1) (ig=1 is
  ! G=0).  For the smallest non-zero |q| we also report the inversion residual
  ! and use eps_M there as a (coarse) epsilon_infinity estimate.
  ! ----------------------------------------------------------------
  allocate(gvec_t4(3, ngmax_t2))
  allocate(gg_t4(ngmax_t2))
  call build_gvectors(system%primitive_b, epsilon_cutoff, ngmax_t2, &
                      ng_t4, gvec_t4, gg_t4)
  allocate(epsinv_t4(ng_t4, ng_t4))
  allocate(epsd_t4(ng_t4))

  nq_t4 = system%nk
  ! locate the smallest non-zero |q| in the q-set (relative to vec_k(:,1))
  qabs_min  = huge(1.0d0)
  ismallest = 0
  do iq_t4 = 1, nq_t4
    qvec_t4(1:3) = system%vec_k(1:3,iq_t4) - system%vec_k(1:3,1)
    qabs = sqrt(qvec_t4(1)**2 + qvec_t4(2)**2 + qvec_t4(3)**2)
    if (qabs > 1.0d-8 .and. qabs < qabs_min) then
      qabs_min  = qabs
      ismallest = iq_t4
    end if
  end do

  if (comm_is_root(nproc_id_global)) then
    write(*,*)
    write(*,*) "[gw][t4] static RPA dielectric matrix (omega=0)"
    write(*,*) "[gw][t4] ng =", ng_t4, "  nq =", nq_t4, &
               "  smallest |q| at iq=", ismallest
    write(*,*) "[gw][t4] eps_M(q) = 1/epsinv(G=0,G=0) per q:"
  end if

  nskip_t4 = 0
  eps_inf  = 0.0d0
  do iq_t4 = 1, nq_t4
    qvec_t4(1:3) = system%vec_k(1:3,iq_t4) - system%vec_k(1:3,1)
    qabs = sqrt(qvec_t4(1)**2 + qvec_t4(2)**2 + qvec_t4(3)**2)

    call calc_epsinv(system, info, mg, lg, spsi, energy%esp, gvec_t4, gg_t4, ng_t4, &
                     iq_t4, qvec_t4, 1, epsinv_t4, eps_diag=epsd_t4, &
                     residual=resid_t4, ok=ok_t4)

    if (comm_is_root(nproc_id_global)) then
      if (.not. ok_t4) then
        nskip_t4 = nskip_t4 + 1
        write(*,'(A,I4,A,ES12.4,A)') "  [gw][t4]   iq=", iq_t4, "  |q|=", qabs, &
          "  SKIPPED (no k+q partner; symmetry-reduced mesh)"
      else
        eps_M = dble( 1.0d0 / epsinv_t4(1,1) )
        write(*,'(A,I4,A,ES12.4,A,F12.5)') "  [gw][t4]   iq=", iq_t4, &
          "  |q|=", qabs, "  eps_M=", eps_M
        if (iq_t4 == ismallest) then
          eps_inf = eps_M
          write(*,'(A,ES12.4)') "  [gw][t4]   inversion residual max|eps.epsinv-I| =", resid_t4
        end if
      end if
    end if
  end do

  if (comm_is_root(nproc_id_global)) then
    write(*,'(A,F12.5)') "  [gw][t4] epsilon_infinity estimate (eps_M at smallest |q|) =", eps_inf
    write(*,*) "[gw][t4] (coarse on a small k-mesh; q->0 head and BZ sampling are limiting factors)"
    if (nskip_t4 > 0) then
      write(*,*) "[gw][t4] NOTE: ", nskip_t4, " q-point(s) skipped (no mesh partner)."
    end if
    write(*,*)
  end if

  deallocate(gvec_t4, gg_t4, epsinv_t4, epsd_t4)
  end if  ! run_sanity_t234

  ! ----------------------------------------------------------------
  ! [gw][sym] self-test: validate the eps^{-1} G-rotation convention used by the
  ! symmetry-reduced BZ sum, in isolation, before it feeds the self-energy.
  ! ----------------------------------------------------------------
  if (run_sanity_sym) then
    call gw_grid_perm_selftest(system, lg)
    ! make the full-mesh orbitals exactly symmetric (orbital rotation), then the
    ! eps G-rotation discrepancy should collapse to machine precision.
    call gw_symmetrize_orbitals(system, info, lg, spsi, energy)
    call gw_sym_selftest(system, info, mg, lg, spsi, energy%esp, epsilon_cutoff)
  end if

  ! ----------------------------------------------------------------
  ! [gw][chi0w] self-test (spec-b1): the omega=0, eta=0 limit of the
  ! real-frequency polarizability must reproduce calc_epsinv's static
  ! per-pair weight (regression gate iii).  One omega=0 point, zero
  ! broadening, smallest non-zero q (exercises the umklapp/imap path).
  ! ----------------------------------------------------------------
  if (run_sanity_chi0w) then
    allocate(gvec_cw(3, ngmax_t2), gg_cw(ngmax_t2))
    call build_gvectors(system%primitive_b, epsilon_cutoff, ngmax_t2, &
                        ng_cw, gvec_cw, gg_cw)
    qabs_min_cw = huge(1.0d0); ismall_cw = 0
    do iq_cw = 1, system%nk
      qvec_cw(1:3) = system%vec_k(1:3,iq_cw) - system%vec_k(1:3,1)
      qabs_cw = sqrt(sum(qvec_cw(1:3)**2))
      if (qabs_cw > 1.0d-8 .and. qabs_cw < qabs_min_cw) then
        qabs_min_cw = qabs_cw; ismall_cw = iq_cw
      end if
    end do
    if (ismall_cw == 0) ismall_cw = 1
    qvec_cw(1:3) = system%vec_k(1:3,ismall_cw) - system%vec_k(1:3,1)
    omg_cw(1) = 0.0d0
    eta_cw    = 0.0d0
    allocate(chi0_cw(ng_cw, ng_cw, 1))
    write(*,'(A,I6,A,I4)') "  [gw][chi0w] ng =", ng_cw, "  q-index =", ismall_cw
    call calc_chi0_freq(system, info, mg, lg, spsi, energy%esp, gvec_cw, ng_cw, &
                        ismall_cw, qvec_cw, 1, 1, omg_cw, eta_cw, chi0_cw, &
                        ok=ok_cw, run_sanity=.true.)
    write(*,'(A,L2)') "  [gw][chi0w] q_ok =", ok_cw
    ! W-level gate (b1-2): eps = 1 - v chi0_w(0) inverted == calc_epsinv static.
    allocate(epsinv_w_cw(ng_cw,ng_cw,1), epsinv_ref_cw(ng_cw,ng_cw), vcoul_cw(ng_cw))
    call calc_w_freq(system, gvec_cw, gg_cw, ng_cw, qvec_cw, 1, chi0_cw, &
                     epsinv_w_cw, vcoul_cw)
    call calc_epsinv(system, info, mg, lg, spsi, energy%esp, gvec_cw, gg_cw, ng_cw, &
                     ismall_cw, qvec_cw, 1, epsinv_ref_cw)
    wmax_cw = 0.0d0
    do jgw_cw = 1, ng_cw
      do igw_cw = 1, ng_cw
        wmax_cw = max(wmax_cw, abs(epsinv_w_cw(igw_cw,jgw_cw,1) - epsinv_ref_cw(igw_cw,jgw_cw)))
      end do
    end do
    write(*,'(A,ES12.4)') "  [gw][wfreq] max|epsinv_w(0)-epsinv_static| = ", wmax_cw
    deallocate(epsinv_w_cw, epsinv_ref_cw, vcoul_cw)
    deallocate(gvec_cw, gg_cw, chi0_cw)
  end if

  ! ----------------------------------------------------------------
  ! Absorption mode (spec-b1, b1-3): full-frequency macroscopic dielectric
  ! eps_M(q->0,w) = 1/epsinv_w(0,0,w) over the real-frequency grid, written to
  ! <sysname>_absorption.data for comparison with the RT-TDDFT Im eps.  This is
  ! a standalone deliverable: it returns before the QP/self-energy path.
  ! StageA head = smallest-|q| proxy for q->0.
  ! ----------------------------------------------------------------
  if (yn_gw_absorption == 'y' .or. &
      (trim(theory) == 'dft_response' .and. yn_out_gw_eps == 'y') .or. &
      (trim(theory) == 'gw_response'  .and. yn_out_gw_eps == 'y')) then
    allocate(gvec_ab(3, ngmax_t2), gg_ab(ngmax_t2))
    call build_gvectors(system%primitive_b, epsilon_cutoff, ngmax_t2, &
                        ng_ab, gvec_ab, gg_ab)
    ! G=0 index (build_gvectors sorts by |G|^2 -> normally first; locate robustly)
    ig0_ab = 1
    do iq_ab = 1, ng_ab
      if (gg_ab(iq_ab) < 1.0d-10) then
        ig0_ab = iq_ab; exit
      end if
    end do
    ! real-frequency grid (a.u.); first point omega=0 -> Re eps_M(0)=eps_inf
    allocate(omega_grid_ab(nomega_gw), eps_macro_ab(nomega_gw))
    do iw_ab = 1, nomega_gw
      omega_grid_ab(iw_ab) = (dble(iw_ab-1)/dble(max(nomega_gw-1,1))) &
                             * omega_max_gw / 27.211386d0
    end do
    eta_ab = eta_gw / 27.211386d0

    ! theory='gw_response': RPA@QP = scissors-corrected RPA.  Shift the conduction
    ! bands (io > n_occ) by gw_scissors so chi0(w) / the velocity head see the QP
    ! gap (the absorption edge blue-shifts from the KS gap to the QP gap).  KS path
    ! (dft_response / yn_gw_absorption) leaves esp_ab = the bare KS energies.
    allocate(esp_ab(system%no, system%nk, system%nspin))
    esp_ab(:,:,:) = energy%esp(:,:,:)
    if (trim(theory) == 'gw_response' .and. yn_gw_qp_inject == 'y') then
      ! b1-3b true G0W0+RPA: inject the per-state QP spectrum (read from a prior
      ! theory='gw' run's _qp_energies.data) into chi0 -- not a rigid scissors.
      allocate(eqp_read_ab(system%no, system%nk, system%nspin))
      call read_qp_energies(system, eqp_read_ab, ok_qp_ab)
      if (ok_qp_ab) then
        esp_ab(:,:,:) = eqp_read_ab(:,:,:)
        if (comm_is_root(nproc_id_global)) &
          write(*,*) "  [gw] gw_response: per-state QP injection from _qp_energies.data (true G0W0+RPA)"
      else if (comm_is_root(nproc_id_global)) then
        write(*,*) "  [gw] gw_response: yn_gw_qp_inject='y' but _qp_energies.data missing -> KS chi0"
      end if
      deallocate(eqp_read_ab)
    else if (trim(theory) == 'gw_response' .and. abs(gw_scissors) > 1.0d-12) then
      n_occ_ab = 0
      do io_ab = 1, system%no
        if (system%rocc(io_ab,1,1) > 1.0d-6) n_occ_ab = io_ab
      end do
      do is_ab = 1, system%nspin
        do iq_ab = 1, system%nk
          do io_ab = n_occ_ab+1, system%no
            esp_ab(io_ab,iq_ab,is_ab) = esp_ab(io_ab,iq_ab,is_ab) + gw_scissors/27.211386d0
          end do
        end do
      end do
      if (comm_is_root(nproc_id_global)) &
        write(*,'(A,F8.4,A,I4)') "  [gw] gw_response: scissors =", gw_scissors, &
          " eV on conduction bands io>", n_occ_ab
    end if

    if (gw_head_mode == 'velocity') then
      ! StageB: q=0 analytic head/wing/body with the momentum matrix element.
      if (comm_is_root(nproc_id_global)) &
        write(*,*) "  [gw] absorption mode: velocity head, full LFE (q=0)"
      ! replicate orbitals so every rank owns all k, then distribute the BZ
      ! k-sum over icomm_k (local_only) and reduce (gpp-8b node-parallel pattern).
      call replicate_orbitals_k(system, info, spsi, spsi_full)
      call calc_absorption_velocity(system, info, mg, lg, stencil, srg, ppg, spsi_full, &
           esp_ab, gvec_ab, gg_ab, ng_ab, ig0_ab, 1, nomega_gw, &
           omega_grid_ab, eta_ab, eps_macro_ab, local_only=.true.)
      if (comm_is_root(nproc_id_global)) then
        open(newunit=fh_ab, file=trim(sysname)//'_absorption.data', status='replace')
        write(fh_ab,'(A)') "# velocity-head full-LFE macroscopic dielectric (RPA@KS)"
        write(fh_ab,'(A)') "# optical constants from eps_M=n^2: (n+ik)=sqrt(eps), &
          &R=((n-1)^2+k^2)/((n+1)^2+k^2), alpha=2 w k / c [1/bohr]"
        write(fh_ab,'(A)') "# 1:omega[eV] 2:Re eps_M 3:Im eps_M 4:n 5:k 6:R 7:alpha[1/bohr]"
        do iw_ab = 1, nomega_gw
          ree_ab = dble (eps_macro_ab(iw_ab))
          ime_ab = aimag(eps_macro_ab(iw_ab))
          mag_ab = sqrt(ree_ab*ree_ab + ime_ab*ime_ab)
          nopt_ab = sqrt(max(0.5d0*(mag_ab + ree_ab), 0.0d0))   ! refractive index
          kopt_ab = sqrt(max(0.5d0*(mag_ab - ree_ab), 0.0d0))   ! extinction coeff
          refl_ab = ( (nopt_ab-1.0d0)**2 + kopt_ab**2 ) &
                  / ( (nopt_ab+1.0d0)**2 + kopt_ab**2 )          ! normal-incidence R
          alph_ab = 2.0d0 * omega_grid_ab(iw_ab) * kopt_ab / 137.035999d0  ! 1/bohr
          write(fh_ab,'(7ES22.12)') omega_grid_ab(iw_ab)*27.211386d0, &
            ree_ab, ime_ab, nopt_ab, kopt_ab, refl_ab, alph_ab
        end do
        close(fh_ab)
      end if
    else
      ! StageA: smallest-|q| proxy (mesh-unstable at small q; kept for reference).
      if (comm_is_root(nproc_id_global)) &
        write(*,*) "  [gw] absorption mode: full-frequency eps_M(w) (q->0 proxy)"
      qabsmin_ab = huge(1.0d0); ismall_ab = 0
      do iq_ab = 1, system%nk
        qvec_ab(1:3) = system%vec_k(1:3,iq_ab) - system%vec_k(1:3,1)
        qabs_ab = sqrt(sum(qvec_ab(1:3)**2))
        if (qabs_ab > 1.0d-8 .and. qabs_ab < qabsmin_ab) then
          qabsmin_ab = qabs_ab; ismall_ab = iq_ab
        end if
      end do
      if (ismall_ab == 0) ismall_ab = 1
      qvec_ab(1:3) = system%vec_k(1:3,ismall_ab) - system%vec_k(1:3,1)
      allocate(chi0_ab(ng_ab,ng_ab,nomega_gw), epsinv_ab(ng_ab,ng_ab,nomega_gw))
      allocate(vcoul_ab(ng_ab))
      call calc_chi0_freq(system, info, mg, lg, spsi, esp_ab, gvec_ab, ng_ab, &
                          ismall_ab, qvec_ab, 1, nomega_gw, omega_grid_ab, eta_ab, &
                          chi0_ab, ok=ok_ab)
      call calc_w_freq(system, gvec_ab, gg_ab, ng_ab, qvec_ab, nomega_gw, chi0_ab, &
                       epsinv_ab, vcoul_ab)
      call calc_absorption(ng_ab, ig0_ab, nomega_gw, omega_grid_ab, epsinv_ab, &
                           eps_macro_ab, comm_is_root(nproc_id_global), sysname)
      deallocate(chi0_ab, epsinv_ab, vcoul_ab)
    end if

    if (comm_is_root(nproc_id_global)) then
      write(*,'(A,I5,A,F7.2,A,F7.4)') "  [gw][abs] nomega=", nomega_gw, &
        "  omega_max[eV]=", omega_max_gw, "  eta[eV]=", eta_gw
      write(*,'(A,2ES14.5)') "  [gw][abs] eps_M(w=0) (Re=eps_inf, Im~0) = ", &
        dble(eps_macro_ab(1)), aimag(eps_macro_ab(1))
    end if
    deallocate(gvec_ab, gg_ab, omega_grid_ab, eps_macro_ab, esp_ab)
    return
  end if

  ! theory='dft_response' (KS-level) and 'gw_response' (RPA@QP via the gw_scissors
  ! input) are analysis modes that never enter the QP / Sigma path below.
  if (trim(theory) == 'dft_response' .or. trim(theory) == 'gw_response') then
    if (comm_is_root(nproc_id_global)) &
      write(*,*) "  [gw] response mode: no analysis requested (set yn_out_gw_eps='y')"
    return
  end if

  ! ----------------------------------------------------------------
  ! Step 5: QP energies — full (no,nk,nspin) arrays for the output file.
  ! Defaults are the KS passthrough; the sigma_type=='sigx' branch below
  ! overwrites the band window [ib_min,ib_max] with the bare-exchange QP solve.
  ! ----------------------------------------------------------------
  allocate(eqp    (system%no, system%nk, system%nspin))
  allocate(zfac   (system%no, system%nk, system%nspin))
  allocate(sigx   (system%no, system%nk, system%nspin))
  allocate(sigc   (system%no, system%nk, system%nspin))
  allocate(vxc_arr(system%no, system%nk, system%nspin))

  eqp     = energy%esp  ! default: QP = KS passthrough
  zfac    = 1.0d0
  sigx    = 0.0d0
  sigc    = 0.0d0
  vxc_arr = 0.0d0

  ! Determine band output window; nband_qp_max=0 means all bands
  if (nband_qp_max > system%no) then
    if (comm_is_root(nproc_id_global)) &
      write(*,*) "  [gw] FATAL: nband_qp_max exceeds the available bands (system%no)"
    stop 'gw: nband_qp_max exceeds system%no'
  end if
  ib_min = max(1, nband_qp_min)
  if (nband_qp_max < 1) then
    ib_max = system%no
  else
    ib_max = min(nband_qp_max, system%no)
  end if

  ! ----------------------------------------------------------------
  ! Bare-exchange (Fock) rung: sigma_type=='sigx'
  !   Sigma_x over [ib_min,ib_max], <Vxc> over the same window, sigc=0, Z=1,
  !   then eps^QP = eps^KS + (Sigma_x - Vxc).  Real columns only.
  ! ----------------------------------------------------------------
  if (trim(sigma_type) == 'sigx') then

    ! G-vector list (rebuilt on every rank; deterministic) for Sigma_x.
    allocate(gvec_t5(3, ngmax_t2))
    allocate(gg_t5(ngmax_t2))
    call build_gvectors(system%primitive_b, epsilon_cutoff, ngmax_t2, &
                        ng_t5, gvec_t5, gg_t5)

    allocate(sigx_w(ib_min:ib_max, system%nk))
    allocate(vxc_w (ib_min:ib_max, system%nk))
    allocate(sigc_w(ib_min:ib_max, system%nk))
    allocate(eqp_w (ib_min:ib_max, system%nk))
    sigc_w(:,:) = 0.0d0

    do is_t5 = 1, system%nspin
      ! Sigma_x and <Vxc> for the band window (both collective inside).
      call calc_sigma_x(system, info, mg, lg, spsi, gvec_t5, gg_t5, ng_t5, &
                        is_t5, ib_min, ib_max, sigx_w)
      call calc_vxc_expect(system, info, mg, spsi, Vxc, is_t5, &
                           ib_min, ib_max, vxc_w)

      ! Linearised QP update (Z=1, sigc=0): eqp = eks + (sigx - vxc).
      call solve_qp(energy%esp(ib_min:ib_max,:,is_t5), sigx_w, sigc_w, vxc_w, &
                    1.0d0, eqp_w)

      ! Scatter the window into the full output arrays for this spin.
      sigx   (ib_min:ib_max,:,is_t5) = sigx_w(:,:)
      vxc_arr(ib_min:ib_max,:,is_t5) = vxc_w (:,:)
      eqp    (ib_min:ib_max,:,is_t5) = eqp_w (:,:)
      ! sigc and zfac stay at their defaults (0 and 1).
    end do

    ! ----------------------------------------------------------------
    ! [gw][t5] sanity block (root only): for spin 1, k=1, report the top
    ! valence and bottom conduction states.  Checks:
    !   - <Sigma_x> negative, O(-10 eV) for valence,
    !   - <Vxc> (eV),
    !   - KS gap vs exchange-corrected (QP) gap -> expect the QP gap WIDER,
    !   - <Sigma_x> invariance to n_empty (Sigma_x sums only occupied v):
    !       recomputed over occupied bands only and compared to the window value.
    ! ----------------------------------------------------------------
    is_t5 = 1
    ik_t5 = 1

    ! locate the top occupied (ivtop) and bottom unoccupied (icbot) bands at k=1
    ivtop = 0
    do io = 1, system%no
      if (system%rocc(io,ik_t5,is_t5) > 1.0d-6) ivtop = io
    end do
    icbot = min(ivtop + 1, system%no)

    if (comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "[gw][t5] bare exchange Sigma_x + QP solve"
      write(*,*) "[gw][t5] ng =", ng_t5, "  band window:", ib_min, "..", ib_max
      write(*,'(A,I4,A,I4)') "  [gw][t5] at k=1: top valence n=", ivtop, &
        "  bottom conduction n=", icbot
    end if

    ! n_empty invariance check: Sigma_x(n,k) sums ONLY occupied bands v, so it
    ! must be identical regardless of how many empty/output bands the run carries.
    ! Demonstrate by recomputing Sigma_x over a window that STOPS at ivtop (no
    ! empty states above it) and comparing the ivtop entry to the full-window
    ! value.  Only meaningful when ivtop is itself inside the output window.
    if (ivtop >= ib_min .and. ivtop <= ib_max) then
      allocate(sigx_w0(ib_min:ivtop, system%nk))
      call calc_sigma_x(system, info, mg, lg, spsi, gvec_t5, gg_t5, ng_t5, &
                        is_t5, ib_min, ivtop, sigx_w0)
    end if

    if (comm_is_root(nproc_id_global)) then
      if (ivtop >= ib_min .and. ivtop <= ib_max) then
        write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t5] top valence: <Sigx>=", &
          sigx(ivtop,ik_t5,is_t5)*hartree2ev, " eV   <Vxc>=", &
          vxc_arr(ivtop,ik_t5,is_t5)*hartree2ev, " eV"
        write(*,'(A,L2,A,ES12.4)') "  [gw][t5] <Sigx> negative? ", &
          (sigx(ivtop,ik_t5,is_t5) < 0.0d0), &
          "   n_empty-invariance |dSigx| (eV)=", &
          abs(sigx_w0(ivtop,ik_t5) - sigx(ivtop,ik_t5,is_t5))*hartree2ev
      end if
      if (icbot >= ib_min .and. icbot <= ib_max) then
        write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t5] bot conduction: <Sigx>=", &
          sigx(icbot,ik_t5,is_t5)*hartree2ev, " eV   <Vxc>=", &
          vxc_arr(icbot,ik_t5,is_t5)*hartree2ev, " eV"
      end if
      ! KS gap vs exchange-corrected gap at k=1 (direct gap proxy)
      if (ivtop >= 1 .and. icbot <= system%no .and. icbot > ivtop) then
        gap_ks = ( energy%esp(icbot,ik_t5,is_t5) - energy%esp(ivtop,ik_t5,is_t5) ) &
                 * hartree2ev
        if (icbot >= ib_min .and. icbot <= ib_max .and. &
            ivtop >= ib_min .and. ivtop <= ib_max) then
          gap_qp = ( eqp(icbot,ik_t5,is_t5) - eqp(ivtop,ik_t5,is_t5) ) * hartree2ev
          write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t5] direct gap @k=1: KS=", &
            gap_ks, " eV   QP(Sigx)=", gap_qp, " eV  (expect QP wider)"
        else
          write(*,'(A,F12.5,A)') "  [gw][t5] direct gap @k=1: KS=", gap_ks, &
            " eV  (QP gap needs both gap states inside the band window)"
        end if
      end if
      write(*,*)
    end if

    if (allocated(sigx_w0)) deallocate(sigx_w0)
    deallocate(sigx_w, vxc_w, sigc_w, eqp_w)
    deallocate(gvec_t5, gg_t5)

    if (comm_is_root(nproc_id_global)) then
      write(*,*) "  [gw] computed bare-exchange QP energies (sigma_type=sigx)"
    end if

  ! ----------------------------------------------------------------
  ! Static COHSEX rung: sigma_type=='cohsex'
  !   Sigma_COHSEX = Sigma_SEX + Sigma_COH over [ib_min,ib_max], <Vxc> over the
  !   same window, Z=1, then eps^QP = eps^KS + (Sigma_COHSEX - Vxc).  The sigx
  !   column is left at the bare-exchange reference (Sigma_x) for comparison; it
  !   is NOT folded into the QP energy (SEX already carries the screened
  !   exchange).  sigc holds Sigma_COHSEX.
  ! ----------------------------------------------------------------
  else if (trim(sigma_type) == 'cohsex') then

    ! G-vector list (eps G-set; rebuilt on every rank, deterministic).
    allocate(gvec_t6(3, ngmax_t2))
    allocate(gg_t6(ngmax_t2))
    call build_gvectors(system%primitive_b, epsilon_cutoff, ngmax_t2, &
                        ng_t6, gvec_t6, gg_t6)

    allocate(sigc_w6(ib_min:ib_max, system%nk))
    allocate(sex_w6 (ib_min:ib_max, system%nk))
    allocate(coh_w6 (ib_min:ib_max, system%nk))
    allocate(vxc_w6 (ib_min:ib_max, system%nk))
    allocate(eqp_w6 (ib_min:ib_max, system%nk))
    allocate(sigx_w6(ib_min:ib_max, system%nk))

    skipfrac6 = 0.0d0

    do is_t6 = 1, system%nspin
      ! Sigma_COHSEX (=SEX+COH), its SEX/COH parts, and <Vxc> for the window.
      call calc_sigma_cohsex(system, info, mg, lg, spsi, energy%esp, &
                             gvec_t6, gg_t6, ng_t6, is_t6, ib_min, ib_max, &
                             sigc_w6, sex_out=sex_w6, coh_out=coh_w6, &
                             skip_frac=skipfrac6)
      call calc_vxc_expect(system, info, mg, spsi, Vxc, is_t6, &
                           ib_min, ib_max, vxc_w6)
      ! Bare-exchange reference (Sigma_x) for the sigx output column only.
      call calc_sigma_x(system, info, mg, lg, spsi, gvec_t6, gg_t6, ng_t6, &
                        is_t6, ib_min, ib_max, sigx_w6)

      ! Linearised QP update (Z=1): eqp = eks + (Sigma_COHSEX - vxc).  Pass
      ! sigc=Sigma_COHSEX and sigx-column=0 so solve_qp adds only (sigc - vxc).
      call solve_qp(energy%esp(ib_min:ib_max,:,is_t6), &
                    0.0d0*sigc_w6, sigc_w6, vxc_w6, 1.0d0, eqp_w6)

      ! Scatter into the full output arrays for this spin.
      sigx   (ib_min:ib_max,:,is_t6) = sigx_w6(:,:)   ! bare-exchange reference
      sigc   (ib_min:ib_max,:,is_t6) = sigc_w6(:,:)   ! Sigma_COHSEX
      vxc_arr(ib_min:ib_max,:,is_t6) = vxc_w6 (:,:)
      eqp    (ib_min:ib_max,:,is_t6) = eqp_w6 (:,:)
      ! zfac stays at its default (1).
    end do

    ! ----------------------------------------------------------------
    ! [gw][t6] sanity block (root only): for spin 1, k=1, report the top
    ! valence and bottom conduction states.  Reports Sigma_SEX, Sigma_COH (eV),
    ! the COHSEX gap, and the skipped-(G,G') fraction.
    ! ----------------------------------------------------------------
    is_t6 = 1
    ik_t6 = 1
    ivtop6 = 0
    do io = 1, system%no
      if (system%rocc(io,ik_t6,is_t6) > 1.0d-6) ivtop6 = io
    end do
    icbot6 = min(ivtop6 + 1, system%no)

    if (comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "[gw][t6] static COHSEX self-energy + QP solve"
      write(*,*) "[gw][t6] ng =", ng_t6, "  band window:", ib_min, "..", ib_max
      write(*,'(A,I4,A,I4)') "  [gw][t6] at k=1: top valence n=", ivtop6, &
        "  bottom conduction n=", icbot6
      write(*,'(A,ES12.4)') "  [gw][t6] skipped (G,G') fraction (COH) =", skipfrac6

      if (ivtop6 >= ib_min .and. ivtop6 <= ib_max) then
        write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t6] top valence: <Sig_SEX>=", &
          sex_w6(ivtop6,ik_t6)*hartree2ev, " eV   <Sig_COH>=", &
          coh_w6(ivtop6,ik_t6)*hartree2ev, " eV"
        write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t6]   <Sig_COHSEX>=", &
          sigc(ivtop6,ik_t6,is_t6)*hartree2ev, " eV   <Vxc>=", &
          vxc_arr(ivtop6,ik_t6,is_t6)*hartree2ev, " eV"
      end if
      if (icbot6 >= ib_min .and. icbot6 <= ib_max) then
        write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t6] bot conduction: <Sig_SEX>=", &
          sex_w6(icbot6,ik_t6)*hartree2ev, " eV   <Sig_COH>=", &
          coh_w6(icbot6,ik_t6)*hartree2ev, " eV"
        write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t6]   <Sig_COHSEX>=", &
          sigc(icbot6,ik_t6,is_t6)*hartree2ev, " eV   <Vxc>=", &
          vxc_arr(icbot6,ik_t6,is_t6)*hartree2ev, " eV"
      end if
      ! KS gap vs COHSEX gap at k=1 (direct gap proxy)
      if (ivtop6 >= 1 .and. icbot6 <= system%no .and. icbot6 > ivtop6) then
        gap_ks6 = ( energy%esp(icbot6,ik_t6,is_t6) - energy%esp(ivtop6,ik_t6,is_t6) ) &
                  * hartree2ev
        if (icbot6 >= ib_min .and. icbot6 <= ib_max .and. &
            ivtop6 >= ib_min .and. ivtop6 <= ib_max) then
          gap_cohsex6 = ( eqp(icbot6,ik_t6,is_t6) - eqp(ivtop6,ik_t6,is_t6) ) * hartree2ev
          write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t6] direct gap @k=1: KS=", &
            gap_ks6, " eV   COHSEX=", gap_cohsex6, " eV  (expect COHSEX wider)"
        else
          write(*,'(A,F12.5,A)') "  [gw][t6] direct gap @k=1: KS=", gap_ks6, &
            " eV  (COHSEX gap needs both gap states inside the band window)"
        end if
      end if
      write(*,*)
    end if

    deallocate(sigc_w6, sex_w6, coh_w6, vxc_w6, eqp_w6, sigx_w6)
    deallocate(gvec_t6, gg_t6)

    if (comm_is_root(nproc_id_global)) then
      write(*,*) "  [gw] computed static COHSEX QP energies (sigma_type=cohsex)"
    end if

  ! ----------------------------------------------------------------
  ! Dynamic G0W0 (generalized plasmon-pole) rung: sigma_type=='gpp'
  !   Sigma_x from calc_sigma_x (bare exchange), Re Sigma_c and Z from
  !   calc_sigma_gpp, <Vxc> over the window, then the linearised one-shot
  !   QP solve  eps^QP = eps^KS + Z (Sigma_x + Re Sigma_c - Vxc).  This is
  !   the PoC goal: the gap must fall BELOW the static COHSEX value.
  ! ----------------------------------------------------------------
  else if (trim(sigma_type) == 'gpp' .or. trim(sigma_type) == 'gpp_real') then

    ! G-vector list (eps G-set; rebuilt on every rank, deterministic).
    allocate(gvec_t7(3, ngmax_t2))
    allocate(gg_t7(ngmax_t2))
    call build_gvectors(system%primitive_b, epsilon_cutoff, ngmax_t2, &
                        ng_t7, gvec_t7, gg_t7)

    allocate(sigc_w7(ib_min:ib_max, system%nk))
    allocate(zfac_w7(ib_min:ib_max, system%nk))
    allocate(vxc_w7 (ib_min:ib_max, system%nk))
    allocate(eqp_w7 (ib_min:ib_max, system%nk))
    allocate(sigx_w7(ib_min:ib_max, system%nk))

    skipfrac7 = 0.0d0

    ! Node scaling: replicate the orbitals to every rank so the output-state
    ! k-loop inside calc_sigma_x / calc_sigma_gpp can be split over info%icomm_k
    ! (each rank owns a k-share; matrix elements / dielectric are then local and
    ! a single reduction assembles the per-k self-energy).  With one rank this is
    ! a copy and reproduces the serial result.  calc_vxc_expect keeps the
    ! distributed spsi (it already guards by k and assembles collectively).
    call replicate_orbitals_k(system, info, spsi, spsi_full)

    ! Optional exact point-group symmetrisation of the (now full-mesh) orbitals
    ! by grid rotation: every star member is overwritten by the rotated
    ! representative so the dielectric matrix is symmetric to machine precision,
    ! which the symmetry-reduced BZ sum relies on.  Needs sym.dat in the run dir.
    if (yn_gw_sym == 'y') then
      call gw_symmetrize_orbitals(system, info, lg, spsi_full, energy)
    end if

    ! Resolve the eps/sigma band caps (0=all, >=nocc enforced inside) once.
    ! nocc = highest occupied band over ALL k and spin (safe for metals/magnetic).
    nocc_t7 = 0
    do is_t7 = 1, system%nspin
      do ik_t7 = 1, system%nk
        do io_t7 = 1, system%no
          if (system%rocc(io_t7,ik_t7,is_t7) > 1.0d-6) nocc_t7 = max(nocc_t7, io_t7)
        end do
      end do
    end do
    call resolve_band_caps(system%no, nocc_t7, nband_eps, nband_sigma, neps_t7, nsig_t7)
    ! chi0 extrapolar: Phase 1 implements 'offset' (constant DeltaE) only; the
    ! 'sumrule' (EET) effective energy is Phase 2 and not wired yet.  The tail is
    ! also a no-op unless nband_eps actually caps the conduction sum (neps<no).
    if (yn_gw_extrapolar == 'y') then
      if (trim(gw_extrapolar_mode) /= 'offset') then
        if (comm_is_root(nproc_id_global)) &
          write(*,*) "  [gw] FATAL: gw_extrapolar_mode='", trim(gw_extrapolar_mode), &
                     "' not implemented (offset only)"
        stop 'gw: extrapolar mode unsupported'
      end if
      if (neps_t7 >= system%no .and. comm_is_root(nproc_id_global)) &
        write(*,*) "  [gw] WARN: yn_gw_extrapolar='y' but nband_eps caps no bands -> tail is a no-op"
      ! Real-axis Sigma_c (gpp_real) integrates the spectral function -Im[W^c]/pi
      ! over [0, omega_max_gw].  The extrapolar lumps the missing empties into one
      ! oscillator at Ebar-eps_v = esp(neps)+DeltaE - eps_v; if that pole sits above
      ! the grid its Im weight is lost while its Re screening still suppresses the
      ! on-grid plasmon -> anti-screening (gap opens).  Auto-extend the grid (keeping
      ! the spacing) so the lumped pole is represented.  The static 'gpp' path uses a
      ! sum-rule plasmon (no grid) and needs no adjustment.
      if (neps_t7 < system%no .and. trim(sigma_type) == 'gpp_real') then
        ebar_max_t7 = -1.0d30;  ev_min_t7 = 1.0d30
        do is_t7 = 1, system%nspin
          do ik_t7 = 1, system%nk
            ebar_max_t7 = max(ebar_max_t7, energy%esp(neps_t7, ik_t7, is_t7))
            do io_t7 = 1, nocc_t7
              if (system%rocc(io_t7,ik_t7,is_t7) > 1.0d-6) &
                ev_min_t7 = min(ev_min_t7, energy%esp(io_t7,ik_t7,is_t7))
            end do
          end do
        end do
        ! needed window (eV): deepest lumped pole + a few broadening widths of margin
        om_need_t7 = (ebar_max_t7 + gw_extrapolar_de - ev_min_t7)*27.211386d0 &
                     + max(10.0d0, 8.0d0*eta_gw)
        if (omega_max_gw > 1.0d-6 .and. omega_max_gw < om_need_t7) then
          ! grid spacing is omega_max_gw/(nomega_gw-1); scale intervals to keep it
          nomega_gw    = ceiling((nomega_gw-1) * om_need_t7 / omega_max_gw) + 1
          if (comm_is_root(nproc_id_global)) &
            write(*,'(A,F7.2,A,F7.2,A,I0,A)') &
              "   [gw] extrapolar: omega_max_gw ", omega_max_gw, &
              " eV too small for the lumped pole; auto-extended to ", om_need_t7, &
              " eV (nomega_gw=", nomega_gw, ")"
          omega_max_gw = om_need_t7
        end if
      end if
    end if
    ! Band separation is only wired on the qcache path (gpp / gpp_real, no sym).
    if ((nband_eps > 0 .or. nband_sigma > 0) .and. &
        (yn_gw_sym == 'y' .or. yn_gw_qcache /= 'y')) then
      if (comm_is_root(nproc_id_global)) &
        write(*,*) "  [gw] FATAL: band separation needs yn_gw_qcache='y', yn_gw_sym='n'"
      stop 'gw: band separation unsupported on this path'
    end if

    do is_t7 = 1, system%nspin
      ! bare exchange Sigma_x (the exchange part of Sigma = Sigma_x + Sigma_c).
      call calc_sigma_x(system, info, mg, lg, spsi_full, gvec_t7, gg_t7, ng_t7, &
                        is_t7, ib_min, ib_max, sigx_w7, local_only=.true.)
      ! dynamic correlation Re Sigma_c(eps^KS) and the renormalization Z.
      ! yn_gw_qcache='y' selects the q-cached BZ sum (O(nk^2), identical result);
      ! both take the replicated orbitals built above.
      if (trim(sigma_type) == 'gpp_real') then
        ! b1-4: full-frequency real-axis Sigma_c from the chi0(w)->W(w) engine
        ! (no plasmon-pole, no pole-skip).  Build the omega' grid (0..omega_max_gw)
        ! that calc_chi0_freq integrates over, the same grid the absorption uses.
        if (.not. allocated(omega_grid_t7)) allocate(omega_grid_t7(nomega_gw))
        do iw_t7 = 1, nomega_gw
          omega_grid_t7(iw_t7) = (dble(iw_t7-1)/dble(max(nomega_gw-1,1))) &
                                 * omega_max_gw / 27.211386d0
        end do
        eta_t7 = eta_gw / 27.211386d0
        ! yn_gw_qcache='y' caches W per distinct q + distributes over icomm_k
        ! (node-scalable, identical result); else the base path (all (ik,iq), -n1).
        if (yn_out_gw_spectral == 'y') then
          ! sp3: scan Sigma_c(in, k=1; w) over a wide window ( band edges +- ~25 eV,
          ! to reach the plasmon satellites) accumulated in the SAME q-cache pass.
          nw_scan_t7 = 400
          if (.not. allocated(w_scan_t7)) allocate(w_scan_t7(nw_scan_t7))
          if (.not. allocated(sigc_scan_t7)) allocate(sigc_scan_t7(ib_min:ib_max, nw_scan_t7))
          do iws_t7 = 1, nw_scan_t7
            w_scan_t7(iws_t7) = ( energy%esp(ib_min,1,is_t7) - 25.0d0/27.211386d0 ) &
              + dble(iws_t7-1)/dble(nw_scan_t7-1) &
              * ( (energy%esp(ib_max,1,is_t7) + 25.0d0/27.211386d0) &
                - (energy%esp(ib_min,1,is_t7) - 25.0d0/27.211386d0) )
          end do
          call calc_sigma_c_real_qcache(system, info, mg, lg, spsi_full, energy%esp, &
                          gvec_t7, gg_t7, ng_t7, is_t7, ib_min, ib_max, &
                          nomega_gw, omega_grid_t7, eta_t7, sigc_w7, zfac_w7, &
                          nw_scan=nw_scan_t7, w_scan=w_scan_t7, k_scan=1, sigc_scan=sigc_scan_t7, &
                          nb_sigma=nsig_t7, nb_eps=neps_t7, &
                          do_remainder=(yn_gw_static_remainder=='y'))
        else if (yn_gw_sym == 'y') then
          ! point-group symmetry-reduced (q-IBZ + output-k IBZ); needs the
          ! symmetrised orbitals built above + sym.dat in the run dir.
          call calc_sigma_c_real_sym(system, info, mg, lg, spsi_full, energy%esp, &
                          gvec_t7, gg_t7, ng_t7, is_t7, ib_min, ib_max, &
                          nomega_gw, omega_grid_t7, eta_t7, sigc_w7, zfac_w7)
        else if (yn_gw_qcache == 'y') then
          call calc_sigma_c_real_qcache(system, info, mg, lg, spsi_full, energy%esp, &
                          gvec_t7, gg_t7, ng_t7, is_t7, ib_min, ib_max, &
                          nomega_gw, omega_grid_t7, eta_t7, sigc_w7, zfac_w7, &
                          nb_sigma=nsig_t7, nb_eps=neps_t7, &
                          do_remainder=(yn_gw_static_remainder=='y'))
        else
          call calc_sigma_c_real(system, info, mg, lg, spsi_full, energy%esp, &
                            gvec_t7, gg_t7, ng_t7, is_t7, ib_min, ib_max, &
                            nomega_gw, omega_grid_t7, eta_t7, &
                            sigc_w7, zfac_w7, local_only=.true.)
        end if
        skipfrac7 = 0.0d0   ! real-axis has no unphysical-pole skip
      else if (yn_gw_sym == 'y') then
        ! point-group symmetry-reduced sum (needs the symmetrised orbitals above).
        call calc_sigma_gpp_sym(system, info, mg, lg, spsi_full, energy%esp, rho, &
                            gvec_t7, gg_t7, ng_t7, is_t7, ib_min, ib_max, &
                            sigc_w7, zfac_w7, skip_frac=skipfrac7)
      else if (yn_gw_qcache == 'y') then
        call calc_sigma_gpp_qcache(system, info, mg, lg, spsi_full, energy%esp, rho, &
                            gvec_t7, gg_t7, ng_t7, is_t7, ib_min, ib_max, &
                            sigc_w7, zfac_w7, skip_frac=skipfrac7, &
                            do_remainder=(yn_gw_static_remainder=='y'), &
                            nb_sigma=nsig_t7, nb_eps=neps_t7)
      else
        call calc_sigma_gpp(system, info, mg, lg, spsi_full, energy%esp, rho, &
                            gvec_t7, gg_t7, ng_t7, is_t7, ib_min, ib_max, &
                            sigc_w7, zfac_w7, skip_frac=skipfrac7, local_only=.true.)
      end if
      call calc_vxc_expect(system, info, mg, spsi, Vxc, is_t7, &
                           ib_min, ib_max, vxc_w7)

      ! Linearised one-shot QP update: eqp = eks + Z (sigx + sigc - vxc).
      call solve_qp(energy%esp(ib_min:ib_max,:,is_t7), sigx_w7, sigc_w7, &
                    vxc_w7, zfac_w7, eqp_w7)

      ! Scatter into the full output arrays for this spin.
      sigx   (ib_min:ib_max,:,is_t7) = sigx_w7(:,:)   ! bare exchange
      sigc   (ib_min:ib_max,:,is_t7) = sigc_w7(:,:)   ! Re Sigma_c
      zfac   (ib_min:ib_max,:,is_t7) = zfac_w7(:,:)
      vxc_arr(ib_min:ib_max,:,is_t7) = vxc_w7 (:,:)
      eqp    (ib_min:ib_max,:,is_t7) = eqp_w7 (:,:)
    end do

    ! ----------------------------------------------------------------
    ! [gw][t7] sanity block (root only): for spin 1, k=1, report the top
    ! valence and bottom conduction states.  Reports Sigma_x, Re Sigma_c, Z,
    ! the direct gap@Gamma (k=1) KS vs G0W0, and the unphysical-pole skip
    ! fraction.  Asserts the PoC ordering KS < G0W0 < COHSEX(5.39 eV).
    ! ----------------------------------------------------------------
    is_t7 = 1
    ik_t7 = 1
    ivtop7 = 0
    do io = 1, system%no
      if (system%rocc(io,ik_t7,is_t7) > 1.0d-6) ivtop7 = io
    end do
    icbot7 = min(ivtop7 + 1, system%no)

    if (comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "[gw][t7] dynamic G0W0 self-energy (generalized plasmon-pole)"
      write(*,*) "[gw][t7] ng =", ng_t7, "  band window:", ib_min, "..", ib_max
      write(*,'(A,I4,A,I4)') "  [gw][t7] at k=1: top valence n=", ivtop7, &
        "  bottom conduction n=", icbot7
      write(*,'(A,ES12.4)') "  [gw][t7] unphysical-pole skip fraction =", skipfrac7

      if (ivtop7 >= ib_min .and. ivtop7 <= ib_max) then
        write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t7] top valence: <Sig_x>=", &
          sigx(ivtop7,ik_t7,is_t7)*hartree2ev, " eV   Re<Sig_c>=", &
          sigc(ivtop7,ik_t7,is_t7)*hartree2ev, " eV"
        write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t7]   Z=", &
          zfac(ivtop7,ik_t7,is_t7), "      <Vxc>=", &
          vxc_arr(ivtop7,ik_t7,is_t7)*hartree2ev, " eV"
      end if
      if (icbot7 >= ib_min .and. icbot7 <= ib_max) then
        write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t7] bot conduction: <Sig_x>=", &
          sigx(icbot7,ik_t7,is_t7)*hartree2ev, " eV   Re<Sig_c>=", &
          sigc(icbot7,ik_t7,is_t7)*hartree2ev, " eV"
        write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t7]   Z=", &
          zfac(icbot7,ik_t7,is_t7), "      <Vxc>=", &
          vxc_arr(icbot7,ik_t7,is_t7)*hartree2ev, " eV"
      end if
      ! direct gap@Gamma (k=1): KS vs one-shot G0W0
      if (ivtop7 >= 1 .and. icbot7 <= system%no .and. icbot7 > ivtop7) then
        gap_ks7 = ( energy%esp(icbot7,ik_t7,is_t7) - energy%esp(ivtop7,ik_t7,is_t7) ) &
                  * hartree2ev
        if (icbot7 >= ib_min .and. icbot7 <= ib_max .and. &
            ivtop7 >= ib_min .and. ivtop7 <= ib_max) then
          gap_g0w0_7 = ( eqp(icbot7,ik_t7,is_t7) - eqp(ivtop7,ik_t7,is_t7) ) * hartree2ev
          write(*,'(A,F12.5,A,F12.5,A)') "  [gw][t7] direct gap @k=1: KS=", &
            gap_ks7, " eV   G0W0=", gap_g0w0_7, " eV"
          write(*,'(A,L2)') "  [gw][t7]   ordering KS < G0W0 ?         ", &
            (gap_g0w0_7 > gap_ks7)
          write(*,'(A,L2,A)') "  [gw][t7]   ordering G0W0 < COHSEX(5.39)? ", &
            (gap_g0w0_7 < 5.39d0), "  (PoC: KS < G0W0 < COHSEX)"
        else
          write(*,'(A,F12.5,A)') "  [gw][t7] direct gap @k=1: KS=", gap_ks7, &
            " eV  (G0W0 gap needs both gap states inside the band window)"
        end if
      end if
      ! sp3: write Sigma_c(n,k=1;w) + spectral function A(n,k;w) for the band edges.
      ! Im Sigma_c(w) is the scattering rate; A peaks at the QP energy and shows
      ! the plasmon satellites -- the research-goal output of the real-axis engine.
      if (yn_out_gw_spectral == 'y' .and. allocated(sigc_scan_t7)) then
        open(newunit=fh_spec, file='Si_sigma_c_spectrum.data', status='replace')
        write(fh_spec,'(A)') '# real-axis Sigma_c(n,k=1;w) and spectral function A(n,k;w)'
        write(fh_spec,'(A)') '# 1:w[eV]  then per edge (VBM,CBM): ReSigc[eV] ImSigc[eV] A[1/eV]'
        do iws_t7 = 1, nw_scan_t7
          write(fh_spec,'(F12.5)',advance='no') w_scan_t7(iws_t7)*hartree2ev
          do io_spec = ivtop7, icbot7, max(icbot7-ivtop7,1)
            ek7   = energy%esp(io_spec, ik_t7, is_t7)
            sx7   = sigx   (io_spec, ik_t7, is_t7)
            vx7   = vxc_arr(io_spec, ik_t7, is_t7)
            resc7 = dble (sigc_scan_t7(io_spec, iws_t7))
            imsc7 = aimag(sigc_scan_t7(io_spec, iws_t7))
            aval7 = (abs(imsc7)/acos(-1.0d0)) &
                  / ( (w_scan_t7(iws_t7) - ek7 - sx7 - resc7 + vx7)**2 + imsc7**2 )
            write(fh_spec,'(3ES15.6)',advance='no') &
              resc7*hartree2ev, imsc7*hartree2ev, aval7/hartree2ev
          end do
          write(fh_spec,*)
        end do
        close(fh_spec)
        write(*,'(A,I4,A,I4)') "  [gw][spectral] wrote Si_sigma_c_spectrum.data: VBM n=", &
          ivtop7, "  CBM n=", icbot7
      end if
      write(*,*)
    end if

    deallocate(sigc_w7, zfac_w7, vxc_w7, eqp_w7, sigx_w7)
    if (allocated(w_scan_t7))    deallocate(w_scan_t7)
    if (allocated(sigc_scan_t7)) deallocate(sigc_scan_t7)
    deallocate(gvec_t7, gg_t7)
    if (allocated(spsi_full%zwf)) deallocate(spsi_full%zwf)
    if (allocated(spsi_full%rwf)) deallocate(spsi_full%rwf)

    if (comm_is_root(nproc_id_global)) then
      write(*,*) "  [gw] computed dynamic G0W0 QP energies (sigma_type=gpp)"
    end if

  else
    if (comm_is_root(nproc_id_global)) then
      write(*,*) "  [gw] sigma_type='"//trim(sigma_type)// &
                 "' not a self-energy rung here: QP = KS passthrough"
    end if
  end if

  ! ----------------------------------------------------------------
  ! Step 6: write output
  ! ----------------------------------------------------------------
  call write_qp_energies(system, energy%esp, eqp, zfac, sigx, sigc, vxc_arr)
  call write_qp_dos(system, energy%esp, eqp)

  if (comm_is_root(nproc_id_global)) then
    write(*,*) "  [gw] wrote QP energies + QP/KS DOS"
    write(*,*) "  [gw] band window: ", ib_min, " to ", ib_max
  end if

  deallocate(eqp, zfac, sigx, sigc, vxc_arr)

end subroutine main_gw
