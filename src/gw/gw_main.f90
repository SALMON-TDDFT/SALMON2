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
  use gw_coulomb_sub,     only: build_gvectors, build_vcoul
  use gw_mtxel_sub,       only: calc_mtxel
  use gw_epsilon_sub,     only: calc_epsinv
  use gw_sigma_x_sub,     only: calc_sigma_x
  use gw_qp_sub,          only: calc_vxc_expect, solve_qp
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
  ! t2/t3/t4 sanity blocks are validated; compile them out of production runs
  ! (set .true. to re-enable the per-stage diagnostics).
  logical, parameter :: run_sanity_t234 = .false.

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

  if (comm_is_root(nproc_id_global)) then
    write(*,*) "  [gw] wrote QP energies"
    write(*,*) "  [gw] band window: ", ib_min, " to ", ib_max
  end if

  deallocate(eqp, zfac, sigx, sigc, vxc_arr)

end subroutine main_gw
