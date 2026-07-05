module ttm3
!
! Three-temperature + carrier model (Stage 1: local rate-equation engine).
!
! State per cell: electron temperature Te, hole temperature Th, lattice
! temperature Tl, electron density Ne, hole density Nh.  Transport is OFF in
! this stage (local ODEs only); the sendrecv-grid argument of ttm3_main is
! reserved for the transport stage.  Maxwell coupling is also a later stage:
! here the heating source and carrier generation are inputs.
!
! Recommended first model (documented for review):
!   * classical carrier heat capacity C = (3/2) N  (atomic units, kB=1),
!   * energy-conserving carrier<->lattice relaxation with time constant tau,
!   * standard Auger recombination R = (A_e Ne + A_h Nh) Ne Nh (removes e-h
!     pairs, so dNe=dNh for recombination).
! The exact reference heat capacities (from Fermi statistics) and the exact
! relaxation/Auger split are a coefficient-level refinement for a later pass;
! they do not change this module's structure or interfaces.
!
  implicit none
  private
  public :: init_ttm3_parameters
  public :: ttm3_set_dt
  public :: init_ttm3_grid
  public :: init_ttm3_alloc
  public :: ttm3_main
  public :: ttm3_step_cell
  public :: ttm3_get_state
  public :: ttm3_get_max
  public :: ttm3_get_front
  public :: ttm3_front_ijk
  public :: ttm3_set_xcomm
  public :: ttm3_write_profile
  public :: ttm3_permittivity
  public :: ttm3_linear_gen
  public :: ttm3_drude_coef
  public :: ttm3_drude_coef_cell
  public :: ttm3_drude_osc_cell
  public :: ttm3_eps_sig
  public :: ttm3_coupling_mode
  public :: ttm3_ninterior
  public :: ttm3_interior_cell
  public :: ttm3_generation
  public :: ttm3_gap

  logical,public :: use_ttm3=.false.

  ! material parameters
  type ttm3_param
     real(8) :: Egap     ! band gap (cold)                 [Hartree]
     real(8) :: mu_e     ! electron effective mass ratio    [-]
     real(8) :: mu_h     ! hole effective mass ratio        [-]
     real(8) :: A_e      ! electron Auger coefficient        [bohr^6 / a.u.time]
     real(8) :: A_h      ! hole Auger coefficient            [bohr^6 / a.u.time]
     real(8) :: tau      ! carrier-lattice relaxation time   [a.u.time]
     real(8) :: Cl       ! lattice heat capacity             [a.u. / (volume * temperature)]
     real(8) :: Tini     ! initial temperature (Te=Th=Tl)    [a.u.]
     ! stage 2-5 (recommended models)
     real(8) :: eps_bg   ! background relative permittivity  [-]            (Stage 2)
     real(8) :: N0       ! saturation/reference density      [1/bohr^3]     (Stage 2,5)
     real(8) :: beta2    ! two-photon generation coefficient [a.u.]         (Stage 3)
     real(8) :: Ddiff    ! ambipolar carrier diffusivity     [bohr^2/a.u.t] (Stage 4)
     real(8) :: kappa_e  ! carrier thermal conductivity      [a.u.]         (Stage 4)
     real(8) :: kappa_l  ! lattice thermal conductivity      [a.u.]         (Stage 4)
     real(8) :: dgap_c   ! gap-renormalisation coefficient   [a.u.]         (Stage 5)
     ! Set A (Fermi transport): carrier mobilities (optional 16th,17th input lines).
     ! When >0 the diffusion/conductivity are computed from the Fermi-Dirac integrals
     ! (mu=mob0*F0/F12, D=mob0*Te*F0/Fm, K=(6 F2/F0-4(F1/F0)^2) N Te mu), superseding
     ! the constant Ddiff/kappa_e.  Input m^2/Vs, converted to a.u. at read.
     real(8) :: mob_e0   ! electron mobility                 [a.u. = m^2/Vs * 2.3505e+5]
     real(8) :: mob_h0   ! hole mobility                     [a.u.]
     ! Layer C-a (ablation): cold bound (interband) conductivity (optional 18th input line).
     ! When >0, carrier generation is restricted to the BOUND absorption fraction
     ! sig_bound/(sig_bound+sig_drude) with sig_bound=sig_cold*(1-Ne/N0) (band filling), so
     ! free-carrier (Drude) absorption heats rather than ionises -> generation self-limits in
     ! the metallisation regime.  <=0 disables the split (generation = full absorbed power).
     real(8) :: sig_cold ! cold interband conductivity       [a.u.]         (Layer C-a)
     integer :: coupling_mode ! Maxwell-coupling method (&ttm ttm_coupling):
                              ! 0=J-update (current-source ADE, bilinear), 1=eps-update (carrier
                              ! eps/sigma folded into the FDTD media coefficients; diverges when
                              ! eps_re -> 0/negative = deep metallization), 2=J-update with ETD
                              ! integrator, 3=exact 2x2 matrix-exponential JE (unconditionally
                              ! stable, energy-conserving -- the DEFAULT).
  end type ttm3_param

  integer,allocatable :: ijk_media_whole(:,:)
  integer,allocatable :: ijk_media_myrnk(:,:)
  integer,allocatable :: ijk_interior(:,:)   ! interior medium cells (all 6 neighbours same medium)
  integer :: nmedia_myrnk
  integer :: ninterior
  integer :: is_array(3), ie_array(3)
  integer :: is_inner(3), ie_inner(3)   ! inner (no-halo) bounds, for the space-charge prefix sum
  integer :: comm
  integer :: nprocs_g = 1                            ! # MPI ranks (diagnostic global reductions)
  integer :: xc_icomm = 0, xc_id = 0, xc_isize = 1   ! x-axis sub-comm (space-charge x-decomposition)
  integer :: gx0 = 0, gx1 = -1                       ! global inner x-range (depth-profile gather)
  ! medium surface faces (medium cell whose neighbour across face dir is not medium):
  ! ttm3_transport mirrors state across these faces = reflecting (zero-flux) boundary,
  ! matching the standalone 3D3TM (a fixed vacuum neighbour would act as a Dirichlet sink).
  integer,allocatable :: ijk_surf(:,:)               ! (4,nsurf): ix,iy,iz, face dir (+-1,+-2,+-3)
  integer :: nsurf = 0
  ! medium mask over the full (halo) box: the space-charge sum must skip vacuum cells,
  ! whose Ne/Nh now hold mirrored (nonzero) values from the reflecting boundary.
  logical,allocatable :: med_mask(:,:,:)
  real(8) :: dcap_worst = 1.0d0                      ! worst D/Dmax seen (const path CFL-cap diagnostic)
  integer :: nsub_worst = 8                          ! worst transport sub-step count reported so far
  logical :: ktw_warned = .false.                    ! one-shot: coefficient kT floor engaged
  real(8) :: hgs(3)
  real(8) :: dt

  type(ttm3_param) :: tp3

  ! state and Runge-Kutta stage buffers
  real(8),allocatable :: Te(:,:,:), Th(:,:,:), Tl(:,:,:)
  real(8),allocatable :: Ne(:,:,:), Nh(:,:,:)
  ! transport increment buffers (Stage 4)
  real(8),allocatable :: rhs_te(:,:,:), rhs_th(:,:,:), rhs_tl(:,:,:)
  real(8),allocatable :: rhs_ne(:,:,:), rhs_nh(:,:,:)
  ! Set A (Fermi transport) per-cell coefficient fields, frozen over the sub-steps of one
  ! FDTD step: thermal conductivity Kc, diffusion Dc, mobility Mc, heat capacity Cf,
  ! heat-per-carrier Cj (cJeh), and the N-normalized cross-advection velocities
  ! bE = Dc*cjE/N (band-edge force) and bT = Dc*cjT/N (thermodiffusion) -- these are
  ! bounded ~Dc/kT, so no N_neighbor/N_local amplification survives at density shoulders
  ! (the 2026-07-06 s2 blow-up mechanism).  Fluxes are built at faces (flux form).
  real(8),allocatable :: Kc_e(:,:,:), Kc_h(:,:,:)
  real(8),allocatable :: Dc_e(:,:,:), Dc_h(:,:,:), Mc_e(:,:,:), Mc_h(:,:,:)
  real(8),allocatable :: Cf_e(:,:,:), Cf_h(:,:,:)
  real(8),allocatable :: bE_e(:,:,:), bT_e(:,:,:), bE_h(:,:,:), bT_h(:,:,:)
  real(8),allocatable :: rgap(:,:,:), Efld(:,:,:)       ! local gap, space-charge field (x, cell centre)
  real(8),allocatable :: Efc(:,:,:)                     ! space-charge field at x-faces (i+1/2, exact prefix)
  real(8),allocatable :: Cj_e(:,:,:), Cj_h(:,:,:)       ! Set A layer 3: heat per carrier (cJeh)

  ! Fermi-Dirac integral table F_j(eta).  Indices 1,2,3 = j = -1/2,1/2,3/2 (heat
  ! capacity / chemical potential); indices 4,5 = j = 1,2 (transport: mobility,
  ! diffusion, thermal conductivity).  F_0(eta)=ln(1+exp(eta)) is analytic (ttm3_F0),
  ! no table.  Built once at init by quadrature; interpolated in the inner loop (the
  ! same table approach as the reference's Table_Fermi).  eta grid: f_eta0 + k*f_deta.
  integer,parameter   :: f_neta = 2400
  real(8),parameter   :: f_eta0 = -60.0d0, f_deta = 0.05d0
  integer,parameter   :: f_nj   = 5
  ! j values and Gamma(j+1) for indices 1..5 = j = -1/2, 1/2, 3/2, 1, 2 (shared by the
  ! table builder and the degenerate (eta > +60) Sommerfeld tail of ttm3_fermi).
  real(8),parameter   :: f_jval(f_nj) = (/ -0.5d0, 0.5d0, 1.5d0, 1.0d0, 2.0d0 /)
  real(8),parameter   :: f_gam1(f_nj) = (/ 1.7724538509055160d0, 0.8862269254527580d0, &
                                           1.3293403881791370d0, 1.0d0, 2.0d0 /)
  real(8)             :: f_tab(0:f_neta,f_nj)
  logical             :: f_built = .false.

  ! density floor to keep the carrier heat capacity well defined.  Set to the
  ! standalone 3D3TM's thermal intrinsic seed density Ns(300 K) so the carrier
  ! population can reach the reference's near-transparent regime; this is
  ! safe only because the carrier temperatures use the exponential integrator
  ! (the old RK4 step blew up at such low Ne -- see ttm3_step_cell).
  real(8),parameter :: N_floor = 3.657d-14  ! [1/bohr^3] (=2.468e11 cm^-3, standalone 3D3TM intrinsic seed Ns)
  real(8),parameter :: T_clamp = 31.67d0   ! [a.u.] (~1e7 K) numerical temperature cap
  ! positivity floor for the temperatures (~1 K).  The explicit transport step can
  ! undershoot at steep (metallization) fronts; a negative temperature would flip the
  ! heat equation backward and turn exp(-1.5*gap/T) into an overflow -> NaN cascade.
  real(8),parameter :: T_floor = 3.17d-6   ! [a.u.] (~1 K)
  real(8),parameter :: pi_ = 3.14159265358979323846d0
  real(8),parameter :: hk_kelvin = 3.1577502480407d5   ! a.u. temperature -> Kelvin
  real(8),parameter :: T_au = 2.48036d5                ! [a.u.time] Auger saturation timescale

  character(12) :: ttm3_file = 'ttm3.inp_3tm'
  logical :: DISPLAY=.false.

contains

  !---------------------------------------------------------------------------
  subroutine init_ttm3_parameters( dt_em )
    use communication, only: comm_get_globalinfo, comm_is_root, comm_bcast
    use salmon_global, only: theory, yn_use_ttm3, ttm_egap, ttm_mu_e, ttm_mu_h, ttm_auger_e, ttm_auger_h, &
                             ttm_tau, ttm_cl, ttm_tini, ttm_eps_bg, ttm_n0, ttm_beta2, &
                             ttm_ddiff, ttm_kappa_e, ttm_kappa_l, ttm_dgap, &
                             ttm_mob_e, ttm_mob_h, ttm_sig_cold, ttm_coupling
    implicit none
    real(8), intent(in) :: dt_em
    integer, parameter :: unit=1112
    integer :: ios
    integer :: npid, nprocs
    logical :: flag
    real(8), parameter :: atomic_unit_of_length = 5.29177210903d-11 ! [m]
    real(8), parameter :: hartree_joule_relationship = 4.3597447222071d-18 ! [J]
    real(8), parameter :: hartree_kelvin_relationship = 3.1577502480407d5 ! [K]
    real(8), parameter :: atomic_unit_of_time = 2.4188843265857d-17 ! [s]
    real(8), parameter :: hartree_ev = 27.211386245988d0 ! [eV]

    call comm_get_globalinfo( comm, npid, nprocs )
    nprocs_g = nprocs
    DISPLAY = comm_is_root(npid)

    if ( DISPLAY ) write(*,'(a60)') repeat("-",30)//" init_ttm3_parameters(start)"

    use_ttm3 = yn_use_ttm3
    if( .not.use_ttm3 )then
       if( DISPLAY )then
          write(*,*) "theory /= 'maxwell-3tm': three-temperature + carrier model is not used."
          write(*,'(a60)') repeat("-",30)//" init_ttm3_parameters(end  )"
       end if
       return
    end if
    if( DISPLAY ) write(*,*) "theory = 'maxwell-3tm': three-temperature + carrier model (&ttm)."

    dt = dt_em

    ! &ttm namelist parameters (read + broadcast by inputoutput; physical units, converted below)
    tp3%Egap   = ttm_egap    ; tp3%mu_e   = ttm_mu_e   ; tp3%mu_h   = ttm_mu_h
    tp3%A_e    = ttm_auger_e ; tp3%A_h    = ttm_auger_h
    tp3%tau    = ttm_tau     ; tp3%Cl     = ttm_cl     ; tp3%Tini   = ttm_tini
    tp3%eps_bg = ttm_eps_bg  ; tp3%N0     = ttm_n0     ; tp3%beta2  = ttm_beta2
    tp3%Ddiff  = ttm_ddiff   ; tp3%kappa_e= ttm_kappa_e; tp3%kappa_l= ttm_kappa_l
    tp3%dgap_c = ttm_dgap    ; tp3%mob_e0 = ttm_mob_e  ; tp3%mob_h0 = ttm_mob_h
    tp3%sig_cold = ttm_sig_cold ; tp3%coupling_mode = ttm_coupling
    if( DISPLAY )then
       write(*,*) "Egap[eV]    =",tp3%Egap
       write(*,*) "mu_e/mu_h   =",tp3%mu_e, tp3%mu_h
       write(*,*) "tau[fs]     =",tp3%tau
       write(*,*) "eps_bg / N0 =",tp3%eps_bg, tp3%N0
       write(*,*) "mob_e/mob_h =",tp3%mob_e0, tp3%mob_h0
       write(*,*) "sig_cold/cm =",tp3%sig_cold, tp3%coupling_mode
    end if

    if( tp3%coupling_mode < 0 .or. tp3%coupling_mode > 3 )then
       if( DISPLAY ) write(*,'(A,I6)') " ERROR: ttm_coupling must be 0..3, got ", tp3%coupling_mode
       error stop 'invalid ttm_coupling'
    end if

! Convert to atomic units
    tp3%Egap = tp3%Egap / hartree_ev
    ! mobility [m^2/Vs] -> a.u.: 1 m^2/Vs = 1e20 A^2/(V s) = 1e5 A^2/(V fs), then x E_h*a_t/a_B^2
    ! => mob_au = mob_SI * 2.3505e5 (8.5e-3 m^2/Vs = 85 cm^2/Vs -> 1998 a.u.).  Cross-check via the
    ! Einstein relation with this file's own Ddiff conversion: mu*kT(300 K) = 1998*9.5e-4 = 1.90
    ! bohr^2/aut = 2.2 cm^2/s, the textbook ambipolar D of Si.  (The standalone reference's init
    ! conversion, 8.5e-3 -> 2.0e-7, is 1e10 too small = its Fermi transport is effectively inert;
    ! the port used to reproduce that number verbatim.)
    tp3%mob_e0 = tp3%mob_e0 * hartree_ev*(atomic_unit_of_time*1.0d15)*1.0d+5/(atomic_unit_of_length*1.0d10)**2
    tp3%mob_h0 = tp3%mob_h0 * hartree_ev*(atomic_unit_of_time*1.0d15)*1.0d+5/(atomic_unit_of_length*1.0d10)**2
    ! Auger coefficient [cm^6/s] -> [bohr^6 / a.u.time]
    tp3%A_e  = tp3%A_e * (1.0d-2/atomic_unit_of_length)**6 * atomic_unit_of_time
    tp3%A_h  = tp3%A_h * (1.0d-2/atomic_unit_of_length)**6 * atomic_unit_of_time
    tp3%tau  = tp3%tau * 1.0d-15 / atomic_unit_of_time
    ! lattice heat capacity [J/(m^3 K)] -> a.u. (same as the two-temperature model)
    tp3%Cl   = tp3%Cl / hartree_joule_relationship &
                       * atomic_unit_of_length**3 &
                       * hartree_kelvin_relationship
    tp3%Tini = tp3%Tini / hartree_kelvin_relationship

    ! Fermi-Dirac integral table for the carrier heat capacity (built once)
    call ttm3_build_fermi_table()
    ! mu_e, mu_h, eps_bg, beta2, dgap_c are taken as given (eps_bg,mu dimensionless; beta2,dgap_c in a.u.)
    tp3%N0    = tp3%N0 * (atomic_unit_of_length*1.0d2)**3                 ! [cm^-3] -> [bohr^-3]
    tp3%Ddiff = tp3%Ddiff * (1.0d-2/atomic_unit_of_length)**2 * atomic_unit_of_time   ! [cm^2/s] -> a.u.
    ! thermal conductivity [W/(m K)] -> a.u. (same form as the two-temperature model)
    tp3%kappa_e = tp3%kappa_e / hartree_joule_relationship*atomic_unit_of_length &
                              * hartree_kelvin_relationship*atomic_unit_of_time
    tp3%kappa_l = tp3%kappa_l / hartree_joule_relationship*atomic_unit_of_length &
                              * hartree_kelvin_relationship*atomic_unit_of_time

    if ( DISPLAY ) write(*,'(a60)') repeat("-",30)//" init_ttm3_parameters(end  )"

  end subroutine init_ttm3_parameters

  !---------------------------------------------------------------------------
  ! Set the time step.  init_ttm3_parameters is called early (so use_ttm3 can
  ! feed flag_save) before dt_em is finalised by the CFL condition; this setter
  ! is called afterwards with the correct dt_em.
  subroutine ttm3_set_dt( dt_em )
    implicit none
    real(8),intent(in) :: dt_em
    dt = dt_em
  end subroutine ttm3_set_dt

  !---------------------------------------------------------------------------
  subroutine init_ttm3_grid( hgs_in, is_a, is, ie, imedia )
    use communication, only: comm_get_min, comm_get_max
    implicit none
    real(8),intent(in) :: hgs_in(3)
    integer,intent(in) :: is_a(3), is(3), ie(3)
    integer,intent(in) :: imedia(is_a(1):,is_a(2):,is_a(3):)
    integer :: ix,iy,iz,m,id,jx,jy,jz,ndup
    integer :: dirs(3,6)
    logical :: allmed
    logical,allocatable :: seen(:,:,:)
    real(8) :: g1(1),g2(1)

    ! Ablation mode (sig_cold>0): the surface metallizes, so the carrier-dependent optics
    ! back-action must reach the illuminated face cell too -> include ALL medium cells, not
    ! only the 6-neighbour interior.  Layer B / regression (sig_cold<=0) keeps interior-only.
    allmed = ( tp3%sig_cold > 0.0d0 )

    hgs(:)      = hgs_in(:)
    is_array(:) = is_a(:)
    ie_array(:) = ie(:)
    is_inner(:) = is(:)
    ie_inner(:) = ie(:)

    ! the temperatures/carriers live only on medium cells (imedia /= 0):
    ! build the local medium-cell index list, mirroring the two-temperature setup
    m = 0
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
       if( imedia(ix,iy,iz) /= 0 ) m = m + 1
    end do
    end do
    end do
    nmedia_myrnk = m
    if( allocated(ijk_media_myrnk) ) deallocate(ijk_media_myrnk)
    allocate( ijk_media_myrnk(3,max(m,1)) )
    m = 0
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
       if( imedia(ix,iy,iz) /= 0 )then
          m = m + 1
          ijk_media_myrnk(1:3,m) = (/ix,iy,iz/)
       end if
    end do
    end do
    end do

    ! interior medium cells (all 6 neighbours the same medium): used by the
    ! permittivity back-action so it never touches fs%imedia (which the FDTD
    ! setup deallocates) and keeps the interface cells' init coefficients.
    m = 0
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
       if( imedia(ix,iy,iz)/=0 .and. ( allmed .or. ( &
           imedia(ix+1,iy,iz)==imedia(ix,iy,iz) .and. imedia(ix-1,iy,iz)==imedia(ix,iy,iz) .and. &
           imedia(ix,iy+1,iz)==imedia(ix,iy,iz) .and. imedia(ix,iy-1,iz)==imedia(ix,iy,iz) .and. &
           imedia(ix,iy,iz+1)==imedia(ix,iy,iz) .and. imedia(ix,iy,iz-1)==imedia(ix,iy,iz) ) ) ) m = m + 1
    end do
    end do
    end do
    ninterior = m
    if( allocated(ijk_interior) ) deallocate(ijk_interior)
    allocate( ijk_interior(3,max(m,1)) )
    m = 0
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
       if( imedia(ix,iy,iz)/=0 .and. ( allmed .or. ( &
           imedia(ix+1,iy,iz)==imedia(ix,iy,iz) .and. imedia(ix-1,iy,iz)==imedia(ix,iy,iz) .and. &
           imedia(ix,iy+1,iz)==imedia(ix,iy,iz) .and. imedia(ix,iy-1,iz)==imedia(ix,iy,iz) .and. &
           imedia(ix,iy,iz+1)==imedia(ix,iy,iz) .and. imedia(ix,iy,iz-1)==imedia(ix,iy,iz) ) ) )then
          m = m + 1
          ijk_interior(1:3,m) = (/ix,iy,iz/)
       end if
    end do
    end do
    end do

    ! medium surface faces: medium cell whose face-d neighbour is not medium.  ttm3_transport
    ! mirrors state across these faces (reflecting, zero-flux) so the never-evolved vacuum
    ! neighbours do not act as a Dirichlet heat/carrier sink.  imedia halos are already
    ! exchanged here, so faces at rank boundaries are classified correctly.
    dirs = reshape( (/ 1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1 /), (/3,6/) )
    m = 0
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
       if( imedia(ix,iy,iz) /= 0 )then
          do id=1,6
             if( imedia(ix+dirs(1,id),iy+dirs(2,id),iz+dirs(3,id)) == 0 ) m = m + 1
          end do
       end if
    end do
    end do
    end do
    nsurf = m
    if( allocated(ijk_surf) ) deallocate(ijk_surf)
    allocate( ijk_surf(4,max(m,1)) )
    allocate( seen(is_a(1):ubound(imedia,1), is_a(2):ubound(imedia,2), is_a(3):ubound(imedia,3)) )
    seen = .false. ; ndup = 0 ; m = 0
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
       if( imedia(ix,iy,iz) /= 0 )then
          do id=1,6
             jx=ix+dirs(1,id); jy=iy+dirs(2,id); jz=iz+dirs(3,id)
             if( imedia(jx,jy,jz) == 0 )then
                m = m + 1
                ijk_surf(1,m)=ix ; ijk_surf(2,m)=iy ; ijk_surf(3,m)=iz
                ijk_surf(4,m) = merge( (id+1)/2, -((id+1)/2), mod(id,2)==1 )   ! +-1,+-2,+-3
                if( seen(jx,jy,jz) ) ndup = ndup + 1
                seen(jx,jy,jz) = .true.
             end if
          end do
       end if
    end do
    end do
    end do
    deallocate( seen )
    if( nprocs_g > 1 )then                       ! duplicates can sit on any rank
       g1(1) = dble(ndup) ; call comm_get_max( g1, g2, 1, comm ) ; ndup = nint(g2(1))
    end if
    if( ndup > 0 .and. DISPLAY ) write(*,'(A,I8,A)') &
       " ttm3 WARNING: non-slab geometry --", ndup, " vacuum cells border multiple medium faces;" &
       // " the reflecting boundary is approximate at those corners."

    ! medium mask (full halo box) for the space-charge sum
    if( allocated(med_mask) ) deallocate(med_mask)
    allocate( med_mask(is_a(1):ubound(imedia,1), is_a(2):ubound(imedia,2), is_a(3):ubound(imedia,3)) )
    med_mask(:,:,:) = ( imedia(:,:,:) /= 0 )

    ! global inner x-range (for the depth-profile gather under x-decomposition)
    if( nprocs_g > 1 )then
       g1(1) = dble(is(1)) ; call comm_get_min( g1(1), comm ) ; gx0 = nint(g1(1))
       g1(1) = dble(ie(1)) ; call comm_get_max( g1, g2, 1, comm ) ; gx1 = nint(g2(1))
    else
       gx0 = is(1) ; gx1 = ie(1)
    end if
  end subroutine init_ttm3_grid

  !---------------------------------------------------------------------------
  ! Register the x-axis sub-communicator (icomm_x) for the space-charge prefix
  ! sum under x-direction domain decomposition (nproc_rgrid(1)>1).  Called once
  ! at init from the FDTD driver.  isz<=1 (default) keeps the single-domain path.
  subroutine ttm3_set_xcomm( ic, idx, isz )
    implicit none
    integer,intent(in) :: ic, idx, isz
    xc_icomm = ic; xc_id = idx; xc_isize = isz
  end subroutine ttm3_set_xcomm

  !---------------------------------------------------------------------------
  subroutine init_ttm3_alloc( srg, rg )
    use structures, only: s_rgrid, s_sendrecv_grid
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_rgrid),         intent(in)    :: rg
    integer :: i1,i2,i3,j1,j2,j3

    i1=rg%is_array(1); j1=rg%ie_array(1)
    i2=rg%is_array(2); j2=rg%ie_array(2)
    i3=rg%is_array(3); j3=rg%ie_array(3)

    allocate( Te(i1:j1,i2:j2,i3:j3), Th(i1:j1,i2:j2,i3:j3), Tl(i1:j1,i2:j2,i3:j3) )
    allocate( Ne(i1:j1,i2:j2,i3:j3), Nh(i1:j1,i2:j2,i3:j3) )
    allocate( rhs_te(i1:j1,i2:j2,i3:j3), rhs_th(i1:j1,i2:j2,i3:j3), rhs_tl(i1:j1,i2:j2,i3:j3) )
    allocate( rhs_ne(i1:j1,i2:j2,i3:j3), rhs_nh(i1:j1,i2:j2,i3:j3) )

    Te = tp3%Tini ; Th = tp3%Tini ; Tl = tp3%Tini
    Ne = N_floor  ; Nh = N_floor

    ! Set A (Fermi transport) coefficient/current arrays: ~18 full-grid arrays accessed only
    ! inside the use_fermi branch of ttm3_transport.  Allocate them only when transport is on
    ! (mobility>0) so a non-transport run does not pay their memory.
    if( tp3%mob_e0>0.0d0 .or. tp3%mob_h0>0.0d0 )then
       allocate( Kc_e(i1:j1,i2:j2,i3:j3), Kc_h(i1:j1,i2:j2,i3:j3) )
       allocate( Dc_e(i1:j1,i2:j2,i3:j3), Dc_h(i1:j1,i2:j2,i3:j3) )
       allocate( Mc_e(i1:j1,i2:j2,i3:j3), Mc_h(i1:j1,i2:j2,i3:j3) )
       allocate( Cf_e(i1:j1,i2:j2,i3:j3), Cf_h(i1:j1,i2:j2,i3:j3) )
       allocate( bE_e(i1:j1,i2:j2,i3:j3), bT_e(i1:j1,i2:j2,i3:j3) )
       allocate( bE_h(i1:j1,i2:j2,i3:j3), bT_h(i1:j1,i2:j2,i3:j3) )
       allocate( rgap(i1:j1,i2:j2,i3:j3), Efld(i1:j1,i2:j2,i3:j3), Efc(i1:j1,i2:j2,i3:j3) )
       allocate( Cj_e(i1:j1,i2:j2,i3:j3), Cj_h(i1:j1,i2:j2,i3:j3) )
       Kc_e = 0.0d0  ; Kc_h = 0.0d0
       Dc_e = 0.0d0  ; Dc_h = 0.0d0 ; Mc_e = 0.0d0 ; Mc_h = 0.0d0
       Cf_e = 0.0d0  ; Cf_h = 0.0d0
       bE_e = 0.0d0  ; bT_e = 0.0d0 ; bE_h = 0.0d0 ; bT_h = 0.0d0
       rgap = tp3%Egap ; Efld = 0.0d0 ; Efc = 0.0d0
       Cj_e = 0.0d0  ; Cj_h = 0.0d0
    end if
  end subroutine init_ttm3_alloc

  !---------------------------------------------------------------------------
  ! Build the Fermi-Dirac integral table F_j(eta), j = -1/2, 1/2, 3/2, once.
  ! F_j(eta) = (2/Gamma(j+1)) * Integral_0^inf u^(2j+1)/(exp(u^2-eta)+1) du
  ! (the substitution t=u^2 removes the t^(-1/2) endpoint singularity).
  subroutine ttm3_build_fermi_table()
    implicit none
    integer,parameter :: nq=4000
    integer :: k,q,jj
    real(8) :: eta,umax,du,u,integ,s
    ! indices 1..5 = j = -1/2, 1/2, 3/2, 1, 2  (f_jval/f_gam1 = module constants)
    do k=0,f_neta
       eta = f_eta0 + dble(k)*f_deta
       umax = sqrt( max(eta,0.0d0) + 50.0d0 )
       du = umax/dble(nq)
       do jj=1,f_nj
          integ = 0.0d0
          do q=0,nq
             u = dble(q)*du
             s = u**(2.0d0*f_jval(jj)+1.0d0)/( exp(u*u-eta) + 1.0d0 )
             if(q==0 .or. q==nq)then       ; integ=integ+s
             elseif(mod(q,2)==1)then        ; integ=integ+4.0d0*s
             else                           ; integ=integ+2.0d0*s
             end if
          end do
          f_tab(k,jj) = 2.0d0*integ*du/3.0d0/f_gam1(jj)
       end do
    end do
    f_built = .true.
    if(DISPLAY)then
       write(*,'(a)') " ttm3 Fermi heat-capacity self-test  Ce/N = 15/4 F3/2/F1/2 - 9/4 F1/2/F-1/2 :"
       write(*,'(a,f9.5,a)') "   eta=-20 (non-degenerate): ", &
            3.75d0*ttm3_fermi(3,-20.0d0)/ttm3_fermi(2,-20.0d0)-2.25d0*ttm3_fermi(2,-20.0d0)/ttm3_fermi(1,-20.0d0), &
            "  (classical limit 1.5)"
       write(*,'(a,f9.5)') "   eta=  0                  : ", &
            3.75d0*ttm3_fermi(3,0.0d0)/ttm3_fermi(2,0.0d0)-2.25d0*ttm3_fermi(2,0.0d0)/ttm3_fermi(1,0.0d0)
       write(*,'(a,f9.5,a)') "   eta=+20 (degenerate)     : ", &
            3.75d0*ttm3_fermi(3,20.0d0)/ttm3_fermi(2,20.0d0)-2.25d0*ttm3_fermi(2,20.0d0)/ttm3_fermi(1,20.0d0), &
            "  (Sommerfeld-suppressed)"
    end if
  end subroutine ttm3_build_fermi_table

  ! F_j(eta) by linear interpolation of the table (jidx: 1=-1/2, 2=1/2, 3=3/2, 4=1, 5=2).
  ! Above the table (eta > +60, deep degeneracy -- reached after metallization) a two-term
  ! Sommerfeld expansion replaces the old hard clamp (which froze all F-ratios and
  ! overestimated the carrier heat capacity by up to ~20x at Ne ~ 1e22 cm^-3):
  !   F_j(eta) = eta^{j+1}/Gamma(j+2) * [ 1 + pi^2/6 j(j+1)/eta^2
  !                                        + 7pi^4/360 (j+1)j(j-1)(j-2)/eta^4 ]
  ! (this table's normalization; at eta=60 it matches the quadrature to <1e-3).
  pure function ttm3_fermi( jidx, eta ) result( F )
    implicit none
    integer,intent(in) :: jidx
    real(8),intent(in) :: eta
    real(8) :: F,x,w,jv
    integer :: k
    x = (eta - f_eta0)/f_deta
    if(x <= 0.0d0)then            ; k=0       ; w=0.0d0
    elseif(x >= dble(f_neta))then
       jv = f_jval(jidx)
       F  = eta**(jv+1.0d0)/( (jv+1.0d0)*f_gam1(jidx) ) * &
            ( 1.0d0 + (pi_*pi_/6.0d0)*jv*(jv+1.0d0)/(eta*eta) &
                    + (7.0d0*pi_**4/360.0d0)*(jv+1.0d0)*jv*(jv-1.0d0)*(jv-2.0d0)/(eta**4) )
       return
    else                          ; k=int(x)  ; w=x-dble(k)
    end if
    F = (1.0d0-w)*f_tab(k,jidx) + w*f_tab(k+1,jidx)
  end function ttm3_fermi

  ! F_0(eta) = ln(1+exp(eta)) (analytic; no table).  Overflow-guarded: -> eta for eta>>0.
  pure function ttm3_F0( eta ) result( F )
    implicit none
    real(8),intent(in) :: eta
    real(8) :: F
    if( eta > 40.0d0 )then ; F = eta
    else                   ; F = log( 1.0d0 + exp(eta) )
    end if
  end function ttm3_F0

  ! Reduced chemical potential eta from (T,mass,N): invert N = Neff*F_{1/2}(eta),
  ! Neff = 2*(mass*T/(2*pi))^(3/2).  F_{1/2} is tabulated (monotone in eta), so invert by a
  ! binary search on the table + linear interpolation.  Because the table is piecewise-linear
  ! in eta, this returns the SAME eta as inverting that piecewise-linear F_{1/2} -- i.e. it is
  ! mathematically identical to the old 60-iteration bisection, but ~10x cheaper (no repeated
  ! Fermi-table evaluation in the loop).
  pure function ttm3_chem_pot( T, mc, N ) result( eta )
    implicit none
    real(8),intent(in) :: T,mc,N
    real(8) :: eta,Neff,y,w,g52
    integer :: lo,hi,mid,it
    Neff = 2.0d0*( mc*max(T,1.0d-12)/(2.0d0*pi_) )**1.5d0
    y = N/Neff                                       ! target value of F_{1/2}(eta)
    if( y <= f_tab(0,2) )then
       eta = f_eta0; return
    elseif( y >= f_tab(f_neta,2) )then
       ! deep degeneracy (eta > +60): invert the two-term Sommerfeld F_{1/2}(eta)
       !   F_{1/2} = ( eta^{3/2} + (pi^2/8) eta^{-1/2} ) / Gamma(5/2)
       ! leading order eta = (Gamma(5/2) y)^{2/3}, then two Newton corrections.
       g52 = f_gam1(3)                                ! Gamma(5/2) = 1.32934...
       eta = ( g52*y )**(2.0d0/3.0d0)
       do it=1,2
          eta = eta - ( (eta**1.5d0 + 0.125d0*pi_*pi_/sqrt(eta))/g52 - y ) &
                     /( (1.5d0*sqrt(eta) - 0.0625d0*pi_*pi_*eta**(-1.5d0))/g52 )
       end do
       return
    end if
    lo = 0; hi = f_neta
    do while( hi-lo > 1 )                            ! binary search: f_tab(lo,2) <= y < f_tab(hi,2)
       mid = (lo+hi)/2
       if( f_tab(mid,2) <= y )then; lo=mid; else; hi=mid; end if
    end do
    w = ( y - f_tab(lo,2) )/( f_tab(lo+1,2) - f_tab(lo,2) )
    eta = f_eta0 + ( dble(lo) + w )*f_deta
  end function ttm3_chem_pot

  ! Fermi-Dirac carrier heat capacity (per volume, atomic units, kB=1):
  ! Ce/N = (15/4) F_{3/2}/F_{1/2} - (9/4) F_{1/2}/F_{-1/2}
  ! (-> 3/2 in the non-degenerate limit, recovering the classical value).
  pure function ttm3_heat_capacity( T, mc, N ) result( Ce )
    implicit none
    real(8),intent(in) :: T,mc,N
    real(8) :: Ce,eta,Fm,F0,Fp
    eta = ttm3_chem_pot(T,mc,N)
    Fm = ttm3_fermi(1,eta); F0 = ttm3_fermi(2,eta); Fp = ttm3_fermi(3,eta)
    Ce = N*( 3.75d0*Fp/F0 - 2.25d0*F0/Fm )
  end function ttm3_heat_capacity

  ! Set A: Fermi-Dirac carrier thermal conductivity (Wiedemann-Franz with degeneracy
  ! corrections), reference Coef_Nc NTc(3)/(10).  K = (6 F_2/F_0 - 4 (F_1/F_0)^2) N Te mu,
  ! mu = mob0 F_0/F_{1/2}.  (T is in energy units, so the reference's extra kbT is absorbed.)
  pure function ttm3_thermal_cond( T, mc, N, mob0 ) result( K )
    implicit none
    real(8),intent(in) :: T,mc,N,mob0
    real(8) :: K,eta,F0,F12,F1,F2,mu
    eta = ttm3_chem_pot(T,mc,N)
    F0  = ttm3_F0(eta)                       ! F_0 = ln(1+exp(eta))
    F12 = ttm3_fermi(2,eta)                  ! F_{1/2}
    F1  = ttm3_fermi(4,eta)                  ! F_1
    F2  = ttm3_fermi(5,eta)                  ! F_2
    mu  = mob0*F0/F12
    K   = ( 6.0d0*F2/F0 - 4.0d0*(F1/F0)**2 )*N*T*mu
  end function ttm3_thermal_cond

  !---------------------------------------------------------------------------
  ! Local rates for one cell (the right-hand side of the five ODEs).
  pure subroutine ttm3_rates( Te_,Th_,Tl_,Ne_,Nh_, source_, gen_, &
                              dTe,dTh,dTl,dNe,dNh )
    implicit none
    real(8),intent(in)  :: Te_,Th_,Tl_,Ne_,Nh_, source_, gen_
    real(8),intent(out) :: dTe,dTh,dTl,dNe,dNh
    real(8) :: Ce,Ch,R,inv_tau,geff,q_heat,red_e,red_h

    ! Fermi-Dirac carrier heat capacities (reduce to the classical 3/2 N when
    ! non-degenerate; suppressed when degenerate)
    Ce = ttm3_heat_capacity( Te_, tp3%mu_e, Ne_+N_floor )
    Ch = ttm3_heat_capacity( Th_, tp3%mu_h, Nh_+N_floor )
    inv_tau = 1.0d0/tp3%tau

    ! generation, saturated at the reference/atomic density N0 (cannot exceed it).
    ! Single bleach: when sig_cold>0 the band-filling factor is already applied once
    ! inside ttm3_linear_gen (same convention as ttm3_step_cell).
    if( tp3%sig_cold > 0.0d0 )then
       geff = gen_
    else
       geff = gen_ * max( 0.0d0, 1.0d0 - Ne_/tp3%N0 )
    end if

    ! standard Auger recombination (removes electron-hole pairs).  This is the
    ! unsaturated limit of the reference R = N*A*N*Ne*Nh/(1+T_au*A*Ne*Nh): at the
    ! carrier densities reached here (T_au*A*Ne*Nh << 1) the two coincide.
    R = ( tp3%A_e*Ne_ + tp3%A_h*Nh_ )*Ne_*Nh_

    ! carrier densities: generation - recombination (charge-neutral)
    dNe = geff - R
    dNh = geff - R

    ! Reference energy partition between the electron and hole subsystems:
    ! red_e = mu_h/(mu_e+mu_h), red_h = mu_e/(mu_e+mu_h) (lighter carrier takes
    ! the larger share).  For Si this is 0.692/0.308, setting the Te:Th ratio.
    red_e = tp3%mu_h/(tp3%mu_e+tp3%mu_h)
    red_h = tp3%mu_e/(tp3%mu_e+tp3%mu_h)

    ! only the absorbed power in excess of the gap cost (Egap per generated pair)
    ! becomes carrier kinetic energy (heating); the rest creates the pair
    q_heat = max( source_ - geff*tp3%Egap, 0.0d0 )

    ! carrier temperatures: relaxation toward the lattice + partitioned heating +
    ! dilution by freshly generated carriers (born near the lattice temperature),
    ! which caps Te near the hot-carrier value (~(hw-Egap)) instead of diverging.
    dTe = -(Te_-Tl_)*inv_tau + (geff/(Ne_+N_floor))*(Tl_-Te_) + red_e*q_heat/Ce
    dTh = -(Th_-Tl_)*inv_tau + (geff/(Nh_+N_floor))*(Tl_-Th_) + red_h*q_heat/Ch

    ! lattice receives the energy lost by the carriers (energy-conserving coupling)
    dTl = ( Ce*(Te_-Tl_) + Ch*(Th_-Tl_) )*inv_tau / tp3%Cl
  end subroutine ttm3_rates

  !---------------------------------------------------------------------------
  ! Advance a single cell by one step: an exponential integrator for the (stiff)
  ! carrier temperatures and explicit Euler for the (non-stiff) densities/lattice.
  ! The carrier-temperature ODE is a relaxation toward the diluted+heated target
  !   dTe/dt = -rate_e*(Te - Te_star),   rate_e = 1/tau + geff/Ne,
  !   Te_star = Tl + (red_e*q_heat/Ce)/rate_e
  ! whose rate geff/Ne (dilution by cold, freshly-generated carriers) can far
  ! exceed 1/dt at low Ne.  The exact update Te <- Te_star+(Te-Te_star)*exp(-rate_e*dt)
  ! is unconditionally stable there (RK4 was not -- it blew up once Ne was small),
  ! and since Te_star -> Tl + (2/3)red_e*(hw-Egap) as geff/Ne dominates, the hot-
  ! carrier fixed point is reached without a clamp.  This removes the stiffness
  ! limit on the carrier-density floor, so Ne can reach the reference's near-
  ! transparent regime.
  subroutine ttm3_step_cell( Te_,Th_,Tl_,Ne_,Nh_, source_, gen_, dt_in )
    implicit none
    real(8),intent(inout) :: Te_,Th_,Tl_,Ne_,Nh_
    real(8),intent(in)    :: source_, gen_, dt_in
    real(8) :: Ce,Ch,R,geff,q_heat,red_e,red_h,inv_tau
    real(8) :: rate_e,rate_h,Te_star,Th_star,dTl
    real(8) :: Cef_e,Cef_h,Re,Rh,gap,TlK,Cl_eff,Col

    Ce = ttm3_heat_capacity( Te_, tp3%mu_e, Ne_+N_floor )
    Ch = ttm3_heat_capacity( Th_, tp3%mu_h, Nh_+N_floor )
    inv_tau = (1.0d0/tp3%tau)/( 1.0d0 + (Ne_*8443.9d0)**2 )      ! density-dependent e-ph relaxation
    if( tp3%sig_cold > 0.0d0 )then
       ! band filling is already applied once inside ttm3_linear_gen (sig_b = sig_cold*bf),
       ! matching the reference's single bleach; a second factor here would give bf^2.
       geff = gen_
    else
       geff = gen_ * max( 0.0d0, 1.0d0 - Ne_/tp3%N0 )
    end if
    Cef_e = tp3%A_e*Ne_*Nh_ ; Re = Ne_*Cef_e/( 1.0d0 + T_au*Cef_e )   ! saturated Auger (e)
    Cef_h = tp3%A_h*Ne_*Nh_ ; Rh = Nh_*Cef_h/( 1.0d0 + T_au*Cef_h )   ! saturated Auger (h)
    R    = ( Re + Rh ) * max( 0.0d0, 1.0d0 - Ne_/tp3%N0 )   ! band-filling (matches standalone Reh*Rate_eh)
    red_e = tp3%mu_h/(tp3%mu_e+tp3%mu_h)
    red_h = tp3%mu_e/(tp3%mu_e+tp3%mu_h)
    gap = ttm3_gap( Ne_, Tl_ )                                   ! carrier + thermal (Varshni) gap
    q_heat = max( source_ - geff*gap, 0.0d0 )
    ! thermal across-gap (collisional) carrier generation, Pauli-blocked
    Col = 3.6d-5*2.419d-2*0.5d0*( Ne_*exp(-1.5d0*gap/Te_) + Nh_*exp(-1.5d0*gap/Th_) ) &
          * max( 0.0d0, 1.0d0 - Ne_/tp3%N0 )

    ! lattice: relaxation heating from the carriers (non-stiff, Cl large), with the
    ! temperature-dependent lattice heat capacity Cl(Tl) (= input value at 300 K).
    ! Floored: the polynomial goes negative below ~1.4 K (and -inf at 0), which would
    ! flip the carrier->lattice coupling sign.
    TlK    = max(Tl_,T_floor)*hk_kelvin
    Cl_eff = tp3%Cl*max( 1.978d0 + 3.54d-4*TlK - 3.68d0/(TlK*TlK), 0.2d0 )/2.084d0
    dTl = ( Ce*(Te_-Tl_) + Ch*(Th_-Tl_) )*inv_tau / Cl_eff

    ! carrier temperatures: exponential integrator (stiff relaxation + dilution).
    ! Carrier-number energy bookkeeping, matching the reference Tnt = -(dN/dt)*Cnteh
    ! with Cnteh = (thermal energy per carrier) + red*gap:
    !  - ALL freshly generated carriers (optical geff AND collisional Col) are born
    !    near Tl and dilute the carrier temperature -> dilution rate (geff+Col)/N;
    !  - impact ionization (Col) draws the carriers' share red*gap from the carrier
    !    kinetic energy to pay for each new pair (an energy SINK).  Col ~ exp(-gap/Te)
    !    grows with Te, so this is a negative feedback that self-limits Te -- the
    !    physical cap, vs the unphysical free-carrier heating runaway when it is left
    !    out (Col was previously fed into dN/dt only, not the energy balance);
    !  - Auger recombination (R) returns red*gap to the carriers (recombination heating).
    rate_e = inv_tau + (geff+Col)/(Ne_+N_floor)
    rate_h = inv_tau + (geff+Col)/(Nh_+N_floor)
    Te_star = Tl_ + ( red_e*( q_heat + gap*(R-Col) )/Ce )/rate_e
    Th_star = Tl_ + ( red_h*( q_heat + gap*(R-Col) )/Ch )/rate_h
    Te_ = Te_star + (Te_-Te_star)*exp( -rate_e*dt_in )
    Th_ = Th_star + (Th_-Th_star)*exp( -rate_h*dt_in )

    ! carrier densities and lattice: explicit Euler (non-stiff); +Col thermal generation
    Ne_ = min( max( Ne_ + dt_in*(geff - R + Col), 0.0d0 ), tp3%N0 )   ! cap at valence density
    Nh_ = min( max( Nh_ + dt_in*(geff - R + Col), 0.0d0 ), tp3%N0 )
    Tl_ = max( Tl_ + dt_in*dTl, T_floor )

    ! numerical guard (should not trigger now that the stiff modes are integrated exactly)
    Te_ = min(max(Te_,T_floor),T_clamp); Th_ = min(max(Th_,T_floor),T_clamp)
  end subroutine ttm3_step_cell

  !---------------------------------------------------------------------------
  ! Advance every medium cell by one time step (transport OFF).
  subroutine ttm3_main( srg, rg, source, gen )
    use structures, only: s_rgrid, s_sendrecv_grid
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg     ! used by ttm3_transport (Stage 4)
    type(s_rgrid),         intent(in)    :: rg
    real(8),intent(in) :: source(rg%is(1):,rg%is(2):,rg%is(3):)
    real(8),intent(in) :: gen   (rg%is(1):,rg%is(2):,rg%is(3):)
    integer :: m,ix,iy,iz

!$omp parallel do private(m,ix,iy,iz)
    do m=1,nmedia_myrnk
       ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
       call ttm3_step_cell( Te(ix,iy,iz),Th(ix,iy,iz),Tl(ix,iy,iz), &
                            Ne(ix,iy,iz),Nh(ix,iy,iz), &
                            source(ix,iy,iz), gen(ix,iy,iz), dt )
    end do
!$omp end parallel do

    ! Stage 4: spatial transport (operator-split explicit diffusion + heat conduction)
    call ttm3_transport( srg, rg )
  end subroutine ttm3_main

  !---------------------------------------------------------------------------
  ! Stage 4: spatial transport (operator-split from the local step).  Two paths:
  !  * constant-coefficient path (Ddiff/kappa inputs, mobility=0): the original
  !    capped explicit Laplacians -- unchanged, regression-stable;
  !  * Set A Fermi path (mobility>0): conservative FLUX-FORM finite volumes with
  !    harmonic interface diffusivities, donor-cell (upwind) cross advection and
  !    drift, compensated convective carrier heat, and adaptive sub-stepping.
  !    The former point-local splitting (K_loc*lapT + gradK.gradT, D_loc*divJ +
  !    gradD.J, Cj*(...)) net-created energy at carrier-density shoulders where
  !    N jumps 20-400x per cell: neighbour-sized fluxes divided by the LOCAL
  !    Ce ~ N_loc amplified the update by N_nb/N_loc (10^2..10^3) -- a
  !    dt-INDEPENDENT blow-up (s2 smoke forensics, 2026-07-06).
  subroutine ttm3_transport( srg, rg )
    use structures,    only: s_rgrid, s_sendrecv_grid
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_rgrid),         intent(in)    :: rg
    if( tp3%mob_e0>0.0d0 .or. tp3%mob_h0>0.0d0 )then
       call ttm3_transport_fermi( srg, rg )
    else if( tp3%Ddiff>0.0d0 .or. tp3%kappa_e>0.0d0 .or. tp3%kappa_l>0.0d0 )then
       call ttm3_transport_const( srg, rg )
    end if
  end subroutine ttm3_transport

  ! harmonic interface mean: bounded by twice the SMALLER argument, so a huge
  ! neighbour coefficient can never drive a small cell (maximum-principle helper).
  pure function harm( a, b ) result( c )
    implicit none
    real(8),intent(in) :: a,b
    real(8) :: c
    c = 2.0d0*a*b/max( a+b, 1.0d-300 )
  end function harm

  !---------------------------------------------------------------------------
  ! Constant-coefficient explicit diffusion (the pre-Set-A path, unchanged).
  subroutine ttm3_transport_const( srg, rg )
    use structures,    only: s_rgrid, s_sendrecv_grid
    use sendrecv_grid, only: update_overlap_real8
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_rgrid),         intent(in)    :: rg
    integer :: m,ix,iy,iz
    real(8) :: cx,cy,cz,Ce,Ch,lNe,lNh,lTe,lTh,lTl,Dmax,Dd
    real(8) :: TlK_t,Cl_t,dcap

    cx=1.0d0/hgs(1)**2; cy=1.0d0/hgs(2)**2; cz=1.0d0/hgs(3)**2
    ! explicit FTCS diffusion stability: cap every diffusivity so dt*D*2*(cx+cy+cz)<=0.5.
    Dmax = 0.5d0/( 2.0d0*dt*(cx+cy+cz) )
    Dd   = min(tp3%Ddiff, Dmax)

    call update_overlap_real8(srg, rg, Ne); call update_overlap_real8(srg, rg, Nh)
    call update_overlap_real8(srg, rg, Te); call update_overlap_real8(srg, rg, Th)
    call update_overlap_real8(srg, rg, Tl)
    ! reflecting (zero-flux) slab boundary: mirror the state into the vacuum face
    ! neighbours (after the halo exchange, so the exchange cannot overwrite it).
    call ttm3_mirror_state( .false. )

    dcap = 0.0d0
!$omp parallel do private(m,ix,iy,iz,Ce,Ch,lNe,lNh,lTe,lTh,lTl,TlK_t,Cl_t) reduction(max:dcap)
    do m=1,nmedia_myrnk
       ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
       lTe=(Te(ix+1,iy,iz)-2*Te(ix,iy,iz)+Te(ix-1,iy,iz))*cx+(Te(ix,iy+1,iz)-2*Te(ix,iy,iz)+Te(ix,iy-1,iz))*cy+(Te(ix,iy,iz+1)-2*Te(ix,iy,iz)+Te(ix,iy,iz-1))*cz
       lTh=(Th(ix+1,iy,iz)-2*Th(ix,iy,iz)+Th(ix-1,iy,iz))*cx+(Th(ix,iy+1,iz)-2*Th(ix,iy,iz)+Th(ix,iy-1,iz))*cy+(Th(ix,iy,iz+1)-2*Th(ix,iy,iz)+Th(ix,iy,iz-1))*cz
       lTl=(Tl(ix+1,iy,iz)-2*Tl(ix,iy,iz)+Tl(ix-1,iy,iz))*cx+(Tl(ix,iy+1,iz)-2*Tl(ix,iy,iz)+Tl(ix,iy-1,iz))*cy+(Tl(ix,iy,iz+1)-2*Tl(ix,iy,iz)+Tl(ix,iy,iz-1))*cz
       Ce=ttm3_heat_capacity(Te(ix,iy,iz),tp3%mu_e,Ne(ix,iy,iz)+N_floor)
       Ch=ttm3_heat_capacity(Th(ix,iy,iz),tp3%mu_h,Nh(ix,iy,iz)+N_floor)
       lNe=(Ne(ix+1,iy,iz)-2*Ne(ix,iy,iz)+Ne(ix-1,iy,iz))*cx+(Ne(ix,iy+1,iz)-2*Ne(ix,iy,iz)+Ne(ix,iy-1,iz))*cy+(Ne(ix,iy,iz+1)-2*Ne(ix,iy,iz)+Ne(ix,iy,iz-1))*cz
       lNh=(Nh(ix+1,iy,iz)-2*Nh(ix,iy,iz)+Nh(ix-1,iy,iz))*cx+(Nh(ix,iy+1,iz)-2*Nh(ix,iy,iz)+Nh(ix,iy-1,iz))*cy+(Nh(ix,iy,iz+1)-2*Nh(ix,iy,iz)+Nh(ix,iy,iz-1))*cz
       rhs_ne(ix,iy,iz)=dt*Dd*lNe
       rhs_nh(ix,iy,iz)=dt*Dd*lNh
       rhs_te(ix,iy,iz)=dt*min(tp3%kappa_e/Ce,Dmax)*lTe
       rhs_th(ix,iy,iz)=dt*min(tp3%kappa_e/Ch,Dmax)*lTh
       ! lattice conduction Kl(Tl)/Cl(Tl): both temperature-dependent, Tl floored before
       ! the power/polynomial (a zero/negative Tl would give inf/NaN before any commit floor)
       TlK_t = max(Tl(ix,iy,iz),T_floor)*hk_kelvin
       Cl_t  = tp3%Cl*max( 1.978d0 + 3.54d-4*TlK_t - 3.68d0/(TlK_t*TlK_t), 0.2d0 )/2.084d0
       dcap  = max( dcap, tp3%kappa_l*(TlK_t/300.0d0)**(-1.23d0)/Cl_t/Dmax )
       rhs_tl(ix,iy,iz)=dt*min( tp3%kappa_l*(TlK_t/300.0d0)**(-1.23d0)/Cl_t, Dmax )*lTl  ! Kl(Tl)
    end do
!$omp end parallel do

    ! warn (rank-local, log-cadence) when the CFL cap is actually rewriting the physics
    if( dcap > 1.0d0 .and. dcap > 1.999d0*dcap_worst )then
       write(*,'(A,ES10.2,A)') " ttm3 WARNING: transport diffusivity CFL-capped (worst D/Dmax =", dcap, &
                               " ) -- conduction/diffusion locally reduced."
    end if
    dcap_worst = max( dcap_worst, dcap )

!$omp parallel do private(m,ix,iy,iz)
    do m=1,nmedia_myrnk
       ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
       Ne(ix,iy,iz)=max(Ne(ix,iy,iz)+rhs_ne(ix,iy,iz),0.0d0)
       Nh(ix,iy,iz)=max(Nh(ix,iy,iz)+rhs_nh(ix,iy,iz),0.0d0)
       ! positivity floors: an explicit-step undershoot to T<=0 would flip the heat
       ! equation backward and blow up exp(-1.5*gap/T) in the next local step
       Te(ix,iy,iz)=max(Te(ix,iy,iz)+rhs_te(ix,iy,iz),T_floor)
       Th(ix,iy,iz)=max(Th(ix,iy,iz)+rhs_th(ix,iy,iz),T_floor)
       Tl(ix,iy,iz)=max(Tl(ix,iy,iz)+rhs_tl(ix,iy,iz),T_floor)
    end do
!$omp end parallel do
  end subroutine ttm3_transport_const

  !---------------------------------------------------------------------------
  ! Set A Fermi transport: conservative flux-form finite volumes + adaptive
  ! sub-stepping.  Per FDTD step: (1) freeze the per-cell coefficient fields
  ! (Kc, Dc, Mc, Cf=Ce, Cj, bE, bT, rgap, space-charge Efld/Efc; the lattice
  ! Kl(Tl)/Cl(Tl) are deliberately LIVE per sub-step -- cheap state functions,
  ! better for the nonlinear lattice diffusion), (2) advance the
  ! state by sub-steps dt_s chosen from an ALL-TERM stability estimate (diffusive
  ! interface CFL + advective CFL + a relative-change guard, re-evaluated every
  ! sub-step; rates are linear in dt_s so no retry pass is needed).
  ! Faces to non-medium cells carry exactly zero flux (reflecting slab boundary),
  ! so no vacuum value is ever read.  Fail-loud: the sub-step count is bounded;
  ! exceeding it stops the run (never freeze-and-continue).
  subroutine ttm3_transport_fermi( srg, rg )
    use structures,    only: s_rgrid, s_sendrecv_grid
    use sendrecv_grid, only: update_overlap_real8
    use communication, only: comm_get_min
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_rgrid),         intent(in)    :: rg
    integer,parameter :: nit_max = 1000
    real(8),parameter :: cfl_saf = 0.5d0, rel_cap = 0.1d0
    integer :: m,ix,iy,iz,jx,jy,jz,d,sg,nit,nktf
    real(8) :: kT,eta,F0,Fm,F12,F1,F2
    real(8) :: h,hi,TlK_t,Cl_i,Kl_i,Kl_j,Ce,Ch
    real(8) :: khe,khh,khl,dhe,dhh,ue,uh,gte,gth,gtl,gne,gnh,ggp,efc_f
    real(8) :: phn_e,phn_h,cjfe,cjfh,rr
    real(8) :: rte,rth,rtl,rne,rnh,rwe,rwh,rmax,relmax,tau,dts,tsc

    ! ---- per-FDTD-step coefficient pass (frozen over the sub-steps) ----------
!$omp parallel do private(m,ix,iy,iz)
    do m=1,nmedia_myrnk
       ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
       rgap(ix,iy,iz)=ttm3_gap( Ne(ix,iy,iz), Tl(ix,iy,iz) )    ! carrier + thermal renormalised gap
    end do
!$omp end parallel do
    call update_overlap_real8(srg, rg, rgap)

    nktf = 0
!$omp parallel do private(m,ix,iy,iz,kT,eta,F0,Fm,F12,F1,F2) reduction(+:nktf)
    do m=1,nmedia_myrnk
       ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
       ! electrons.  The coefficient kT is floored at a PHYSICAL scale: the state floor
       ! (~1 K) would let the 1/kT thermodiffusion coefficient amplify ~300x (s2 crash).
       kT = max( Te(ix,iy,iz), 0.5d0*max(Tl(ix,iy,iz),tp3%Tini) )
       if( kT > Te(ix,iy,iz) ) nktf = nktf + 1
       eta=ttm3_chem_pot(kT,tp3%mu_e,Ne(ix,iy,iz)+N_floor)
       F0=ttm3_F0(eta); Fm=ttm3_fermi(1,eta); F12=ttm3_fermi(2,eta)
       F1=ttm3_fermi(4,eta); F2=ttm3_fermi(5,eta)
       Mc_e(ix,iy,iz)=tp3%mob_e0*F0/F12
       Dc_e(ix,iy,iz)=tp3%mob_e0*kT*F0/Fm
       Kc_e(ix,iy,iz)=(6.0d0*F2/F0-4.0d0*(F1/F0)**2)*Ne(ix,iy,iz)*kT*Mc_e(ix,iy,iz)
       ! thermodynamic divisor (Ce*dT = dW - Cj*dN): LIVE temperature, not the transport
       ! kT floor -- flooring the EOS would misconvert energy whenever the floor engages
       Cf_e(ix,iy,iz)=ttm3_heat_capacity(max(Te(ix,iy,iz),T_floor),tp3%mu_e,Ne(ix,iy,iz)+N_floor)
       ! N-normalized cross-advection velocities per unit gradient (bounded ~Dc/kT;
       ! the N_neighbour/N_local amplification cancels exactly):
       !   bE = Dc*cjE/(Ne+floor) ;  bT = Dc*cjT/(Ne+floor)
       bE_e(ix,iy,iz)=Dc_e(ix,iy,iz)*(tp3%mu_h/(tp3%mu_e+tp3%mu_h))*Fm/F12/kT
       bT_e(ix,iy,iz)=Dc_e(ix,iy,iz)*(2.0d0*Fm*F1/(F12*F0)-1.5d0)/kT
       Cj_e(ix,iy,iz)=2.0d0*kT*F1/F0 + (tp3%mu_h/(tp3%mu_e+tp3%mu_h))*rgap(ix,iy,iz)   ! heat per carrier
       ! holes
       kT = max( Th(ix,iy,iz), 0.5d0*max(Tl(ix,iy,iz),tp3%Tini) )
       if( kT > Th(ix,iy,iz) ) nktf = nktf + 1
       eta=ttm3_chem_pot(kT,tp3%mu_h,Nh(ix,iy,iz)+N_floor)
       F0=ttm3_F0(eta); Fm=ttm3_fermi(1,eta); F12=ttm3_fermi(2,eta)
       F1=ttm3_fermi(4,eta); F2=ttm3_fermi(5,eta)
       Mc_h(ix,iy,iz)=tp3%mob_h0*F0/F12
       Dc_h(ix,iy,iz)=tp3%mob_h0*kT*F0/Fm
       Kc_h(ix,iy,iz)=(6.0d0*F2/F0-4.0d0*(F1/F0)**2)*Nh(ix,iy,iz)*kT*Mc_h(ix,iy,iz)
       Cf_h(ix,iy,iz)=ttm3_heat_capacity(max(Th(ix,iy,iz),T_floor),tp3%mu_h,Nh(ix,iy,iz)+N_floor)
       bE_h(ix,iy,iz)=Dc_h(ix,iy,iz)*(tp3%mu_e/(tp3%mu_e+tp3%mu_h))*Fm/F12/kT
       bT_h(ix,iy,iz)=Dc_h(ix,iy,iz)*(2.0d0*Fm*F1/(F12*F0)-1.5d0)/kT
       Cj_h(ix,iy,iz)=2.0d0*kT*F1/F0 + (tp3%mu_e/(tp3%mu_e+tp3%mu_h))*rgap(ix,iy,iz)   ! heat per carrier
    end do
!$omp end parallel do
    if( nktf > 0 .and. .not.ktw_warned )then
       write(*,'(A)') " ttm3 NOTE: transport-coefficient kT floored at 0.5*max(Tl,Tini) on this rank."
       ktw_warned = .true.
    end if

    call update_overlap_real8(srg, rg, Kc_e); call update_overlap_real8(srg, rg, Kc_h)
    call update_overlap_real8(srg, rg, Dc_e); call update_overlap_real8(srg, rg, Dc_h)
    call update_overlap_real8(srg, rg, Mc_e); call update_overlap_real8(srg, rg, Mc_h)
    call update_overlap_real8(srg, rg, Cf_e); call update_overlap_real8(srg, rg, Cf_h)
    call update_overlap_real8(srg, rg, bE_e); call update_overlap_real8(srg, rg, bT_e)
    call update_overlap_real8(srg, rg, bE_h); call update_overlap_real8(srg, rg, bT_h)
    call update_overlap_real8(srg, rg, Cj_e); call update_overlap_real8(srg, rg, Cj_h)
    call ttm3_space_charge()   ! Efld (cell centres) + Efc (x-faces, exact prefix); frozen over sub-steps
    ! make rank-boundary faces bit-identical on both owners: the halo copy of Efc
    ! replaces the locally-recomputed offset value with the owning rank's exact one
    call update_overlap_real8(srg, rg, Efc)

    call update_overlap_real8(srg, rg, Ne); call update_overlap_real8(srg, rg, Nh)
    call update_overlap_real8(srg, rg, Te); call update_overlap_real8(srg, rg, Th)
    call update_overlap_real8(srg, rg, Tl)

    ! ---- adaptive sub-step loop ----------------------------------------------
    tau = dt ; nit = 0
    do while( tau > 0.0d0 )
       nit = nit + 1
       if( nit > nit_max )then
          write(*,'(A,I8,A)') " ttm3 ERROR: transport sub-stepping exceeded", nit_max, &
                              " iterations in one FDTD step (unresolvable stiffness)."
          error stop 'ttm3_transport_fermi: sub-step limit'
       end if

       ! one flux pass: per-cell RATES (per unit time) + stability measures.
       ! Gather form: each cell evaluates its own faces; the face formulas are
       ! symmetric under (i,j) swap, so both owners compute bit-identical fluxes
       ! (conservative to machine precision, incl. across MPI ranks via halos).
       rmax = 0.0d0 ; relmax = 0.0d0
!$omp parallel do private(m,ix,iy,iz,jx,jy,jz,d,sg,h,hi,kT,TlK_t,Cl_i,Kl_i,Kl_j,Ce,Ch, &
!$omp&                    khe,khh,khl,dhe,dhh,ue,uh,gte,gth,gtl,gne,gnh,ggp,efc_f, &
!$omp&                    phn_e,phn_h,cjfe,cjfh,rr,rte,rth,rtl,rne,rnh,rwe,rwh,tsc) &
!$omp&         reduction(max:rmax,relmax)
       do m=1,nmedia_myrnk
          ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
          Ce = Cf_e(ix,iy,iz) ; Ch = Cf_h(ix,iy,iz)
          TlK_t = max(Tl(ix,iy,iz),T_floor)*hk_kelvin
          Cl_i  = tp3%Cl*max( 1.978d0 + 3.54d-4*TlK_t - 3.68d0/(TlK_t*TlK_t), 0.2d0 )/2.084d0
          Kl_i  = tp3%kappa_l*(TlK_t/300.0d0)**(-1.23d0)
          rte=0.0d0; rth=0.0d0; rtl=0.0d0; rne=0.0d0; rnh=0.0d0
          rwe=0.0d0; rwh=0.0d0; rr=0.0d0
          do d=1,3
             h = hgs(d) ; hi = 1.0d0/h
             do sg=-1,1,2
                jx=ix; jy=iy; jz=iz
                select case(d)
                case(1); jx=ix+sg
                case(2); jy=iy+sg
                case(3); jz=iz+sg
                end select
                if( .not.med_mask(jx,jy,jz) ) cycle      ! reflecting boundary: zero face flux
                ! +d-oriented face gradients (sg folds the orientation)
                gte=dble(sg)*(Te(jx,jy,jz)-Te(ix,iy,iz))*hi
                gth=dble(sg)*(Th(jx,jy,jz)-Th(ix,iy,iz))*hi
                gtl=dble(sg)*(Tl(jx,jy,jz)-Tl(ix,iy,iz))*hi
                gne=dble(sg)*(Ne(jx,jy,jz)-Ne(ix,iy,iz))*hi
                gnh=dble(sg)*(Nh(jx,jy,jz)-Nh(ix,iy,iz))*hi
                ggp=dble(sg)*(rgap(jx,jy,jz)-rgap(ix,iy,iz))*hi
                ! interface diffusivities: harmonic mean (bounded by the small side)
                khe = harm( Kc_e(ix,iy,iz), Kc_e(jx,jy,jz) )
                khh = harm( Kc_h(ix,iy,iz), Kc_h(jx,jy,jz) )
                TlK_t = max(Tl(jx,jy,jz),T_floor)*hk_kelvin
                Kl_j  = tp3%kappa_l*(TlK_t/300.0d0)**(-1.23d0)
                khl   = harm( Kl_i, Kl_j )
                dhe = harm( Dc_e(ix,iy,iz), Dc_e(jx,jy,jz) )
                dhh = harm( Dc_h(ix,iy,iz), Dc_h(jx,jy,jz) )
                ! cross advection (band-edge force + thermodiffusion) + x drift = face velocity
                ue = -0.5d0*(bE_e(ix,iy,iz)+bE_e(jx,jy,jz))*ggp - 0.5d0*(bT_e(ix,iy,iz)+bT_e(jx,jy,jz))*gte
                uh = -0.5d0*(bE_h(ix,iy,iz)+bE_h(jx,jy,jz))*ggp - 0.5d0*(bT_h(ix,iy,iz)+bT_h(jx,jy,jz))*gth
                if( d==1 )then
                   efc_f = Efc( min(ix,jx), iy, iz )     ! exact prefix face field (x only)
                   ue = ue - harm(Mc_e(ix,iy,iz),Mc_e(jx,jy,jz))*efc_f   ! electrons: v = -mu E
                   uh = uh + harm(Mc_h(ix,iy,iz),Mc_h(jx,jy,jz))*efc_f   ! holes:     v = +mu E
                end if
                ! +d-oriented carrier fluxes: Fickian diffusion + donor-cell advection
                if( sg==1 )then
                   phn_e = -dhe*gne + ue*merge( Ne(ix,iy,iz), Ne(jx,jy,jz), ue>0.0d0 )
                   phn_h = -dhh*gnh + uh*merge( Nh(ix,iy,iz), Nh(jx,jy,jz), uh>0.0d0 )
                   cjfe  = merge( Cj_e(ix,iy,iz), Cj_e(jx,jy,jz), phn_e>0.0d0 )
                   cjfh  = merge( Cj_h(ix,iy,iz), Cj_h(jx,jy,jz), phn_h>0.0d0 )
                else
                   phn_e = -dhe*gne + ue*merge( Ne(jx,jy,jz), Ne(ix,iy,iz), ue>0.0d0 )
                   phn_h = -dhh*gnh + uh*merge( Nh(jx,jy,jz), Nh(ix,iy,iz), uh>0.0d0 )
                   cjfe  = merge( Cj_e(jx,jy,jz), Cj_e(ix,iy,iz), phn_e>0.0d0 )
                   cjfh  = merge( Cj_h(jx,jy,jz), Cj_h(ix,iy,iz), phn_h>0.0d0 )
                end if
                ! cell-i contribution of a +d-oriented face flux is -sg*PHI/h
                rne = rne - dble(sg)*phn_e*hi
                rnh = rnh - dble(sg)*phn_h*hi
                rwe = rwe - dble(sg)*cjfe*phn_e*hi       ! advected carrier heat (donor Cj)
                rwh = rwh - dble(sg)*cjfh*phn_h*hi
                rte = rte + dble(sg)*khe*gte*hi          ! conduction: -sg*(-K dT)/h
                rth = rth + dble(sg)*khh*gth*hi
                rtl = rtl + dble(sg)*khl*gtl*hi
                ! all-term stability accumulator (diffusive + advective)
                rr  = rr + max( khe/Ce, khh/Ch, khl/Cl_i, dhe, dhh )*hi*hi &
                         + max( abs(ue), abs(uh) )*hi
             end do
          end do
          ! temperature rates: conduction + COMPENSATED convective heat.  Continuous
          ! identity div(Cj*F) - Cj*div(F) = F.grad(Cj): advection alone must not
          ! change T (exactly zero for uniform Cj since rwe and rne use the SAME fluxes).
          rhs_te(ix,iy,iz) = ( rte + rwe - Cj_e(ix,iy,iz)*rne )/Ce
          rhs_th(ix,iy,iz) = ( rth + rwh - Cj_h(ix,iy,iz)*rnh )/Ch
          rhs_tl(ix,iy,iz) = rtl/Cl_i
          rhs_ne(ix,iy,iz) = rne
          rhs_nh(ix,iy,iz) = rnh
          rmax = max( rmax, rr )
          ! relative-change guard.  Channels that sit AT their floor and keep draining
          ! must not throttle dt_s: the commit clamp makes the downward dynamics inert
          ! there, and holding dt_s hostage to an unresolvable decay deadlocks the
          ! sub-stepper at nit_max (s8b: front-face Th collapse, 2026-07-06).  Carrier
          ! scales are floored at a ppb of the saturation density for the same reason.
          if( .not.( Te(ix,iy,iz) <= 1.01d0*T_floor .and. rhs_te(ix,iy,iz) < 0.0d0 ) )then
             tsc = max( Te(ix,iy,iz), 0.05d0*tp3%Tini )
             relmax = max( relmax, abs(rhs_te(ix,iy,iz))/tsc )
          end if
          if( .not.( Th(ix,iy,iz) <= 1.01d0*T_floor .and. rhs_th(ix,iy,iz) < 0.0d0 ) )then
             tsc = max( Th(ix,iy,iz), 0.05d0*tp3%Tini )
             relmax = max( relmax, abs(rhs_th(ix,iy,iz))/tsc )
          end if
          if( .not.( Tl(ix,iy,iz) <= 1.01d0*T_floor .and. rhs_tl(ix,iy,iz) < 0.0d0 ) )then
             tsc = max( Tl(ix,iy,iz), 0.05d0*tp3%Tini )
             relmax = max( relmax, abs(rhs_tl(ix,iy,iz))/tsc )
          end if
          if( .not.( Ne(ix,iy,iz) <= 0.0d0 .and. rhs_ne(ix,iy,iz) < 0.0d0 ) ) &
             relmax = max( relmax, abs(rhs_ne(ix,iy,iz))/max(Ne(ix,iy,iz),N_floor,1.0d-9*tp3%N0) )
          if( .not.( Nh(ix,iy,iz) <= 0.0d0 .and. rhs_nh(ix,iy,iz) < 0.0d0 ) ) &
             relmax = max( relmax, abs(rhs_nh(ix,iy,iz))/max(Nh(ix,iy,iz),N_floor,1.0d-9*tp3%N0) )
       end do
!$omp end parallel do

       ! sub-step size: interface CFL + relative-change cap; rates are linear in
       ! dt_s so no retry pass is needed.  Global min keeps all ranks in lockstep
       ! (identical sub-step count => matched halo exchanges).
       dts = tau
       if( rmax   > 0.0d0 ) dts = min( dts, cfl_saf/rmax )
       if( relmax > 0.0d0 ) dts = min( dts, rel_cap/relmax )
       if( nprocs_g > 1 ) call comm_get_min( dts, comm )
       ! underflow guard: only when the STABILITY caps (not a small final remainder tau)
       ! forced the collapse
       if( dts < dt/dble(nit_max) .and. tau > dt/dble(nit_max) )then
          write(*,'(A,ES10.2)') " ttm3 ERROR: transport sub-step collapsed, dt_s/dt =", dts/dt
          error stop 'ttm3_transport_fermi: sub-step underflow'
       end if

!$omp parallel do private(m,ix,iy,iz)
       do m=1,nmedia_myrnk
          ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
          Ne(ix,iy,iz)=max(Ne(ix,iy,iz)+dts*rhs_ne(ix,iy,iz),0.0d0)
          Nh(ix,iy,iz)=max(Nh(ix,iy,iz)+dts*rhs_nh(ix,iy,iz),0.0d0)
          ! positivity floors: last-resort guards (the sub-step control should keep
          ! the update far from them; see the rel_cap bound above)
          Te(ix,iy,iz)=max(Te(ix,iy,iz)+dts*rhs_te(ix,iy,iz),T_floor)
          Th(ix,iy,iz)=max(Th(ix,iy,iz)+dts*rhs_th(ix,iy,iz),T_floor)
          Tl(ix,iy,iz)=max(Tl(ix,iy,iz)+dts*rhs_tl(ix,iy,iz),T_floor)
       end do
!$omp end parallel do

       tau = tau - dts
       if( tau > 0.0d0 )then
          if( tau < 1.0d-12*dt )then                    ! FP crumbs: fold into done
             tau = 0.0d0
          else                                          ! state changed: refresh halos
             call update_overlap_real8(srg, rg, Ne); call update_overlap_real8(srg, rg, Nh)
             call update_overlap_real8(srg, rg, Te); call update_overlap_real8(srg, rg, Th)
             call update_overlap_real8(srg, rg, Tl)
          end if
       end if
    end do

    ! rank-local, log-cadence sub-step report (visibility without spam)
    if( nit > nsub_worst )then
       write(*,'(A,I6,A)') " ttm3 NOTE: transport used", nit, " sub-steps in one FDTD step."
       nsub_worst = max( 2*nsub_worst, nit )
    end if
  end subroutine ttm3_transport_fermi

  !---------------------------------------------------------------------------
  ! Reflecting (zero-flux) closure at the medium surface: copy each surface cell's
  ! state into its vacuum face neighbour, so the central differences see a zero
  ! normal gradient.  Serial loops (nsurf is small; deterministic on duplicate
  ! targets, which only occur for non-slab corner geometries -- warned at init).
  subroutine ttm3_mirror_state( use_fermi_ )
    implicit none
    logical,intent(in) :: use_fermi_
    integer :: m,ix,iy,iz,id,jx,jy,jz
    do m=1,nsurf
       ix=ijk_surf(1,m); iy=ijk_surf(2,m); iz=ijk_surf(3,m); id=ijk_surf(4,m)
       jx=ix; jy=iy; jz=iz
       select case( abs(id) )
       case(1); jx = ix + sign(1,id)
       case(2); jy = iy + sign(1,id)
       case(3); jz = iz + sign(1,id)
       end select
       Ne(jx,jy,jz)=Ne(ix,iy,iz); Nh(jx,jy,jz)=Nh(ix,iy,iz)
       Te(jx,jy,jz)=Te(ix,iy,iz); Th(jx,jy,jz)=Th(ix,iy,iz); Tl(jx,jy,jz)=Tl(ix,iy,iz)
       if( use_fermi_ ) rgap(jx,jy,jz)=rgap(ix,iy,iz)
    end do
  end subroutine ttm3_mirror_state


  !---------------------------------------------------------------------------
  ! Space-charge field Efld along x: Ef(ix)=2*pi*inveps*[sum_{i<ix} - sum_{i>ix}] rho*dx.
  ! A sheet at x' with areal charge sigma = rho*dx contributes E = +-2*pi*sigma/eps, so the
  ! prefix weight is the x cell size hgs(1) ONLY (the old dV=hx*hy*hz overscaled Efld by the
  ! transverse cell area and made it dimensionally inconsistent with the local Gauss term).
  ! rho=(Nh-Ne) on MEDIUM cells only (vacuum cells hold mirrored values from the reflecting
  ! boundary and must not contribute).  Per (iy,iz) line, prefix sum over the inner x-range.
  subroutine ttm3_space_charge()
    use communication, only: comm_summation
    implicit none
    integer :: ix,iy,iz,p
    real(8) :: dx,c2,total,pre,rho
    real(8),allocatable :: ain(:,:,:), aout(:,:,:)
    dx = hgs(1)
    c2 = 2.0d0*pi_/tp3%eps_bg
    if( xc_isize <= 1 )then
       ! single x-domain (nproc_rgrid(1)=1): local exclusive prefix sum per (iy,iz) line
!$omp parallel do collapse(2) private(ix,iy,iz,total,pre,rho)
       do iz=is_inner(3),ie_inner(3)
       do iy=is_inner(2),ie_inner(2)
          total=0.0d0
          do ix=is_inner(1),ie_inner(1)
             if( med_mask(ix,iy,iz) ) total=total+(Nh(ix,iy,iz)-Ne(ix,iy,iz))
          end do
          pre=0.0d0                                  ! exclusive prefix sum (i<ix)
          Efc(is_inner(1)-1,iy,iz)=c2*dx*( -total )  ! face left of the first inner cell
          do ix=is_inner(1),ie_inner(1)
             rho=merge( Nh(ix,iy,iz)-Ne(ix,iy,iz), 0.0d0, med_mask(ix,iy,iz) )
             ! Ef = (left - right)*dx*c2 ;  right = total - pre - rho
             Efld(ix,iy,iz)=c2*dx*( 2.0d0*pre + rho - total )
             pre=pre+rho
             ! exact field at face ix+1/2: all charge through ix is now to the left
             Efc(ix,iy,iz)=c2*dx*( 2.0d0*pre - total )
          end do
       end do
       end do
!$omp end parallel do
    else
       ! x-decomposed (nproc_rgrid(1)>1): each x-rank deposits its per-(iy,iz)-line
       ! charge total into its own slot; an allreduce over icomm_x gathers all x-ranks
       ! (= allgather), giving the left-offset and global total per line for the prefix.
       allocate( ain (0:xc_isize-1, is_inner(2):ie_inner(2), is_inner(3):ie_inner(3)) )
       allocate( aout(0:xc_isize-1, is_inner(2):ie_inner(2), is_inner(3):ie_inner(3)) )
       ain = 0.0d0
       do iz=is_inner(3),ie_inner(3)
       do iy=is_inner(2),ie_inner(2)
          rho=0.0d0
          do ix=is_inner(1),ie_inner(1)
             if( med_mask(ix,iy,iz) ) rho=rho+(Nh(ix,iy,iz)-Ne(ix,iy,iz))
          end do
          ain(xc_id,iy,iz)=rho
       end do
       end do
       call comm_summation( ain, aout, size(ain), xc_icomm )
       do iz=is_inner(3),ie_inner(3)
       do iy=is_inner(2),ie_inner(2)
          total=0.0d0; pre=0.0d0
          do p=0,xc_isize-1
             total=total+aout(p,iy,iz)
             if( p<xc_id ) pre=pre+aout(p,iy,iz)   ! charge in x-domains left of this one
          end do
          Efc(is_inner(1)-1,iy,iz)=c2*dx*( 2.0d0*pre - total )   ! face into this x-domain
          do ix=is_inner(1),ie_inner(1)
             rho=merge( Nh(ix,iy,iz)-Ne(ix,iy,iz), 0.0d0, med_mask(ix,iy,iz) )
             Efld(ix,iy,iz)=c2*dx*( 2.0d0*pre + rho - total )
             pre=pre+rho
             Efc(ix,iy,iz)=c2*dx*( 2.0d0*pre - total )           ! exact face ix+1/2
          end do
       end do
       end do
       deallocate( ain, aout )
    end if
  end subroutine ttm3_space_charge

  !---------------------------------------------------------------------------
  ! Return the five fields at a probe cell, converted to output units.
  subroutine ttm3_get_state( a, Te_o, Th_o, Tl_o, Ne_o, Nh_o )
    implicit none
    integer,intent(in)  :: a(3)
    real(8),intent(out) :: Te_o, Th_o, Tl_o, Ne_o, Nh_o   ! [K], [K], [K], [1/cm^3], [1/cm^3]
    real(8), parameter :: hartree_kelvin_relationship = 3.1577502480407d5 ! [K]
    real(8), parameter :: atomic_unit_of_length = 5.29177210903d-11 ! [m]
    real(8) :: cm3_per_bohr3

    cm3_per_bohr3 = (atomic_unit_of_length*1.0d2)**3   ! (bohr in cm)^3
    Te_o = Te(a(1),a(2),a(3)) * hartree_kelvin_relationship
    Th_o = Th(a(1),a(2),a(3)) * hartree_kelvin_relationship
    Tl_o = Tl(a(1),a(2),a(3)) * hartree_kelvin_relationship
    Ne_o = Ne(a(1),a(2),a(3)) / cm3_per_bohr3
    Nh_o = Nh(a(1),a(2),a(3)) / cm3_per_bohr3
  end subroutine ttm3_get_state

  !---------------------------------------------------------------------------
  ! Peak values over this rank's medium cells (output units), for diagnostics.
  subroutine ttm3_get_max( Te_o, Th_o, Tl_o, Ne_o, Nh_o )
    use communication, only: comm_get_max
    implicit none
    real(8),intent(out) :: Te_o, Th_o, Tl_o, Ne_o, Nh_o   ! [K],[K],[K],[1/cm^3],[1/cm^3]
    real(8), parameter :: hartree_kelvin_relationship = 3.1577502480407d5
    real(8), parameter :: atomic_unit_of_length = 5.29177210903d-11
    real(8) :: cm3_per_bohr3, vin(5), vout(5)
    integer :: m,ix,iy,iz
    Te_o=0d0; Th_o=0d0; Tl_o=0d0; Ne_o=0d0; Nh_o=0d0
    cm3_per_bohr3 = (atomic_unit_of_length*1.0d2)**3
    do m=1,nmedia_myrnk
       ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
       Te_o=max(Te_o,Te(ix,iy,iz)); Th_o=max(Th_o,Th(ix,iy,iz)); Tl_o=max(Tl_o,Tl(ix,iy,iz))
       Ne_o=max(Ne_o,Ne(ix,iy,iz)); Nh_o=max(Nh_o,Nh(ix,iy,iz))
    end do
    Te_o=Te_o*hartree_kelvin_relationship; Th_o=Th_o*hartree_kelvin_relationship
    Tl_o=Tl_o*hartree_kelvin_relationship
    Ne_o=Ne_o/cm3_per_bohr3; Nh_o=Nh_o/cm3_per_bohr3
    if( nprocs_g>1 )then            ! global max across the spatial (r-space) decomposition
       vin(1)=Te_o; vin(2)=Th_o; vin(3)=Tl_o; vin(4)=Ne_o; vin(5)=Nh_o
       call comm_get_max( vin, vout, 5, comm )
       Te_o=vout(1); Th_o=vout(2); Tl_o=vout(3); Ne_o=vout(4); Nh_o=vout(5)
    end if
  end subroutine ttm3_get_max

  !---------------------------------------------------------------------------
  ! State at the illuminated front face: the medium cell on the GLOBAL minimum-ix
  ! plane nearest the GLOBAL transverse centroid -- the analogue of the reference's
  ! x=0 surface cell and, unlike ttm3_get_max, not biased by internal Fabry-Perot
  ! antinodes.  COLLECTIVE: must be called by every rank.  The winner cell is chosen
  ! by a decomposition-independent key that is unique per cell (centroid distance,
  ! then iy,iz), so ALL outputs -- state, env/power/gen, indices -- come from the
  ! same single cell for any MPI layout (the old per-rank centroid + elementwise max
  ! could mix fields from different cells and report a vacuum cell on rank 0).
  subroutine ttm3_get_front( env, pw, gn, Te_o, Th_o, Tl_o, Ne_o, Nh_o, &
                             env_o, pow_o, gen_o, gjx, gjy, gjz )
    use communication, only: comm_get_min, comm_get_max, comm_summation
    implicit none
    real(8),intent(in)  :: env(is_array(1):,is_array(2):,is_array(3):)  ! field envelope
    real(8),intent(in)  :: pw (is_inner(1):,is_inner(2):,is_inner(3):)  ! heating source
    real(8),intent(in)  :: gn (is_inner(1):,is_inner(2):,is_inner(3):)  ! generation rate
    real(8),intent(out) :: Te_o, Th_o, Tl_o, Ne_o, Nh_o   ! [K],[K],[K],[1/cm^3],[1/cm^3]
    real(8),intent(out) :: env_o, pow_o, gen_o            ! [a.u.] at the front cell
    integer,intent(out) :: gjx, gjy, gjz                  ! front cell (global grid indices)
    real(8), parameter :: hartree_kelvin_relationship = 3.1577502480407d5
    real(8), parameter :: atomic_unit_of_length = 5.29177210903d-11
    real(8) :: cm3_per_bohr3, cy, cz, gfix, dkey, dkeyg, cyz(3), cyzg(3), vin(11), vout(11)
    integer :: m,ix,iy,iz,ixmin,a(3),dbest,d
    logical :: have

    Te_o=tp3%Tini*hartree_kelvin_relationship; Th_o=Te_o; Tl_o=Te_o
    Ne_o=0d0; Nh_o=0d0; env_o=0d0; pow_o=0d0; gen_o=0d0
    gjx=is_inner(1); gjy=is_inner(2); gjz=is_inner(3)
    cm3_per_bohr3 = (atomic_unit_of_length*1.0d2)**3

    ! (1) global front plane = min ix over all medium cells
    gfix = dble(huge(1))
    do m=1,nmedia_myrnk
       gfix = min( gfix, dble(ijk_media_myrnk(1,m)) )
    end do
    if( nprocs_g>1 ) call comm_get_min( gfix, comm )
    if( gfix >= dble(huge(1)) ) return                  ! no medium cell anywhere
    ixmin = nint(gfix)

    ! (2) global transverse centroid of the front plane
    cyz = 0.0d0
    do m=1,nmedia_myrnk
       if( ijk_media_myrnk(1,m)==ixmin )then
          cyz(1)=cyz(1)+dble(ijk_media_myrnk(2,m)); cyz(2)=cyz(2)+dble(ijk_media_myrnk(3,m))
          cyz(3)=cyz(3)+1.0d0
       end if
    end do
    if( nprocs_g>1 )then
       call comm_summation( cyz, cyzg, 3, comm )
    else
       cyzg = cyz
    end if
    cy = cyzg(1)/max(cyzg(3),1.0d0); cz = cyzg(2)/max(cyzg(3),1.0d0)

    ! (3) local candidate on the plane: nearest the centroid, ties broken by (iy,iz)
    have=.false.; dbest=huge(0); a=(/ixmin,0,0/)
    do m=1,nmedia_myrnk
       if( ijk_media_myrnk(1,m)==ixmin )then
          iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
          d = (iy-nint(cy))**2 + (iz-nint(cz))**2
          if( .not.have .or. d<dbest .or. ( d==dbest .and. ( iy<a(2) .or. &
              ( iy==a(2) .and. iz<a(3) ) ) ) )then
             dbest=d; a(2)=iy; a(3)=iz; have=.true.
          end if
       end if
    end do

    ! (4) global winner by staged exact min-reductions on (distance, iy, iz):
    !     integer values carried in doubles are exact, so the winner cell is unique
    !     (one owner rank) for any grid size -- no packed-key collisions.
    dkey = merge( dble(dbest), huge(1.0d0), have )
    dkeyg = dkey
    if( nprocs_g>1 ) call comm_get_min( dkeyg, comm )
    if( dkeyg >= huge(1.0d0) ) return                   ! no candidate anywhere (paranoia)
    have = have .and. ( dkey == dkeyg )
    dkey = merge( dble(a(2)), huge(1.0d0), have )
    dkeyg = dkey
    if( nprocs_g>1 ) call comm_get_min( dkeyg, comm )
    have = have .and. ( dkey == dkeyg )
    dkey = merge( dble(a(3)), huge(1.0d0), have )
    dkeyg = dkey
    if( nprocs_g>1 ) call comm_get_min( dkeyg, comm )
    have = have .and. ( dkey == dkeyg )

    ! (5) the winner packs every output from its candidate cell; elementwise max
    !     is then a pure broadcast (all non-winners contribute -huge)
    if( have )then
       ix=a(1); iy=a(2); iz=a(3)
       vin(1)=Te(ix,iy,iz); vin(2)=Th(ix,iy,iz); vin(3)=Tl(ix,iy,iz)
       vin(4)=Ne(ix,iy,iz); vin(5)=Nh(ix,iy,iz)
       vin(6)=env(ix,iy,iz); vin(7)=pw(ix,iy,iz); vin(8)=gn(ix,iy,iz)
       vin(9)=dble(ix); vin(10)=dble(iy); vin(11)=dble(iz)
    else
       vin(:) = -huge(1.0d0)
    end if
    if( nprocs_g>1 )then
       call comm_get_max( vin, vout, 11, comm )
    else
       vout = vin
    end if
    Te_o=vout(1)*hartree_kelvin_relationship; Th_o=vout(2)*hartree_kelvin_relationship
    Tl_o=vout(3)*hartree_kelvin_relationship
    Ne_o=vout(4)/cm3_per_bohr3; Nh_o=vout(5)/cm3_per_bohr3
    env_o=vout(6); pow_o=vout(7); gen_o=vout(8)
    gjx=nint(vout(9)); gjy=nint(vout(10)); gjz=nint(vout(11))
  end subroutine ttm3_get_front

  ! Index of the illuminated front-face medium cell (same selection as ttm3_get_front).
  ! Diagnostic: lets the caller read the field envelope there to compare the surface
  ! intensity against the reference's out_all(6).
  subroutine ttm3_front_ijk( jx, jy, jz )
    implicit none
    integer,intent(out) :: jx,jy,jz
    integer :: m,iy,iz,ixmin,nmin,a(3),dbest,d
    real(8) :: cy,cz
    jx=is_inner(1); jy=is_inner(2); jz=is_inner(3)
    if( nmedia_myrnk<=0 ) return
    ixmin = ijk_media_myrnk(1,1)
    do m=2,nmedia_myrnk
       if( ijk_media_myrnk(1,m) < ixmin ) ixmin = ijk_media_myrnk(1,m)
    end do
    cy=0d0; cz=0d0; nmin=0
    do m=1,nmedia_myrnk
       if( ijk_media_myrnk(1,m)==ixmin )then
          cy=cy+ijk_media_myrnk(2,m); cz=cz+ijk_media_myrnk(3,m); nmin=nmin+1
       end if
    end do
    cy=cy/max(nmin,1); cz=cz/max(nmin,1)
    a = ijk_media_myrnk(1:3,1); dbest=huge(0)
    do m=1,nmedia_myrnk
       if( ijk_media_myrnk(1,m)==ixmin )then
          iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
          d = (iy-nint(cy))**2 + (iz-nint(cz))**2
          if( d<dbest )then; dbest=d; a=ijk_media_myrnk(1:3,m); end if
       end if
    end do
    jx=a(1); jy=a(2); jz=a(3)
  end subroutine ttm3_front_ijk

  ! Depth profile along x at the GLOBAL front-cell transverse line (gjy,gjz from
  ! ttm3_get_front): each rank fills its slice of the global x-line, a sum-gather
  ! assembles it, and root writes one block per call (time, ix, Te[K], Ne[cm^-3],
  ! Tl[K], |E|^2).  COLLECTIVE: must be called by every rank (the old version wrote
  ! only root's own x-slice = truncated profile under x-decomposition).
  subroutine ttm3_write_profile( fname, t_fs, u, field, gjy, gjz )
    use communication, only: comm_summation
    implicit none
    character(*),intent(in) :: fname
    real(8),     intent(in) :: t_fs
    integer,     intent(in) :: u, gjy, gjz
    real(8),     intent(in) :: field(is_array(1):,is_array(2):,is_array(3):)  ! |E|^2 envelope [a.u.]
    real(8),parameter :: hk = 3.1577502480407d5, aul = 5.29177210903d-11
    real(8) :: cm3
    integer :: m,ix,iy,iz
    real(8),allocatable :: buf(:,:), bufg(:,:)
    if( gx1 < gx0 ) return
    allocate( buf(gx0:gx1,0:4), bufg(gx0:gx1,0:4) ); buf=0.0d0
    cm3 = (aul*1.0d2)**3
    do m=1,nmedia_myrnk
       ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
       if( iy==gjy .and. iz==gjz )then                ! my cells on the global centre line
          buf(ix,0)=1.0d0                             ! occupancy (skips vacuum rows on write)
          buf(ix,1)=Te(ix,iy,iz)*hk
          buf(ix,2)=Ne(ix,iy,iz)/cm3
          buf(ix,3)=Tl(ix,iy,iz)*hk
          buf(ix,4)=field(ix,iy,iz)
       end if
    end do
    if( nprocs_g>1 )then
       call comm_summation( buf, bufg, size(buf), comm )
    else
       bufg = buf
    end if
    if( DISPLAY )then
       open(u, file=fname, status='unknown', position='append')
       do ix=gx0,gx1                                  ! sorted by depth
          if( bufg(ix,0) > 0.5d0 ) &
             write(u,'(F12.4,1X,I6,4(1X,E16.8))') t_fs, ix, bufg(ix,1), bufg(ix,2), bufg(ix,3), bufg(ix,4)
       end do
       write(u,*) ' '                                 ! blank line separates time blocks
       close(u)
    end if
    deallocate( buf, bufg )
  end subroutine ttm3_write_profile

  !---------------------------------------------------------------------------
  ! Stage 2: carrier-dependent permittivity (recommended Drude model).
  ! eps(omega) = eps_bg - omega_p^2/omega^2,  omega_p^2 = 4*pi*(Ne/m_e + Nh/m_h);
  ! the Drude damping (rate 1/tau) gives the conductivity sigma. All atomic units.
  ! Layer C-b1: strong carrier-temperature-dependent Drude (reference TOTAL_DIELE + Diel_Func).
  ! The electron-hole scattering rate gamma_eh rises with carrier density and temperature, so
  ! the free-carrier absorption climbs sharply as the material metallizes -- the T_e "second
  ! rise" that the fixed-tau Drude misses.  Ported verbatim with the reference constants
  ! (kbT, a_t, E_h, optical masses mu_eo/mu_ho); ttm3 stores T_e = k_B T [Ha] so the reference's
  ! Kelvin temperature is T/kbT.  Ne in bohr^-3 matches the reference.
  pure subroutine ttm3_drude( Ne_, Te_, Th_, Tl_, omega, eps_re, sig )
    implicit none
    real(8),intent(in)  :: Ne_, Te_, Th_, Tl_, omega
    real(8),intent(out) :: eps_re, sig
    real(8),parameter :: mu_eo=0.26d0, mu_ho=0.37d0                 ! Drude optical masses
    real(8),parameter :: kbT_r=3.1668d-6, a_t_=2.419d-2, E_h_=27.2114d0
    real(8) :: red,red2,Cg,Tryd,Pe,Tf,aa,gsp,gph,g1e,g1h,ge,gh,TeK,ThK,TlK,xE,xH,w2,wpe,wph,argl
    TeK=max(Te_,1.0d-12)/kbT_r; ThK=max(Th_,1.0d-12)/kbT_r; TlK=Tl_/kbT_r
    red  = 1.0d0/(1.0d0/mu_eo+1.0d0/mu_ho)
    red2 = 1.0d0/(1.0d0/tp3%mu_e+1.0d0/tp3%mu_h)                    ! band masses
    Cg   = 16.0d0/9.0d0*pi_**(-1.5d0)
    Tryd = red*0.5d0/(11.7d0**2)*(16.0d0*pi_*pi_)
    Pe   = (3.0d0*pi_*pi_*max(Ne_,1.0d-30))**(2.0d0/3.0d0)
    Tf   = 0.5d0*Pe/(red2*kbT_r)                                   ! Fermi temperature [K]
    argl = sqrt(Tryd*Tf)*Tf
    aa   = 1.0d0/sqrt( sqrt(kbT_r)/argl )*exp(2.0d0/3.0d0)
    if( TeK > aa )then                                             ! Spitzer e-h scattering
       gsp = abs( Cg*Tryd*(Tf/TeK)**1.5d0 * log( TeK*TeK*sqrt(kbT_r)/argl ) )
    else
       gsp = Cg*Tryd*(Tf/aa)**1.5d0 * log( aa*aa*sqrt(kbT_r)/argl )
    end if
    gph = 4.0d-4*TlK*a_t_                                          ! phonon
    g1e = a_t_/( 0.98d0 + 0.2d0*(kbT_r*TeK*E_h_)**(-3.5d0) )       ! correction (electron temp)
    g1h = a_t_/( 0.98d0 + 0.2d0*(kbT_r*ThK*E_h_)**(-3.5d0) )       ! correction (hole temp)
    ge  = max( gsp+gph+g1e, 1.0d-30 ); gh = max( gsp+gph+g1h, 1.0d-30 )
    w2  = omega*omega; xE = ge/omega; xH = gh/omega
    wpe = 4.0d0*pi_*Ne_/mu_eo; wph = 4.0d0*pi_*Ne_/mu_ho           ! ambipolar: Nh=Ne
    eps_re = -wpe/w2/(1.0d0+xE*xE) - wph/w2/(1.0d0+xH*xH)          ! Re(eps_Drude), e+h
    sig    = ( wpe/w2*xE/(1.0d0+xE*xE) + wph/w2*xH/(1.0d0+xH*xH) )*omega/(4.0d0*pi_)
  end subroutine ttm3_drude

  ! ADE-FDTD coefficients for the carrier Drude polarization current J (the STABLE,
  ! SALMON-native current-source coupling): J_new = c1*J + c2*E, with the current injected
  ! via eh_add_curr.  Same gamma_eh as ttm3_drude; combined e+h plasma frequency
  ! wp^2 = 4*pi*Ne*(1/mu_eo+1/mu_ho).  c1=(1-g dt/2)/(1+g dt/2), c2=wp^2 dt/(8 pi (1+g dt/2))
  ! (the same form as the validated Lorentz-Drude pole, c2_jx_ld).  This replaces the
  ! eps/sigma back-action (which is unstable in the real-field FDTD as sigma grows fast).
  pure subroutine ttm3_gamma_wp2( Ne_, Te_, Tl_, ge, wp2 )
    implicit none
    real(8),intent(in)  :: Ne_, Te_, Tl_
    real(8),intent(out) :: ge, wp2
    real(8),parameter :: mu_eo=0.26d0, mu_ho=0.37d0, kbT_r=3.1668d-6, a_t_=2.419d-2, E_h_=27.2114d0
    real(8) :: red,red2,Cg,Tryd,Pe,Tf,aa,gsp,gph,g1e,TeK,TlK,argl
    TeK=max(Te_,1.0d-12)/kbT_r; TlK=Tl_/kbT_r
    red  = 1.0d0/(1.0d0/mu_eo+1.0d0/mu_ho)
    red2 = 1.0d0/(1.0d0/tp3%mu_e+1.0d0/tp3%mu_h)
    Cg   = 16.0d0/9.0d0*pi_**(-1.5d0)
    Tryd = red*0.5d0/(11.7d0**2)*(16.0d0*pi_*pi_)
    Pe   = (3.0d0*pi_*pi_*max(Ne_,1.0d-30))**(2.0d0/3.0d0)
    Tf   = 0.5d0*Pe/(red2*kbT_r)
    argl = sqrt(Tryd*Tf)*Tf
    aa   = 1.0d0/sqrt( sqrt(kbT_r)/argl )*exp(2.0d0/3.0d0)
    if( TeK > aa )then
       gsp = abs( Cg*Tryd*(Tf/TeK)**1.5d0 * log( TeK*TeK*sqrt(kbT_r)/argl ) )
    else
       gsp = Cg*Tryd*(Tf/aa)**1.5d0 * log( aa*aa*sqrt(kbT_r)/argl )
    end if
    gph = 4.0d-4*TlK*a_t_
    g1e = a_t_/( 0.98d0 + 0.2d0*(kbT_r*TeK*E_h_)**(-3.5d0) )
    ge  = max( gsp+gph+g1e, 1.0d-30 )
    wp2 = 4.0d0*pi_*Ne_*( 1.0d0/mu_eo + 1.0d0/mu_ho )
  end subroutine ttm3_gamma_wp2

  pure subroutine ttm3_drude_coef( Ne_, Te_, Tl_, omega, dt, c1, c2 )
    implicit none
    real(8),intent(in)  :: Ne_, Te_, Tl_, omega, dt
    real(8),intent(out) :: c1, c2
    real(8) :: ge, wp2, c0
    call ttm3_gamma_wp2( Ne_, Te_, Tl_, ge, wp2 )
    if( tp3%coupling_mode == 2 )then              ! ETD (exponential) Drude current integrator
       c1  = exp( -ge*dt )
       if( ge*dt > 1.0d-6 )then
          c2 = ( 1.0d0 - c1 )/ge * wp2/( 4.0d0*pi_ )
       else                                       ! series form: (1-exp(-x))/x -> catastrophic
          c2 = dt*( 1.0d0 - 0.5d0*ge*dt ) * wp2/( 4.0d0*pi_ )   ! cancellation for x <~ 1e-8
       end if
    else                                          ! bilinear (Crank-Nicolson), mode 0
       c0  = 1.0d0 + ge*dt/2.0d0
       c1  = ( 1.0d0 - ge*dt/2.0d0 )/c0
       c2  = wp2*dt/( 4.0d0*pi_*c0 )
    end if
  end subroutine ttm3_drude_coef

  ! Cell wrapper: ADE Drude coefficients from the cell's carrier state (for the FDTD coupling).
  subroutine ttm3_drude_coef_cell( ix, iy, iz, omega, dt, c1, c2 )
    implicit none
    integer,intent(in)  :: ix,iy,iz
    real(8),intent(in)  :: omega, dt
    real(8),intent(out) :: c1, c2
    call ttm3_drude_coef( Ne(ix,iy,iz)+N_floor, Te(ix,iy,iz), Tl(ix,iy,iz), omega, dt, c1, c2 )
  end subroutine ttm3_drude_coef_cell

  ! Cell wrapper: raw Drude-oscillator parameters (eps_bg, gamma, wp^2) for the exact
  ! exponential sub-step (coupling_mode=3).
  subroutine ttm3_drude_osc_cell( ix, iy, iz, eps_bg, gamma, wp2 )
    implicit none
    integer,intent(in)  :: ix,iy,iz
    real(8),intent(out) :: eps_bg, gamma, wp2
    call ttm3_gamma_wp2( Ne(ix,iy,iz)+N_floor, Te(ix,iy,iz), Tl(ix,iy,iz), gamma, wp2 )
    eps_bg = tp3%eps_bg
  end subroutine ttm3_drude_osc_cell

  pure subroutine ttm3_permittivity( Ne_, Nh_, Te_, Th_, Tl_, omega, eps_re, sig )
    implicit none
    real(8),intent(in)  :: Ne_, Nh_, Te_, Th_, Tl_, omega
    real(8),intent(out) :: eps_re, sig
    real(8) :: wp2, bf, eps_d, sig_d
    if( tp3%sig_cold > 0.0d0 )then
       ! Layer C-a/b: band-filled bound interband part + strong carrier-T Drude.
       call ttm3_drude( Ne_, Te_, Th_, Tl_, omega, eps_d, sig_d )
       bf     = max( 1.0d0 - Ne_/tp3%N0, 0.0d0 )
       ! band filling depletes only the interband SUSCEPTIBILITY (eps_bg-1), not the
       ! vacuum contribution: full depletion must leave eps -> 1 + eps_Drude, not eps_Drude.
       eps_re = 1.0d0 + (tp3%eps_bg-1.0d0)*bf + eps_d
       sig    = tp3%sig_cold*bf + sig_d
    else
       wp2    = 4.0d0*pi_*( Ne_/tp3%mu_e + Nh_/tp3%mu_h )          ! legacy simple Drude
       eps_re = tp3%eps_bg - wp2/(omega*omega)
       sig    = wp2/(4.0d0*pi_)*(1.0d0/tp3%tau)/(omega*omega)
    end if
  end subroutine ttm3_permittivity

  !---------------------------------------------------------------------------
  ! Layer C-a carrier generation: only the BOUND (interband) fraction of the absorbed
  ! power creates electron-hole pairs.  sig_bound = sig_cold*(1-Ne/N0) depletes by band
  ! filling; sig_drude grows with carriers.  As the material metallises the bound fraction
  ! -> 0, so Drude absorption heats existing carriers instead of ionising -> generation
  ! self-limits (the reference's metallisation behaviour).  sig_cold<=0 reproduces the old
  ! "all absorbed power ionises" behaviour exactly (gen = max(power,0)/omega).
  pure function ttm3_linear_gen( ix, iy, iz, omega, power ) result( gen )
    implicit none
    integer,intent(in) :: ix,iy,iz
    real(8),intent(in) :: omega, power
    real(8) :: gen, bf, sig_b, sig_d, eps_d
    if( power <= 0.0d0 )then
       gen = 0.0d0; return
    end if
    if( tp3%sig_cold <= 0.0d0 )then
       gen = power/omega; return                                ! old behaviour (no split)
    end if
    bf    = max( 1.0d0 - Ne(ix,iy,iz)/tp3%N0, 0.0d0 )
    sig_b = tp3%sig_cold*bf
    call ttm3_drude( Ne(ix,iy,iz), Te(ix,iy,iz), Th(ix,iy,iz), Tl(ix,iy,iz), omega, eps_d, sig_d )
    gen   = power/omega * sig_b/( sig_b + sig_d + 1.0d-50 )      ! bound fraction only
  end function ttm3_linear_gen

  !---------------------------------------------------------------------------
  ! Stage 2 back-action: carrier-dependent permittivity/conductivity at a cell
  ! (drives the FDTD media coefficients; single-frequency static eps at omega,
  ! the same approximation the reference uses).
  subroutine ttm3_eps_sig( ix, iy, iz, omega, eps_re, sig )
    implicit none
    integer,intent(in)  :: ix,iy,iz
    real(8),intent(in)  :: omega
    real(8),intent(out) :: eps_re, sig
    call ttm3_permittivity( Ne(ix,iy,iz), Nh(ix,iy,iz), Te(ix,iy,iz), Th(ix,iy,iz), Tl(ix,iy,iz), omega, eps_re, sig )
  end subroutine ttm3_eps_sig

  pure function ttm3_coupling_mode() result( m )
    implicit none
    integer :: m
    m = tp3%coupling_mode
  end function ttm3_coupling_mode

  integer function ttm3_ninterior()
    implicit none
    ttm3_ninterior = ninterior
  end function ttm3_ninterior

  subroutine ttm3_interior_cell( m, ix, iy, iz )
    implicit none
    integer,intent(in)  :: m
    integer,intent(out) :: ix,iy,iz
    ix=ijk_interior(1,m); iy=ijk_interior(2,m); iz=ijk_interior(3,m)
  end subroutine ttm3_interior_cell

  !---------------------------------------------------------------------------
  ! Stage 3: carrier generation rate from the field-intensity envelope.
  ! Recommended two-photon form gen = beta2 * I^2, with I = |E|^2 + |E_g|^2
  ! (the cycle-averaged intensity from eh_get_field_envelope), atomic units.
  pure function ttm3_generation( intensity ) result( gen )
    implicit none
    real(8),intent(in) :: intensity
    real(8) :: gen
    gen = tp3%beta2 * intensity*intensity
  end function ttm3_generation

  !---------------------------------------------------------------------------
  ! Band-gap renormalisation: carrier-density shrinkage (dgap_c*Ne^(1/3)) plus the
  ! thermal (Varshni) shift dET = 7.02e-4*TlK^2/(TlK+1108)/E_h [Hartree], TlK = Tl in K.
  pure function ttm3_gap( Ne_, Tl_ ) result( Eg )
    implicit none
    real(8),intent(in) :: Ne_, Tl_
    real(8) :: Eg, TlK
    TlK = Tl_*hk_kelvin
    Eg = tp3%Egap - tp3%dgap_c*( max(Ne_,0.0d0) )**(1.0d0/3.0d0) &
                  - 7.02d-4*TlK*TlK/(TlK+1108.0d0)/27.2114d0
    Eg = max( Eg, 0.0d0 )   ! gap collapse: clamp at 0 (Mott handoff), matches standalone
  end function ttm3_gap

end module ttm3
