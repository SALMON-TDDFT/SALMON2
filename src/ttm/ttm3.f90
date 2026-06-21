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
  public :: ttm3_write_profile
  public :: ttm3_permittivity
  public :: ttm3_linear_gen
  public :: ttm3_eps_sig
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
     real(8) :: mob_e0   ! electron mobility                 [a.u. = m^2/Vs*2.353e-5]
     real(8) :: mob_h0   ! hole mobility                     [a.u.]
     ! Layer C-a (ablation): cold bound (interband) conductivity (optional 18th input line).
     ! When >0, carrier generation is restricted to the BOUND absorption fraction
     ! sig_bound/(sig_bound+sig_drude) with sig_bound=sig_cold*(1-Ne/N0) (band filling), so
     ! free-carrier (Drude) absorption heats rather than ionises -> generation self-limits in
     ! the metallisation regime.  <=0 disables the split (generation = full absorbed power).
     real(8) :: sig_cold ! cold interband conductivity       [a.u.]         (Layer C-a)
  end type ttm3_param

  integer,allocatable :: ijk_media_whole(:,:)
  integer,allocatable :: ijk_media_myrnk(:,:)
  integer,allocatable :: ijk_interior(:,:)   ! interior medium cells (all 6 neighbours same medium)
  integer :: nmedia_myrnk
  integer :: ninterior
  integer :: is_array(3), ie_array(3)
  integer :: is_inner(3), ie_inner(3)   ! inner (no-halo) bounds, for the space-charge prefix sum
  integer :: comm
  real(8) :: hgs(3)
  real(8) :: dt

  type(ttm3_param) :: tp3

  ! state and Runge-Kutta stage buffers
  real(8),allocatable :: Te(:,:,:), Th(:,:,:), Tl(:,:,:)
  real(8),allocatable :: Ne(:,:,:), Nh(:,:,:)
  ! transport increment buffers (Stage 4)
  real(8),allocatable :: rhs_te(:,:,:), rhs_th(:,:,:), rhs_tl(:,:,:)
  real(8),allocatable :: rhs_ne(:,:,:), rhs_nh(:,:,:)
  ! Set A (Fermi transport): per-cell thermal conductivity (electron/hole)
  real(8),allocatable :: Kc_e(:,:,:), Kc_h(:,:,:)
  ! Set A layer 2 (carrier transport): diffusion Dc, mobility Mc, current components
  ! Jc(species,dir)=grad N + coefjE grad gap + coefjT grad T, band gap, space-charge field.
  real(8),allocatable :: Dc_e(:,:,:), Dc_h(:,:,:), Mc_e(:,:,:), Mc_h(:,:,:)
  real(8),allocatable :: Jc_e(:,:,:,:), Jc_h(:,:,:,:)   ! (i,j,k,dir) -- dir last for contiguous halo
  real(8),allocatable :: rgap(:,:,:), Efld(:,:,:)       ! local gap, space-charge field (x)
  real(8),allocatable :: Cj_e(:,:,:), Cj_h(:,:,:)       ! Set A layer 3: heat per carrier (cJeh)

  ! Fermi-Dirac integral table F_j(eta).  Indices 1,2,3 = j = -1/2,1/2,3/2 (heat
  ! capacity / chemical potential); indices 4,5 = j = 1,2 (transport: mobility,
  ! diffusion, thermal conductivity).  F_0(eta)=ln(1+exp(eta)) is analytic (ttm3_F0),
  ! no table.  Built once at init by quadrature; interpolated in the inner loop (the
  ! same table approach as the reference's Table_Fermi).  eta grid: f_eta0 + k*f_deta.
  integer,parameter   :: f_neta = 2400
  real(8),parameter   :: f_eta0 = -60.0d0, f_deta = 0.05d0
  integer,parameter   :: f_nj   = 5
  real(8)             :: f_tab(0:f_neta,f_nj)
  logical             :: f_built = .false.

  ! density floor to keep the carrier heat capacity well defined.  Set to the
  ! reference's intrinsic seed density Ns ~ 5.4e-15 bohr^-3 (~3.4e10 cm^-3) so the
  ! carrier population can reach the reference's near-transparent regime; this is
  ! safe only because the carrier temperatures use the exponential integrator
  ! (the old RK4 step blew up at such low Ne -- see ttm3_step_cell).
  real(8),parameter :: N_floor = 5.4d-15   ! [1/bohr^3] (~3.4e10 cm^-3, reference seed)
  real(8),parameter :: T_clamp = 31.67d0   ! [a.u.] (~1e7 K) numerical temperature cap
  real(8),parameter :: pi_ = 3.14159265358979323846d0
  real(8),parameter :: hk_kelvin = 3.1577502480407d5   ! a.u. temperature -> Kelvin
  real(8),parameter :: T_au = 2.48036d5                ! [a.u.time] Auger saturation timescale

  character(12) :: ttm3_file = 'ttm3.inp_3tm'
  logical :: DISPLAY=.false.

contains

  !---------------------------------------------------------------------------
  subroutine init_ttm3_parameters( dt_em )
    use communication, only: comm_get_globalinfo, comm_is_root, comm_bcast
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
    DISPLAY = comm_is_root(npid)

    if ( DISPLAY ) write(*,'(a60)') repeat("-",30)//" init_ttm3_parameters(start)"

    inquire( FILE=ttm3_file, EXIST=flag )
    if( .not.flag )then
       if( DISPLAY )then
          write(*,*) "No three-temperature file ("//ttm3_file//")."
          write(*,*) "Three-temperature + carrier model is not used."
          write(*,'(a60)') repeat("-",30)//" init_ttm3_parameters(end  )"
       end if
       return
    else
       if( DISPLAY )then
          write(*,*) "Three-temperature file ("//ttm3_file//") is found."
          write(*,*) "Calculation uses the three-temperature + carrier model."
       end if
       use_ttm3 = .true.
    end if

    dt = dt_em

! Input parameters (one value per line):
!   Egap [eV], mu_e [-], mu_h [-], A_e [cm^6/s], A_h [cm^6/s],
!   tau [fs], Cl [J/(m^3 K)], Tini [K]
    if( comm_is_root(npid) )then
       open(unit, file=ttm3_file, status='old')
       read(unit,*) tp3%Egap
       read(unit,*) tp3%mu_e
       read(unit,*) tp3%mu_h
       read(unit,*) tp3%A_e
       read(unit,*) tp3%A_h
       read(unit,*) tp3%tau
       read(unit,*) tp3%Cl
       read(unit,*) tp3%Tini
       read(unit,*) tp3%eps_bg
       read(unit,*) tp3%N0
       read(unit,*) tp3%beta2
       read(unit,*) tp3%Ddiff
       read(unit,*) tp3%kappa_e
       read(unit,*) tp3%kappa_l
       read(unit,*) tp3%dgap_c
       tp3%mob_e0 = 0.0d0; tp3%mob_h0 = 0.0d0     ! optional (Set A Fermi transport)
       read(unit,*,iostat=ios) tp3%mob_e0
       read(unit,*,iostat=ios) tp3%mob_h0
       tp3%sig_cold = 0.0d0                        ! optional (Layer C-a bound/Drude split); <=0 disables
       read(unit,*,iostat=ios) tp3%sig_cold
       close(unit)
       write(*,*) "Egap[eV]    =",tp3%Egap
       write(*,*) "mu_e        =",tp3%mu_e
       write(*,*) "mu_h        =",tp3%mu_h
       write(*,*) "A_e[cm6/s]  =",tp3%A_e
       write(*,*) "A_h[cm6/s]  =",tp3%A_h
       write(*,*) "tau[fs]     =",tp3%tau
       write(*,*) "Cl[J/m3K]   =",tp3%Cl
       write(*,*) "Tini[K]     =",tp3%Tini
       write(*,*) "eps_bg      =",tp3%eps_bg
       write(*,*) "N0[cm-3]    =",tp3%N0
       write(*,*) "beta2[a.u.] =",tp3%beta2
       write(*,*) "Ddiff[cm2/s]=",tp3%Ddiff
       write(*,*) "kappa_e[W/mK]=",tp3%kappa_e
       write(*,*) "kappa_l[W/mK]=",tp3%kappa_l
       write(*,*) "dgap_c[a.u.]=",tp3%dgap_c
    end if

    call comm_bcast(tp3%Egap,comm,0)
    call comm_bcast(tp3%mu_e,comm,0)
    call comm_bcast(tp3%mu_h,comm,0)
    call comm_bcast(tp3%A_e ,comm,0)
    call comm_bcast(tp3%A_h ,comm,0)
    call comm_bcast(tp3%tau ,comm,0)
    call comm_bcast(tp3%Cl  ,comm,0)
    call comm_bcast(tp3%Tini,comm,0)
    call comm_bcast(tp3%eps_bg ,comm,0)
    call comm_bcast(tp3%N0     ,comm,0)
    call comm_bcast(tp3%beta2  ,comm,0)
    call comm_bcast(tp3%Ddiff  ,comm,0)
    call comm_bcast(tp3%kappa_e,comm,0)
    call comm_bcast(tp3%kappa_l,comm,0)
    call comm_bcast(tp3%dgap_c ,comm,0)
    call comm_bcast(tp3%mob_e0 ,comm,0)
    call comm_bcast(tp3%mob_h0 ,comm,0)
    call comm_bcast(tp3%sig_cold,comm,0)

! Convert to atomic units
    tp3%Egap = tp3%Egap / hartree_ev
    ! mobility [m^2/Vs] -> a.u.  (mob_au = mob_SI * E_h * a_t * 1e-5 / a_B^2, the reference's
    ! conversion: 8.5e-3 -> 2.0e-7).  Kept in SI-derived form via the atomic-unit constants.
    tp3%mob_e0 = tp3%mob_e0 * hartree_ev*(atomic_unit_of_time*1.0d15)*1.0d-5/(atomic_unit_of_length*1.0d10)**2
    tp3%mob_h0 = tp3%mob_h0 * hartree_ev*(atomic_unit_of_time*1.0d15)*1.0d-5/(atomic_unit_of_length*1.0d10)**2
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
    implicit none
    real(8),intent(in) :: hgs_in(3)
    integer,intent(in) :: is_a(3), is(3), ie(3)
    integer,intent(in) :: imedia(is_a(1):,is_a(2):,is_a(3):)
    integer :: ix,iy,iz,m
    logical :: allmed

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
  end subroutine init_ttm3_grid

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
       allocate( Jc_e(i1:j1,i2:j2,i3:j3,3), Jc_h(i1:j1,i2:j2,i3:j3,3) )
       allocate( rgap(i1:j1,i2:j2,i3:j3), Efld(i1:j1,i2:j2,i3:j3) )
       allocate( Cj_e(i1:j1,i2:j2,i3:j3), Cj_h(i1:j1,i2:j2,i3:j3) )
       Kc_e = 0.0d0  ; Kc_h = 0.0d0
       Dc_e = 0.0d0  ; Dc_h = 0.0d0 ; Mc_e = 0.0d0 ; Mc_h = 0.0d0
       Jc_e = 0.0d0  ; Jc_h = 0.0d0 ; rgap = tp3%Egap ; Efld = 0.0d0
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
    real(8) :: eta,umax,du,u,integ,s,gam(f_nj),jval(f_nj)
    ! indices 1..5 = j = -1/2, 1/2, 3/2, 1, 2  (gam = Gamma(j+1))
    jval = (/ -0.5d0, 0.5d0, 1.5d0, 1.0d0, 2.0d0 /)
    gam  = (/ 1.7724538509055160d0, 0.8862269254527580d0, 1.3293403881791370d0, &
              1.0d0, 2.0d0 /)                                            ! Gamma(2)=1, Gamma(3)=2
    do k=0,f_neta
       eta = f_eta0 + dble(k)*f_deta
       umax = sqrt( max(eta,0.0d0) + 50.0d0 )
       du = umax/dble(nq)
       do jj=1,f_nj
          integ = 0.0d0
          do q=0,nq
             u = dble(q)*du
             s = u**(2.0d0*jval(jj)+1.0d0)/( exp(u*u-eta) + 1.0d0 )
             if(q==0 .or. q==nq)then       ; integ=integ+s
             elseif(mod(q,2)==1)then        ; integ=integ+4.0d0*s
             else                           ; integ=integ+2.0d0*s
             end if
          end do
          f_tab(k,jj) = 2.0d0*integ*du/3.0d0/gam(jj)
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
  pure function ttm3_fermi( jidx, eta ) result( F )
    implicit none
    integer,intent(in) :: jidx
    real(8),intent(in) :: eta
    real(8) :: F,x,w
    integer :: k
    x = (eta - f_eta0)/f_deta
    if(x <= 0.0d0)then            ; k=0       ; w=0.0d0
    elseif(x >= dble(f_neta))then ; k=f_neta-1; w=1.0d0
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
    real(8) :: eta,Neff,y,w
    integer :: lo,hi,mid
    Neff = 2.0d0*( mc*max(T,1.0d-12)/(2.0d0*pi_) )**1.5d0
    y = N/Neff                                       ! target value of F_{1/2}(eta)
    if( y <= f_tab(0,2) )then
       eta = f_eta0; return
    elseif( y >= f_tab(f_neta,2) )then
       eta = f_eta0 + dble(f_neta)*f_deta; return
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

    ! generation, saturated at the reference/atomic density N0 (cannot exceed it)
    geff = gen_ * max( 0.0d0, 1.0d0 - Ne_/tp3%N0 )

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
    geff = gen_ * max( 0.0d0, 1.0d0 - Ne_/tp3%N0 )
    Cef_e = tp3%A_e*Ne_*Nh_ ; Re = Ne_*Cef_e/( 1.0d0 + T_au*Cef_e )   ! saturated Auger (e)
    Cef_h = tp3%A_h*Ne_*Nh_ ; Rh = Nh_*Cef_h/( 1.0d0 + T_au*Cef_h )   ! saturated Auger (h)
    R    = Re + Rh
    red_e = tp3%mu_h/(tp3%mu_e+tp3%mu_h)
    red_h = tp3%mu_e/(tp3%mu_e+tp3%mu_h)
    gap = ttm3_gap( Ne_, Tl_ )                                   ! carrier + thermal (Varshni) gap
    q_heat = max( source_ - geff*gap, 0.0d0 )
    ! thermal across-gap (collisional) carrier generation, Pauli-blocked
    Col = 3.6d-5*2.419d-2*0.5d0*( Ne_*exp(-1.5d0*gap/Te_) + Nh_*exp(-1.5d0*gap/Th_) ) &
          * max( 0.0d0, 1.0d0 - Ne_/tp3%N0 )

    ! lattice: relaxation heating from the carriers (non-stiff, Cl large), with the
    ! temperature-dependent lattice heat capacity Cl(Tl) (= input value at 300 K)
    TlK    = Tl_*hk_kelvin
    Cl_eff = tp3%Cl*( 1.978d0 + 3.54d-4*TlK - 3.68d0/(TlK*TlK) )/2.084d0
    dTl = ( Ce*(Te_-Tl_) + Ch*(Th_-Tl_) )*inv_tau / Cl_eff

    ! carrier temperatures: exponential integrator (stiff relaxation + dilution)
    rate_e = inv_tau + geff/(Ne_+N_floor)
    rate_h = inv_tau + geff/(Nh_+N_floor)
    Te_star = Tl_ + (red_e*q_heat/Ce)/rate_e
    Th_star = Tl_ + (red_h*q_heat/Ch)/rate_h
    Te_ = Te_star + (Te_-Te_star)*exp( -rate_e*dt_in )
    Th_ = Th_star + (Th_-Th_star)*exp( -rate_h*dt_in )

    ! carrier densities and lattice: explicit Euler (non-stiff); +Col thermal generation
    Ne_ = max( Ne_ + dt_in*(geff - R + Col), 0.0d0 )
    Nh_ = max( Nh_ + dt_in*(geff - R + Col), 0.0d0 )
    Tl_ = max( Tl_ + dt_in*dTl, 0.0d0 )

    ! numerical guard (should not trigger now that the stiff modes are integrated exactly)
    Te_ = min(max(Te_,0.0d0),T_clamp); Th_ = min(max(Th_,0.0d0),T_clamp)
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
  ! Stage 4: carrier diffusion + heat conduction (recommended explicit scheme,
  ! operator-split from the local step).  Skipped when all coefficients are zero
  ! (then ttm3_main reproduces the verified Stage-1 local engine exactly).
  subroutine ttm3_transport( srg, rg )
    use structures,    only: s_rgrid, s_sendrecv_grid
    use sendrecv_grid, only: update_overlap_real8
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_rgrid),         intent(in)    :: rg
    integer :: m,ix,iy,iz
    real(8) :: cx,cy,cz,Ce,Ch,lNe,lNh,lTe,lTh,lTl,Dmax,Dd
    real(8) :: gx,gy,gz,gKe,gKh
    real(8) :: eta,F0,Fm,F12,F1,F2,kT,cjE,cjT,fp4i
    real(8) :: dJe,dJh,gDJe,gDJh,gMxe,gMxh,gNxe,gNxh,dEf,nabNe,nabNh,gCJe,gCJh
    logical :: use_fermi

    use_fermi = ( tp3%mob_e0>0.0d0 .or. tp3%mob_h0>0.0d0 )
    if( tp3%Ddiff<=0.0d0 .and. tp3%kappa_e<=0.0d0 .and. tp3%kappa_l<=0.0d0 &
        .and. .not.use_fermi ) return

    cx=1.0d0/hgs(1)**2; cy=1.0d0/hgs(2)**2; cz=1.0d0/hgs(3)**2
    gx=0.5d0/hgs(1); gy=0.5d0/hgs(2); gz=0.5d0/hgs(3)   ! central first-difference
    ! explicit FTCS diffusion stability: cap every diffusivity so dt*D*2*(cx+cy+cz)<=0.5.
    Dmax = 0.5d0/( 2.0d0*dt*(cx+cy+cz) )
    Dd   = min(tp3%Ddiff, Dmax)
    fp4i = 4.0d0*pi_/tp3%eps_bg                ! 4*pi*inveps (local Gauss term, dEf)

    ! Set A: per-cell Fermi-Dirac transport coefficients + carrier current components.
    if( use_fermi )then
!$omp parallel do private(m,ix,iy,iz)
       do m=1,nmedia_myrnk
          ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
          rgap(ix,iy,iz)=ttm3_gap( Ne(ix,iy,iz), Tl(ix,iy,iz) )    ! carrier + thermal renormalised gap
       end do
!$omp end parallel do
       call update_overlap_real8(srg, rg, rgap)
    end if

    call update_overlap_real8(srg, rg, Ne); call update_overlap_real8(srg, rg, Nh)
    call update_overlap_real8(srg, rg, Te); call update_overlap_real8(srg, rg, Th)
    call update_overlap_real8(srg, rg, Tl)

    if( use_fermi )then
       ! pass 1: Fermi coefficients (mobility Mc, diffusion Dc, conductivity Kc) and the
       ! current components Jc(dir)=grad N + coefjE grad gap + coefjT grad T (per species).
!$omp parallel do private(m,ix,iy,iz,eta,F0,Fm,F12,F1,F2,kT,cjE,cjT)
       do m=1,nmedia_myrnk
          ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
          ! electrons
          kT=Te(ix,iy,iz); eta=ttm3_chem_pot(kT,tp3%mu_e,Ne(ix,iy,iz)+N_floor)
          F0=ttm3_F0(eta); Fm=ttm3_fermi(1,eta); F12=ttm3_fermi(2,eta)
          F1=ttm3_fermi(4,eta); F2=ttm3_fermi(5,eta)
          Mc_e(ix,iy,iz)=tp3%mob_e0*F0/F12
          Dc_e(ix,iy,iz)=tp3%mob_e0*kT*F0/Fm
          Kc_e(ix,iy,iz)=(6.0d0*F2/F0-4.0d0*(F1/F0)**2)*Ne(ix,iy,iz)*kT*Mc_e(ix,iy,iz)
          cjE=(tp3%mu_h/(tp3%mu_e+tp3%mu_h))*(Ne(ix,iy,iz)+N_floor)*Fm/F12/kT
          cjT=(Ne(ix,iy,iz)+N_floor)*(2.0d0*Fm*F1/(F12*F0)-1.5d0)/kT
          Jc_e(ix,iy,iz,1)=(Ne(ix+1,iy,iz)-Ne(ix-1,iy,iz))*gx + cjE*(rgap(ix+1,iy,iz)-rgap(ix-1,iy,iz))*gx + cjT*(Te(ix+1,iy,iz)-Te(ix-1,iy,iz))*gx
          Jc_e(ix,iy,iz,2)=(Ne(ix,iy+1,iz)-Ne(ix,iy-1,iz))*gy + cjE*(rgap(ix,iy+1,iz)-rgap(ix,iy-1,iz))*gy + cjT*(Te(ix,iy+1,iz)-Te(ix,iy-1,iz))*gy
          Jc_e(ix,iy,iz,3)=(Ne(ix,iy,iz+1)-Ne(ix,iy,iz-1))*gz + cjE*(rgap(ix,iy,iz+1)-rgap(ix,iy,iz-1))*gz + cjT*(Te(ix,iy,iz+1)-Te(ix,iy,iz-1))*gz
          Cj_e(ix,iy,iz)=2.0d0*kT*F1/F0 + (tp3%mu_h/(tp3%mu_e+tp3%mu_h))*rgap(ix,iy,iz)   ! heat per carrier
          ! holes
          kT=Th(ix,iy,iz); eta=ttm3_chem_pot(kT,tp3%mu_h,Nh(ix,iy,iz)+N_floor)
          F0=ttm3_F0(eta); Fm=ttm3_fermi(1,eta); F12=ttm3_fermi(2,eta)
          F1=ttm3_fermi(4,eta); F2=ttm3_fermi(5,eta)
          Mc_h(ix,iy,iz)=tp3%mob_h0*F0/F12
          Dc_h(ix,iy,iz)=tp3%mob_h0*kT*F0/Fm
          Kc_h(ix,iy,iz)=(6.0d0*F2/F0-4.0d0*(F1/F0)**2)*Nh(ix,iy,iz)*kT*Mc_h(ix,iy,iz)
          cjE=(tp3%mu_e/(tp3%mu_e+tp3%mu_h))*(Nh(ix,iy,iz)+N_floor)*Fm/F12/kT
          cjT=(Nh(ix,iy,iz)+N_floor)*(2.0d0*Fm*F1/(F12*F0)-1.5d0)/kT
          Jc_h(ix,iy,iz,1)=(Nh(ix+1,iy,iz)-Nh(ix-1,iy,iz))*gx + cjE*(rgap(ix+1,iy,iz)-rgap(ix-1,iy,iz))*gx + cjT*(Th(ix+1,iy,iz)-Th(ix-1,iy,iz))*gx
          Jc_h(ix,iy,iz,2)=(Nh(ix,iy+1,iz)-Nh(ix,iy-1,iz))*gy + cjE*(rgap(ix,iy+1,iz)-rgap(ix,iy-1,iz))*gy + cjT*(Th(ix,iy+1,iz)-Th(ix,iy-1,iz))*gy
          Jc_h(ix,iy,iz,3)=(Nh(ix,iy,iz+1)-Nh(ix,iy,iz-1))*gz + cjE*(rgap(ix,iy,iz+1)-rgap(ix,iy,iz-1))*gz + cjT*(Th(ix,iy,iz+1)-Th(ix,iy,iz-1))*gz
          Cj_h(ix,iy,iz)=2.0d0*kT*F1/F0 + (tp3%mu_e/(tp3%mu_e+tp3%mu_h))*rgap(ix,iy,iz)   ! heat per carrier
       end do
!$omp end parallel do
       call update_overlap_real8(srg, rg, Dc_e); call update_overlap_real8(srg, rg, Dc_h)
       call update_overlap_real8(srg, rg, Mc_e); call update_overlap_real8(srg, rg, Mc_h)
       call update_overlap_real8(srg, rg, Kc_e); call update_overlap_real8(srg, rg, Kc_h)
       call update_overlap_real8(srg, rg, Cj_e); call update_overlap_real8(srg, rg, Cj_h)
       call update_overlap_real8(srg, rg, Jc_e(:,:,:,1)); call update_overlap_real8(srg, rg, Jc_e(:,:,:,2)); call update_overlap_real8(srg, rg, Jc_e(:,:,:,3))
       call update_overlap_real8(srg, rg, Jc_h(:,:,:,1)); call update_overlap_real8(srg, rg, Jc_h(:,:,:,2)); call update_overlap_real8(srg, rg, Jc_h(:,:,:,3))
       call ttm3_space_charge()                 ! Efld = space-charge field (prefix sum along x)
    end if

!$omp parallel do private(m,ix,iy,iz,Ce,Ch,lNe,lNh,lTe,lTh,lTl,gKe,gKh, &
!$omp&                    dJe,dJh,gDJe,gDJh,gMxe,gMxh,gNxe,gNxh,dEf,nabNe,nabNh,gCJe,gCJh)
    do m=1,nmedia_myrnk
       ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
       lTe=(Te(ix+1,iy,iz)-2*Te(ix,iy,iz)+Te(ix-1,iy,iz))*cx+(Te(ix,iy+1,iz)-2*Te(ix,iy,iz)+Te(ix,iy-1,iz))*cy+(Te(ix,iy,iz+1)-2*Te(ix,iy,iz)+Te(ix,iy,iz-1))*cz
       lTh=(Th(ix+1,iy,iz)-2*Th(ix,iy,iz)+Th(ix-1,iy,iz))*cx+(Th(ix,iy+1,iz)-2*Th(ix,iy,iz)+Th(ix,iy-1,iz))*cy+(Th(ix,iy,iz+1)-2*Th(ix,iy,iz)+Th(ix,iy,iz-1))*cz
       lTl=(Tl(ix+1,iy,iz)-2*Tl(ix,iy,iz)+Tl(ix-1,iy,iz))*cx+(Tl(ix,iy+1,iz)-2*Tl(ix,iy,iz)+Tl(ix,iy-1,iz))*cy+(Tl(ix,iy,iz+1)-2*Tl(ix,iy,iz)+Tl(ix,iy,iz-1))*cz
       Ce=ttm3_heat_capacity(Te(ix,iy,iz),tp3%mu_e,Ne(ix,iy,iz)+N_floor)
       Ch=ttm3_heat_capacity(Th(ix,iy,iz),tp3%mu_h,Nh(ix,iy,iz)+N_floor)
       if( use_fermi )then
          ! heat conduction nab.(K nabT) = nabK.nabT + K lap(T), /Ce (K/Ce CFL-capped)
          gKe=(Kc_e(ix+1,iy,iz)-Kc_e(ix-1,iy,iz))*gx*(Te(ix+1,iy,iz)-Te(ix-1,iy,iz))*gx &
             +(Kc_e(ix,iy+1,iz)-Kc_e(ix,iy-1,iz))*gy*(Te(ix,iy+1,iz)-Te(ix,iy-1,iz))*gy &
             +(Kc_e(ix,iy,iz+1)-Kc_e(ix,iy,iz-1))*gz*(Te(ix,iy,iz+1)-Te(ix,iy,iz-1))*gz
          gKh=(Kc_h(ix+1,iy,iz)-Kc_h(ix-1,iy,iz))*gx*(Th(ix+1,iy,iz)-Th(ix-1,iy,iz))*gx &
             +(Kc_h(ix,iy+1,iz)-Kc_h(ix,iy-1,iz))*gy*(Th(ix,iy+1,iz)-Th(ix,iy-1,iz))*gy &
             +(Kc_h(ix,iy,iz+1)-Kc_h(ix,iy,iz-1))*gz*(Th(ix,iy,iz+1)-Th(ix,iy,iz-1))*gz
          rhs_te(ix,iy,iz)=dt*( min(Kc_e(ix,iy,iz)/Ce,Dmax)*lTe + gKe/Ce )
          rhs_th(ix,iy,iz)=dt*( min(Kc_h(ix,iy,iz)/Ch,Dmax)*lTh + gKh/Ch )
          ! carrier transport: nab.(Dc Jc) - drift(space charge).  reference Evolve_N_T.
          dJe=(Jc_e(ix+1,iy,iz,1)-Jc_e(ix-1,iy,iz,1))*gx+(Jc_e(ix,iy+1,iz,2)-Jc_e(ix,iy-1,iz,2))*gy+(Jc_e(ix,iy,iz+1,3)-Jc_e(ix,iy,iz-1,3))*gz
          dJh=(Jc_h(ix+1,iy,iz,1)-Jc_h(ix-1,iy,iz,1))*gx+(Jc_h(ix,iy+1,iz,2)-Jc_h(ix,iy-1,iz,2))*gy+(Jc_h(ix,iy,iz+1,3)-Jc_h(ix,iy,iz-1,3))*gz
          gDJe=(Dc_e(ix+1,iy,iz)-Dc_e(ix-1,iy,iz))*gx*Jc_e(ix,iy,iz,1)+(Dc_e(ix,iy+1,iz)-Dc_e(ix,iy-1,iz))*gy*Jc_e(ix,iy,iz,2)+(Dc_e(ix,iy,iz+1)-Dc_e(ix,iy,iz-1))*gz*Jc_e(ix,iy,iz,3)
          gDJh=(Dc_h(ix+1,iy,iz)-Dc_h(ix-1,iy,iz))*gx*Jc_h(ix,iy,iz,1)+(Dc_h(ix,iy+1,iz)-Dc_h(ix,iy-1,iz))*gy*Jc_h(ix,iy,iz,2)+(Dc_h(ix,iy,iz+1)-Dc_h(ix,iy,iz-1))*gz*Jc_h(ix,iy,iz,3)
          ! Layer 3: convective heat -- the carrier current carries energy cJeh per carrier:
          ! nabW += cJeh*(Dc div Jc + grad Dc . Jc) + Dc*(grad cJeh . Jc).  Added to the heat RHS.
          gCJe=(Cj_e(ix+1,iy,iz)-Cj_e(ix-1,iy,iz))*gx*Jc_e(ix,iy,iz,1)+(Cj_e(ix,iy+1,iz)-Cj_e(ix,iy-1,iz))*gy*Jc_e(ix,iy,iz,2)+(Cj_e(ix,iy,iz+1)-Cj_e(ix,iy,iz-1))*gz*Jc_e(ix,iy,iz,3)
          gCJh=(Cj_h(ix+1,iy,iz)-Cj_h(ix-1,iy,iz))*gx*Jc_h(ix,iy,iz,1)+(Cj_h(ix,iy+1,iz)-Cj_h(ix,iy-1,iz))*gy*Jc_h(ix,iy,iz,2)+(Cj_h(ix,iy,iz+1)-Cj_h(ix,iy,iz-1))*gz*Jc_h(ix,iy,iz,3)
          rhs_te(ix,iy,iz)=rhs_te(ix,iy,iz)+dt*( Cj_e(ix,iy,iz)*(Dc_e(ix,iy,iz)*dJe+gDJe) + Dc_e(ix,iy,iz)*gCJe )/Ce
          rhs_th(ix,iy,iz)=rhs_th(ix,iy,iz)+dt*( Cj_h(ix,iy,iz)*(Dc_h(ix,iy,iz)*dJh+gDJh) + Dc_h(ix,iy,iz)*gCJh )/Ch
          gMxe=(Mc_e(ix+1,iy,iz)-Mc_e(ix-1,iy,iz))*gx; gNxe=(Ne(ix+1,iy,iz)-Ne(ix-1,iy,iz))*gx
          gMxh=(Mc_h(ix+1,iy,iz)-Mc_h(ix-1,iy,iz))*gx; gNxh=(Nh(ix+1,iy,iz)-Nh(ix-1,iy,iz))*gx
          dEf=fp4i*(Ne(ix,iy,iz)-Nh(ix,iy,iz))               ! 4*pi*inveps*(Ne-Nh)
          nabNe=min(Dc_e(ix,iy,iz),Dmax)*dJe + gDJe - (gMxe*Ne(ix,iy,iz)+Mc_e(ix,iy,iz)*gNxe)*Efld(ix,iy,iz) - Mc_e(ix,iy,iz)*Ne(ix,iy,iz)*dEf
          nabNh=min(Dc_h(ix,iy,iz),Dmax)*dJh + gDJh + (gMxh*Nh(ix,iy,iz)+Mc_h(ix,iy,iz)*gNxh)*Efld(ix,iy,iz) + Mc_h(ix,iy,iz)*Nh(ix,iy,iz)*dEf
          rhs_ne(ix,iy,iz)=dt*nabNe
          rhs_nh(ix,iy,iz)=dt*nabNh
       else
          lNe=(Ne(ix+1,iy,iz)-2*Ne(ix,iy,iz)+Ne(ix-1,iy,iz))*cx+(Ne(ix,iy+1,iz)-2*Ne(ix,iy,iz)+Ne(ix,iy-1,iz))*cy+(Ne(ix,iy,iz+1)-2*Ne(ix,iy,iz)+Ne(ix,iy,iz-1))*cz
          lNh=(Nh(ix+1,iy,iz)-2*Nh(ix,iy,iz)+Nh(ix-1,iy,iz))*cx+(Nh(ix,iy+1,iz)-2*Nh(ix,iy,iz)+Nh(ix,iy-1,iz))*cy+(Nh(ix,iy,iz+1)-2*Nh(ix,iy,iz)+Nh(ix,iy,iz-1))*cz
          rhs_ne(ix,iy,iz)=dt*Dd*lNe
          rhs_nh(ix,iy,iz)=dt*Dd*lNh
          rhs_te(ix,iy,iz)=dt*min(tp3%kappa_e/Ce,Dmax)*lTe
          rhs_th(ix,iy,iz)=dt*min(tp3%kappa_e/Ch,Dmax)*lTh
       end if
       rhs_tl(ix,iy,iz)=dt*min( tp3%kappa_l*(Tl(ix,iy,iz)*hk_kelvin/300.0d0)**(-1.23d0)/tp3%Cl, Dmax )*lTl  ! Kl(Tl)
    end do
!$omp end parallel do

!$omp parallel do private(m,ix,iy,iz)
    do m=1,nmedia_myrnk
       ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
       Ne(ix,iy,iz)=max(Ne(ix,iy,iz)+rhs_ne(ix,iy,iz),0.0d0)
       Nh(ix,iy,iz)=max(Nh(ix,iy,iz)+rhs_nh(ix,iy,iz),0.0d0)
       Te(ix,iy,iz)=Te(ix,iy,iz)+rhs_te(ix,iy,iz)
       Th(ix,iy,iz)=Th(ix,iy,iz)+rhs_th(ix,iy,iz)
       Tl(ix,iy,iz)=Tl(ix,iy,iz)+rhs_tl(ix,iy,iz)
    end do
!$omp end parallel do
  end subroutine ttm3_transport

  !---------------------------------------------------------------------------
  ! Space-charge field Efld along x: Ef(ix)=2*pi*inveps*[sum_{i<ix} - sum_{i>ix}] rho dV,
  ! rho=(Nh-Ne) (charge density; zero in non-medium cells, Ne=Nh=floor).  Per (iy,iz)
  ! line, a prefix sum over the local inner x-range (single-domain; nproc_rgrid=1).
  subroutine ttm3_space_charge()
    implicit none
    integer :: ix,iy,iz
    real(8) :: dV,c2,total,pre,rho
    dV = hgs(1)*hgs(2)*hgs(3)
    c2 = 2.0d0*pi_/tp3%eps_bg
!$omp parallel do collapse(2) private(ix,iy,iz,total,pre,rho)
    do iz=is_inner(3),ie_inner(3)
    do iy=is_inner(2),ie_inner(2)
       total=0.0d0
       do ix=is_inner(1),ie_inner(1)
          total=total+(Nh(ix,iy,iz)-Ne(ix,iy,iz))
       end do
       pre=0.0d0                                  ! exclusive prefix sum (i<ix)
       do ix=is_inner(1),ie_inner(1)
          rho=Nh(ix,iy,iz)-Ne(ix,iy,iz)
          ! Ef = (left - right)*dV*c2 ;  right = total - pre - rho
          Efld(ix,iy,iz)=c2*dV*( 2.0d0*pre + rho - total )
          pre=pre+rho
       end do
    end do
    end do
!$omp end parallel do
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
    implicit none
    real(8),intent(out) :: Te_o, Th_o, Tl_o, Ne_o, Nh_o   ! [K],[K],[K],[1/cm^3],[1/cm^3]
    real(8), parameter :: hartree_kelvin_relationship = 3.1577502480407d5
    real(8), parameter :: atomic_unit_of_length = 5.29177210903d-11
    real(8) :: cm3_per_bohr3
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
  end subroutine ttm3_get_max

  !---------------------------------------------------------------------------
  ! State at the illuminated front face: the minimum-ix medium cell nearest the
  ! transverse centroid.  This is the analogue of the reference's x=0 surface cell
  ! and, unlike ttm3_get_max, is not biased by the Fabry-Perot field antinodes that
  ! form inside a finite slab (which sit deeper than the front face and read high).
  subroutine ttm3_get_front( Te_o, Th_o, Tl_o, Ne_o, Nh_o )
    implicit none
    real(8),intent(out) :: Te_o, Th_o, Tl_o, Ne_o, Nh_o   ! [K],[K],[K],[1/cm^3],[1/cm^3]
    real(8), parameter :: hartree_kelvin_relationship = 3.1577502480407d5
    real(8), parameter :: atomic_unit_of_length = 5.29177210903d-11
    real(8) :: cm3_per_bohr3, cy, cz
    integer :: m,ix,iy,iz,ixmin,nmin,a(3),dbest,d
    Te_o=tp3%Tini*hartree_kelvin_relationship; Th_o=Te_o; Tl_o=Te_o; Ne_o=0d0; Nh_o=0d0
    if( nmedia_myrnk<=0 ) return
    cm3_per_bohr3 = (atomic_unit_of_length*1.0d2)**3
    ! minimum ix among medium cells = the illuminated face
    ixmin = ijk_media_myrnk(1,1)
    do m=2,nmedia_myrnk
       if( ijk_media_myrnk(1,m) < ixmin ) ixmin = ijk_media_myrnk(1,m)
    end do
    ! transverse centroid of that face
    cy=0d0; cz=0d0; nmin=0
    do m=1,nmedia_myrnk
       if( ijk_media_myrnk(1,m)==ixmin )then
          cy=cy+ijk_media_myrnk(2,m); cz=cz+ijk_media_myrnk(3,m); nmin=nmin+1
       end if
    end do
    cy=cy/max(nmin,1); cz=cz/max(nmin,1)
    ! face cell nearest the centroid
    a = ijk_media_myrnk(1:3,1); dbest=huge(0)
    do m=1,nmedia_myrnk
       if( ijk_media_myrnk(1,m)==ixmin )then
          iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
          d = (iy-nint(cy))**2 + (iz-nint(cz))**2
          if( d<dbest )then; dbest=d; a=ijk_media_myrnk(1:3,m); end if
       end if
    end do
    Te_o=Te(a(1),a(2),a(3))*hartree_kelvin_relationship
    Th_o=Th(a(1),a(2),a(3))*hartree_kelvin_relationship
    Tl_o=Tl(a(1),a(2),a(3))*hartree_kelvin_relationship
    Ne_o=Ne(a(1),a(2),a(3))/cm3_per_bohr3; Nh_o=Nh(a(1),a(2),a(3))/cm3_per_bohr3
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

  ! Depth profile along x at the transverse centre line (front cell's iy,iz): writes
  ! one block per call (time, ix, Te[K], Ne[cm^-3], Tl[K]) for spatial-distribution diagnostics.
  subroutine ttm3_write_profile( fname, t_fs, u, field )
    implicit none
    character(*),intent(in) :: fname
    real(8),     intent(in) :: t_fs
    integer,     intent(in) :: u
    real(8),     intent(in) :: field(is_array(1):,is_array(2):,is_array(3):)  ! |E|^2 envelope [a.u.]
    real(8),parameter :: hk = 3.1577502480407d5, aul = 5.29177210903d-11
    real(8) :: cm3
    integer :: m,ix,iy,iz,fx,fy,fz
    if( nmedia_myrnk<=0 ) return
    call ttm3_front_ijk( fx, fy, fz )                 ! transverse centre = front cell line
    cm3 = (aul*1.0d2)**3
    open(u, file=fname, status='unknown', position='append')
    do ix=is_inner(1),ie_inner(1)                     ! sorted by depth
       do m=1,nmedia_myrnk
          if( ijk_media_myrnk(1,m)==ix .and. ijk_media_myrnk(2,m)==fy .and. ijk_media_myrnk(3,m)==fz )then
             iy=fy; iz=fz
             write(u,'(F12.4,1X,I6,4(1X,E16.8))') t_fs, ix, &
                  Te(ix,iy,iz)*hk, Ne(ix,iy,iz)/cm3, Tl(ix,iy,iz)*hk, field(ix,iy,iz)
          end if
       end do
    end do
    write(u,*) ' '                                    ! blank line separates time blocks
    close(u)
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

  pure subroutine ttm3_permittivity( Ne_, Nh_, Te_, Th_, Tl_, omega, eps_re, sig )
    implicit none
    real(8),intent(in)  :: Ne_, Nh_, Te_, Th_, Tl_, omega
    real(8),intent(out) :: eps_re, sig
    real(8) :: wp2, bf, eps_d, sig_d
    if( tp3%sig_cold > 0.0d0 )then
       ! Layer C-a/b: band-filled bound interband part + strong carrier-T Drude.
       call ttm3_drude( Ne_, Te_, Th_, Tl_, omega, eps_d, sig_d )
       bf     = max( 1.0d0 - Ne_/tp3%N0, 0.0d0 )
       eps_re = tp3%eps_bg*bf + eps_d
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
  end function ttm3_gap

end module ttm3
