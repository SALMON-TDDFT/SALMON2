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
  public :: ttm3_permittivity
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
  end type ttm3_param

  integer,allocatable :: ijk_media_whole(:,:)
  integer,allocatable :: ijk_media_myrnk(:,:)
  integer,allocatable :: ijk_interior(:,:)   ! interior medium cells (all 6 neighbours same medium)
  integer :: nmedia_myrnk
  integer :: ninterior
  integer :: is_array(3), ie_array(3)
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

  ! Fermi-Dirac integral table F_j(eta) for j = -1/2, 1/2, 3/2 (indices 1,2,3),
  ! built once at init by quadrature, interpolated in the inner loop (the same
  ! table approach as the reference's Table_Fermi).  eta grid: f_eta0 + k*f_deta.
  integer,parameter   :: f_neta = 2400
  real(8),parameter   :: f_eta0 = -60.0d0, f_deta = 0.05d0
  real(8)             :: f_tab(0:f_neta,3)
  logical             :: f_built = .false.

  ! density floor to keep the carrier heat capacity well defined.  Set to the
  ! reference's intrinsic seed density Ns ~ 5.4e-15 bohr^-3 (~3.4e10 cm^-3) so the
  ! carrier population can reach the reference's near-transparent regime; this is
  ! safe only because the carrier temperatures use the exponential integrator
  ! (the old RK4 step blew up at such low Ne -- see ttm3_step_cell).
  real(8),parameter :: N_floor = 5.4d-15   ! [1/bohr^3] (~3.4e10 cm^-3, reference seed)
  real(8),parameter :: T_clamp = 31.67d0   ! [a.u.] (~1e7 K) numerical temperature cap
  real(8),parameter :: pi_ = 3.14159265358979323846d0

  character(12) :: ttm3_file = 'ttm3.inp_3tm'
  logical :: DISPLAY=.false.

contains

  !---------------------------------------------------------------------------
  subroutine init_ttm3_parameters( dt_em )
    use communication, only: comm_get_globalinfo, comm_is_root, comm_bcast
    implicit none
    real(8), intent(in) :: dt_em
    integer, parameter :: unit=1112
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

! Convert to atomic units
    tp3%Egap = tp3%Egap / hartree_ev
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

    hgs(:)      = hgs_in(:)
    is_array(:) = is_a(:)
    ie_array(:) = ie(:)

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
       if( imedia(ix,iy,iz)/=0 .and. &
           imedia(ix+1,iy,iz)==imedia(ix,iy,iz) .and. imedia(ix-1,iy,iz)==imedia(ix,iy,iz) .and. &
           imedia(ix,iy+1,iz)==imedia(ix,iy,iz) .and. imedia(ix,iy-1,iz)==imedia(ix,iy,iz) .and. &
           imedia(ix,iy,iz+1)==imedia(ix,iy,iz) .and. imedia(ix,iy,iz-1)==imedia(ix,iy,iz) ) m = m + 1
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
       if( imedia(ix,iy,iz)/=0 .and. &
           imedia(ix+1,iy,iz)==imedia(ix,iy,iz) .and. imedia(ix-1,iy,iz)==imedia(ix,iy,iz) .and. &
           imedia(ix,iy+1,iz)==imedia(ix,iy,iz) .and. imedia(ix,iy-1,iz)==imedia(ix,iy,iz) .and. &
           imedia(ix,iy,iz+1)==imedia(ix,iy,iz) .and. imedia(ix,iy,iz-1)==imedia(ix,iy,iz) )then
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
  end subroutine init_ttm3_alloc

  !---------------------------------------------------------------------------
  ! Build the Fermi-Dirac integral table F_j(eta), j = -1/2, 1/2, 3/2, once.
  ! F_j(eta) = (2/Gamma(j+1)) * Integral_0^inf u^(2j+1)/(exp(u^2-eta)+1) du
  ! (the substitution t=u^2 removes the t^(-1/2) endpoint singularity).
  subroutine ttm3_build_fermi_table()
    implicit none
    integer,parameter :: nq=4000
    integer :: k,q,jj
    real(8) :: eta,umax,du,u,integ,s,gam(3),jval(3)
    jval = (/ -0.5d0, 0.5d0, 1.5d0 /)
    gam  = (/ 1.7724538509055160d0, 0.8862269254527580d0, 1.3293403881791370d0 /) ! Gamma(j+1)
    do k=0,f_neta
       eta = f_eta0 + dble(k)*f_deta
       umax = sqrt( max(eta,0.0d0) + 50.0d0 )
       du = umax/dble(nq)
       do jj=1,3
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

  ! F_j(eta) by linear interpolation of the table (jidx: 1=-1/2, 2=1/2, 3=3/2).
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

  ! Reduced chemical potential eta from (T,mass,N) by bisecting N = Neff*F_{1/2}(eta),
  ! Neff = 2*(mass*T/(2*pi))^(3/2)  (atomic units).
  pure function ttm3_chem_pot( T, mc, N ) result( eta )
    implicit none
    real(8),intent(in) :: T,mc,N
    real(8) :: eta,Neff,lo,hi,mid,f
    integer :: it
    Neff = 2.0d0*( mc*max(T,1.0d-12)/(2.0d0*pi_) )**1.5d0
    lo = f_eta0; hi = f_eta0 + dble(f_neta)*f_deta
    do it=1,60
       mid = 0.5d0*(lo+hi)
       f = N - Neff*ttm3_fermi(2,mid)
       if(f > 0.0d0)then; lo=mid; else; hi=mid; end if
    end do
    eta = 0.5d0*(lo+hi)
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

    Ce = ttm3_heat_capacity( Te_, tp3%mu_e, Ne_+N_floor )
    Ch = ttm3_heat_capacity( Th_, tp3%mu_h, Nh_+N_floor )
    inv_tau = 1.0d0/tp3%tau
    geff = gen_ * max( 0.0d0, 1.0d0 - Ne_/tp3%N0 )
    R    = ( tp3%A_e*Ne_ + tp3%A_h*Nh_ )*Ne_*Nh_
    red_e = tp3%mu_h/(tp3%mu_e+tp3%mu_h)
    red_h = tp3%mu_e/(tp3%mu_e+tp3%mu_h)
    q_heat = max( source_ - geff*tp3%Egap, 0.0d0 )

    ! lattice: relaxation heating from the carriers (non-stiff, Cl large), using
    ! the start-of-step carrier temperatures (energy-conserving relaxation channel)
    dTl = ( Ce*(Te_-Tl_) + Ch*(Th_-Tl_) )*inv_tau / tp3%Cl

    ! carrier temperatures: exponential integrator (stiff relaxation + dilution)
    rate_e = inv_tau + geff/(Ne_+N_floor)
    rate_h = inv_tau + geff/(Nh_+N_floor)
    Te_star = Tl_ + (red_e*q_heat/Ce)/rate_e
    Th_star = Tl_ + (red_h*q_heat/Ch)/rate_h
    Te_ = Te_star + (Te_-Te_star)*exp( -rate_e*dt_in )
    Th_ = Th_star + (Th_-Th_star)*exp( -rate_h*dt_in )

    ! carrier densities and lattice: explicit Euler (non-stiff)
    Ne_ = max( Ne_ + dt_in*(geff - R), 0.0d0 )
    Nh_ = max( Nh_ + dt_in*(geff - R), 0.0d0 )
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

    if( tp3%Ddiff<=0.0d0 .and. tp3%kappa_e<=0.0d0 .and. tp3%kappa_l<=0.0d0 ) return

    cx=1.0d0/hgs(1)**2; cy=1.0d0/hgs(2)**2; cz=1.0d0/hgs(3)**2
    ! explicit FTCS diffusion stability: the update u += dt*D*lap with the 7-point
    ! Laplacian amplifies unless dt*D*2*(cx+cy+cz) <= 1.  Cap every diffusivity so
    ! this factor stays at 0.5 (safety).  [The previous 0.4*min(h)^2/dt cap was the
    ! 1-D form; in 3-D it gives dt*D*6/h^2 = 2.4 > 1, i.e. unstable -> Te blow-up.]
    Dmax = 0.5d0/( 2.0d0*dt*(cx+cy+cz) )
    Dd   = min(tp3%Ddiff, Dmax)
    call update_overlap_real8(srg, rg, Ne); call update_overlap_real8(srg, rg, Nh)
    call update_overlap_real8(srg, rg, Te); call update_overlap_real8(srg, rg, Th)
    call update_overlap_real8(srg, rg, Tl)

!$omp parallel do private(m,ix,iy,iz,Ce,Ch,lNe,lNh,lTe,lTh,lTl)
    do m=1,nmedia_myrnk
       ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
       lNe=(Ne(ix+1,iy,iz)-2*Ne(ix,iy,iz)+Ne(ix-1,iy,iz))*cx &
          +(Ne(ix,iy+1,iz)-2*Ne(ix,iy,iz)+Ne(ix,iy-1,iz))*cy &
          +(Ne(ix,iy,iz+1)-2*Ne(ix,iy,iz)+Ne(ix,iy,iz-1))*cz
       lNh=(Nh(ix+1,iy,iz)-2*Nh(ix,iy,iz)+Nh(ix-1,iy,iz))*cx &
          +(Nh(ix,iy+1,iz)-2*Nh(ix,iy,iz)+Nh(ix,iy-1,iz))*cy &
          +(Nh(ix,iy,iz+1)-2*Nh(ix,iy,iz)+Nh(ix,iy,iz-1))*cz
       lTe=(Te(ix+1,iy,iz)-2*Te(ix,iy,iz)+Te(ix-1,iy,iz))*cx &
          +(Te(ix,iy+1,iz)-2*Te(ix,iy,iz)+Te(ix,iy-1,iz))*cy &
          +(Te(ix,iy,iz+1)-2*Te(ix,iy,iz)+Te(ix,iy,iz-1))*cz
       lTh=(Th(ix+1,iy,iz)-2*Th(ix,iy,iz)+Th(ix-1,iy,iz))*cx &
          +(Th(ix,iy+1,iz)-2*Th(ix,iy,iz)+Th(ix,iy-1,iz))*cy &
          +(Th(ix,iy,iz+1)-2*Th(ix,iy,iz)+Th(ix,iy,iz-1))*cz
       lTl=(Tl(ix+1,iy,iz)-2*Tl(ix,iy,iz)+Tl(ix-1,iy,iz))*cx &
          +(Tl(ix,iy+1,iz)-2*Tl(ix,iy,iz)+Tl(ix,iy-1,iz))*cy &
          +(Tl(ix,iy,iz+1)-2*Tl(ix,iy,iz)+Tl(ix,iy,iz-1))*cz
       Ce=1.5d0*(Ne(ix,iy,iz)+N_floor); Ch=1.5d0*(Nh(ix,iy,iz)+N_floor)
       rhs_ne(ix,iy,iz)=dt*Dd*lNe
       rhs_nh(ix,iy,iz)=dt*Dd*lNh
       rhs_te(ix,iy,iz)=dt*min(tp3%kappa_e/Ce,Dmax)*lTe
       rhs_th(ix,iy,iz)=dt*min(tp3%kappa_e/Ch,Dmax)*lTh
       rhs_tl(ix,iy,iz)=dt*min(tp3%kappa_l/tp3%Cl,Dmax)*lTl
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

  !---------------------------------------------------------------------------
  ! Stage 2: carrier-dependent permittivity (recommended Drude model).
  ! eps(omega) = eps_bg - omega_p^2/omega^2,  omega_p^2 = 4*pi*(Ne/m_e + Nh/m_h);
  ! the Drude damping (rate 1/tau) gives the conductivity sigma. All atomic units.
  pure subroutine ttm3_permittivity( Ne_, Nh_, omega, eps_re, sig )
    implicit none
    real(8),intent(in)  :: Ne_, Nh_, omega
    real(8),intent(out) :: eps_re, sig
    real(8) :: wp2
    wp2    = 4.0d0*pi_*( Ne_/tp3%mu_e + Nh_/tp3%mu_h )
    eps_re = tp3%eps_bg - wp2/(omega*omega)
    sig    = wp2/(4.0d0*pi_)*(1.0d0/tp3%tau)/(omega*omega)
  end subroutine ttm3_permittivity

  !---------------------------------------------------------------------------
  ! Stage 2 back-action: carrier-dependent permittivity/conductivity at a cell
  ! (drives the FDTD media coefficients; single-frequency static eps at omega,
  ! the same approximation the reference uses).
  subroutine ttm3_eps_sig( ix, iy, iz, omega, eps_re, sig )
    implicit none
    integer,intent(in)  :: ix,iy,iz
    real(8),intent(in)  :: omega
    real(8),intent(out) :: eps_re, sig
    call ttm3_permittivity( Ne(ix,iy,iz), Nh(ix,iy,iz), omega, eps_re, sig )
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
  ! Stage 5: band-gap renormalisation (recommended form).
  ! Eg(Ne) = Egap - dgap_c * Ne^(1/3)  (band-gap shrinkage with carrier density).
  pure function ttm3_gap( Ne_ ) result( Eg )
    implicit none
    real(8),intent(in) :: Ne_
    real(8) :: Eg
    Eg = tp3%Egap - tp3%dgap_c*( max(Ne_,0.0d0) )**(1.0d0/3.0d0)
  end function ttm3_gap

end module ttm3
