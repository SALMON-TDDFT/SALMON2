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
  public :: init_ttm3_grid
  public :: init_ttm3_alloc
  public :: ttm3_main
  public :: ttm3_step_cell
  public :: ttm3_get_state

  logical,public :: use_ttm3=.false.

  ! material parameters
  type ttm3_param
     real(8) :: Egap     ! band gap                        [Hartree]
     real(8) :: mu_e     ! electron effective mass ratio    [-]
     real(8) :: mu_h     ! hole effective mass ratio        [-]
     real(8) :: A_e      ! electron Auger coefficient        [bohr^6 / a.u.time]
     real(8) :: A_h      ! hole Auger coefficient            [bohr^6 / a.u.time]
     real(8) :: tau      ! carrier-lattice relaxation time   [a.u.time]
     real(8) :: Cl       ! lattice heat capacity             [a.u. / (volume * temperature)]
     real(8) :: Tini     ! initial temperature (Te=Th=Tl)    [a.u.]
  end type ttm3_param

  integer,allocatable :: ijk_media_whole(:,:)
  integer,allocatable :: ijk_media_myrnk(:,:)
  integer :: nmedia_myrnk
  integer :: is_array(3), ie_array(3)
  integer :: comm
  real(8) :: hgs(3)
  real(8) :: dt

  type(ttm3_param) :: tp3

  ! state and Runge-Kutta stage buffers
  real(8),allocatable :: Te(:,:,:), Th(:,:,:), Tl(:,:,:)
  real(8),allocatable :: Ne(:,:,:), Nh(:,:,:)

  ! density floor to keep the carrier heat capacity well defined
  real(8),parameter :: N_floor = 1.0d-24   ! [1/bohr^3]

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
       close(unit)
       write(*,*) "Egap[eV]    =",tp3%Egap
       write(*,*) "mu_e        =",tp3%mu_e
       write(*,*) "mu_h        =",tp3%mu_h
       write(*,*) "A_e[cm6/s]  =",tp3%A_e
       write(*,*) "A_h[cm6/s]  =",tp3%A_h
       write(*,*) "tau[fs]     =",tp3%tau
       write(*,*) "Cl[J/m3K]   =",tp3%Cl
       write(*,*) "Tini[K]     =",tp3%Tini
    end if

    call comm_bcast(tp3%Egap,comm,0)
    call comm_bcast(tp3%mu_e,comm,0)
    call comm_bcast(tp3%mu_h,comm,0)
    call comm_bcast(tp3%A_e ,comm,0)
    call comm_bcast(tp3%A_h ,comm,0)
    call comm_bcast(tp3%tau ,comm,0)
    call comm_bcast(tp3%Cl  ,comm,0)
    call comm_bcast(tp3%Tini,comm,0)

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
    ! mu_e, mu_h are dimensionless

    if ( DISPLAY ) write(*,'(a60)') repeat("-",30)//" init_ttm3_parameters(end  )"

  end subroutine init_ttm3_parameters

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

    Te = tp3%Tini ; Th = tp3%Tini ; Tl = tp3%Tini
    Ne = N_floor  ; Nh = N_floor
  end subroutine init_ttm3_alloc

  !---------------------------------------------------------------------------
  ! Local rates for one cell (the right-hand side of the five ODEs).
  pure subroutine ttm3_rates( Te_,Th_,Tl_,Ne_,Nh_, source_, gen_, &
                              dTe,dTh,dTl,dNe,dNh )
    implicit none
    real(8),intent(in)  :: Te_,Th_,Tl_,Ne_,Nh_, source_, gen_
    real(8),intent(out) :: dTe,dTh,dTl,dNe,dNh
    real(8) :: Ce,Ch,R,inv_tau

    ! classical carrier heat capacities (floored to stay well defined)
    Ce = 1.5d0*(Ne_ + N_floor)
    Ch = 1.5d0*(Nh_ + N_floor)
    inv_tau = 1.0d0/tp3%tau

    ! standard Auger recombination (removes electron-hole pairs)
    R = ( tp3%A_e*Ne_ + tp3%A_h*Nh_ )*Ne_*Nh_

    ! carrier densities: generation - recombination (charge-neutral)
    dNe = gen_ - R
    dNh = gen_ - R

    ! temperatures: relaxation toward the lattice + heating source (split e/h).
    ! The lattice term is energy-conserving: the energy lost by the carriers
    ! (Ce*(Te-Tl) + Ch*(Th-Tl))/tau is deposited into the lattice / Cl.
    dTe = -(Te_-Tl_)*inv_tau + 0.5d0*source_/Ce
    dTh = -(Th_-Tl_)*inv_tau + 0.5d0*source_/Ch
    dTl = ( Ce*(Te_-Tl_) + Ch*(Th_-Tl_) )*inv_tau / tp3%Cl
  end subroutine ttm3_rates

  !---------------------------------------------------------------------------
  ! One classical Runge-Kutta (RK4) step for a single cell.  Pure scalar core
  ! with no grid/MPI dependence, so it can be exercised standalone.
  subroutine ttm3_step_cell( Te_,Th_,Tl_,Ne_,Nh_, source_, gen_, dt_in )
    implicit none
    real(8),intent(inout) :: Te_,Th_,Tl_,Ne_,Nh_
    real(8),intent(in)    :: source_, gen_, dt_in
    real(8) :: k1(5),k2(5),k3(5),k4(5)
    real(8) :: y(5), yt(5)

    y = (/Te_,Th_,Tl_,Ne_,Nh_/)

    call ttm3_rates( y(1),y(2),y(3),y(4),y(5), source_,gen_, k1(1),k1(2),k1(3),k1(4),k1(5) )
    yt = y + 0.5d0*dt_in*k1
    call ttm3_rates( yt(1),yt(2),yt(3),yt(4),yt(5), source_,gen_, k2(1),k2(2),k2(3),k2(4),k2(5) )
    yt = y + 0.5d0*dt_in*k2
    call ttm3_rates( yt(1),yt(2),yt(3),yt(4),yt(5), source_,gen_, k3(1),k3(2),k3(3),k3(4),k3(5) )
    yt = y + dt_in*k3
    call ttm3_rates( yt(1),yt(2),yt(3),yt(4),yt(5), source_,gen_, k4(1),k4(2),k4(3),k4(4),k4(5) )

    y = y + (dt_in/6.0d0)*( k1 + 2.0d0*k2 + 2.0d0*k3 + k4 )

    Te_=y(1); Th_=y(2); Tl_=y(3); Ne_=max(y(4),0.0d0); Nh_=max(y(5),0.0d0)
  end subroutine ttm3_step_cell

  !---------------------------------------------------------------------------
  ! Advance every medium cell by one time step (transport OFF).
  subroutine ttm3_main( srg, rg, source, gen )
    use structures, only: s_rgrid, s_sendrecv_grid
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg     ! reserved for the transport stage
    type(s_rgrid),         intent(in)    :: rg
    real(8),intent(in) :: source(rg%is_array(1):,rg%is_array(2):,rg%is_array(3):)
    real(8),intent(in) :: gen   (rg%is_array(1):,rg%is_array(2):,rg%is_array(3):)
    integer :: m,ix,iy,iz

!$omp parallel do private(m,ix,iy,iz)
    do m=1,nmedia_myrnk
       ix=ijk_media_myrnk(1,m); iy=ijk_media_myrnk(2,m); iz=ijk_media_myrnk(3,m)
       call ttm3_step_cell( Te(ix,iy,iz),Th(ix,iy,iz),Tl(ix,iy,iz), &
                            Ne(ix,iy,iz),Nh(ix,iy,iz), &
                            source(ix,iy,iz), gen(ix,iy,iz), dt )
    end do
!$omp end parallel do
  end subroutine ttm3_main

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

end module ttm3
