module ttm

  implicit none
  private
  public :: init_ttm_parameters
  public :: init_ttm_grid
  public :: init_ttm_alloc
  public :: ttm_penetration
  public :: ttm_main
  public :: ttm_get_temperatures

  logical,public :: use_ttm=.false.

  type ttm_param
     real(8) :: g
     real(8) :: Cl
     real(8) :: kappa_e0 ! Thermal conductivity for Te==Tl
     real(8) :: Ce_prime ! a constant depending on the material
     real(8) :: zeta_bal ! ballistic penetration depth
     real(8) :: Tini     ! Initial temperature (Te==Tl==Tini)
  end type ttm_param
  
  integer,allocatable :: ijk_media_whole(:,:)
  integer,allocatable :: ijk_media_myrnk(:,:)
  integer :: is_array(3), ie_array(3)
  integer :: lb_array(3), ub_array(3)
  integer :: comm
  real(8) :: hgs(3)
  real(8) :: dt

  type(ttm_param) :: tp

  real(8),allocatable :: Te(:,:,:), NABLA_Te(:,:,:)
  real(8),allocatable :: Tl(:,:,:)
  real(8),allocatable :: work(:,:,:), NABLA_work(:,:,:)
  real(8),allocatable :: rhs_e(:,:,:)
  real(8),allocatable :: rhs_l(:,:,:)

  character(7) :: ttm_file = 'ttm.inp_ttm'
  logical :: DISPLAY=.false.

contains

  subroutine init_ttm_parameters( dt_em )
    use communication, only: comm_get_globalinfo, comm_is_root, comm_bcast
    implicit none
    real(8), intent(in) :: dt_em
    integer, parameter :: unit=1111
    integer :: npid, nprocs
    logical :: flag
    real(8), parameter :: atomic_unit_of_length = 5.29177210903d-11 ! [m] (+/- 8.0e-21)
    real(8), parameter :: hartree_joule_relationship = 4.3597447222071d-18 ! [J] (+/- 8.5e-30)
    real(8), parameter :: hartree_kelvin_relationship = 3.1577502480407d5 ! [K] (+/- 6.1e-07)
    real(8), parameter :: atomic_unit_of_time = 2.4188843265857d-17 ! [s] (+/- 4.7e-29)

    call comm_get_globalinfo( comm, npid, nprocs )
    DISPLAY = comm_is_root(npid)

    if ( DISPLAY ) write(*,'(a60)') repeat("-",33)//" init_ttm_parameters(start)"

    inquire( FILE=ttm_file, EXIST=flag )
    if( .not.flag )then
       if( DISPLAY )then
          write(*,*) "No TTM file ("//ttm_file//")."
          write(*,*) "Two-Temperature Model is not used."
       end if
       return
    else
       if( DISPLAY )then
          write(*,*) "TTM file ("//ttm_file//") is found."
          write(*,*) "Calculation is performed with Two-Temperature Model."
       end if
       use_ttm = .true.
    end if

    dt = dt_em

! Input parameters
! Symbols and units are followed by the paper Phys.Rev.B84,033405(2011).
!    tp%g        = 24.5 ! J/m^3
!    tp%Cl       = 2.44 ! J/(m^3K)
!    tp%kappa_e0 = 235  ! J/(msK)
!    tp%Ce_prime = 135  ! J/(m^3K^2)
!    tp%zeta_bal = 46   ! nm
!    tp%Tini     = 80   ! K

    if( comm_is_root(npid) )then
       open(unit, file=ttm_file, status='old')
       read(unit,*) tp%g
       read(unit,*) tp%Cl
       read(unit,*) tp%kappa_e0
       read(unit,*) tp%Ce_prime
       read(unit,*) tp%zeta_bal
       read(unit,*) tp%Tini
       close(unit)
       write(*,*) "g       =",tp%g
       write(*,*) "Cl      =",tp%Cl
       write(*,*) "kappa_e0=",tp%kappa_e0
       write(*,*) "Ce'     =",tp%Ce_prime
       write(*,*) "zeta_bal=",tp%zeta_bal
       write(*,*) "Tini    =",tp%Tini
    end if

    call comm_bcast(tp%g       ,comm,0)
    call comm_bcast(tp%Cl      ,comm,0)
    call comm_bcast(tp%kappa_e0,comm,0)
    call comm_bcast(tp%Ce_prime,comm,0)
    call comm_bcast(tp%zeta_bal,comm,0)
    call comm_bcast(tp%Tini    ,comm,0)

! Convert to the atomic unit
!
    tp%g = tp%g /hartree_joule_relationship &
                *atomic_unit_of_length**3
    tp%Cl = tp%Cl /hartree_joule_relationship &
                  *atomic_unit_of_length**3 &
                  *hartree_kelvin_relationship
    tp%kappa_e0 = tp%kappa_e0 /hartree_joule_relationship &
                              *atomic_unit_of_length &
                              *hartree_kelvin_relationship &
                              *atomic_unit_of_time
    tp%Ce_prime = tp%Ce_prime /hartree_joule_relationship*atomic_unit_of_length**3 &
                              *hartree_kelvin_relationship**2
    tp%zeta_bal = tp%zeta_bal /(atomic_unit_of_length*1.0d9)
    tp%Tini = tp%Tini /hartree_kelvin_relationship

  end subroutine init_ttm_parameters

  subroutine init_ttm_grid( hgs_in, is_a, is, ie, imedia )
    use parallelization, only: nproc_id_global, nproc_size_global, nproc_group_global
    use communication, only: comm_summation
    implicit none
    real(8), intent(in) :: hgs_in(3)
    integer, intent(in) :: is_a(3), is(3), ie(3)
    integer, intent(in) :: imedia(is_a(1):,is_a(2):,is_a(3):)
    integer :: ii,ij,ix,iy,iz,i
    integer,allocatable :: ircnt(:),idisp(:),ijk_tmp(:,:)

    if ( DISPLAY ) write(*,'(a60)') repeat("-",39)//" init_ttm_grid(start)"
    hgs(1:3) = hgs_in(1:3)

    ii=count(imedia(is(1):ie(1),is(2):ie(2),is(3):ie(3))/=0)
    allocate( ijk_media_myrnk(3,ii) ); ijk_media_myrnk=0

    ii=0
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
       if ( imedia(ix,iy,iz) /= 0 ) then
          ii=ii+1
          ijk_media_myrnk(1,ii)=ix
          ijk_media_myrnk(2,ii)=iy
          ijk_media_myrnk(3,ii)=iz
       end if
    end do
    end do
    end do

! ---

    allocate( ircnt(0:nproc_size_global-1) ); ircnt=0
    allocate( idisp(0:nproc_size_global-1) ); idisp=0

    idisp(nproc_id_global) = ii
    call comm_summation( idisp, ircnt, nproc_size_global, nproc_group_global )

    idisp=0
    do i=0,nproc_size_global-1
       idisp(i) = sum(ircnt(0:i)) - ircnt(i)
    end do

    ij=sum(ircnt)
    allocate( ijk_media_whole(3,ij) ); ijk_media_whole=0
    allocate( ijk_tmp(3,ij) ); ijk_tmp=0

    do ii=1,ircnt(nproc_id_global)
       ij=idisp(nproc_id_global)+ii
       ijk_tmp(:,ij) = ijk_media_myrnk(:,ii)
    end do
    call comm_summation( ijk_tmp, ijk_media_whole, size(ijk_tmp), nproc_group_global )

! ---

    do i=1,3
       lb_array(i) = minval( ijk_media_whole(i,:) )
       ub_array(i) = maxval( ijk_media_whole(i,:) )
    end do

    do i=1,3
       is_array(i) = is_a(i)
       ie_array(i) = is_a(i) + size(imedia,i) - 1
    end do

  end subroutine init_ttm_grid


  subroutine init_ttm_alloc( srg, rg )
    use structures, only: s_rgrid, s_sendrecv_grid
    use sendrecv_grid, only: update_overlap_real8
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_rgrid), intent(in) :: rg
    integer :: i,ix,iy,iz

    if ( DISPLAY ) write(*,'(a60)') repeat("-",38)//" init_ttm_alloc(start)"

! ---

    allocate( Te(rg%is_array(1):rg%ie_array(1), &
                 rg%is_array(2):rg%ie_array(2), &
                 rg%is_array(3):rg%ie_array(3)) ); Te=0.0d0
    allocate( Tl(rg%is_array(1):rg%ie_array(1), &
                 rg%is_array(2):rg%ie_array(2), &
                 rg%is_array(3):rg%ie_array(3)) ); Tl=0.0d0

    do i = 1, size(ijk_media_myrnk,2)
       ix = ijk_media_myrnk(1,i)
       iy = ijk_media_myrnk(2,i)
       iz = ijk_media_myrnk(3,i)
       Te(ix,iy,iz) = tp%Tini
       Tl(ix,iy,iz) = tp%Tini
    end do

    call update_overlap_real8( srg, rg, Te )
    call update_overlap_real8( srg, rg, Tl )

    allocate( work(rg%is_array(1):rg%ie_array(1), &
                   rg%is_array(2):rg%ie_array(2), &
                   rg%is_array(3):rg%ie_array(3)) ); work=0.0d0
    allocate( NABLA_Te(rg%is_array(1):rg%ie_array(1), &
                       rg%is_array(2):rg%ie_array(2), &
                       rg%is_array(3):rg%ie_array(3)) ); NABLA_Te=0.0d0
    allocate( NABLA_work(rg%is_array(1):rg%ie_array(1), &
                         rg%is_array(2):rg%ie_array(2), &
                         rg%is_array(3):rg%ie_array(3)) ); NABLA_work=0.0d0
    allocate( rhs_e(rg%is_array(1):rg%ie_array(1), &
                    rg%is_array(2):rg%ie_array(2), &
                    rg%is_array(3):rg%ie_array(3)) ); rhs_e=0.0d0
    allocate( rhs_l(rg%is_array(1):rg%ie_array(1), &
                    rg%is_array(2):rg%ie_array(2), &
                    rg%is_array(3):rg%ie_array(3)) ); rhs_e=0.0d0

  end subroutine init_ttm_alloc


  subroutine ttm_main( srg, rg, source )
    use structures, only: s_rgrid, s_sendrecv_grid
    use sendrecv_grid, only: update_overlap_real8
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_rgrid), intent(in) :: rg
    real(8), intent(in) :: source(rg%is(1):,rg%is(2):,rg%is(3):)
    integer :: i,ii,ix,iy,iz
    real(8) :: const_e, const_l

    rhs_e(:,:,:) = 0.0d0
    rhs_l(:,:,:) = 0.0d0
    const_e = dt/tp%Ce_prime
    const_l = dt/tp%Cl*tp%g

    call update_overlap_real8( srg, rg, Te )

    do i = 1, 3
       call calc_nabla( Te, NABLA_Te, i )
       where( Tl /= 0.0d0 )
          work(:,:,:) = tp%kappa_e0*( Te(:,:,:)/Tl(:,:,:) ) * NABLA_Te(:,:,:)
       end where
       call update_overlap_real8( srg, rg, work )
       call calc_nabla( work, NABLA_work, i )
       rhs_e = rhs_e + NABLA_work
    end do

    do ii = 1, size(ijk_media_myrnk,2)
       ix = ijk_media_myrnk(1,ii)
       iy = ijk_media_myrnk(2,ii)
       iz = ijk_media_myrnk(3,ii)
       rhs_e(ix,iy,iz) = rhs_e(ix,iy,iz) &
            - tp%g*(Te(ix,iy,iz)-Tl(ix,iy,iz)) + source(ix,iy,iz)
       rhs_e(ix,iy,iz) = rhs_e(ix,iy,iz)*const_e/Te(ix,iy,iz)
       rhs_l(ix,iy,iz) = const_l*(Te(ix,iy,iz)-Tl(ix,iy,iz))
    end do

    do ii = 1, size(ijk_media_myrnk,2)
       ix = ijk_media_myrnk(1,ii)
       iy = ijk_media_myrnk(2,ii)
       iz = ijk_media_myrnk(3,ii)
       Te(ix,iy,iz) = Te(ix,iy,iz) + rhs_e(ix,iy,iz)
       Tl(ix,iy,iz) = Tl(ix,iy,iz) + rhs_l(ix,iy,iz)
    end do

!---
!   Ce(Te)*dTe(r,t)/dt = NABLA{ kappa_e[Te(r,t),Tl(r,t)]*NABLA[Te(r,t)] }
!                      - g*{Te(r,t)-Tl(r,t)} + S(r,t)
!
!   Cl*dTl/dt = g*{Te(r,t)-Tl(r,t)}
!
!   kappa_e[Te(r,t),Tl(r,t)] = kappa_e0*{Te(r,t)/Tl(r,t)}
!
!   Ce(Te) = Ce_prime * Te(r,t)
!
!   Initial condition
!    Qe(0,t)=Qe(L,t)=0
!    Te(r,0)=Tl(r,0)=80K
!---

  contains

    subroutine calc_nabla( f, NABLA_f, idir )
      implicit none
      real(8), intent(in) :: f(is_array(1):,is_array(2):,is_array(3):)
      real(8), intent(inout) :: NABLA_f(is_array(1):,is_array(2):,is_array(3):)
      integer, intent(in) :: idir
      real(8) :: c1,c2,c3
      integer :: ii,ix,iy,iz
      c1 = 0.5d0/hgs(1)
      c2 = 0.5d0/hgs(2)
      c3 = 0.5d0/hgs(3)
      select case(idir)
      case( 1 )
         do ii = 1, size(ijk_media_myrnk,2)
            ix = ijk_media_myrnk(1,ii)
            iy = ijk_media_myrnk(2,ii)
            iz = ijk_media_myrnk(3,ii)
            NABLA_f(ix,iy,iz) = ( f(ix+1,iy,iz) - f(ix-1,iy,iz) )*c1
         end do
      case( 2 )
         do ii = 1, size(ijk_media_myrnk,2)
            ix = ijk_media_myrnk(1,ii)
            iy = ijk_media_myrnk(2,ii)
            iz = ijk_media_myrnk(3,ii)
            NABLA_f(ix,iy,iz) = ( f(ix,iy+1,iz) - f(ix,iy-1,iz) )*c2
         end do
      case( 3 )
         do ii = 1, size(ijk_media_myrnk,2)
            ix = ijk_media_myrnk(1,ii)
            iy = ijk_media_myrnk(2,ii)
            iz = ijk_media_myrnk(3,ii)
            NABLA_f(ix,iy,iz) = ( f(ix,iy,iz+1) - f(ix,iy,iz-1) )*c3
         end do
      end select
    end subroutine calc_nabla

  end subroutine ttm_main


  subroutine ttm_penetration( is, f )
    use communication, only: comm_summation
    implicit none
    integer, intent(in) :: is(3)   ! lower bounds of f 
    real(8), intent(inout) :: f(is(1):,is(2):,is(3):)
    integer :: num_points_in_myrnk,num_points_in_whole
    integer :: ii,jj,ix,iy,iz,jx,jy,jz
    real(8) :: dx,dy,dz,r,rr,factor
    real(8) :: sum_f,sum_fnew,c,const,tmp
    real(8),allocatable :: fnew(:,:,:), work(:,:,:)

    const = 1.0d0/tp%zeta_bal

    num_points_in_myrnk = size( ijk_media_myrnk,2 )
    num_points_in_whole = size( ijk_media_whole,2 ) 

    tmp=0.0d0
    do ii = 1, num_points_in_myrnk
       ix = ijk_media_myrnk(1,ii)
       iy = ijk_media_myrnk(2,ii)
       iz = ijk_media_myrnk(3,ii)
       tmp = tmp + f(ix,iy,iz)
    end do
    call comm_summation( tmp, sum_f, comm )

    if ( sum_f == 0.0d0 ) return

    allocate( fnew(lb_array(1):ub_array(1),lb_array(2):ub_array(2),lb_array(3):ub_array(3)) )
    fnew=0.0d0
    allocate( work(lb_array(1):ub_array(1),lb_array(2):ub_array(2),lb_array(3):ub_array(3)) )
    work=0.0d0

    do ii = 1, num_points_in_myrnk
       ix = ijk_media_myrnk(1,ii)
       iy = ijk_media_myrnk(2,ii)
       iz = ijk_media_myrnk(3,ii)
       do jj = 1, num_points_in_whole
          jx = ijk_media_whole(1,jj)
          jy = ijk_media_whole(2,jj)
          jz = ijk_media_whole(3,jj)
          dx = (ix-jx)*hgs(1)
          dy = (iy-jy)*hgs(2)
          dz = (iz-jz)*hgs(3)
          rr = dx*dx + dy*dy + dz*dz
          r = sqrt(rr)
          factor=exp(-r*const)
          work(jx,jy,jz) = work(jx,jy,jz) + factor*f(ix,iy,iz)
       end do !jj
    end do !ii

    call comm_summation( work, fnew, size(fnew), comm )

    sum_fnew = sum(fnew)
    if ( sum_fnew == 0.0d0 ) return

    c = sum_f/sum_fnew

    f=0.0d0
    tmp=0.0d0
    do ii = 1, num_points_in_myrnk
       ix = ijk_media_myrnk(1,ii)
       iy = ijk_media_myrnk(2,ii)
       iz = ijk_media_myrnk(3,ii)
       f(ix,iy,iz) = c*fnew(ix,iy,iz)
       tmp = tmp + f(ix,iy,iz)
    end do
    call comm_summation( tmp, sum_fnew, comm )

  end subroutine ttm_penetration


  subroutine ttm_get_temperatures( a, t1, t2 )
    implicit none
    integer, intent(in) :: a(3)
    real(8), intent(inout) :: t1(a(1):,a(2):,a(3):), t2(a(1):,a(2):,a(3):)
    integer :: ii, ix, iy, iz
    do ii = 1, size(ijk_media_myrnk,2)
       ix = ijk_media_myrnk(1,ii)
       iy = ijk_media_myrnk(2,ii)
       iz = ijk_media_myrnk(3,ii)
       t1(ix,iy,iz)=Te(ix,iy,iz)
       t2(ix,iy,iz)=Tl(ix,iy,iz)
    end do
  end subroutine ttm_get_temperatures

end module ttm
