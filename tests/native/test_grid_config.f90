program probe
 use structures,only:s_rgrid,s_rgrid_config
 use initialization_sub,only:init_grid_whole,init_grid_whole_config,grid_config_from_input
 use salmon_global,only:iperiodic,num_rgrid,dl,theory,al_em,dl_em,num_rgrid_em
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 type(s_rgrid)::a,b,reference
 type(s_rgrid_config)::ca,cb,bad
 integer::status,period,mode,variant
 real(8)::lengths(3),spacing(3),saved_nan
 do period=0,3,3
 do mode=1,2
 do variant=1,2
  lengths=[8d0,12d0,16d0];spacing=1d0
  if(variant==2)lengths=[9d0,11d0,13d0]
  iperiodic=period;num_rgrid=int(lengths);theory='tddft_response';dl=0d0
  if(mode==2)dl=spacing
  call grid_config_from_input(lengths,ca)
  call frozen_grid_whole(lengths,spacing,reference)
  ! Global input no longer chooses the explicit grid.
  iperiodic=3-period;num_rgrid=3;dl=7d0;theory='maxwell'
  call init_grid_whole_config(lengths,spacing,ca,a,status)
  if(status/=0)stop 1
  if(any(a%is/=reference%is).or.any(a%ie/=reference%ie).or.any(a%idx/=reference%idx))stop 2
  cb=ca;cb%points=[12,8,20];cb%lengths=[12d0,8d0,20d0]
  call init_grid_whole_config(cb%lengths,spacing,cb,b,status)
  if(status/=0.or.any(b%num/=cb%points))stop 3
  if(any(a%num/=ca%points))stop 4
  if(period==3)then
    if(any(a%is/=1).or.any(a%ie/=ca%points))stop 9
  else
    if(any(a%ie/=ca%points/2))stop 10
    if(variant==1.and.any(a%is/=1-ca%points/2))stop 11
    if(variant==2.and.any(a%is/=-ca%points/2))stop 12
  endif
  if(any(a%is_overlap/=a%is-4).or.any(a%ie_overlap/=a%ie+4))stop 13
  bad=ca;bad%points=1;bad%spacing=0d0
  call init_grid_whole_config(lengths,spacing,bad,a,status)
  if(status==0.or.any(a%idx/=reference%idx))stop 5
  saved_nan=ieee_value(0d0,ieee_quiet_nan)
  call init_grid_whole_config(lengths,[saved_nan,1d0,1d0],ca,a,status)
  if(status==0.or.any(a%ie/=reference%ie))stop 6
  call init_grid_whole_config([huge(1d0),12d0,16d0],spacing,ca,a,status)
  if(status==0.or.any(a%idx/=reference%idx))stop 8
 enddo
 enddo
 enddo
 ! Maxwell's consistency lengths are al_em, preserving the legacy rule.
 iperiodic=3;theory='maxwell';al_em=[8d0,12d0,16d0];dl_em=1d0;num_rgrid_em=[8,12,16]
 call grid_config_from_input(al_em,ca)
 call frozen_grid_whole(al_em,dl_em,reference)
 call init_grid_whole_config(al_em,dl_em,ca,a,status)
 if(status/=0.or.any(a%idz/=reference%idz))stop 7
 ! Freeze the original epsilon-boundary behavior independently of the new wrapper.
 do period=0,3,3
  lengths=[8d0-5d-11,12d0,16d0];spacing=1d0
  iperiodic=period;theory='tddft_response';dl=0d0;num_rgrid=[8,12,16]
  if(period==0)num_rgrid(1)=9
  call grid_config_from_input(lengths,ca)
  call frozen_grid_whole(lengths,spacing,reference)
  call init_grid_whole_config(lengths,spacing,ca,a,status)
  if(status/=0.or.any(a%is/=reference%is).or.any(a%idx/=reference%idx))stop 14
  if(period==3.and.a%num(1)/=8)stop 15
  if(period==0.and.a%num(1)/=9)stop 16
 enddo
 print *,'Explicit real-space grid configuration PASS'
contains
! Frozen original SALMON init_grid_whole, Apache-2.0, SALMON developers 2019-2020.
! Source: runtime/sources/hse-fast-reference/src/common/initialization.f90
! Original commit c9e8ce431f43eefd4a88998134adb4796f8ec837; renamed only, Nd declared locally.
subroutine frozen_grid_whole(rsize,hgs,lg)
  use structures, only: s_rgrid
  use salmon_global, only: iperiodic,dl,num_rgrid,theory,al_em,dl_em,num_rgrid_em
  implicit none
  integer,parameter::Nd=4
  real(8),intent(in) :: rsize(3),hgs(3)
  type(s_rgrid)      :: lg
  !
  real(8),parameter :: epsilon=1.d-10
  integer :: j

  lg%ndir = 3 ! high symmetry nonorthogonal lattice is not implemented
  lg%nd = Nd

  select case(iperiodic)
    case(0)
      lg%ie(:)=int((rsize(:)+epsilon)/2.d0/Hgs(:))
      do j=1,3
        if(mod(int(rsize(j)/Hgs(j)+1.d-12),2)==1)then
          lg%is(j)=-(int((rsize(j)+epsilon)/2.d0/Hgs(j)))
        else
          lg%is(j)=-(int((rsize(j)+epsilon)/2.d0/Hgs(j)))+1
        end if
      end do
    case(3)
      lg%is(:)=1
      lg%ie(:)=int((rsize(:)+epsilon)/Hgs(:))
  end select
  lg%num(:)=lg%ie(:)-lg%is(:)+1

  lg%is_overlap(1:3) = lg%is(1:3) - Nd
  lg%ie_overlap(1:3) = lg%ie(1:3) + Nd
  lg%is_array(1:3)   = lg%is_overlap(1:3)
  lg%ie_array(1:3)   = lg%ie_overlap(1:3)

  if ( allocated(lg%idx) ) deallocate(lg%idx)
  if ( allocated(lg%idy) ) deallocate(lg%idy)
  if ( allocated(lg%idz) ) deallocate(lg%idz)
  allocate(lg%idx(lg%is_overlap(1):lg%ie_overlap(1)) &
          ,lg%idy(lg%is_overlap(2):lg%ie_overlap(2)) &
          ,lg%idz(lg%is_overlap(3):lg%ie_overlap(3)))

  do j=lg%is_overlap(1),lg%ie_overlap(1)
    lg%idx(j) = j
  end do
  do j=lg%is_overlap(2),lg%ie_overlap(2)
    lg%idy(j) = j
  end do
  do j=lg%is_overlap(3),lg%ie_overlap(3)
    lg%idz(j) = j
  end do

  select case(theory)
  case('maxwell')
    if(sum(abs(dl_em)) <= 1d-12) then
      if( maxval(abs(num_rgrid_em-lg%num)) > 0) stop "error: num_rgrid_em /= lg%num"
    else
      if( maxval(abs((al_em/dl_em)-dble(lg%num))) > 1d-4 ) stop "error: abs((al_em/dl_em)-dble(lg%num)) is too large"
    end if
  case default
    if(sum(abs(dl)) <= 1d-12) then
      if( maxval(abs(num_rgrid-lg%num)) > 0) stop "error: num_rgrid /= lg%num"
    else
      if( maxval(abs((rsize/dl)-dble(lg%num))) > 1d-4 ) stop "error: abs((rsize/dl)-dble(lg%num)) is too large"
    end if
  end select

  return
end subroutine frozen_grid_whole
end program
