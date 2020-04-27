!
!  Copyright 2019-2020 SALMON developers
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
module em_field
  implicit none

contains

!===================================================================================================================================

subroutine calc_emfields(itt,nspin,curr_in,rt)
  use structures, only : s_rt
  use math_constants, only : pi
  use salmon_global, only : dt,trans_longi
  implicit none
  integer   ,intent(in)    :: itt,nspin
  real(8)   ,intent(in)    :: curr_in(3,nspin)
  type(s_rt),intent(inout) :: rt

  if(nspin==1) then
    rt%curr(1:3,itt) = curr_in(1:3,1)
  else if(nspin==2) then
    rt%curr(1:3,itt) = curr_in(1:3,1) + curr_in(1:3,2)
  end if

  if(trans_longi=="lo")then
    rt%Ac_ind(:,itt+1) = 2d0*rt%Ac_ind(:,itt) -rt%Ac_ind(:,itt-1) -4d0*Pi*rt%curr(:,itt)*dt**2
  else if(trans_longi=="tr")then
    rt%Ac_ind(:,itt+1) = 0d0
  end if

  rt%Ac_tot(:,itt+1) = rt%Ac_ext(:,itt+1) + rt%Ac_ind(:,itt+1)

  rt%E_ext(:,itt) = -(rt%Ac_ext(:,itt+1) - rt%Ac_ext(:,itt-1))/(2d0*dt)
  rt%E_ind(:,itt) = -(rt%Ac_ind(:,itt+1) - rt%Ac_ind(:,itt-1))/(2d0*dt)
  rt%E_tot(:,itt) = -(rt%Ac_tot(:,itt+1) - rt%Ac_tot(:,itt-1))/(2d0*dt)

end subroutine calc_emfields

!===================================================================================================================================

Subroutine calc_Ac_ext(t,Ac_ext)
  real(8),intent(in) :: t
  real(8)            :: Ac_ext(1:3)
  real(8)            :: Ac_ext_t(1:3, 0:0)

  call calc_Ac_ext_t(t, 0d0, 0, 0, Ac_ext_t)
  Ac_ext(1:3) = Ac_ext_t(1:3, 0)
  return
end subroutine calc_Ac_ext

!===================================================================================================================================

Subroutine calc_Ac_ext_t(t0, delta_t, is, ie, Ac_ext_t)
  use math_constants,only: zi,pi
  use salmon_global, only: I_wcm2_1,I_wcm2_2,E_amplitude1,E_amplitude2,ae_shape1,ae_shape2, &
                         & epdir_re1,epdir_re2,epdir_im1,epdir_im2,tw1,tw2,t1_start,omega1,omega2, &
                         & phi_CEP1,phi_CEP2,T1_T2,e_impulse,file_input1
  implicit none
  real(8),intent(in) :: t0
  real(8),intent(in) :: delta_t
  integer,intent(in) :: is
  integer,intent(in) :: ie
  real(8),intent(out) :: Ac_ext_t(1:3, is:ie)
  !
  integer :: i,npower
  real(8) :: t,f0_1,f0_2,tt,T1_T2_tmp
  
  Ac_ext_t(:,:) = 0d0
  if(t < 0d0) return

  if(I_wcm2_1 < 0d0)then
    f0_1 = E_amplitude1
  else
    f0_1=5.338d-9*sqrt(I_wcm2_1)      ! electric field in a.u.
  end if
  if(I_wcm2_2 < 0d0)then
    f0_2 = E_amplitude2
  else
    f0_2=5.338d-9*sqrt(I_wcm2_2)      ! electric field in a.u.
  end if
  
  T1_T2_tmp = T1_T2

  select case(AE_shape1)
  case('impulse')
    do i=is, ie
      t=t0+i*delta_t
      if (0d0 <= t) then
        Ac_ext_t(1,i) = epdir_re1(1)*e_impulse
        Ac_ext_t(2,i) = epdir_re1(2)*e_impulse
        Ac_ext_t(3,i) = epdir_re1(3)*e_impulse
      end if
    end do
    
  case('Acos2','Acos3','Acos4','Acos6','Acos8')
  
    select case(ae_shape1)
    case('Acos2'); npower = 2
    case('Acos3'); npower = 3
    case('Acos4'); npower = 4
    case('Acos6'); npower = 6
    case('Acos8'); npower = 8
    case default
      stop 'Error in init_rt.f90'
    end select

    do i=is, ie
      t=t0+i*delta_t
      tt = t - 0.5d0*tw1 - t1_start
      if (abs(tt)<0.5d0*tw1) then
        Ac_ext_t(:,i) = -f0_1/omega1*(cos(pi*tt/tw1))**npower &
          *aimag( (epdir_re1(:) + zI*epdir_im1(:)) &
          *exp(zI*(omega1*tt+phi_CEP1*2d0*pi))  &
          )
      end if
    end do
    T1_T2_tmp = T1_T2 + t1_start

  case('Ecos2')
  
    if(phi_CEP1 /= 0.75d0)then
      stop "Error: phi_cep1 should be 0.75 when ae_shape1 is 'Ecos2'."
    end if
    if(sum(abs(epdir_im1(:)))>1.0d-8)then
      stop "Error: ae_shape1 should be 'Acos2' when epdir_im1 is used."
    end if
    do i=is, ie
      t=t0+i*delta_t
      tt = t - 0.5d0*tw1 - t1_start
      if (abs(tt)<0.5d0*tw1) then
        Ac_ext_t(:,i) = -epdir_re1(:)*f0_1/(8d0*pi**2*omega1 - 2d0*tw1**2*omega1**3) &
          *( &
          (-4d0*pi**2+tw1**2*omega1**2 + tw1**2*omega1**2*cos(2d0*pi*tt/tw1))*cos(omega1*tt) &
          +2d0*pi*(2d0*pi*cos(tw1*omega1/2d0) &
          +tw1*omega1*sin(2d0*pi*tt/tw1)*sin(omega1*tt)))
      end if
    end do
    T1_T2_tmp = T1_T2 + t1_start

  case('Esin2sin')
  
    stop "Esin2sin is not implemented"
    
  case('Asin2cos')
  
    ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
    ! pump laser
    do i=is, ie
      t=t0+i*delta_t
      tt = t
      if (tt<tw1) then
        Ac_ext_t(:,i) = -epdir_re1(:)*f0_1/omega1*(sin(pi*tt/tw1))**2*cos(omega1*tt+phi_CEP1*2d0*pi)
      end if
    end do
    
  case('input')
     
    call add_Ac_from_file(trim(file_input1), t0, delta_t, is, ie, Ac_ext_t)
    !  stop "ae_shape1='input' is not implemented"
    
  case('Asin2_cw')
  
    ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
    ! pump laser
    do i=is, ie
      t=t0+i*delta_t
      tt = t
      if (tt<tw1*0.5d0) then
        Ac_ext_t(:,i) = -Epdir_re1(:)*f0_1/omega1*(sin(pi*tt/tw1))**2*cos(omega1*tt+phi_CEP1*2d0*pi)
      else
        Ac_ext_t(:,i) = -Epdir_re1(:)*f0_1/omega1*cos(omega1*tt+phi_CEP1*2d0*pi)
      end if
    end do
    
  case('none')
  
    Ac_ext_t(:,:) = 0d0
    
  case default
  
    stop "Invalid pulse_shape_1 parameter!"
    
  end select

! Probe
  select case(ae_shape2)
  case('impulse')

    do i=is, ie
      t=t0+i*delta_t
      tt = t - 0.5d0*tw1 - T1_T2_tmp
      if(tt > 0d0)then
        Ac_ext_t(1,i) = Ac_ext_t(1,i) + epdir_re2(1)*e_impulse
        Ac_ext_t(2,i) = Ac_ext_t(2,i) + epdir_re2(2)*e_impulse
        Ac_ext_t(3,i) = Ac_ext_t(3,i) + epdir_re2(3)*e_impulse
      end if
    end do
    
  case('Acos2','Acos3','Acos4','Acos6','Acos8')
  
    select case(ae_shape2)
    case('Acos2'); npower = 2
    case('Acos3'); npower = 3
    case('Acos4'); npower = 4
    case('Acos6'); npower = 6
    case('Acos8'); npower = 8
    case default
      stop 'Error in init_Ac.f90'
    end select

    do i=is, ie
      t=t0+i*delta_t
      tt = t - 0.5d0*tw1 - T1_T2_tmp
      if (abs(tt)<0.5d0*tw2) then
        Ac_ext_t(:,i)=Ac_ext_t(:,i) &
          -f0_2/omega2*(cos(pi*tt/tw2))**npower &
          *aimag( (epdir_re2(:) + zI*epdir_im2(:)) &
          *exp(zI*(omega2*tt+phi_CEP2*2d0*pi))  &
          )
      end if
    end do

  case('Ecos2')
  
    if(phi_CEP2 /= 0.75d0)then
      stop "Error: phi_cep2 should be 0.75 when ae_shape2 is 'Ecos2'."
    end if
    if(sum(abs(epdir_im2(:)))>1.0d-8)then
      stop "Error: ae_shape2 should be 'Acos2' when epdir_im2 is used."
    end if

    do i=is, ie
      t=t0+i*delta_t
      tt = t - 0.5d0*tw1 - T1_T2_tmp
      if (abs(tt)<0.5d0*tw2) then
        Ac_ext_t(:,i)=Ac_ext_t(:,i) &
          -epdir_re2(:)*f0_2/(8d0*pi**2*omega2 - 2d0*tw2**2*omega2**3) &
          *( &
          (-4d0*pi**2+tw2**2*omega2**2 + tw2**2*omega2**2*cos(2d0*pi*tt/tw2))*cos(omega2*tt) &
          +2d0*pi*(2d0*pi*cos(tw2*omega2/2d0) &
          +tw2*omega2*sin(2d0*pi*tt/tw2)*sin(omega2*tt)))
      end if
    end do

  case('Esin2sin')
      
    stop "Esin2sin is not implemented"
    
  case('Asin2cos')
  
      ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
    ! probe laser

    do i=is, ie
      t=t0+i*delta_t
      tt = t
      if ( (tt-T1_T2_tmp>0d0) .and. (tt-T1_T2_tmp<tw2) ) then
        Ac_ext_t(:,i) = Ac_ext_t(:,i) &
          &-Epdir_re2(:)*f0_2/omega2*(sin(pi*(tt-T1_T2_tmp)/tw2))**2*cos(omega2*(tt-T1_T2_tmp)+phi_CEP2*2d0*pi)
      endif
    end do

  case('input')
  case('Asin2_cw')
  case('none')
  case default
  
    stop "Invalid pulse_shape_2 parameter!"
    
  end select

  return
End Subroutine calc_Ac_ext_t

!===================================================================================================================================

Subroutine calc_E_ext(ipulse,t,E_ext,amp_flag)
  use salmon_global,  only: ae_shape1,ae_shape2,tw1,tw2,omega1,omega2,&
                            phi_cep1,phi_cep2,E_amplitude1,E_amplitude2,t1_t2
  use math_constants, only: Pi
  implicit none
  integer,intent(in)      :: ipulse  ! 1: first pulse, 2: second pulse
  real(8),intent(in)      :: t
  real(8),intent(out)     :: E_ext
  character(1),intent(in) :: amp_flag
  real(8)                 :: alpha,beta,theta1,theta2,amp
  character(16)           :: tae_shape
  
  if(ipulse==1)then
    ! cos(theta1)**2
    theta1 = Pi/tw1*(t-0.5d0*tw1)
    alpha  = Pi/tw1
    ! cos(theta2)
    theta2 = omega1*(t-0.5d0*tw1)+phi_cep1*2d0*pi
    beta   = omega1
    !other
    tae_shape = ae_shape1
    amp       = E_amplitude1
  else if(ipulse==2)then
    ! cos(theta1)**2
    theta1 = Pi/tw2*(t-t1_t2-0.5d0*tw1)
    alpha  = Pi/tw2
    ! cos(theta2)
    theta2 = omega2*(t-t1_t2-0.5d0*tw1)+phi_cep2*2d0*pi
    beta   = omega2
    !other
    tae_shape = ae_shape2
    amp       = E_amplitude2
  end if
  
  if(tae_shape=='Ecos2')then
    E_ext = cos(theta1)**2*cos(theta2)
  else if(tae_shape=='Acos2')then
    E_ext = -( -alpha*sin(2d0*theta1)*cos(theta2)   &
               -beta*cos(theta1)**2*sin(theta2) )/beta
  end if
  
  if(amp_flag=='y') E_ext = amp*E_ext
  
  return
End Subroutine calc_E_ext

!===================================================================================================================================

subroutine add_Ac_from_file(filename, t0, delta_t, is, ie, Ac_ext_t)
  use read_rtdata_file, only: count_rows_from_rtdata_file, load_Ac_from_rtdata_file
  implicit none
  character(256), intent(in) :: filename
  real(8),intent(in) :: t0
  real(8),intent(in) :: delta_t
  integer,intent(in) :: is, ie
  real(8),intent(inout) :: Ac_ext_t(1:3, is:ie)
  
  integer :: n_dat, i
  real(8), allocatable :: t_dat(:), Ac_dat(:, :)
  real(8) :: t, dt_dat, Ac(1:3)

  n_dat = count_rows_from_rtdata_file(trim(filename))

  allocate(t_dat(1:n_dat), Ac_dat(1:3, 1:n_dat))
  call load_Ac_from_rtdata_file(trim(filename), n_dat, t_dat, Ac_dat)

  do i = is, ie
      call interp_Ac(t0 + delta_t * i, Ac)
      Ac_ext_t(1:3, i) = Ac_ext_t(1:3, i) + Ac(1:3)
  end do

  deallocate(t_dat, Ac_dat)

  return

contains

  subroutine interp_Ac(tt, Ac)
    ! Linear interpolation of `Ac` at time `tt`.
    ! If a value outside the range of the data is requested, 
    ! the value at the nearest end is returned.
    implicit none
    real(8), intent(in) :: tt
    real(8), intent(out) :: Ac(1:3)
    integer :: ii
    real(8) :: x, rii
    if (tt < t_dat(1)) then
      Ac(1:3) = Ac_dat(1:3, 1)
    else if (t_dat(n_dat) <= tt) then
      Ac(1:3) = Ac_dat(1:3, n_dat)
    else
      ! Linear estimation of row index: ii
      ii = int(floor((tt - t_dat(1)) / (t_dat(n_dat) - t_dat(1)) * (n_dat - 1))) + 1
      if (.not. (t_dat(ii) <= tt .and. tt < t_dat(ii+1))) then
        ! Bluteforce search of row index: ii
        do ii = 1, n_dat-1
          if ((t_dat(ii) <= tt) .and. (tt < t_dat(ii+1))) exit
        end do
      end if
      x = (tt - t_dat(ii)) / (t_dat(ii+1) - t_dat(ii))
      Ac(1:3) = (1d0 - x) * Ac_dat(1:3, ii) + x * Ac_dat(1:3, ii+1)
    end if
    return
  end subroutine interp_Ac


end subroutine add_Ac_from_file

!===================================================================================================================================

subroutine calcVbox(mg,lg,itt_t,system,rt,Vbox)
  use structures, only: s_rgrid, s_dft_system, s_rt, s_scalar
  use communication, only: comm_is_root
  use misc_routines, only: get_wtime
  use inputoutput
  implicit none
  type(s_rgrid),intent(in) :: mg,lg
  integer :: itt_t
  type(s_dft_system),intent(inout) :: system
  type(s_scalar),intent(inout)     :: Vbox
  type(s_rt),intent(in) :: rt
  integer :: ix,iy,iz,jj
  integer :: ix_sta_Vbox(3),ix_end_Vbox(3)
  integer :: ipulse
  real(8) :: env_trigon_1,env_trigon_2
  integer,parameter :: Nd = 4

  if(iperiodic==0)then
    if(yn_md=='y' .or. yn_out_rvf_rt=='y')then
      do jj=1,3
        if(lg%is(jj)==mg%is(jj))then
          ix_sta_Vbox(jj)=mg%is(jj)
        else
          ix_sta_Vbox(jj)=mg%is(jj)-Nd
        end if
        if(lg%ie(jj)==mg%ie(jj))then
          ix_end_Vbox(jj)=mg%ie(jj)
        else
          ix_end_Vbox(jj)=mg%ie(jj)+Nd
        end if
      end do
    else
      ix_sta_Vbox(1:3)=mg%is(1:3)
      ix_end_Vbox(1:3)=mg%ie(1:3)
    end if
  end if
 
  !$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    Vbox%f(ix,iy,iz)=0.d0
  end do
  end do
  end do

  system%vec_E_ext(1:3)=0.d0
   
  if(iperiodic==0)then
    if(ae_shape1 == 'impulse')then
      continue
    else
        if(dt*dble(itt_t) <= tw1)then
          ipulse=1
          call calc_E_ext(ipulse,dt*dble(itt_t),env_trigon_1,'y')
        !$OMP parallel do collapse(2) private(ix,iy,iz)
          do iz=ix_sta_Vbox(3),ix_end_Vbox(3)
          do iy=ix_sta_Vbox(2),ix_end_Vbox(2)
          do ix=ix_sta_Vbox(1),ix_end_Vbox(1)
            Vbox%f(ix,iy,iz)=Vbox%f(ix,iy,iz)  &
                             +( epdir_re1(1)*lg%coordinate(ix,1)   &
                               +epdir_re1(2)*lg%coordinate(iy,2)   &
                               +epdir_re1(3)*lg%coordinate(iz,3) )*env_trigon_1  &
                             +( epdir_im1(1)*lg%coordinate(ix,1)   &
                               +epdir_im1(2)*lg%coordinate(iy,2)   &
                               +epdir_im1(3)*lg%coordinate(iz,3) )*env_trigon_1
          end do
          end do
          end do
          system%vec_E_ext(1:3)=system%vec_E_ext(1:3)+E_amplitude1*epdir_re1(1:3)*env_trigon_1   &
                                                       +E_amplitude1*epdir_im1(1:3)*env_trigon_1
        end if
        if(abs(dt*dble(itt_t)-0.5d0*tw1-t1_t2) < 0.5d0*tw2)then
          ipulse=2
          call calc_E_ext(ipulse,dt*dble(itt_t),env_trigon_2,'y')
          !$OMP parallel do collapse(2) private(ix,iy,iz)
          do iz=ix_sta_Vbox(3),ix_end_Vbox(3)
          do iy=ix_sta_Vbox(2),ix_end_Vbox(2)
          do ix=ix_sta_Vbox(1),ix_end_Vbox(1)
            Vbox%f(ix,iy,iz)=Vbox%f(ix,iy,iz)   &
                             +( epdir_re2(1)*lg%coordinate(ix,1)   &
                               +epdir_re2(2)*lg%coordinate(iy,2)   &
                               +epdir_re2(3)*lg%coordinate(iz,3) )*env_trigon_2  &
                             +( epdir_im2(1)*lg%coordinate(ix,1)   &
                               +epdir_im2(2)*lg%coordinate(iy,2)   &
                               +epdir_im2(3)*lg%coordinate(iz,3) )*env_trigon_2
          end do
          end do
          end do
          system%vec_E_ext(1:3)=system%vec_E_ext(1:3)+E_amplitude2*epdir_re2(1:3)*env_trigon_2   &
                                                       +E_amplitude2*epdir_im2(1:3)*env_trigon_2
        end if
    end if
  end if
  system%vec_E(1:3)=system%vec_E_ext(1:3) 
  system%vec_Ac(1:3)=system%vec_Ac(1:3)-system%vec_E(1:3)*dt
  system%vec_Ac_ext(1:3)=system%vec_Ac(1:3) 
   
  if(num_dipole_source>=1)then
    if(dt*dble(itt_t) <= tw1)then
      ipulse=1
      call calc_E_ext(ipulse,dt*dble(itt_t),env_trigon_1,'n')
!$OMP parallel do collapse(2) private(ix,iy,iz)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        Vbox%f(ix,iy,iz)=Vbox%f(ix,iy,iz)+rt%vonf%f(ix,iy,iz)*env_trigon_1
      end do
      end do
      end do
    end if
    if(abs(dt*dble(itt_t)-0.5d0*tw1-t1_t2) < 0.5d0*tw2)then
      ipulse=2
      call calc_E_ext(ipulse,dt*dble(itt_t),env_trigon_2,'n')
!$OMP parallel do collapse(2) private(ix,iy,iz)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        Vbox%f(ix,iy,iz)=Vbox%f(ix,iy,iz)+rt%vonf%f(ix,iy,iz)*env_trigon_2
      end do
      end do
      end do
    end if
  end if

  return
end subroutine calcVbox

!===================================================================================================================================

subroutine set_vonf(mg,lg,Hgs,rt)
!$ use omp_lib
  use structures, only: s_rgrid,s_rt
  use communication, only: comm_is_root, comm_summation
  use salmon_global
  implicit none
  real(8),intent(in) :: Hgs(3)
  type(s_rgrid) :: mg,lg
  type(s_rt) :: rt
  integer :: i
  integer :: ix,iy,iz,iix,iiy,iiz
  integer :: max_icell
  real(8) :: rr

  if(iperiodic==0)then
    max_icell=0
  else if(iperiodic==3)then
    max_icell=2
  end if

  rt%vonf%f=0.d0

  do i=1,num_dipole_source
    do iiz=-max_icell,max_icell
    do iiy=-max_icell,max_icell
    do iix=-max_icell,max_icell
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        rr=sqrt((lg%coordinate(ix,1)-(cood_dipole_source(1,i)+dble(iix*lg%num(1))*Hgs(1)))**2 &
               +(lg%coordinate(iy,2)-(cood_dipole_source(2,i)+dble(iiy*lg%num(2))*Hgs(2)))**2 &
               +(lg%coordinate(iz,3)-(cood_dipole_source(3,i)+dble(iiz*lg%num(3))*Hgs(3)))**2)
        if(rr>rad_dipole_source)then
          rt%vonf%f(ix,iy,iz)=rt%vonf%f(ix,iy,iz)  &
                        -(vec_dipole_source(1,i)*(lg%coordinate(ix,1)-(cood_dipole_source(1,i)+dble(iix*lg%num(1))*Hgs(1))) &
                         +vec_dipole_source(2,i)*(lg%coordinate(iy,2)-(cood_dipole_source(2,i)+dble(iiy*lg%num(2))*Hgs(2))) &
                         +vec_dipole_source(3,i)*(lg%coordinate(iz,3)-(cood_dipole_source(3,i)+dble(iiz*lg%num(3))*Hgs(3))))/rr**3
        else
          rt%vonf%f(ix,iy,iz)=rt%vonf%f(ix,iy,iz)  &
                        -(vec_dipole_source(1,i)*(lg%coordinate(ix,1)-(cood_dipole_source(1,i)+dble(iix*lg%num(1))*Hgs(1))) &
                         +vec_dipole_source(2,i)*(lg%coordinate(iy,2)-(cood_dipole_source(2,i)+dble(iiy*lg%num(2))*Hgs(2))) &
                         +vec_dipole_source(3,i)*(lg%coordinate(iz,3)-(cood_dipole_source(3,i)+dble(iiz*lg%num(3))*Hgs(3))))&
                                                                                                        /rad_dipole_source**3
        end if
      end do
      end do
      end do
    end do
    end do
    end do
  end do

  return
end subroutine set_vonf

end module em_field
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
