!
!  Copyright 2019 SALMON developers
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

Subroutine calc_Ac_ext(t,Ac_ext)
  use math_constants,only: zi,pi
  use salmon_global, only: I_wcm2_1,I_wcm2_2,E_amplitude1,E_amplitude2,ae_shape1,ae_shape2, &
                         & epdir_re1,epdir_re2,epdir_im1,epdir_im2,tw1,tw2,t1_start,omega1,omega2, &
                         & phi_CEP1,phi_CEP2,T1_T2,e_impulse
  implicit none
  real(8),intent(in) :: t
  real(8)            :: Ac_ext(3)
  !
  integer :: npower
  real(8) :: f0_1,f0_2,tt,T1_T2_tmp
  
  Ac_ext = 0d0
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
  
    Ac_ext(1) = epdir_re1(1)*e_impulse
    Ac_ext(2) = epdir_re1(2)*e_impulse
    Ac_ext(3) = epdir_re1(3)*e_impulse
    
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

    tt = t - 0.5d0*tw1 - t1_start
    if (abs(tt)<0.5d0*tw1) then
      Ac_ext(:) = -f0_1/omega1*(cos(pi*tt/tw1))**npower &
        *aimag( (epdir_re1(:) + zI*epdir_im1(:)) &
        *exp(zI*(omega1*tt+phi_CEP1*2d0*pi))  &
        )
    end if
    T1_T2_tmp = T1_T2 + t1_start

  case('Ecos2')
  
    if(phi_CEP1 /= 0.75d0)then
      stop "Error: phi_cep1 should be 0.75 when ae_shape1 is 'Ecos2'."
    end if
    if(sum(abs(epdir_im1(:)))>1.0d-8)then
      stop "Error: ae_shape1 should be 'Acos2' when epdir_im1 is used."
    end if
    tt = t - 0.5d0*tw1 - t1_start
    if (abs(tt)<0.5d0*tw1) then
      Ac_ext(:) = -epdir_re1(:)*f0_1/(8d0*pi**2*omega1 - 2d0*tw1**2*omega1**3) &
        *( &
        (-4d0*pi**2+tw1**2*omega1**2 + tw1**2*omega1**2*cos(2d0*pi*tt/tw1))*cos(omega1*tt) &
        +2d0*pi*(2d0*pi*cos(tw1*omega1/2d0) &
        +tw1*omega1*sin(2d0*pi*tt/tw1)*sin(omega1*tt)))
    end if
    T1_T2_tmp = T1_T2 + t1_start

  case('Esin2sin')
  
    stop "Esin2sin is not implemented"
    
  case('Asin2cos')
  
    ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
    ! pump laser
    tt = t
    if (tt<tw1) then
      Ac_ext(:) = -epdir_re1(:)*f0_1/omega1*(sin(pi*tt/tw1))**2*cos(omega1*tt+phi_CEP1*2d0*pi)
    end if
    
  case('input')
     
     stop "ae_shape1='input' is not implemented"
    
  case('Asin2_cw')
  
    ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
    ! pump laser
    tt = t
    if (tt<tw1*0.5d0) then
      Ac_ext(:) = -Epdir_re1(:)*f0_1/omega1*(sin(pi*tt/tw1))**2*cos(omega1*tt+phi_CEP1*2d0*pi)
    else
      Ac_ext(:) = -Epdir_re1(:)*f0_1/omega1*cos(omega1*tt+phi_CEP1*2d0*pi)
    end if
    
  case('none')
  
    Ac_ext = 0d0
    
  case default
  
    stop "Invalid pulse_shape_1 parameter!"
    
  end select

! Probe
  select case(ae_shape2)
  case('impulse')
  
    tt = t - 0.5d0*tw1 - T1_T2_tmp
    if(tt > 0d0)then
      Ac_ext(1) = Ac_ext(1) + epdir_re2(1)*e_impulse
      Ac_ext(2) = Ac_ext(2) + epdir_re2(2)*e_impulse
      Ac_ext(3) = Ac_ext(3) + epdir_re2(3)*e_impulse
    end if
    
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

    tt = t - 0.5d0*tw1 - T1_T2_tmp
    if (abs(tt)<0.5d0*tw2) then
      Ac_ext(:)=Ac_ext(:) &
        -f0_2/omega2*(cos(pi*tt/tw2))**npower &
        *aimag( (epdir_re2(:) + zI*epdir_im2(:)) &
        *exp(zI*(omega2*tt+phi_CEP2*2d0*pi))  &
        )
    end if

  case('Ecos2')
  
    if(phi_CEP2 /= 0.75d0)then
      stop "Error: phi_cep2 should be 0.75 when ae_shape2 is 'Ecos2'."
    end if
    if(sum(abs(epdir_im2(:)))>1.0d-8)then
      stop "Error: ae_shape2 should be 'Acos2' when epdir_im2 is used."
    end if
    tt = t - 0.5d0*tw1 - T1_T2_tmp
    if (abs(tt)<0.5d0*tw2) then
      Ac_ext(:)=Ac_ext(:) &
        -epdir_re2(:)*f0_2/(8d0*pi**2*omega2 - 2d0*tw2**2*omega2**3) &
        *( &
        (-4d0*pi**2+tw2**2*omega2**2 + tw2**2*omega2**2*cos(2d0*pi*tt/tw2))*cos(omega2*tt) &
        +2d0*pi*(2d0*pi*cos(tw2*omega2/2d0) &
        +tw2*omega2*sin(2d0*pi*tt/tw2)*sin(omega2*tt)))
    end if

  case('Esin2sin')
      
    stop "Esin2sin is not implemented"
    
  case('Asin2cos')
  
      ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
    ! probe laser
    tt = t
    if ( (tt-T1_T2_tmp>0d0) .and. (tt-T1_T2_tmp<tw2) ) then
      Ac_ext(:) = Ac_ext(:) &
        &-Epdir_re2(:)*f0_2/omega2*(sin(pi*(tt-T1_T2_tmp)/tw2))**2*cos(omega2*(tt-T1_T2_tmp)+phi_CEP2*2d0*pi)
    endif

  case('input')
  case('Asin2_cw')
  case('none')
  case default
  
    stop "Invalid pulse_shape_2 parameter!"
    
  end select

  return
End Subroutine calc_Ac_ext

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

end module em_field
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
