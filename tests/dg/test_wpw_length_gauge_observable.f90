program test_wpw_length_gauge_observable
  use rt_dg_wpw_length_gauge
  implicit none
  complex(8) :: s(2,2), z(2,2), h0(2,2), hplus(2,2), hminus(2,2), c(2,1)
  real(8) :: pz, delta_pz, jz, residual
  logical :: ok

  s = reshape([cmplx(1d0,0d0,8),cmplx(0.2d0,0d0,8), &
               cmplx(0.2d0,0d0,8),cmplx(1.1d0,0d0,8)],[2,2])
  z = reshape([cmplx(0.3d0,0d0,8),cmplx(0d0,-0.1d0,8), &
               cmplx(0d0,0.1d0,8),cmplx(-0.2d0,0d0,8)],[2,2])
  h0 = reshape([cmplx(-1d0,0d0,8),cmplx(0.05d0,0d0,8), &
                cmplx(0.05d0,0d0,8),cmplx(-0.4d0,0d0,8)],[2,2])
  c(:,1) = [cmplx(1d0,0d0,8),cmplx(0d0,0d0,8)]

  call validate_wpw_position_operator(s,z,1d-12,residual,ok)
  if(.not.ok .or. residual > 1d-12) error stop 'valid position rejected'
  call build_wpw_length_gauge_hamiltonian(h0,z,0.4d0,hplus,ok)
  if(.not.ok) error stop 'positive field rejected'
  call build_wpw_length_gauge_hamiltonian(h0,z,-0.4d0,hminus,ok)
  if(.not.ok .or. maxval(abs((hplus-h0)+(hminus-h0))) > 1d-14) error stop 'field oddness'
  call evaluate_wpw_polarization(z,c,[2d0],pz,ok)
  if(.not.ok .or. abs(pz-0.6d0)>1d-14) error stop 'polarization expectation'
  call update_wpw_polarization_branch(-4.8d0,4.9d0,10d0,0.2d0,delta_pz,jz,ok)
  if(.not.ok .or. abs(delta_pz-5.2d0)>1d-14 .or. abs(jz-1.5d0)>1d-14) &
    error stop 'continuous branch/current'

  z(1,2)=cmplx(0.3d0,0d0,8)
  call validate_wpw_position_operator(s,z,1d-12,residual,ok)
  if(ok) error stop 'non-Hermitian position accepted'
end program
