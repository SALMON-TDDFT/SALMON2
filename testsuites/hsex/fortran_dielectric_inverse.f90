program test_dielectric_inverse
  use hsex_screen_math, only: screen_inverse,screen_finite
  use, intrinsic :: ieee_exceptions, only: ieee_set_halting_mode,ieee_set_flag,ieee_invalid,ieee_divide_by_zero,ieee_overflow
  implicit none
  complex(8) :: matrix(3,3),inverse(3,3),residual(3,3)
  integer :: i,status
  matrix=0d0
  matrix(1,1)=1.0181490145761602d0
  matrix(2,2)=1.0178518737677067d0
  matrix(3,3)=1.0178518851325662d0
  matrix(2,1)=cmplx(.011785982628343955d0,-.0003799027353899603d0,8)
  matrix(3,1)=cmplx(-.011785986379939342d0,-.000379902856316776d0,8)
  matrix(3,2)=cmplx(-.01710915638237355d0,.0005305164769728233d0,8)
  do i=1,3
    matrix(i,i+1:)=conjg(matrix(i+1:,i))
  enddo
  ! Only the production inversion is under test here. Reference arithmetic
  ! is outside the trap region: some compilers relax complex division/ABS.
  call ieee_set_flag(ieee_invalid,.false.)
  call ieee_set_flag(ieee_divide_by_zero,.false.)
  call ieee_set_flag(ieee_overflow,.false.)
  call ieee_set_halting_mode(ieee_invalid,.true.)
  call ieee_set_halting_mode(ieee_divide_by_zero,.true.)
  call ieee_set_halting_mode(ieee_overflow,.true.)
  call screen_inverse(matrix,inverse,status)
  call ieee_set_halting_mode(ieee_invalid,.false.)
  call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
  call ieee_set_halting_mode(ieee_overflow,.false.)
  if(status/=0.or..not.screen_finite(inverse))error stop 'finite dielectric inversion rejected'
  residual=matmul(matrix,inverse)
  do i=1,3
    residual(i,i)=residual(i,i)-1d0
  enddo
  if(.not.screen_finite(residual))error stop 'nonfinite dielectric inverse residual'
  if(maxval(abs(residual))>5d-14)error stop 'dielectric inverse residual'
  print *, 'DIELECTRIC_INVERSE_IEEE_PASS',maxval(abs(residual))
end program
