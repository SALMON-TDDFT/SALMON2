program test_ld_api
  use lorentz_drude, only: s_ld_state, ld_update, ld_store_old, ld_coefficients
  implicit none
  type(s_ld_state) :: ld
  print *, 'LD module API available'
end program
