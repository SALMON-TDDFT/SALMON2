program occupation_probe
 use ptcn_core,only:pt_uniform_occupation
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 real(8)::a(16,2,1)
 a=2d0
 if(.not.pt_uniform_occupation(a,2d0))error stop 'occupied-only rejected'
 a(16,2,1)=0d0
 if(pt_uniform_occupation(a,2d0))error stop 'empty accepted'
 a(16,2,1)=1d0
 if(pt_uniform_occupation(a,2d0))error stop 'unequal accepted'
 a(16,2,1)=ieee_value(0d0,ieee_quiet_nan)
 if(pt_uniform_occupation(a,2d0))error stop 'NaN accepted'
 a(16,2,1)=ieee_value(0d0,ieee_positive_inf)
 if(pt_uniform_occupation(a,2d0))error stop 'Inf accepted'
 print *,'OCCUPIED_ONLY_PASS'
end program
