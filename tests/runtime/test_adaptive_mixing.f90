program test_adaptive_mixing
 use structures,only:s_mixing
 use mixing_sub,only:check_mixing_half
 use parallelization,only:setup_parallel,end_parallel
 implicit none
 type(s_mixing)::a,b,saved
 integer::i
 call setup_parallel
 if(a%input%mixrate/=0.3d0.or.a%input%alpha_mb/=0.3d0.or.a%input%beta_p/=0.3d0)error stop 'mixing defaults differ'
 a%input%method='broyden';a%input%auto_mixing='y';a%alpha_mb=0.4d0
 a%adaptive%min_rate=0.05d0;a%adaptive%growth=1.2d0
 a%adaptive%patience=3;a%adaptive%cooldown=2
 call check_mixing_half(1,1d0,a)
 call check_mixing_half(2,2d0,a)
 call check_mixing_half(3,1d0,a)
 if(a%alpha_mb/=0.4d0)error stop 'isolated rise caused reduction'
 call check_mixing_half(4,2d0,a)
 call check_mixing_half(5,2d0,a)
 saved=a;b=a
 call check_mixing_half(6,2d0,a)
 if(a%alpha_mb/=0.2d0)error stop 'persistent growth not reduced'
 if(b%alpha_mb/=0.4d0)error stop 'sibling modified'
 b=saved;call check_mixing_half(6,2d0,b)
 if(a%alpha_mb/=b%alpha_mb.or.a%adaptive%changes/=b%adaptive%changes)error stop 'resume mismatch'
 call check_mixing_half(7,8d0,a);call check_mixing_half(8,8d0,a)
 if(a%alpha_mb/=0.2d0)error stop 'cooldown ignored'
 do i=9,80
  call check_mixing_half(i,2d0**min(i-7,30),a)
 enddo
 if(a%alpha_mb/=0.05d0)error stop 'floor violated or not reached'
 b=saved;b%input%auto_mixing='n'
 call check_mixing_half(6,100d0,b)
 if(b%alpha_mb/=saved%alpha_mb)error stop 'fixed mode changed'
 print *, 'ADAPTIVE_MIXING_PASS'
 call end_parallel
end program
