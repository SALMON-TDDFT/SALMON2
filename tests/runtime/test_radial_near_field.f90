! Independent analytic fields at positions not used for fitting; no GS or SALMON launch.
program test_radial_near_field
 use iso_c_binding
 use salmonpy_solid_harmonics,only:sp_near_fit,sp_near_evaluate
 implicit none
 integer,parameter::n=400,nv=53
 real(c_double)::xyz(3,n),points(3,nv),fields(3,n,2),geo(7),ce(3,165),ch(3,165),stats(2)
 real(c_double)::E(3,nv),H(3,nv),potential(nv),Ac(3,nv),r(3),r2,radial,expected(3),cross(3)
 real(c_double)::ev(3),hv(3),field_error,potential_error,magnetic_error
 integer::seedsize,i,degree,status,p
 integer,allocatable::seed(:)
 call random_seed(size=seedsize);allocate(seed(seedsize));seed=952;call random_seed(put=seed)
 call random_number(xyz);call random_number(points)
 geo=[.17d0,-.23d0,.31d0,2d0,0d0,0d0,0d0]
 do i=1,n
  xyz(:,i)=2d0*xyz(:,i)-1d0+geo(1:3)
 enddo
 do i=1,nv
  points(:,i)=1.6d0*points(:,i)-.8d0+geo(1:3)
 enddo
 ev=[.01d0,.02d0,-.03d0];hv=[.003d0,-.004d0,.002d0]
 do degree=2,8,2
  p=(degree+1)*(degree+2)*(degree+3)/6
  do i=1,n
   r=xyz(:,i)-geo(1:3);radial=dot_product(r,r)**(degree/2)
   fields(:,i,1)=radial*ev;fields(:,i,2)=radial*hv
  enddo
  ce=0d0;ch=0d0;stats=0d0
  status=sp_near_fit(n,degree,xyz,fields,geo,0,0,ce,ch,stats)
  if(status/=0)error stop 'radial fit failed'
  status=sp_near_evaluate(nv,degree,points,geo,ce,ch,E,H,potential,Ac)
  if(status/=0)error stop 'radial evaluation failed'
  field_error=0d0;potential_error=0d0;magnetic_error=0d0
  do i=1,nv
   r=points(:,i)-geo(1:3);r2=dot_product(r,r);radial=r2**(degree/2)
   field_error=max(field_error,maxval(abs(E(:,i)-radial*ev)),maxval(abs(H(:,i)-radial*hv)))
   potential_error=max(potential_error,abs(potential(i)-dot_product(r,radial*ev)/(degree+1)))
   expected=radial*hv
   cross=[r(2)*expected(3)-r(3)*expected(2),r(3)*expected(1)-r(1)*expected(3), &
          r(1)*expected(2)-r(2)*expected(1)]
   magnetic_error=max(magnetic_error,maxval(abs(Ac(:,i)+cross/((degree+2)*137.035999084d0))))
  enddo
  if(field_error>2d-14.or.potential_error>2d-14.or.magnetic_error>2d-16) &
   error stop 'radial off-sample field/potential mismatch'
  print '(a,i2,a,i3,a,3es24.15)','degree=',degree,' coefficients=',p,' absolute errors E/H,V,Ac=', &
       field_error,potential_error,magnetic_error
 enddo
 print *, 'PASS: radial E/H and analytic radial-gauge potentials at independent off-grid positions'
end program
