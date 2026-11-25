! Direct-grid interpolation and Iwasa path integrals: no GS or TDDFT startup.
program test_direct_grid_field
 use iso_c_binding
 use salmonpy_grid_field,only:sp_grid_fields,sp_grid_potentials,grid_line_phase
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 integer,parameter::nq=97
 real(8),parameter::gauss(3)=[-.774596669241483377d0,0d0,.774596669241483377d0]
 real(8),parameter::gw(3)=[5d0/9d0,8d0/9d0,5d0/9d0]
 integer::counts(3),i,j,k,d,which,n,level,status,seedsize
 integer,allocatable::seed(:)
 real(8),allocatable::data(:,:,:,:,:)
 real(8)::origins(3,3,2),spacing(3),field_center(3),center(3),points(3,nq),relative(3,nq)
 real(8)::E(3,nq),H(3,nq),V(nq),Ac(3,nq),expectedE(3,nq),expectedH(3,nq),expectedV(nq)
 real(8)::x(3),r(3),G(3,3),GH(3,3),ev(3),hv(3),eint(3),hint(3),cross(3),error
 real(8)::phase,reference_phase,t,vec(3),Eerror(3),Verror(3),previous,scalar
 call random_seed(size=seedsize);allocate(seed(seedsize));seed=3952;call random_seed(put=seed)
 call random_number(relative);relative=4d0*(relative-.5d0)
 field_center=[.17d0,-.23d0,.31d0];center=[.2d0,-.3d0,.1d0]
 do i=1,nq
  points(:,i)=center+relative(:,i)
 enddo
 counts=9;spacing=1d0
 ev=[.01d0,.02d0,-.03d0];hv=[.003d0,-.004d0,.002d0]
 G=reshape([.01d0,.002d0,0d0,.002d0,-.02d0,.001d0,0d0,.001d0,.01d0],[3,3])
 GH=reshape([.002d0,.001d0,0d0,0d0,-.001d0,0d0,0d0,0d0,-.001d0],[3,3])
 do which=1,2
 do d=1,3
  origins(:,d,which)=-4d0
  do k=1,3
   if((which==1.and.k==d).or.(which==2.and.k/=d))origins(k,d,which)=-3.5d0
  enddo
 enddo
 enddo
 allocate(data(3,9,9,9,2))
 do which=1,2
 do k=1,9
 do j=1,9
 do i=1,9
 do d=1,3
  x=origins(:,d,which)+[i-1,j-1,k-1]*spacing
  if(which==1)then
   vec=ev+matmul(G,x)
  else
   vec=hv+matmul(GH,x)
  endif
  data(d,i,j,k,which)=vec(d)
 enddo
 enddo
 enddo
 enddo
 enddo
 status=sp_grid_fields(nq,counts,origins,spacing,field_center,center,points,data,E,H)
 if(status/=0)error stop 'staggered interpolation failed'
 status=sp_grid_potentials(nq,counts,origins,spacing,field_center,center,points,data,V,Ac)
 if(status/=0)error stop 'grid path integration failed'
 error=0d0
 do i=1,nq
  r=relative(:,i);x=field_center+r
  expectedE(:,i)=ev+matmul(G,x);expectedH(:,i)=hv+matmul(GH,x)
  eint=ev+matmul(G,field_center)+matmul(G,r)/2d0
  hint=(hv+matmul(GH,field_center))/2d0+matmul(GH,r)/3d0
  cross=[r(2)*hint(3)-r(3)*hint(2),r(3)*hint(1)-r(1)*hint(3),r(1)*hint(2)-r(2)*hint(1)]
  error=max(error,maxval(abs(E(:,i)-expectedE(:,i))),maxval(abs(H(:,i)-expectedH(:,i))), &
    abs(V(i)-dot_product(r,eint)),maxval(abs(Ac(:,i)+cross/137.035999084d0)))
 enddo
 if(error>3d-15)error stop 'staggered affine field or cell-crossing integral is inaccurate'
 print *, 'staggered off-grid affine E/H/potential max absolute error=',error
 ! Independent three-point reference is exact for the quadratic A from affine H.
 call grid_line_phase(counts,origins,spacing,data,field_center,center, &
   center+[-.7d0,.3d0,-.4d0],center+[.9d0,-.2d0,1.1d0],phase,status)
 if(status/=0)error stop 'grid magnetic line failed'
 reference_phase=0d0
 do i=1,3
  t=.5d0*(1d0+gauss(i))
  r=[-.7d0,.3d0,-.4d0]+t*[1.6d0,-.5d0,1.5d0]
  hint=(hv+matmul(GH,field_center))/2d0+matmul(GH,r)/3d0
  cross=[r(2)*hint(3)-r(3)*hint(2),r(3)*hint(1)-r(1)*hint(3),r(1)*hint(2)-r(2)*hint(1)]
  reference_phase=reference_phase-dot_product(cross,[1.6d0,-.5d0,1.5d0])/137.035999084d0* &
    .5d0*gw(i)
 enddo
 if(abs(phase-reference_phase)>1d-17)error stop 'magnetic line disagrees with analytic polynomial'
 ! All source values are checked, even outside the interpolation stencil; outputs stay intact.
 do i=1,2
  scalar=data(1,1,1,1,1)
  if(i==1)data(1,1,1,1,1)=ieee_value(0d0,ieee_quiet_nan)
  if(i==2)data(1,1,1,1,1)=ieee_value(0d0,ieee_positive_inf)
  E=7d0;H=8d0
  status=sp_grid_fields(nq,counts,origins,spacing,field_center,center,points,data,E,H)
  if(status==0.or.any(E/=7d0).or.any(H/=8d0))error stop 'nonfinite source changed outputs'
  data(1,1,1,1,1)=scalar
 enddo
 points(3,nq)=100d0;E=7d0;H=8d0
 status=sp_grid_fields(nq,counts,origins,spacing,field_center,center,points,data,E,H)
 if(status==0.or.any(E/=7d0).or.any(H/=8d0))error stop 'out-of-grid query changed outputs'
 deallocate(data)
 ! Smooth source-free electrostatic field: interpolation error, not a fitted expansion.
 call random_number(relative);relative=1.4d0*(relative-.5d0)
 field_center=[.07d0,-.13d0,.11d0]
 do i=1,nq
  points(:,i)=center+relative(:,i);x=field_center+relative(:,i)
  expectedE(:,i)=exp(-2d0*x(3))*[sin(2d0*x(1)),0d0,cos(2d0*x(1))]
  expectedV(i)=.5d0*(exp(-2d0*field_center(3))*cos(2d0*field_center(1))-exp(-2d0*x(3))*cos(2d0*x(1)))
 enddo
 do level=1,3
  n=2**(level+2)+1;counts=n;spacing=4d0/(n-1);origins=-2d0
  allocate(data(3,n,n,n,2));data=0d0
  do k=1,n
  do j=1,n
  do i=1,n
   x=-2d0+[i-1,j-1,k-1]*spacing
   data(:,i,j,k,1)=exp(-2d0*x(3))*[sin(2d0*x(1)),0d0,cos(2d0*x(1))]
  enddo
  enddo
  enddo
  status=sp_grid_fields(nq,counts,origins,spacing,field_center,center,points,data,E,H)
  if(status/=0)error stop 'refined grid interpolation failed'
  status=sp_grid_potentials(nq,counts,origins,spacing,field_center,center,points,data,V,Ac)
  if(status/=0)error stop 'refined grid path failed'
  Eerror(level)=sqrt(sum((E-expectedE)**2)/sum(expectedE**2))
  Verror(level)=sqrt(sum((V-expectedV)**2)/sum(expectedV**2))
  print *, 'spacing, E/V relative RMS error=',spacing(1),Eerror(level),Verror(level)
  deallocate(data)
 enddo
 if(Eerror(3)>Eerror(1)/10d0.or.Verror(3)>Verror(1)/10d0)error stop 'direct-grid refinement did not converge'
 print *, 'PASS: staggered interpolation, cell faces, orbital magnetic phase, finite guards, grid refinement'
end program
