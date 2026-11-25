program test_electric_path_response
 use iso_c_binding
 use salmonpy_grid_response,only:sp_grid_electric_response
 use salmonpy_grid_field,only:grid_potentials
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 integer,parameter::n=19
 integer::counts(3),i,j,k,d,a,status
 real(8)::origins(3,3,2),spacing(3),position(3),center(3),points(3,n),old(n),new(n),volumes(n),r(3),x(3)
 real(8)::data(3,7,8,9,2),P(3,7,8,9),Jgrid(3,7,8,9),V(n),Ac(3),dt,dv,work,reference,total(3),expected(3)
 counts=[7,8,9];spacing=[.5d0,.4d0,.3d0];origins=-1.2d0
 origins(1,1,1)=-1.17d0;origins(2,2,1)=-1.24d0;origins(3,3,1)=-1.18d0
 position=[.17d0,-.13d0,.21d0];center=[.2d0,-.3d0,.1d0];dt=.02d0;dv=product(spacing)
 data=0d0
 do k=1,9
 do j=1,8
 do i=1,7
 do d=1,3
  x=origins(:,d,1)+[i-1,j-1,k-1]*spacing
  data(d,i,j,k,1)=.01d0*(sin(x(1)*d)+x(2)*x(3)+d)
 enddo
 enddo
 enddo
 enddo
 do a=1,n
  points(:,a)=center+.7d0*[sin(real(a,8)),cos(real(2*a,8)),sin(real(3*a,8))]
  old(a)=-.1d0-.01d0*a;new(a)=old(a)+.001d0*cos(real(a,8));volumes(a)=.001d0+.0001d0*a
  r=points(:,a)-center
  call grid_potentials(counts,origins,spacing,data,position,r,V(a),Ac,status)
  if(status/=0)error stop 'forward potential failed'
 enddo
 status=sp_grid_electric_response(0,n,counts,origins(:,:,1),spacing,position,center,points,volumes,old,new,dt,P)
 if(status/=0)error stop 'polarization failed'
 status=sp_grid_electric_response(1,n,counts,origins(:,:,1),spacing,position,center,points,volumes,old,new,dt,Jgrid)
 if(status/=0)error stop 'current failed'
 work=sum(data(:,:,:,:,1)*Jgrid)*dv*dt;reference=sum((new-old)*volumes*V)
 total=0d0;expected=0d0
 do k=1,9
 do j=1,8
 do i=1,7
  total=total+P(:,i,j,k)*dv
 enddo
 enddo
 enddo
 do a=1,n
  expected=expected+new(a)*volumes(a)*(points(:,a)-center)
 enddo
 print *,'spatial work absolute/relative error=',abs(work-reference),abs(work-reference)/max(abs(reference),tiny(1d0))
 print *,'integrated dipole absolute error=',maxval(abs(total-expected))
 if(abs(work-reference)>1d-17.or.maxval(abs(total-expected))>1d-17)error stop 'transpose identity failed'
 P=7d0;points(1,n)=100d0
 status=sp_grid_electric_response(1,n,counts,origins(:,:,1),spacing,position,center,points,volumes,old,new,dt,P)
 if(status/=1.or.any(P/=7d0))error stop 'domain rejection changed output'
 points(:,n)=center
 new(n)=ieee_value(1d0,ieee_quiet_nan)
 status=sp_grid_electric_response(1,n,counts,origins(:,:,1),spacing,position,center,points,volumes,old,new,dt,P)
 if(status/=1.or.any(P/=7d0))error stop 'NaN rejection changed output'
 new(n)=ieee_value(1d0,ieee_positive_inf)
 status=sp_grid_electric_response(1,n,counts,origins(:,:,1),spacing,position,center,points,volumes,old,new,dt,P)
 if(status/=1.or.any(P/=7d0))error stop 'Infinity rejection changed output'
 new(n)=huge(1d0);volumes(n)=huge(1d0)
 status=sp_grid_electric_response(1,n,counts,origins(:,:,1),spacing,position,center,points,volumes,old,new,dt,P)
 if(status/=4.or.any(P/=7d0))error stop 'overflow rejection changed output'
 print *,'PASS: staggered E, off-grid molecule, cell-crossing paths, work/dipole identities and atomic failure'
end program
