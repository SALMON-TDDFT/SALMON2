program test_magnetic_path_response
 use salmonpy_grid_field,only:grid_potentials,grid_magnetic_add,grid_line_phase,grid_phase_magnetic_add
 use salmonpy_grid_response,only:sp_grid_magnetization_current
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 integer::counts(3),i,j,k,d,status
 real(8)::origins(3,3,2),spacing(3),position(3),r(3),data(3,9,10,11,2),M(3,9,10,11),Jm(3,9,10,11)
 real(8)::V,A(3),q(3),work,reference,start(3),finish(3),phase,dv,divergence
 counts=[9,10,11];spacing=[.5d0,.4d0,.3d0];origins=-1.2d0
 do d=1,3
  do k=1,3
   if(k==d)then
    origins(k,d,1)=origins(k,d,1)+.5d0*spacing(k)
   else
    origins(k,d,2)=origins(k,d,2)+.5d0*spacing(k)
   endif
  enddo
 enddo
 data=0d0
 do k=1,11
 do j=1,10
 do i=1,9
 do d=1,3
  data(d,i,j,k,2)=sin(real(i+2*j+3*k+d,8))*.01d0
 enddo
 enddo
 enddo
 enddo
 position=[.17d0,-.13d0,.21d0];r=[.7d0,.6d0,.5d0];q=[.2d0,-.3d0,.4d0];dv=product(spacing)
 call grid_potentials(counts,origins,spacing,data,position,r,V,A,status)
 if(status/=0)error stop 'forward'
 M=0d0
 call grid_magnetic_add(counts,origins(:,:,2),spacing,position,r,q,M,status)
 if(status/=0)error stop 'adjoint'
 work=sum(data(:,:,:,:,2)*M)*dv;reference=-dot_product(q,A)
 if(abs(work-reference)>1d-17)error stop 'magnetic transpose work'
 start=[.1d0,.2d0,.3d0];finish=[.7d0,.8d0,.9d0]
 call grid_line_phase(counts,origins,spacing,data,position,[0d0,0d0,0d0],start,finish,phase,status)
 if(status/=0)error stop 'phase'
 M=0d0
 call grid_phase_magnetic_add(counts,origins(:,:,2),spacing,position,[0d0,0d0,0d0],start,finish,2d0,M,status)
 if(status/=0)error stop 'phase adjoint'
 work=sum(data(:,:,:,:,2)*M)*dv
 if(abs(work+2d0*phase)>1d-17)error stop 'nonlocal phase transpose'
 status=sp_grid_magnetization_current(counts,spacing,M,Jm)
 if(status/=0)error stop 'curl'
 divergence=0d0
 do k=2,11
 do j=2,10
 do i=2,9
  divergence=max(divergence,abs((Jm(1,i,j,k)-Jm(1,i-1,j,k))/spacing(1)+ &
   (Jm(2,i,j,k)-Jm(2,i,j-1,k))/spacing(2)+(Jm(3,i,j,k)-Jm(3,i,j,k-1))/spacing(3)))
 enddo
 enddo
 enddo
 if(divergence>1d-14)error stop 'curl creates charge'
 Jm=7d0;M=0d0;M(1,3,3,3)=ieee_value(1d0,ieee_quiet_nan)
 status=sp_grid_magnetization_current(counts,spacing,M,Jm)
 if(status/=1.or.any(Jm/=7d0))error stop 'NaN changed current output'
 M(1,3,3,3)=ieee_value(1d0,ieee_positive_inf)
 status=sp_grid_magnetization_current(counts,spacing,M,Jm)
 if(status/=1.or.any(Jm/=7d0))error stop 'Infinity changed current output'
 M(1,3,3,3)=huge(1d0)*.5d0
 status=sp_grid_magnetization_current(counts,spacing,M,Jm)
 if(status/=4.or.any(Jm/=7d0))error stop 'overflow changed current output'
 print *,'PASS: magnetic and projector-path work, Yee divergence=',divergence
end program
