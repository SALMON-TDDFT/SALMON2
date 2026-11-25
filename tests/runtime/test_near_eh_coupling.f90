! Analytic E1/E2 and orbital M1 checks, without GS data.
program test_near_eh_coupling
 use salmonpy_solid_harmonics,only:near_potentials,near_line_integral
 use stencil_sub,only:zstencil_microAc_symmetrized
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 integer,parameter::a=-5,b=5
 real(8),parameter::c=137.035999084d0,h=.2d0,mag=.04d0,pi=acos(-1d0)
 real(8)::ce(3,1),ch(3,1),Ac(3,a:b,a:b,a:b),potential,vector(3),r(3),geometry(4)
 real(8)::from(3),finish(3),reference,integral,lapt(4,3),nabt(4,3),v(-1:1,-1:1,-1:1),error
 complex(8)::psi(a:b,a:b,a:b),out(a:b,a:b,a:b),expected
 integer::idx(a:b),i,j,k
 ce(:,1)=[.01d0,.02d0,-.03d0]*sqrt(4d0*pi)
 ch(:,1)=[0d0,0d0,mag]*sqrt(4d0*pi)
 geometry=[.3d0,-.4d0,.2d0,2d0]
 from=[1d0,.2d0,-.1d0];finish=[-.2d0,.8d0,.1d0]
 integral=near_line_integral(from,finish,geometry,0,ch)
 from=from-geometry(1:3);finish=finish-geometry(1:3)
 reference=mag*(from(1)*finish(2)-from(2)*finish(1))/(2*c)
 if(abs(integral-reference)>1d-17)error stop 'magnetic projector phase sign/normalization'
 do k=a,b
 do j=a,b
 do i=a,b
   r=[h*i+.35d0,h*j+.23d0,h*k]
   call near_potentials(r,0,2d0,ce,ch,potential,vector)
   reference=dot_product(r,[.01d0,.02d0,-.03d0])
   if(abs(potential-reference)>1d-16)error stop 'electric dipole normalization'
   Ac(:,i,j,k)=vector;psi(i,j,k)=cmplx(r(1),r(2),kind=8)
 enddo
 enddo
 enddo
 do i=a,b
  idx(i)=i
 enddo
 lapt=0d0;nabt=0d0;lapt(1,:)=1d0/h**2;nabt(1,:)=.5d0/h;v=0d0;out=0d0
 call zstencil_microAc_symmetrized([a,a,a],[b,b,b],[-1,-1,-1],[1,1,1],idx,idx,idx, &
                                  psi,out,v,Ac,3d0/h**2,lapt,nabt,[0d0,0d0,0d0])
 error=0d0
 do k=-1,1
 do j=-1,1
 do i=-1,1
   r=[h*i+.35d0,h*j+.23d0,h*k]
   expected=(mag/(2*c)+mag**2*(r(1)**2+r(2)**2)/(8*c*c))*psi(i,j,k)
   error=max(error,abs(out(i,j,k)-expected))
   if(.not.ieee_is_finite(real(out(i,j,k))))error stop 'nonfinite magnetic kinetic real part'
   if(.not.ieee_is_finite(aimag(out(i,j,k))))error stop 'nonfinite magnetic kinetic imaginary part'
 enddo
 enddo
 enddo
 if(error>2d-14)error stop 'orbital M1 or A squared kinetic action'
 print *, 'PASS: electric dipole, magnetic line phase, orbital Lz coupling and A squared',error
end program
