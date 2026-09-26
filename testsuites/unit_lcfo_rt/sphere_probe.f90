program sphere_probe
 use lcfo_rt_wannier
 use lcfo_rt_basis
 implicit none
 complex(8) :: c(1,1)
 complex(8),allocatable :: source(:,:,:)
 real(8) :: d(3),pos(3),center(3),r,expected,err,norm
 integer :: x,y,z,g
 character(64) :: arg
 logical :: weak
 call get_command_argument(1,arg);read(arg,*)r
 call get_command_argument(2,arg);weak=trim(arg)=='weak'
 lcfo_grid=[8,8,8];lcfo_core=lcfo_grid;lcfo_h=[1d0,1.25d0,1.5d0]
 lcfo_dv=product(lcfo_h);center=[7d0,7.5d0,1.5d0]
 allocate(lcfo_basis(512,1),lcfo_counts(1),lcfo_offsets(2),lcfo_origins(3,1))
 lcfo_counts=1;lcfo_offsets=[0,1];lcfo_origins=0;g=0
 do z=0,7;do y=0,7;do x=0,7
  g=g+1;pos=[x,y,z]*lcfo_h
  d=modulo(pos-center+.5d0*lcfo_grid*lcfo_h,lcfo_grid*lcfo_h)-.5d0*lcfo_grid*lcfo_h
  if(weak)d(2)=0d0
  lcfo_basis(g,1)=exp(-sum(d*d)/8d0)
  if(weak)lcfo_basis(g,1)=lcfo_basis(g,1)*sqrt(1d0+.1d0*cos(2*acos(-1d0)*(pos(2)-center(2))/10d0))
 enddo;enddo;enddo
 norm=sqrt(sum(abs(lcfo_basis)**2)*lcfo_dv);lcfo_basis=lcfo_basis/norm;c=1d0
 call lcfo_mlwf_configure();call lcfo_mlwf_source(c,lcfo_basis,[1],source)
 err=0;g=0
 do z=0,7;do y=0,7;do x=0,7
  g=g+1;pos=[x,y,z]*lcfo_h
  d=modulo(pos-center+.5d0*lcfo_grid*lcfo_h,lcfo_grid*lcfo_h)-.5d0*lcfo_grid*lcfo_h
  expected=abs(lcfo_basis(g,1))
  if(r>0d0.and.sum(d*d)>r*r.and..not.weak)expected=0d0
  err=max(err,abs(abs(source(g,1,1))-expected))
 enddo;enddo;enddo
 if(err>1d-12)error stop '3D periodic spherical mask mismatch'
 write(*,*)'3D source mask passed: radius, protected, max error',r,weak,err
end program
