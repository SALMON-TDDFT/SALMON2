program probe
 use hse_symmetry
 use, intrinsic :: ieee_arithmetic
 implicit none
 type(hse_symmetry_map):: m
 real(8):: a(3,4,8), b(3,4,8), k(3,1),w(1),pi,atoms(3,2)
 complex(8),allocatable:: s(:,:,:),t(:,:,:),es(:,:,:),et(:,:,:)
 complex(8):: transformed(64,1)
 complex(8):: oracle
 integer:: op,d,g,x,y,z,ierr,j,q,qj,xj,yj,zj
 pi=acos(-1d0); a=0; b=0
 do op=1,8
 do d=1,3
 a(d,d,op)=1-2*ibits(op-1,d-1,1)
 enddo
 enddo
 ! Put identity last to ensure representative targets never acquire a gauge transform.
 a(:,:,1)=-a(:,:,1);a(:,:,8)=-a(:,:,8)
 b=a;k(:,1)=pi/8;w=1
 call symmetry_init(m,[4,4,4],[1d0,1d0,1d0],k,w,a,b,2,ierr)
 if(ierr/=0)stop 1
 if(m%nfull/=8.or.m%max_little/=1.or.m%first(1)/=1.or.m%first(2)/=9)stop 2
 allocate(s(64,1,1),t(64,1,1))
 do g=1,64
 x=mod(g-1,4); y=mod((g-1)/4,4);z=(g-1)/16
 s(g,1,1)=exp(cmplx(0d0,pi/2*(x+2*y+z),8))
 enddo
 t=s
 call symmetry_expand(m,1,s,t,es,et,ierr)
 if(ierr/=0)stop 3
 if(maxval(abs(et(:,:,1)-t(:,:,1)))>1d-12)stop 4
 do j=1,8
 op=m%operations(1,j)
 do g=1,64
 x=mod(g-1,4);y=mod((g-1)/4,4);z=(g-1)/16
 if(abs(et(g,1,j)-exp(cmplx(0d0,pi/2*(a(1,1,op)*x+2*a(2,2,op)*y+a(3,3,op)*z),8)))>1d-12)stop 5
 enddo
 enddo
 ! Gamma: all 8 operations belong to the little group; source is projector average.
 k=0
 call symmetry_init(m,[4,4,4],[1d0,1d0,1d0],k,w,a,b,1,ierr)
 if(ierr/=0.or.m%max_little/=8.or.m%multiplicity(1)/=8)stop 6
 call symmetry_expand(m,1,s,t,es,et,ierr)
 if(ierr/=0)stop 7
 if(abs(sum(abs(es)**2)-sum(abs(s)**2))>1d-11)stop 8
 ! Independent real-space projector orbit average for arbitrary grid pairs.
 do g=1,64
 x=mod(g-1,4);y=mod((g-1)/4,4);z=(g-1)/16
 do j=1,64
 xj=mod(j-1,4);yj=mod((j-1)/4,4);zj=(j-1)/16
 oracle=0
 do op=1,8
 q=1+modulo(nint(a(1,1,op))*x,4)+4*modulo(nint(a(2,2,op))*y,4)+16*modulo(nint(a(3,3,op))*z,4)
 qj=1+modulo(nint(a(1,1,op))*xj,4)+4*modulo(nint(a(2,2,op))*yj,4)+16*modulo(nint(a(3,3,op))*zj,4)
 oracle=oracle+s(q,1,1)*conjg(s(qj,1,1))/8
 enddo
 if(abs(sum(es(g,:,1)*conjg(es(j,:,1)))-oracle)>1d-12)stop 13
 enddo
 enddo
 call symmetry_expand(m,1,s,t(:63,:,:),es,et,ierr)
 if(ierr==0)stop 14
 t(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
 call symmetry_expand(m,1,s,t,es,et,ierr)
 if(ierr==0)stop 15
 t=s
 ! Half translation with reciprocal wrapping, k outside canonical cell retained first.
 a=0;b=0
 do d=1,3
 a(d,d,1)=1;a(d,d,2)=-1
 enddo
 a(1,4,2)=0.5d0;b=a;k(:,1)=pi/2
 call symmetry_init(m,[4,4,4],[1d0,1d0,1d0],k,w,a(:,:,1:2),b(:,:,1:2),1,ierr)
 if(ierr/=0.or.m%multiplicity(1)/=2)stop 9
 call symmetry_expand(m,1,s,t,es,et,ierr)
 if(ierr/=0)stop 10
 call symmetry_transform(m,1,2,s(:,:,1),transformed,ierr)
 if(ierr/=0)stop 21
 if(maxval(abs(transformed/sqrt(2d0)-es(:,2:2,1)))>1d-12)stop 22
 call symmetry_transform(m,1,3,s(:,:,1),transformed,ierr)
 if(ierr==0)stop 23
 do g=1,64
 x=mod(g-1,4);y=mod((g-1)/4,4);z=(g-1)/16
 q=1+modulo(2-x,4)+4*modulo(-y,4)+16*modulo(-z,4)
 if(abs(es(g,2,1)-s(q,1,1)*exp(cmplx(0d0,-pi*(x+y+z)+pi,8))/sqrt(2d0))>1d-12)stop 11
 enddo
 atoms=0;atoms(1,2)=2d0
 call symmetry_validate_atoms(m,atoms,[1,1],ierr)
 if(ierr/=0)stop 16
 call symmetry_validate_atoms(m,atoms,[1,2],ierr)
 if(ierr==0)stop 17
 atoms(1,2)=1d0
 call symmetry_validate_atoms(m,atoms,[1,1],ierr)
 if(ierr==0)stop 18
 w=0.5d0
 call symmetry_init(m,[4,4,4],[1d0,1d0,1d0],k,w,a(:,:,1:2),b(:,:,1:2),1,ierr)
 if(ierr==0)stop 12
 ! Non-group operation subset (I and 90 degree rotation, missing its square).
 a=0;b=0
 do d=1,3
 a(d,d,1)=1
 enddo
 a(1,2,2)=-1;a(2,1,2)=1;a(3,3,2)=1;b=a;w=1;k=0
 call symmetry_init(m,[4,4,4],[1d0,1d0,1d0],k,w,a(:,:,1:2),b(:,:,1:2),1,ierr)
 if(ierr==0)stop 19
 ! Translation set must also close modulo lattice vectors.
 a(:,:,2)=a(:,:,1);a(1,4,2)=0.25d0;b=a
 call symmetry_init(m,[4,4,4],[1d0,1d0,1d0],k,w,a(:,:,1:2),b(:,:,1:2),1,ierr)
 if(ierr==0)stop 20
 print *, 'symmetry probe passed'
end program
