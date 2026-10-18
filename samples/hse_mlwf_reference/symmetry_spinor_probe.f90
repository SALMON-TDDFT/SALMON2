program probe
 use exx_symmetry
 use exx_k_exchange
 use exx_spinor
 implicit none
 type(s_exx_symmetry_map)::m
 real(8)::a(3,4,4),k(3,1),w(1),h(3)
 complex(8)::u(2,2),v(2,2),psi(96,1),out(96,1)
 type(exx_k_kernel)::op
 complex(8),allocatable::es(:,:,:),et(:,:,:),act(:,:,:)
 complex(8)::ref(96,1,1),dm
 integer::i,j,e,g,q,b,t,delta(3)
 a=0; k=0;w=1;h=[1d0,2d0,3d0]
 do j=1,4
 do i=1,3
 a(i,i,j)=1
 enddo
 enddo
 a(1,1,2)=-1;a(3,3,3)=-1;a(1,1,4)=-1;a(3,3,4)=-1
 call symmetry_init(m,[4,6,2],h,k,w,a,a,[1,1,1],e)
 if(e/=0)stop 1
 call symmetry_spin_rotation(a(:,1:3,2),u,e)
 if(e/=0)stop 2
 v=matmul(u,conjg(transpose(u)))
 if(abs(v(1,1)-1)>1d-14.or.abs(v(2,2)-1)>1d-14.or.abs(v(1,2))>1d-14)stop 3
 ! Mx is axial pi rotation around x, mixing both spin components.
 psi(:48,1)=1;psi(49:,1)=2
 call symmetry_transform(m,1,2,psi,out,e)
 if(e/=0)stop 4
 if(maxval(abs(out(:48,1)-(0d0,-2d0)))>1d-14)stop 5
 if(maxval(abs(out(49:,1)-(0d0,-1d0)))>1d-14)stop 6
 call symmetry_spin_rotation(a(:,1:3,3),u,e)
 if(e/=0.or.abs(u(1,1)-(0d0,-1d0))>1d-14.or.abs(u(2,2)-(0d0,1d0))>1d-14)stop 15
 ! U sigma_x Udag = -sigma_x for horizontal mirror (axial spin).
 v=0;v(1,2)=1;v(2,1)=1
 v=matmul(u,matmul(v,conjg(transpose(u))))
 if(abs(v(1,2)+1)>1d-14.or.abs(v(2,1)+1)>1d-14)stop 16
 ! Little-group projector averaging must preserve cross-spin density blocks.
 do i=1,96
 psi(i,1)=cmplx(sin(real(i,8)),cos(real(2*i,8)),8)/7d0
 enddo
 call symmetry_expand(m,1,reshape(psi,[96,1,1]),reshape(psi,[96,1,1]),es,et,e)
 if(e/=0)stop 7
 allocate(act(96,1,1))
 call exx_k_kernel_init(op,[4,6,2],[1,1,1],h,k,.11d0,3,e,fft_layout='strided')
 if(e/=0)stop 8
 call exx_spinor_apply_distributed(op,es,et,act,[1],[1],0,0,transpose_tiles,e)
 if(e/=0)stop 9
 ref=0
 do i=1,48;do j=1,48
 delta=modulo(op%point(:,i)-op%point(:,j),op%n)
 do b=1,2;do t=1,2
 dm=sum(es(i+(b-1)*48,:,1)*conjg(es(j+(t-1)*48,:,1)))
 ref(i+(b-1)*48,1,1)=ref(i+(b-1)*48,1,1)-dm*op%kernel(delta(1),delta(2),delta(3))*et(j+(t-1)*48,1,1)
 enddo;enddo
 enddo;enddo
 if(maxval(abs(ref-act))>1d-11)stop 10
 ! An incompatible fractional translation is rejected before any FFT work.
 a(:,4,2)=[.25d0,0d0,0d0]
 call symmetry_init(m,[4,6,2],h,k,w,a,a,[1,1,1],e)
 if(e==0)stop 11
 block
 type(s_exx_symmetry_map)::stars
 real(8)::ops(3,4,4),ks(3,2),ws(2),pi
 complex(8)::input(96,1),output(96,1),phase
 integer::r,c(3)
 ops=0;pi=acos(-1d0);ks=0;ks(1,2)=-2*pi/(3*4*h(1));ws=[1d0/3,2d0/3]
 do r=1,4
 do i=1,3
 ops(i,i,r)=1
 enddo
 enddo
 ops(1,1,2)=-1;ops(1,1,4)=-1
 ops(:,4,3)=[.5d0,.5d0,0d0];ops(:,4,4)=ops(:,4,3)
 call symmetry_init(stars,[4,6,2],h,ks,ws,ops,ops,[3,1,1],e)
 if(e/=0.or.stars%nfull/=3.or.stars%max_little/=4)stop 12
 input=1
 ! Locate the mirror+centering transform for the -1/3 representative.
 do r=stars%first(2),stars%first(3)-1
 do q=1,stars%multiplicity(r)
 if(stars%operations(q,r)/=4)cycle
 call symmetry_transform(stars,r,q,input,output,e)
 if(e/=0)stop 13
 ! Mx swaps spin and the nonsymmorphic translation adds exp(-i Rk.t).
 phase=exp(cmplx(0d0,-pi/3,8))*(0d0,-1d0)
 if(maxval(abs(output-phase))>1d-12)stop 14
 enddo
 enddo
 end block
 print *, 'spinor symmetry PASS' 
contains
 subroutine transpose_tiles(send,recv,count,comm)
 complex(8),intent(in)::send(:)
 complex(8),intent(out)::recv(:)
 integer,intent(in)::count,comm
 recv=send
 end subroutine
end program
