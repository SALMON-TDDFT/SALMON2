program probe
 use exx_k_exchange
 use exx_spinor, only: exx_spinor_apply_distributed
 use exx_ace
 implicit none
 type(exx_k_kernel)::op
 type(s_exx_ace)::ace
 integer::err,i,j,a,b,s,t,delta(3)
 complex(8)::u(16,2,1),x(16,2,1),w(16,2,1),ref(16,2,1),wx(16,2,1),wa(16,2,1),ur(16,2,1),xr(16,2,1),wr(16,2,1)
 complex(8)::rot(2,2),dm,scalar_action(8,2,1)
 real(8)::k(3,1)
 k=0d0
 do a=1,2;do i=1,16
  u(i,a,1)=cmplx(sin(real(i+3*a,8)),cos(real(2*i-a,8)),8)/5d0
  x(i,a,1)=cmplx(cos(real(i-2*a,8)),sin(real(3*i+a,8)),8)/4d0
 enddo;enddo
 call exx_k_kernel_init(op,2,1,.7d0,k,.11d0,3,err,fft_layout='strided')
 if(err/=0)stop 1
 call exx_spinor_apply_distributed(op,u,x,w,[1],[1],0,0,transpose_tiles,err)
 if(err/=0)stop 2
 ref=0d0
 do i=1,8;do j=1,8
  delta=modulo(op%point(:,i)-op%point(:,j),op%n)
  do s=1,2;do t=1,2
   dm=sum(u(i+(s-1)*8,:,1)*conjg(u(j+(t-1)*8,:,1)))
   ref(i+(s-1)*8,:,1)=ref(i+(s-1)*8,:,1)-dm*op%kernel(delta(1),delta(2),delta(3))*x(j+(t-1)*8,:,1)
  enddo;enddo
 enddo;enddo
 if(maxval(abs(ref-w))>1d-11)stop 3
 rot(1,:)=[cmplx(.6d0,0d0,8),cmplx(0d0,.8d0,8)]
 rot(2,:)=[cmplx(0d0,.8d0,8),cmplx(.6d0,0d0,8)]
 do s=1,2
  ur((s-1)*8+1:s*8,:,1)=rot(s,1)*u(1:8,:,1)+rot(s,2)*u(9:16,:,1)
  xr((s-1)*8+1:s*8,:,1)=rot(s,1)*x(1:8,:,1)+rot(s,2)*x(9:16,:,1)
 enddo
 call exx_spinor_apply_distributed(op,ur,xr,wr,[1],[1],0,0,transpose_tiles,err)
 if(err/=0)stop 4
 do s=1,2
  ref((s-1)*8+1:s*8,:,1)=rot(s,1)*w(1:8,:,1)+rot(s,2)*w(9:16,:,1)
 enddo
 if(maxval(abs(ref-wr))>1d-11)stop 5
 call exx_spinor_apply_distributed(op,u,u,wx,[1],[1],0,0,transpose_tiles,err)
 if(err/=0)stop 6
 call exx_ace_build(ace,u,wx,.7d0**3,err)
 if(err/=0)stop 7
 call exx_ace_apply(ace,u,wa,err)
 if(err/=0.or.maxval(abs(wa-wx))>1d-10)stop 8
 ! Closed-shell spin pairs reduce to the scalar one-spin kernel exactly.
 u=0d0
 do i=1,8
  u(i,1,1)=cmplx(sin(real(i,8)),cos(real(i,8)),8)/4d0
  u(i+8,2,1)=u(i,1,1)
 enddo
 call exx_spinor_apply_distributed(op,u,x,w,[1],[1],0,0,transpose_tiles,err)
 if(err/=0)stop 9
 call exx_k_kernel_apply(op,u(1:8,1:1,:),x(1:8,:,:),scalar_action,0,1,err)
 if(err/=0.or.maxval(abs(w(1:8,:,:)-scalar_action))>1d-11)stop 10
 call exx_k_kernel_apply(op,u(9:16,2:2,:),x(9:16,:,:),scalar_action,0,1,err)
 if(err/=0.or.maxval(abs(w(9:16,:,:)-scalar_action))>1d-11)stop 11
 print *, 'SPINOR_FOCK_OFFDIAGONAL_SU2_ACE_CLOSED_SHELL_PASS'
contains
 subroutine transpose_tiles(send,recv,count,comm)
  complex(8),intent(in)::send(:)
  complex(8),intent(out)::recv(:)
  integer,intent(in)::count,comm
  recv=send
 end subroutine
end program
