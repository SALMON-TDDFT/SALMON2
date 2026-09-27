program spatial_driver
 use mpi
 use hse_ace
 use, intrinsic :: ieee_arithmetic
 implicit none
 integer,parameter :: ng=11,no=3,nt=2,nk=2
 complex(8) :: u(ng,no,nk),w(ng,no,nk),t(ng,nt,nk),ref(ng,nt,nk),serial(ng,nt,nk),source_action(ng,no,nk),flag(1,1)
 complex(8),allocatable :: a(:,:,:),ul(:,:,:),wl(:,:,:),tl(:,:,:)
 type(hse_ace_state) :: ace,other,mid,dense
 integer :: rank,np,i,j,k,lo,hi,n,status,err
 real(8),parameter :: dv=.37d0
 real(8) :: difference
 call MPI_Init(err)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,err)
 call MPI_Comm_size(MPI_COMM_WORLD,np,err)
 do k=1,nk
 do j=1,no
 do i=1,ng
 u(i,j,k)=cmplx(sin(real(i*j+k,8)),cos(real(2*i+j*k,8)),8)
 enddo
 enddo
 do j=1,nt
 do i=1,ng
 t(i,j,k)=cmplx(cos(real(i+j,8)),sin(real(i*j+k,8)),8)
 enddo
 enddo
 ! Independent rank-no negative operator K=-U U^dagger, dv=.37.
 w(:,:,k)=-dv*matmul(u(:,:,k),matmul(conjg(transpose(u(:,:,k))),u(:,:,k)))
 ref(:,:,k)=-dv*matmul(u(:,:,k),matmul(conjg(transpose(u(:,:,k))),t(:,:,k)))
 enddo
 call hse_ace_build(dense,u,w,dv,status)
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,12,err)
 call hse_ace_apply(dense,t,serial,status)
 if(status/=0.or.maxval(abs(serial-ref))>1d-10)call MPI_Abort(MPI_COMM_WORLD,13,err)
 call hse_ace_apply(dense,u,source_action,status)
 if(status/=0.or.maxval(abs(source_action-w))>1d-10)call MPI_Abort(MPI_COMM_WORLD,14,err)
 ! With four ranks, last rank deliberately owns no rows.
 n=np
 if(np==4)n=3
 lo=ng*min(rank,n)/n+1;hi=ng*min(rank+1,n)/n
 ul=u(lo:hi,:,:);wl=w(lo:hi,:,:);tl=t(lo:hi,:,:)
 allocate(a(hi-lo+1,nt,nk))
 call hse_ace_build(ace,ul,wl,dv,status,sum_grid)
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,1,err)
 call hse_ace_apply(ace,tl,a,status,sum_grid)
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,2,err)
 difference=0d0
 if(size(a)>0)difference=maxval(abs(a-ref(lo:hi,:,:)))
 flag=cmplx(difference,0d0,8);call sum_grid(flag)
 if(real(flag(1,1))>1d-10)call MPI_Abort(MPI_COMM_WORLD,3,err)
 call hse_ace_build(other,ul,2*wl,dv,status,sum_grid)
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,4,err)
 call hse_ace_average(ace,other,mid,status)
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,5,err)
 call hse_ace_apply(mid,tl,a,status,sum_grid)
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,6,err)
 if(size(a)>0)then
 if(maxval(abs(a-1.5d0*ref(lo:hi,:,:)))>1d-10)call MPI_Abort(MPI_COMM_WORLD,7,err)
 endif
 if(rank==0)tl(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
 call hse_ace_apply(ace,tl,a,status,sum_grid)
 if(status==0)call MPI_Abort(MPI_COMM_WORLD,8,err)
 tl=t(lo:hi,:,:)
 wl=0d0
 call hse_ace_build(ace,ul,wl,dv,status,sum_grid)
 if(status/=0)call MPI_Abort(MPI_COMM_WORLD,9,err)
 call hse_ace_apply(ace,tl,a,status,sum_grid)
 if(status/=0.or.any(a/=(0d0,0d0)))call MPI_Abort(MPI_COMM_WORLD,10,err)
 if(rank==0)ul(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
 call hse_ace_build(ace,ul,wl,dv,status,sum_grid)
 if(status==0)call MPI_Abort(MPI_COMM_WORLD,11,err)
 if(rank==0)print *, 'PASS spatial ACE ranks=',np,' error=',real(flag(1,1))
 call MPI_Finalize(err)
contains
 subroutine sum_grid(x)
 complex(8),intent(inout) :: x(:,:)
 integer :: code
 call MPI_Allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code)
 end subroutine
end program
