program probe
 use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
 use mpi
 use exx_k_exchange
 use exx_spinor
 use exx_symmetry
 implicit none
 type(exx_k_kernel)::op,full_op
 type(s_exx_symmetry_map)::map
 real(8)::a(3,4,2),kr(3,2),wt(2)
 complex(8),allocatable::es(:,:,:),et(:,:,:)
 integer::rank,np,err,status,first,nlocal,i,j,q
 integer,allocatable::starts(:),counts(:)
 complex(8)::u(16,4,3),x(16,2,3),reference(16,2,3),assembled(16,2,3),local_full(16,2,3)
 complex(8),allocatable::result(:,:,:)
 real(8)::k(3,3),pi
 character(16)::invalid
 call get_command_argument(1,invalid)
 call MPI_Init(status)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
 call MPI_Comm_size(MPI_COMM_WORLD,np,status)
 if(np/=2)call MPI_Abort(MPI_COMM_WORLD,1,status)
 allocate(starts(np),counts(np));counts=[1,2];starts=[1,2]
 first=starts(rank+1);nlocal=counts(rank+1)
 allocate(result(16,2,nlocal))
 pi=acos(-1d0);a=0;kr=0;kr(1,2)=-2*pi/(3*2*.7d0);wt=[1d0/3,2d0/3]
 do q=1,2
 do i=1,3
 a(i,i,q)=1
 enddo
 enddo
 a(1,1,2)=-1
 call symmetry_init(map,[2,2,2],[.7d0,.7d0,.7d0],kr,wt,a,a,[3,1,1],err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,8,status)
 k=map%full_k
 do q=1,2;do j=1,2;do i=1,16
  u(i,j,q)=cmplx(sin(real(i+3*j+q,8)),cos(real(2*i-j+q,8)),8)/5d0
  x(i,j,q)=cmplx(cos(real(i-2*j+q,8)),sin(real(3*i+j+q,8)),8)/4d0
 enddo;enddo;enddo
 call symmetry_expand(map,1,u(:,:2,:2),x(:,:,:2),es,et,err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,9,status)
 ! Explicit complete star inputs carry little-group projector weights.
 u=es;x=et
 if(rank==0)then
  call exx_k_kernel_init(full_op,[2,2,2],[3,1,1],[.7d0,.7d0,.7d0],k,.11d0,3,err,fft_layout='strided')
  if(err/=0)call MPI_Abort(MPI_COMM_WORLD,2,status)
  call exx_spinor_apply_distributed(full_op,u,x,reference,[1],[3],0,0,transpose_local,err)
  if(err/=0)call MPI_Abort(MPI_COMM_WORLD,3,status)
 endif
 call exx_k_kernel_init(op,[2,2,2],[3,1,1],[.7d0,.7d0,.7d0],k,.11d0,3,err, &
  phase_start=first,phase_count=nlocal,fft_layout='strided')
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,4,status)
 if(trim(invalid)=='nan'.and.rank==0)u(9,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
 call exx_spinor_apply_distributed(op,u(:,:,first:first+nlocal-1),x(:,:,first:first+nlocal-1), &
  result,starts,counts,rank,MPI_COMM_WORLD,transpose_tiles,err)
 if(trim(invalid)=='nan')then
  if(err==0)call MPI_Abort(MPI_COMM_WORLD,7,status)
  if(rank==0)print *, 'SPINOR_INVALID_CROSS_SPIN_COLLECTIVE_PASS'
  call MPI_Finalize(status)
  stop
 endif
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,5,status)
 local_full=0d0;local_full(:,:,first:first+nlocal-1)=result
 call MPI_Reduce(local_full,assembled,size(assembled),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,status)
 if(rank==0)then
  if(maxval(abs(assembled-reference))>1d-11)call MPI_Abort(MPI_COMM_WORLD,6,status)
  print *, 'SPINOR_SYMMETRY_UNEVEN_DISTRIBUTED_PASS',np,maxval(abs(assembled-reference))
 endif
 call MPI_Finalize(status)
contains
 subroutine transpose_tiles(send,recv,count,comm)
  complex(8),intent(in)::send(:)
  complex(8),intent(out)::recv(:)
  integer,intent(in)::count,comm
  integer::stat
  call MPI_Alltoall(send,count,MPI_DOUBLE_COMPLEX,recv,count,MPI_DOUBLE_COMPLEX,comm,stat)
 end subroutine
 subroutine transpose_local(send,recv,count,comm)
  complex(8),intent(in)::send(:)
  complex(8),intent(out)::recv(:)
  integer,intent(in)::count,comm
  recv=send
 end subroutine
end program
