program probe
 use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
 use mpi, only: MPI_Init,MPI_Comm_rank,MPI_Comm_size,MPI_Abort,MPI_Finalize,MPI_Alltoall, &
  MPI_Reduce,MPI_COMM_WORLD,MPI_DOUBLE_COMPLEX,MPI_SUM
 use exx_k_exchange, only: exx_k_kernel,exx_k_kernel_init
 use exx_spinor, only: exx_spinor_apply_distributed
 implicit none
 type(exx_k_kernel)::op,full_op
 integer::rank,np,err,status,first,nlocal,i,j,q
 integer,allocatable::starts(:),counts(:)
 complex(8)::u(24,2,6),x(24,2,6),reference(24,2,6),assembled(24,2,6),local_full(24,2,6)
 complex(8),allocatable::result(:,:,:)
 real(8)::k(3,6),pi
 character(16)::invalid,layout_name
 call get_command_argument(1,invalid)
 layout_name='strided'
 if(trim(invalid)=='contiguous')layout_name='contiguous'
 call MPI_Init(status)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
 call MPI_Comm_size(MPI_COMM_WORLD,np,status)
 if(np>8)call MPI_Abort(MPI_COMM_WORLD,1,status)
 allocate(starts(np),counts(np));counts=6/np
 do i=1,mod(6,np);counts(i)=counts(i)+1;enddo
 do i=1,np;starts(i)=1+sum(counts(:i-1));enddo
 first=starts(rank+1);nlocal=counts(rank+1)
 allocate(result(24,2,nlocal))
 pi=acos(-1d0);k=0d0
 do q=1,6
  k(1,q)=2*pi*(mod(q-1,3)+.5d0)/(3*2*.7d0)
  k(2,q)=2*pi*((q-1)/3+.5d0)/(2*3*.7d0)
 enddo
 do q=1,6;do j=1,2;do i=1,24
  u(i,j,q)=cmplx(sin(real(i+3*j+q,8)),cos(real(2*i-j+q,8)),8)/5d0
  x(i,j,q)=cmplx(cos(real(i-2*j+q,8)),sin(real(3*i+j+q,8)),8)/4d0
 enddo;enddo;enddo
 if(rank==0)then
  call exx_k_kernel_init(full_op,[2,3,2],[3,2,1],[.7d0,.7d0,.7d0],k,.11d0,1,err,fft_layout='strided')
  if(err/=0)call MPI_Abort(MPI_COMM_WORLD,2,status)
  call exx_spinor_apply_distributed(full_op,u,x,reference,[1],[6],0,0,transpose_local,err)
  if(err/=0)call MPI_Abort(MPI_COMM_WORLD,3,status)
 endif
 call exx_k_kernel_init(op,[2,3,2],[3,2,1],[.7d0,.7d0,.7d0],k,.11d0,1,err, &
  phase_start=first,phase_count=nlocal,fft_layout=trim(layout_name),defer_fft=.true.)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,4,status)
 if(trim(invalid)=='overflow'.and.rank==0)op%block=huge(0)
 if(trim(invalid)=='nan'.and.rank==0)u(9,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
 call exx_spinor_apply_distributed(op,u(:,:,first:first+nlocal-1),x(:,:,first:first+nlocal-1), &
  result,starts,counts,rank,MPI_COMM_WORLD,transpose_tiles,err,pair_tiles=trim(invalid)/='mixed'.or.rank==0)
 if(trim(invalid)=='nan'.or.trim(invalid)=='mixed'.or.trim(invalid)=='overflow')then
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
  print *, 'SPINOR_PAIR_DISTRIBUTED_PASS',np,maxval(abs(assembled-reference))
 endif
 call MPI_Finalize(status)
contains
 subroutine transpose_tiles(send,recv,count,comm)
  implicit none
  complex(8),intent(in)::send(:)
  complex(8),intent(out)::recv(:)
  integer,intent(in)::count,comm
  integer::stat
  call MPI_Alltoall(send,count,MPI_DOUBLE_COMPLEX,recv,count,MPI_DOUBLE_COMPLEX,comm,stat)
 end subroutine
 subroutine transpose_local(send,recv,count,comm)
  implicit none
  complex(8),intent(in)::send(:)
  complex(8),intent(out)::recv(:)
  integer,intent(in)::count,comm
  recv=send
 end subroutine
end program
