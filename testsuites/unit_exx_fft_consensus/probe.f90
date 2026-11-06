program probe
 use mpi
 use ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_is_finite
 use exx_fft_consensus,only:agree_fft_layout
 use exx_k_exchange,only:exx_k_kernel,exx_k_kernel_init,exx_k_kernel_apply,exx_k_kernel_apply_distributed
 implicit none
 type(exx_k_kernel)::op,reference_op
 integer::rank,np,e,stat,mode,i,j,q,first,nlocal,low,high,flag,reduce(3)
 integer,allocatable::starts(:),counts(:)
 real(8)::k(3,8),pi,err,global_err,h(3)
 complex(8)::u(8,2,8),reference(8,2,8),assembled(8,2,8),local(8,2,8)
 complex(8),allocatable::action(:,:,:)
 character(20)::layout_name
 call MPI_Init(e);call MPI_Comm_rank(MPI_COMM_WORLD,rank,e);call MPI_Comm_size(MPI_COMM_WORLD,np,e)
 allocate(starts(np),counts(np));counts=8/np
 do i=1,mod(8,np);counts(i)=counts(i)+1;enddo
 do i=1,np;starts(i)=1+sum(counts(:i-1));enddo
 first=starts(rank+1);nlocal=counts(rank+1);allocate(action(8,2,nlocal))
 pi=acos(-1d0)
 do q=1,8
 k(:,q)=pi*real([mod(q-1,2),mod((q-1)/2,2),(q-1)/4],8)/1.4d0
 do j=1,2;do i=1,8
 u(i,j,q)=cmplx(sin(.13d0*(i+j+q)),cos(.17d0*(i+2*j+q)),8)
 enddo;enddo
 enddo
 call exx_k_kernel_init(reference_op,2,2,.7d0,k,.11d0,1,stat,fft_layout='strided')
 if(stat/=0)call MPI_Abort(MPI_COMM_WORLD,1,e)
 call exx_k_kernel_apply(reference_op,u,u,reference,0,1,stat)
 if(stat/=0)call MPI_Abort(MPI_COMM_WORLD,2,e)
 do mode=0,8
 layout_name="auto";reduce=1
 if(mode==6.and.rank==np-1)reduce(1)=0
 if(mode==7.and.rank==np-1)layout_name="contiguous"
 if(mode==8.and.rank==np-1)layout_name="invalid"
 h=.7d0
 if(mode==5.and.rank==np-1)h(1)=-.7d0
 call exx_k_kernel_init(op,[2,2,2],[2,2,2],h,k,.11d0,1,stat, &
  phase_start=first,phase_count=nlocal,fft_layout=layout_name,q_reduction=reduce,fft_layout_sync=sync_layout)
 if(mode==2.or.mode==3.or.mode>=5)then
 if(stat==0)call MPI_Abort(MPI_COMM_WORLD,3,e)
 cycle
 endif
 if(stat/=0)call MPI_Abort(MPI_COMM_WORLD,4,e)
 flag=merge(1,0,op%contiguous_fft)
 call MPI_Allreduce(flag,low,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,e)
 call MPI_Allreduce(flag,high,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,e)
 if(low/=high.or.(flag==1.neqv.mode==1))call MPI_Abort(MPI_COMM_WORLD,5,e)
 call exx_k_kernel_apply_distributed(op,u(:,:,first:first+nlocal-1),u(:,:,first:first+nlocal-1),action, &
  starts,counts,rank,MPI_COMM_WORLD,transpose_tiles,stat)
 if(stat/=0)call MPI_Abort(MPI_COMM_WORLD,6,e)
 local=0;local(:,:,first:first+nlocal-1)=action
 call MPI_Allreduce(local,assembled,size(local),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,e)
 err=maxval(abs(assembled-reference));call MPI_Allreduce(err,global_err,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,e)
 if(global_err>1d-11)call MPI_Abort(MPI_COMM_WORLD,7,e)
 if(rank==0)print *,'FFT_CONSENSUS_ORACLE_PASS',mode,global_err
 enddo
 mode=0
 call exx_k_kernel_init(op,2,2,.7d0,k,.11d0,1,stat,fft_layout='auto',fft_layout_sync=sync_layout)
 if(stat/=0.or.op%contiguous_fft)call MPI_Abort(MPI_COMM_WORLD,8,e)
 call exx_k_kernel_init(op,2,2,.7d0,k,.11d0,1,stat,fft_layout='contiguous',fft_layout_sync=sync_layout)
 if(stat/=0.or..not.op%contiguous_fft)call MPI_Abort(MPI_COMM_WORLD,9,e)
 call exx_k_kernel_init(op,2,2,.7d0,k,.11d0,1,stat,defer_fft=.true.,fft_layout_sync=sync_layout)
 if(stat/=0.or.op%auto_fft)call MPI_Abort(MPI_COMM_WORLD,10,e)
 if(rank==0)print *,'FFT_CONSENSUS_COLLECTIVE_REJECTION_PASS'
 call MPI_Finalize(e)
contains
 subroutine sync_layout(times,contiguous,status,automatic)
 real(8),intent(inout)::times(2)
 logical,intent(inout)::contiguous
 logical,intent(in)::automatic
 integer,intent(inout)::status
 ! Opposing preferences model different startup timings without changing physical data.
 times=[2d0,1d0]
 if(rank/=0)times=[1d0,2d0]
 if(mode==1)times=[2d0+.1d0*rank,.5d0]
 if(mode==2.and.rank==np-1)times(2)=ieee_value(0d0,ieee_quiet_nan)
 if(mode==3.and.rank==np-1)status=1
 if(mode==4)times=[1d-8,1d-10]
 call agree_fft_layout(times,contiguous,status,MPI_COMM_WORLD,automatic)
 end subroutine
 subroutine transpose_tiles(send,recv,count,comm)
 complex(8),intent(in)::send(:)
 complex(8),intent(out)::recv(:)
 integer,intent(in)::count,comm
 integer::z
 call MPI_Alltoall(send,count,MPI_DOUBLE_COMPLEX,recv,count,MPI_DOUBLE_COMPLEX,comm,z)
 end subroutine
end program
