program scalar_worker_probe
 use iso_fortran_env,only:int64
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf,ieee_is_finite
 use mpi
 use exx_k_exchange
 use exx_symmetry
 use exx_spinor_wave_workers
 !$ use omp_lib,only:omp_set_num_threads
 implicit none
 type(exx_k_kernel)::op
 type(s_exx_symmetry_map)::map
 real(8)::a(3,4,2),kr(3,2),wt(2),pi,diff,totaldiff
 complex(8)::u(8,2,2),x(8,2,2)
 complex(8),allocatable::es(:,:,:),et(:,:,:),local_es(:,:,:),target(:,:,:),ref(:,:,:),action(:,:,:),direct(:,:,:)
 integer::rank,np,err,rc,provided,i,j,q,rep,repcount,phasefirst,nlocal,mode,group
 integer,allocatable::starts(:),counts(:)
 integer(int64)::metrics(5),reference_metrics(5)
 call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,rc)
 if(provided<MPI_THREAD_FUNNELED)call MPI_Abort(MPI_COMM_WORLD,1,rc)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,rc)
 call MPI_Comm_size(MPI_COMM_WORLD,np,rc)
 if(np<2)call MPI_Abort(MPI_COMM_WORLD,2,rc)
 allocate(starts(np),counts(np));starts=4;counts=0;starts(:2)=[1,2];counts(:2)=[1,2]
 rep=min(rank+1,3);repcount=merge(1,0,rank<2);phasefirst=starts(rank+1);nlocal=counts(rank+1)
 allocate(local_es(8,2,nlocal),target(8,2,repcount),ref(8,2,repcount),action(8,2,repcount))
 pi=acos(-1d0);a=0;kr=0;kr(1,2)=-2*pi/(3*2*.7d0);wt=[1d0/3,2d0/3]
 do q=1,2;do i=1,3
 a(i,i,q)=1
 enddo;enddo
 a(1,1,2)=-1
 call symmetry_init(map,[2,2,2],[.7d0,.7d0,.7d0],kr,wt,a,a,[3,1,1],err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,3,rc)
 do q=1,2;do j=1,2;do i=1,8
 u(i,j,q)=cmplx(sin(real(i+3*j+q,8)),cos(real(2*i-j+q,8)),8)/5d0
 x(i,j,q)=cmplx(cos(real(i-2*j+q,8)),sin(real(3*i+j+q,8)),8)/4d0
 enddo;enddo;enddo
 call symmetry_expand(map,1,u,x,es,et,err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,4,rc)
 local_es=es(:,:,phasefirst:phasefirst+nlocal-1)
 if(repcount>0)target(:,:,1)=et(:,:,map%first(rep))
 call exx_k_kernel_init(op,[2,2,2],[3,1,1],[.7d0,.7d0,.7d0],map%full_k,.11d0,2,err, &
 phase_start=phasefirst,phase_count=nlocal,fft_layout='strided')
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,5,rc)
 !$ call omp_set_num_threads(1)
 call apply_wave_workers(op,map,local_es,target,ref,starts,counts,rep,rank,MPI_COMM_WORLD,err, &
 reference_metrics(1),reference_metrics(2),reference_metrics(3),reference_metrics(4),reference_metrics(5),rprime_group_count=1)
 if(err/=0)then
 if(rank==0)print *,'SCALAR_WORKER_REJECTED',err
 call MPI_Abort(MPI_COMM_WORLD,6,rc)
 endif
 !$ call omp_set_num_threads(4)
 call apply_wave_workers(op,map,local_es,target,action,starts,counts,rep,rank,MPI_COMM_WORLD,err, &
 metrics(1),metrics(2),metrics(3),metrics(4),metrics(5),rprime_group_count=1)
 if(err/=0.or.any(metrics/=reference_metrics))call MPI_Abort(MPI_COMM_WORLD,7,rc)
 diff=maxval(abs(action-ref));if(repcount==0)diff=0
 call MPI_Allreduce(diff,totaldiff,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,rc)
 if(.not.ieee_is_finite(totaldiff).or.totaldiff>1d-11)call MPI_Abort(MPI_COMM_WORLD,8,rc)
 allocate(direct(8,2,nlocal))
 call exx_k_kernel_apply_distributed(op,local_es,et(:,:,phasefirst:phasefirst+nlocal-1),direct, &
 starts,counts,rank,MPI_COMM_WORLD,transpose_tiles,err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,10,rc)
 diff=0
 if(repcount>0)diff=maxval(abs(ref(:,:,1)-direct(:,:,1)))
 call MPI_Allreduce(diff,totaldiff,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,rc)
 if(.not.ieee_is_finite(totaldiff).or.totaldiff>1d-10)then
 print *,'SCALAR_DIRECT_DIFF',rank,totaldiff
 call MPI_Abort(MPI_COMM_WORLD,11,rc)
 endif
 do group=2,4,2
 if(modulo(np,group)/=0)cycle
 call apply_wave_workers(op,map,local_es,target,action,starts,counts,rep,rank,MPI_COMM_WORLD,err, &
 metrics(1),metrics(2),metrics(3),metrics(4),metrics(5),rprime_group_count=group)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,12,rc)
 diff=0
 if(repcount>0)diff=maxval(abs(action-ref))
 call MPI_Allreduce(diff,totaldiff,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,rc)
 if(.not.ieee_is_finite(totaldiff).or.totaldiff>1d-10)call MPI_Abort(MPI_COMM_WORLD,13,rc)
 enddo
 ! Existing worker collective finite guards must still reject source NaN and Inf.
 do mode=1,2
 if(rank==0)then
 if(mode==1)local_es(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
 if(mode==2)local_es(1,1,1)=cmplx(0d0,ieee_value(0d0,ieee_positive_inf),8)
 endif
 call apply_wave_workers(op,map,local_es,target,action,starts,counts,rep,rank,MPI_COMM_WORLD,err, &
 metrics(1),metrics(2),metrics(3),metrics(4),metrics(5))
 if(err==0)call MPI_Abort(MPI_COMM_WORLD,9,rc)
 enddo
 if(rank==0)print *,'SCALAR_WORKER_OMP_ACTION_PASS',np,totaldiff,' FINITE_REJECTION_PASS'
 call exx_k_kernel_destroy(op)
 call MPI_Finalize(rc)
contains
 subroutine transpose_tiles(send,recv,count,comm)
 complex(8),intent(in)::send(:)
 complex(8),intent(out)::recv(:)
 integer,intent(in)::count,comm
 integer::e
 call MPI_Alltoall(send,count,MPI_DOUBLE_COMPLEX,recv,count,MPI_DOUBLE_COMPLEX,comm,e)
 end subroutine
end program
