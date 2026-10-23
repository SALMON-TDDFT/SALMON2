program action_probe
 use iso_fortran_env,only:int64
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_value,ieee_quiet_nan,ieee_positive_inf
 use mpi
 use exx_k_exchange
 use exx_k_pairs
 use exx_k_orbits
 use exx_symmetry
 !$ use omp_lib,only:omp_set_num_threads
 implicit none
 type(exx_k_kernel)::op
 type(s_exx_symmetry_map)::map
 real(8)::ops(3,4,16),k(3,8),weights(8),pi,diff,globaldiff,thread_diff
 complex(8)::u(8,2,8),x(8,3,8)
 complex(8),allocatable::es(:,:,:),et(:,:,:),su(:,:,:),tx(:,:,:),ref(:,:,:),pair(:,:,:),orbit(:,:,:),old(:,:,:)
 integer,allocatable::starts(:),counts(:)
 integer::rank,np,rc,provided,s,d,j,g,i,p,first,nloc,err,pass,layout,threads,wanted,mode,badrank
 integer(int64)::expected_pairs,total_pairs
 character(10)::layout_name
 call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,rc)
 if(provided<MPI_THREAD_FUNNELED)call MPI_Abort(MPI_COMM_WORLD,1,rc)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,rc);call MPI_Comm_size(MPI_COMM_WORLD,np,rc)
 allocate(starts(np),counts(np))
 do p=0,np-1
 starts(p+1)=1+(p*8)/np;counts(p+1)=((p+1)*8)/np-(p*8)/np
 enddo
 first=starts(rank+1);nloc=counts(rank+1);badrank=0
 do p=1,np
 if(counts(p)>0)then
 badrank=p-1;exit
 endif
 enddo
 pi=acos(-1d0);ops=0;weights=1d0/8
 do s=1,16
 do d=1,3
 ops(d,d,s)=1-2*ibits(s-1,d-1,1)
 enddo
 ops(1,4,s)=.5d0*ibits(s-1,3,1)
 enddo
 do j=1,8
 do d=1,3
 k(d,j)=pi/2*ibits(j-1,d-1,1)
 enddo
 do i=1,2;do g=1,8
 u(g,i,j)=cmplx(sin(real(g+3*i+7*j,8)),cos(real(3*g+7*i+j,8)),8)/5
 enddo;enddo
 do i=1,3;do g=1,8
 x(g,i,j)=cmplx(cos(real(g+5*i+2*j,8)),sin(real(2*g+i+3*j,8)),8)/4
 enddo;enddo
 enddo
 call symmetry_init(map,[2,2,2],[1d0,1d0,1d0],k,weights,ops,ops,2,err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,2,rc)
 call symmetry_expand(map,1,u,x,es,et,err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,3,rc)
 su=es(:,:,first:first+nloc-1);tx=et(:,:,first:first+nloc-1)
 allocate(ref(8,3,nloc),pair(8,3,nloc),orbit(8,3,nloc),old(8,3,nloc))
 globaldiff=0;thread_diff=0
 do layout=1,2
 layout_name='strided'
 if(layout==2)layout_name='contiguous'
 do pass=1,2
 wanted=1
 if(pass==2)wanted=4
 !$ call omp_set_num_threads(wanted)
 call exx_k_kernel_init(op,2,2,1d0,map%full_k,.11d0,2,err,phase_start=first,phase_count=nloc,fft_layout=layout_name)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,4,rc)
 call exx_k_kernel_apply_distributed(op,su,tx,ref,starts,counts,rank,MPI_COMM_WORLD,transpose_tiles,err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,5,rc)
 call apply_k_pairs(op,su,tx,pair,starts,counts,rank,MPI_COMM_WORLD,transpose_tiles,err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,6,rc)
 if(op%pair_fft_pairs>0.and.op%threads_used/=wanted)call MPI_Abort(MPI_COMM_WORLD,6,rc)
 expected_pairs=op%pair_fft_pairs
 call MPI_Allreduce(expected_pairs,total_pairs,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,rc)
 if(total_pairs/=36_int64)call MPI_Abort(MPI_COMM_WORLD,7,rc)
 call apply_k_orbits(op,map,tx,orbit,starts,counts,first,rank,MPI_COMM_WORLD,transpose_tiles,fill_density,err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,8,rc)
 diff=0
 if(nloc>0)diff=max(maxval(abs(ref-pair)),maxval(abs(ref-orbit)))
 call MPI_Allreduce(diff,globaldiff,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,rc)
 if(.not.ieee_is_finite(globaldiff).or.globaldiff>1d-10)then
 print *, 'ACTION_DIFF',layout,pass,globaldiff
 call MPI_Abort(MPI_COMM_WORLD,9,rc)
 endif
 if(pass==1)then
 old=pair
 else
 diff=0
 if(nloc>0)diff=maxval(abs(pair-old))
 call MPI_Allreduce(diff,thread_diff,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,rc)
 if(.not.ieee_is_finite(thread_diff).or.thread_diff>1d-10)call MPI_Abort(MPI_COMM_WORLD,10,rc)
 endif
 do mode=1,2
 if(rank==badrank)then
 if(mode==1)su(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
 if(mode==2)su(1,1,1)=cmplx(0d0,ieee_value(0d0,ieee_positive_inf),8)
 endif
 call apply_k_pairs(op,su,tx,pair,starts,counts,rank,MPI_COMM_WORLD,transpose_tiles,err)
 if(err==0)call MPI_Abort(MPI_COMM_WORLD,11,rc)
 enddo
 su=es(:,:,first:first+nloc-1)
 call exx_k_kernel_destroy(op)
 enddo
 ! Reduced-q residue classes use the same storage, with independent smaller FFTs.
 do pass=1,2
 wanted=1
 if(pass==2)wanted=4
 !$ call omp_set_num_threads(wanted)
 call exx_k_kernel_init(op,[2,2,2],[2,2,2],[1d0,1d0,1d0],map%full_k,.11d0,2,err, &
 phase_start=first,phase_count=nloc,fft_layout=layout_name,q_reduction=[2,1,1])
 if(err/=0.or.op%nq/=4)call MPI_Abort(MPI_COMM_WORLD,12,rc)
 call exx_k_kernel_apply_distributed(op,su,tx,ref,starts,counts,rank,MPI_COMM_WORLD,transpose_tiles,err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,13,rc)
 if(pass==1)then
 old=ref
 else
 diff=0
 if(nloc>0)diff=maxval(abs(ref-old))
 call MPI_Allreduce(diff,thread_diff,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,rc)
 if(.not.ieee_is_finite(thread_diff).or.thread_diff>1d-10)call MPI_Abort(MPI_COMM_WORLD,14,rc)
 endif
 call exx_k_kernel_destroy(op)
 enddo
 enddo
 if(rank==0)print *,'COLLINEAR_ACTION_OMP_PASS',np,globaldiff,thread_diff,' PAIRS36 FINITE_REJECTION_PASS'
 call MPI_Finalize(rc)
contains
 subroutine transpose_tiles(send,recv,count,comm)
 complex(8),intent(in)::send(:)
 complex(8),intent(out)::recv(:)
 integer,intent(in)::count,comm
 integer::status
 call MPI_Alltoall(send,count,MPI_DOUBLE_COMPLEX,recv,count,MPI_DOUBLE_COMPLEX,comm,status)
 end subroutine
 subroutine fill_density(j,points,density)
 integer,intent(in)::j,points(:)
 complex(8),intent(out)::density(:,:)
 integer::r,g
 do g=1,8;do r=1,size(points)
 density(r,g)=sum(su(points(r),:,j)*conjg(su(g,:,j)))*op%phase(points(r),j)*conjg(op%phase(g,j))
 enddo;enddo
 end subroutine
end program
