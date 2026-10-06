program benchmark
 use iso_fortran_env,only:int64
 use exx_wannier
 implicit none
 integer,parameter :: n(3)=[32,24,20],ng=15360,no=4,nt=8
 type(s_exx_wannier) :: op
 complex(8),allocatable :: target(:,:,:),full(:,:,:),local(:,:,:)
 real(8) :: h(3)=[.7d0,.8d0,.9d0],k(3,1),d(3),t0,t1,t2,error
 integer(int64) :: tick0,tick1,tick2,rate
 integer :: g,j,trial,status,center(3)
 k=0
 call wannier_init(op,n,[1,1,1],h,k,0d0,status)
 if(status/=0)error stop 'init'
 allocate(op%source(ng,no),target(ng,nt,1),full(ng,nt,1),local(ng,nt,1))
 op%source=0
 do j=1,no
 center=modulo([7*j,5*j,3*j],n)
 if(j==1)center=n-1 ! support crosses all three periodic boundaries
 do g=1,ng
 d=modulo(real(op%point(:,g)-center,8)+n/2d0,real(n,8))-n/2d0
 if(any(abs(d)>2))cycle
 op%source(g,j)=exp(-sum(d*d)/3d0)*cmplx(1d0,.1d0*j,8)
 enddo
 enddo
 do j=1,nt;do g=1,ng
 target(g,j,1)=cmplx(sin(.003d0*g*j),cos(.007d0*(g+j)),8)
 enddo;enddo
 target(:,nt,1)=0 ! exact zero pairs remain skipped
 call wannier_apply(op,target,full,status)
 if(status/=0)error stop 'full'
 op%use_local_fft=.true.
 call wannier_apply(op,target,local,status)
 if(status/=0)error stop 'local'
 error=maxval(abs(full-local))/maxval(abs(full))
 if(error>2d-12)error stop 'full/local mismatch'
 if(op%local_fft_pairs_executed/=no*(nt-1))error stop 'local/zero-pair count'
 if(op%fft_pair_grid_points>=op%fft_pairs_executed*op%ngs/10)error stop 'insufficient volume reduction'
 print *, 'relative_action_error',error
 print *, 'local_pairs',op%local_fft_pairs_executed
 print *, 'local_pair_grid_points',op%fft_pair_grid_points
 print *, 'full_pair_grid_points',op%fft_pairs_executed*op%ngs
 call system_clock(tick0,rate)
 op%use_local_fft=.false.
 do trial=1,5
 call wannier_apply(op,target,full,status)
 if(status/=0)error stop 'timed full'
 enddo
 call system_clock(tick1)
 op%use_local_fft=.true.
 do trial=1,5
 call wannier_apply(op,target,local,status)
 if(status/=0)error stop 'timed local'
 enddo
 call system_clock(tick2)
 print *, 'full_wall_seconds_per_action',real(tick1-tick0,8)/rate/5
 print *, 'local_wall_seconds_per_action',real(tick2-tick1,8)/rate/5
 call wannier_destroy(op)
end program
