program test_q_batch
  use hsex_q_response, only: hsex_q_solve,hsex_q_shared,hsex_q_workspace, &
    hsex_q_prepare,hsex_q_release,hsex_q_worker_prepare,hsex_q_worker_release,hsex_q_one
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: n=384,nb=7,nq=3,mg=80,nu=3,nc=1
  complex(8) :: adj(n,n),u(n,nu),psi(n,nb,nq),phase(n,mg,nq),dummy(1,1)
  complex(8) :: actual(mg,mg,nc,nq),expected(mg,mg,nc,nq),one(mg,mg,nc,nq),b(n,mg),rot(n,mg)
  real(8) :: d(n),e(n-1),levels(nb,nq),occ(nb,nq),eta(nc),residual(nq),pi,error,slow(3),fast(3),tmp
  real(8) :: saved_occ(nb,nq)
  integer :: counts(nq),status(nq),team,i,j,k,v,trial,repetition,pump,threads
  integer(8) :: start,finish,rate
  character(32) :: argument
  pi=acos(-1d0);dummy=0d0;eta=0d0
  do j=1,n
    do i=1,n
      adj(i,j)=sqrt(2d0/n)*cos(pi*(i-.5d0)*(j-1)/n)
    enddo
    adj(:,1)=1d0/sqrt(dble(n))
    d(j)=2d0+.01d0*j
  enddo
  e=0d0;u=0d0
  do j=1,nu
    u(j,j)=1d0
  enddo
  do k=1,nq
    do v=1,nb
      levels(v,k)=-.7d0+.1d0*v+.01d0*k
      occ(v,k)=2d0
      do i=1,n
        psi(i,v,k)=cmplx(cos(dble(i*v+k)),sin(dble(i+v*k)),8)/sqrt(dble(n))
      enddo
    enddo
    do j=1,mg
      do i=1,n
        phase(i,j,k)=exp(cmplx(0d0,2d0*pi*i*j/n,8))
      enddo
    enddo
  enddo
  counts=[mg,mg-7,mg/2];occ(2,2)=0d0;occ(4,3)=.6d0
  do pump=0,1
    eta=0d0
    if(pump==1)eta=.02d0
    expected=0d0
    do v=1,nb
      call hsex_q_solve(n,1,nq,mg,nu,nc,1,pump,adj,d,e,u,psi(:,v:v,:),phase,counts, &
        levels(v:v,:),occ(v:v,:),eta,one,residual,status,team,0,dummy)
      if(any(status/=0).or.team/=1)error stop 'single-band reference failed'
      expected=expected+one
    enddo
    do threads=1,2
      call hsex_q_solve(n,nb,nq,mg,nu,nc,threads,pump,adj,d,e,u,psi,phase,counts, &
        levels,occ,eta,actual,residual,status,team,0,dummy)
      if(any(status/=0).or.team/=threads)error stop 'batched response rejected'
      error=0d0
      do k=1,nq
        do j=1,mg
          do i=1,mg
            if(.not.ieee_is_finite(real(actual(i,j,1,k),8)))error stop 'nonfinite batched real response'
            if(.not.ieee_is_finite(aimag(actual(i,j,1,k))))error stop 'nonfinite batched imaginary response'
            error=max(error,abs(actual(i,j,1,k)-expected(i,j,1,k)))
          enddo
        enddo
      enddo
      print *, 'Q_BATCH_ERROR pump/threads/error',pump,threads,error
      if(error>2d-12)error stop 'band batching changes response'
    enddo
  enddo
  ! Independent explicit diagonal resolvent checks the static projection,
  ! occupation prefactor and target transform, rather than only self-comparison.
  expected=0d0
  do k=1,nq
    do v=1,nb
      b=phase(:,:,k)*spread(psi(:,v,k),2,mg);b(:nu,:)=0d0
      rot=matmul(adj,b)
      do j=1,counts(k)
        do i=1,counts(k)
          expected(i,j,1,k)=expected(i,j,1,k)+.5d0*occ(v,k)* &
            sum(conjg(rot(:,i))*rot(:,j)/(d-levels(v,k)))
        enddo
      enddo
    enddo
  enddo
  eta=0d0
  call hsex_q_solve(n,nb,nq,mg,nu,nc,1,0,adj,d,e,u,psi,phase,counts, &
    levels,occ,eta,actual,residual,status,team,0,dummy)
  if(any(status/=0).or.maxval(abs(actual-expected))>2d-12)error stop 'independent static resolvent differs'
  block
    type(hsex_q_shared) :: shared
    type(hsex_q_workspace) :: workspace
    integer :: iteration,q_index,code,worker_error,workers
    expected=actual
    call hsex_q_prepare(n,nb,mg,nu,0,0,eta,dummy,shared,code)
    if(code/=0)error stop 'persistent target setup failed'
    worker_error=0;workers=0
    !$omp parallel num_threads(2) private(workspace,iteration,q_index,code) reduction(max:worker_error) &
    !$omp& reduction(+:workers)
    call hsex_q_worker_prepare(n,nb,mg,nu,0,0,shared,workspace)
    workers=workers+1
    do iteration=1,3
      !$omp do schedule(dynamic,1)
      do q_index=1,nq
        call hsex_q_one(n,nb,mg,nu,nc,0,0,adj,d,e,u,psi(:,:,q_index),phase(:,:,q_index),counts(q_index), &
          levels(:,q_index),occ(:,q_index),eta,shared,workspace,actual(:,:,:,q_index),residual(q_index),status(q_index))
      enddo
      !$omp end do
    enddo
    call hsex_q_worker_release(workspace,code)
    worker_error=max(worker_error,code)
    !$omp end parallel
    call hsex_q_release(shared)
    if(worker_error/=0.or.workers/=2.or.any(status/=0))error stop 'persistent workers rejected'
    if(maxval(abs(actual-expected))>2d-12)error stop 'persistent workers change response'
    print *, 'Q_PERSISTENT_WORKERS_PASS batches/workspaces',3,workers
  end block
  saved_occ=occ;occ=0d0
  call hsex_q_solve(n,nb,nq,mg,nu,nc,2,0,adj,d,e,u,psi,phase,counts, &
    levels,occ,eta,actual,residual,status,team,0,dummy)
  if(any(status/=0).or.team/=2)error stop 'empty source set rejected'
  if(any(actual/=cmplx(0d0,0d0,8)).or.any(residual/=0d0))error stop 'empty source set gives nonzero response'
  occ=saved_occ
  call get_command_argument(1,argument)
  if(argument=='--benchmark')then
    do trial=1,3
      call system_clock(start,rate)
      do repetition=1,5
        do v=1,nb
          call hsex_q_solve(n,1,nq,mg,nu,nc,1,0,adj,d,e,u,psi(:,v:v,:),phase,counts, &
            levels(v:v,:),occ(v:v,:),eta,one,residual,status,team,0,dummy)
          if(any(status/=0))error stop 'reference benchmark rejected'
        enddo
      enddo
      call system_clock(finish);slow(trial)=dble(finish-start)/rate
      call system_clock(start)
      do repetition=1,5
        call hsex_q_solve(n,nb,nq,mg,nu,nc,1,0,adj,d,e,u,psi,phase,counts, &
          levels,occ,eta,actual,residual,status,team,0,dummy)
        if(any(status/=0))error stop 'batch benchmark rejected'
      enddo
      call system_clock(finish);fast(trial)=dble(finish-start)/rate
    enddo
    do i=1,2
      do j=i+1,3
        if(slow(j)<slow(i))then
          tmp=slow(i);slow(i)=slow(j);slow(j)=tmp
        endif
        if(fast(j)<fast(i))then
          tmp=fast(i);fast(i)=fast(j);fast(j)=tmp
        endif
      enddo
    enddo
    print *, 'Q_BATCH_MEDIAN single/batched/speedup',slow(2),fast(2),slow(2)/fast(2)
  endif
  print *, 'Q_BATCH_PASS'
end program test_q_batch
