program probe
  use hse_wannier
  use hse_wannier_gauge
  use hse_ace
  implicit none
  type(s_hse_wannier) :: op
  type(hse_ace_state) :: ace
  integer :: n(3),mesh(3),no,nt,ng,nk,iu,ierr,i,j,ik,expected_workers
  real(8) :: h(3),omega,checks(5),spread,gradient,smin
  real(8),allocatable :: k(:,:),occ(:,:),trial_occ(:,:),eval(:),rwork(:)
  complex(8),allocatable :: psi(:,:,:),target(:,:,:),action(:,:,:),w(:,:,:),wa(:,:,:),back(:,:,:)
  complex(8),allocatable :: transported_action(:,:,:)
  complex(8),allocatable :: gauge(:,:,:),rot(:,:,:),previous(:,:,:),transported(:,:,:),metric(:,:),work(:)
  character(1024) :: path,out,worker_setting
  call check_gauge_minimizer()
  call get_command_argument(1,path)
  call get_command_argument(2,out)
  open(newunit=iu,file=trim(path),form='unformatted',access='stream')
  read(iu)n,mesh,no,nt
  ng=product(n);nk=product(mesh)
  allocate(k(3,nk),occ(no,nk),psi(ng,no,nk),target(ng,nt,nk),action(ng,nt,nk),w(ng,no,nk),wa(ng,no,nk))
  read(iu)h,omega,k,occ,psi,target
  close(iu)
  call wannier_init(op,n,mesh,h,k,omega,ierr)
  if(ierr/=0)error stop 'init'
  call get_environment_variable('WANNIER_TEST_BATCH',worker_setting,status=ierr)
  if(ierr==0)read(worker_setting,*)op%fft_batch_size
  allocate(gauge(no,no,nk),rot(no,no,nk),previous(ng,no,nk),transported(no,no,nk))
  gauge=0d0
  do ik=1,nk
    do i=1,no
      gauge(i,i,ik)=exp(cmplx(0d0,0.17d0*i*ik,8))
    enddo
    previous(:,:,ik)=matmul(psi(:,:,ik),gauge(:,:,ik))
    rot(:,:,ik)=0d0
    do i=1,no
      rot(i,no+1-i,ik)=exp(cmplx(0d0,0.29d0*i*ik,8))
    enddo
    w(:,:,ik)=matmul(psi(:,:,ik),rot(:,:,ik))
  enddo
  call gauge_transport(w,previous,product(h),transported,smin,ierr)
  if(ierr/=0)error stop 'transport'
  checks=0d0
  do ik=1,nk
    checks(1)=max(checks(1),maxval(abs(matmul(w(:,:,ik),transported(:,:,ik))-previous(:,:,ik))))
  enddo
  call wannier_set_source(op,psi,occ,gauge,ierr)
  if(ierr/=0)error stop 'source'
  call wannier_apply(op,target,action,ierr)
  if(ierr/=0)error stop 'action'
  call get_environment_variable('WANNIER_EXPECT_WORKERS',worker_setting,status=ierr)
  if(ierr==0)then
    read(worker_setting,*)expected_workers
    if(op%workers/=expected_workers)error stop 'wrong FFT worker count'
  endif
  call wannier_apply(op,psi,w,ierr)
  if(ierr/=0)error stop 'construction action'
  allocate(metric(no,no),eval(no),work(4*no),rwork(3*no))
  do ik=1,nk
    metric=matmul(conjg(transpose(psi(:,:,ik))),w(:,:,ik))*product(h)
    checks(2)=max(checks(2),maxval(abs(metric-conjg(transpose(metric)))))
    call zheev('V','U',no,metric,no,eval,work,size(work),rwork,ierr)
    if(ierr/=0)error stop 'metric'
    checks(5)=max(checks(5),maxval(eval))
  enddo
  call hse_ace_build(ace,psi,w,product(h),ierr)
  if(ierr/=0)error stop 'ACE build'
  call hse_ace_apply(ace,psi,wa,ierr)
  if(ierr/=0)error stop 'ACE apply'
  checks(3)=maxval(abs(wa-w))
  allocate(back(ng*nk,no,1))
  call wannier_forward(op,psi,back(:,:,1))
  call wannier_backward(op,back(:,:,1),wa)
  checks(4)=maxval(abs(wa-psi))
  open(newunit=iu,file=trim(out),form='unformatted',access='stream',status='replace')
  write(iu)action
  close(iu)
  open(newunit=iu,file=trim(out)//'.checks',status='replace')
  write(iu,'(es24.15)')checks
  close(iu)
  ! Localization is an acceleration gauge: full-support exchange is invariant.
  call wannier_localize(op,psi,100,1d-7,ierr)
  if(ierr/=0)error stop 'localization setup'
  call wannier_set_source(op,psi,occ,op%gauge,ierr)
  if(ierr/=0)error stop 'localized source'
  call wannier_apply(op,target,action,ierr)
  if(ierr/=0)error stop 'localized action'
  open(newunit=iu,file=trim(out)//'.localized',form='unformatted',access='stream',status='replace')
  write(iu)action
  close(iu)
  call wannier_snapshot(op,occ,omega,0d0,1d-9,7,.true.,trim(out)//'.snapshot',ierr)
  if(ierr/=0)error stop 'snapshot write'
  call wannier_localize(op,psi,0,1d-7,ierr)
  if(ierr/=0.or.op%localization_status/=2)error stop 'transport-only gauge status'
  call wannier_set_source(op,psi,occ,op%gauge,ierr)
  if(ierr/=0)error stop 'transport-only source'
  allocate(transported_action,mold=action)
  call wannier_apply(op,target(:,1:0,:),transported_action(:,1:0,:),ierr)
  if(ierr==0)error stop 'empty target accepted'
  call wannier_apply(op,target,transported_action,ierr)
  if(ierr/=0.or.maxval(abs(transported_action-action))>1d-10)error stop 'transport-only action changed'
  call wannier_refresh_source(op,psi,occ,0,1d-7,ierr)
  if(ierr/=0)error stop 'active source refresh'
  if(size(op%source,2)/=max(1,count(any(occ>0d0,dim=2))))error stop 'active source count'
  call wannier_apply(op,target,transported_action,ierr)
  if(ierr/=0.or.maxval(abs(transported_action-action))>1d-10)error stop 'active source action changed'
  ! Exercise source-set changes with and without a change in source count.
  trial_occ=occ
  do j=1,no
    trial_occ=0d0;trial_occ(j,:)=1d0
    call wannier_refresh_source(op,psi,trial_occ,0,1d-7,ierr)
    if(ierr/=0.or.size(op%source,2)/=1)error stop 'changed source set'
    call wannier_apply(op,target,transported_action,ierr)
    if(ierr/=0)error stop 'changed source action'
    call wannier_set_source(op,psi,trial_occ,gauge,ierr)
    if(ierr/=0)error stop 'full reference source'
    call wannier_apply(op,target,action,ierr)
    if(ierr/=0.or.maxval(abs(transported_action-action))>1d-10)error stop 'changed source mismatch'
  enddo
  trial_occ=0d0
  call wannier_refresh_source(op,psi,trial_occ,0,1d-7,ierr)
  if(ierr/=0)error stop 'zero source refresh'
  call wannier_apply(op,target,transported_action,ierr)
  if(ierr/=0.or.maxval(abs(transported_action))>0d0)error stop 'zero source action'
  call wannier_refresh_source(op,psi,occ,0,1d-7,ierr)
  if(ierr/=0)error stop 'source reactivation'
  call wannier_apply(op,target,transported_action,ierr)
  call wannier_set_source(op,psi,occ,gauge,ierr)
  call wannier_apply(op,target,action,ierr)
  if(ierr/=0.or.maxval(abs(transported_action-action))>1d-10)error stop 'reactivated source mismatch'
  w=0d0
  call hse_ace_build(ace,psi,w,product(h),ierr)
  if(ierr/=0)error stop 'zero exchange ACE build'
  call hse_ace_apply(ace,psi,wa,ierr)
  if(ierr/=0.or.maxval(abs(wa))>0d0)error stop 'zero exchange ACE apply'
  call wannier_destroy(op)
contains
  subroutine check_gauge_minimizer()
    implicit none
    complex(8) :: u(2,2,3),raw(2,2,2,3),d(2,2,3),dp(2,2,3),up(2,2,3)
    real(8) :: b(3,2),weights(2),f,variable,fp,fm,dummy,eps,norm,initial
    integer :: neighbors(2,3),k,i,status,iterations
    neighbors(:,1)=[2,3];neighbors(:,2)=[3,1];neighbors(:,3)=[1,2]
    b=0d0;b(1,:)=[1d0,-1d0];weights=.5d0;raw=0d0;u=0d0
    do k=1,3
      u(1,1,k)=cos(.21d0*k);u(2,2,k)=u(1,1,k)
      u(1,2,k)=sin(.21d0*k);u(2,1,k)=-u(1,2,k)
      do i=1,2
        raw(i,i,:,k)=.8d0
      enddo
    enddo
    call gauge_functional(u,raw,neighbors,b,weights,f,variable,d,status)
    if(status/=0)error stop 'gauge functional'
    initial=f;eps=1d-6
    do k=1,3
      up(:,:,k)=u(:,:,k)+eps*matmul(u(:,:,k),d(:,:,k))
    enddo
    call gauge_functional(up,raw,neighbors,b,weights,fp,dummy,dp,status)
    do k=1,3
      up(:,:,k)=u(:,:,k)-eps*matmul(u(:,:,k),d(:,:,k))
    enddo
    call gauge_functional(up,raw,neighbors,b,weights,fm,dummy,dp,status)
    if(abs((fp-fm)/(2*eps)+sum(abs(d)**2))>1d-7)error stop 'MV gradient'
    call gauge_minimize(u,raw,neighbors,b,weights,500,1d-7,f,norm,iterations,status)
    if(status/=0.or.f>=initial.or.norm>1d-7)error stop 'MV convergence'
  end subroutine
end program
