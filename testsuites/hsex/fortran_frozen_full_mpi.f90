program test_frozen_full
  use mpi
  use hsex_frozen_full, only: frozen_full_w,s_frozen_window_cache,frozen_window_anchor,frozen_window_reset
  use hsex_screen_math, only: screen_geometry,screen_assemble,screen_transform
  use hsex_longwave, only: longwave_cell
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  integer,parameter :: n(3)=[6,2,2],ng=24,nk=2,na=8
  integer :: rank,np,ierr,provided,ks,ke,first,count,status,j,p,boxes,action_calls,nqproc
  real(8) :: cell(3,3),k(3,nk),closure_tol,eta
  real(8),allocatable :: w(:,:),part(:,:),total(:,:),reference_w(:,:)
  real(8) :: expected(na*na,nk)
  complex(8),allocatable :: psi(:,:,:)
  logical :: fail_window=.false.
  integer :: filter,previous_calls
  type(s_frozen_window_cache) :: cache
  character(32) :: mode
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  ks=1+nk*rank/np;ke=nk*(rank+1)/np
  first=1+na*na*rank/np;count=na*na*(rank+1)/np-first+1
  allocate(psi(ng,1,max(0,ke-ks+1)),w(count,nk),part(na*na,nk),total(na*na,nk))
  psi=1d0/sqrt(dble(ng));cell=0d0;cell(1,1)=12d0;cell(2,2)=3d0;cell(3,3)=3d0
  k=0d0;k(1,:)=[-.25d0,.25d0]*2d0*acos(-1d0)/12d0
  mode='valid'
  if(command_argument_count()>0)call get_command_argument(1,mode)
  p=ng;action_calls=0
  if(size(psi,3)>0)then
    if(mode=='nonfinite')psi(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
    if(mode=='infinite')psi(1,1,1)=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
  endif
  closure_tol=.05d0
  if(mode=='hermitian-coherent'.or.mode=='hermitian-q')then
    do j=1,size(psi,3)
      do ierr=0,ng-1
        psi(ierr+1,1,j)=(cos(.01d0)+cmplx(0d0,sin(.01d0),8)* &
          exp(cmplx(0d0,4d0*acos(-1d0)*mod(ierr,n(1))/n(1),8)))/sqrt(dble(ng))
      enddo
    enddo
  endif
  if(mode=='coherent-reject'.or.mode=='optical-coherent')then
    closure_tol=1d-12
    do j=1,size(psi,3)
      do ierr=0,ng-1
        psi(ierr+1,1,j)=(cos(.3d0)+cmplx(0d0,sin(.3d0),8)* &
          exp(cmplx(0d0,4d0*acos(-1d0)*mod(ierr,n(1))/n(1),8)))/sqrt(dble(ng))
      enddo
    enddo
  endif
  boxes=2
  if(mode=='empty-pairs'.or.mode=='optical-empty'.or.mode=='window-empty')then
    boxes=1;first=rank+1;count=merge(1,0,rank==0)
    deallocate(w);allocate(w(count,nk))
  endif
  if(mode=='cap')p=1
  if(mode=='q2-short'.or.mode=='q2-conflict')then
    if(mode=='q2-conflict')then
      call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
        first,count,p,1d-8,.01d0,closure_tol,action,w,status,velocity,density_gamma=.true.)
    else
      call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
        first,count,p,1d-8,.01d0,closure_tol,action,w,status,density_gamma=.true.)
    endif
    if(status/=11.or.action_calls/=0)error stop 'q2 Gamma preflight failed or performed H actions'
    if(rank==0)print *,'FROZEN_DENSITY_GAMMA_PREFLIGHT_PASS'
    call MPI_Finalize(ierr)
    stop
  endif
  if(index(mode,'window')==1)then
    allocate(reference_w(count,nk))
    nqproc=min(2,np);filter=1
    if(mode=='window-exp')filter=2
    if(mode=='window-missing')then
      call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
        first,count,p,1d-8,1d-4,closure_tol,action,w,status,closure='hermitian',prepare=prepare,nproc_q=nqproc, &
        window_kind=filter,window_time=10d0,window_cache=cache)
      if(status/=10.or.action_calls/=0)error stop 'uninitialized window cache accepted'
      if(rank==0)print *,'FROZEN_WINDOW_MISSING_CACHE_PASS'
      call MPI_Finalize(ierr)
      stop
    endif
    if(mode=='window-failed')then
      fail_window=.true.
    call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
      first,count,p,1d-8,1d-4,closure_tol,action,w,status,closure='hermitian',prepare=prepare,nproc_q=nqproc, &
      window_kind=filter,window_time=10d0,window_cache=cache,initialize_window=.true.)
      if(status/=23)error stop 'failed reference initialization accepted'
      fail_window=.false.;previous_calls=action_calls
    call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
      first,count,p,1d-8,1d-4,closure_tol,action,w,status,closure='hermitian',prepare=prepare,nproc_q=nqproc, &
      window_kind=filter,window_time=10d0,window_cache=cache,initialize_window=.true.)
      if(status/=10.or.action_calls/=previous_calls)error stop 'partial cache retried without reset'
      call frozen_window_reset(cache)
    endif
    call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
      first,count,p,1d-8,1d-4,closure_tol,action,w,status,closure='hermitian',prepare=prepare,nproc_q=nqproc, &
      window_kind=filter,window_time=10d0,window_cache=cache,initialize_window=.true.)
    if(status/=0)error stop 'window reference failed'
    reference_w=w
    call frozen_window_anchor(cache,2d0*w,1d-8,MPI_COMM_WORLD,status)
    if(status==0)error stop 'mismatched equilibrium W accepted'
    call frozen_window_anchor(cache,w,1d-8,MPI_COMM_WORLD,status)
    if(status/=0)error stop 'matching equilibrium W rejected'
    if(mode=='window-key')then
      previous_calls=action_calls
      call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
        first,count,p,1d-8,1d-4,closure_tol,action,w,status,closure='hermitian',prepare=prepare,nproc_q=nqproc, &
        window_kind=filter,window_time=11d0,window_cache=cache)
      if(status/=10.or.action_calls/=previous_calls)error stop 'incompatible window cache accepted'
    endif
    if(mode=='window-one-reset')then
      if(rank==np-1)call frozen_window_reset(cache)
      previous_calls=action_calls
      call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
        first,count,p,1d-8,1d-4,closure_tol,action,w,status,closure='hermitian',prepare=prepare,nproc_q=nqproc, &
        window_kind=filter,window_time=10d0,window_cache=cache)
      if(status/=10.or.action_calls/=previous_calls)error stop 'rank-local missing cache not collective'
      if(rank==0)print *,'FROZEN_WINDOW_COLLECTIVE_CACHE_REJECT_PASS'
      call MPI_Finalize(ierr)
      stop
    endif

    call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
      first,count,p,1d-8,1d-4,closure_tol,action,w,status,closure='hermitian',prepare=prepare,nproc_q=nqproc, &
      window_kind=filter,window_time=10d0,window_cache=cache)
    if(status/=0.or.maxval(abs(w-reference_w))>1d-10)error stop 'window equilibrium subtraction failed'
    if(rank==0)print *,'FROZEN_WINDOW_REFERENCE_ANCHOR_ZERO_PASS'
    call MPI_Finalize(ierr)
    stop
  else if(index(mode,'hermitian')==1)then
    eta=1d-4
    if(mode=='hermitian-invalid')then
      call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
        first,count,p,1d-8,eta,closure_tol,action,w,status,closure='invalid')
      if(status/=10.or.action_calls/=0)error stop 'invalid closure did not reject before H actions'
      if(rank==0)print *,'FROZEN_STATIC_INVALID_PASS'
      call MPI_Finalize(ierr)
      stop
    endif
    if(mode=='hermitian-coherent'.or.mode=='hermitian-q')then
      call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
        first,count,p,1d-8,eta,closure_tol,action,w,status)
      if(status/=0.and.status/=32)error stop 'strict coherent reference failed before static closure'
      allocate(reference_w(count,nk))
      call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
        first,count,p,1d-8,eta,closure_tol,action,reference_w,status,closure='hermitian')
      if(status/=0)error stop 'projected coherent reference rejected'
      nqproc=1
      if(mode=='hermitian-q')nqproc=min(2,np)
      call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
        first,count,p,1d-8,eta,closure_tol,action,w,status,closure='hermitian',prepare=prepare,nproc_q=nqproc)
    else
      if(mode=='hermitian-eta')eta=.1d0
      call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
        first,count,p,1d-8,eta,closure_tol,action,w,status,closure='hermitian')
    endif
    if(mode=='hermitian-eta')then
      if(status/=34)error stop 'projected mode skipped eta convergence guard'
      if(rank==0)print *,'FROZEN_STATIC_ETA_REJECT_PASS'
      call MPI_Finalize(ierr)
      stop
    endif
  else if(index(mode,'pairs')==1)then
    if(mode=='pairs-coherent')then
      do j=1,size(psi,3)
        do ierr=0,ng-1
          psi(ierr+1,1,j)=(cos(1d-4)+cmplx(0d0,sin(1d-4),8)* &
            exp(cmplx(0d0,4d0*acos(-1d0)*mod(ierr,n(1))/n(1),8)))/sqrt(dble(ng))
        enddo
      enddo
      allocate(reference_w(count,nk))
      call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
        first,count,p,1d-8,.01d0,closure_tol,action,reference_w,status)
      if(status/=0)error stop 'coherent reference rejected'
    endif
    nqproc=1
    if(mode=='pairs-invalid')nqproc=np+1
    if(mode=='pairs-q'.or.mode=='pairs-bad')nqproc=min(2,np)
    call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
      first,count,p,1d-8,.01d0,closure_tol,action,w,status,prepare=prepare,nproc_q=nqproc)
  else if(index(mode,'optical')==1)then
    eta=1d-4
    if(mode=='optical-eta-reject')eta=.1d0
    call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
      first,count,p,1d-8,eta,closure_tol,action,w,status,velocity)
  else
    call frozen_full_w(n,[2,1,1],[2,1,1],boxes,.2d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
      first,count,p,1d-8,.01d0,closure_tol,action,w,status)
  endif
  if(mode=='pairs-prepare'.or.mode=='pairs-invalid')then
    if(mode=='pairs-prepare'.and.status/=25)error stop 'frame prepare failure not collective'
    if(mode=='pairs-invalid'.and.status/=10)error stop 'invalid q layout not rejected'
    if(rank==0)print *,'FROZEN_PAIR_PREFLIGHT_FAILURE_PASS'
    call MPI_Finalize(ierr)
    stop
  endif
  if(mode=='pairs-bad')then
    if(status/=23)error stop 'distributed action failure not collective'
    if(rank==0)print *,'FROZEN_PAIR_FAILURE_PASS'
    call MPI_Finalize(ierr)
    stop
  endif
  if(mode=='bad-action'.or.mode=='cap'.or.mode=='coherent-reject'.or.mode=='nonfinite'.or.mode=='infinite'.or. &
      mode=='optical-failure'.or.mode=='optical-coherent'.or.mode=='optical-nan'.or.mode=='optical-eta-reject')then
    if(status==0)error stop 'failed full response accepted'
    if((mode=='nonfinite'.or.mode=='infinite').and.status/=10)error stop 'nonfinite response rejected after input guard'
    if(mode=='coherent-reject'.and.status/=32)error stop 'coherent response rejected before static closure check'
    if(index(mode,'optical')==1.and.mode/='optical-eta-reject'.and.status/=50) &
      error stop 'optical input rejected at wrong stage'
    if(mode=='optical-eta-reject'.and.status/=34)error stop 'eta limit rejected at wrong stage'
    if(rank==0)print *,'FROZEN_FULL_FAILURE_REJECTED'
  else
    if(status/=0)error stop 'full W rejected'
    part=0d0;part(first:first+count-1,:)=w
    call MPI_Allreduce(part,total,size(part),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    do j=1,nk
      if(maxval(abs(reshape(total(:,j),[na,na])-transpose(reshape(total(:,j),[na,na]))))>1d-8) &
        error stop 'nonreciprocal W accepted'
    enddo
    if(maxval(abs(total))<1d-3)error stop 'empty W'
    if(mode=='pairs-coherent'.or.mode=='hermitian-coherent'.or.mode=='hermitian-q')then
      if(maxval(abs(w-reference_w))>1d-10)error stop 'coherent MPI reference mismatch'
    else if(mode/='empty-pairs'.and.mode/='optical-empty')then
      call analytic_reference(expected)
      if(maxval(abs(total-expected))>1d-10)error stop 'analytic spin/transition/RPA normalization mismatch'
    endif
    if(rank==0)then
      print '(a,2es22.12)','FROZEN_FULL_W_CHECKSUM ',sum(total),sum(total**2)
      print *,'FROZEN_FULL_ALL_Q_LOCAL_FIELDS_PASS'
    endif
  endif
  call MPI_Finalize(ierr)
contains
  subroutine prepare(source,target,worker,source_owner,target_owner,code)
    implicit none
    integer,intent(in) :: source,target,worker,source_owner,target_owner
    integer,intent(out) :: code
    code=0
    if(mode=='pairs-prepare'.and.rank==np-1)code=1
    if(min(source,target)<1.or.max(source,target)>nk)code=1
    if(min(worker,source_owner,target_owner)<0.or.max(worker,source_owner,target_owner)>=np)code=1
  end subroutine
  subroutine analytic_reference(result)
    implicit none
    real(8),intent(out) :: result(na*na,nk)
    real(8) :: rec(3,3),vol,q(3),qp(3,3),head,qhead,einv,qc,fk(2),qt,source_energy,target_energy,difference,macro(3,3)
    complex(8) :: chi(3,3),value(na,na),values(na*na,nk),z
    complex(8) :: optical(5,5),point(na,na),origin(na,na)
    integer :: optical_g(3,2)
    integer :: iq,source,target,g,mg,lo,hi,fold,code
    fk=[-.25d0,.25d0];z=cmplx(0d0,.005d0,8);einv=0d0
    if(index(mode,'optical')==1.or.index(mode,'hermitian')==1)z=0d0
    call screen_geometry(cell,rec,vol,code)
    if(code/=0)error stop 'reference geometry'
    do iq=1,nk
      qt=merge(0d0,-.5d0,iq==1);q=[qt,0d0,0d0]
      lo=merge(-1,0,iq==1);hi=1;mg=hi-lo+1;chi=0d0;qp=0d0
      do g=lo,hi
        qp(:,g-lo+1)=[qt+g,0d0,0d0]
        do source=1,nk
          target=source
          if(iq==2)target=3-source
          fold=nint(fk(source)+qt-fk(target))
          source_energy=.2d0-.2d0*cos(2d0*acos(-1d0)*fk(source)/n(1))
          target_energy=.2d0-.2d0*cos(2d0*acos(-1d0)*(g+fold+fk(target))/n(1))
          difference=target_energy-source_energy
          if(abs(z-difference)>1d-14)chi(g-lo+1,g-lo+1)=chi(g-lo+1,g-lo+1)+2d0/(z-difference)
          source_energy=.2d0-.2d0*cos(2d0*acos(-1d0)*(-g-fold+fk(source))/n(1))
          target_energy=.2d0-.2d0*cos(2d0*acos(-1d0)*fk(target)/n(1))
          difference=target_energy-source_energy
          if(abs(z-difference)>1d-14)chi(g-lo+1,g-lo+1)=chi(g-lo+1,g-lo+1)-2d0/(z-difference)
        enddo
      enddo
      chi=.5d0*(chi+conjg(transpose(chi)))
      call screen_assemble(n,2,qp(:,:mg),rec,vol,nk,-.5d0*chi(:mg,:mg),.true.,value,head,qhead,code)
      if(code/=0)error stop 'analytic RPA reference'
      if(iq==2)einv=head
      values(:,iq)=reshape(value,[na*na])
      if(iq==1.and.index(mode,'optical')==1)then
        optical_g=0;optical_g(1,:)=[-1,1];optical=0d0
        optical(4,4)=-.5d0*chi(1,1);optical(5,5)=-.5d0*chi(3,3)
        call longwave_cell(n,2,optical_g,rec,vol,nk,optical,macro,point,origin,code)
        if(code/=0)error stop 'analytic optical reference failed'
        values(:,iq)=reshape(origin,[na*na])
      endif
    enddo
    call screen_transform(values,[2,1,1])
    qc=(6d0*acos(-1d0)**2/(vol*nk))**(1d0/3d0)
    result=real(values,8)
    if(index(mode,'optical')/=1)result=result+12d0*acos(-1d0)*einv/(qc*qc*vol*nk)
  end subroutine
  subroutine velocity(index,x,y,code)
    implicit none
    integer,intent(in) :: index
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:,:)
    integer,intent(out) :: code
    code=1
    if(mode(:5)/='pairs'.and.mode(:9)/='hermitian'.and.mode(:6)/='window')then
      if(index<ks.or.index>ke)return
    endif
    if(mode=='pairs-bad'.and.index==nk)return
    if(mode=='optical-failure'.and.index==nk)return
    ! The plane-wave occupied state has no interband velocity elements.
    ! Consequently its exact optical epsilon is one, while finite-q epsilon
    ! is not: this catches retention/double counting of the old Gamma head.
    y=0d0;code=0
    if(mode=='optical-nan'.and.index==nk)y(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  end subroutine
  subroutine action(index,x,y,code)
    implicit none
    integer,intent(in) :: index
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    integer :: a,left,right
    real(8) :: angle
    action_calls=action_calls+1
    code=1
    if(fail_window)return
    if(mode(:5)/='pairs'.and.mode(:9)/='hermitian'.and.mode(:6)/='window')then
      if(index<ks.or.index>ke)return
    endif
    if(mode=='pairs-bad'.and.index==nk)return
    if(mode=='bad-action'.and.index==nk)return
    angle=2d0*acos(-1d0)*merge(-.25d0,.25d0,index==1)/n(1)
    do a=0,ng-1
      left=a-mod(a,n(1))+modulo(mod(a,n(1))-1,n(1))+1
      right=a-mod(a,n(1))+mod(mod(a,n(1))+1,n(1))+1
      y(a+1)=.2d0*x(a+1)-.1d0*(exp(cmplx(0d0,angle,8))*x(right)+exp(cmplx(0d0,-angle,8))*x(left))
    enddo
    code=0
  end subroutine
end program
