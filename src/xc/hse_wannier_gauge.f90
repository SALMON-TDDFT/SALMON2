! Occupied/retained-subspace MV localization and polar temporal transport.
! Port of the tested TDCDFT reference, using explicit state and grid-weighted overlaps.
module hse_wannier_gauge
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: gauge_transport,gauge_functional,gauge_minimize,gauge_seed,gauge_minimize_gamma
contains
  subroutine gauge_seed(psi,position,k,u,status)
    implicit none
    complex(8),intent(in) :: psi(:,:,:)
    real(8),intent(in) :: position(:,:),k(:,:)
    complex(8),intent(out) :: u(:,:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: columns(:,:),tau(:),work(:),overlap(:,:),left(:,:),right(:,:)
    real(8),allocatable :: rwork(:),singular(:)
    integer,allocatable :: pivots(:)
    integer :: n,ng,ik,j,point
    n=size(psi,2);ng=size(psi,1);status=1
    if(n>ng)return
    allocate(columns(n,ng),pivots(ng),tau(n),work(max(8*n,34*ng+32)))
    allocate(overlap(n,n),left(n,n),right(n,n),rwork(max(2*ng,5*n)),singular(n))
    ! Pivoted columns of the retained projector select spatial delta trials.
    columns=conjg(transpose(psi(:,:,1)));pivots=0
    call zgeqp3(n,ng,columns,n,pivots,tau,work,size(work),rwork,status)
    if(status/=0)return
    do ik=1,size(psi,3)
      do j=1,n
        point=pivots(j)
        overlap(:,j)=conjg(psi(point,:,ik))*exp(cmplx(0d0,-sum(k(:,ik)*position(:,point)),8))
      enddo
      call zgesvd('A','A',n,n,overlap,n,singular,left,n,right,n,work,size(work),rwork,status)
      if(status/=0)return
      if(minval(singular)<1d-10*maxval(singular))then
        status=1;return
      endif
      u(:,:,ik)=matmul(left,right)
    enddo
    status=0
  end subroutine
  subroutine gauge_transport(current,previous,dv,u,min_singular,status)
    implicit none
    complex(8),intent(in) :: current(:,:,:),previous(:,:,:)
    real(8),intent(in) :: dv
    complex(8),intent(out) :: u(:,:,:)
    real(8),intent(out) :: min_singular
    integer,intent(out) :: status
    complex(8),allocatable :: overlap(:,:),left(:,:),right(:,:),work(:)
    real(8),allocatable :: singular(:),rwork(:)
    integer :: n,ik
    status=1;min_singular=0d0
    if(any(shape(current)/=shape(previous)).or.dv<=0d0)return
    n=size(current,2)
    if(any(shape(u)/=[n,n,size(current,3)]))return
    allocate(overlap(n,n),left(n,n),right(n,n),work(8*n),singular(n),rwork(5*n))
    min_singular=huge(1d0)
    do ik=1,size(current,3)
      overlap=matmul(conjg(transpose(current(:,:,ik))),previous(:,:,ik))*dv
      call zgesvd('A','A',n,n,overlap,n,singular,left,n,right,n,work,size(work),rwork,status)
      if(status/=0)return
      min_singular=min(min_singular,minval(singular))
      if(.not.all(ieee_is_finite(singular)).or.min_singular<1d-8)then
        status=1;return
      endif
      u(:,:,ik)=matmul(left,right)
    enddo
    status=0
  end subroutine

  subroutine gauge_functional(u,raw,neighbors,b,weights,spread,variable,gradient,status,phase_reference,phase_out)
    implicit none
    complex(8),intent(in) :: u(:,:,:),raw(:,:,:,:)
    integer,intent(in) :: neighbors(:,:)
    real(8),intent(in) :: b(:,:),weights(:)
    real(8),intent(out) :: spread,variable
    real(8),intent(in),optional :: phase_reference(:,:,:)
    real(8),intent(out),optional :: phase_out(:,:,:)
    complex(8),intent(out) :: gradient(:,:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: m(:,:,:,:),back(:,:,:),q(:)
    real(8),allocatable :: theta(:,:,:),center(:,:)
    real(8) :: residual,weight
    integer :: n,nk,nb,k,l,j,i,next
    n=size(u,1);nk=size(u,3);nb=size(weights)
    allocate(m(n,n,nb,nk),back(n,n,nk),q(n),theta(n,nb,nk),center(3,n))
    center=0d0;back=0d0;spread=0d0;variable=0d0;gradient=0d0;status=1
    do k=1,nk
      do l=1,nb
        next=neighbors(l,k)
        m(:,:,l,k)=matmul(conjg(transpose(u(:,:,k))),matmul(raw(:,:,l,k),u(:,:,next)))
        do i=1,n
          if(abs(m(i,i,l,k))<1d-12)return
          theta(i,l,k)=atan2(aimag(m(i,i,l,k)),real(m(i,i,l,k),8))
          ! Follow the accepted link branch through the principal atan2 cut.
          ! Reverse links inherit opposite branches from their accepted phases.
          if(present(phase_reference))theta(i,l,k)=theta(i,l,k)+2*acos(-1d0)* &
            anint((phase_reference(i,l,k)-theta(i,l,k))/(2*acos(-1d0)))
          center(:,i)=center(:,i)-weights(l)*b(:,l)*theta(i,l,k)/nk
        enddo
      enddo
    enddo
    do k=1,nk
      do l=1,nb
        next=neighbors(l,k);weight=weights(l)/nk
        do i=1,n
          residual=theta(i,l,k)+sum(b(:,l)*center(:,i))
          spread=spread+weight*(1d0-abs(m(i,i,l,k))**2+residual**2)
          variable=variable+weight*residual**2
          q(i)=weight*(-2d0*conjg(m(i,i,l,k))-cmplx(0d0,2d0*residual,8)/m(i,i,l,k))
        enddo
        do j=1,n
          do i=1,n
            if(i/=j)variable=variable+weight*abs(m(i,j,l,k))**2
            back(i,j,k)=back(i,j,k)-m(i,j,l,k)*q(j)
            back(i,j,next)=back(i,j,next)+q(i)*m(i,j,l,k)
          enddo
        enddo
      enddo
    enddo
    do k=1,nk
      gradient(:,:,k)=0.5d0*(back(:,:,k)-conjg(transpose(back(:,:,k))))
    enddo
    if(present(phase_out))phase_out=theta
    status=0
  end subroutine

  subroutine gauge_minimize(u,raw,neighbors,b,weights,maxiter,tolerance,spread,gradnorm,iterations,status)
    implicit none
    complex(8),intent(inout) :: u(:,:,:)
    complex(8),intent(in) :: raw(:,:,:,:)
    integer,intent(in) :: neighbors(:,:),maxiter
    real(8),intent(in) :: b(:,:),weights(:),tolerance
    real(8),intent(out) :: spread,gradnorm
    integer,intent(out) :: iterations,status
    complex(8),allocatable :: direction(:,:,:),candidate(:,:,:),dtrial(:,:,:),vectors(:,:),work(:),rotation(:,:)
    real(8),allocatable :: eigen(:),rwork(:),phases(:,:,:),trialphases(:,:,:)
    real(8) :: variable,trial,trialspread,step
    integer :: n,nk,k,j,backtrack,stat
    n=size(u,1);nk=size(u,3)
    allocate(direction(n,n,nk),candidate(n,n,nk),dtrial(n,n,nk),vectors(n,n),rotation(n,n))
    allocate(eigen(n),work(4*n),rwork(3*n))
    allocate(phases(n,size(weights),nk),trialphases(n,size(weights),nk))
    step=1d0;status=1;gradnorm=huge(1d0);spread=huge(1d0)
    do iterations=0,maxiter
      if(iterations==0)then
        call gauge_functional(u,raw,neighbors,b,weights,spread,variable,direction,stat,phase_out=phases)
      else
        call gauge_functional(u,raw,neighbors,b,weights,spread,variable,direction,stat, &
                              phase_reference=phases,phase_out=trialphases)
      endif
      if(stat/=0)return
      gradnorm=sqrt(sum(abs(direction)**2))
      if(gradnorm<tolerance)then
        status=0;return
      endif
      if(iterations==maxiter)return
      do backtrack=1,35
        do k=1,nk
          vectors=cmplx(0d0,-1d0,8)*direction(:,:,k)
          call zheev('V','U',n,vectors,n,eigen,work,size(work),rwork,stat)
          if(stat/=0)return
          rotation=vectors
          do j=1,n
            rotation(:,j)=rotation(:,j)*exp(cmplx(0d0,step*eigen(j),8))
          enddo
          candidate(:,:,k)=matmul(u(:,:,k),matmul(rotation,conjg(transpose(vectors))))
        enddo
        call gauge_functional(candidate,raw,neighbors,b,weights,trialspread,trial,dtrial,stat, &
                              phase_reference=phases,phase_out=trialphases)
        if(stat==0)then
          if(trial<=variable-1d-4*step*gradnorm**2)exit
        endif
        step=0.5d0*step
      enddo
      if(backtrack>35)return
      u=candidate;phases=trialphases;step=min(1.5d0*step,1000d0)
    enddo
  end subroutine
  subroutine gauge_minimize_gamma(u,raw,b,weights,maxiter,tolerance,spread,gradnorm,iterations,status)
    ! At Gamma, the +/- link center residual cancels exactly. Minimizing MV
    ! spread is simultaneous maximization of the squared link diagonals.
    ! Each complex two-column SU(2) rotation maximizes n^T G n on |n|=1,
    ! G_ab=sum_link weight*Re(conjg(v_a)*v_b), A=a0 I+v.sigma.
    ! Unlike global steepest descent, distant centers do not limit every step.
    complex(8),intent(inout) :: u(:,:,:)
    complex(8),intent(in) :: raw(:,:,:,:)
    real(8),intent(in) :: b(:,:),weights(:),tolerance
    integer,intent(in) :: maxiter
    real(8),intent(out) :: spread,gradnorm
    integer,intent(out) :: iterations,status
    complex(8),allocatable :: links(:,:,:),tmp(:),direction(:,:,:)
    complex(8) :: v(3),sine
    real(8) :: metric(3,3),original(3,3),eval(3),work(32),axis(3),cosine,variable,scale
    integer :: n,nb,l,i,j,a,c,istat,neighbors(size(weights),1)
    n=size(u,1);nb=size(weights);status=1;iterations=0
    spread=huge(1d0);gradnorm=huge(1d0)
    if(size(u,3)/=1.or.size(raw,4)/=1.or.nb/=6)return
    if(maxval(abs(b(:,1:3)+b(:,4:6)))>1d-12)return
    if(maxval(abs(weights(1:3)-weights(4:6)))>1d-12)return
    allocate(links(n,n,nb),tmp(n),direction(n,n,1));neighbors=1
    do l=1,nb
      links(:,:,l)=matmul(conjg(transpose(u(:,:,1))),matmul(raw(:,:,l,1),u(:,:,1)))
    enddo
    do iterations=0,maxiter
      call gauge_functional(u,raw,neighbors,b,weights,spread,variable,direction,istat)
      if(istat/=0)return
      gradnorm=sqrt(sum(abs(direction)**2))
      if(gradnorm<tolerance.and.iterations>0)then
        status=0;return
      endif
      if(iterations==maxiter)return
      do j=2,n;do i=1,j-1
        metric=0d0
        do l=1,nb
          v(1)=.5d0*(links(i,j,l)+links(j,i,l))
          v(2)=cmplx(0d0,.5d0,8)*(links(i,j,l)-links(j,i,l))
          v(3)=.5d0*(links(i,i,l)-links(j,j,l))
          do c=1,3;do a=1,3
            metric(a,c)=metric(a,c)+weights(l)*real(conjg(v(a))*v(c),8)
          enddo;enddo
        enddo
        original=metric
        call dsyev('V','U',3,metric,3,eval,work,size(work),istat)
        if(istat/=0)return
        scale=maxval(abs(original))
        if(eval(3)<=original(3,3)+32*epsilon(1d0)*scale.and. &
           maxval(abs(original(1:2,3)))<tolerance/(16*n))cycle
        axis=metric(:,3)
        if(axis(3)<0d0)axis=-axis
        cosine=sqrt(.5d0*(1d0+axis(3)))
        sine=cmplx(axis(1),axis(2),8)/(2*cosine)
        if(abs(sine)<1d-15)cycle
        do l=1,nb
          tmp=links(:,i,l)
          links(:,i,l)=cosine*tmp+sine*links(:,j,l)
          links(:,j,l)=-conjg(sine)*tmp+cosine*links(:,j,l)
          tmp=links(i,:,l)
          links(i,:,l)=cosine*tmp+conjg(sine)*links(j,:,l)
          links(j,:,l)=-sine*tmp+cosine*links(j,:,l)
        enddo
        tmp=u(:,i,1)
        u(:,i,1)=cosine*tmp+sine*u(:,j,1)
        u(:,j,1)=-conjg(sine)*tmp+cosine*u(:,j,1)
      enddo;enddo
    enddo
  end subroutine
end module
