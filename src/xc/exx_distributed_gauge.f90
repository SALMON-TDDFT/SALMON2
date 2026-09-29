#include "config.h"
! Gamma MLWF matrices in unique block-cyclic tiles. Mesh functions stay local.
module exx_distributed_gauge
  use communication, only: comm_get_groupinfo,comm_get_max,comm_summation,comm_bcast
  use exx_orbitals, only: orbital_layout,orbital_check
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: s_exx_gauge,gauge_tiles_clear,gauge_tiles_refresh,gauge_tiles_rotate
  type s_exx_gauge
    integer :: n=0,comm=0
    integer,allocatable :: rows(:),cols(:)
    complex(8),allocatable :: matrix(:,:)
  end type
contains
  subroutine gauge_tiles_clear(gauge)
    implicit none
    type(s_exx_gauge),intent(inout) :: gauge
    type(s_exx_gauge) :: empty
    gauge=empty
  end subroutine

  subroutine gauge_tiles_refresh(gauge,psi,previous,dv,phase,b,weights,comm_r,comm_o,comm, &
      maxiter,tolerance,seed,needed,retain,last_status,retained,minimum,spread,gradient,iterations,loc_status,status)
#ifdef USE_SCALAPACK
    use exx_blas_threads, only: exx_blas_thread_control
!$  use omp_lib, only: omp_get_max_threads,omp_in_parallel
#endif
    implicit none
    type(s_exx_gauge),intent(inout) :: gauge
    complex(8),intent(in) :: psi(:,:),phase(:,:)
    complex(8),allocatable,intent(in) :: previous(:,:,:)
    real(8),intent(in) :: dv,b(3,6),weights(6),tolerance
    integer,intent(in) :: comm_r,comm_o,comm,maxiter
    logical,intent(in) :: seed,retain
    logical,intent(inout) :: needed
    integer,intent(inout) :: last_status
    logical,intent(out) :: retained
    real(8),intent(out) :: minimum,spread,gradient
    integer,intent(out) :: iterations,loc_status,status
#ifdef USE_SCALAPACK
    integer :: workers,thread_state(3)
    workers=1;thread_state=0
!$  workers=omp_get_max_threads()
!$  if(omp_in_parallel())workers=1
    if(workers>1)call exx_blas_thread_control(workers,thread_state)
#endif
    ! Keep one cleanup path for normal completion and every early return below.
    call refresh_core()
#ifdef USE_SCALAPACK
    call exx_blas_thread_control(0,thread_state)
#endif
  contains
    subroutine refresh_core()
      implicit none
#ifdef USE_SCALAPACK
    integer :: context,desc(9),nr,nc,n,first,bad,axis
    integer,allocatable :: counts(:)
    complex(8),allocatable :: links(:,:,:),a(:,:),tmp(:,:),transported(:,:)
    logical :: accepted,initialized
    status=1;retained=.false.;minimum=0d0
    spread=-1d0;gradient=-1d0;iterations=0;loc_status=2
    call orbital_layout(size(psi,2),comm_r,comm_o,counts,first,bad)
    if(bad/=0.or.sum(counts)<1)return
    n=sum(counts);initialized=allocated(gauge%matrix)
    bad=merge(1,0,initialized)
    call comm_get_max(bad,comm)
    bad=abs(merge(1,0,initialized)-bad)
    call comm_get_max(bad,comm)
    if(bad/=0)return
    if(initialized)then
      bad=0
      if(gauge%n/=n.or.gauge%comm/=comm.or..not.allocated(previous))bad=1
      if(allocated(previous))then
        if(any(shape(previous)/=[size(psi,1),size(psi,2),1]))bad=1
        if(size(previous,3)==1)then
          if(.not.finite_matrix(previous(:,:,1)))bad=1
        endif
      endif
      call comm_get_max(bad,comm)
      if(bad/=0)return
    endif
    call setup(gauge,n,comm,context,desc,bad)
    if(bad/=0)then
      call blacs_gridexit(context)
      return
    endif
    nr=size(gauge%rows);nc=size(gauge%cols)
    if(initialized)then
      allocate(a(nr,nc))
      call overlap(gauge,psi,previous(:,:,1),dv,comm_o,counts,first,a)
      call polar(gauge,a,desc,minimum,bad)
      deallocate(a)
      if(bad/=0)then
        call identity(gauge)
        last_status=1;needed=.true.
      endif
    else
      call identity(gauge)
    endif
    accepted=retain.and.seed.and.last_status==0
    if(accepted)transported=gauge%matrix
    if(maxiter>0)then
      allocate(links(nr,nc,6),tmp(nr,nc),a(nr,nc))
      do axis=1,3
        call overlap(gauge,psi,psi,dv,comm_o,counts,first,links(:,:,axis),phase(:,axis))
        call pztranc(n,n,(1d0,0d0),links(:,:,axis),1,1,desc,(0d0,0d0),links(:,:,axis+3),1,1,desc)
      enddo
      if(seed.and.needed)then
        call position_seed(gauge,links,desc,bad)
        if(bad/=0)then
          call blacs_gridexit(context)
          return
        endif
        needed=.false.
      endif
      do axis=1,6
        call pzgemm('N','N',n,n,n,(1d0,0d0),links(:,:,axis),1,1,desc,gauge%matrix,1,1,desc, &
          (0d0,0d0),tmp,1,1,desc)
        call pzgemm('C','N',n,n,n,(1d0,0d0),gauge%matrix,1,1,desc,tmp,1,1,desc,(0d0,0d0),a,1,1,desc)
        links(:,:,axis)=a
      enddo
      deallocate(tmp,a)
      call minimize(gauge,links,desc,b,weights,maxiter,tolerance,spread,gradient,iterations,loc_status)
      if(loc_status/=0.and.accepted)then
        gauge%matrix=transported;retained=.true.
      else
        last_status=loc_status
      endif
    endif
    bad=0
    if(.not.finite_matrix(gauge%matrix))bad=1
    call comm_get_max(bad,comm)
    call blacs_gridexit(context)
    if(bad==0)status=0
#else
    status=1;retained=.false.;minimum=0d0
    spread=-1d0;gradient=-1d0;iterations=0;loc_status=2
#endif
    end subroutine refresh_core
  end subroutine gauge_tiles_refresh

#ifdef USE_SCALAPACK
  subroutine setup(gauge,n,comm,context,desc,status)
    implicit none
    type(s_exx_gauge),intent(inout) :: gauge
    integer,intent(in) :: n,comm
    integer,intent(out) :: context,desc(9),status
    integer :: rank,np,nprow,npcol,myrow,mycol,block,nr,nc,i,k
    integer,external :: sys2blacs_handle,numroc
    call comm_get_groupinfo(comm,rank,np)
    nprow=int(sqrt(dble(np)))
    do while(mod(np,nprow)/=0)
      nprow=nprow-1
    enddo
    npcol=np/nprow
    context=sys2blacs_handle(comm)
    call blacs_gridinit(context,'R',nprow,npcol)
    call blacs_gridinfo(context,nprow,npcol,myrow,mycol)
    block=min(32,max(1,n/max(nprow,npcol)))
    nr=max(1,numroc(n,block,myrow,0,nprow));nc=max(1,numroc(n,block,mycol,0,npcol))
    call descinit(desc,n,n,block,block,0,0,context,nr,status)
    status=abs(status);call comm_get_max(status,comm)
    if(status/=0)return
    if(allocated(gauge%matrix))then
      if(any(shape(gauge%matrix)/=[nr,nc]))status=1
      call comm_get_max(status,comm)
      return
    endif
    allocate(gauge%matrix(nr,nc),gauge%rows(nr),gauge%cols(nc))
    gauge%n=n;gauge%comm=comm;gauge%matrix=0d0;gauge%rows=0;gauge%cols=0
    do i=1,n
      if(mod((i-1)/block,nprow)==myrow)then
        k=((i-1)/(block*nprow))*block+mod(i-1,block)+1;gauge%rows(k)=i
      endif
      if(mod((i-1)/block,npcol)==mycol)then
        k=((i-1)/(block*npcol))*block+mod(i-1,block)+1;gauge%cols(k)=i
      endif
    enddo
  end subroutine

  subroutine identity(gauge)
    implicit none
    type(s_exx_gauge),intent(inout) :: gauge
    integer :: i,j
    gauge%matrix=0d0
    do j=1,size(gauge%cols)
      do i=1,size(gauge%rows)
        if(gauge%rows(i)>0.and.gauge%rows(i)==gauge%cols(j))gauge%matrix(i,j)=1d0
      enddo
    enddo
  end subroutine

  subroutine overlap(gauge,left,right,dv,comm_o,counts,first,a,phase)
    implicit none
    type(s_exx_gauge),intent(in) :: gauge
    complex(8),intent(in) :: left(:,:),right(:,:)
    complex(8),intent(in),optional :: phase(:)
    real(8),intent(in) :: dv
    integer,intent(in) :: comm_o,counts(0:),first
    complex(8),intent(out) :: a(:,:)
    complex(8),allocatable :: column(:),partial(:),total(:)
    integer :: rank,np,owner,k,col,i,j
    call comm_get_groupinfo(comm_o,rank,np)
    allocate(column(size(right,1)),partial(gauge%n),total(gauge%n))
    a=0d0;k=0
    do owner=0,np-1
      do col=1,counts(owner)
        k=k+1
        if(rank==owner)column=right(:,col)
        call comm_bcast(column,comm_o,owner)
        if(present(phase))column=column*phase
        partial=0d0
        do i=1,size(left,2)
          partial(first+i-1)=dv*dot_product(left(:,i),column)
        enddo
        call comm_summation(partial,total,gauge%n,gauge%comm)
        do j=1,size(gauge%cols)
          if(gauge%cols(j)/=k)cycle
          do i=1,size(gauge%rows)
            if(gauge%rows(i)>0)a(i,j)=total(gauge%rows(i))
          enddo
        enddo
      enddo
    enddo
  end subroutine

  subroutine polar(gauge,a,desc,minimum,status)
    implicit none
    type(s_exx_gauge),intent(inout) :: gauge
    complex(8),intent(inout) :: a(:,:)
    integer,intent(in) :: desc(9)
    real(8),intent(out) :: minimum
    integer,intent(out) :: status
    complex(8),allocatable :: left(:,:),right(:,:),work(:)
    real(8),allocatable :: singular(:),rwork(:)
    complex(8) :: query(1)
    integer :: n,i,heterogeneous_status
    logical :: heterogeneous
    n=gauge%n;minimum=0d0;status=0
    if(.not.finite_matrix(a))status=1
    call comm_get_max(status,gauge%comm)
    if(status/=0)return
    allocate(left(size(a,1),size(a,2)),right(size(a,1),size(a,2)),singular(n),rwork(1+4*n))
    left=0d0;right=0d0
    gauge%matrix=a ! Keep the original overlap for a heterogeneous-SVD retry.
    call pzgesvd('V','V',n,n,a,1,1,desc,singular,left,1,1,desc,right,1,1,desc,query,-1,rwork,status)
    status=abs(status);call comm_get_max(status,gauge%comm)
    if(status/=0)return
    allocate(work(max(1,int(real(query(1))))))
    call pzgesvd('V','V',n,n,a,1,1,desc,singular,left,1,1,desc,right,1,1,desc,work,size(work),rwork,status)
    ! Preserve INFO's sign: a negative argument error must never be mistaken
    ! for the positive N+1 notification, nor hidden by another rank's N+1.
    heterogeneous_status=merge(1,0,status==n+1)
    if(status==n+1)status=0
    status=abs(status);call comm_get_max(status,gauge%comm)
    if(status/=0)return
    call comm_get_max(heterogeneous_status,gauge%comm)
    heterogeneous=heterogeneous_status==1
    minimum=minval(singular)
    do i=1,n
      if(.not.ieee_is_finite(singular(i)))status=1
    enddo
    if(minimum<1d-8)status=1
    call comm_get_max(status,gauge%comm)
    if(status/=0)return
    if(heterogeneous)then
      deallocate(left,right,work,rwork)
      call polar_iteration(gauge,desc,status)
    else
      call pzgemm('N','N',n,n,n,(1d0,0d0),left,1,1,desc,right,1,1,desc,(0d0,0d0),gauge%matrix,1,1,desc)
    endif
  end subroutine

  subroutine polar_iteration(gauge,desc,status)
    ! Newton-Schulz polar iteration. Frobenius scaling puts all singular values
    ! in (0,1]; the polar factor is unchanged. Used only after INFO=N+1, with
    ! the original SVD rank-loss threshold already checked on every rank.
    implicit none
    type(s_exx_gauge),intent(inout) :: gauge
    integer,intent(in) :: desc(9)
    integer,intent(out) :: status
    complex(8),allocatable :: gram(:,:),next(:,:)
    real(8) :: local_norm,total_norm,error,total_error
    integer :: n,i,j,iteration
    n=gauge%n;status=1
    local_norm=sum(abs(gauge%matrix)**2)
    call comm_summation(local_norm,total_norm,gauge%comm)
    if(.not.ieee_is_finite(total_norm).or.total_norm<=0d0)return
    gauge%matrix=gauge%matrix/sqrt(total_norm)
    allocate(gram(size(gauge%matrix,1),size(gauge%matrix,2)),next(size(gauge%matrix,1),size(gauge%matrix,2)))
    do iteration=1,100
      gram=0d0
      call pzgemm('C','N',n,n,n,(1d0,0d0),gauge%matrix,1,1,desc,gauge%matrix,1,1,desc, &
        (0d0,0d0),gram,1,1,desc)
      do j=1,size(gauge%cols)
        do i=1,size(gauge%rows)
          if(gauge%rows(i)>0.and.gauge%rows(i)==gauge%cols(j))gram(i,j)=gram(i,j)-1d0
        enddo
      enddo
      error=sum(abs(gram)**2)
      call comm_summation(error,total_error,gauge%comm)
      if(.not.ieee_is_finite(total_error))return
      if(sqrt(total_error)<=1d-12*sqrt(dble(n)))then
        status=0;return
      endif
      ! X <- X - X (X^H X-I)/2, using distinct GEMM input/output storage.
      next=0d0
      call pzgemm('N','N',n,n,n,(-.5d0,0d0),gauge%matrix,1,1,desc,gram,1,1,desc,(0d0,0d0),next,1,1,desc)
      gauge%matrix=gauge%matrix+next
    enddo
  end subroutine

  subroutine position_seed(gauge,links,desc,status)
    implicit none
    type(s_exx_gauge),intent(inout) :: gauge
    complex(8),intent(in) :: links(:,:,:)
    integer,intent(in) :: desc(9)
    integer,intent(out) :: status
    complex(8),allocatable :: a(:,:),tmp(:,:),work(:),partial(:),column(:)
    real(8),allocatable :: eigen(:),rwork(:)
    complex(8) :: coefficient,query(1),phase
    real(8) :: rquery(1),cw(3),sw(3)
    integer :: n,i,j,k,axis,pivot
    n=gauge%n;cw=sqrt([2d0,3d0,5d0]);sw=sqrt([7d0,11d0,13d0])
    allocate(a(size(links,1),size(links,2)),tmp(size(links,1),size(links,2)),eigen(n))
    a=0d0;tmp=0d0
    do axis=1,3
      coefficient=cmplx(cw(axis),-sw(axis),8)/2d0
      a=a+coefficient*links(:,:,axis)+conjg(coefficient)*links(:,:,axis+3)
    enddo
    call pztranc(n,n,(1d0,0d0),a,1,1,desc,(0d0,0d0),tmp,1,1,desc)
    a=.5d0*(a+tmp)
    call pzheev('V','U',n,a,1,1,desc,eigen,gauge%matrix,1,1,desc,query,-1,rquery,-1,status)
    status=abs(status);call comm_get_max(status,gauge%comm)
    if(status/=0)return
    allocate(work(max(1,int(real(query(1))))),rwork(max(1,int(rquery(1)))))
    call pzheev('V','U',n,a,1,1,desc,eigen,gauge%matrix,1,1,desc,work,size(work),rwork,size(rwork),status)
    status=abs(status)
    if(.not.finite_matrix(gauge%matrix))status=1
    call comm_get_max(status,gauge%comm)
    if(status/=0)return
    allocate(partial(n),column(n))
    do k=1,n
      partial=0d0
      do j=1,size(gauge%cols)
        if(gauge%cols(j)/=k)cycle
        do i=1,size(gauge%rows)
          if(gauge%rows(i)>0)partial(gauge%rows(i))=gauge%matrix(i,j)
        enddo
      enddo
      call comm_summation(partial,column,n,gauge%comm)
      pivot=maxloc(abs(column),dim=1);phase=column(pivot)/abs(column(pivot))
      do j=1,size(gauge%cols)
        if(gauge%cols(j)==k)gauge%matrix(:,j)=gauge%matrix(:,j)*conjg(phase)
      enddo
    enddo
  end subroutine

  subroutine functional(gauge,links,desc,b,weights,spread,gradient,status)
    implicit none
    type(s_exx_gauge),intent(in) :: gauge
    complex(8),intent(in) :: links(:,:,:)
    integer,intent(in) :: desc(9)
    real(8),intent(in) :: b(3,6),weights(6)
    real(8),intent(out) :: spread,gradient
    integer,intent(out) :: status
    complex(8),allocatable :: partial(:,:),diagonal(:,:),q(:,:),back(:,:),adjoint(:,:)
    real(8),allocatable :: theta(:,:),center(:,:)
    real(8) :: residual,norm,total
    integer :: n,i,j,l,ir,jc
    n=gauge%n;status=0;spread=0d0;gradient=huge(1d0)
    allocate(partial(n,6),diagonal(n,6),theta(n,6),center(3,n),q(n,6))
    allocate(back(size(links,1),size(links,2)),adjoint(size(links,1),size(links,2)))
    partial=0d0;center=0d0;back=0d0;adjoint=0d0
    do j=1,size(gauge%cols)
      do i=1,size(gauge%rows)
        ir=gauge%rows(i)
        if(ir>0.and.ir==gauge%cols(j))partial(ir,:)=links(i,j,:)
      enddo
    enddo
    call comm_summation(partial,diagonal,size(partial),gauge%comm)
    if(.not.finite_matrix(diagonal))status=1
    if(any(abs(diagonal)<1d-12))status=1
    if(status/=0)return
    do l=1,6
      do i=1,n
        theta(i,l)=atan2(aimag(diagonal(i,l)),real(diagonal(i,l),8))
        center(:,i)=center(:,i)-weights(l)*b(:,l)*theta(i,l)
      enddo
    enddo
    do l=1,6
      do i=1,n
        residual=theta(i,l)+sum(b(:,l)*center(:,i))
        spread=spread+weights(l)*(1d0-abs(diagonal(i,l))**2+residual**2)
        q(i,l)=weights(l)*(-2d0*conjg(diagonal(i,l))-cmplx(0d0,2d0*residual,8)/diagonal(i,l))
      enddo
      do j=1,size(gauge%cols)
        jc=gauge%cols(j)
        if(jc==0)cycle
        do i=1,size(gauge%rows)
          ir=gauge%rows(i)
          if(ir>0)back(i,j)=back(i,j)-links(i,j,l)*q(jc,l)+q(ir,l)*links(i,j,l)
        enddo
      enddo
    enddo
    call pztranc(n,n,(1d0,0d0),back,1,1,desc,(0d0,0d0),adjoint,1,1,desc)
    norm=sum(abs(.5d0*(back-adjoint))**2)
    call comm_summation(norm,total,gauge%comm)
    gradient=sqrt(total)
    if(.not.ieee_is_finite(gradient).or..not.ieee_is_finite(spread))status=1
  end subroutine

  subroutine minimize(gauge,links,desc,b,weights,maxiter,tolerance,spread,gradient,iterations,status)
    implicit none
    type(s_exx_gauge),intent(inout) :: gauge
    complex(8),intent(inout) :: links(:,:,:)
    integer,intent(in) :: desc(9),maxiter
    real(8),intent(in) :: b(3,6),weights(6),tolerance
    real(8),intent(out) :: spread,gradient
    integer,intent(out) :: iterations,status
    complex(8),allocatable :: partial(:,:,:),pair(:,:,:)
    complex(8) :: sine,local_block(2,2,6),block(2,2,6)
    real(8) :: cosine,rotation(4)
    integer :: n,i,j,a,c,ir,jc,r,col,istat,rank,np
    n=gauge%n;status=1
    call comm_get_groupinfo(gauge%comm,rank,np)
    allocate(partial(n,2,7),pair(n,2,7))
    do iterations=0,maxiter
      call functional(gauge,links,desc,b,weights,spread,gradient,istat)
      if(istat/=0)return
      if(gradient<tolerance.and.iterations>0)then
        status=0;return
      endif
      if(iterations==maxiter)return
      do j=2,n
        do i=1,j-1
          ! Test a pair using only its 2x2 blocks. Skip full row/column traffic
          ! when the existing Jacobi criterion says no rotation is needed.
          local_block=0d0
          do col=1,size(gauge%cols)
            jc=gauge%cols(col)
            if(jc/=i.and.jc/=j)cycle
            c=1
            if(jc==j)c=2
            do r=1,size(gauge%rows)
              ir=gauge%rows(r)
              if(ir/=i.and.ir/=j)cycle
              a=1
              if(ir==j)a=2
              local_block(a,c,:)=links(r,col,:)
            enddo
          enddo
          call comm_summation(local_block,block,size(block),gauge%comm)
          if(rank==0)call jacobi_rotation(block,weights,n,tolerance,rotation)
          call comm_bcast(rotation,gauge%comm,0)
          if(rotation(4)/=0d0)return
          cosine=rotation(1);sine=cmplx(rotation(2),rotation(3),8)
          if(abs(sine)<1d-15)cycle
          ! Columns i/j of six links and U; only O(N) buffers are replicated.
          partial=0d0
          do col=1,size(gauge%cols)
            jc=gauge%cols(col)
            if(jc/=i.and.jc/=j)cycle
            c=1
            if(jc==j)c=2
            do r=1,size(gauge%rows)
              ir=gauge%rows(r)
              if(ir==0)cycle
              partial(ir,c,:6)=links(r,col,:)
              partial(ir,c,7)=gauge%matrix(r,col)
            enddo
          enddo
          call comm_summation(partial,pair,size(pair),gauge%comm)
          do col=1,size(gauge%cols)
            jc=gauge%cols(col)
            if(jc/=i.and.jc/=j)cycle
            do r=1,size(gauge%rows)
              ir=gauge%rows(r)
              if(ir==0)cycle
              if(jc==i)then
                links(r,col,:)=cosine*pair(ir,1,:6)+sine*pair(ir,2,:6)
                gauge%matrix(r,col)=cosine*pair(ir,1,7)+sine*pair(ir,2,7)
              else
                links(r,col,:)=-conjg(sine)*pair(ir,1,:6)+cosine*pair(ir,2,:6)
                gauge%matrix(r,col)=-conjg(sine)*pair(ir,1,7)+cosine*pair(ir,2,7)
              endif
            enddo
          enddo
          ! Left multiplication uses rows after the right rotation above.
          partial=0d0
          do r=1,size(gauge%rows)
            ir=gauge%rows(r)
            if(ir/=i.and.ir/=j)cycle
            a=1
            if(ir==j)a=2
            do col=1,size(gauge%cols)
              jc=gauge%cols(col)
              if(jc>0)partial(jc,a,:6)=links(r,col,:)
            enddo
          enddo
          call comm_summation(partial,pair,size(pair),gauge%comm)
          do r=1,size(gauge%rows)
            ir=gauge%rows(r)
            if(ir/=i.and.ir/=j)cycle
            do col=1,size(gauge%cols)
              jc=gauge%cols(col)
              if(jc==0)cycle
              if(ir==i)then
                links(r,col,:)=cosine*pair(jc,1,:6)+conjg(sine)*pair(jc,2,:6)
              else
                links(r,col,:)=-sine*pair(jc,1,:6)+cosine*pair(jc,2,:6)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine
  subroutine jacobi_rotation(block,weights,n,tolerance,rotation)
    ! One rank chooses the angle so all peers follow the same collective path.
    implicit none
    complex(8),intent(in) :: block(2,2,6)
    real(8),intent(in) :: weights(6),tolerance
    integer,intent(in) :: n
    real(8),intent(out) :: rotation(4)
    complex(8) :: v(3),sine
    real(8) :: metric(3,3),original(3,3),eval(3),work(32),axis(3),scale,cosine
    integer :: l,a,c,status
    metric=0d0;rotation=[1d0,0d0,0d0,0d0]
    do l=1,6
      v(1)=.5d0*(block(1,2,l)+block(2,1,l))
      v(2)=cmplx(0d0,.5d0,8)*(block(1,2,l)-block(2,1,l))
      v(3)=.5d0*(block(1,1,l)-block(2,2,l))
      do c=1,3
        do a=1,3
          metric(a,c)=metric(a,c)+weights(l)*real(conjg(v(a))*v(c),8)
        enddo
      enddo
    enddo
    original=metric
    call dsyev('V','U',3,metric,3,eval,work,size(work),status)
    rotation(4)=dble(status)
    if(status/=0)return
    scale=maxval(abs(original))
    if(eval(3)<=original(3,3)+32*epsilon(1d0)*scale.and. &
       maxval(abs(original(1:2,3)))<tolerance/(16*n))return
    axis=metric(:,3)
    if(axis(3)<0d0)axis=-axis
    cosine=sqrt(.5d0*(1d0+axis(3)));sine=cmplx(axis(1),axis(2),8)/(2*cosine)
    rotation(:3)=[cosine,real(sine,8),aimag(sine)]
  end subroutine
#endif

  subroutine gauge_tiles_rotate(gauge,input,comm_r,comm_o,output,status,adjoint,weights)
    implicit none
    type(s_exx_gauge),intent(in) :: gauge
    complex(8),intent(in) :: input(:,:)
    complex(8),intent(out) :: output(:,:)
    integer,intent(in) :: comm_r,comm_o
    integer,intent(out) :: status
    logical,intent(in),optional :: adjoint
    real(8),intent(in),optional :: weights(:)
    integer,allocatable :: counts(:)
    complex(8),allocatable :: column(:),partial(:),total(:)
    integer :: first,rank,np,owner,k,col,i,j,bad
    logical :: reverse
    reverse=.false.
    if(present(adjoint))reverse=adjoint
    status=1;bad=0;output=0d0
    if(.not.allocated(gauge%matrix))bad=1
    if(any(shape(input)/=shape(output)))bad=1
    if(present(weights))then
      if(size(weights)/=size(input,2))bad=1
    endif
    call orbital_check(bad,comm_r,comm_o)
    if(bad/=0)return
    call orbital_layout(size(input,2),comm_r,comm_o,counts,first,bad)
    if(bad/=0.or.sum(counts)/=gauge%n)return
    call comm_get_groupinfo(comm_o,rank,np)
    allocate(column(size(input,1)),partial(gauge%n),total(gauge%n))
    k=0
    do owner=0,np-1
      do col=1,counts(owner)
        k=k+1
        if(rank==owner)then
          column=input(:,col)
          if(present(weights))column=column*weights(col)
        endif
        call comm_bcast(column,comm_o,owner)
        partial=0d0
        if(reverse)then
          do j=1,size(gauge%cols)
            if(gauge%cols(j)/=k)cycle
            do i=1,size(gauge%rows)
              if(gauge%rows(i)>0)partial(gauge%rows(i))=conjg(gauge%matrix(i,j))
            enddo
          enddo
        else
          do i=1,size(gauge%rows)
            if(gauge%rows(i)/=k)cycle
            do j=1,size(gauge%cols)
              if(gauge%cols(j)>0)partial(gauge%cols(j))=gauge%matrix(i,j)
            enddo
          enddo
        endif
        call comm_summation(partial,total,size(total),gauge%comm)
        do j=1,size(output,2)
          output(:,j)=output(:,j)+column*total(first+j-1)
        enddo
      enddo
    enddo
    status=0
  end subroutine

  logical function finite_matrix(values) result(finite)
    implicit none
    complex(8),intent(in) :: values(:,:)
    real(8) :: component
    integer :: i,j
    finite=.false.
    do j=1,size(values,2)
      do i=1,size(values,1)
        component=real(values(i,j),8)
        if(.not.ieee_is_finite(component))return
        component=aimag(values(i,j))
        if(.not.ieee_is_finite(component))return
      enddo
    enddo
    finite=.true.
  end function
end module
