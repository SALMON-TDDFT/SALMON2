#include "config.h"
! Gamma MLWF matrices in unique block-cyclic tiles. Mesh functions stay local.
module exx_distributed_gauge
  use memory_trace,only:trace_memory
  use small_symmetric3,only:symmetric3_eigen
  use communication, only: comm_get_groupinfo,comm_get_max,comm_summation,comm_bcast
  use exx_orbitals, only: orbital_layout,orbital_check
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: s_exx_gauge,gauge_tiles_clear,gauge_tiles_refresh,gauge_tiles_rotate
#ifdef USE_SCALAPACK
  public :: jacobi_matching,gauge_jacobi_round,jacobi_rotation
#endif
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
    call trace_memory('gauge.enter')
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
    call trace_memory('gauge.setup.done')
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
      call trace_memory('links.allocate.before')
      allocate(links(nr,nc,6),tmp(nr,nc),a(nr,nc))
      ! BLACS does not write padding cells on ranks with partial/empty tiles.
      links=0d0;tmp=0d0;a=0d0
      call trace_memory('links.allocate.after')
      do axis=1,3
        call trace_memory('overlap.before',axis)
        call overlap(gauge,psi,psi,dv,comm_o,counts,first,links(:,:,axis),phase(:,axis))
        call trace_memory('overlap.after',axis)
        call pztranc(n,n,(1d0,0d0),links(:,:,axis),1,1,desc,(0d0,0d0),links(:,:,axis+3),1,1,desc)
      enddo
      if(seed.and.needed)then
        call trace_memory('seed.before')
        call position_seed(gauge,links,desc,bad)
        call trace_memory('seed.after')
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
      call trace_memory('minimize.before')
      call minimize(gauge,links,desc,b,weights,maxiter,tolerance,spread,gradient,iterations,loc_status)
      call trace_memory('minimize.after',iterations)
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
    integer,allocatable :: iwork(:)
    complex(8) :: coefficient,query(1),phase
    real(8) :: rquery(1),cw(3),sw(3)
    integer :: n,i,j,k,axis,pivot,lwork,lrwork,liwork,iquery(1),nprow,npcol,myrow,mycol
    n=gauge%n;cw=sqrt([2d0,3d0,5d0]);sw=sqrt([7d0,11d0,13d0])
    allocate(a(size(links,1),size(links,2)),tmp(size(links,1),size(links,2)),eigen(n))
    a=0d0;tmp=0d0
    do axis=1,3
      coefficient=cmplx(cw(axis),-sw(axis),8)/2d0
      a=a+coefficient*links(:,:,axis)+conjg(coefficient)*links(:,:,axis+3)
    enddo
    call pztranc(n,n,(1d0,0d0),a,1,1,desc,(0d0,0d0),tmp,1,1,desc)
    a=.5d0*(a+tmp)
    call pzheevd('V','U',n,a,1,1,desc,eigen,gauge%matrix,1,1,desc,query,-1,rquery,-1,iquery,-1,status)
    status=abs(status);call comm_get_max(status,gauge%comm)
    if(status/=0)return
    ! PZHEEV's redistribution path corrupts the heap on the failing Fujitsu run.
    ! Use the independently validated divide-and-conquer driver and retain its
    ! queried sizes, with the same conservative local workspace floors as LCFO.
    status=0
    if(.not.ieee_is_finite(real(query(1),8)).or..not.ieee_is_finite(rquery(1)))status=1
    if(real(query(1),8)<1d0.or.rquery(1)<1d0.or.iquery(1)<1)status=1
    if(real(query(1),8)>dble(huge(lwork)-1).or.rquery(1)>dble(huge(lrwork)-1))status=1
    call comm_get_max(status,gauge%comm)
    if(status/=0)return
    lwork=max(1,ceiling(real(query(1),8)))
    lrwork=max(ceiling(rquery(1)),2*(1+8*n+2*size(a,1)*size(a,2)))
    call blacs_gridinfo(desc(2),nprow,npcol,myrow,mycol)
    liwork=max(iquery(1),2+7*n+8*npcol)
    allocate(work(lwork),rwork(lrwork),iwork(liwork))
    call pzheevd('V','U',n,a,1,1,desc,eigen,gauge%matrix,1,1,desc, &
      work,lwork,rwork,lrwork,iwork,liwork,status)
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

  subroutine jacobi_matching(n,round,mate,slot)
    integer,intent(in) :: n,round
    integer,intent(out) :: mate(n),slot(n)
    integer :: m,k,a,b,ia,ib
    m=n+mod(n,2);slot=0
    do k=1,n;mate(k)=k;enddo
    do k=1,m/2
      ia=k;ib=m+1-k
      a=1;b=1
      if(ia>1)a=2+modulo(ia-2-(round-1),m-1)
      if(ib>1)b=2+modulo(ib-2-(round-1),m-1)
      if(a>n.or.b>n)cycle
      mate(a)=b;mate(b)=a;slot(a)=k;slot(b)=k
    enddo
  end subroutine

  subroutine gauge_jacobi_round(gauge,links,desc,weights,tolerance,round,status,times,accepted,row_comm,col_comm)
    use mpi
!$  use omp_lib,only:omp_get_num_threads
    type(s_exx_gauge),intent(inout) :: gauge
    complex(8),intent(inout) :: links(:,:,:)
    integer,intent(in) :: desc(9),round
    real(8),intent(in) :: weights(6),tolerance
    integer,intent(out) :: status
    real(8),intent(out) :: times(4)
    integer(kind=8),intent(out) :: accepted
    integer,intent(in),optional :: row_comm,col_comm
    integer :: n,m,nr,nc,nprow,npcol,myrow,mycol,rank,np,cr,cc,e,k,i,j,a,b,si,sj,istat,l,team
    integer,allocatable :: mate(:),slot(:),sc(:),rc(:),sd(:),rd(:)
    complex(8),allocatable :: local_block(:,:,:,:),block(:,:,:,:),send(:),recv(:)
    real(8),allocatable :: local_rot(:,:),rot(:,:)
    real(8) :: t,cosine
    complex(8) :: sine,factor
    logical :: own_comms
    n=gauge%n;m=n+mod(n,2);nr=size(links,1);nc=size(links,2)
    status=0;times=0;accepted=0;own_comms=.not.present(row_comm)
    call comm_get_groupinfo(gauge%comm,rank,np)
    call blacs_gridinfo(desc(2),nprow,npcol,myrow,mycol)
    if(present(row_comm).neqv.present(col_comm))status=1
    if(n<1.or.round<1.or.round>m-1.or.size(links,3)/=6)status=1
    if(any(shape(gauge%matrix)/=[nr,nc]))status=1
    if(nprow*npcol/=np.or.myrow<0.or.mycol<0.or.desc(5)/=desc(6))status=1
    if(desc(7)/=0.or.desc(8)/=0)status=1
    if(.not.all(ieee_is_finite(weights)).or..not.ieee_is_finite(tolerance))status=1
    if(tolerance<=0d0.or..not.finite_matrix(gauge%matrix))status=1
    do l=1,size(links,3)
      if(.not.finite_matrix(links(:,:,l)))status=1
    enddo
    call comm_get_max(status,gauge%comm)
    if(status/=0)return
    if(own_comms)then
      call MPI_Comm_split(gauge%comm,myrow,mycol,cr,e)
      call jacobi_mpi_fail_stop(e,gauge%comm)
      call MPI_Comm_split(gauge%comm,mycol,myrow,cc,e)
      call jacobi_mpi_fail_stop(e,gauge%comm)
    else
      cr=row_comm;cc=col_comm
    endif
    allocate(mate(n),slot(n),local_block(2,2,6,m/2),block(2,2,6,m/2), &
      local_rot(4,m/2),rot(4,m/2),sc(0:max(nprow,npcol)-1),rc(0:max(nprow,npcol)-1), &
      sd(0:max(nprow,npcol)-1),rd(0:max(nprow,npcol)-1), &
      send(max(1,nr*nc*7)),recv(max(1,nr*nc*7)),stat=istat)
    status=abs(istat);call comm_get_max(status,gauge%comm)
    if(status/=0)goto 900
    call jacobi_matching(n,round,mate,slot)
    t=MPI_Wtime();local_block=0d0;team=1
    !$omp parallel do default(shared) private(j,b,k,sj,i,a,si) reduction(max:team) schedule(static)
    do j=1,nc
!$    team=max(team,omp_get_num_threads())
      b=gauge%cols(j);if(b<1)cycle
      k=slot(b);if(k==0)cycle
      sj=merge(1,2,b<mate(b))
      do i=1,nr
        a=gauge%rows(i);if(a<1)cycle
        if(a/=b.and.a/=mate(b))cycle
        si=merge(1,2,a<mate(a))
        local_block(si,sj,:,k)=links(i,j,:)
      enddo
    enddo
    !$omp end parallel do
    if(rank==0.and.round==1)print *,"MLWF_JACOBI_THREADS",team
    call MPI_Allreduce(local_block,block,size(block),MPI_DOUBLE_COMPLEX,MPI_SUM,gauge%comm,e)
    call jacobi_mpi_fail_stop(e,gauge%comm)
    status=abs(e);call comm_get_max(status,gauge%comm)
    if(status/=0)goto 900
    times(1)=MPI_Wtime()-t;t=MPI_Wtime();local_rot=0d0
    !$omp parallel do default(shared) private(k) schedule(static)
    do k=rank+1,m/2,np
      call jacobi_rotation(block(:,:,:,k),weights,n,tolerance,local_rot(:,k))
      if(sqrt(sum(local_rot(2:3,k)**2))<1d-15)local_rot(:3,k)=[1d0,0d0,0d0]
    enddo
    !$omp end parallel do
    call MPI_Allreduce(local_rot,rot,size(rot),MPI_DOUBLE_PRECISION,MPI_SUM,gauge%comm,e)
    call jacobi_mpi_fail_stop(e,gauge%comm)
    status=abs(e)
    if(any(rot(4,:)/=0d0).or..not.all(ieee_is_finite(rot)))status=1
    call comm_get_max(status,gauge%comm)
    if(status/=0)goto 900
    accepted=int(count(sqrt(sum(rot(2:3,:)**2,dim=1))>=1d-15),8)
    times(2)=MPI_Wtime()-t;t=MPI_Wtime()
    call exchange(1)
    if(status/=0)goto 900
    times(3)=MPI_Wtime()-t;t=MPI_Wtime()
    call exchange(2)
    times(4)=MPI_Wtime()-t
    if(.not.finite_matrix(gauge%matrix))status=1
    do l=1,6
      if(.not.finite_matrix(links(:,:,l)))status=1
    enddo
    call comm_get_max(status,gauge%comm)
900 continue
    if(own_comms)then
      call MPI_Comm_free(cr,e);call jacobi_mpi_fail_stop(e,gauge%comm)
      call MPI_Comm_free(cc,e);call jacobi_mpi_fail_stop(e,gauge%comm)
    endif
  contains
    integer function owner(g,peers)
      integer,intent(in) :: g,peers
      owner=mod((g-1)/desc(5),peers)
    end function
    integer function index_local(g,peers)
      integer,intent(in) :: g,peers
      index_local=((g-1)/(desc(5)*peers))*desc(5)+mod(g-1,desc(5))+1
    end function
    logical function active(g)
      integer,intent(in) :: g
      active=.false.
      if(slot(g)==0)return
      active=abs(cmplx(rot(2,slot(g)),rot(3,slot(g)),8))>=1d-15
    end function
    subroutine exchange(axis)
      integer,intent(in) :: axis
      integer :: peers,coord,comm,stride,ng,ns,d,g,src,dst,idx,ll,v,ii,jj,pos,kk
      integer :: pack_offset(n),unpack_offset(n)
      real(8) :: c_local
      complex(8) :: s_local,f_local
      peers=npcol;coord=mycol;comm=cr;stride=nr*7
      if(axis==2)then;peers=nprow;coord=myrow;comm=cc;stride=nc*6;endif
      sc=0;rc=0;sd=0;rd=0
      do g=1,n
        if(.not.active(g))cycle
        src=owner(mate(g),peers);dst=owner(g,peers)
        if(src==coord)sc(dst)=sc(dst)+stride
        if(dst==coord)rc(src)=rc(src)+stride
      enddo
      do d=1,peers-1
        sd(d)=sd(d-1)+sc(d-1);rd(d)=rd(d-1)+rc(d-1)
      enddo
      status=0
      if(sum(sc(:peers-1))>size(send).or.sum(rc(:peers-1))>size(recv))status=1
      call comm_get_max(status,gauge%comm)
      if(status/=0)return
      pack_offset=-1;unpack_offset=-1;pos=0
      do d=0,peers-1
        do g=1,n
          if(.not.active(g))cycle
          if(owner(g,peers)/=d.or.owner(mate(g),peers)/=coord)cycle
          pack_offset(g)=pos;pos=pos+stride
        enddo
      enddo
      pos=0
      do d=0,peers-1
        do g=1,n
          if(.not.active(g))cycle
          if(owner(g,peers)/=coord.or.owner(mate(g),peers)/=d)cycle
          unpack_offset(g)=pos;pos=pos+stride
        enddo
      enddo
      !$omp parallel do default(shared) private(g,idx,ll,ii,jj,pos) schedule(static)
      do g=1,n
          if(pack_offset(g)<0)cycle
          pos=pack_offset(g)
          idx=index_local(mate(g),peers)
          if(axis==1)then
            do ll=1,7;do ii=1,nr
              pos=pos+1
              if(ll<=6)then;send(pos)=links(ii,idx,ll)
              else;send(pos)=gauge%matrix(ii,idx);endif
            enddo;enddo
          else
            do ll=1,6;do jj=1,nc
              pos=pos+1;send(pos)=links(idx,jj,ll)
            enddo;enddo
          endif
      enddo
      !$omp end parallel do
      call MPI_Alltoallv(send,sc(:peers-1),sd(:peers-1),MPI_DOUBLE_COMPLEX, &
        recv,rc(:peers-1),rd(:peers-1),MPI_DOUBLE_COMPLEX,comm,e)
      call jacobi_mpi_fail_stop(e,gauge%comm)
      status=abs(e);call comm_get_max(status,gauge%comm)
      if(status/=0)return
      !$omp parallel do default(shared) private(g,idx,kk,c_local,s_local,f_local,ll,ii,jj,pos) schedule(static)
      do g=1,n
          if(unpack_offset(g)<0)cycle
          pos=unpack_offset(g)
          idx=index_local(g,peers);kk=slot(g);c_local=rot(1,kk)
          s_local=cmplx(rot(2,kk),rot(3,kk),8)
          if(axis==1)then
            f_local=-conjg(s_local);if(g<mate(g))f_local=s_local
            do ll=1,7;do ii=1,nr
              pos=pos+1
              if(ll<=6)then;links(ii,idx,ll)=c_local*links(ii,idx,ll)+f_local*recv(pos)
              else;gauge%matrix(ii,idx)=c_local*gauge%matrix(ii,idx)+f_local*recv(pos);endif
            enddo;enddo
          else
            f_local=-s_local;if(g<mate(g))f_local=conjg(s_local)
            do ll=1,6;do jj=1,nc
              pos=pos+1;links(idx,jj,ll)=c_local*links(idx,jj,ll)+f_local*recv(pos)
            enddo;enddo
          endif
      enddo
      !$omp end parallel do
    end subroutine
  end subroutine

  subroutine jacobi_mpi_fail_stop(error,comm)
    use mpi
    integer,intent(in) :: error,comm
    integer :: ignored
    ! Also fail closed if the caller selected MPI_ERRORS_RETURN. Never use an
    ! invalid communicator or an undefined collective output after an error.
    if(error/=MPI_SUCCESS)call MPI_Abort(comm,max(1,abs(error)),ignored)
  end subroutine

  subroutine minimize(gauge,links,desc,b,weights,maxiter,tolerance,spread,gradient,iterations,status)
    use mpi
    type(s_exx_gauge),intent(inout) :: gauge
    complex(8),intent(inout) :: links(:,:,:)
    integer,intent(in) :: desc(9),maxiter
    real(8),intent(in) :: b(3,6),weights(6),tolerance
    real(8),intent(out) :: spread,gradient
    integer,intent(out) :: iterations,status
    integer :: rd,istat,nprow,npcol,myrow,mycol,cr,cc,e,rank,np
    integer(kind=8) :: accepted,total_accepted
    real(8) :: times(4),totals(4),maximum(4),t,wall,maxwall,orth,maxorth
    complex(8),allocatable :: gram(:,:)
    integer :: i,j
    call comm_get_groupinfo(gauge%comm,rank,np)
    call blacs_gridinfo(desc(2),nprow,npcol,myrow,mycol)
    call MPI_Comm_split(gauge%comm,myrow,mycol,cr,e)
    call jacobi_mpi_fail_stop(e,gauge%comm)
    call MPI_Comm_split(gauge%comm,mycol,myrow,cc,e)
    call jacobi_mpi_fail_stop(e,gauge%comm)
    status=1
    do iterations=0,maxiter
      call functional(gauge,links,desc,b,weights,spread,gradient,istat)
      if(rank==0)then
        print *,'MLWF_PARALLEL_FUNCTIONAL',iterations,spread,gradient,tolerance
        flush(6)
      endif
      if(istat/=0)exit
      if(gradient<tolerance.and.iterations>0)then;status=0;exit;endif
      if(iterations==maxiter)exit
      t=MPI_Wtime();totals=0;total_accepted=0
      do rd=1,gauge%n+mod(gauge%n,2)-1
        call gauge_jacobi_round(gauge,links,desc,weights,tolerance,rd,istat,times,accepted,cr,cc)
        if(istat/=0)exit
        totals=totals+times;total_accepted=total_accepted+accepted
      enddo
      if(istat/=0)exit
      wall=MPI_Wtime()-t
      call MPI_Allreduce(totals,maximum,4,MPI_DOUBLE_PRECISION,MPI_MAX,gauge%comm,e)
      call jacobi_mpi_fail_stop(e,gauge%comm)
      call MPI_Allreduce(wall,maxwall,1,MPI_DOUBLE_PRECISION,MPI_MAX,gauge%comm,e)
      call jacobi_mpi_fail_stop(e,gauge%comm)
      if(rank==0)then
        print *,'MLWF_PARALLEL_SWEEP',iterations,total_accepted,maxwall,maximum
        flush(6)
      endif
    enddo
    if(status==0)then
      allocate(gram(size(gauge%rows),size(gauge%cols)))
      gram=0d0
      call pzgemm('C','N',gauge%n,gauge%n,gauge%n,(1d0,0d0),gauge%matrix,1,1,desc, &
        gauge%matrix,1,1,desc,(0d0,0d0),gram,1,1,desc)
      do j=1,size(gauge%cols);do i=1,size(gauge%rows)
        if(gauge%rows(i)>0.and.gauge%rows(i)==gauge%cols(j))gram(i,j)=gram(i,j)-1d0
      enddo;enddo
      orth=maxval(abs(gram));call MPI_Allreduce(orth,maxorth,1,MPI_DOUBLE_PRECISION,MPI_MAX,gauge%comm,e)
      call jacobi_mpi_fail_stop(e,gauge%comm)
      if(rank==0)print *,'MLWF_PARALLEL_ORTHOGONAL',maxorth
      if(.not.ieee_is_finite(maxorth).or.maxorth>1d-10)status=1
    endif
    call MPI_Comm_free(cr,e);call jacobi_mpi_fail_stop(e,gauge%comm)
    call MPI_Comm_free(cc,e);call jacobi_mpi_fail_stop(e,gauge%comm)
  end subroutine

  recursive subroutine jacobi_rotation(block,weights,n,tolerance,rotation)
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
    call symmetric3_eigen(original,metric,eval,status)
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
    integer :: first,rank,np,owner,k,col,i,j,g,bad
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
        ! Preserve the input-column summation order and the original mesh ownership.
!$omp parallel do collapse(2) default(none) schedule(static) &
!$omp shared(output,column,total,first) private(j,g)
        do j=1,size(output,2)
          do g=1,size(output,1)
            output(g,j)=output(g,j)+column(g)*total(first+j-1)
          enddo
        enddo
!$omp end parallel do
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
