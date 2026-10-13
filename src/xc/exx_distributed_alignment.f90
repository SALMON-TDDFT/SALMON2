#include "config.h"
! Align real-space ACE columns without replicating an orbital-space matrix.
module exx_distributed_alignment
  use communication, only: comm_get_groupinfo,comm_bcast,comm_summation,comm_get_max
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: distributed_factor_align
contains
  subroutine distributed_factor_align(x,ref,comm_r,comm_o,comm,status)
    use exx_blas_threads, only: exx_blas_thread_control
!$  use omp_lib, only: omp_get_max_threads,omp_in_parallel
    implicit none
    complex(8),intent(inout) :: x(:,:)
    complex(8),intent(in) :: ref(:,:)
    integer,intent(in) :: comm_r,comm_o,comm
    integer,intent(out) :: status
    integer :: workers,state(3)
    workers=1;state=0
!$  workers=omp_get_max_threads()
!$  if(omp_in_parallel())workers=1
    if(workers>1)call exx_blas_thread_control(workers,state)
    call factor_align_core(x,ref,comm_r,comm_o,comm,status)
    call exx_blas_thread_control(0,state)
  end subroutine

  subroutine factor_align_core(x,ref,comm_r,comm_o,comm,status)
    implicit none
    complex(8),intent(inout) :: x(:,:)
    complex(8),intent(in) :: ref(:,:)
    integer,intent(in) :: comm_r,comm_o,comm
    integer,intent(out) :: status
#ifdef USE_SCALAPACK
    integer :: rank,np,rank_o,np_o,rank_r,np_r,n,first,owner,col,k,i,j,ir,jc,bad
    integer :: context,nprow,npcol,myrow,mycol,block,nrows,ncols,desc(9),lwork,svd_status
    integer,external :: sys2blacs_handle,numroc
    integer,allocatable :: counts(:),local_counts(:),rows(:),cols(:)
    complex(8),allocatable :: a(:,:),u(:,:),vt(:,:),rotation(:,:),work(:),column(:),partial(:),total(:),aligned(:,:),saved(:,:)
    real(8),allocatable :: singular(:),rwork(:)
    real(8) :: local_norm(4),global_norm(4)
    complex(8) :: query(1)
    status=1
    call comm_get_groupinfo(comm,rank,np)
    call comm_get_groupinfo(comm_o,rank_o,np_o)
    call comm_get_groupinfo(comm_r,rank_r,np_r)
    bad=0
    if(any(shape(x)/=shape(ref)))bad=1
    ! Reject collectively before accessing ref with x's bounds.
    call comm_get_max(bad,comm)
    if(bad/=0)return
    do j=1,size(x,2)
      do i=1,size(x,1)
        if(.not.ieee_is_finite(real(x(i,j),8)).or..not.ieee_is_finite(aimag(x(i,j))))bad=1
        if(.not.ieee_is_finite(real(ref(i,j),8)).or..not.ieee_is_finite(aimag(ref(i,j))))bad=1
      enddo
    enddo
    call comm_get_max(bad,comm)
    if(bad/=0)return
    allocate(counts(0:np_o-1),local_counts(0:np_o-1));local_counts=0
    ! Each orbital column has identical ownership on all spatial ranks.
    if(rank_r==0)local_counts(rank_o)=size(x,2)
    call comm_summation(local_counts,counts,np_o,comm)
    n=sum(counts)
    if(counts(rank_o)/=size(x,2))bad=1
    call comm_get_max(bad,comm)
    if(n<1.or.bad/=0)return
    first=1+sum(counts(:rank_o-1))
    nprow=int(sqrt(dble(np)))
    do while(mod(np,nprow)/=0)
      nprow=nprow-1
    enddo
    npcol=np/nprow
    context=sys2blacs_handle(comm)
    call blacs_gridinit(context,'R',nprow,npcol)
    call blacs_gridinfo(context,nprow,npcol,myrow,mycol)
    block=min(32,max(1,n/max(nprow,npcol)))
    nrows=max(1,numroc(n,block,myrow,0,nprow));ncols=max(1,numroc(n,block,mycol,0,npcol))
    call descinit(desc,n,n,block,block,0,0,context,nrows,bad)
    bad=abs(bad);call comm_get_max(bad,comm)
    if(bad/=0)then
      call blacs_gridexit(context)
      return
    endif
    allocate(rows(nrows),cols(ncols));rows=0;cols=0
    do i=1,n
      if(mod((i-1)/block,nprow)==myrow)then
        ir=((i-1)/(block*nprow))*block+mod(i-1,block)+1
        rows(ir)=i
      endif
      if(mod((i-1)/block,npcol)==mycol)then
        jc=((i-1)/(block*npcol))*block+mod(i-1,block)+1
        cols(jc)=i
      endif
    enddo
    allocate(a(nrows,ncols),u(nrows,ncols),vt(nrows,ncols),rotation(nrows,ncols))
    allocate(column(size(x,1)),partial(n),total(n),singular(n),rwork(1+4*n))
    a=0;u=0;vt=0;rotation=0;k=0
    do owner=0,np_o-1
      do col=1,counts(owner)
        k=k+1
        if(rank_o==owner)column=ref(:,col)
        call comm_bcast(column,comm_o,owner)
        partial=0
        do i=1,size(x,2)
          partial(first+i-1)=dot_product(x(:,i),column)
        enddo
        call comm_summation(partial,total,n,comm)
        do j=1,ncols
          if(cols(j)/=k)cycle
          do i=1,nrows
            if(rows(i)>0)a(i,j)=total(rows(i))
          enddo
        enddo
      enddo
    enddo
    allocate(saved(nrows,ncols));saved=a
    call pzgesvd('V','V',n,n,a,1,1,desc,singular,u,1,1,desc,vt,1,1,desc,query,-1,rwork,bad)
    bad=abs(bad)
    if(.not.ieee_is_finite(real(query(1),8)))bad=1
    call comm_get_max(bad,comm)
    if(bad==0)then
      lwork=max(1,int(real(query(1),8)));allocate(work(lwork))
      call pzgesvd('V','V',n,n,a,1,1,desc,singular,u,1,1,desc,vt,1,1,desc,work,lwork,rwork,bad)
      svd_status=bad
      ! PZGESVD flags any bitwise difference in singular values as n+1.
      ! Such output is accepted only after explicit residual/unitarity checks.
      if(svd_status==n+1)bad=0
      bad=abs(bad)
      ! A local solver error must stop every rank before the next collective.
      call comm_get_max(bad,comm)
      if(bad/=0)then
        call blacs_gridexit(context)
        return
      endif
      if(bad==0)then
        call comm_bcast(singular,comm,0)
        a=u
        do j=1,ncols
          if(cols(j)>0)a(:,j)=a(:,j)*singular(cols(j))
        enddo
        call pzgemm('N','N',n,n,n,(1d0,0d0),a,1,1,desc,vt,1,1,desc,(0d0,0d0),rotation,1,1,desc)
        local_norm(1)=sum(abs(saved)**2)
        local_norm(2)=sum(abs(saved-rotation)**2)
        call pzgemm('C','N',n,n,n,(1d0,0d0),u,1,1,desc,u,1,1,desc,(0d0,0d0),a,1,1,desc)
        call pzgemm('N','C',n,n,n,(1d0,0d0),vt,1,1,desc,vt,1,1,desc,(0d0,0d0),rotation,1,1,desc)
        do j=1,ncols
          do i=1,nrows
            if(rows(i)>0.and.rows(i)==cols(j))then
              a(i,j)=a(i,j)-1d0
              rotation(i,j)=rotation(i,j)-1d0
            endif
          enddo
        enddo
        local_norm(3)=sum(abs(a)**2)
        local_norm(4)=sum(abs(rotation)**2)
        call comm_summation(local_norm,global_norm,4,comm)
        do i=1,4
          if(.not.ieee_is_finite(global_norm(i)))bad=1
        enddo
        if(global_norm(1)<=0d0.or.sqrt(global_norm(2))>1d-10*sqrt(global_norm(1)))bad=1
        if(max(global_norm(3),global_norm(4))>1d-20*dble(n))bad=1
        if(svd_status==n+1.and.rank==0)write(*,'(a,i6,3es15.6)') &
          'ACE alignment SVD heterogeneity/residual/unitarity: ',bad, &
          sqrt(global_norm(2)/max(global_norm(1),tiny(1d0))),sqrt(global_norm(3)),sqrt(global_norm(4))
      endif
      do i=1,n
        if(.not.ieee_is_finite(singular(i)))bad=1
      enddo
      if(minval(singular)<=maxval(singular)*1d-12)then
        if(rank==0)write(*,'(a,2es20.10)')'ACE distributed alignment singular range: ', &
          minval(singular),maxval(singular)
        bad=1
      endif
      call comm_get_max(bad,comm)
      if(bad==0)then
        call pzgemm('N','N',n,n,n,(1d0,0d0),u,1,1,desc,vt,1,1,desc,(0d0,0d0),rotation,1,1,desc)
        allocate(aligned(size(x,1),size(x,2)));aligned=0;k=0
        ! Rotation rows are unique across the BLACS grid. Only an O(n) row is reduced.
        do owner=0,np_o-1
          do col=1,counts(owner)
            k=k+1
            if(rank_o==owner)column=x(:,col)
            call comm_bcast(column,comm_o,owner)
            partial=0
            do i=1,nrows
              if(rows(i)/=k)cycle
              do j=1,ncols
                if(cols(j)>0)partial(cols(j))=rotation(i,j)
              enddo
            enddo
            call comm_summation(partial,total,n,comm)
            do j=1,size(x,2)
              aligned(:,j)=aligned(:,j)+column*total(first+j-1)
            enddo
          enddo
        enddo
        do j=1,size(aligned,2)
          do i=1,size(aligned,1)
            if(.not.ieee_is_finite(real(aligned(i,j),8)).or..not.ieee_is_finite(aimag(aligned(i,j))))bad=1
          enddo
        enddo
        call comm_get_max(bad,comm)
        if(bad==0)x=aligned
      endif
    endif
    call blacs_gridexit(context)
    if(bad==0)status=0
#else
    status=1
#endif
  end subroutine
end module
