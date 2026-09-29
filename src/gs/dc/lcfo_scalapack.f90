! Full complex LCFO eigensolver: global dense matrices stay block-cyclic.
module lcfo_scalapack
  use communication, only: comm_get_groupinfo,comm_summation,comm_get_max
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: lcfo_dense_state,lcfo_dense_init,lcfo_dense_add,lcfo_dense_solve,lcfo_dense_rows,lcfo_dense_free
  type lcfo_dense_state
    integer :: context=-1,comm,n,nb,nprow,npcol,myrow,mycol,nr,nc,desc(9)
    integer,allocatable :: rows(:),cols(:)
    complex(8),allocatable :: h(:,:),vectors(:,:)
    real(8),allocatable :: values(:)
  end type
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_1d,finite_real_2d
  interface salmon_all_finite
    module procedure finite_real_1d,finite_real_2d
  end interface
contains
  subroutine max_real(value,comm)
    real(8),intent(inout) :: value
    integer,intent(in) :: comm
    real(8) :: result(1)
    call comm_get_max([value],result,1,comm)
    value=result(1)
  end subroutine
  subroutine lcfo_dense_init(s,n,comm,status)
    type(lcfo_dense_state),intent(out) :: s
    integer,intent(in) :: n,comm
    integer,intent(out) :: status
    integer :: rank,np,iam,nworld,i,j,numroc
    integer,allocatable :: ranks(:),global_ranks(:),map(:,:)
    external :: numroc
    s%comm=comm;s%n=n;status=0
    call comm_get_groupinfo(comm,rank,np)
    call blacs_pinfo(iam,nworld)
    if(n<1.or.nworld<np)status=1
    call comm_get_max(status,comm)
    if(status/=0)return
    s%nprow=int(sqrt(real(np,8)))
    do while(mod(np,s%nprow)/=0)
      s%nprow=s%nprow-1
    enddo
    s%npcol=np/s%nprow
    allocate(ranks(np),global_ranks(np),map(s%nprow,s%npcol));ranks=0
    ranks(rank+1)=iam+1
    call comm_summation(ranks,global_ranks,np,comm)
    map=reshape(global_ranks-1,shape(map))
    call blacs_get(0,0,s%context)
    call blacs_gridmap(s%context,map,s%nprow,s%nprow,s%npcol)
    call blacs_gridinfo(s%context,s%nprow,s%npcol,s%myrow,s%mycol)
    s%nb=min(32,max(1,n/max(s%nprow,s%npcol)))
    s%nr=numroc(n,s%nb,s%myrow,0,s%nprow);s%nc=numroc(n,s%nb,s%mycol,0,s%npcol)
    call descinit(s%desc,n,n,s%nb,s%nb,0,0,s%context,max(1,s%nr),status)
    status=abs(status)
    call comm_get_max(status,comm)
    if(status/=0)return
    allocate(s%h(max(1,s%nr),max(1,s%nc)),s%vectors(max(1,s%nr),max(1,s%nc)),s%values(n))
    allocate(s%rows(s%nr),s%cols(s%nc))
    do i=1,s%nr
      s%rows(i)=((i-1)/s%nb*s%nprow+s%myrow)*s%nb+mod(i-1,s%nb)+1
    enddo
    do j=1,s%nc
      s%cols(j)=((j-1)/s%nb*s%npcol+s%mycol)*s%nb+mod(j-1,s%nb)+1
    enddo
    s%h=0d0;s%vectors=0d0;s%values=0d0
    write(*,'(a,4i12)')'LCFO_DENSE rank/local_rows/local_cols/global: ',rank,s%nr,s%nc,n
  end subroutine

  subroutine lcfo_dense_add(s,first_row,first_col,block)
    type(lcfo_dense_state),intent(inout) :: s
    integer,intent(in) :: first_row,first_col
    complex(8),intent(in) :: block(:,:)
    integer :: i,j,ib,jb
    do j=1,s%nc
      jb=s%cols(j)-first_col+1
      if(jb<1.or.jb>size(block,2))cycle
      do i=1,s%nr
        ib=s%rows(i)-first_row+1
        if(ib<1.or.ib>size(block,1))cycle
        s%h(i,j)=s%h(i,j)+block(ib,jb)
      enddo
    enddo
  end subroutine

  subroutine lcfo_dense_solve(s,nt,hermitian,orthogonal,residual,status)
    type(lcfo_dense_state),intent(inout) :: s
    integer,intent(in) :: nt
    real(8),intent(out) :: hermitian,orthogonal,residual
    integer,intent(out) :: status
    complex(8),allocatable :: original(:,:),work(:)
    real(8),allocatable :: rwork(:),local_res(:),total_res(:)
    complex(8) :: query(1)
    real(8) :: rquery(1),scale,local_norm,norm
    integer :: i,j,n,lwork,lrwork
    complex(8),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    n=s%n;status=0;hermitian=huge(1d0);orthogonal=huge(1d0);residual=huge(1d0)
    if(nt<1.or.nt>n)status=1
    if(.not.salmon_all_finite(real(s%h)).or..not.salmon_all_finite(aimag(s%h)))status=1
    call comm_get_max(status,s%comm)
    if(status/=0)return
    call pztranc(n,n,one,s%h,1,1,s%desc,zero,s%vectors,1,1,s%desc)
    hermitian=maxval(abs(s%h-s%vectors));scale=maxval(abs(s%h))
    call max_real(hermitian,s%comm);call max_real(scale,s%comm)
    hermitian=hermitian/max(1d0,scale)
    if(hermitian>1d-10.or..not.ieee_is_finite(hermitian))then
      status=1;return
    endif
    s%h=.5d0*(s%h+s%vectors)
    if(n==1)then
      ! Exact scalar eigensystem, including ranks with empty local tiles.
      ! Avoid the vendor PZHEEV scalar shortcut returning a nonunit vector.
      local_norm=0d0;s%vectors=0d0
      if(s%nr==1.and.s%nc==1)then
        local_norm=real(s%h(1,1),8);s%vectors(1,1)=one
      endif
      call comm_summation(local_norm,norm,s%comm)
      s%values(1)=norm;orthogonal=0d0;residual=0d0;return
    endif
    original=s%h
    call pzheev('V','U',n,s%h,1,1,s%desc,s%values,s%vectors,1,1,s%desc,query,-1,rquery,-1,status)
    status=abs(status)
    call comm_get_max(status,s%comm)
    if(status/=0)return
    lwork=max(1,nint(real(query(1),8)));lrwork=max(1,nint(rquery(1)))
    allocate(work(lwork),rwork(lrwork))
    call pzheev('V','U',n,s%h,1,1,s%desc,s%values,s%vectors,1,1,s%desc,work,lwork,rwork,lrwork,status)
    status=abs(status)
    call comm_get_max(status,s%comm)
    if(status/=0)return
    deallocate(work,rwork)
    ! Validate all returned eigenvectors, matching the LAPACK reference contract.
    call pzgemm('C','N',n,n,n,one,s%vectors,1,1,s%desc,s%vectors,1,1,s%desc,zero,s%h,1,1,s%desc)
    orthogonal=0d0
    do j=1,s%nc;do i=1,s%nr
      if(s%rows(i)==s%cols(j))s%h(i,j)=s%h(i,j)-one
      orthogonal=max(orthogonal,abs(s%h(i,j)))
    enddo;enddo
    call max_real(orthogonal,s%comm)
    call pzgemm('N','N',n,nt,n,one,original,1,1,s%desc,s%vectors,1,1,s%desc,zero,s%h,1,1,s%desc)
    allocate(local_res(nt),total_res(nt));local_res=0d0
    do j=1,s%nc
      if(s%cols(j)>nt)cycle
      do i=1,s%nr
        local_res(s%cols(j))=local_res(s%cols(j))+ &
          abs(s%h(i,j)-s%values(s%cols(j))*s%vectors(i,j))**2
      enddo
    enddo
    call comm_summation(local_res,total_res,nt,s%comm)
    local_norm=sum(abs(original(:s%nr,:s%nc))**2)
    call comm_summation(local_norm,norm,s%comm)
    residual=maxval(sqrt(total_res)/max(1d0,sqrt(norm),abs(s%values(:nt))))
    if(.not.salmon_all_finite(s%values).or..not.ieee_is_finite(residual).or. &
       .not.ieee_is_finite(orthogonal).or.max(orthogonal,residual)>1d-10)status=1
    call comm_get_max(status,s%comm)
  end subroutine

  subroutine lcfo_dense_rows(s,first,count,nt,root,output)
    type(lcfo_dense_state),intent(in) :: s
    integer,intent(in) :: first,count,nt,root
    complex(8),intent(out) :: output(:,:)
    complex(8),allocatable :: local(:,:)
    integer :: i,j,row
    allocate(local(count,nt));local=0d0;output=0d0
    do j=1,s%nc
      if(s%cols(j)>nt)cycle
      do i=1,s%nr
        row=s%rows(i)-first+1
        if(row>=1.and.row<=count)local(row,s%cols(j))=s%vectors(i,j)
      enddo
    enddo
    call comm_summation(local,output,size(local),s%comm,root)
  end subroutine

  subroutine lcfo_dense_free(s)
    type(lcfo_dense_state),intent(inout) :: s
    if(allocated(s%h))deallocate(s%h,s%vectors,s%values,s%rows,s%cols)
    if(s%context>=0)call blacs_gridexit(s%context)
    s%context=-1
  end subroutine

  pure logical function finite_real_1d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:)
    real(8) :: value
    integer :: i
    finite=.false.
    do i=1,size(values,1)
      value=values(i)
      if(.not.ieee_is_finite(value))return
    enddo
    finite=.true.
  end function

  pure logical function finite_real_2d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:)
    real(8) :: value
    integer :: i,j
    finite=.false.
    do j=1,size(values,2)
      do i=1,size(values,1)
        value=values(i,j)
        if(.not.ieee_is_finite(value))return
      enddo
    enddo
    finite=.true.
  end function

end module
