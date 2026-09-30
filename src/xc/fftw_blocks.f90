! Axis transforms return every field to its original Cartesian block.
! Only axis-line segments move; no rank holds a replicated three-dimensional box.
module fftw_blocks
  use iso_c_binding
  use iso_fortran_env, only: int64
  use communication, only: comm_alltoall,comm_get_max,comm_get_groupinfo
  implicit none
  private
  include 'fftw3.f03'
  public :: block_layout,block_transform,block_clear,mesh_transform
  integer,parameter :: local_batch=4
  type axis_workspace
    integer :: length=0,segment=0,lines=0,peers=0
    complex(c_double_complex),allocatable :: send(:),recv(:),work(:)
    type(c_ptr) :: forward=c_null_ptr,backward=c_null_ptr
    type(c_ptr) :: batch_forward=c_null_ptr,batch_backward=c_null_ptr
  end type
  type(axis_workspace),save :: cache(3)
  real(8),save,public :: fftw_block_seconds(4)=0d0 ! setup, FFT, MPI, packing/copies
  real(8),save,public :: fftw_block_fused_seconds=0d0 ! local gather + FFT + scatter
contains
  real(8) function stamp()
    implicit none
    integer(int64) :: count,rate
    call system_clock(count,rate)
    stamp=real(count,8)/real(rate,8)
  end function
  subroutine mesh_transform(n,dims,coords,comm,input,output,sign,status,spectral_z)
    use fftw_pencils, only: pencil_transform
    implicit none
    integer,intent(in) :: n(3),dims(:),coords(:),comm(:),sign
    complex(8),intent(in) :: input(:,:)
    complex(8),intent(out) :: output(:,:)
    integer,intent(out) :: status
    logical,intent(in),optional :: spectral_z
    status=1
    if(size(dims)/=size(coords).or.size(dims)/=size(comm))return
    select case(size(dims))
    case(2)
      call pencil_transform(n,dims,coords,comm,input,output,sign,status,spectral_z)
    case(3)
      ! Spectra share Cartesian ownership and xyz order with real-space fields.
      call block_transform(n,dims,coords,comm,input,output,sign,status)
    end select
  end subroutine
  subroutine block_layout(n,dims,coords,m,lo,status)
    implicit none
    integer,intent(in) :: n(3),dims(:),coords(:)
    integer,intent(out) :: m(3),lo(3),status
    integer :: d(3),c(3)
    status=1;m=0;lo=0
    if(size(dims)/=size(coords))return
    select case(size(dims))
    case(2)
      d=[1,dims];c=[0,coords]
    case(3)
      d=dims;c=coords
    case default
      return
    end select
    if(any(n<1).or.any(d<1))return
    if(any(modulo(n,d)/=0).or.any(c<0).or.any(c>=d))return
    m=n/d;lo=c*m;status=0
  end subroutine
  subroutine release(p)
    implicit none
    type(axis_workspace),intent(inout) :: p
    if(c_associated(p%forward))call fftw_destroy_plan(p%forward)
    if(c_associated(p%backward))call fftw_destroy_plan(p%backward)
    if(c_associated(p%batch_forward))call fftw_destroy_plan(p%batch_forward)
    if(c_associated(p%batch_backward))call fftw_destroy_plan(p%batch_backward)
    p%batch_forward=c_null_ptr;p%batch_backward=c_null_ptr
    p%forward=c_null_ptr;p%backward=c_null_ptr
    if(allocated(p%work))deallocate(p%work,p%send,p%recv)
    p%length=0;p%segment=0;p%lines=0;p%peers=0
  end subroutine
  subroutine block_clear()
    implicit none
    integer :: a
    do a=1,3
      call release(cache(a))
    enddo
  end subroutine
  subroutine prepare(p,n,m,lines,peers,status)
    implicit none
    type(axis_workspace),intent(inout) :: p
    integer,intent(in) :: n,m,lines,peers
    integer,intent(out) :: status
    integer(int64) :: count
    status=0
    if(p%length==n.and.p%segment==m.and.p%lines==lines.and.p%peers==peers)return
    call release(p)
    count=int(n,int64)*lines
    if(count>int(huge(0),int64))then
      status=1;return
    endif
    allocate(p%work(int(count)),p%send(int(count)),p%recv(int(count)))
    p%work=0d0
    p%forward=fftw_plan_dft_1d(n,p%work,p%work,FFTW_FORWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    p%backward=fftw_plan_dft_1d(n,p%work,p%work,FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    if(.not.c_associated(p%forward).or..not.c_associated(p%backward))then
      call release(p);status=1;return
    endif
    if(peers==1.and.lines>=local_batch)then
      p%batch_forward=fftw_plan_many_dft(1,[n],local_batch,p%work,[n],1,n,p%work,[n],1,n, &
        FFTW_FORWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
      p%batch_backward=fftw_plan_many_dft(1,[n],local_batch,p%work,[n],1,n,p%work,[n],1,n, &
        FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
      if(.not.c_associated(p%batch_forward).or..not.c_associated(p%batch_backward))then
        call release(p);status=1;return
      endif
    endif
    p%length=n;p%segment=m;p%lines=lines;p%peers=peers
  end subroutine
  subroutine block_transform(n,dims,coords,comm,input,output,sign,status)
    implicit none
    integer,intent(in) :: n(3),dims(3),coords(3),comm(3),sign
    complex(8),intent(in) :: input(:,:)
    complex(8),intent(out) :: output(:,:)
    integer,intent(out) :: status
    integer :: m(3),lo(3),a,b,bad,rank,peers,nt,ntmax
    call block_layout(n,dims,coords,m,lo,bad)
    if(sign/=-1.and.sign/=1)bad=1
    if(any(shape(input)/=shape(output)))bad=1
    if(size(input,1)/=product(m))bad=1
    nt=size(input,2)
    do a=1,3
      call comm_get_groupinfo(comm(a),rank,peers)
      if(rank/=coords(a).or.peers/=dims(a))bad=1
      ntmax=nt;call comm_get_max(ntmax,comm(a))
      if(ntmax/=nt)bad=1
    enddo
    do a=1,3
      call comm_get_max(bad,comm(a))
    enddo
    status=bad
    if(status/=0)return
    output=input
    if(nt==0)return
    do a=1,3
      call axis_transform(cache(a),n,m,a,dims(a),comm(a),output,sign,bad)
      do b=1,3
        call comm_get_max(bad,comm(b))
      enddo
      if(bad/=0)then
        status=bad;return
      endif
    enddo
    if(sign==1)output=output/real(product(int(n,int64)),8)
  end subroutine
  subroutine axis_transform(p,n,m,axis,peers,comm,field,sign,status)
    implicit none
    type(axis_workspace),intent(inout) :: p
    integer,intent(in) :: n(3),m(3),axis,peers,comm,sign
    complex(8),intent(inout) :: field(:,:)
    integer,intent(out) :: status
    integer :: cross,lines,nowned,count,ell,owner,slot,i,j,t,b,g,point(3),other(2),s,base,full_slots,tail_lines,q,stride
    integer :: first,last,batch_base
    type(c_ptr) :: plan,batch_plan
    integer(int64) :: total_lines
    real(8) :: started
    other=pack([1,2,3],[1,2,3]/=axis)
    total_lines=int(m(other(1)),int64)*m(other(2))*size(field,2)
    if(total_lines>int(huge(0),int64))then
      status=1;return
    endif
    cross=m(other(1))*m(other(2));lines=int(total_lines)
    full_slots=lines/peers;tail_lines=lines-full_slots*peers
    nowned=full_slots+min(1,tail_lines)
    started=stamp()
    call prepare(p,n(axis),m(axis),nowned,peers,status)
    call comm_get_max(status,comm)
    if(status/=0)return
    fftw_block_seconds(1)=fftw_block_seconds(1)+stamp()-started
    started=stamp()
    count=nowned*m(axis)
    if(peers==1.and.axis==1)then
      ! The contiguous, unpartitioned axis can transform the field in place.
      plan=p%forward
      if(sign==1)plan=p%backward
!$omp parallel do default(none) private(ell,b,t,base,q) shared(lines,cross,n,plan,field)
      do ell=1,lines
        q=(ell-1)/cross
        b=1+q;t=(ell-1)-q*cross
        base=1+t*n(1)
        call fftw_execute_dft(plan,field(base:base+n(1)-1,b),field(base:base+n(1)-1,b))
      enddo
!$omp end parallel do
      fftw_block_seconds(2)=fftw_block_seconds(2)+stamp()-started
      return
    endif
    if(peers==1)then
      ! Batch adjacent lines in existing disjoint workspace, including a scalar tail.
      plan=p%forward;batch_plan=p%batch_forward
      if(sign==1)then
        plan=p%backward;batch_plan=p%batch_backward
      endif
      stride=product(m(:axis-1))
!$omp parallel do default(none) private(first,last,ell,b,t,q,point,g,base,i,batch_base) &
!$omp shared(lines,cross,other,m,n,axis,stride,p,plan,batch_plan,field)
      do first=1,lines,local_batch
        last=min(lines,first+local_batch-1)
        do ell=first,last
          q=(ell-1)/cross
          b=1+q;t=(ell-1)-q*cross
          q=t/m(other(1))
          point=0;point(other(2))=q;point(other(1))=t-q*m(other(1))
          g=1+point(1)+m(1)*(point(2)+m(2)*point(3))
          base=1+(ell-1)*n(axis)
          do i=0,n(axis)-1
            p%work(base+i)=field(g+i*stride,b)
          enddo
        enddo
        batch_base=1+(first-1)*n(axis)
        if(last-first+1==local_batch)then
          call fftw_execute_dft(batch_plan,p%work(batch_base:),p%work(batch_base:))
        else
          do ell=first,last
            base=1+(ell-1)*n(axis)
            call fftw_execute_dft(plan,p%work(base:),p%work(base:))
          enddo
        endif
        do ell=first,last
          q=(ell-1)/cross
          b=1+q;t=(ell-1)-q*cross
          q=t/m(other(1))
          point=0;point(other(2))=q;point(other(1))=t-q*m(other(1))
          g=1+point(1)+m(1)*(point(2)+m(2)*point(3))
          base=1+(ell-1)*n(axis)
          do i=0,n(axis)-1
            field(g+i*stride,b)=p%work(base+i)
          enddo
        enddo
      enddo
!$omp end parallel do
      fftw_block_fused_seconds=fftw_block_fused_seconds+stamp()-started
      return
    endif
    plan=p%forward
    if(sign==1)plan=p%backward
    ! Reuse one team across packing/FFT/restoration; only the master calls MPI.
!$omp parallel default(none) private(ell,b,t,point,owner,slot,i,g,j,s,base,q) &
!$omp shared(peers,lines,nowned,cross,other,m,n,axis,count,p,field,plan,comm,started,fftw_block_seconds) &
!$omp shared(full_slots,tail_lines)
      ! Rectangular full slots preserve OMP parallelism even with few peers.
!$omp do collapse(2) schedule(static)
      do owner=0,peers-1
        do slot=0,full_slots-1
          ell=1+owner+peers*slot
          q=(ell-1)/cross
          b=1+q;t=(ell-1)-q*cross
          q=t/m(other(1))
          point=0;point(other(2))=q;point(other(1))=t-q*m(other(1))
          do i=0,m(axis)-1
            point(axis)=i;g=1+point(1)+m(1)*(point(2)+m(2)*point(3))
            p%send(1+owner*count+slot*m(axis)+i)=field(g,b)
          enddo
        enddo
      enddo
!$omp end do
!$omp master
      ! The final partial slot has fewer than peers data lines.
      slot=full_slots
      do owner=0,tail_lines-1
        ell=1+owner+peers*slot
        q=(ell-1)/cross
        b=1+q;t=(ell-1)-q*cross
        q=t/m(other(1))
        point=0;point(other(2))=q;point(other(1))=t-q*m(other(1))
        do i=0,m(axis)-1
          point(axis)=i;g=1+point(1)+m(1)*(point(2)+m(2)*point(3))
          p%send(1+owner*count+slot*m(axis)+i)=field(g,b)
        enddo
      enddo
      do owner=lines-(nowned-1)*peers,peers-1
        p%send(1+owner*count+full_slots*m(axis):owner*count+nowned*m(axis))=0d0
      enddo
      fftw_block_seconds(4)=fftw_block_seconds(4)+stamp()-started
      started=stamp()
      call comm_alltoall(p%send,p%recv,comm,count)
      fftw_block_seconds(3)=fftw_block_seconds(3)+stamp()-started
      started=stamp()
!$omp end master
!$omp barrier
!$omp do collapse(2) schedule(static)
      do j=0,peers-1
        do s=0,nowned-1
          do i=0,m(axis)-1
            p%work(1+s*n(axis)+j*m(axis)+i)=p%recv(1+j*count+s*m(axis)+i)
          enddo
        enddo
      enddo
!$omp end do
!$omp master
      fftw_block_seconds(4)=fftw_block_seconds(4)+stamp()-started
!$omp end master
!$omp master
    started=stamp()
!$omp end master
    ! FFTW new-array execution is thread safe; plans are created outside OMP.
!$omp do schedule(static)
    do s=0,nowned-1
      base=1+s*n(axis)
      call fftw_execute_dft(plan,p%work(base:base+n(axis)-1),p%work(base:base+n(axis)-1))
    enddo
!$omp end do
!$omp master
    fftw_block_seconds(2)=fftw_block_seconds(2)+stamp()-started
    started=stamp()
!$omp end master
!$omp do collapse(2) schedule(static)
      do j=0,peers-1
        do s=0,nowned-1
          do i=0,m(axis)-1
            p%send(1+j*count+s*m(axis)+i)=p%work(1+s*n(axis)+j*m(axis)+i)
          enddo
        enddo
      enddo
!$omp end do
!$omp master
      fftw_block_seconds(4)=fftw_block_seconds(4)+stamp()-started
      started=stamp()
      call comm_alltoall(p%send,p%recv,comm,count)
      fftw_block_seconds(3)=fftw_block_seconds(3)+stamp()-started
      started=stamp()
!$omp end master
!$omp barrier
!$omp do collapse(2) schedule(static)
      do owner=0,peers-1
        do slot=0,full_slots-1
          ell=1+owner+peers*slot
          q=(ell-1)/cross
          b=1+q;t=(ell-1)-q*cross
          q=t/m(other(1))
          point=0;point(other(2))=q;point(other(1))=t-q*m(other(1))
          do i=0,m(axis)-1
            point(axis)=i;g=1+point(1)+m(1)*(point(2)+m(2)*point(3))
            field(g,b)=p%recv(1+owner*count+slot*m(axis)+i)
          enddo
        enddo
      enddo
!$omp end do
!$omp master
      slot=full_slots
      do owner=0,tail_lines-1
        ell=1+owner+peers*slot
        q=(ell-1)/cross
        b=1+q;t=(ell-1)-q*cross
        q=t/m(other(1))
        point=0;point(other(2))=q;point(other(1))=t-q*m(other(1))
        do i=0,m(axis)-1
          point(axis)=i;g=1+point(1)+m(1)*(point(2)+m(2)*point(3))
          field(g,b)=p%recv(1+owner*count+slot*m(axis)+i)
        enddo
      enddo
!$omp end master
!$omp master
    fftw_block_seconds(4)=fftw_block_seconds(4)+stamp()-started
!$omp end master
!$omp end parallel
  end subroutine
end module
