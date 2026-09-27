! Reusable FFTW transforms with native y/z pencil communicators.
! Local FFTW is serial within each MPI rank; channel batches bound workspace.
! Collective caller contract: identical grid, process dimensions and channel count;
! axis communicators ordered by coordinates. Saved caches require serial entry per rank.
module fftw_pencils
  use iso_c_binding
  use iso_fortran_env, only: int64
  use communication, only: comm_alltoall,comm_get_max
  implicit none
  private
  include 'fftw3.f03'
  integer,parameter :: max_batch=4
  type redistribution
    integer :: peers=0,axis=0
    integer,allocatable :: pack(:),unpack(:)
  end type
  type pencil_cache
    integer :: n(3)=0,dims(2)=0,coords(2)=0,batch=0,nt=0
    complex(c_double_complex),allocatable :: work(:)
    complex(c_double_complex),allocatable :: send(:),recv(:)
    type(c_ptr) :: plans(3,2)=c_null_ptr
    type(redistribution) :: steps(4)
  end type
  type(pencil_cache),save :: cache(max_batch)
  integer,save,public :: fftw_pencil_plans_created=0
  real(8),save,public :: fftw_pencil_seconds(4)=0d0 ! setup, FFT, MPI, packing/copies
  public :: pencil_transform,pencil_clear
contains
  real(8) function stamp()
    integer(int64) :: count,rate
    call system_clock(count,rate)
    stamp=real(count,8)/real(rate,8)
  end function
  subroutine destroy(p)
    type(pencil_cache),intent(inout) :: p
    integer :: a,b
    do b=1,2;do a=1,3
      if(c_associated(p%plans(a,b)))call fftw_destroy_plan(p%plans(a,b))
    enddo;enddo
    p%plans=c_null_ptr
    if(allocated(p%work))deallocate(p%work)
    if(allocated(p%send))deallocate(p%send,p%recv)
    do a=1,4
      if(allocated(p%steps(a)%pack))deallocate(p%steps(a)%pack,p%steps(a)%unpack)
    enddo
    p%n=0;p%batch=0;p%nt=0
  end subroutine
  subroutine pencil_clear()
    integer :: a
    do a=1,max_batch
      call destroy(cache(a))
    enddo
  end subroutine
  subroutine layout(n,dims,coords,which,shape,origin,order)
    integer,intent(in) :: n(3),dims(2),coords(2),which
    integer,intent(out) :: shape(3),origin(3),order(3)
    select case(which)
    case(1)
      shape=[n(1),n(2)/dims(1),n(3)/dims(2)]
      origin=[0,coords(1)*shape(2),coords(2)*shape(3)];order=[1,2,3]
    case(2)
      shape=[n(1)/dims(1),n(2),n(3)/dims(2)]
      origin=[coords(1)*shape(1),0,coords(2)*shape(3)];order=[2,1,3]
    case(3)
      shape=[n(1)/dims(1),n(2)/dims(2),n(3)]
      origin=[coords(1)*shape(1),coords(2)*shape(2),0];order=[3,1,2]
    end select
  end subroutine
  subroutine mapping(p,step,from,to,axis,split_to,split_from)
    type(pencil_cache),intent(inout) :: p
    integer,intent(in) :: step,from,to,axis,split_to,split_from
    integer :: shape(3),origin(3),order(3),point(3),global(3),counter(p%dims(axis))
    integer :: x,y,z,index,peer,block,pass,which,split
    p%steps(step)%peers=p%dims(axis);p%steps(step)%axis=axis
    allocate(p%steps(step)%pack(p%nt),p%steps(step)%unpack(p%nt))
    block=p%nt/p%dims(axis)
    do pass=1,2
      which=from;split=split_to
      if(pass==2)then
        which=to;split=split_from
      endif
      call layout(p%n,p%dims,p%coords,which,shape,origin,order)
      counter=0
      ! Identical physical xyz ordering on both sides of each peer intersection.
      do z=0,shape(3)-1;do y=0,shape(2)-1;do x=0,shape(1)-1
        point=[x,y,z];global=point+origin
        peer=global(split)/(p%n(split)/p%dims(axis))
        counter(peer+1)=counter(peer+1)+1
        index=1+point(order(1))+shape(order(1))*(point(order(2))+shape(order(2))*point(order(3)))
        if(pass==1)then
          p%steps(step)%pack(index)=peer*block*p%batch+counter(peer+1)
        else
          p%steps(step)%unpack(index)=peer*block*p%batch+counter(peer+1)
        endif
      enddo;enddo;enddo
    enddo
  end subroutine
  subroutine prepare(p,n,dims,coords,batch,status)
    type(pencil_cache),intent(inout) :: p
    integer,intent(in) :: n(3),dims(2),coords(2),batch
    integer,intent(out) :: status
    integer :: a,b,sgn
    real(8) :: start
    status=0
    if(all(p%n==n).and.all(p%dims==dims).and.all(p%coords==coords).and.p%batch==batch)return
    start=stamp();call destroy(p)
    p%n=n;p%dims=dims;p%coords=coords;p%batch=batch;p%nt=n(1)*(n(2)/dims(1))*(n(3)/dims(2))
    allocate(p%work(p%nt*batch),p%send(p%nt*batch),p%recv(p%nt*batch))
    p%work=0d0
    call mapping(p,1,1,2,1,1,2)
    call mapping(p,2,2,3,2,2,3)
    call mapping(p,3,3,2,2,3,2)
    call mapping(p,4,2,1,1,2,1)
    do b=1,2
      sgn=FFTW_FORWARD
      if(b==2)sgn=FFTW_BACKWARD
      do a=1,3
        p%plans(a,b)=fftw_plan_many_dft(1,[n(a)],p%nt*batch/n(a),p%work,[n(a)],1,n(a), &
          p%work,[n(a)],1,n(a),sgn,FFTW_MEASURE)
        if(.not.c_associated(p%plans(a,b)))status=1
        fftw_pencil_plans_created=fftw_pencil_plans_created+1
      enddo
    enddo
    fftw_pencil_seconds(1)=fftw_pencil_seconds(1)+stamp()-start
    if(status/=0)call destroy(p)
  end subroutine
  subroutine redistribute(p,step,comm)
    type(pencil_cache),intent(inout) :: p
    integer,intent(in) :: step,comm(2)
    integer :: i,q,block
    real(8) :: start
    block=p%nt/p%steps(step)%peers;start=stamp()
    do q=1,p%batch;do i=1,p%nt
      p%send(p%steps(step)%pack(i)+(q-1)*block)=p%work(i+(q-1)*p%nt)
    enddo;enddo
    fftw_pencil_seconds(4)=fftw_pencil_seconds(4)+stamp()-start;start=stamp()
    call comm_alltoall(p%send,p%recv,comm(p%steps(step)%axis),block*p%batch)
    fftw_pencil_seconds(3)=fftw_pencil_seconds(3)+stamp()-start;start=stamp()
    do q=1,p%batch;do i=1,p%nt
      p%work(i+(q-1)*p%nt)=p%recv(p%steps(step)%unpack(i)+(q-1)*block)
    enddo;enddo
    fftw_pencil_seconds(4)=fftw_pencil_seconds(4)+stamp()-start
  end subroutine
  subroutine pencil_transform(n,dims,coords,comm,input,output,sgn,status)
    integer,intent(in) :: n(3),dims(2),coords(2),comm(2),sgn
    complex(8),intent(in) :: input(:,:)
    complex(8),intent(out) :: output(:,:)
    integer,intent(out) :: status
    integer :: first,batch,a,direction,nt,bad
    integer(int64) :: wide_nt
    real(8) :: start
    status=0
    if(any(n<1).or.any(dims<1))status=1
    if(status==0)then
      if(modulo(n(1),dims(1))/=0.or.modulo(n(2),dims(1))/=0.or. &
         modulo(n(2),dims(2))/=0.or.modulo(n(3),dims(2))/=0)status=1
      if(any(coords<0).or.any(coords>=dims).or.abs(sgn)/=1)status=1
      wide_nt=int(n(1),int64)*(n(2)/dims(1))*(n(3)/dims(2))
      if(wide_nt>huge(nt)/max_batch)status=1
      nt=int(min(wide_nt,int(huge(nt)/max_batch,int64)))
      if(size(input,1)/=nt.or.any(shape(input)/=shape(output)).or.size(input,2)<1)status=1
    endif
    start=stamp()
    call comm_get_max(status,comm(1));call comm_get_max(status,comm(2))
    fftw_pencil_seconds(3)=fftw_pencil_seconds(3)+stamp()-start
    if(status/=0)return
    direction=1
    if(sgn==1)direction=2
    batch=min(max_batch,size(input,2))
    call prepare(cache(batch),n,dims,coords,batch,status)
    batch=modulo(size(input,2),max_batch)
    if(batch>0.and.size(input,2)>max_batch)then
      call prepare(cache(batch),n,dims,coords,batch,bad)
      status=max(status,bad)
    endif
    start=stamp()
    call comm_get_max(status,comm(1));call comm_get_max(status,comm(2))
    fftw_pencil_seconds(3)=fftw_pencil_seconds(3)+stamp()-start
    if(status/=0)return
    do first=1,size(input,2),max_batch
      batch=min(max_batch,size(input,2)-first+1)
      start=stamp();cache(batch)%work=reshape(input(:,first:first+batch-1),[nt*batch])
      fftw_pencil_seconds(4)=fftw_pencil_seconds(4)+stamp()-start
      do a=1,3
        start=stamp()
        call fftw_execute_dft(cache(batch)%plans(a,direction),cache(batch)%work,cache(batch)%work)
        fftw_pencil_seconds(2)=fftw_pencil_seconds(2)+stamp()-start
        if(a<3)call redistribute(cache(batch),a,comm)
      enddo
      call redistribute(cache(batch),3,comm);call redistribute(cache(batch),4,comm)
      start=stamp()
      if(sgn==1)cache(batch)%work=cache(batch)%work/product(real(n,8))
      output(:,first:first+batch-1)=reshape(cache(batch)%work,[nt,batch])
      fftw_pencil_seconds(4)=fftw_pencil_seconds(4)+stamp()-start
    enddo
  end subroutine
end module
