! Demand-sized payload IO for validated complex DC files (interleaved real/imaginary).
! Storage follows the requested domain, coefficient column and contiguous runs.
module lcfo_mesh_stream
  use iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: lcfo_stream_contract
contains
  subroutine read_values(unit,position,values,status)
    integer,intent(in) :: unit
    integer(int64),intent(in) :: position
    complex(8),intent(out) :: values(:)
    integer,intent(out) :: status
    real(8),allocatable :: wire(:)
    integer :: i,n
    status=1;n=size(values)
    if(position<1_int64)return
    allocate(wire(2*n))
    read(unit,pos=position,iostat=status)wire(:2*n)
    if(status/=0)return
    if(.not.all(ieee_is_finite(wire(:2*n))))then
      status=1;return
    endif
    do i=1,n
      values(i)=cmplx(wire(2*i-1),wire(2*i),8)
    enddo
  end subroutine

  subroutine lcfo_stream_contract(ub,uc,basis_pos,coef_pos,core,nb,io,jxyz,lo,m,first,tile,status)
    integer,intent(in) :: ub,uc,core(3),nb,io,jxyz(:,:),lo(3),m(3),first
    integer(int64),intent(in) :: basis_pos,coef_pos
    complex(8),intent(inout) :: tile(:)
    integer,intent(out) :: status
    integer(int64),allocatable :: source(:)
    integer,allocatable :: target(:)
    complex(8),allocatable :: values(:),coeff(:)
    integer(int64) :: ng,volume,index,offset,limit
    integer :: x,y,z,g,np,b,p,count,j,max_run
    status=1
    if(any(core<1).or.any(m<1).or.nb<0.or.io<1.or.first<1)return
    if(size(jxyz,1)<maxval(core).or.size(jxyz,2)<3)return
    ng=int(core(1),int64)*core(2)*core(3)
    volume=int(m(1),int64)*m(2)*m(3)
    if(int(first,int64)+size(tile)-1>volume)return
    if(basis_pos<1.or.coef_pos<1)return
    limit=(huge(ng)-max(basis_pos,coef_pos))/16_int64
    if(ng>limit/max(1,nb).or.int(io,int64)>limit/max(1,nb))return
    if(nb==0)then
      status=0;return
    endif
    allocate(source(size(tile)),target(size(tile)));np=0
    ! File-order traversal produces ordered runs even with wrapped global maps.
    do z=1,core(3)
      if(jxyz(z,3)<lo(3).or.jxyz(z,3)>=lo(3)+m(3))cycle
      do y=1,core(2)
        if(jxyz(y,2)<lo(2).or.jxyz(y,2)>=lo(2)+m(2))cycle
        do x=1,core(1)
          if(jxyz(x,1)<lo(1).or.jxyz(x,1)>=lo(1)+m(1))cycle
          index=1_int64+jxyz(x,1)-lo(1)+int(m(1),int64)* &
            (jxyz(y,2)-lo(2)+int(m(2),int64)*(jxyz(z,3)-lo(3)))
          index=index-first+1
          if(index<1.or.index>size(tile))cycle
          np=np+1
          if(np>size(tile))return ! preflight excludes overlapping core points
          target(np)=int(index)
          source(np)=x+int(core(1),int64)*(y-1+int(core(2),int64)*(z-1))
        enddo
      enddo
    enddo
    if(np==0)then
      status=0;return
    endif
    ! A single requested orbital's coefficient column, never the full matrix.
    allocate(coeff(nb))
    offset=coef_pos+16_int64*int(io-1,int64)*nb
    call read_values(uc,offset,coeff,status)
    if(status/=0)return
    max_run=0;p=1
    do while(p<=np)
      count=1
      do while(p+count<=np)
        if(source(p+count)/=source(p)+count)exit
        count=count+1
      enddo
      max_run=max(max_run,count);p=p+count
    enddo
    allocate(values(max_run))
    do b=1,nb
      p=1
      do while(p<=np)
        count=1
        do while(p+count<=np)
          if(source(p+count)/=source(p)+count)exit
          count=count+1
        enddo
        offset=basis_pos+16_int64*(ng*(b-1)+source(p)-1)
        call read_values(ub,offset,values(:count),status)
        if(status/=0)return
        do j=1,count
          g=target(p+j-1)
          tile(g)=tile(g)+values(j)*coeff(b)
        enddo
        p=p+count
      enddo
    enddo
    status=0
  end subroutine
end module
