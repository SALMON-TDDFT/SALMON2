! Cached payload IO for validated complex DC files (interleaved real/imaginary).
! Fragment payloads are read once and contracted without repeated file access.
module lcfo_mesh_stream
  use iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: s_lcfo_cache, lcfo_cache_load, lcfo_cache_contract
  type :: s_lcfo_cache
    complex(8),allocatable :: basis(:,:),coeff(:,:)
  end type
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_1d
  interface salmon_all_finite
    module procedure finite_real_1d
  end interface
contains
  subroutine lcfo_cache_load(ub,uc,basis_pos,coef_pos,core,nb,norb,cache,status)
    implicit none
    integer,intent(in)::ub,uc,core(3),nb,norb
    integer(int64),intent(in)::basis_pos,coef_pos
    type(s_lcfo_cache),intent(inout)::cache
    integer,intent(out)::status
    integer::b,io,ng,ios
    status=1
    if(any(core<1).or.nb<1.or.norb<1)return
    if(allocated(cache%basis))deallocate(cache%basis)
    if(allocated(cache%coeff))deallocate(cache%coeff)
    ng=product(core)
    allocate(cache%basis(ng,nb),cache%coeff(nb,norb),stat=ios)
    if(ios/=0)return
    do b=1,nb
      call read_values(ub,basis_pos+16_int64*int(ng,int64)*(b-1),cache%basis(:,b),status)
      if(status/=0)return
    enddo
    do io=1,norb
      call read_values(uc,coef_pos+16_int64*int(nb,int64)*(io-1),cache%coeff(:,io),status)
      if(status/=0)return
    enddo
    status=0
  end subroutine

  subroutine lcfo_cache_contract(cache,core,io,jxyz,lo,m,first,tile,status)
    implicit none
    type(s_lcfo_cache),intent(in)::cache
    integer,intent(in)::core(3),io,jxyz(:,:),lo(3),m(3),first
    complex(8),intent(inout)::tile(:)
    integer,intent(out)::status
    integer::b,x,y,z,g,source
    integer(int64)::index
    status=1
    if(.not.allocated(cache%basis).or..not.allocated(cache%coeff))return
    if(io<1.or.io>size(cache%coeff,2).or.any(core<1).or.any(m<1).or.first<1)return
    if(size(jxyz,1)<maxval(core).or.size(jxyz,2)<3)return
    if(size(cache%basis,1)/=product(core).or.size(cache%basis,2)/=size(cache%coeff,1))return
    if(int(first,int64)+size(tile)-1>int(m(1),int64)*m(2)*m(3))return
    ! Keep the scalar path's basis accumulation order; reuse payloads, not approximations.
    do b=1,size(cache%coeff,1)
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
            g=int(index);source=x+core(1)*(y-1+core(2)*(z-1))
            tile(g)=tile(g)+cache%basis(source,b)*cache%coeff(b,io)
          enddo
        enddo
      enddo
    enddo
    status=0
  end subroutine

  subroutine read_values(unit,position,values,status)
    implicit none
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
    if(.not.salmon_all_finite(wire(:2*n)))then
      status=1;return
    endif
    do i=1,n
      values(i)=cmplx(wire(2*i-1),wire(2*i),8)
    enddo
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

end module
