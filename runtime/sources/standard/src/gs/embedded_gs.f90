! Retained standard GS C handles. Host owns MPI; all calls are subgroup collective.
module salmonpy_embedded_gs
use iso_c_binding
use salmonpy_gs_storage, only: s_gs_storage,gs_storage_clear
use salmonpy_gs_model, only: s_gs_model
use communication, only: comm_is_attached,comm_get_globalinfo,comm_bcast,comm_get_max,comm_child_counts
implicit none
private
public :: borrow_gs
type :: gs_owner
  type(s_gs_model) :: model
  type(s_gs_storage) :: state
end type
type :: gs_slot
  type(gs_owner),pointer :: owner=>null()
  integer(c_int64_t) :: token=0
end type
type(gs_slot),allocatable,save :: slots(:)
integer(c_int64_t),save :: serial=0
integer(c_int64_t),parameter :: prefix=int(z'4700000000000000',c_int64_t)
logical,save :: initialized=.false.
contains
integer(c_int) function sp_gs_protocol() bind(C)
  sp_gs_protocol=1
end function

subroutine lookup(token,index,status)
  integer(c_int64_t),intent(in) :: token
  integer,intent(out) :: index,status
  integer :: group,rank,nproc,i,parts(2),reference(2)
  index=0;status=10
  if(.not.comm_is_attached())return
  call comm_get_globalinfo(group,rank,nproc)
  parts=transfer(token,parts)
  reference=parts
  call comm_bcast(reference,group)
  if(all(reference==parts).and.token>prefix.and.allocated(slots))then
    do i=1,size(slots)
      if(slots(i)%token/=token.or..not.associated(slots(i)%owner))cycle
      index=i;status=0;exit
    enddo
  endif
  call comm_get_max(status,group)
  if(status/=0)index=0
end subroutine

subroutine borrow_gs(token,gs,model,status)
  integer(c_int64_t),intent(in) :: token
  type(s_gs_storage),pointer,intent(out) :: gs
  type(s_gs_model),pointer,intent(out) :: model
  integer,intent(out) :: status
  integer :: index
  nullify(gs,model)
  call lookup(token,index,status)
  if(status/=0)return
  gs=>slots(index)%owner%state
  model=>slots(index)%owner%model
end subroutine

integer(c_int) function sp_gs_create(filename,token) bind(C)
  character(c_char),intent(in) :: filename(*)
  integer(c_int64_t),intent(out) :: token
  sp_gs_create=create_gs(filename,token,.false.)
end function

integer(c_int) function sp_gs_import(filename,token) bind(C)
  character(c_char),intent(in) :: filename(*)
  integer(c_int64_t),intent(out) :: token
  sp_gs_import=create_gs(filename,token,.true.)
end function

integer(c_int) function create_gs(filename,token,import_only)
  use math_constants, only: set_math_constants
  use parallelization, only: setup_parallel
  use timer, only: timer_initialize
  character(c_char),intent(in) :: filename(*)
  integer(c_int64_t),intent(out) :: token
  type(gs_slot),allocatable :: enlarged(:)
  character(4096) :: path,expected
  integer :: group,rank,nproc,i,index,status,n
  interface
    integer(c_int) function sp_tddft_park_for_gs() bind(C)
      import c_int
    end function
    subroutine sp_tddft_resume_after_gs() bind(C)
    end subroutine
  end interface
  logical,intent(in) :: import_only
  token=0;create_gs=10
  if(.not.comm_is_attached())return
  call comm_get_globalinfo(group,rank,nproc)
  path='';status=0
  do n=1,len(path)
    if(filename(n)==c_null_char)exit
    path(n:n)=filename(n)
  enddo
  if(n>len(path).or.len_trim(path)==0)status=10
  expected=path
  call comm_bcast(expected,group)
  if(expected/=path)status=10
  if(serial>=huge(serial)-prefix)status=14
  call comm_get_max(status,group)
  if(status/=0)then
    create_gs=status;return
  endif
  status=sp_tddft_park_for_gs()
  call comm_get_max(status,group)
  if(status/=0)then
    create_gs=status;return
  endif
  if(.not.initialized)then
    call set_math_constants
    call setup_parallel
    call timer_initialize
    initialized=.true.
  endif
  if(.not.allocated(slots))allocate(slots(8))
  index=0
  do i=1,size(slots)
    if(associated(slots(i)%owner))cycle
    index=i;exit
  enddo
  if(index==0)then
    index=size(slots)+1
    allocate(enlarged(size(slots)*2))
    enlarged(:size(slots))=slots
    call move_alloc(enlarged,slots)
  endif
  allocate(slots(index)%owner)
  call slots(index)%owner%model%load(trim(path),status)
  call comm_get_max(status,group)
  if(status==0)call slots(index)%owner%model%run(status,slots(index)%owner%state,import_only)
  call comm_get_max(status,group)
  call sp_tddft_resume_after_gs()
  if(status/=0)then
    call gs_storage_clear(slots(index)%owner%state)
    call slots(index)%owner%model%clear()
    deallocate(slots(index)%owner)
    create_gs=status;return
  endif
  serial=serial+1
  token=prefix+serial
  slots(index)%token=token
  create_gs=0
end function

integer(c_int) function sp_gs_info(token,stats,values) bind(C)
  integer(c_int64_t),value :: token
  integer(c_int),intent(out) :: stats(4)
  real(c_double),intent(out) :: values(3)
  integer :: index,status
  stats=0;values=0d0
  call lookup(token,index,status)
  sp_gs_info=status
  if(status/=0)return
  associate(s=>slots(index)%owner%state)
    stats=[s%iterations,merge(1,0,s%converged),s%system%nion,s%system%no]
    if(s%imported)stats(2)=-1
    values=[s%energy%E_tot,s%residual,s%threshold]
  end associate
end function

subroutine array_layout(s,quantity,ndim,shape_out,lower,status)
  type(s_gs_storage),intent(in) :: s
  integer,intent(in) :: quantity
  integer(c_int),intent(out) :: ndim
  integer(c_int64_t),intent(out) :: shape_out(7),lower(7)
  integer,intent(out) :: status
  ndim=0;shape_out=1;lower=1;status=0
  select case(quantity)
  case(1)
    ndim=3;shape_out(:3)=shape(s%rho%f,kind=c_int64_t);lower(:3)=lbound(s%rho%f,kind=c_int64_t)
  case(2)
    ndim=3;shape_out(:3)=shape(s%system%rocc,kind=c_int64_t);lower(:3)=lbound(s%system%rocc,kind=c_int64_t)
  case(3)
    ndim=3;shape_out(:3)=shape(s%energy%esp,kind=c_int64_t);lower(:3)=lbound(s%energy%esp,kind=c_int64_t)
  case(4)
    if(s%imported)then
      status=11;return
    endif
    ndim=2;shape_out(:2)=shape(s%system%Force,kind=c_int64_t)
  case(5)
    ndim=7
    if(allocated(s%spsi%rwf))then
      shape_out=shape(s%spsi%rwf,kind=c_int64_t);lower=lbound(s%spsi%rwf,kind=c_int64_t)
    else
      shape_out=shape(s%spsi%zwf,kind=c_int64_t);lower=lbound(s%spsi%zwf,kind=c_int64_t)
    endif
  case default
    status=11
  end select
end subroutine

integer(c_int) function sp_gs_array_info(token,quantity,ndim,shape_out,lower) bind(C)
  integer(c_int64_t),value :: token
  integer(c_int),value :: quantity
  integer(c_int),intent(out) :: ndim
  integer(c_int64_t),intent(out) :: shape_out(7),lower(7)
  integer :: index,status
  call lookup(token,index,status)
  if(status==0)call array_layout(slots(index)%owner%state,quantity,ndim,shape_out,lower,status)
  sp_gs_array_info=status
end function

integer(c_int) function sp_gs_read(token,quantity,real_ptr,imag_ptr,length) bind(C)
  integer(c_int64_t),value :: token,length
  integer(c_int),value :: quantity
  type(c_ptr),value :: real_ptr,imag_ptr
  real(c_double),pointer :: re(:),im(:)
  integer(c_int64_t) :: shape_out(7),lower(7),k
  integer(c_int) :: ndim
  integer :: index,status,i1,i2,i3,i4,i5,i6,i7,group,rank,nproc
  call lookup(token,index,status)
  if(status/=0)then
    sp_gs_read=status;return
  endif
  call array_layout(slots(index)%owner%state,quantity,ndim,shape_out,lower,status)
  if(status==0)then
    if(length/=product(shape_out).or..not.c_associated(real_ptr))status=11
    if(quantity==5.and..not.c_associated(imag_ptr))status=11
  endif
  call comm_get_globalinfo(group,rank,nproc)
  call comm_get_max(status,group)
  sp_gs_read=status
  if(status/=0)return
  call c_f_pointer(real_ptr,re,[length])
  associate(s=>slots(index)%owner%state)
    select case(quantity)
    case(1)
      re=reshape(s%rho%f,[length])
    case(2)
      re=reshape(s%system%rocc,[length])
    case(3)
      re=reshape(s%energy%esp,[length])
    case(4)
      re=reshape(s%system%Force,[length])
    case(5)
      call c_f_pointer(imag_ptr,im,[length])
      k=0
      do i7=int(lower(7)),int(lower(7)+shape_out(7)-1)
      do i6=int(lower(6)),int(lower(6)+shape_out(6)-1)
      do i5=int(lower(5)),int(lower(5)+shape_out(5)-1)
      do i4=int(lower(4)),int(lower(4)+shape_out(4)-1)
      do i3=int(lower(3)),int(lower(3)+shape_out(3)-1)
      do i2=int(lower(2)),int(lower(2)+shape_out(2)-1)
      do i1=int(lower(1)),int(lower(1)+shape_out(1)-1)
        k=k+1
        if(allocated(s%spsi%rwf))then
          re(k)=s%spsi%rwf(i1,i2,i3,i4,i5,i6,i7);im(k)=0d0
        else
          re(k)=real(s%spsi%zwf(i1,i2,i3,i4,i5,i6,i7),c_double)
          im(k)=aimag(s%spsi%zwf(i1,i2,i3,i4,i5,i6,i7))
        endif
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
    end select
  end associate
end function

integer(c_int) function sp_gs_close(token) bind(C)
  integer(c_int64_t),value :: token
  integer :: index,status
  call lookup(token,index,status)
  sp_gs_close=status
  if(status/=0)return
  call gs_storage_clear(slots(index)%owner%state)
  call slots(index)%owner%model%clear()
  deallocate(slots(index)%owner)
  slots(index)%token=0
end function

integer(c_int) function sp_gs_resources(counts) bind(C)
  integer(c_int),intent(out) :: counts(3)
  integer :: i,live,nslots
  counts=0
  if(allocated(slots))then
    do i=1,size(slots)
      if(associated(slots(i)%owner))counts(1)=counts(1)+1
    enddo
  endif
  call comm_child_counts(live,nslots)
  counts(2:3)=[live,nslots]
  sp_gs_resources=0
end function
end module
