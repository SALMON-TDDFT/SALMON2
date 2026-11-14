program host
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  use salmonpy_material_iso
  implicit none
  type(salmonpy_material_api_v1) :: api
  procedure(sp_model_create),pointer :: model_create
  procedure(sp_destroy),pointer :: model_destroy,state_destroy
  procedure(sp_state_create),pointer :: state_create
  procedure(sp_state_observe),pointer :: observe
  procedure(sp_state_step),pointer :: step
  real(c_double),target :: config(2)=[2d0,.8d0]
  real(c_double) :: A(3),zero(3)=0d0,current(3),reference(3)=0d0,last(3)=0d0,time,bad(3)
  integer(c_int64_t) :: model,one,two
  integer(c_int32_t) :: status
  character(kind=c_char) :: error(256)
  integer :: n,k
  call require(salmonpy_material_query(SALMONPY_MATERIAL_ABI,c_sizeof(api),api)==SALMONPY_OK)
  call require(api%abi_version==SALMONPY_MATERIAL_ABI .and. api%struct_size==c_sizeof(api))
  call c_f_procpointer(api%model_create,model_create)
  call c_f_procpointer(api%model_destroy,model_destroy)
  call c_f_procpointer(api%state_create,state_create)
  call c_f_procpointer(api%state_observe,observe)
  call c_f_procpointer(api%state_step,step)
  call c_f_procpointer(api%state_destroy,state_destroy)
  call require(model_create(c_loc(config),c_sizeof(config),model,error,size(error,kind=c_size_t))==SALMONPY_OK)
  config(1)=999d0
  call require(state_create(model,.02d0,0d0,zero,one,error,size(error,kind=c_size_t))==SALMONPY_OK)
  call require(state_create(model,.02d0,0d0,zero,two,error,size(error,kind=c_size_t))==SALMONPY_OK)
  call require(one/=two)
  call require(model_destroy(model,error,size(error,kind=c_size_t))==SALMONPY_BUSY)
  do n=1,8
    A=[n*.01d0,n*.02d0,-n*.005d0]
    call require(step(one,A,time,current,error,size(error,kind=c_size_t))==SALMONPY_OK)
    do k=1,3
      reference(k)=.8d0*reference(k)+2d0*(A(k)-last(k));last(k)=A(k)
      call require(ieee_is_finite(current(k)))
      call require(abs(current(k)-reference(k))<1d-14)
    end do
    call require(time==n*.02d0)
    call require(observe(two,time,current,error,size(error,kind=c_size_t))==SALMONPY_OK)
    call require(time==0d0)
    do k=1,3
      call require(current(k)==0d0)
    end do
  end do
  bad=0d0;bad(1)=ieee_value(0d0,ieee_quiet_nan)
  call require(step(one,bad,time,current,error,size(error,kind=c_size_t))==SALMONPY_NUMERICAL)
  call require(observe(one,time,current,error,size(error,kind=c_size_t))==SALMONPY_OK)
  call require(time==8*.02d0)
  do k=1,3
    call require(current(k)==reference(k))
  end do
  call require(state_destroy(one,error,size(error,kind=c_size_t))==SALMONPY_OK)
  call require(observe(one,time,current,error,size(error,kind=c_size_t))==SALMONPY_HANDLE)
  call require(state_destroy(two,error,size(error,kind=c_size_t))==SALMONPY_OK)
  call require(model_destroy(model,error,size(error,kind=c_size_t))==SALMONPY_OK)
  print *,'Fortran material ABI, independent handles, ownership and atomic failure PASS'
contains
  subroutine require(ok)
    logical,intent(in) :: ok
    if(.not.ok) error stop 'Fortran ABI test failed'
  end subroutine
end program
