! Preview ABI1 declarations for a serial endpoint-A material plugin.
module salmonpy_material_iso
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: salmonpy_material_api_v1,salmonpy_material_query
  public :: sp_model_create,sp_destroy,sp_state_create,sp_state_observe,sp_state_step
  public :: SALMONPY_MATERIAL_ABI,SALMONPY_OK,SALMONPY_ARGUMENT,SALMONPY_HANDLE
  public :: SALMONPY_MEMORY,SALMONPY_NUMERICAL,SALMONPY_BUSY,SALMONPY_ABI,SALMONPY_INTERNAL
  integer(c_int32_t),parameter :: SALMONPY_MATERIAL_ABI=1,SALMONPY_OK=0
  integer(c_int32_t),parameter :: SALMONPY_ARGUMENT=1,SALMONPY_HANDLE=2,SALMONPY_MEMORY=3
  integer(c_int32_t),parameter :: SALMONPY_NUMERICAL=4,SALMONPY_BUSY=5,SALMONPY_ABI=6,SALMONPY_INTERNAL=7
  type, bind(C) :: salmonpy_material_api_v1
    integer(c_int32_t) :: abi_version,struct_size
    type(c_ptr) :: backend
    integer(c_int32_t) :: input_lag,uses_initial_A
    type(c_funptr) :: model_create,model_destroy,state_create,state_observe,state_step,state_destroy
  end type
  interface
    function salmonpy_material_query(version,capacity,api) bind(C,name='salmonpy_material_query') result(status)
      import
      integer(c_int32_t),value :: version
      integer(c_size_t),value :: capacity
      type(salmonpy_material_api_v1),intent(out) :: api
      integer(c_int32_t) :: status
    end function
  end interface
  abstract interface
    function sp_model_create(config,size,handle,error,capacity) bind(C) result(status)
      import
      type(c_ptr),value :: config
      integer(c_size_t),value :: size,capacity
      integer(c_int64_t),intent(out) :: handle
      character(kind=c_char),intent(out) :: error(*)
      integer(c_int32_t) :: status
    end function
    function sp_destroy(handle,error,capacity) bind(C) result(status)
      import
      integer(c_int64_t),value :: handle
      character(kind=c_char),intent(out) :: error(*)
      integer(c_size_t),value :: capacity
      integer(c_int32_t) :: status
    end function
    function sp_state_create(model,dt,t0,A,handle,error,capacity) bind(C) result(status)
      import
      integer(c_int64_t),value :: model
      real(c_double),value :: dt,t0
      real(c_double),intent(in) :: A(3)
      integer(c_int64_t),intent(out) :: handle
      character(kind=c_char),intent(out) :: error(*)
      integer(c_size_t),value :: capacity
      integer(c_int32_t) :: status
    end function
    function sp_state_observe(handle,time,current,error,capacity) bind(C) result(status)
      import
      integer(c_int64_t),value :: handle
      real(c_double),intent(out) :: time,current(3)
      character(kind=c_char),intent(out) :: error(*)
      integer(c_size_t),value :: capacity
      integer(c_int32_t) :: status
    end function
    function sp_state_step(handle,A,time,current,error,capacity) bind(C) result(status)
      import
      integer(c_int64_t),value :: handle
      real(c_double),intent(in) :: A(3)
      real(c_double),intent(out) :: time,current(3)
      character(kind=c_char),intent(out) :: error(*)
      integer(c_size_t),value :: capacity
      integer(c_int32_t) :: status
    end function
  end interface
end module
