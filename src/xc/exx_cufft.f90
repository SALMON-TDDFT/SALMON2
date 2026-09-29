#include "config.h"
! Optional resident-batch backend for the existing compact discrete kernel.
! Assign ACC_DEVICE_NUM (or CUDA_VISIBLE_DEVICES) per MPI rank. This owning
! object must not be copied: use the factory, prepare/apply, then release.
module exx_cufft
  use exx_batch_backend, only: s_exx_batch_backend
  use iso_c_binding, only: c_intptr_t
  implicit none
  private
  public :: s_exx_cufft,exx_cufft_create,exx_cufft_apply
  type,extends(s_exx_batch_backend) :: s_exx_cufft
    private
    integer :: padded(3)=0,ns=0,npoints=0,capacity=0,device=-1,plan=0
    integer(c_intptr_t) :: stream=0
    logical :: ready=.false.,resident=.false.,plan_created=.false.
    integer,allocatable :: indices(:)
    complex(8),allocatable :: filter(:,:,:),source(:),work(:),targets(:,:),output(:,:)
    ! Actual operations, accumulated over the lifetime of this backend object.
    integer,public :: plan_builds=0,filter_uploads=0,index_uploads=0,source_uploads=0,batch_uploads=0
  contains
    procedure :: prepare => cufft_prepare
    procedure :: apply => cufft_apply_batch
    procedure :: release => cufft_release
    final :: cufft_finalize
  end type
#ifdef USE_EXX_CUFFT
  interface
    subroutine set_acc_error_routine(callback) bind(C,name='acc_set_error_routine')
      use iso_c_binding, only: c_funptr
      implicit none
      type(c_funptr),value,intent(in) :: callback
    end subroutine
  end interface
#endif
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_1d,finite_real_2d,finite_real_3d
  interface salmon_all_finite
    module procedure finite_real_1d,finite_real_2d,finite_real_3d
  end interface
contains
  subroutine exx_cufft_create(backend)
    implicit none
    class(s_exx_batch_backend),allocatable,intent(out) :: backend
    allocate(s_exx_cufft::backend)
  end subroutine

  subroutine validate_source(padded,indices,filter,source,capacity,npoints,status)
    use iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,intent(in) :: padded(3),indices(:),capacity
    complex(8),intent(in) :: filter(:,:,:),source(:)
    integer,intent(out) :: npoints,status
    integer :: ns,axis,i,ios
    integer(int64) :: points
    logical,allocatable :: seen(:)
    status=-2;npoints=0
    if(any(padded<1).or.capacity<1)return
    if(any(shape(filter)/=padded))return
    if(size(source,kind=int64)>int(huge(0),int64))return
    ns=size(source)
    if(size(indices)/=ns)return
    points=1_int64
    do axis=1,3
      if(points>int(huge(0),int64)/int(padded(axis),int64))return
      points=points*int(padded(axis),int64)
    enddo
    if(points>int(huge(0),int64)/int(capacity,int64))return
    npoints=int(points)
    if(ns>npoints)return
    if(any(indices<1).or.any(indices>npoints))return
    if(ns==0)then
      status=0;return
    endif
    if(.not.salmon_all_finite(real(filter)).or..not.salmon_all_finite(aimag(filter)))return
    if(.not.salmon_all_finite(real(source)).or..not.salmon_all_finite(aimag(source)))return
    allocate(seen(npoints),stat=ios)
    if(ios/=0)then
      status=-3;return
    endif
    seen=.false.
    do i=1,ns
      if(seen(indices(i)))return
      seen(indices(i))=.true.
    enddo
    status=0
  end subroutine

  subroutine cufft_prepare(self,padded,indices,filter,source,capacity,status)
#ifdef USE_EXX_CUFFT
    use iso_c_binding, only: c_size_t,c_funloc
    use cufft, only: cufftCreate,cufftMakePlanMany,cufftSetStream,CUFFT_SUCCESS,CUFFT_Z2Z
    use openacc, only: acc_init,acc_device_nvidia,acc_get_num_devices,acc_get_device_num, &
      acc_get_cuda_stream,acc_async_sync
#endif
    implicit none
    class(s_exx_cufft),target,intent(inout) :: self
    integer,intent(in) :: padded(3),indices(:),capacity
    complex(8),intent(in) :: filter(:,:,:),source(:)
    integer,intent(out) :: status
    integer :: npoints,ns,cleanup
#ifdef USE_EXX_CUFFT
    integer :: device,dims(3),ios
    integer(c_size_t) :: workspace_bytes(1)
    integer,pointer,contiguous :: index_map(:)
    complex(8),pointer,contiguous :: kernel(:,:,:),source_buffer(:),work(:),target_buffer(:,:),output_buffer(:,:)
    logical :: rebuild
#endif
    call validate_source(padded,indices,filter,source,capacity,npoints,status)
    if(status/=0)then
      call self%release(cleanup)
      return
    endif
    ns=size(source)
    if(ns==0)then
      call self%release(status)
      if(status/=0)return
      self%padded=padded;self%ns=0;self%npoints=npoints;self%capacity=capacity;self%ready=.true.
      return
    endif
#ifndef USE_EXX_CUFFT
    call self%release(cleanup)
    status=-1
#else
    ! Fatal OpenACC errors abort the complete MPI job; registration precedes all runtime calls.
    call set_acc_error_routine(c_funloc(exx_cufft_fatal_error))
    if(acc_get_num_devices(acc_device_nvidia)<1)then
      call self%release(cleanup)
      status=-1;return
    endif
    call acc_init(acc_device_nvidia)
    device=acc_get_device_num(acc_device_nvidia)
    rebuild=.not.self%ready.or.any(self%padded/=padded).or.self%ns/=ns.or.self%capacity/=capacity.or.self%device/=device
    if(rebuild)then
      call self%release(status)
      if(status/=0)return
      self%padded=padded;self%ns=ns;self%npoints=npoints;self%capacity=capacity;self%device=device
      self%stream=acc_get_cuda_stream(acc_async_sync)
      allocate(self%indices(ns),self%filter(padded(1),padded(2),padded(3)),self%source(ns), &
        self%work(npoints*capacity),self%targets(ns,capacity),self%output(ns,capacity),stat=ios)
      if(ios/=0)then
        call self%release(cleanup)
        status=-3;return
      endif
      self%indices(:)=indices;self%filter(:,:,:)=filter;self%source(:)=source
      dims=padded(3:1:-1)
      status=cufftCreate(self%plan)
      if(status==CUFFT_SUCCESS)then
        self%plan_created=.true.
        status=cufftMakePlanMany(self%plan,3,dims,dims,1,npoints,dims,1,npoints,CUFFT_Z2Z,capacity,workspace_bytes)
        if(status==CUFFT_SUCCESS)then
          self%plan_builds=self%plan_builds+1
          status=cufftSetStream(self%plan,self%stream)
        endif
      endif
      if(status/=CUFFT_SUCCESS)then
        call self%release(cleanup)
        return
      endif
      ! Map only owned array storage. Kernels use aliases, never the polymorphic host object.
      index_map=>self%indices;kernel=>self%filter;source_buffer=>self%source
      work=>self%work;target_buffer=>self%targets;output_buffer=>self%output
      !$acc enter data copyin(index_map,kernel,source_buffer) create(work,target_buffer,output_buffer)
      self%resident=.true.
      self%index_uploads=self%index_uploads+1
      self%filter_uploads=self%filter_uploads+1
      self%source_uploads=self%source_uploads+1
    else
      index_map=>self%indices;kernel=>self%filter;source_buffer=>self%source
      if(any(self%indices/=indices))then
        self%indices(:)=indices
        !$acc update device(index_map)
        self%index_uploads=self%index_uploads+1
      endif
      if(any(self%filter/=filter))then
        self%filter(:,:,:)=filter
        !$acc update device(kernel)
        self%filter_uploads=self%filter_uploads+1
      endif
      if(any(self%source/=source))then
        self%source(:)=source
        !$acc update device(source_buffer)
        self%source_uploads=self%source_uploads+1
      endif
    endif
    self%ready=.true.;status=0
#endif
  end subroutine

  subroutine cufft_apply_batch(self,targets,action,status)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_EXX_CUFFT
    use cufft, only: cufftExecZ2Z,CUFFT_SUCCESS,CUFFT_FORWARD,CUFFT_INVERSE
    use openacc, only: acc_get_device_num,acc_device_nvidia
    use cudafor, only: cudaStreamSynchronize,cudaSuccess
#endif
    implicit none
    class(s_exx_cufft),target,intent(inout) :: self
    complex(8),intent(in) :: targets(:,:)
    complex(8),intent(out) :: action(:,:)
    integer,intent(out) :: status
    integer :: ns,nb,cleanup
#ifdef USE_EXX_CUFFT
    integer :: i,j,x,y,z,row,nx,ny,nz,npoints,total,capacity,ierr
    integer,pointer,contiguous :: index_map(:)
    complex(8),pointer,contiguous :: kernel(:,:,:),source_buffer(:),work(:),target_buffer(:,:),output_buffer(:,:)
    real(8) :: scale
#endif
    ! -1 unavailable; -2 invalid; -3 host allocation; -4 nonfinite output;
    ! -5 device changed (prepare again); <= -1000 CUDA error; positive cuFFT error.
    status=-2;action=0d0;ns=self%ns;nb=size(targets,2)
    if(.not.self%ready.or.size(targets,1)/=ns.or.nb>self%capacity.or.any(shape(action)/=[ns,nb]))then
      call self%release(cleanup)
      return
    endif
    if(ns==0.or.nb==0)then
      status=0;return
    endif
    if(.not.salmon_all_finite(real(targets)).or..not.salmon_all_finite(aimag(targets)))then
      call self%release(cleanup)
      return
    endif
#ifndef USE_EXX_CUFFT
    status=-1
#else
    if(acc_get_device_num(acc_device_nvidia)/=self%device)then
      call self%release(cleanup)
      status=-5;return
    endif
    index_map=>self%indices;kernel=>self%filter;source_buffer=>self%source
    work=>self%work;target_buffer=>self%targets;output_buffer=>self%output
    nx=self%padded(1);ny=self%padded(2);nz=self%padded(3)
    npoints=self%npoints;capacity=self%capacity;total=npoints*capacity;scale=1d0/dble(npoints)
    target_buffer(:,1:nb)=targets
    !$acc update device(target_buffer(:,1:nb))
    self%batch_uploads=self%batch_uploads+1
    ! Always clear inactive columns: the plan retains capacity, including the tail batch.
    !$acc parallel loop present(work)
    do i=1,total
      work(i)=(0d0,0d0)
    enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(2) present(index_map,source_buffer,target_buffer,work)
    do j=1,nb
      do i=1,ns
        work(index_map(i)+(j-1)*npoints)=conjg(source_buffer(i))*target_buffer(i,j)
      enddo
    enddo
    !$acc end parallel loop
    !$acc host_data use_device(work)
    status=cufftExecZ2Z(self%plan,work,work,CUFFT_FORWARD)
    !$acc end host_data
    ierr=cudaStreamSynchronize(self%stream)
    if(ierr/=cudaSuccess.and.status==CUFFT_SUCCESS)status=-1000-abs(ierr)
    if(status==CUFFT_SUCCESS)then
      !$acc parallel loop collapse(4) present(work,kernel) private(row)
      do j=1,capacity
        do z=1,nz
          do y=1,ny
            do x=1,nx
              row=x+nx*((y-1)+ny*(z-1))
              work(row+(j-1)*npoints)=work(row+(j-1)*npoints)*kernel(x,y,z)
            enddo
          enddo
        enddo
      enddo
      !$acc end parallel loop
      !$acc host_data use_device(work)
      status=cufftExecZ2Z(self%plan,work,work,CUFFT_INVERSE)
      !$acc end host_data
      ierr=cudaStreamSynchronize(self%stream)
      if(ierr/=cudaSuccess.and.status==CUFFT_SUCCESS)status=-1000-abs(ierr)
    endif
    if(status==CUFFT_SUCCESS)then
      !$acc parallel loop collapse(2) present(index_map,source_buffer,work,output_buffer)
      do j=1,nb
        do i=1,ns
          output_buffer(i,j)=-source_buffer(i)*work(index_map(i)+(j-1)*npoints)*scale
        enddo
      enddo
      !$acc end parallel loop
      ierr=cudaStreamSynchronize(self%stream)
      if(ierr/=cudaSuccess)status=-1000-abs(ierr)
    endif
    if(status==CUFFT_SUCCESS)then
      !$acc update self(output_buffer(:,1:nb))
      action=output_buffer(:,1:nb)
      if(.not.salmon_all_finite(real(action)).or..not.salmon_all_finite(aimag(action)))status=-4
    endif
    if(status/=CUFFT_SUCCESS)then
      action=0d0
      call self%release(cleanup)
    endif
#endif
  end subroutine

  subroutine cufft_release(self,status)
#ifdef USE_EXX_CUFFT
    use cufft, only: cufftDestroy,CUFFT_SUCCESS
    use openacc, only: acc_get_device_num,acc_set_device_num,acc_device_nvidia
    use cudafor, only: cudaStreamSynchronize,cudaSuccess
#endif
    implicit none
    class(s_exx_cufft),target,intent(inout) :: self
    integer,intent(out) :: status
#ifdef USE_EXX_CUFFT
    integer :: device,ierr
    integer,pointer,contiguous :: index_map(:)
    complex(8),pointer,contiguous :: kernel(:,:,:),source_buffer(:),work(:),target_buffer(:,:),output_buffer(:,:)
#endif
    status=0
#ifdef USE_EXX_CUFFT
    if(self%resident.or.self%plan_created)then
      device=acc_get_device_num(acc_device_nvidia)
      if(device/=self%device)call acc_set_device_num(self%device,acc_device_nvidia)
      ierr=cudaStreamSynchronize(self%stream)
      if(ierr/=cudaSuccess)status=-1000-abs(ierr)
      if(self%plan_created)then
        ierr=cufftDestroy(self%plan)
        if(status==0.and.ierr/=CUFFT_SUCCESS)status=ierr
      endif
      if(self%resident)then
        index_map=>self%indices;kernel=>self%filter;source_buffer=>self%source
        work=>self%work;target_buffer=>self%targets;output_buffer=>self%output
        !$acc exit data delete(index_map,kernel,source_buffer,work,target_buffer,output_buffer)
      endif
      if(device/=self%device)call acc_set_device_num(device,acc_device_nvidia)
    endif
#endif
    self%ready=.false.;self%resident=.false.;self%plan_created=.false.
    self%padded=0;self%ns=0;self%npoints=0;self%capacity=0;self%device=-1;self%plan=0;self%stream=0
    if(allocated(self%indices))deallocate(self%indices)
    if(allocated(self%filter))deallocate(self%filter)
    if(allocated(self%source))deallocate(self%source)
    if(allocated(self%work))deallocate(self%work)
    if(allocated(self%targets))deallocate(self%targets)
    if(allocated(self%output))deallocate(self%output)
  end subroutine

  subroutine cufft_finalize(self)
    implicit none
    type(s_exx_cufft),intent(inout) :: self
    integer :: status
    ! Normal callers use release to inspect its status; finalization prevents orphaned maps.
    call self%release(status)
  end subroutine

  subroutine exx_cufft_apply(padded,indices,filter,source,targets,action,status)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,intent(in) :: padded(3),indices(:)
    complex(8),intent(in) :: filter(:,:,:),source(:),targets(:,:)
    complex(8),intent(out) :: action(:,:)
    integer,intent(out) :: status
    type(s_exx_cufft),target :: backend
    integer :: nb,npoints,cleanup
    action=0d0;status=-2;nb=size(targets,2)
    if(size(targets,1)/=size(source).or.any(shape(action)/=shape(targets)))return
    if(nb==0.or.size(source)==0)then
      call validate_source(padded,indices,filter,source,max(1,nb),npoints,status)
      return
    endif
    if(.not.salmon_all_finite(real(targets)).or..not.salmon_all_finite(aimag(targets)))return
    call backend%prepare(padded,indices,filter,source,nb,status)
    if(status==0)call backend%apply(targets,action,status)
    call backend%release(cleanup)
    if(status==0)status=cleanup
    if(status/=0)action=0d0
  end subroutine
  ! Keep this host callback available for CPU-only MPI fault-injection tests.
  subroutine exx_cufft_fatal_error(message) bind(C,name='salmon_exx_cufft_abort')
    use iso_c_binding, only: c_ptr
    use iso_fortran_env, only: error_unit
#ifdef USE_MPI
    use mpi, only: MPI_Initialized,MPI_Finalized,MPI_Abort,MPI_COMM_WORLD,MPI_SUCCESS
#endif
    implicit none
    type(c_ptr),value,intent(in) :: message
    integer :: ios
#ifdef USE_MPI
    integer :: ierr
    logical :: initialized,finalized
#endif
    ! The runtime prints the detailed C error message. Do not allocate memory
    ! or access the failed device from this nonreturning host callback.
    write(error_unit,'(a)',iostat=ios)'EXX cuFFT: fatal OpenACC runtime error; terminating the job'
    flush(error_unit,iostat=ios)
#ifdef USE_MPI
    initialized=.false.;finalized=.true.
    call MPI_Initialized(initialized,ierr)
    if(ierr==MPI_SUCCESS.and.initialized)then
      call MPI_Finalized(finalized,ierr)
      if(ierr==MPI_SUCCESS.and..not.finalized)call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
#endif
    ! Returning from an OpenACC error callback has undefined behavior.
    error stop 'EXX cuFFT: fatal OpenACC runtime error'
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

  pure logical function finite_real_3d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:,:)
    real(8) :: value
    integer :: i,j,k
    finite=.false.
    do k=1,size(values,3)
      do j=1,size(values,2)
        do i=1,size(values,1)
          value=values(i,j,k)
          if(.not.ieee_is_finite(value))return
        enddo
      enddo
    enddo
    finite=.true.
  end function
end module exx_cufft
