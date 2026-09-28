#include "config.h"
! Optional batch backend for the existing compact, discrete exchange kernel.
! Each MPI rank uses its OpenACC-selected NVIDIA device. Assign ACC_DEVICE_NUM
! (or restrict CUDA_VISIBLE_DEVICES) per rank; this routine does not map ranks.
module exx_cufft
  implicit none
  private
  public :: exx_cufft_apply
#ifdef USE_EXX_CUFFT
  interface
    subroutine set_acc_error_routine(callback) bind(C,name='acc_set_error_routine')
      use iso_c_binding, only: c_funptr
      implicit none
      type(c_funptr),value,intent(in) :: callback
    end subroutine
  end interface
#endif
contains
  subroutine exx_cufft_apply(padded,indices,filter,source,targets,action,status)
    use iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_EXX_CUFFT
    use iso_c_binding, only: c_size_t,c_intptr_t,c_funloc
    use cufft, only: cufftCreate,cufftMakePlanMany,cufftSetStream,cufftExecZ2Z,cufftDestroy, &
      CUFFT_SUCCESS,CUFFT_Z2Z,CUFFT_FORWARD,CUFFT_INVERSE
    use openacc, only: acc_init,acc_device_nvidia,acc_get_num_devices,acc_get_cuda_stream,acc_async_sync
    use cudafor, only: cudaStreamSynchronize,cudaSuccess
#endif
    implicit none
    integer,intent(in) :: padded(3),indices(:)
    complex(8),intent(in) :: filter(:,:,:),source(:),targets(:,:)
    complex(8),intent(out) :: action(:,:)
    integer,intent(out) :: status
    integer :: ns,nb,npoints,axis,i,ios
    integer(int64) :: points
    logical,allocatable :: seen(:)
#ifdef USE_EXX_CUFFT
    complex(8),allocatable :: work(:)
    integer :: plan,dims(3),ierr,cleanup,total,j,x,y,z,row,nx,ny,nz
    integer(c_size_t) :: workspace_bytes(1)
    integer(c_intptr_t) :: stream
    real(8) :: scale
#endif
    ! Negative codes are local: -1 unavailable, -2 invalid input, -3 allocation,
    ! -4 nonfinite result, <= -1000 CUDA runtime errors. cuFFT codes are positive.
    status=-2;action=0d0
    if(any(padded<1))return
    if(any(shape(filter)/=padded))return
    if(size(source,kind=int64)>int(huge(0),int64))return
    if(size(targets,2,kind=int64)>int(huge(0),int64))return
    ns=size(source);nb=size(targets,2)
    if(size(indices)/=ns.or.size(targets,1)/=ns)return
    if(any(shape(action)/=[ns,nb]))return
    ! The 32-bit cuFFT batch layout and device indexing must not overflow.
    points=1_int64
    do axis=1,3
      if(points>int(huge(0),int64)/int(padded(axis),int64))return
      points=points*int(padded(axis),int64)
    enddo
    npoints=int(points)
    if(ns>npoints)return
    if(nb>0)then
      if(points>int(huge(0),int64)/int(nb,int64))return
    endif
    if(any(indices<1).or.any(indices>npoints))return
    if(ns==0.or.nb==0)then
      status=0;return
    endif
    if(.not.all(ieee_is_finite(real(filter))).or..not.all(ieee_is_finite(aimag(filter))))return
    if(.not.all(ieee_is_finite(real(source))).or..not.all(ieee_is_finite(aimag(source))))return
    if(.not.all(ieee_is_finite(real(targets))).or..not.all(ieee_is_finite(aimag(targets))))return
    allocate(seen(npoints),stat=ios)
    if(ios/=0)then
      status=-3;return
    endif
    seen=.false.
    do i=1,ns
      if(seen(indices(i)))return
      seen(indices(i))=.true.
    enddo
    deallocate(seen)
#ifndef USE_EXX_CUFFT
    status=-1
#else
    status=-1
    ! A fatal OpenACC error cannot return a status. Terminate the complete MPI
    ! job rather than leaving other ranks waiting in exchange collectives.
    call set_acc_error_routine(c_funloc(exx_cufft_fatal_error))
    if(acc_get_num_devices(acc_device_nvidia)<1)return
    ! Initialize the device selected by the OpenACC environment, without changing its number.
    call acc_init(acc_device_nvidia)
    stream=acc_get_cuda_stream(acc_async_sync)
    total=npoints*nb
    allocate(work(total),stat=ios)
    if(ios/=0)then
      status=-3;return
    endif
    ! Reverse Fortran dimensions for cuFFT's row-major rank-3 layout.
    dims=padded(3:1:-1)
    status=cufftCreate(plan)
    if(status/=CUFFT_SUCCESS)return
    status=cufftMakePlanMany(plan,3,dims,dims,1,npoints,dims,1,npoints,CUFFT_Z2Z,nb,workspace_bytes)
    if(status==CUFFT_SUCCESS)status=cufftSetStream(plan,stream)
    if(status==CUFFT_SUCCESS)then
      scale=1d0/dble(npoints)
      nx=padded(1);ny=padded(2);nz=padded(3)
      ! Inputs and the filter are uploaded once for the whole batch. Intermediate
      ! pair densities, Fourier coefficients and potentials remain on the device.
      !$acc data copyin(indices,filter,source,targets) create(work) copyout(action)
      !$acc parallel loop present(work)
      do i=1,total
        work(i)=(0d0,0d0)
      enddo
      !$acc end parallel loop
      !$acc parallel loop collapse(2) present(indices,source,targets,work,action)
      do j=1,nb
        do i=1,ns
          action(i,j)=(0d0,0d0)
          work(indices(i)+(j-1)*npoints)=conjg(source(i))*targets(i,j)
        enddo
      enddo
      !$acc end parallel loop
      ! Library launches use the same stream as synchronous OpenACC kernels.
      !$acc host_data use_device(work)
      status=cufftExecZ2Z(plan,work,work,CUFFT_FORWARD)
      !$acc end host_data
      ierr=cudaStreamSynchronize(stream)
      if(ierr/=cudaSuccess.and.status==CUFFT_SUCCESS)status=-1000-abs(ierr)
      if(status==CUFFT_SUCCESS)then
        !$acc parallel loop collapse(4) present(work,filter) private(row)
        do j=1,nb
          do z=1,nz
            do y=1,ny
              do x=1,nx
                row=x+nx*((y-1)+ny*(z-1))
                work(row+(j-1)*npoints)=work(row+(j-1)*npoints)*filter(x,y,z)
              enddo
            enddo
          enddo
        enddo
        !$acc end parallel loop
        !$acc host_data use_device(work)
        status=cufftExecZ2Z(plan,work,work,CUFFT_INVERSE)
        !$acc end host_data
        ierr=cudaStreamSynchronize(stream)
        if(ierr/=cudaSuccess.and.status==CUFFT_SUCCESS)status=-1000-abs(ierr)
      endif
      if(status==CUFFT_SUCCESS)then
        !$acc parallel loop collapse(2) present(indices,source,work,action)
        do j=1,nb
          do i=1,ns
            action(i,j)=-source(i)*work(indices(i)+(j-1)*npoints)*scale
          enddo
        enddo
        !$acc end parallel loop
      endif
      ierr=cudaStreamSynchronize(stream)
      if(ierr/=cudaSuccess.and.status==CUFFT_SUCCESS)status=-1000-abs(ierr)
      !$acc end data
    endif
    cleanup=cufftDestroy(plan)
    if(status==CUFFT_SUCCESS.and.cleanup/=CUFFT_SUCCESS)status=cleanup
    if(status==CUFFT_SUCCESS)then
      if(.not.all(ieee_is_finite(real(action))).or..not.all(ieee_is_finite(aimag(action))))status=-4
    endif
    if(status/=CUFFT_SUCCESS)action=0d0
#endif
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
end module exx_cufft
