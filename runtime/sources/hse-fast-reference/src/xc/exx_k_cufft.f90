#include "config.h"
! Resident k-mesh FFT convolution of MPI density tiles. The caller performs
! orbital contractions and the 1/nk normalization; this backend does neither.
! Device selection follows OpenACC (ACC_DEVICE_NUM/CUDA_VISIBLE_DEVICES per rank).
! This type owns its device mappings and plan and must not be copied.
module exx_k_cufft
  use exx_k_backend, only: s_exx_k_backend
  use iso_c_binding, only: c_intptr_t
  implicit none
  private
  public :: s_exx_k_cufft,exx_k_cufft_create
  type,extends(s_exx_k_backend) :: s_exx_k_cufft
    private
    integer :: n(3)=0,mesh(3)=0,block=0,ng=0,nk=0,ns(3)=0,nslots=0,device=-1,plan=0
    integer(c_intptr_t) :: stream=0
    logical :: ready=.false.,resident=.false.,plan_created=.false.
    integer,allocatable :: point(:,:),shift(:,:),slot(:)
    real(8),allocatable :: kernel(:,:,:)
    complex(8),allocatable :: work(:),buffer(:,:,:),output(:,:,:)
    ! Lifetime totals. constant_uploads counts arrays, initially four; tiles count apply calls.
    integer,public :: plan_builds=0,constant_uploads=0,tile_uploads=0
  contains
    procedure :: prepare => k_cufft_prepare
    procedure :: apply => k_cufft_apply
    procedure :: release => k_cufft_release
    final :: k_cufft_finalize
  end type
#ifdef USE_EXX_CUFFT
  interface
    subroutine set_acc_error_routine(callback) bind(C,name='acc_set_error_routine')
      use iso_c_binding, only: c_funptr
      implicit none
      type(c_funptr),value,intent(in) :: callback
    end subroutine
    subroutine shared_fatal_error(message) bind(C,name='salmon_exx_cufft_abort')
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr),value,intent(in) :: message
    end subroutine
  end interface
#endif
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_3d
  interface salmon_all_finite
    module procedure finite_real_3d
  end interface
contains
  subroutine exx_k_cufft_create(backend)
    implicit none
    class(s_exx_k_backend),allocatable,intent(out) :: backend
    allocate(s_exx_k_cufft::backend)
  end subroutine

  subroutine validate_geometry(n,mesh,block,kernel,point,shift,slot,nslots,ng,nk,ns,status)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer,intent(in) :: n(3),mesh(3),block,point(:,:),shift(:,:),slot(:),nslots
    real(8),intent(in) :: kernel(0:,0:,0:)
    integer,intent(out) :: ng,nk,ns(3),status
    integer :: axis,batches,i,ios
    logical,allocatable :: used(:)
    status=-2;ng=1;nk=1;ns=0
    if(any(n<1).or.any(mesh<1).or.block<1.or.nslots<1)return
    ! Check before each multiplication, including the full MPI tile extent.
    do axis=1,3
      if(ng>huge(ng)/n(axis).or.nk>huge(nk)/mesh(axis))return
      ng=ng*n(axis);nk=nk*mesh(axis)
    enddo
    if(any(n>huge(1)/mesh))return
    ns=n*mesh
    if(nslots<nk)return
    if(block>huge(batches)/ng)return
    batches=block*ng
    if(nslots>huge(batches)/batches)return
    ! Since nslots >= nk, this also bounds the contiguous cuFFT work size.
    if(any(shape(kernel)/=ns))return
    if(any(shape(point)/=[3,ng]).or.any(shape(shift)/=[3,nk]).or.size(slot)/=nk)return
    do axis=1,3
      if(any(point(axis,:)<0).or.any(point(axis,:)>=n(axis)))return
      if(any(shift(axis,:)<0).or.any(shift(axis,:)>ns(axis)-n(axis)))return
      if(any(modulo(shift(axis,:),n(axis))/=0))return
    enddo
    if(any(slot<1).or.any(slot>nslots))return
    if(.not.salmon_all_finite(kernel))return
    allocate(used(nslots),stat=ios)
    if(ios/=0)then
      status=-3;return
    endif
    used=.false.
    do i=1,nk
      if(used(slot(i)))return
      used(slot(i))=.true.
    enddo
    status=0
  end subroutine

  subroutine k_cufft_prepare(self,n,mesh,block,kernel,point,shift,slot,nslots,status)
#ifdef USE_EXX_CUFFT
    use iso_c_binding, only: c_size_t,c_funloc
    use cufft, only: cufftCreate,cufftMakePlanMany,cufftSetStream,CUFFT_SUCCESS,CUFFT_Z2Z
    use openacc, only: acc_init,acc_device_nvidia,acc_get_num_devices,acc_get_device_num, &
      acc_get_cuda_stream,acc_async_sync
#endif
    implicit none
    class(s_exx_k_cufft),target,intent(inout) :: self
    integer,intent(in) :: n(3),mesh(3),block,point(:,:),shift(:,:),slot(:),nslots
    real(8),intent(in) :: kernel(0:,0:,0:)
    integer,intent(out) :: status
    integer :: ng,nk,ns(3),cleanup
#ifdef USE_EXX_CUFFT
    integer :: device,dims(3),ios
    integer(c_size_t) :: workspace_bytes(1)
    integer,pointer,contiguous :: point_map(:,:),shift_map(:,:),slot_map(:)
    real(8),pointer,contiguous :: kernel_gpu(:,:,:)
    complex(8),pointer,contiguous :: work(:),tile(:,:,:),output(:,:,:)
    logical :: rebuild
#endif
    call validate_geometry(n,mesh,block,kernel,point,shift,slot,nslots,ng,nk,ns,status)
    if(status/=0)then
      call self%release(cleanup)
      return
    endif
#ifndef USE_EXX_CUFFT
    ! Retain validated geometry so zero-row calls and argument checks stay usable in the stub.
    call self%release(cleanup)
    self%n=n;self%mesh=mesh;self%block=block;self%ng=ng;self%nk=nk;self%ns=ns;self%nslots=nslots
    self%ready=.true.;status=-1
#else
    call set_acc_error_routine(c_funloc(shared_fatal_error))
    if(acc_get_num_devices(acc_device_nvidia)<1)then
      call self%release(cleanup)
      status=-1;return
    endif
    call acc_init(acc_device_nvidia)
    device=acc_get_device_num(acc_device_nvidia)
    rebuild=.not.self%ready.or.any(self%n/=n).or.any(self%mesh/=mesh).or.self%block/=block.or. &
      self%nslots/=nslots.or.self%device/=device
    if(rebuild)then
      call self%release(status)
      if(status/=0)return
      self%n=n;self%mesh=mesh;self%block=block;self%ng=ng;self%nk=nk;self%ns=ns;self%nslots=nslots
      self%device=device;self%stream=acc_get_cuda_stream(acc_async_sync)
      allocate(self%kernel(0:ns(1)-1,0:ns(2)-1,0:ns(3)-1),self%point(3,ng),self%shift(3,nk),self%slot(nk), &
        self%work(nk*block*ng),self%buffer(block,ng,nslots),self%output(block,ng,nslots),stat=ios)
      if(ios/=0)then
        call self%release(cleanup)
        status=-3;return
      endif
      self%kernel(:,:,:)=kernel;self%point(:,:)=point;self%shift(:,:)=shift;self%slot(:)=slot
      dims=mesh(3:1:-1)
      status=cufftCreate(self%plan)
      if(status==CUFFT_SUCCESS)then
        self%plan_created=.true.
        ! Contiguous (nk,block,ng): one mesh^3 transform for every (row,grid point).
        status=cufftMakePlanMany(self%plan,3,dims,dims,1,nk,dims,1,nk,CUFFT_Z2Z,block*ng,workspace_bytes)
        if(status==CUFFT_SUCCESS)then
          self%plan_builds=self%plan_builds+1
          status=cufftSetStream(self%plan,self%stream)
        endif
      endif
      if(status/=CUFFT_SUCCESS)then
        call self%release(cleanup)
        return
      endif
      point_map=>self%point;shift_map=>self%shift;slot_map=>self%slot;kernel_gpu=>self%kernel
      work=>self%work;tile=>self%buffer;output=>self%output
      ! Only owned array storage is mapped; GPU kernels do not reference the host object.
      !$acc enter data copyin(point_map,shift_map,slot_map,kernel_gpu) create(work,tile,output)
      self%resident=.true.;self%constant_uploads=self%constant_uploads+4
    else
      point_map=>self%point;shift_map=>self%shift;slot_map=>self%slot;kernel_gpu=>self%kernel
      if(any(self%point/=point))then
        self%point(:,:)=point
        !$acc update device(point_map)
        self%constant_uploads=self%constant_uploads+1
      endif
      if(any(self%shift/=shift))then
        self%shift(:,:)=shift
        !$acc update device(shift_map)
        self%constant_uploads=self%constant_uploads+1
      endif
      if(any(self%slot/=slot))then
        self%slot(:)=slot
        !$acc update device(slot_map)
        self%constant_uploads=self%constant_uploads+1
      endif
      if(any(self%kernel/=kernel))then
        self%kernel(:,:,:)=kernel
        !$acc update device(kernel_gpu)
        self%constant_uploads=self%constant_uploads+1
      endif
    endif
    self%ready=.true.;status=0
#endif
  end subroutine

  subroutine k_cufft_apply(self,lo,rows,buffer,action,status)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_EXX_CUFFT
    use cufft, only: cufftExecZ2Z,CUFFT_SUCCESS,CUFFT_FORWARD,CUFFT_INVERSE
    use openacc, only: acc_get_device_num,acc_device_nvidia
    use cudafor, only: cudaStreamSynchronize,cudaSuccess
#endif
    implicit none
    class(s_exx_k_cufft),target,intent(inout) :: self
    integer,intent(in) :: lo,rows
    complex(8),intent(in) :: buffer(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    integer,intent(out) :: status
    integer :: cleanup
#ifdef USE_EXX_CUFFT
    integer :: b,ng,nk,nsx,nsy,nsz,nslots,i,r,g,ki,total,offset_x,offset_y,offset_z,flat,ierr
    integer,pointer,contiguous :: point_map(:,:),shift_map(:,:),slot_map(:)
    real(8),pointer,contiguous :: kernel_gpu(:,:,:)
    complex(8),pointer,contiguous :: work(:),tile(:,:,:),output(:,:,:)
#endif
    ! -1 unavailable; -2 invalid; -3 allocation; -4 nonfinite result; -5 device changed;
    ! <= -1000 CUDA runtime errors; positive values are cuFFT errors.
    status=-2;action=0d0
    if(.not.self%ready.or.any(shape(buffer)/=[self%block,self%ng,self%nslots]).or. &
      any(shape(action)/=shape(buffer)))then
      call self%release(cleanup)
      return
    endif
    if(lo<1.or.rows<0.or.rows>self%block)then
      call self%release(cleanup)
      return
    endif
    if(rows>0)then
      if(lo>self%ng)then
        call self%release(cleanup)
        return
      endif
      if(rows>self%ng-lo+1)then
        call self%release(cleanup)
        return
      endif
    endif
    if(.not.salmon_all_finite(real(buffer)).or..not.salmon_all_finite(aimag(buffer)))then
      call self%release(cleanup)
      return
    endif
    if(rows==0)then
      status=0;return
    endif
#ifndef USE_EXX_CUFFT
    status=-1
#else
    if(acc_get_device_num(acc_device_nvidia)/=self%device)then
      call self%release(cleanup)
      status=-5;return
    endif
    point_map=>self%point;shift_map=>self%shift;slot_map=>self%slot;kernel_gpu=>self%kernel
    work=>self%work;tile=>self%buffer;output=>self%output
    b=self%block;ng=self%ng;nk=self%nk;nsx=self%ns(1);nsy=self%ns(2);nsz=self%ns(3);nslots=self%nslots;total=nk*b*ng
    tile=buffer
    !$acc update device(tile)
    self%tile_uploads=self%tile_uploads+1
    !$acc parallel loop present(work)
    do i=1,total
      work(i)=(0d0,0d0)
    enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(3) present(output)
    do ki=1,nslots
      do g=1,ng
        do r=1,b
          output(r,g,ki)=(0d0,0d0)
        enddo
      enddo
    enddo
    !$acc end parallel loop
    !$acc parallel loop collapse(3) present(work,tile,slot_map)
    do g=1,ng
      do r=1,rows
        do ki=1,nk
          work(ki+nk*((r-1)+b*(g-1)))=tile(r,g,slot_map(ki))
        enddo
      enddo
    enddo
    !$acc end parallel loop
    !$acc host_data use_device(work)
    status=cufftExecZ2Z(self%plan,work,work,CUFFT_FORWARD)
    !$acc end host_data
    ierr=cudaStreamSynchronize(self%stream)
    if(ierr/=cudaSuccess.and.status==CUFFT_SUCCESS)status=-1000-abs(ierr)
    if(status==CUFFT_SUCCESS)then
      !$acc parallel loop collapse(3) present(work,kernel_gpu,point_map,shift_map) &
      !$acc& private(offset_x,offset_y,offset_z,flat)
      do g=1,ng
        do r=1,rows
          do ki=1,nk
            offset_x=modulo(point_map(1,lo+r-1)-point_map(1,g)-shift_map(1,ki),nsx)
            offset_y=modulo(point_map(2,lo+r-1)-point_map(2,g)-shift_map(2,ki),nsy)
            offset_z=modulo(point_map(3,lo+r-1)-point_map(3,g)-shift_map(3,ki),nsz)
            flat=ki+nk*((r-1)+b*(g-1))
            work(flat)=work(flat)*kernel_gpu(offset_x,offset_y,offset_z)
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
      !$acc parallel loop collapse(3) present(work,output,slot_map)
      do g=1,ng
        do r=1,rows
          do ki=1,nk
            ! Deliberately unnormalized, matching FFTW before the caller's -1/nk BLAS contraction.
            output(r,g,slot_map(ki))=work(ki+nk*((r-1)+b*(g-1)))
          enddo
        enddo
      enddo
      !$acc end parallel loop
      ierr=cudaStreamSynchronize(self%stream)
      if(ierr/=cudaSuccess)status=-1000-abs(ierr)
    endif
    if(status==CUFFT_SUCCESS)then
      !$acc update self(output)
      action=output
      if(.not.salmon_all_finite(real(action)).or..not.salmon_all_finite(aimag(action)))status=-4
    endif
    if(status/=CUFFT_SUCCESS)then
      action=0d0
      call self%release(cleanup)
    endif
#endif
  end subroutine

  subroutine k_cufft_release(self,status)
#ifdef USE_EXX_CUFFT
    use cufft, only: cufftDestroy,CUFFT_SUCCESS
    use openacc, only: acc_get_device_num,acc_set_device_num,acc_device_nvidia
    use cudafor, only: cudaStreamSynchronize,cudaSuccess
#endif
    implicit none
    class(s_exx_k_cufft),target,intent(inout) :: self
    integer,intent(out) :: status
#ifdef USE_EXX_CUFFT
    integer :: device,ierr
    integer,pointer,contiguous :: point_map(:,:),shift_map(:,:),slot_map(:)
    real(8),pointer,contiguous :: kernel_gpu(:,:,:)
    complex(8),pointer,contiguous :: work(:),tile(:,:,:),output(:,:,:)
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
        point_map=>self%point;shift_map=>self%shift;slot_map=>self%slot;kernel_gpu=>self%kernel
        work=>self%work;tile=>self%buffer;output=>self%output
        !$acc exit data delete(point_map,shift_map,slot_map,kernel_gpu,work,tile,output)
      endif
      if(device/=self%device)call acc_set_device_num(device,acc_device_nvidia)
    endif
#endif
    self%ready=.false.;self%resident=.false.;self%plan_created=.false.
    self%n=0;self%mesh=0;self%block=0;self%ng=0;self%nk=0;self%ns=0;self%nslots=0
    self%device=-1;self%plan=0;self%stream=0
    if(allocated(self%kernel))deallocate(self%kernel)
    if(allocated(self%point))deallocate(self%point)
    if(allocated(self%shift))deallocate(self%shift)
    if(allocated(self%slot))deallocate(self%slot)
    if(allocated(self%work))deallocate(self%work)
    if(allocated(self%buffer))deallocate(self%buffer)
    if(allocated(self%output))deallocate(self%output)
  end subroutine

  subroutine k_cufft_finalize(self)
    implicit none
    type(s_exx_k_cufft),intent(inout) :: self
    integer :: status
    call self%release(status)
  end subroutine

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
end module exx_k_cufft
