program gpu_probe
 use iso_c_binding
 use exx_k_backend, only: s_exx_k_backend
 use exx_k_cufft, only: s_exx_k_cufft,exx_k_cufft_create
 implicit none
 include 'fftw3.f03'
 integer,parameter :: n=2,mesh=3,ng=n**3,nk=mesh**3,ns=n*mesh
 class(s_exx_k_backend),allocatable :: backend
 integer :: point(3,ng),shift(3,nk),slot(nk),expected(3),status,x,y,z,i,block,nslots,swap(3),slot_swap
 real(8) :: kernel(0:ns-1,0:ns-1,0:ns-1)
 complex(8),allocatable :: buffer(:,:,:),action(:,:,:),reference(:,:,:)
 block=3;nslots=nk+3;i=0
 do z=0,n-1;do y=0,n-1;do x=0,n-1
  i=i+1;point(:,i)=[x,y,z]
 enddo;enddo;enddo
 i=0
 do z=0,mesh-1;do y=0,mesh-1;do x=0,mesh-1
  i=i+1;shift(:,i)=n*[x,y,z]
 enddo;enddo;enddo
 ! An odd 3^3 mesh distinguishes forward from inverse FFT signs; padded slots are permuted.
 do i=1,nk
  slot(i)=1+modulo(7*(i-1)+4,nslots)
 enddo
 do z=0,ns-1;do y=0,ns-1;do x=0,ns-1
  kernel(x,y,z)=.2d0+.017d0*x+.031d0*y+.047d0*z+.006d0*x*y*z
 enddo;enddo;enddo
 call reset_buffers()
 call exx_k_cufft_create(backend)
 expected=0;call check_counters()
 call prepare()
 expected=[1,4,0];call check_counters()
 call compare_apply(1,3)
 call compare_apply(7,2) ! A non-full final row tile must clear the unused row.
 call compare_apply(99,0) ! Idle MPI workers can have a lower bound past the mesh.
 call prepare()
 call check_counters() ! Identical prepare must leave all resident data untouched.
 ! The density tile contains the contracted source: changing it uploads only a tile.
 buffer(1,2,slot(3))=buffer(1,2,slot(3))*cmplx(1.3d0,.4d0,8)
 call compare_apply(3,3)
 kernel=kernel*.8d0+.017d0
 call prepare()
 expected(2)=expected(2)+1;call check_counters()
 call compare_apply(2,3)
 swap=point(:,1);point(:,1)=point(:,3);point(:,3)=swap
 call prepare()
 expected(2)=expected(2)+1;call check_counters()
 call compare_apply(1,3)
 swap=shift(:,2);shift(:,2)=shift(:,7);shift(:,7)=swap
 call prepare()
 expected(2)=expected(2)+1;call check_counters()
 call compare_apply(6,3)
 slot_swap=slot(1);slot(1)=slot(8);slot(8)=slot_swap
 call prepare()
 expected(2)=expected(2)+1;call check_counters()
 call compare_apply(4,3)
 block=2
 call reset_buffers()
 call prepare()
 expected(1:2)=expected(1:2)+[1,4];call check_counters()
 call compare_apply(8,1)
 nslots=nslots+1
 call reset_buffers()
 call prepare()
 expected(1:2)=expected(1:2)+[1,4];call check_counters()
 call compare_apply(2,2)
 call compare_apply(99,0)
 call backend%release(status)
 if(status/=0)error stop 'k-cuFFT explicit release failed'
 call check_counters()
 call backend%release(status)
 if(status/=0)error stop 'k-cuFFT repeated release failed'
 call check_counters()
 call prepare()
 expected(1:2)=expected(1:2)+[1,4];call check_counters()
 call compare_apply(5,2)
 deallocate(backend) ! Finalization of a prepared resident object.
 call exx_k_cufft_create(backend)
 expected=0;call check_counters()
 call prepare()
 expected=[1,4,0];call check_counters()
 call compare_apply(7,2)
 deallocate(backend)
 print *, 'PASS k-cuFFT/FFTW parity: padded slots, row tails, zero rows, uploads, rebuilds and finalization'
contains
 subroutine prepare()
  implicit none
  call backend%prepare([n,n,n],[mesh,mesh,mesh],block,kernel,point,shift,slot,nslots,status)
  if(status/=0)error stop 'k-cuFFT prepare failed (NVHPC and NVIDIA device required)'
 end subroutine

 subroutine reset_buffers()
  implicit none
  integer :: r,g,s
  if(allocated(buffer))deallocate(buffer,action,reference)
  allocate(buffer(block,ng,nslots),action(block,ng,nslots),reference(block,ng,nslots))
  ! Padding deliberately contains large nonzero values: only slot(:) may enter the transform.
  buffer=cmplx(99d0,-77d0,8)
  do s=1,nk;do g=1,ng;do r=1,block
   buffer(r,g,slot(s))=cmplx(sin(.19d0*r+.31d0*g+.27d0*s),cos(.17d0*r-.13d0*g+.23d0*s),8)
  enddo;enddo;enddo
 end subroutine

 subroutine check_counters()
  implicit none
  integer :: actual(3)
  select type(backend)
  type is(s_exx_k_cufft)
   actual=[backend%plan_builds,backend%constant_uploads,backend%tile_uploads]
  class default
   error stop 'Wrong k-cuFFT factory dynamic type'
  end select
  if(any(actual/=expected))then
   print *, 'k-cuFFT counters actual/expected: ',actual,expected
   error stop 'k-cuFFT resident counter mismatch'
  endif
 end subroutine

 subroutine compare_apply(lo,rows)
  implicit none
  integer,intent(in) :: lo,rows
  integer :: r,g,k,ix,iy,iz,offset(3),s
  complex(c_double_complex) :: work(mesh,mesh,mesh)
  type(c_ptr) :: forward,backward
  real(8) :: scale
  reference=0d0
  forward=fftw_plan_dft_3d(mesh,mesh,mesh,work,work,FFTW_FORWARD,FFTW_ESTIMATE)
  backward=fftw_plan_dft_3d(mesh,mesh,mesh,work,work,FFTW_BACKWARD,FFTW_ESTIMATE)
  if(.not.c_associated(forward).or..not.c_associated(backward))error stop 'FFTW oracle plan failed'
  do g=1,ng;do r=1,rows
   work=reshape(buffer(r,g,slot),[mesh,mesh,mesh])
   call fftw_execute_dft(forward,work,work)
   k=0
   do iz=1,mesh;do iy=1,mesh;do ix=1,mesh
    k=k+1
    offset=modulo(point(:,lo+r-1)-point(:,g)-shift(:,k),ns)
    work(ix,iy,iz)=work(ix,iy,iz)*kernel(offset(1),offset(2),offset(3))
   enddo;enddo;enddo
   call fftw_execute_dft(backward,work,work)
   ! Deliberately no division by mesh**3: the caller owns this normalization.
   reference(r,g,slot)=reshape(work,[nk])
  enddo;enddo
  call fftw_destroy_plan(forward);call fftw_destroy_plan(backward)
  action=cmplx(999d0,-999d0,8)
  call backend%apply(lo,rows,buffer,action,status)
  if(status/=0)error stop 'k-cuFFT resident apply failed'
  if(rows>0)expected(3)=expected(3)+1
  call check_counters()
  scale=max(1d0,maxval(abs(reference)))
  if(maxval(abs(action-reference))/scale>3d-12)error stop 'k-cuFFT/FFTW unnormalized action mismatch'
  if(rows<block)then
   if(any(action(rows+1:block,:,:)/=(0d0,0d0)))error stop 'Inactive row was not cleared'
  endif
  do s=1,nslots
   if(any(slot==s))cycle
   if(any(action(:,:,s)/=(0d0,0d0)))error stop 'Padding slot was not cleared'
  enddo
 end subroutine
end program
