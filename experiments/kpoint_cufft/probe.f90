module fftw_backend
 use iso_c_binding
 use exx_k_backend, only: s_exx_k_backend
 implicit none
 include 'fftw3.f03'
 logical :: fail_prepare=.false.,fail_apply=.false.
 integer :: release_count=0
 type,extends(s_exx_k_backend) :: oracle
  integer :: mesh(3),n(3),prepares=0,calls=0
  real(8),allocatable :: kernel(:,:,:)
  integer,allocatable :: point(:,:),shift(:,:),slot(:)
 contains
  procedure :: prepare=>prepare_oracle
  procedure :: apply=>apply_oracle
  procedure :: release=>release_oracle
 end type
contains
 subroutine create_oracle(backend)
  implicit none
  class(s_exx_k_backend),allocatable,intent(out) :: backend
  allocate(oracle::backend)
 end subroutine
 subroutine prepare_oracle(self,n,mesh,block,kernel,point,shift,slot,nslots,status)
  implicit none
  class(oracle),target,intent(inout) :: self
  integer,intent(in) :: n(3),mesh(3),block,point(:,:),shift(:,:),slot(:),nslots
  real(8),intent(in) :: kernel(0:,0:,0:)
  integer,intent(out) :: status
  self%prepares=self%prepares+1
  status=-101
  if(fail_prepare)return
  self%n=n;self%mesh=mesh;self%kernel=kernel
  self%point=point;self%shift=shift;self%slot=slot
  status=0
 end subroutine
 subroutine apply_oracle(self,lo,rows,buffer,action,status)
  implicit none
  class(oracle),target,intent(inout) :: self
  integer,intent(in) :: lo,rows
  complex(8),intent(in) :: buffer(:,:,:)
  complex(8),intent(out) :: action(:,:,:)
  integer,intent(out) :: status
  complex(8),allocatable :: work(:,:,:)
  type(c_ptr) :: forward,backward
  integer :: m(3),r,g,k,x,y,z,offset(3)
  self%calls=self%calls+1
  status=-102;action=0d0
  if(fail_apply)return
  m=self%mesh
  allocate(work(m(1),m(2),m(3)))
  forward=fftw_plan_dft_3d(m(3),m(2),m(1),work,work,FFTW_FORWARD,FFTW_ESTIMATE)
  backward=fftw_plan_dft_3d(m(3),m(2),m(1),work,work,FFTW_BACKWARD,FFTW_ESTIMATE)
  if(.not.c_associated(forward).or..not.c_associated(backward))error stop 'oracle plans'
  do g=1,size(buffer,2);do r=1,rows
   work=reshape(buffer(r,g,self%slot),m)
   call fftw_execute_dft(forward,work,work)
   k=0
   do z=1,m(3);do y=1,m(2);do x=1,m(1)
    k=k+1
    offset=modulo(self%point(:,lo+r-1)-self%point(:,g)-self%shift(:,k),self%n*m)
    work(x,y,z)=work(x,y,z)*self%kernel(offset(1)+lbound(self%kernel,1), &
      offset(2)+lbound(self%kernel,2),offset(3)+lbound(self%kernel,3))
   enddo;enddo;enddo
   call fftw_execute_dft(backward,work,work)
   action(r,g,self%slot)=reshape(work,[product(m)])
  enddo;enddo
  call fftw_destroy_plan(forward);call fftw_destroy_plan(backward)
  status=0
 end subroutine
 subroutine release_oracle(self,status)
  implicit none
  class(oracle),target,intent(inout) :: self
  integer,intent(out) :: status
  if(allocated(self%kernel))deallocate(self%kernel,self%point,self%shift,self%slot)
  release_count=release_count+1;status=0
 end subroutine
end module

program probe
 use mpi
 use exx_k_backend, only: k_backend_factory
 use exx_k_exchange, only: exx_k_kernel,exx_k_kernel_init,exx_k_kernel_apply_distributed,exx_k_kernel_destroy
 use fftw_backend, only: oracle,create_oracle,fail_prepare,fail_apply,release_count
#ifdef USE_EXX_CUFFT
 use exx_k_cufft, only: exx_k_cufft_create
#endif
 implicit none
 type(exx_k_kernel) :: cpu,device
 procedure(k_backend_factory),pointer :: factory
 integer :: rank,np,ierr,status,test_comm,n(3),m(3),nk,ng,p,j,g,pass,mode,block,x,y,z,i
 integer,allocatable :: starts(:),counts(:)
 real(8),allocatable :: k(:,:)
 complex(8),allocatable :: source(:,:,:),target(:,:,:),a(:,:,:),b(:,:,:)
 real(8) :: err,scale,h(3),omega
 call MPI_Init(ierr)
 call MPI_Comm_dup(MPI_COMM_WORLD,test_comm,ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 factory=>create_oracle
#ifdef USE_EXX_CUFFT
 factory=>exx_k_cufft_create
#endif
 allocate(starts(np),counts(np))
 do mode=1,4
  n=2;m=2;h=.7d0
  if(modulo(mode,2)==0)then
   n=[3,2,4];m=[2,3,1];h=[.6d0,.8d0,.7d0]
  endif
  ng=product(n);nk=product(m);block=5
  omega=.11d0
  if(mode>2)omega=0d0
  allocate(k(3,nk))
  i=0
  do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
   i=i+1;k(:,nk+1-i)=2d0*acos(-1d0)*[x,y,z]/(m*n*h)+[.03d0,-.02d0,.01d0]
  enddo;enddo;enddo
  do p=0,np-1
   starts(p+1)=1+p*nk/np;counts(p+1)=(p+1)*nk/np-p*nk/np
  enddo
  allocate(source(ng,2,counts(rank+1)),target(ng,3,counts(rank+1)), &
    a(ng,3,counts(rank+1)),b(ng,3,counts(rank+1)))
  do j=1,counts(rank+1);do g=1,ng
   i=starts(rank+1)+j-1
   source(g,1,j)=cmplx(sin(.31d0*g+.2d0*i),cos(.21d0*g-.13d0*i),8)
   source(g,2,j)=cmplx(cos(.17d0*g+.3d0*i),sin(.41d0*g),8)
   target(g,1:2,j)=source(g,:,j)*cmplx(.8d0,.1d0,8)
   target(g,3,j)=0d0
  enddo;enddo
  call exx_k_kernel_init(cpu,n,m,h,k,omega,block,status,starts(rank+1),counts(rank+1),fft_layout='strided')
  if(status/=0)error stop 'cpu init'
  call exx_k_kernel_init(device,n,m,h,k,omega,block,status,starts(rank+1),counts(rank+1),create_backend=factory)
  if(status/=0.or..not.allocated(device%accelerator))error stop 'device init'
  if(allocated(device%work))error stop 'unused CPU FFT workspace allocated'
  do pass=1,2
   call exx_k_kernel_apply_distributed(cpu,source,target,a,starts,counts,rank,test_comm,transpose_tiles,status)
   if(status/=0)error stop 'cpu apply'
   call exx_k_kernel_apply_distributed(device,source,target,b,starts,counts,rank,test_comm,transpose_tiles,status)
   if(status/=0)error stop 'device apply'
   scale=max(1d0,maxval(abs(a)));err=maxval(abs(a-b))/scale
   if(err>3d-12)error stop 'distributed backend parity'
   if(any(b(:,3,:)/=(0d0,0d0)))error stop 'zero target'
   source=source*cmplx(.7d0,.2d0,8);target=target*cmplx(.9d0,-.1d0,8)
  enddo
#ifndef USE_EXX_CUFFT
  select type(backend=>device%accelerator)
  type is(oracle)
   if(backend%prepares/=2)error stop 'prepare not once per action'
   if(backend%calls/=2*((ng+np*block-1)/(np*block)))error stop 'tile calls'
  end select
  fail_prepare=rank==0
  call exx_k_kernel_apply_distributed(device,source,target,b,starts,counts,rank,test_comm,transpose_tiles,status)
  if(status==0)error stop 'prepare failure not collective'
  fail_prepare=.false.;fail_apply=rank==0
  call exx_k_kernel_apply_distributed(device,source,target,b,starts,counts,rank,test_comm,transpose_tiles,status)
  if(status==0)error stop 'apply failure not collective'
  fail_apply=.false.
#endif
  if(np>1)then
   ! Mismatched backend selection must fail in the first fixed-size handshake.
   if(rank==0)then
    call exx_k_kernel_apply_distributed(cpu,source,target,a,starts,counts,rank,test_comm,transpose_tiles,status)
   else
    call exx_k_kernel_apply_distributed(device,source,target,b,starts,counts,rank,test_comm,transpose_tiles,status)
   endif
   if(status==0)error stop 'mixed backend accepted'
  endif
  call exx_k_kernel_destroy(cpu,status)
  if(status/=0)error stop 'cpu release'
  call exx_k_kernel_destroy(device,status)
  if(status/=0)error stop 'device release'
  deallocate(k,source,target,a,b)
 enddo
#ifndef USE_EXX_CUFFT
 if(release_count/=4)error stop 'release lifecycle'
#endif
 if(rank==0)print *, 'PASS distributed k backend parity, tails and lifecycle',np
 call MPI_Comm_free(test_comm,ierr)
 call MPI_Finalize(ierr)
contains
 subroutine transpose_tiles(send,recv,count,comm)
  implicit none
  complex(8),intent(in) :: send(:)
  complex(8),intent(out) :: recv(:)
  integer,intent(in) :: count,comm
  integer :: e
  if(comm/=test_comm.or.comm==MPI_COMM_NULL)error stop 'wrong explicit transpose communicator'
  call MPI_Alltoall(send,count,MPI_DOUBLE_COMPLEX,recv,count,MPI_DOUBLE_COMPLEX,comm,e)
  if(e/=MPI_SUCCESS)call MPI_Abort(MPI_COMM_WORLD,4,e)
 end subroutine
end program
