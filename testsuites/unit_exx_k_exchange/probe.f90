program probe
  use mpi
  use exx_k_exchange, only: exx_k_kernel,exx_k_kernel_init,exx_k_kernel_apply_distributed,exx_k_kernel_destroy
  use exx_wannier, only: s_exx_wannier,wannier_init,wannier_refresh_source,wannier_apply,wannier_destroy
  implicit none
  type(exx_k_kernel) :: kernel
  type(s_exx_wannier) :: reference
  integer :: rank,np,ierr,status,n(3),mesh(3),nk,ng,no,nt,i,j,g,x,y,z,mode,first,last
  integer,allocatable :: starts(:),counts(:)
  real(8) :: h(3),omega,radius,err,scale
  real(8),allocatable :: k(:,:),occ(:,:)
  complex(8),allocatable :: source(:,:,:),target(:,:,:),expected(:,:,:),actual(:,:,:)
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  n=[3,2,4];mesh=[2,3,1];h=[.6d0,.8d0,.7d0];ng=product(n);nk=product(mesh);no=2;nt=3
  allocate(k(3,nk),occ(no,nk),source(ng,no,nk),target(ng,nt,nk),expected(ng,nt,nk))
  allocate(starts(np),counts(np))
  do i=0,np-1
    starts(i+1)=1+i*nk/np;counts(i+1)=(i+1)*nk/np-i*nk/np
  enddo
  first=starts(rank+1);last=first+counts(rank+1)-1
  allocate(actual(ng,nt,counts(rank+1)))
  i=0
  do z=0,mesh(3)-1;do y=0,mesh(2)-1;do x=0,mesh(1)-1
    i=i+1;k(:,nk+1-i)=2*acos(-1d0)*[x,y,z]/(mesh*n*h)+[.03d0,-.02d0,.01d0]
  enddo;enddo;enddo
  do j=1,nk;do g=1,ng
    source(g,1,j)=cmplx(sin(.31d0*g+.2d0*j),cos(.21d0*g-.13d0*j),8)
    source(g,2,j)=cmplx(cos(.17d0*g+.3d0*j),sin(.41d0*g),8)
    target(g,1:2,j)=source(g,:,j)*cmplx(.8d0,.1d0,8)
    target(g,3,j)=cmplx(sin(.11d0*g),cos(.7d0*j),8)
  enddo;enddo
  occ=2d0
  do mode=1,3
    omega=.11d0;radius=0d0
    if(mode>=2)omega=0d0
    if(mode==3)radius=.9d0
    call wannier_init(reference,n,mesh,h,k,omega,status,radius)
    if(status/=0)error stop 'reference init'
    call wannier_refresh_source(reference,source,occ,1,1d-7,status,localize=.false.)
    if(status/=0)error stop 'reference update'
    call wannier_apply(reference,target,expected,status)
    if(status/=0)error stop 'reference action'
    call exx_k_kernel_init(kernel,n,mesh,h,k,omega,5,status,first,counts(rank+1),coulomb_radius=radius)
    if(status/=0)error stop 'distributed init'
    call exx_k_kernel_apply_distributed(kernel,source(:,:,first:last),target(:,:,first:last),actual, &
      starts,counts,rank,transpose_tiles,status)
    if(status/=0)error stop 'distributed action'
    err=maxval(abs(actual-expected(:,:,first:last)));scale=max(1d0,maxval(abs(expected)))
    if(err>2d-11*scale)then
      print *, 'mismatch',mode,rank,err,scale
      error stop 'exchange action differs from full-support Wannier oracle'
    endif
    if(rank==0)print *, 'PASS kernel/np/error:',mode,np,err
    call exx_k_kernel_destroy(kernel)
    call wannier_destroy(reference)
  enddo
  call MPI_Finalize(ierr)
contains
  subroutine transpose_tiles(send,recv,count)
    implicit none
    complex(8),intent(in) :: send(:)
    complex(8),intent(out) :: recv(:)
    integer,intent(in) :: count
    integer :: stat
    call MPI_Alltoall(send,count,MPI_DOUBLE_COMPLEX,recv,count,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,stat)
    if(stat/=MPI_SUCCESS)error stop 'test transpose'
  end subroutine
end program
