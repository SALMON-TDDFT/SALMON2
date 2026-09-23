! Full periodic sampled HSE kernel. Unit one-spin source occupations;
! neither the hybrid mixing fraction nor a second spin factor is included.
module hse_exchange
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  include 'fftw3.f03'
  public :: hse_kernel, hse_kernel_init, hse_kernel_apply, hse_kernel_destroy
  type hse_kernel
    integer :: n=0, mesh=0, ng=0, nk=0, block=0
    integer, allocatable :: order(:), point(:,:),shift(:,:)
    complex(c_double_complex), allocatable :: phase(:,:),work(:,:,:)
    real(c_double), allocatable :: kernel(:,:,:)
    type(c_ptr) :: forward=c_null_ptr, backward=c_null_ptr
  end type
contains
  subroutine hse_kernel_init(op,n,mesh,h,k,omega,block,ierr)
    type(hse_kernel),intent(inout) :: op
    integer,intent(in) :: n,mesh,block
    real(c_double),intent(in) :: h,omega,k(:,:)
    integer,intent(out) :: ierr
    integer :: ns,nk,ng,i,j,x,y,z,ix,iy,iz,index,idx(3),dims(3)
    real(c_double) :: pi,q2,q(3),scaled(3)
    complex(c_double_complex),allocatable :: spectrum(:,:,:)
    type(c_ptr) :: plan
    ierr=1
    call hse_kernel_destroy(op)
    if(n<1.or.mesh<1.or.block<1.or.h<=0.or.omega<=0) return
    if(.not.ieee_is_finite(h).or..not.ieee_is_finite(omega))return
    nk=mesh**3;ng=n**3;ns=n*mesh;pi=acos(-1d0)
    if(size(k,1)/=3.or.size(k,2)/=nk.or..not.all(ieee_is_finite(k)))return
    op%n=n;op%mesh=mesh;op%ng=ng;op%nk=nk;op%block=min(block,ng)
    allocate(op%order(nk),op%point(3,ng),op%shift(3,nk),op%phase(ng,nk),op%kernel(0:ns-1,0:ns-1,0:ns-1))
    op%order=0
    do i=1,nk
      scaled=(k(:,i)-k(:,1))*real(n*mesh,c_double)*h/(2*pi)
      if(maxval(abs(scaled-anint(scaled)))>1d-8)goto 900
      idx=modulo(nint(scaled),mesh);index=1+idx(1)+mesh*idx(2)+mesh**2*idx(3)
      if(op%order(index)/=0)goto 900
      op%order(index)=i
    enddo
    i=0
    do z=0,n-1;do y=0,n-1;do x=0,n-1
      i=i+1;op%point(:,i)=[x,y,z]
      do j=1,nk
        op%phase(i,j)=exp(cmplx(0d0,sum((k(:,j)-k(:,1))*op%point(:,i))*h,c_double))
      enddo
    enddo;enddo;enddo
    i=0
    do z=0,mesh-1;do y=0,mesh-1;do x=0,mesh-1
      i=i+1;op%shift(:,i)=[x,y,z]*n
    enddo;enddo;enddo
    allocate(spectrum(ns,ns,ns))
    do z=0,ns-1;do y=0,ns-1;do x=0,ns-1
      ix=x;if(x>=(ns+1)/2)ix=x-ns
      iy=y;if(y>=(ns+1)/2)iy=y-ns
      iz=z;if(z>=(ns+1)/2)iz=z-ns
      q=2*pi*real([ix,iy,iz],c_double)/(ns*h);q2=sum(q*q)
      if(q2<1d-24)then
        spectrum(x+1,y+1,z+1)=pi/omega**2
      else
        spectrum(x+1,y+1,z+1)=4*pi*(1-exp(-q2/(4*omega**2)))/q2
      endif
    enddo;enddo;enddo
    plan=fftw_plan_dft_3d(ns,ns,ns,spectrum,spectrum,FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    if(.not.c_associated(plan))goto 900
    call fftw_execute_dft(plan,spectrum,spectrum)
    call fftw_destroy_plan(plan)
    op%kernel=real(spectrum,c_double)/real(ns,c_double)**3
    allocate(op%work(op%block,ng,nk));dims=mesh
    op%forward=fftw_plan_many_dft(3,dims,op%block*ng,op%work,dims,op%block*ng,1, &
      op%work,dims,op%block*ng,1,FFTW_FORWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    op%backward=fftw_plan_many_dft(3,dims,op%block*ng,op%work,dims,op%block*ng,1, &
      op%work,dims,op%block*ng,1,FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    if(.not.c_associated(op%forward).or..not.c_associated(op%backward))goto 900
    ierr=0;return
900 call hse_kernel_destroy(op)
  end subroutine

  subroutine hse_kernel_apply(op,source,target,action,rank,nproc,ierr)
    type(hse_kernel),intent(inout) :: op
    complex(c_double_complex),intent(in) :: source(:,:,:),target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(in) :: rank,nproc
    integer,intent(out) :: ierr
    complex(c_double_complex),allocatable :: s(:,:,:),t(:,:,:)
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    integer :: lo,rows,ik,ki,b,j,r,offset(3),ng,nk,no,nt,ns
    external :: zgemm
    ierr=1;action=zero
    if(.not.c_associated(op%forward))return
    ng=op%ng;nk=op%nk;ns=op%n*op%mesh;b=op%block
    if(nproc<1.or.rank<0.or.rank>=nproc)return
    if(size(source,1)/=ng.or.size(target,1)/=ng.or.size(source,3)/=nk.or.size(target,3)/=nk)return
    if(any(shape(action)/=shape(target)))return
    no=size(source,2);nt=size(target,2)
    if(no<1.or.nt<1)return
    if(.not.all(ieee_is_finite(real(source))).or..not.all(ieee_is_finite(aimag(source))))return
    if(.not.all(ieee_is_finite(real(target))).or..not.all(ieee_is_finite(aimag(target))))return
    allocate(s(ng,no,nk),t(ng,nt,nk))
    do ik=1,nk
      do j=1,no;s(:,j,ik)=source(:,j,ik)*op%phase(:,ik);enddo
      do j=1,nt;t(:,j,ik)=target(:,j,ik)*op%phase(:,ik);enddo
    enddo
    do lo=1+rank*b,ng,nproc*b
      rows=min(b,ng-lo+1);op%work=zero
      do ki=1,nk
        ik=op%order(ki)
        call zgemm('N','C',rows,ng,no,one,s(lo,1,ik),ng,s(1,1,ik),ng,zero,op%work(1,1,ki),b)
      enddo
      call fftw_execute_dft(op%forward,op%work,op%work)
      do ki=1,nk;do j=1,ng;do r=1,rows
        offset=modulo(op%point(:,lo+r-1)-op%point(:,j)-op%shift(:,ki),ns)
        op%work(r,j,ki)=op%work(r,j,ki)*op%kernel(offset(1),offset(2),offset(3))
      enddo;enddo;enddo
      call fftw_execute_dft(op%backward,op%work,op%work)
      do ki=1,nk
        ik=op%order(ki)
        call zgemm('N','N',rows,nt,ng,-one/real(nk,c_double),op%work(1,1,ki),b, &
          t(1,1,ik),ng,zero,action(lo,1,ik),ng)
        do j=1,nt
          action(lo:lo+rows-1,j,ik)=action(lo:lo+rows-1,j,ik)*conjg(op%phase(lo:lo+rows-1,ik))
        enddo
      enddo
    enddo
    ierr=0
  end subroutine

  subroutine hse_kernel_destroy(op)
    type(hse_kernel),intent(inout) :: op
    if(c_associated(op%forward))call fftw_destroy_plan(op%forward)
    if(c_associated(op%backward))call fftw_destroy_plan(op%backward)
    op%forward=c_null_ptr;op%backward=c_null_ptr
    if(allocated(op%work))deallocate(op%work)
    if(allocated(op%order))deallocate(op%order,op%point,op%shift,op%phase,op%kernel)
    op%n=0;op%mesh=0;op%ng=0;op%nk=0;op%block=0
  end subroutine
end module
