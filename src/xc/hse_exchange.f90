! Full periodic sampled HSE kernel. Unit one-spin source occupations;
! neither the hybrid mixing fraction nor a second spin factor is included.
module hse_exchange
!$ use omp_lib, only: omp_get_num_threads
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  include 'fftw3.f03'
  public :: hse_kernel, hse_kernel_init, hse_kernel_apply, hse_kernel_destroy
  public :: hse_kernel_apply_distributed
  type hse_kernel
    integer :: n=0, mesh=0, ng=0, nk=0, block=0, phase_start=1, threads_used=1
    integer, allocatable :: order(:), point(:,:),shift(:,:)
    complex(c_double_complex), allocatable :: phase(:,:),work(:,:,:)
    real(c_double), allocatable :: kernel(:,:,:)
    type(c_ptr) :: forward=c_null_ptr, backward=c_null_ptr
  end type
contains
  subroutine hse_kernel_init(op,n,mesh,h,k,omega,block,ierr,phase_start,phase_count)
    type(hse_kernel),intent(inout) :: op
    integer,intent(in) :: n,mesh,block
    real(c_double),intent(in) :: h,omega,k(:,:)
    integer,intent(out) :: ierr
    integer,optional,intent(in) :: phase_start,phase_count
    integer :: ns,nk,ng,i,j,x,y,z,ix,iy,iz,index,idx(3),dims(3),first,nphase
    real(c_double) :: pi,q2,q(3),scaled(3)
    complex(c_double_complex),allocatable :: spectrum(:,:,:)
    type(c_ptr) :: plan
    ierr=1
    call hse_kernel_destroy(op)
    if(n<1.or.mesh<1.or.block<1.or.h<=0.or.omega<=0) return
    if(.not.ieee_is_finite(h).or..not.ieee_is_finite(omega))return
    nk=mesh**3;ng=n**3;ns=n*mesh;pi=acos(-1d0)
    if(size(k,1)/=3.or.size(k,2)/=nk.or..not.all(ieee_is_finite(k)))return
    first=1;nphase=nk
    if(present(phase_start))first=phase_start
    if(present(phase_count))nphase=phase_count
    if(first<1.or.nphase<0.or.first+nphase-1>nk)return
    op%phase_start=first
    op%n=n;op%mesh=mesh;op%ng=ng;op%nk=nk;op%block=min(block,ng)
    allocate(op%order(nk),op%point(3,ng),op%shift(3,nk),op%phase(ng,nphase),op%kernel(0:ns-1,0:ns-1,0:ns-1))
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
      do j=1,nphase
        op%phase(i,j)=exp(cmplx(0d0,sum((k(:,first+j-1)-k(:,1))*op%point(:,i))*h,c_double))
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
    if(op%phase_start/=1.or.size(op%phase,2)/=op%nk)return
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

  ! K-distributed source/target/action; transpose density tiles, never orbitals.
  ! Caller supplies identical layout/kernel metadata and communicator size on all ranks.
  subroutine hse_kernel_apply_distributed(op,source,target,action,starts,counts,rank,transpose_tiles,ierr)
    type(hse_kernel),intent(inout) :: op
    complex(c_double_complex),intent(in) :: source(:,:,:),target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(in) :: starts(:),counts(:),rank
    integer,intent(out) :: ierr
    interface
      subroutine transpose_tiles(send,recv,count)
        import c_double_complex
        complex(c_double_complex),intent(in) :: send(:)
        complex(c_double_complex),intent(out) :: recv(:)
        integer,intent(in) :: count
      end subroutine
    end interface
    complex(c_double_complex),allocatable :: s(:,:,:),t(:,:,:)
    complex(c_double_complex),allocatable,target :: send(:),recv(:)
    complex(c_double_complex),pointer :: sb(:,:,:,:),rb(:,:,:,:)
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    integer :: ng,nk,np,nlocal,no,nt,b,km,nmsg,p,j,ki,ik,base,lo,rows,r,g,offset(3),ns
    integer,allocatable :: inverse(:)
    complex(c_double_complex) :: valid_send(size(counts)),valid_recv(size(counts))
    logical :: valid
    external :: zgemm
    ierr=1;action=zero
    np=size(counts);ng=op%ng;nk=op%nk;b=op%block;ns=op%n*op%mesh
    if(np<1.or.rank<0.or.rank>=np.or.size(starts)/=np)return
    ! All communicator members must participate even if one has invalid data.
    valid=c_associated(op%forward)
    valid=valid.and..not.any(counts<0).and.sum(counts)==nk.and.starts(1)==1
    do p=2,np
      valid=valid.and.starts(p)==starts(p-1)+counts(p-1)
    enddo
    nlocal=counts(rank+1);no=size(source,2);nt=size(target,2)
    valid=valid.and.min(no,nt)>0.and.size(source,1)==ng.and.size(target,1)==ng
    valid=valid.and.size(source,3)==nlocal.and.size(target,3)==nlocal
    valid=valid.and.all(shape(action)==shape(target))
    if(allocated(op%phase))then
      valid=valid.and.op%phase_start==starts(rank+1).and.size(op%phase,2)==nlocal
    else
      valid=.false.
    endif
    valid=valid.and.all(ieee_is_finite(real(source))).and.all(ieee_is_finite(aimag(source)))
    valid=valid.and.all(ieee_is_finite(real(target))).and.all(ieee_is_finite(aimag(target)))
    valid_send=zero
    if(.not.valid)valid_send=one
    call transpose_tiles(valid_send,valid_recv,1)
    if(any(valid_recv/=zero))return
    km=maxval(counts);nmsg=b*ng*km
    allocate(send(nmsg*np),recv(nmsg*np),s(ng,no,nlocal),t(ng,nt,nlocal),inverse(nk))
    sb(1:b,1:ng,1:km,1:np)=>send;rb(1:b,1:ng,1:km,1:np)=>recv
    do ki=1,nk;inverse(op%order(ki))=ki;enddo
    op%threads_used=1
    !$omp parallel private(j,g)
    !$omp single
!$  op%threads_used=omp_get_num_threads()
    !$omp end single
    !$omp do schedule(static)
    do j=1,nlocal
      do g=1,no;s(:,g,j)=source(:,g,j)*op%phase(:,j);enddo
      do g=1,nt;t(:,g,j)=target(:,g,j)*op%phase(:,j);enddo
    enddo
    !$omp end do
    !$omp end parallel
    do base=1,ng,np*b
      send=zero
      ! BLAS owns threading here; call only outside application OpenMP regions.
      do p=0,np-1
        lo=base+p*b;rows=min(b,ng-lo+1)
        if(rows<=0)cycle
        do j=1,nlocal
          call zgemm('N','C',rows,ng,no,one,s(lo,1,j),ng,s(1,1,j),ng,zero,sb(1,1,j,p+1),b)
        enddo
      enddo
      call transpose_tiles(send,recv,nmsg)
      op%work=zero
      !$omp parallel do private(j,ki) schedule(static)
      do p=1,np;do j=1,counts(p)
        ki=inverse(starts(p)+j-1);op%work(:,:,ki)=rb(:,:,j,p)
      enddo;enddo
      !$omp end parallel do
      call fftw_execute_dft(op%forward,op%work,op%work)
      lo=base+rank*b;rows=min(b,ng-lo+1)
      !$omp parallel do collapse(2) private(r,offset) schedule(static)
      do ki=1,nk;do g=1,ng;do r=1,max(0,rows)
        offset=modulo(op%point(:,lo+r-1)-op%point(:,g)-op%shift(:,ki),ns)
        op%work(r,g,ki)=op%work(r,g,ki)*op%kernel(offset(1),offset(2),offset(3))
      enddo;enddo;enddo
      !$omp end parallel do
      call fftw_execute_dft(op%backward,op%work,op%work)
      send=zero
      !$omp parallel do private(j,ki) schedule(static)
      do p=1,np;do j=1,counts(p)
        ki=inverse(starts(p)+j-1);sb(:,:,j,p)=op%work(:,:,ki)
      enddo;enddo
      !$omp end parallel do
      call transpose_tiles(send,recv,nmsg)
      ! BLAS owns threading here; call only outside application OpenMP regions.
      do p=0,np-1
        lo=base+p*b;rows=min(b,ng-lo+1)
        if(rows<=0)cycle
        do j=1,nlocal
          call zgemm('N','N',rows,nt,ng,-one/real(nk,c_double),rb(1,1,j,p+1),b, &
            t(1,1,j),ng,zero,action(lo,1,j),ng)
          do g=1,nt
            action(lo:lo+rows-1,g,j)=action(lo:lo+rows-1,g,j)*conjg(op%phase(lo:lo+rows-1,j))
          enddo
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
