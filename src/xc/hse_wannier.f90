! Rectangular fragment-periodic screened exchange in a Wannier representation.
! Q = Psi sqrt(f/2) U preserves fractional occupations; Phi = Psi U is the
! orthonormal localization frame. Neither the HSE fraction nor spin doubling
! is included in this module's action. Full periodic support is retained.
module hse_wannier
  use iso_c_binding
  use iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use hse_wannier_gauge, only: gauge_transport,gauge_minimize,gauge_seed
  !$ use omp_lib, only: omp_get_max_threads,omp_get_thread_num
  implicit none
  private
  include 'fftw3.f03'
  public :: wannier_snapshot,wannier_refresh_source
  public :: s_hse_wannier,wannier_init,wannier_destroy,wannier_localize
  public :: wannier_set_source,wannier_apply,wannier_forward,wannier_backward
  type s_hse_wannier
    integer :: n(3)=0,mesh(3)=0,ns(3)=0,ng=0,nk=0,ngs=0,updates=0
    integer(int64) :: fft_pairs_total=0_int64,fft_pairs_executed=0_int64,fft_batches_executed=0_int64
    integer :: fft_batch_size=1,worker_batch=0
    integer :: localization_iterations=0,localization_status=1,workers=0
    real(8) :: h(3)=0d0,dv=0d0,spread=0d0,gradient=huge(1d0),min_singular=0d0
    real(8),allocatable :: k(:,:),position(:,:),multiplier(:,:,:),source_occupation(:,:)
    integer,allocatable :: source_indices(:)
    integer,allocatable :: point(:,:),primitive_point(:),cell(:,:),neighbors(:,:)
    complex(8),allocatable :: phase(:,:),source(:,:),previous(:,:,:),gauge(:,:,:),work(:,:,:)
    complex(8),allocatable :: worker_work(:,:,:,:,:)
    type(c_ptr),allocatable :: worker_forward(:,:),worker_backward(:,:)
    type(c_ptr) :: forward=c_null_ptr,backward=c_null_ptr
  end type
contains
  subroutine wannier_snapshot(op,occupation,omega,exchange,residual,iteration,converged,path,status)
    use iso_fortran_env, only: int32
    implicit none
    type(s_hse_wannier),intent(in) :: op
    real(8),intent(in) :: occupation(:,:),omega,exchange,residual
    integer,intent(in) :: iteration
    logical,intent(in) :: converged
    character(*),intent(in) :: path
    integer,intent(out) :: status
    integer :: iu,close_status
    status=1
    if(.not.allocated(op%source).or..not.allocated(op%previous).or..not.allocated(op%gauge))return
    if(any(shape(occupation)/=[size(op%gauge,1),op%nk]))return
    open(newunit=iu,file=path,status='replace',access='stream',form='unformatted',iostat=status)
    if(status/=0)return
    write(iu,iostat=status)int([16909060,1,op%n,op%mesh,size(op%gauge,1),op%updates, &
      op%localization_iterations,op%localization_status,iteration,merge(1,0,converged)],int32)
    if(status==0)write(iu,iostat=status)op%h,omega,op%spread,op%gradient,op%min_singular,exchange,residual
    if(status==0)write(iu,iostat=status)occupation,op%gauge,op%previous,op%source
    close(iu,iostat=close_status)
    if(status==0)status=close_status
  end subroutine
  subroutine wannier_destroy(op)
    implicit none
    type(s_hse_wannier),intent(inout) :: op
    type(s_hse_wannier) :: empty
    call clear_workers(op)
    if(c_associated(op%forward))call fftw_destroy_plan(op%forward)
    if(c_associated(op%backward))call fftw_destroy_plan(op%backward)
    op=empty
  end subroutine

  subroutine wannier_init(op,n,mesh,h,k,omega,status)
    implicit none
    type(s_hse_wannier),intent(inout) :: op
    integer,intent(in) :: n(3),mesh(3)
    real(8),intent(in) :: h(3),k(:,:),omega
    integer,intent(out) :: status
    integer :: x,y,z,g,ik,j,axis,offset(3),index(3),flat,ns(3),p(3),r(3),ic
    integer,allocatable :: order(:)
    real(8) :: pi,scaled(3),length(3),q(3),q2,delta(3)
    call wannier_destroy(op)
    status=1
    if(any(n<1).or.any(mesh<1).or.any(h<=0d0).or.omega<=0d0)return
    if(.not.all(ieee_is_finite(h)).or..not.ieee_is_finite(omega))return
    if(size(k,1)/=3.or.size(k,2)/=product(mesh).or..not.all(ieee_is_finite(k)))return
    op%n=n;op%mesh=mesh;op%ns=n*mesh;op%ng=product(n);op%nk=product(mesh);op%ngs=product(op%ns)
    op%h=h;op%dv=product(h);op%k=k;ns=op%ns;pi=acos(-1d0);length=n*h
    allocate(order(op%nk),op%neighbors(6,op%nk));order=0
    do ik=1,op%nk
      scaled=(k(:,ik)-k(:,1))*length*mesh/(2*pi)
      if(maxval(abs(scaled-anint(scaled)))>1d-8)goto 900
      index=modulo(nint(scaled),mesh)
      flat=1+index(1)+mesh(1)*(index(2)+mesh(2)*index(3))
      if(order(flat)/=0)goto 900
      order(flat)=ik
    enddo
    do ik=1,op%nk
      index=modulo(nint((k(:,ik)-k(:,1))*length*mesh/(2*pi)),mesh)
      do axis=1,3
        offset=0;offset(axis)=1
        p=modulo(index+offset,mesh)
        op%neighbors(axis,ik)=order(1+p(1)+mesh(1)*(p(2)+mesh(2)*p(3)))
        p=modulo(index-offset,mesh)
        op%neighbors(axis+3,ik)=order(1+p(1)+mesh(1)*(p(2)+mesh(2)*p(3)))
      enddo
    enddo
    allocate(op%point(3,op%ngs),op%primitive_point(op%ngs),op%position(3,op%ng))
    allocate(op%cell(3,op%nk),op%phase(op%ngs,op%nk))
    g=0
    do z=0,n(3)-1;do y=0,n(2)-1;do x=0,n(1)-1
      g=g+1;op%position(:,g)=[x,y,z]*h
    enddo;enddo;enddo
    g=0
    do z=0,ns(3)-1;do y=0,ns(2)-1;do x=0,ns(1)-1
      g=g+1;op%point(:,g)=[x,y,z];p=modulo([x,y,z],n)
      op%primitive_point(g)=1+p(1)+n(1)*(p(2)+n(2)*p(3))
      do ik=1,op%nk
        op%phase(g,ik)=exp(cmplx(0d0,sum((k(:,ik)-k(:,1))*[x,y,z]*h),8))
      enddo
    enddo;enddo;enddo
    ic=0
    do z=0,mesh(3)-1;do y=0,mesh(2)-1;do x=0,mesh(1)-1
      ic=ic+1;op%cell(:,ic)=[x,y,z]*n
    enddo;enddo;enddo
    allocate(op%multiplier(ns(1),ns(2),ns(3)),op%work(ns(1),ns(2),ns(3)))
    do z=0,ns(3)-1;do y=0,ns(2)-1;do x=0,ns(1)-1
      p=[x,y,z]
      where(p>=(ns+1)/2)p=p-ns
      q=2*pi*p/(ns*h);q2=sum(q*q)
      if(q2<1d-24)then
        op%multiplier(x+1,y+1,z+1)=pi/omega**2
      else
        op%multiplier(x+1,y+1,z+1)=4*pi*(1-exp(-q2/(4*omega**2)))/q2
      endif
    enddo;enddo;enddo
    op%forward=fftw_plan_dft_3d(ns(3),ns(2),ns(1),op%work,op%work,FFTW_FORWARD,FFTW_ESTIMATE)
    op%backward=fftw_plan_dft_3d(ns(3),ns(2),ns(1),op%work,op%work,FFTW_BACKWARD,FFTW_ESTIMATE)
    if(.not.c_associated(op%forward).or..not.c_associated(op%backward))goto 900
    status=0;return
900 call wannier_destroy(op)
  end subroutine

  subroutine wannier_forward(op,bloch,home)
    implicit none
    type(s_hse_wannier),intent(in) :: op
    complex(8),intent(in) :: bloch(:,:,:)
    complex(8),intent(out) :: home(:,:)
    integer :: ik,j,g,p
    home=0d0
    do ik=1,op%nk
      do j=1,size(bloch,2)
        do g=1,op%ngs
          p=op%primitive_point(g)
          home(g,j)=home(g,j)+bloch(p,j,ik)*op%phase(g,ik)/op%nk
        enddo
      enddo
    enddo
  end subroutine

  subroutine wannier_backward(op,home,bloch)
    implicit none
    type(s_hse_wannier),intent(in) :: op
    complex(8),intent(in) :: home(:,:)
    complex(8),intent(out) :: bloch(:,:,:)
    integer :: ik,j,g,p
    bloch=0d0
    do ik=1,op%nk
      do j=1,size(home,2)
        do g=1,op%ngs
          p=op%primitive_point(g)
          bloch(p,j,ik)=bloch(p,j,ik)+home(g,j)*conjg(op%phase(g,ik))
        enddo
      enddo
    enddo
  end subroutine

  subroutine wannier_refresh_source(op,psi,occupation,maxiter,tolerance,status)
    type(s_hse_wannier),intent(inout) :: op
    complex(8),intent(in) :: psi(:,:,:)
    real(8),intent(in) :: occupation(:,:),tolerance
    integer,intent(in) :: maxiter
    integer,intent(out) :: status
    integer,allocatable :: indices(:)
    integer :: j
    logical :: reset
    status=1
    if(size(psi,2)<1.or.size(psi,1)/=op%ng.or.size(psi,3)/=op%nk)return
    if(any(shape(occupation)/=[size(psi,2),op%nk]))return
    if(any(occupation<0d0).or.any(occupation>2d0).or..not.all(ieee_is_finite(occupation)))return
    ! Keep a common band set over k; never discard a positive occupation.
    indices=pack([(j,j=1,size(psi,2))],any(occupation>0d0,dim=2))
    if(size(indices)==0)indices=[1] ! zero-density operator still has a valid frame
    reset=.true.
    if(allocated(op%source_indices))then
      if(size(op%source_indices)==size(indices))reset=any(op%source_indices/=indices)
    endif
    if(reset)then
      if(allocated(op%gauge))deallocate(op%gauge)
      if(allocated(op%previous))deallocate(op%previous)
      op%min_singular=0d0
    endif
    op%source_indices=indices
    call wannier_localize(op,psi(:,indices,:),maxiter,tolerance,status)
    if(status==0)call wannier_set_source(op,psi(:,indices,:),occupation(indices,:),op%gauge,status)
  end subroutine

  subroutine wannier_localize(op,psi,maxiter,tolerance,status)
    implicit none
    type(s_hse_wannier),intent(inout) :: op
    complex(8),intent(in) :: psi(:,:,:)
    integer,intent(in) :: maxiter
    real(8),intent(in) :: tolerance
    integer,intent(out) :: status
    complex(8),allocatable :: raw(:,:,:,:),shifted(:,:)
    real(8) :: b(3,6),weights(6),length(3),gvec(3),pi,delta
    integer :: n,ik,axis,next,j,g,transport_status
    status=1;n=size(psi,2)
    if(size(psi,1)/=op%ng.or.size(psi,3)/=op%nk)return
    if(allocated(op%gauge))then
      if(size(op%gauge,1)/=n)then
        deallocate(op%gauge)
        if(allocated(op%previous))deallocate(op%previous)
        op%min_singular=0d0
      endif
    endif
    if(.not.allocated(op%gauge))then
      allocate(op%gauge(n,n,op%nk));op%gauge=0d0
      do ik=1,op%nk
        do j=1,n
          op%gauge(j,j,ik)=1d0
        enddo
      enddo
    endif
    if(.not.allocated(op%previous))then
      call gauge_seed(psi,op%position,op%k,op%gauge,transport_status)
      if(transport_status/=0)then
        op%gauge=0d0
        do ik=1,op%nk
          do j=1,n
            op%gauge(j,j,ik)=1d0
          enddo
        enddo
      endif
    endif
    if(allocated(op%previous))then
      call gauge_transport(psi,op%previous,op%dv,op%gauge,op%min_singular,transport_status)
      if(transport_status/=0)then
        ! Overlap loss: restart the gauge, never project away current states.
        op%gauge=0d0
        do ik=1,op%nk
          do j=1,n
            op%gauge(j,j,ik)=1d0
          enddo
        enddo
      endif
    endif
    if(maxiter==0)then
      ! Polar transport already supplies U. Avoid six grid-by-band overlap
      ! products when no spread minimization is requested for this refresh.
      op%localization_iterations=0;op%localization_status=2
      op%spread=-1d0;op%gradient=-1d0 ! not evaluated for the current source
      if(.not.allocated(op%previous))allocate(op%previous(op%ng,n,op%nk))
      do ik=1,op%nk
        op%previous(:,:,ik)=matmul(psi(:,:,ik),op%gauge(:,:,ik))
      enddo
      op%updates=op%updates+1
      status=0;return
    endif
    allocate(raw(n,n,6,op%nk),shifted(op%ng,n))
    b=0d0;pi=acos(-1d0);length=op%n*op%h
    do axis=1,3
      delta=2*pi/(length(axis)*op%mesh(axis))
      b(axis,axis)=delta;b(axis,axis+3)=-delta
      weights(axis)=1d0/(2*delta**2);weights(axis+3)=weights(axis)
    enddo
    do ik=1,op%nk
      do axis=1,3
        next=op%neighbors(axis,ik)
        gvec=op%k(:,ik)+b(:,axis)-op%k(:,next)
        do j=1,n
          do g=1,op%ng
            shifted(g,j)=psi(g,j,next)*exp(cmplx(0d0,-sum(op%position(:,g)*gvec),8))
          enddo
        enddo
        raw(:,:,axis,ik)=matmul(conjg(transpose(psi(:,:,ik))),shifted)*op%dv
        raw(:,:,axis+3,next)=conjg(transpose(raw(:,:,axis,ik)))
      enddo
    enddo
    call gauge_minimize(op%gauge,raw,op%neighbors,b,weights,maxiter,tolerance,op%spread,op%gradient, &
                        op%localization_iterations,op%localization_status)
    if(.not.allocated(op%previous))allocate(op%previous(op%ng,n,op%nk))
    do ik=1,op%nk
      op%previous(:,:,ik)=matmul(psi(:,:,ik),op%gauge(:,:,ik))
    enddo
    op%updates=op%updates+1
    ! A valid but not converged gauge changes cost/locality, not full-support EXX.
    status=0
  end subroutine

  subroutine wannier_set_source(op,psi,occupation,gauge,status)
    implicit none
    type(s_hse_wannier),intent(inout) :: op
    complex(8),intent(in) :: psi(:,:,:),gauge(:,:,:)
    real(8),intent(in) :: occupation(:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: weighted(:,:,:),tmp(:,:)
    integer :: n,j,ik
    status=1;n=size(psi,2)
    if(size(psi,1)/=op%ng.or.size(psi,3)/=op%nk)return
    if(any(shape(occupation)/=[n,op%nk]).or.any(shape(gauge)/=[n,n,op%nk]))return
    if(any(occupation<0d0).or.any(occupation>2d0).or..not.all(ieee_is_finite(occupation)))return
    allocate(weighted(op%ng,n,op%nk),tmp(op%ng,n))
    do ik=1,op%nk
      do j=1,n
        tmp(:,j)=psi(:,j,ik)*sqrt(occupation(j,ik)/2d0)
      enddo
      weighted(:,:,ik)=matmul(tmp,gauge(:,:,ik))
    enddo
    if(allocated(op%source))deallocate(op%source)
    allocate(op%source(op%ngs,n))
    call wannier_forward(op,weighted,op%source)
    op%source_occupation=occupation
    status=0
  end subroutine

  subroutine clear_workers(op)
    type(s_hse_wannier),intent(inout) :: op
    integer :: t,b
    if(allocated(op%worker_forward))then
      do t=1,size(op%worker_forward,2);do b=1,size(op%worker_forward,1)
        if(c_associated(op%worker_forward(b,t)))call fftw_destroy_plan(op%worker_forward(b,t))
        if(c_associated(op%worker_backward(b,t)))call fftw_destroy_plan(op%worker_backward(b,t))
      enddo;enddo
      deallocate(op%worker_forward,op%worker_backward,op%worker_work)
    endif
    op%workers=0;op%worker_batch=0
  end subroutine

  subroutine prepare_workers(op,count,batch,status)
    type(s_hse_wannier),intent(inout) :: op
    integer,intent(in) :: count,batch
    integer,intent(out) :: status
    integer :: t,b,dims(3)
    status=0
    if(op%workers==count.and.op%worker_batch==batch)return
    call clear_workers(op)
    allocate(op%worker_work(op%ns(1),op%ns(2),op%ns(3),batch,count))
    allocate(op%worker_forward(batch,count),op%worker_backward(batch,count))
    op%worker_forward=c_null_ptr;op%worker_backward=c_null_ptr
    dims=op%ns(3:1:-1)
    ! Serial planning, one buffer per worker. Cache every possible compact tail
    ! so a partly filled tile never performs padded zero-density FFTs.
    do t=1,count;do b=1,batch
      op%worker_forward(b,t)=fftw_plan_many_dft(3,dims,b,op%worker_work(:,:,:,:,t),dims,1,op%ngs, &
        op%worker_work(:,:,:,:,t),dims,1,op%ngs,FFTW_FORWARD,FFTW_ESTIMATE)
      op%worker_backward(b,t)=fftw_plan_many_dft(3,dims,b,op%worker_work(:,:,:,:,t),dims,1,op%ngs, &
        op%worker_work(:,:,:,:,t),dims,1,op%ngs,FFTW_BACKWARD,FFTW_ESTIMATE)
      if(.not.c_associated(op%worker_forward(b,t)).or..not.c_associated(op%worker_backward(b,t)))then
        status=1;call clear_workers(op);return
      endif
    enddo;enddo
    op%workers=count;op%worker_batch=batch
  end subroutine

  subroutine wannier_apply(op,target,action,status)
    implicit none
    type(s_hse_wannier),intent(inout) :: op
    complex(8),intent(in) :: target(:,:,:)
    complex(8),intent(out) :: action(:,:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: home(:,:),result(:,:),source(:)
    integer :: i,j,ic,g,p(3),index,nt,t,nworkers,batch,lo,nb,k
    integer :: columns(32)
    integer(int64) :: executed,batches
    status=1;action=0d0
    op%fft_pairs_total=0_int64;op%fft_pairs_executed=0_int64;op%fft_batches_executed=0_int64
    if(.not.allocated(op%source))return
    if(size(target,1)/=op%ng.or.size(target,3)/=op%nk.or.any(shape(action)/=shape(target)))return
    nt=size(target,2)
    if(nt<1)return
    op%fft_pairs_total=int(op%nk,int64)*int(size(op%source,2),int64)*int(nt,int64)
    executed=0_int64;batches=0_int64
    if(op%fft_batch_size<1.or.op%fft_batch_size>size(columns))return
    allocate(home(op%ngs,nt),result(op%ngs,nt),source(op%ngs))
    nworkers=1
    !$ nworkers=min(nt,omp_get_max_threads())
    batch=min(op%fft_batch_size,max(1,nt/nworkers))
    call prepare_workers(op,nworkers,batch,status)
    if(status/=0)return
    call wannier_forward(op,target,home);result=0d0
    do ic=1,op%nk
      do i=1,size(op%source,2)
        do g=1,op%ngs
          p=modulo(op%point(:,g)-op%cell(:,ic),op%ns)
          index=1+p(1)+op%ns(1)*(p(2)+op%ns(2)*p(3))
          source(g)=op%source(index,i)
        enddo
        if(all(source==(0d0,0d0)))cycle
        !$omp parallel do default(none) schedule(static) num_threads(nworkers) &
        !$omp shared(op,home,result,source,nt,batch) private(lo,j,t,nb,k,columns) reduction(+:executed,batches)
        do lo=1,nt,batch
          t=1
          !$ t=omp_get_thread_num()+1
          nb=0
          do j=lo,min(nt,lo+batch-1)
            op%worker_work(:,:,:,nb+1,t)=reshape(conjg(source)*home(:,j),op%ns)
            ! Compact only exact nonzero pair densities; no magnitude threshold.
            if(all(op%worker_work(:,:,:,nb+1,t)==(0d0,0d0)))cycle
            nb=nb+1;columns(nb)=j
          enddo
          if(nb==0)cycle
          executed=executed+int(nb,int64);batches=batches+1_int64
          call fftw_execute_dft(op%worker_forward(nb,t),op%worker_work(:,:,:,:,t),op%worker_work(:,:,:,:,t))
          do k=1,nb
            op%worker_work(:,:,:,k,t)=op%worker_work(:,:,:,k,t)*op%multiplier
          enddo
          call fftw_execute_dft(op%worker_backward(nb,t),op%worker_work(:,:,:,:,t),op%worker_work(:,:,:,:,t))
          do k=1,nb
            j=columns(k)
            result(:,j)=result(:,j)-source*reshape(op%worker_work(:,:,:,k,t),[op%ngs])/op%ngs
          enddo
        enddo
        !$omp end parallel do
      enddo
    enddo
    op%fft_pairs_executed=executed;op%fft_batches_executed=batches
    call wannier_backward(op,result,action)
    status=0
  end subroutine
end module
