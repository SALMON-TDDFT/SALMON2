! Full periodic sampled hybrid exchange kernel. Unit one-spin source occupations;
! neither the hybrid mixing fraction nor a second spin factor is included.
module exx_k_exchange
!$ use omp_lib, only: omp_get_num_threads,omp_get_max_threads
  use exx_k_backend, only: s_exx_k_backend,k_backend_factory
  use iso_fortran_env, only: int64
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  include 'fftw3.f03'
  public :: exx_k_kernel, exx_k_kernel_init, exx_k_kernel_apply, exx_k_kernel_destroy
  public :: exx_k_kernel_apply_distributed
  type exx_k_kernel
    class(s_exx_k_backend),allocatable :: accelerator
    integer :: n(3)=0, mesh(3)=0, ng=0, nk=0, block=0, phase_start=1, threads_used=1
    integer :: fft_nk=0,cosets=1
    logical :: profile=.false.,contiguous_fft=.false.,auto_fft=.true.
    real(c_double) :: seconds(6)=0d0,fft_trial_seconds(2)=0d0
    integer, allocatable :: order(:), point(:,:),shift(:,:),distance_index(:,:,:,:)
    complex(c_double_complex), allocatable :: phase(:,:),work(:,:,:)
    real(c_double), allocatable :: kernel(:,:,:)
    type(c_ptr) :: forward=c_null_ptr, backward=c_null_ptr
  end type
interface exx_k_kernel_init
    module procedure init_cubic,init_rectangular
  end interface
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_1d,finite_real_2d,finite_real_3d
  interface salmon_all_finite
    module procedure finite_real_1d,finite_real_2d,finite_real_3d
  end interface
contains
  subroutine init_cubic(op,n,mesh,h,k,omega,block,ierr,phase_start,phase_count, &
                        block_rows,profile,fft_layout,create_backend,coulomb_radius,source_reduction)
    implicit none
    type(exx_k_kernel),intent(inout) :: op
    integer,intent(in) :: n,mesh,block
    real(c_double),intent(in) :: h,k(:,:),omega
    integer,intent(out) :: ierr
    integer,optional,intent(in) :: phase_start,phase_count,block_rows
    logical,optional,intent(in) :: profile
    character(*),optional,intent(in) :: fft_layout
    procedure(k_backend_factory),pointer,optional :: create_backend
    integer,optional,intent(in) :: source_reduction(3)
    real(c_double),optional,intent(in) :: coulomb_radius
    call init_rectangular(op,[n,n,n],[mesh,mesh,mesh],[h,h,h],k,omega,block,ierr, &
      phase_start,phase_count,block_rows,profile,fft_layout,create_backend,coulomb_radius,source_reduction)
  end subroutine

  subroutine init_rectangular(op,n,mesh,h,k,omega,block,ierr,phase_start,phase_count, &
                             block_rows,profile,fft_layout,create_backend,coulomb_radius,source_reduction)
    implicit none
    type(exx_k_kernel),intent(inout) :: op
    integer,intent(in) :: n(3),mesh(3),block
    real(c_double),intent(in) :: h(3),omega,k(:,:)
    integer,intent(out) :: ierr
    integer,optional,intent(in) :: phase_start,phase_count
    integer,optional,intent(in) :: block_rows
    logical,optional,intent(in) :: profile
    character(*),optional,intent(in) :: fft_layout
    procedure(k_backend_factory),pointer,optional :: create_backend
    integer,optional,intent(in) :: source_reduction(3)
    real(c_double),optional,intent(in) :: coulomb_radius
    integer :: ns(3),nk,ng,i,j,x,y,z,ix,iy,iz,index,idx(3),dims(3),first,nphase,chosen,axis
    integer :: reduction(3),coarse(3),coarse_ns(3),residue(3),cell(3),ax,ay,az
    real(c_double) :: pi,q2,q(3),scaled(3),radius
    complex(c_double_complex),allocatable :: spectrum(:,:,:)
    real(c_double),allocatable :: folded_kernel(:,:,:)
    type(c_ptr) :: plan
    ierr=1
    call exx_k_kernel_destroy(op,ierr)
    if(ierr/=0)return
    ierr=1
    if(any(n<1).or.any(mesh<1).or.block<1.or.any(h<=0).or.omega<0) return
    if(.not.salmon_all_finite(h).or..not.ieee_is_finite(omega))return
    radius=.5d0*minval(n*mesh*h)
    if(present(coulomb_radius))then
      if(.not.ieee_is_finite(coulomb_radius).or.coulomb_radius<0d0)return
      if(coulomb_radius>0d0)radius=coulomb_radius
    endif
    if(omega==0d0.and.radius>.5d0*minval(n*mesh*h)*(1d0+1d-12))return
    reduction=1
    if(present(source_reduction))reduction=source_reduction
    if(any(reduction<1))return
    if(any(modulo(mesh,reduction)/=0))return
    if(any(reduction>1).and.omega<=0d0)return
    if(any(reduction>1).and.present(create_backend))then
      if(associated(create_backend))return
    endif
    coarse=mesh/reduction
    coarse_ns=n*coarse
    op%fft_nk=product(coarse);op%cosets=product(reduction)
    chosen=block
    if(present(block_rows))then
      if(block_rows<0.or.block_rows>product(n))return
      if(block_rows>0)chosen=block_rows
    endif
    op%profile=.false.
    if(present(profile))op%profile=profile
    op%seconds=0d0
    op%contiguous_fft=.false.;op%auto_fft=.true.;op%fft_trial_seconds=0d0
    if(present(fft_layout))then
      select case(trim(fft_layout))
      case('auto')
      case('strided')
        op%auto_fft=.false.
      case('contiguous')
        op%auto_fft=.false.
        op%contiguous_fft=.true.
      case default
        return
      end select
    endif
    if(op%cosets>1)then
      if(op%contiguous_fft)return
      op%auto_fft=.false.
    endif
    nk=product(mesh);ng=product(n);ns=n*mesh;pi=acos(-1d0)
    if(size(k,1)/=3.or.size(k,2)/=nk.or..not.salmon_all_finite(k))return
    first=1;nphase=nk
    if(present(phase_start))first=phase_start
    if(present(phase_count))nphase=phase_count
    if(first<1.or.nphase<0.or.first+nphase-1>nk)return
    op%phase_start=first
    op%n=n;op%mesh=coarse;op%ng=ng;op%nk=nk;op%block=min(chosen,ng)
    allocate(op%order(nk),op%point(3,ng),op%shift(3,nk),op%phase(ng,nphase),op%kernel(0:ns(1)-1,0:ns(2)-1,0:ns(3)-1))
    ! Separable periodic coordinate differences: O(n*n*mesh) integers, not
    ! O(ng*ng*nk). Avoid integer modulo in every exchange-kernel multiply.
    allocate(op%distance_index(0:maxval(n)-1,0:maxval(n)-1,0:maxval(mesh)-1,3))
    do axis=1,3
      do z=0,coarse(axis)-1;do y=0,n(axis)-1;do x=0,n(axis)-1
        op%distance_index(x,y,z,axis)=modulo(x-y-z*n(axis),coarse_ns(axis))
      enddo;enddo;enddo
    enddo
    op%order=0
    do i=1,nk
      scaled=(k(:,i)-k(:,1))*real(n*mesh,c_double)*h/(2*pi)
      if(maxval(abs(scaled-anint(scaled)))>1d-8)goto 900
      idx=modulo(nint(scaled),mesh)
      cell=idx/reduction;residue=modulo(idx,reduction)
      index=1+cell(1)+coarse(1)*(cell(2)+coarse(2)*cell(3))+ &
        op%fft_nk*(residue(1)+reduction(1)*(residue(2)+reduction(2)*residue(3)))
      if(op%order(index)/=0)goto 900
      op%order(index)=i
    enddo
    i=0
    do z=0,n(3)-1;do y=0,n(2)-1;do x=0,n(1)-1
      i=i+1;op%point(:,i)=[x,y,z]
      do j=1,nphase
        op%phase(i,j)=exp(cmplx(0d0,sum((k(:,first+j-1)-k(:,1))*op%point(:,i)*h),c_double))
      enddo
    enddo;enddo;enddo
    do i=1,nk
      index=modulo(i-1,op%fft_nk)
      op%shift(:,i)=[modulo(index,coarse(1)),modulo(index/coarse(1),coarse(2)), &
        index/(coarse(1)*coarse(2))]*n
    enddo
    allocate(spectrum(ns(1),ns(2),ns(3)))
    do z=0,ns(3)-1;do y=0,ns(2)-1;do x=0,ns(1)-1
      ix=x;if(x>=(ns(1)+1)/2)ix=x-ns(1)
      iy=y;if(y>=(ns(2)+1)/2)iy=y-ns(2)
      iz=z;if(z>=(ns(3)+1)/2)iz=z-ns(3)
      q=2*pi*real([ix,iy,iz],c_double)/(ns*h);q2=sum(q*q)
      if(omega==0d0)then
        if(q2<1d-24)then
          spectrum(x+1,y+1,z+1)=2*pi*radius**2
        else
          spectrum(x+1,y+1,z+1)=8*pi*sin(.5d0*sqrt(q2)*radius)**2/q2
        endif
      else if(q2<1d-24)then
        spectrum(x+1,y+1,z+1)=pi/omega**2
      else
        spectrum(x+1,y+1,z+1)=4*pi*(1-exp(-q2/(4*omega**2)))/q2
      endif
    enddo;enddo;enddo
    plan=fftw_plan_dft_3d(ns(3),ns(2),ns(1),spectrum,spectrum,FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    if(.not.c_associated(plan))goto 900
    call fftw_execute_dft(plan,spectrum,spectrum)
    call fftw_destroy_plan(plan)
    op%kernel=real(spectrum,c_double)/real(product(ns),c_double)
    if(op%cosets>1)then
      ! Fold the full real-space kernel into the smaller Born-von Karman cell.
      ! Each k coset then uses only sources differing by multiples of reduction.
      allocate(folded_kernel(0:coarse_ns(1)-1,0:coarse_ns(2)-1,0:coarse_ns(3)-1))
      folded_kernel=0d0
      do az=0,reduction(3)-1;do ay=0,reduction(2)-1;do ax=0,reduction(1)-1
        folded_kernel=folded_kernel+op%kernel( &
          ax*coarse_ns(1):(ax+1)*coarse_ns(1)-1, &
          ay*coarse_ns(2):(ay+1)*coarse_ns(2)-1, &
          az*coarse_ns(3):(az+1)*coarse_ns(3)-1)
      enddo;enddo;enddo
      call move_alloc(folded_kernel,op%kernel)
    endif
    if(present(create_backend))then
      if(associated(create_backend))call create_backend(op%accelerator)
    endif
    if(allocated(op%accelerator))then
      ! Device plans are prepared once the distributed slot layout is known.
      op%auto_fft=.false.;op%contiguous_fft=.true.
      ierr=0;return
    endif
    allocate(op%work(op%block,ng,nk));dims=coarse(3:1:-1)
    if(op%auto_fft)then
      call choose_fft_layout(op,ierr)
      if(ierr/=0)goto 900
      ierr=1
    endif
    if(op%contiguous_fft)then
      ! Same allocation volume, interpreted as (nk,block,ng) during FFT work.
      op%forward=fftw_plan_many_dft(3,dims,op%block*ng,op%work,dims,1,nk, &
        op%work,dims,1,nk,FFTW_FORWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
      op%backward=fftw_plan_many_dft(3,dims,op%block*ng,op%work,dims,1,nk, &
        op%work,dims,1,nk,FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    else
    op%forward=fftw_plan_many_dft(3,dims,op%block*ng,op%work,dims,op%block*ng,1, &
      op%work,dims,op%block*ng,1,FFTW_FORWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    op%backward=fftw_plan_many_dft(3,dims,op%block*ng,op%work,dims,op%block*ng,1, &
      op%work,dims,op%block*ng,1,FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    endif
    if(.not.c_associated(op%forward).or..not.c_associated(op%backward))goto 900
    ierr=0;return
900 call exx_k_kernel_destroy(op)
  end subroutine

  ! Synthetic scratch only: no physical density or orbital is read here.
  ! Preserve the actual strided leading dimension while sampling FFT batches.
  subroutine choose_fft_layout(op,ierr)
    implicit none
    type(exx_k_kernel),target,intent(inout) :: op
    integer,intent(out) :: ierr
    complex(c_double_complex),allocatable :: tile(:,:,:)
    complex(c_double_complex),pointer :: kw(:,:,:)
    type(c_ptr) :: f(2),back(2)
    integer :: cols,b,nk,dims(3),flags,mode,pass,order,j,g0,k0,g,r,k,stat,threads
    real(c_double) :: stamp,elapsed(2,3),value
    ierr=1;f=c_null_ptr;back=c_null_ptr
    b=op%block;nk=op%nk;dims=op%mesh(3:1:-1)
    ! At least one g tile per OpenMP worker where the64MiB bound permits.
    ! One column is the minimum even if it exceeds the scratch budget.
    threads=1
    !$ threads=omp_get_max_threads()
    cols=max(1,min(op%ng,max(32,8*threads),int(4194304_int64/(int(b,int64)*nk))))
    allocate(tile(b,cols,nk),stat=stat)
    if(stat/=0)return
    kw(1:nk,1:b,1:op%ng)=>op%work
    flags=ior(FFTW_ESTIMATE,FFTW_UNALIGNED)
    f(1)=fftw_plan_many_dft(3,dims,b*cols,op%work,dims,b*op%ng,1, &
      op%work,dims,b*op%ng,1,FFTW_FORWARD,flags)
    back(1)=fftw_plan_many_dft(3,dims,b*cols,op%work,dims,b*op%ng,1, &
      op%work,dims,b*op%ng,1,FFTW_BACKWARD,flags)
    f(2)=fftw_plan_many_dft(3,dims,b*cols,op%work,dims,1,nk, &
      op%work,dims,1,nk,FFTW_FORWARD,flags)
    back(2)=fftw_plan_many_dft(3,dims,b*cols,op%work,dims,1,nk, &
      op%work,dims,1,nk,FFTW_BACKWARD,flags)
    do mode=1,2
      if(.not.c_associated(f(mode)).or..not.c_associated(back(mode)))goto 800
    enddo
    do pass=0,3
      do order=1,2
        mode=1+mod(order+pass,2)
        tile=(0.25d0,0.125d0)
        stamp=kernel_walltime()
        if(mode==1)then
          !$omp parallel do schedule(static)
          do k=1,nk
            op%work(:,1:cols,k)=tile(:,:,k)
          enddo
          !$omp end parallel do
        else
          !$omp parallel do private(k0,g,r,k) schedule(static)
          do g0=1,cols,8;do k0=1,nk,32
            do g=g0,min(cols,g0+7);do r=1,b;do k=k0,min(nk,k0+31)
              kw(k,r,g)=tile(r,g,k)
            enddo;enddo;enddo
          enddo;enddo
          !$omp end parallel do
        endif
        call fftw_execute_dft(f(mode),op%work,op%work)
        call fftw_execute_dft(back(mode),op%work,op%work)
        if(mode==1)then
          !$omp parallel do schedule(static)
          do k=1,nk
            tile(:,:,k)=op%work(:,1:cols,k)
          enddo
          !$omp end parallel do
        else
          !$omp parallel do private(k0,g,r,k) schedule(static)
          do g0=1,cols,8;do k0=1,nk,32
            do k=k0,min(nk,k0+31);do g=g0,min(cols,g0+7);do r=1,b
              tile(r,g,k)=kw(k,r,g)
            enddo;enddo;enddo
          enddo;enddo
          !$omp end parallel do
        endif
        value=kernel_walltime()-stamp
        if(pass>0)elapsed(mode,pass)=value
        if(maxval(abs(tile/real(nk,c_double)-(0.25d0,0.125d0)))>1d-10)goto 800
      enddo
    enddo
    do mode=1,2
      ! Median of three, excluding warm-up and planning.
      op%fft_trial_seconds(mode)=sum(elapsed(mode,:))-minval(elapsed(mode,:))-maxval(elapsed(mode,:))
    enddo
    ! Sub-resolution/tiny trials do not support a reliable speed decision.
    op%contiguous_fft=minval(op%fft_trial_seconds)>1d-5.and. &
      op%fft_trial_seconds(2)<0.95d0*op%fft_trial_seconds(1)
    ierr=0
800 do j=1,2
      if(c_associated(f(j)))call fftw_destroy_plan(f(j))
      if(c_associated(back(j)))call fftw_destroy_plan(back(j))
    enddo
  end subroutine

  subroutine exx_k_kernel_apply(op,source,target,action,rank,nproc,ierr)
    implicit none
    type(exx_k_kernel),target,intent(inout) :: op
    complex(c_double_complex),intent(in) :: source(:,:,:),target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(in) :: rank,nproc
    integer,intent(out) :: ierr
    complex(c_double_complex),allocatable :: s(:,:,:),t(:,:,:),tile(:,:)
    complex(c_double_complex),pointer :: kw(:,:,:)
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    integer :: lo,rows,ik,ki,b,j,r,offset(3),ng,nk,no,nt
    external :: zgemm
    ierr=1;action=zero
    if(.not.c_associated(op%forward))return
    if(op%phase_start/=1.or.size(op%phase,2)/=op%nk)return
    ng=op%ng;nk=op%nk;b=op%block
    if(nproc<1.or.rank<0.or.rank>=nproc)return
    if(size(source,1)/=ng.or.size(target,1)/=ng.or.size(source,3)/=nk.or.size(target,3)/=nk)return
    if(any(shape(action)/=shape(target)))return
    no=size(source,2);nt=size(target,2)
    if(no<1.or.nt<1)return
    if(.not.salmon_all_finite(real(source)).or..not.salmon_all_finite(aimag(source)))return
    if(.not.salmon_all_finite(real(target)).or..not.salmon_all_finite(aimag(target)))return
    allocate(s(ng,no,nk),t(ng,nt,nk))
    if(op%contiguous_fft)then
      allocate(tile(b,ng))
      kw(1:nk,1:b,1:ng)=>op%work
    endif
    do ik=1,nk
      do j=1,no;s(:,j,ik)=source(:,j,ik)*op%phase(:,ik);enddo
      do j=1,nt;t(:,j,ik)=target(:,j,ik)*op%phase(:,ik);enddo
    enddo
    do lo=1+rank*b,ng,nproc*b
      rows=min(b,ng-lo+1);op%work=zero
      do ki=1,nk
        ik=op%order(ki)
        if(op%contiguous_fft)then
          call zgemm('N','C',rows,ng,no,one,s(lo,1,ik),ng,s(1,1,ik),ng,zero,tile(1,1),b)
          kw(ki,1:rows,:)=tile(1:rows,:)
        else
          call zgemm('N','C',rows,ng,no,one,s(lo,1,ik),ng,s(1,1,ik),ng,zero,op%work(1,1,ki),b)
        endif
      enddo
      call execute_k_fft(op,op%forward)
      call multiply_kernel(op,lo,rows)
      call execute_k_fft(op,op%backward)
      do ki=1,nk
        ik=op%order(ki)
        if(op%contiguous_fft)then
          tile=kw(ki,:,:)
          call zgemm('N','N',rows,nt,ng,-one/real(op%fft_nk,c_double),tile(1,1),b, &
            t(1,1,ik),ng,zero,action(lo,1,ik),ng)
        else
        call zgemm('N','N',rows,nt,ng,-one/real(op%fft_nk,c_double),op%work(1,1,ki),b, &
          t(1,1,ik),ng,zero,action(lo,1,ik),ng)
        endif
        do j=1,nt
          action(lo:lo+rows-1,j,ik)=action(lo:lo+rows-1,j,ik)*conjg(op%phase(lo:lo+rows-1,ik))
        enddo
      enddo
    enddo
    ierr=0
  end subroutine

  ! K-distributed source/target/action; transpose density tiles, never orbitals.
  ! Caller supplies identical layout/kernel metadata and communicator size on all ranks.
  subroutine exx_k_kernel_apply_distributed(op,source,target,action,starts,counts,rank,comm,transpose_tiles,ierr,fill_density)
    implicit none
    type(exx_k_kernel),target,intent(inout) :: op
    complex(c_double_complex),intent(in) :: source(:,:,:),target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(in) :: starts(:),counts(:),rank,comm
    integer,intent(out) :: ierr
    interface
      subroutine transpose_tiles(send,recv,count,comm)
        import c_double_complex
        implicit none
        complex(c_double_complex),intent(in) :: send(:)
        complex(c_double_complex),intent(out) :: recv(:)
        integer,intent(in) :: count,comm
      end subroutine
      subroutine fill_density(j,lo,rows,density)
        import c_double_complex
        implicit none
        integer,intent(in) :: j,lo,rows
        complex(c_double_complex),intent(out) :: density(:,:)
      end subroutine
    end interface
    optional :: fill_density
    complex(c_double_complex),allocatable :: s(:,:,:),t(:,:,:),density_batch(:,:)
    complex(c_double_complex),allocatable,target :: send(:),recv(:)
    complex(c_double_complex),pointer :: sb(:,:,:,:),rb(:,:,:,:),flat_send(:,:,:),flat_recv(:,:,:)
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    integer :: ng,nk,np,nlocal,no,nt,b,km,nmsg,p,j,ki,ik,base,lo,rows,r,g,offset(3)
    integer,allocatable :: inverse(:),slot(:)
    complex(c_double_complex) :: valid_send(size(counts)),valid_recv(size(counts))
    logical :: valid
    integer :: batch_rows,active_rows
    real(c_double) :: stamp
    external :: zgemm
    ierr=1;action=zero
    op%seconds=0d0
    stamp=0d0
    if(op%profile)stamp=kernel_walltime()
    np=size(counts);ng=op%ng;nk=op%nk;b=op%block
    if(np<1.or.rank<0.or.rank>=np.or.size(starts)/=np)return
    ! All communicator members must participate even if one has invalid data.
    valid=c_associated(op%forward).or.allocated(op%accelerator)
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
    valid=valid.and.salmon_all_finite(real(source)).and.salmon_all_finite(aimag(source))
    valid=valid.and.salmon_all_finite(real(target)).and.salmon_all_finite(aimag(target))
    ! Fixed-size handshake must agree before variable-size tile collectives.
    valid_send=cmplx(b,merge(2,0,allocated(op%accelerator)),c_double)
    if(.not.valid)valid_send=cmplx(b,-1,c_double)
    call transpose_tiles(valid_send,valid_recv,1,comm)
    if(any(real(valid_recv)/=real(b,c_double)).or. &
       any(aimag(valid_recv)/=real(merge(2,0,allocated(op%accelerator)),c_double)))return
    km=maxval(counts);nmsg=b*ng*km
    allocate(send(nmsg*np),recv(nmsg*np),t(ng,nt,nlocal),inverse(nk))
    if(present(fill_density))then
      allocate(s(ng,0,nlocal))
    else
      allocate(s(ng,no,nlocal))
    endif
    sb(1:b,1:ng,1:km,1:np)=>send;rb(1:b,1:ng,1:km,1:np)=>recv
    do ki=1,nk;inverse(op%order(ki))=ki;enddo
    flat_send(1:b,1:ng,1:km*np)=>send
    flat_recv(1:b,1:ng,1:km*np)=>recv
    allocate(slot(nk))
    do p=1,np;do j=1,counts(p)
      slot(inverse(starts(p)+j-1))=j+(p-1)*km
    enddo;enddo
    if(allocated(op%accelerator))then
      call op%accelerator%prepare(op%n,op%mesh,b,op%kernel,op%point,op%shift,slot,km*np,ierr)
      valid_send=cmplx(merge(1,0,ierr/=0),0,c_double)
      call transpose_tiles(valid_send,valid_recv,1,comm)
      ierr=0
      if(any(real(valid_recv)/=0d0))ierr=1
      if(ierr/=0)return
    endif
    batch_rows=min(ng,np*b)
    if(present(fill_density))allocate(density_batch(batch_rows,ng))
    op%threads_used=1
    !$omp parallel private(j,g)
    !$omp single
!$  op%threads_used=omp_get_num_threads()
    !$omp end single
    !$omp do schedule(static)
    do j=1,nlocal
      if(.not.present(fill_density))then
        do g=1,no;s(:,g,j)=source(:,g,j)*op%phase(:,j);enddo
      endif
      do g=1,nt;t(:,g,j)=target(:,g,j)*op%phase(:,j);enddo
    enddo
    !$omp end do
    !$omp end parallel
    if(op%profile)call mark_stage(op,5,stamp)
    do base=1,ng,np*b
      send=zero
      ! Reconstruct each source only once for all destination rows this round.
      ! The explicit-source path retains its measured small-GEMM layout.
      if(present(fill_density))then
        active_rows=min(batch_rows,ng-base+1)
        do j=1,nlocal
          call fill_density(j,base,active_rows,density_batch)
          do p=0,np-1
            lo=p*b+1;rows=min(b,active_rows-lo+1)
            if(rows<=0)cycle
            sb(1:rows,:,j,p+1)=density_batch(lo:lo+rows-1,:)
          enddo
        enddo
      else
        ! BLAS owns threading; calls stay outside application OpenMP regions.
        do p=0,np-1
          lo=base+p*b;rows=min(b,ng-lo+1)
          if(rows<=0)cycle
          do j=1,nlocal
            call zgemm('N','C',rows,ng,no,one,s(lo,1,j),ng,s(1,1,j),ng,zero,sb(1,1,j,p+1),b)
          enddo
        enddo
      endif
      if(op%profile)call mark_stage(op,1,stamp)
      call transpose_tiles(send,recv,nmsg,comm)
      if(op%profile)call mark_stage(op,2,stamp)
      if(allocated(op%accelerator))then
        lo=base+rank*b;rows=max(0,min(b,ng-lo+1))
        call op%accelerator%apply(lo,rows,flat_recv,flat_send,ierr)
        ! A failed rank must not leave peers entering the next tile transpose.
        valid_send=cmplx(merge(1,0,ierr/=0),0,c_double)
        call transpose_tiles(valid_send,valid_recv,1,comm)
        ierr=0
        if(any(real(valid_recv)/=0d0))ierr=1
        if(ierr/=0)return
        ! Device packing/FFT/kernel/transfers are a combined measured stage.
        if(op%profile)call mark_stage(op,3,stamp)
      else
      ! Every k tile (including padded rows) is overwritten below.
      if(op%contiguous_fft)then
        call transpose_contiguous(op,flat_recv,slot,.false.)
      else
      !$omp parallel do schedule(static)
      do ki=1,nk
        op%work(:,:,ki)=flat_recv(:,:,slot(ki))
      enddo
      !$omp end parallel do
      endif
      if(op%profile)call mark_stage(op,5,stamp)
      call execute_k_fft(op,op%forward)
      if(op%profile)call mark_stage(op,3,stamp)
      lo=base+rank*b;rows=min(b,ng-lo+1)
      call multiply_kernel(op,lo,rows)
      if(op%profile)call mark_stage(op,4,stamp)
      call execute_k_fft(op,op%backward)
      if(op%profile)call mark_stage(op,3,stamp)
      ! Reuse send storage: useful entries are all overwritten; padded k
      ! entries remain defined from the density stage and are never consumed.
      if(op%contiguous_fft)then
        call transpose_contiguous(op,flat_send,slot,.true.)
      else
      !$omp parallel do schedule(static)
      do ki=1,nk
        flat_send(:,:,slot(ki))=op%work(:,:,ki)
      enddo
      !$omp end parallel do
      endif
      if(op%profile)call mark_stage(op,5,stamp)
      endif
      call transpose_tiles(send,recv,nmsg,comm)
      if(op%profile)call mark_stage(op,2,stamp)
      ! BLAS owns threading here; call only outside application OpenMP regions.
      do p=0,np-1
        lo=base+p*b;rows=min(b,ng-lo+1)
        if(rows<=0)cycle
        do j=1,nlocal
          call zgemm('N','N',rows,nt,ng,-one/real(op%fft_nk,c_double),rb(1,1,j,p+1),b, &
            t(1,1,j),ng,zero,action(lo,1,j),ng)
          do g=1,nt
            action(lo:lo+rows-1,g,j)=action(lo:lo+rows-1,g,j)*conjg(op%phase(lo:lo+rows-1,j))
          enddo
        enddo
      enddo
      if(op%profile)call mark_stage(op,6,stamp)
    enddo
    ierr=0
  end subroutine

  subroutine transpose_contiguous(op,buffer,slot,unpack)
    implicit none
    type(exx_k_kernel),target,intent(inout) :: op
    complex(c_double_complex),intent(inout) :: buffer(:,:,:)
    integer,intent(in) :: slot(:)
    logical,intent(in) :: unpack
    complex(c_double_complex),pointer :: kw(:,:,:)
    integer :: g0,k0,g,r,k
    kw(1:op%nk,1:op%block,1:op%ng)=>op%work
    ! Direction is fixed for the entire pack/unpack operation.
    if(unpack)then
      !$omp parallel do private(k0,g,r,k) schedule(static)
      do g0=1,op%ng,8
        do k0=1,op%nk,32
          do k=k0,min(op%nk,k0+31);do g=g0,min(op%ng,g0+7);do r=1,op%block
            buffer(r,g,slot(k))=kw(k,r,g)
          enddo;enddo;enddo
        enddo
      enddo
      !$omp end parallel do
    else
      !$omp parallel do private(k0,g,r,k) schedule(static)
      do g0=1,op%ng,8
        do k0=1,op%nk,32
          do g=g0,min(op%ng,g0+7);do r=1,op%block;do k=k0,min(op%nk,k0+31)
            kw(k,r,g)=buffer(r,g,slot(k))
          enddo;enddo;enddo
        enddo
      enddo
      !$omp end parallel do
    endif
  end subroutine

  subroutine execute_k_fft(op,plan)
    implicit none
    type(exx_k_kernel),intent(inout) :: op
    type(c_ptr),intent(in) :: plan
    integer :: first,last
    do first=1,op%nk,op%fft_nk
      last=first+op%fft_nk-1
      call fftw_execute_dft(plan,op%work(:,:,first:last),op%work(:,:,first:last))
    enddo
  end subroutine

  subroutine multiply_kernel(op,lo,rows)
    implicit none
    type(exx_k_kernel),target,intent(inout) :: op
    integer,intent(in) :: lo,rows
    complex(c_double_complex),pointer :: kw(:,:,:)
    integer :: ki,g,r,offset(3)
    if(op%contiguous_fft)then
      kw(1:op%nk,1:op%block,1:op%ng)=>op%work
      !$omp parallel do collapse(2) private(ki,offset) schedule(static)
      do g=1,op%ng;do r=1,max(0,rows)
        do ki=1,op%nk
          offset(1)=op%distance_index(op%point(1,lo+r-1),op%point(1,g),op%shift(1,ki)/op%n(1),1)
          offset(2)=op%distance_index(op%point(2,lo+r-1),op%point(2,g),op%shift(2,ki)/op%n(2),2)
          offset(3)=op%distance_index(op%point(3,lo+r-1),op%point(3,g),op%shift(3,ki)/op%n(3),3)
          kw(ki,r,g)=kw(ki,r,g)*op%kernel(offset(1),offset(2),offset(3))
        enddo
      enddo;enddo
      !$omp end parallel do
    else
      !$omp parallel do collapse(2) private(r,offset) schedule(static)
      do ki=1,op%nk;do g=1,op%ng;do r=1,max(0,rows)
        offset(1)=op%distance_index(op%point(1,lo+r-1),op%point(1,g),op%shift(1,ki)/op%n(1),1)
        offset(2)=op%distance_index(op%point(2,lo+r-1),op%point(2,g),op%shift(2,ki)/op%n(2),2)
        offset(3)=op%distance_index(op%point(3,lo+r-1),op%point(3,g),op%shift(3,ki)/op%n(3),3)
        op%work(r,g,ki)=op%work(r,g,ki)*op%kernel(offset(1),offset(2),offset(3))
      enddo;enddo;enddo
      !$omp end parallel do
    endif
  end subroutine

  real(c_double) function kernel_walltime()
    implicit none
    integer(int64) :: count,rate
    call system_clock(count,rate)
    kernel_walltime=real(count,c_double)/real(rate,c_double)
  end function

  ! Callers guard this helper so disabled profiling makes no timing calls.
  subroutine mark_stage(op,stage,stamp)
    implicit none
    type(exx_k_kernel),intent(inout) :: op
    integer,intent(in) :: stage
    real(c_double),intent(inout) :: stamp
    real(c_double) :: now
    now=kernel_walltime()
    op%seconds(stage)=op%seconds(stage)+now-stamp
    stamp=now
  end subroutine

  subroutine exx_k_kernel_destroy(op,status)
    implicit none
    type(exx_k_kernel),intent(inout) :: op
    integer,intent(out),optional :: status
    integer :: cleanup
    cleanup=0
    if(allocated(op%accelerator))then
      call op%accelerator%release(cleanup)
      deallocate(op%accelerator)
    endif
    if(present(status))status=cleanup
    if(c_associated(op%forward))call fftw_destroy_plan(op%forward)
    if(c_associated(op%backward))call fftw_destroy_plan(op%backward)
    op%forward=c_null_ptr;op%backward=c_null_ptr
    if(allocated(op%work))deallocate(op%work)
    if(allocated(op%distance_index))deallocate(op%distance_index)
    if(allocated(op%order))deallocate(op%order,op%point,op%shift,op%phase,op%kernel)
    op%n=0;op%mesh=0;op%ng=0;op%nk=0;op%block=0
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
end module
