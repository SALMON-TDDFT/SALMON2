! Full periodic sampled HSE kernel. Unit one-spin source occupations;
! neither the hybrid mixing fraction nor a second spin factor is included.
module hse_exchange
!$ use omp_lib, only: omp_get_num_threads,omp_get_max_threads
  use iso_fortran_env, only: int64
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  include 'fftw3.f03'
  public :: hse_kernel, hse_kernel_init, hse_kernel_apply, hse_kernel_destroy
  public :: hse_kernel_apply_distributed
  type hse_kernel
    integer :: n=0, mesh=0, ng=0, nk=0, block=0, phase_start=1, threads_used=1
    logical :: profile=.false.,contiguous_fft=.false.,auto_fft=.true.
    real(c_double) :: seconds(6)=0d0,fft_trial_seconds(2)=0d0
    integer, allocatable :: order(:), point(:,:),shift(:,:),distance_index(:,:,:)
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
    integer :: ns,nk,ng,i,j,x,y,z,ix,iy,iz,index,idx(3),dims(3),first,nphase,chosen,env_status,env_length
    character(64) :: setting
    real(c_double) :: pi,q2,q(3),scaled(3)
    complex(c_double_complex),allocatable :: spectrum(:,:,:)
    type(c_ptr) :: plan
    ierr=1
    call hse_kernel_destroy(op)
    if(n<1.or.mesh<1.or.block<1.or.h<=0.or.omega<=0) return
    if(.not.ieee_is_finite(h).or..not.ieee_is_finite(omega))return
    chosen=block
    call get_environment_variable('SALMON_HSE_BLOCK_ROWS',setting,length=env_length,status=env_status)
    if(env_status/=1)then
      if(env_status/=0.or.env_length<1.or.env_length>len(setting))return
      if(verify(trim(setting),'0123456789')/=0)return
      read(setting,*,iostat=env_status)chosen
      if(env_status/=0.or.chosen<1.or.chosen>n**3)return
    endif
    call get_environment_variable('SALMON_HSE_PROFILE',setting,status=env_status)
    op%profile=env_status==0.and.trim(setting)=='1'
    op%seconds=0d0
    op%contiguous_fft=.false.;op%auto_fft=.true.;op%fft_trial_seconds=0d0
    call get_environment_variable('SALMON_HSE_FFT_LAYOUT',setting,status=env_status)
    if(env_status/=1)then
      if(env_status/=0)return
      select case(trim(setting))
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
    nk=mesh**3;ng=n**3;ns=n*mesh;pi=acos(-1d0)
    if(size(k,1)/=3.or.size(k,2)/=nk.or..not.all(ieee_is_finite(k)))return
    first=1;nphase=nk
    if(present(phase_start))first=phase_start
    if(present(phase_count))nphase=phase_count
    if(first<1.or.nphase<0.or.first+nphase-1>nk)return
    op%phase_start=first
    op%n=n;op%mesh=mesh;op%ng=ng;op%nk=nk;op%block=min(chosen,ng)
    allocate(op%order(nk),op%point(3,ng),op%shift(3,nk),op%phase(ng,nphase),op%kernel(0:ns-1,0:ns-1,0:ns-1))
    ! Separable periodic coordinate differences: O(n*n*mesh) integers, not
    ! O(ng*ng*nk). Avoid integer modulo in every exchange-kernel multiply.
    allocate(op%distance_index(0:n-1,0:n-1,0:mesh-1))
    do z=0,mesh-1;do y=0,n-1;do x=0,n-1
      op%distance_index(x,y,z)=modulo(x-y-z*n,ns)
    enddo;enddo;enddo
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
900 call hse_kernel_destroy(op)
  end subroutine

  ! Synthetic scratch only: no physical density or orbital is read here.
  ! Preserve the actual strided leading dimension while sampling FFT batches.
  subroutine choose_fft_layout(op,ierr)
    type(hse_kernel),target,intent(inout) :: op
    integer,intent(out) :: ierr
    complex(c_double_complex),allocatable :: tile(:,:,:)
    complex(c_double_complex),pointer :: kw(:,:,:)
    type(c_ptr) :: f(2),back(2)
    integer :: cols,b,nk,dims(3),flags,mode,pass,order,j,g0,k0,g,r,k,stat,threads
    real(c_double) :: stamp,elapsed(2,3),value
    ierr=1;f=c_null_ptr;back=c_null_ptr
    b=op%block;nk=op%nk;dims=op%mesh
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

  subroutine hse_kernel_apply(op,source,target,action,rank,nproc,ierr)
    type(hse_kernel),target,intent(inout) :: op
    complex(c_double_complex),intent(in) :: source(:,:,:),target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(in) :: rank,nproc
    integer,intent(out) :: ierr
    complex(c_double_complex),allocatable :: s(:,:,:),t(:,:,:),tile(:,:)
    complex(c_double_complex),pointer :: kw(:,:,:)
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
      call fftw_execute_dft(op%forward,op%work,op%work)
      call multiply_kernel(op,lo,rows)
      call fftw_execute_dft(op%backward,op%work,op%work)
      do ki=1,nk
        ik=op%order(ki)
        if(op%contiguous_fft)then
          tile=kw(ki,:,:)
          call zgemm('N','N',rows,nt,ng,-one/real(nk,c_double),tile(1,1),b, &
            t(1,1,ik),ng,zero,action(lo,1,ik),ng)
        else
        call zgemm('N','N',rows,nt,ng,-one/real(nk,c_double),op%work(1,1,ki),b, &
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
  subroutine hse_kernel_apply_distributed(op,source,target,action,starts,counts,rank,transpose_tiles,ierr,fill_density)
    type(hse_kernel),target,intent(inout) :: op
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
      subroutine fill_density(j,lo,rows,density)
        import c_double_complex
        integer,intent(in) :: j,lo,rows
        complex(c_double_complex),intent(out) :: density(:,:)
      end subroutine
    end interface
    optional :: fill_density
    complex(c_double_complex),allocatable :: s(:,:,:),t(:,:,:),density_batch(:,:)
    complex(c_double_complex),allocatable,target :: send(:),recv(:)
    complex(c_double_complex),pointer :: sb(:,:,:,:),rb(:,:,:,:),flat_send(:,:,:),flat_recv(:,:,:)
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    integer :: ng,nk,np,nlocal,no,nt,b,km,nmsg,p,j,ki,ik,base,lo,rows,r,g,offset(3),ns
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
    ! Fixed-size handshake must agree before variable-size tile collectives.
    valid_send=cmplx(b,0,c_double)
    if(.not.valid)valid_send=cmplx(b,1,c_double)
    call transpose_tiles(valid_send,valid_recv,1)
    if(any(real(valid_recv)/=real(b,c_double)).or.any(aimag(valid_recv)/=0d0))return
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
      call transpose_tiles(send,recv,nmsg)
      if(op%profile)call mark_stage(op,2,stamp)
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
      call fftw_execute_dft(op%forward,op%work,op%work)
      if(op%profile)call mark_stage(op,3,stamp)
      lo=base+rank*b;rows=min(b,ng-lo+1)
      call multiply_kernel(op,lo,rows)
      if(op%profile)call mark_stage(op,4,stamp)
      call fftw_execute_dft(op%backward,op%work,op%work)
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
      call transpose_tiles(send,recv,nmsg)
      if(op%profile)call mark_stage(op,2,stamp)
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
      if(op%profile)call mark_stage(op,6,stamp)
    enddo
    ierr=0
  end subroutine

  subroutine transpose_contiguous(op,buffer,slot,unpack)
    type(hse_kernel),target,intent(inout) :: op
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

  subroutine multiply_kernel(op,lo,rows)
    type(hse_kernel),target,intent(inout) :: op
    integer,intent(in) :: lo,rows
    complex(c_double_complex),pointer :: kw(:,:,:)
    integer :: ki,g,r,offset(3)
    if(op%contiguous_fft)then
      kw(1:op%nk,1:op%block,1:op%ng)=>op%work
      !$omp parallel do collapse(2) private(ki,offset) schedule(static)
      do g=1,op%ng;do r=1,max(0,rows)
        do ki=1,op%nk
          offset(1)=op%distance_index(op%point(1,lo+r-1),op%point(1,g),op%shift(1,ki)/op%n)
          offset(2)=op%distance_index(op%point(2,lo+r-1),op%point(2,g),op%shift(2,ki)/op%n)
          offset(3)=op%distance_index(op%point(3,lo+r-1),op%point(3,g),op%shift(3,ki)/op%n)
          kw(ki,r,g)=kw(ki,r,g)*op%kernel(offset(1),offset(2),offset(3))
        enddo
      enddo;enddo
      !$omp end parallel do
    else
      !$omp parallel do collapse(2) private(r,offset) schedule(static)
      do ki=1,op%nk;do g=1,op%ng;do r=1,max(0,rows)
        offset(1)=op%distance_index(op%point(1,lo+r-1),op%point(1,g),op%shift(1,ki)/op%n)
        offset(2)=op%distance_index(op%point(2,lo+r-1),op%point(2,g),op%shift(2,ki)/op%n)
        offset(3)=op%distance_index(op%point(3,lo+r-1),op%point(3,g),op%shift(3,ki)/op%n)
        op%work(r,g,ki)=op%work(r,g,ki)*op%kernel(offset(1),offset(2),offset(3))
      enddo;enddo;enddo
      !$omp end parallel do
    endif
  end subroutine

  real(c_double) function kernel_walltime()
    integer(int64) :: count,rate
    call system_clock(count,rate)
    kernel_walltime=real(count,c_double)/real(rate,c_double)
  end function

  ! Callers guard this helper so disabled profiling makes no timing calls.
  subroutine mark_stage(op,stage,stamp)
    type(hse_kernel),intent(inout) :: op
    integer,intent(in) :: stage
    real(c_double),intent(inout) :: stamp
    real(c_double) :: now
    now=kernel_walltime()
    op%seconds(stage)=op%seconds(stage)+now-stamp
    stamp=now
  end subroutine

  subroutine hse_kernel_destroy(op)
    type(hse_kernel),intent(inout) :: op
    if(c_associated(op%forward))call fftw_destroy_plan(op%forward)
    if(c_associated(op%backward))call fftw_destroy_plan(op%backward)
    op%forward=c_null_ptr;op%backward=c_null_ptr
    if(allocated(op%work))deallocate(op%work)
    if(allocated(op%distance_index))deallocate(op%distance_index)
    if(allocated(op%order))deallocate(op%order,op%point,op%shift,op%phase,op%kernel)
    op%n=0;op%mesh=0;op%ng=0;op%nk=0;op%block=0
  end subroutine
end module
