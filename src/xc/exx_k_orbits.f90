#include "config.h"
! CPU space-group reuse of full-q exchange FFTs, in bounded first-coordinate tiles.
! This path uses representative target k points and the existing
! occupied-projector callback. MPI remains outside threaded local loops.
module exx_k_orbits
  use communication, only: comm_summation,comm_get_max
  use iso_c_binding
  use exx_k_exchange, only: exx_k_kernel,execute_kernel_fft
  use exx_collinear_fft, only: collinear_fft_plans,collinear_fft_create
  use exx_symmetry, only: s_exx_symmetry_map,symmetry_transform
  use exx_blas_threads, only: exx_blas_thread_control
!$  use omp_lib, only: omp_get_max_threads,omp_get_thread_num,omp_get_num_threads
  implicit none
  private
  include 'fftw3.f03'
  public :: apply_k_orbits,apply_spinor_k_orbits
contains
  subroutine execute_strided_lanes(op,fr,fs,slot,pair_delta,rows,ns,done)
    type(exx_k_kernel),target,intent(inout) :: op
    complex(c_double_complex),intent(in) :: fr(:,:,:)
    complex(c_double_complex),intent(inout) :: fs(:,:,:)
    integer,intent(in) :: slot(:),pair_delta(:,:,:),rows,ns(3)
    logical,intent(out) :: done
    complex(c_double_complex),pointer,contiguous :: work(:)
    done=.false.
    if(op%contiguous_fft)return
    if(.not.allocated(op%omp_fft%forward))return
    if(op%omp_fft%groups/=1.or.size(op%omp_fft%forward)<2)return
    work(1:size(op%work))=>op%work
    call execute_orbit_lanes(op%omp_fft,work,fr,fs,slot,pair_delta,rows,size(fr,2),ns, &
      op%shift,op%kernel,.false.,op%threads_used)
    done=.true.
  end subroutine

  subroutine execute_orbit_lanes(plans,work,fr,fs,slot,pair_delta,rows,ncols,ns,shift,kernel,contiguous,threads)
    type(collinear_fft_plans),intent(in) :: plans
    complex(c_double_complex),contiguous,intent(inout) :: work(:)
    complex(c_double_complex),intent(in) :: fr(:,:,:)
    complex(c_double_complex),intent(inout) :: fs(:,:,:)
    integer,intent(in) :: slot(:),pair_delta(:,:,:),rows,ncols,ns(3),shift(:,:)
    real(c_double),intent(in) :: kernel(0:,0:,0:)
    logical,intent(in) :: contiguous
    integer,intent(out) :: threads
    integer :: b,lanes,nk,workers,w,lo,hi,ki,lane,r,g,offset(3),at,lane_stride,k_stride
    b=size(fr,1);lanes=b*size(fr,2);nk=size(slot);workers=size(plans%forward)
    lane_stride=1;k_stride=lanes
    if(contiguous)then
      lane_stride=nk;k_stride=1
    endif
    threads=1
    !$omp parallel do schedule(static) num_threads(workers) private(lo,hi,ki,lane,r,g,offset,at)
    do w=1,workers
!$    if(w==1)threads=omp_get_num_threads()
      lo=(plans%offset(w)-1)/lane_stride+1
      hi=lanes
      if(w<workers)hi=(plans%offset(w+1)-1)/lane_stride
      do ki=1,nk
        do lane=lo,hi
          r=1+modulo(lane-1,b);g=1+(lane-1)/b
          work(1+(lane-1)*lane_stride+(ki-1)*k_stride)=fr(r,g,slot(ki))
        enddo
      enddo
      call fftw_execute_dft(plans%forward(w),work(plans%offset(w):),work(plans%offset(w):))
      do ki=1,nk
        do lane=lo,hi
          r=1+modulo(lane-1,b);g=1+(lane-1)/b
          if(r>rows.or.g>ncols)cycle
          offset=modulo(pair_delta(:,g,r)-shift(:,ki),ns)
          at=1+(lane-1)*lane_stride+(ki-1)*k_stride
          work(at)=work(at)*kernel(offset(1),offset(2),offset(3))
        enddo
      enddo
      call fftw_execute_dft(plans%backward(w),work(plans%offset(w):),work(plans%offset(w):))
      do ki=1,nk
        do lane=lo,hi
          r=1+modulo(lane-1,b);g=1+(lane-1)/b
          fs(r,g,slot(ki))=work(1+(lane-1)*lane_stride+(ki-1)*k_stride)
        enddo
      enddo
    enddo
    !$omp end parallel do
  end subroutine

  subroutine apply_k_orbits(op,map,source,target,action,starts,counts,rep_start,rank,comm,transpose_tiles,ierr)
    implicit none
    type(exx_k_kernel),target,intent(inout) :: op
    type(s_exx_symmetry_map),intent(in) :: map
    complex(c_double_complex),intent(in) :: source(:,:,:),target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(in) :: starts(:),counts(:),rep_start,rank,comm
    integer,intent(out) :: ierr
    interface
      subroutine transpose_tiles(send,recv,sendcounts,recvcounts,displs,comm)
        import c_double_complex
        implicit none
        complex(c_double_complex),intent(in) :: send(:)
        complex(c_double_complex),intent(out) :: recv(:)
        integer,intent(in) :: sendcounts(:),recvcounts(:),displs(:),comm
      end subroutine
      subroutine zgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        import c_double_complex
        implicit none
        character,intent(in) :: ta,tb
        integer,intent(in) :: m,n,k,lda,ldb,ldc
        complex(c_double_complex),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
        complex(c_double_complex),intent(inout) :: c(ldc,*)
      end subroutine
    end interface
    integer :: ng,nk,b,np,km,nmsg,nfirst,nsym,nlocal,nt,nrep,ns(3),i,j,g,s,p,r,base,nactive,lo,rows,ki
    integer :: local_error,total_error
    integer :: sendcounts(size(starts)),recvcounts(size(starts)),displs(size(starts))
    integer :: workers,worker,blas_state(3)
    integer :: density_workers,no,nops,symop,full_index,rep,stat,band,active_blocks
    complex(c_double_complex) :: weight,beta
    logical :: fused_fft
    integer :: a,dest,source_g,first_id,global_k,local_k,op_id,ip,offset(3),delta(3),c(3),idx(3)
    integer,allocatable :: forward(:,:),wrap(:,:,:),first(:),owner(:),operation(:),kperm(:,:),lookup(:)
    integer,allocatable :: inverse(:),slot(:),points(:)
    integer,allocatable :: pair_delta(:,:,:),target_k(:,:),reference_k(:)
    complex(c_double_complex),allocatable,target :: send(:),recv(:)
    complex(c_double_complex),pointer :: sb(:,:,:,:),rb(:,:,:,:),fs(:,:,:),fr(:,:,:),kw(:,:,:)
    complex(c_double_complex),allocatable :: rotated(:,:,:),selected(:,:,:),t(:,:,:),row(:,:),contracted(:,:),gauge(:,:),phase(:,:)
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    real(c_double) :: pi,scaled(3),v(3),length(3),angle,k0(3)
    real(c_double) :: stamp,stage_seconds(4),slowest(4)
    ierr=1;action=zero;stage_seconds=0d0
    ng=op%ng;nk=op%nk;b=op%block;np=size(starts);nsym=map%nsym
    nlocal=counts(rank+1);nt=size(target,2);nrep=size(target,3);ns=op%n*op%mesh
    local_error=0
    if(allocated(op%accelerator).or.op%nq/=nk)local_error=1
    if(size(target,1)/=ng.or.sum(counts)/=nk.or.any(shape(target)/=shape(action)))local_error=1
    if(any(op%n/=map%n).or.any(op%mesh/=op%mesh(1)))local_error=1
    call comm_summation(local_error,total_error,comm)
    if(total_error/=0)return
    allocate(forward(ng,nsym),wrap(3,ng,nsym),first(ng),owner(ng),operation(ng),gauge(ng,nsym))
    length=map%n*map%h;pi=acos(-1d0);k0=map%full_k(:,1)
    do s=1,nsym
      do g=1,ng
        v=matmul(map%rotation(:,:,s),real(op%point(:,g),c_double)/op%n)+map%translation(:,s)
        c=nint(v*op%n);wrap(:,g,s)=floor(real(c,c_double)/op%n);c=modulo(c,op%n)
        forward(g,s)=1+c(1)+op%n(1)*(c(2)+op%n(2)*c(3))
      enddo
      do g=1,ng
        source_g=map%grid(g,s)
        angle=sum(k0*map%h*real(op%point(:,source_g)-op%point(:,g),c_double))
        gauge(g,s)=exp(cmplx(0d0,angle,c_double))
      enddo
    enddo
    owner=0;operation=0;nfirst=0
    do g=1,ng
      if(owner(g)/=0)cycle
      nfirst=nfirst+1;first(nfirst)=g
      do s=1,nsym
        dest=forward(g,s)
        if(owner(dest)/=0)cycle
        owner(dest)=nfirst;operation(dest)=s
      enddo
    enddo
    deallocate(forward)
    allocate(kperm(nk,nsym),lookup(nk))
    lookup=0
    do i=1,nk
      scaled=(map%full_k(:,i)-k0)*length*op%mesh/(2*pi)
      idx=modulo(nint(scaled),op%mesh)
      j=1+idx(1)+op%mesh(1)*(idx(2)+op%mesh(2)*idx(3));lookup(j)=i
    enddo
    do s=1,nsym;do i=1,nk
      v=matmul(transpose(map%rotation(:,:,s)),map%full_k(:,i))
      scaled=(v-k0)*length*op%mesh/(2*pi)
      if(maxval(abs(scaled-anint(scaled)))>1d-8)return
      idx=modulo(nint(scaled),op%mesh)
      j=1+idx(1)+op%mesh(1)*(idx(2)+op%mesh(2)*idx(3));kperm(i,s)=lookup(j)
      if(kperm(i,s)<1)return
    enddo;enddo
    deallocate(lookup)
    ! A star remains on one k owner, including all inverse-rotated targets.
    allocate(target_k(nsym,nrep),reference_k(nrep))
    local_error=0
    do j=1,nrep
      global_k=map%first(rep_start+j-1)
      reference_k(j)=global_k-starts(rank+1)+1
      do s=1,nsym
        local_k=kperm(global_k,s)-starts(rank+1)+1
        target_k(s,j)=local_k
        if(local_k<1.or.local_k>nlocal)local_error=1
      enddo
    enddo
    call comm_summation(local_error,total_error,comm)
    if(total_error/=0)return
    km=maxval(counts);nmsg=b*ng*km
    displs=[((p-1)*nmsg,p=1,np)]
    workers=1
#if defined(HAVE_EXX_FUJITSU_THREADS) || defined(HAVE_EXX_NVPL_THREADS) || defined(HAVE_EXX_OPENBLAS_THREADS)
!$  if(nrep>0)workers=min(ng,omp_get_max_threads())
#endif
    no=size(source,2);nops=0
    do j=1,nlocal
      nops=max(nops,map%multiplicity(starts(rank+1)+j-1))
    enddo
    density_workers=1
#if defined(HAVE_EXX_FUJITSU_THREADS) || defined(HAVE_EXX_NVPL_THREADS) || defined(HAVE_EXX_OPENBLAS_THREADS)
!$  density_workers=min(max(1,nlocal*min(np,(nfirst+b-1)/b)),omp_get_max_threads())
#endif
    allocate(send(nmsg*np),recv(nmsg*np),points(np*b),t(ng,nt,nrep),phase(27,nrep))
    allocate(rotated(ng,no,nlocal),selected(b,no,density_workers))
    allocate(row(ng,workers),contracted(nt,workers))
    allocate(pair_delta(3,ng,b))
    sb(1:b,1:ng,1:km,1:np)=>send;rb(1:b,1:ng,1:km,1:np)=>recv
    fs(1:b,1:ng,1:km*np)=>send;fr(1:b,1:ng,1:km*np)=>recv
    allocate(inverse(nk),slot(nk))
    do i=1,nk;inverse(op%order(i))=i;enddo
    do p=1,np;do j=1,counts(p)
      slot(inverse(starts(p)+j-1))=j+(p-1)*km
    enddo;enddo
    deallocate(inverse)
    do j=1,nrep
      global_k=map%first(rep_start+j-1);local_k=reference_k(j)
      do i=1,nt;t(:,i,j)=target(:,i,j)*op%phase(:,local_k);enddo
      do i=1,27
        delta=[mod(i-1,3)-1,mod((i-1)/3,3)-1,(i-1)/9-1]
        angle=sum(map%full_k(:,global_k)*length*delta)
        phase(i,j)=exp(cmplx(0d0,angle,c_double))
      enddo
    enddo
    if(rank==0)write(*,'(a,2i12)')'EXX_FFT_ORBITS first points/full: ',nfirst,ng
    do base=1,nfirst,np*b
      if(op%profile)stamp=clock_now()
      nactive=min(np*b,nfirst-base+1);points(:nactive)=first(base:base+nactive-1)
      active_blocks=(nactive+b-1)/b
      sendcounts=0;recvcounts=0
      sendcounts(:active_blocks)=b*ng*nlocal
      if(rank<active_blocks)recvcounts=b*ng*counts
#ifdef HAVE_EXX_OPENBLAS_THREADS
      blas_state=0
      if(density_workers>1)then
        call exx_blas_thread_control(1,blas_state)
        if(blas_state(1)==0)density_workers=1
      endif
#endif
      !$omp parallel if(density_workers>1) num_threads(density_workers) &
      !$omp private(worker,symop,j,p,full_index,rep,stat,band,lo,rows,r,weight,beta) &
#if defined(HAVE_EXX_FUJITSU_THREADS) || defined(HAVE_EXX_NVPL_THREADS)
      !$omp private(blas_state)
      blas_state=0
      if(density_workers>1)call exx_blas_thread_control(1,blas_state)
#else
      !$omp shared(blas_state)
#endif
      worker=1
!$  worker=omp_get_thread_num()+1
      !$omp do collapse(2) schedule(static)
      do p=1,active_blocks
        do j=1,nlocal
          if(no==0)then
            sb(:,:,j,p)=zero
          else if(p==active_blocks)then
            rows=nactive-(p-1)*b
            sb(rows+1:,:,j,p)=zero
          endif
        enddo
      enddo
      !$omp end do
      ! Preserve the little-group sum order while sharing one transformed source across its row blocks.
      do symop=1,nops
        !$omp do schedule(static)
        do j=1,nlocal
          full_index=starts(rank+1)+j-1
          if(symop>map%multiplicity(full_index))cycle
          rep=map%owner(full_index)-rep_start+1
          call symmetry_transform(map,full_index,symop,source(:,:,rep),rotated(:,:,j),stat)
          if(stat/=0)error stop 'HSE symmetry source transformation failed'
          do band=1,no
            rotated(:,band,j)=rotated(:,band,j)*op%phase(:,j)
          enddo
        enddo
        !$omp end do
        beta=one
        if(symop==1)beta=zero
        !$omp do collapse(2) schedule(static)
        do p=1,active_blocks
          do j=1,nlocal
            full_index=starts(rank+1)+j-1
            if(symop>map%multiplicity(full_index).or.no==0)cycle
            lo=(p-1)*b+1;rows=min(b,nactive-lo+1)
            weight=cmplx(1d0/map%multiplicity(full_index),0d0,c_double)
            do band=1,no
              do r=1,rows
                selected(r,band,worker)=rotated(points(lo+r-1),band,j)
              enddo
            enddo
            call zgemm('N','C',rows,ng,no,weight,selected(1,1,worker),b,rotated(1,1,j),ng,beta,sb(:,:,j,p),b)
          enddo
        enddo
        !$omp end do
      enddo
#if defined(HAVE_EXX_FUJITSU_THREADS) || defined(HAVE_EXX_NVPL_THREADS)
      call exx_blas_thread_control(0,blas_state)
#endif
      !$omp end parallel
#ifdef HAVE_EXX_OPENBLAS_THREADS
      call exx_blas_thread_control(0,blas_state)
#endif
      if(op%profile)call mark(1)
      call transpose_tiles(send,recv,sendcounts,recvcounts,displs,comm)
      if(op%profile)call mark(2)
      rows=max(0,min(b,nfirst-base-rank*b+1))
      if(rows>0)then
      ! The point-pair displacement is independent of the FFT k index.
      do r=1,rows
        a=first(base+rank*b+r-1)
        do g=1,ng
          pair_delta(:,g,r)=op%point(:,a)-op%point(:,g)
        enddo
      enddo
      call execute_strided_lanes(op,fr,fs,slot,pair_delta,rows,ns,fused_fft)
      if(.not.fused_fft)then
      if(op%contiguous_fft)then
        kw(1:nk,1:b,1:ng)=>op%work
        !$omp parallel do private(ki,r) schedule(static)
        do g=1,ng
          do r=1,b;do ki=1,nk
            kw(ki,r,g)=fr(r,g,slot(ki))
          enddo;enddo
        enddo
        !$omp end parallel do
      else
        !$omp parallel do schedule(static)
        do ki=1,nk;op%work(:,:,ki)=fr(:,:,slot(ki));enddo
        !$omp end parallel do
      endif
      call execute_kernel_fft(op,op%forward)
      !$omp parallel do collapse(2) private(r,offset) schedule(static)
      do ki=1,nk;do g=1,ng;do r=1,rows
        offset=modulo(pair_delta(:,g,r)-op%shift(:,ki),ns)
        if(op%contiguous_fft)then
          kw(ki,r,g)=kw(ki,r,g)*op%kernel(offset(1),offset(2),offset(3))
        else
          op%work(r,g,ki)=op%work(r,g,ki)*op%kernel(offset(1),offset(2),offset(3))
        endif
      enddo;enddo;enddo
      !$omp end parallel do
      call execute_kernel_fft(op,op%backward)
      if(op%contiguous_fft)then
        !$omp parallel do private(g,r) schedule(static)
        do ki=1,nk
          do g=1,ng;do r=1,b
            fs(r,g,slot(ki))=kw(ki,r,g)
          enddo;enddo
        enddo
        !$omp end parallel do
      else
        !$omp parallel do schedule(static)
        do ki=1,nk;fs(:,:,slot(ki))=op%work(:,:,ki);enddo
        !$omp end parallel do
      endif
      endif
      endif
      if(op%profile)call mark(3)
      call transpose_tiles(send,recv,recvcounts,sendcounts,displs,comm)
      if(op%profile)call mark(2)
      ! Reconstruct only representative-target rows; no expanded target/action.
      if(nrep==0)cycle
#ifdef HAVE_EXX_OPENBLAS_THREADS
      blas_state=0
      if(workers>1)then
        call exx_blas_thread_control(1,blas_state)
        if(blas_state(1)==0)workers=1
      endif
#endif
      !$omp parallel if(workers>1) num_threads(workers) &
      !$omp private(worker,first_id,p,r,a,op_id,j,local_k,g,source_g,delta,ip) &
#if defined(HAVE_EXX_FUJITSU_THREADS) || defined(HAVE_EXX_NVPL_THREADS)
      !$omp private(blas_state)
      blas_state=0
      if(workers>1)call exx_blas_thread_control(1,blas_state)
#else
      !$omp shared(blas_state)
#endif
      worker=1
!$    worker=omp_get_thread_num()+1
      !$omp do schedule(static)
      do dest=1,ng
        first_id=owner(dest)
        if(first_id<base.or.first_id>=base+nactive)cycle
        p=(first_id-base)/b+1;r=mod(first_id-base,b)+1
        a=first(first_id);op_id=operation(dest)
        do j=1,nrep
          local_k=target_k(op_id,j)
          do g=1,ng
            source_g=map%grid(g,op_id)
            delta=wrap(:,source_g,op_id)-wrap(:,a,op_id)
            if(any(abs(delta)>1))error stop 'EXX orbit wrap exceeds one cell'
            ip=1+delta(1)+1+3*(delta(2)+1+3*(delta(3)+1))
            row(g,worker)=rb(r,source_g,local_k,p)*phase(ip,j)*gauge(dest,op_id)*conjg(gauge(g,op_id))
          enddo
          call zgemm('N','N',1,nt,ng,-one/real(nk,c_double),row(1,worker),1,t(1,1,j),ng,zero,contracted(1,worker),1)
          action(dest,:,j)=contracted(:,worker)
          local_k=reference_k(j)
          action(dest,:,j)=action(dest,:,j)*conjg(op%phase(dest,local_k))
        enddo
      enddo
      !$omp end do nowait
#if defined(HAVE_EXX_FUJITSU_THREADS) || defined(HAVE_EXX_NVPL_THREADS)
      call exx_blas_thread_control(0,blas_state)
#endif
      !$omp end parallel
#ifdef HAVE_EXX_OPENBLAS_THREADS
      call exx_blas_thread_control(0,blas_state)
#endif
      if(op%profile)call mark(4)
    enddo
    if(op%profile)then
      call comm_get_max(stage_seconds,slowest,4,comm)
      if(rank==0)write(*,"(a,4es14.5)") "HSE_ORBIT_PROFILE max-rank seconds density transpose fft_pack_kernel action=",slowest
    endif
    ierr=0
  contains
    subroutine mark(stage)
      integer,intent(in) :: stage
      real(c_double) :: now
      now=clock_now()
      stage_seconds(stage)=stage_seconds(stage)+now-stamp
      stamp=now
    end subroutine
  end subroutine

  subroutine apply_spinor_k_orbits(op,map,source,target,action,starts,counts, &
                                  rep_start,rank,comm,transpose_tiles,ierr)
    use communication, only: comm_summation,comm_bcast
    use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
    use iso_fortran_env,only:int64
    use exx_k_exchange, only: exx_k_kernel
    use exx_symmetry, only: s_exx_symmetry_map,symmetry_spin_rotation
    implicit none
    type(exx_k_kernel),target,intent(inout) :: op
    type(s_exx_symmetry_map),intent(in) :: map
    ! The expanded source includes the k-dependent phase for both spin components.
    complex(c_double_complex),contiguous,target,intent(in) :: source(:,:,:)
    complex(c_double_complex),intent(in) :: target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(in) :: starts(:),counts(:),rep_start,rank,comm
    integer,intent(out) :: ierr
    interface
      subroutine transpose_tiles(send,recv,sendcounts,recvcounts,displs,comm)
        import c_double_complex
        implicit none
        complex(c_double_complex),intent(in) :: send(:)
        complex(c_double_complex),intent(out) :: recv(:)
        integer,intent(in) :: sendcounts(:),recvcounts(:),displs(:),comm
      end subroutine
      subroutine zgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        import c_double_complex
        implicit none
        character,intent(in) :: ta,tb
        integer,intent(in) :: m,n,k,lda,ldb,ldc
        complex(c_double_complex),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
        complex(c_double_complex),intent(inout) :: c(ldc,*)
      end subroutine
    end interface
    integer :: ng,nk,b,np,km,nmsg,nfirst,nsym,nlocal,nt,nrep,ns(3),i,j,g,s,p,r,base,nactive,lo,rows,dest_rows
    integer :: local_error,total_error,component,spin_left,spin_right,spin_out,spin_target
    integer :: settings(13),expected(13),divisor,chunk,colbase,ncols,stat,nchars
    integer :: nsource,density_workers,right_start,band,active_blocks
    integer :: sendcounts(size(starts)),recvcounts(size(starts)),displs(size(starts))
    character(32) :: setting
    integer :: workers,worker,fft_workers,blas_state(3)
    integer,allocatable :: pair_delta(:,:,:)
    integer :: a,dest,source_g,first_id,global_k,source_k,local_k,op_id,ip,delta(3),c(3),idx(3)
    integer,allocatable :: forward(:,:),wrap(:,:,:),first(:),owner(:),operation(:),kperm(:,:),lookup(:)
    integer,allocatable :: inverse(:),slot(:),points(:)
    complex(c_double_complex),allocatable,target :: send(:),recv(:),flat(:)
    complex(c_double_complex),allocatable :: spin_rotation(:,:,:)
    complex(c_double_complex) :: coefficient
    complex(c_double_complex),pointer,contiguous :: source_flat(:)
    complex(c_double_complex),pointer :: sb(:,:,:,:),rb(:,:,:,:),fs(:,:,:),fr(:,:,:)
    complex(c_double_complex),allocatable :: left(:,:,:),t(:,:,:),row(:,:),contracted(:,:),gauge(:,:),phase(:,:)
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    real(c_double) :: pi,scaled(3),v(3),angle,k0(3),stamp
    ierr=1;action=zero;op%seconds=0d0;op%threads_used=1;op%pair_fft_pairs=0_int64
    stamp=clock_now()
    ng=op%ng;nk=op%nk;b=op%block;np=size(starts);nsym=map%nsym
    nt=size(target,2);nrep=size(target,3);ns=op%n*op%mesh
    local_error=0
    if(size(counts)/=np.or.rank<0.or.rank>=np.or.np<1)local_error=1
    if(min(ng,nk,b,nsym)<1.or.any(counts<0).or.any(starts<1))local_error=1
    if(allocated(op%accelerator).or.op%nq/=nk)local_error=1
    if(.not.allocated(op%phase))local_error=1
    if(size(target,1)/=2*ng.or.sum(counts)/=nk.or.any(shape(target)/=shape(action)))local_error=1
    if(any(op%n/=map%n))local_error=1
    call comm_summation(local_error,total_error,comm)
    if(total_error/=0)return
    nlocal=counts(rank+1);nsource=size(source,2)
    if(size(source,1)/=2*ng.or.size(source,3)/=nlocal)local_error=1
    if(starts(1)/=1.or.op%phase_start/=starts(rank+1))local_error=1
    do p=2,np
      if(starts(p)/=starts(p-1)+counts(p-1))local_error=1
    enddo
    if(size(op%phase,2)/=nlocal.or.rep_start<1.or.rep_start+nrep>map%nrep+1)local_error=1
    if(maxval(counts)>0)then
      if(int(b,int64)*ng>int(huge(nmsg),int64)/np/maxval(counts))local_error=1
    endif
    divisor=4;setting=' '
    call get_environment_variable('SALMON_EXX_ORBIT_TILE_DIVISOR',setting,length=nchars,status=stat)
    if(stat/=1)then
      if(stat/=0.or.nchars<1.or.nchars>len(setting))then
        local_error=1
      else if(verify(setting(:nchars),'0123456789')/=0)then
        local_error=1
      else
        read(setting(:nchars),*,iostat=stat)divisor
        if(stat/=0.or.divisor<1)local_error=1
      endif
    endif
    settings=[op%n,op%mesh,ng,nk,b,nsym,nt,map%nrep,divisor];expected=settings
    call comm_bcast(expected,comm,0)
    if(any(settings/=expected))local_error=1
    do j=1,nrep;do i=1,nt;do g=1,2*ng
      if(.not.ieee_is_finite(real(target(g,i,j))).or..not.ieee_is_finite(aimag(target(g,i,j))))local_error=1
    enddo;enddo;enddo
    call comm_summation(local_error,total_error,comm)
    if(total_error/=0)return
    allocate(forward(ng,nsym),wrap(3,ng,nsym),first(ng),owner(ng),operation(ng),gauge(ng,nsym))
    allocate(spin_rotation(2,2,nsym))
    do s=1,nsym
      call symmetry_spin_rotation(map%cart_rotation(:,:,s),spin_rotation(:,:,s),local_error)
      if(local_error/=0)return
    enddo
    pi=acos(-1d0);k0=map%full_k(:,1)
    do s=1,nsym
      do g=1,ng
        v=matmul(map%rotation(:,:,s),real(op%point(:,g),c_double)/op%n)+map%translation(:,s)
        c=nint(v*op%n);wrap(:,g,s)=floor(real(c,c_double)/op%n);c=modulo(c,op%n)
        forward(g,s)=1+c(1)+op%n(1)*(c(2)+op%n(2)*c(3))
      enddo
      do g=1,ng
        source_g=map%grid(g,s)
        angle=dot_product(k0,matmul(map%cell,real(op%point(:,source_g)-op%point(:,g),c_double)/op%n))
        gauge(g,s)=exp(cmplx(0d0,angle,c_double))
      enddo
    enddo
    owner=0;operation=0;nfirst=0
    do g=1,ng
      if(owner(g)/=0)cycle
      nfirst=nfirst+1;first(nfirst)=g
      do s=1,nsym
        dest=forward(g,s)
        if(owner(dest)/=0)cycle
        owner(dest)=nfirst;operation(dest)=s
      enddo
    enddo
    ! Ownership is complete; do not overlap this map with exchange scratch.
    ! Retain the inverse grid map to stream second-coordinate tiles.
    allocate(kperm(nk,nsym),lookup(nk))
    lookup=0
    do i=1,nk
      scaled=matmul(transpose(map%cell),map%full_k(:,i)-k0)*op%mesh/(2*pi)
      idx=modulo(nint(scaled),op%mesh)
      j=1+idx(1)+op%mesh(1)*(idx(2)+op%mesh(2)*idx(3));lookup(j)=i
    enddo
    do s=1,nsym;do i=1,nk
      v=matmul(transpose(map%cart_rotation(:,:,s)),map%full_k(:,i))
      scaled=matmul(transpose(map%cell),v-k0)*op%mesh/(2*pi)
      if(maxval(abs(scaled-anint(scaled)))>1d-8)return
      idx=modulo(nint(scaled),op%mesh)
      j=1+idx(1)+op%mesh(1)*(idx(2)+op%mesh(2)*idx(3));kperm(i,s)=lookup(j)
      if(kperm(i,s)<1)return
    enddo;enddo
    deallocate(lookup)
    ! A star remains on one k owner, including all inverse-rotated targets.
    local_error=0
    do j=1,nrep
      global_k=map%first(rep_start+j-1)
      do s=1,nsym
        local_k=kperm(global_k,s)-starts(rank+1)+1
        if(local_k<1.or.local_k>nlocal)local_error=1
      enddo
    enddo
    call comm_summation(local_error,total_error,comm)
    if(total_error/=0)return
    km=maxval(counts)
    allocate(inverse(nk),slot(nk))
    do i=1,nk;inverse(op%order(i))=i;enddo
    do p=1,np;do j=1,counts(p)
      slot(inverse(starts(p)+j-1))=j+(p-1)*km
    enddo;enddo
    deallocate(inverse)
    chunk=1+(ng-1)/divisor;nmsg=b*chunk*km
    displs=[((p-1)*nmsg,p=1,np)]
    workers=1
#if defined(HAVE_EXX_FUJITSU_THREADS) || defined(HAVE_EXX_NVPL_THREADS) || defined(HAVE_EXX_OPENBLAS_THREADS)
!$  if(nrep>0)workers=min(ng,omp_get_max_threads())
#endif
    density_workers=1
#if defined(HAVE_EXX_FUJITSU_THREADS) || defined(HAVE_EXX_NVPL_THREADS) || defined(HAVE_EXX_OPENBLAS_THREADS)
!$  density_workers=min(max(1,nlocal*min(np,(nfirst+b-1)/b)),omp_get_max_threads())
#endif
    source_flat(1:size(source))=>source
    allocate(send(nmsg*np),recv(nmsg*np),flat(b*chunk*nk),points(np*b),phase(27,nrep))
    allocate(left(b,nsource,density_workers))
    allocate(t(chunk,nt,workers),row(chunk,workers),contracted(nt,workers),pair_delta(3,chunk,b))
    sb(1:b,1:chunk,1:km,1:np)=>send;rb(1:b,1:chunk,1:km,1:np)=>recv
    fs(1:b,1:chunk,1:km*np)=>send;fr(1:b,1:chunk,1:km*np)=>recv
    do j=1,nrep
      global_k=map%first(rep_start+j-1);local_k=global_k-starts(rank+1)+1
      do i=1,27
        delta=[mod(i-1,3)-1,mod((i-1)/3,3)-1,(i-1)/9-1]
        angle=dot_product(map%full_k(:,global_k),matmul(map%cell,real(delta,c_double)))
        phase(i,j)=exp(cmplx(0d0,angle,c_double))
      enddo
    enddo
    fft_workers=1
!$  fft_workers=min(b*chunk,omp_get_max_threads())
    local_error=0
    if(rank*b<nfirst)then
      if(op%orbit_tile_lanes/=b*chunk.or.op%orbit_tile_workers/=fft_workers.or. &
         (op%orbit_tile_contiguous.neqv.op%contiguous_fft))then
        call collinear_fft_create(op%omp_orbit_tile,op%mesh,nk,b*chunk,b*chunk,op%contiguous_fft,flat,local_error)
        op%orbit_tile_lanes=0
        if(local_error==0)op%orbit_tile_lanes=b*chunk
        op%orbit_tile_workers=fft_workers
        op%orbit_tile_contiguous=op%contiguous_fft
      endif
    endif
    call comm_summation(local_error,total_error,comm)
    if(total_error/=0)return
    if(rank==0)write(*,'(a,2i12)')'EXX_SPINOR_ORBIT_TILE divisor/columns: ',divisor,chunk
    if(rank==0)write(*,'(a,2i12)')'EXX_SPINOR_FFT_ORBITS first points/full: ',nfirst,ng
    if(op%profile)call mark(5)
    do base=1,nfirst,np*b
      nactive=min(np*b,nfirst-base+1);points(:nactive)=first(base:base+nactive-1)
      active_blocks=(nactive+b-1)/b
      sendcounts=0;recvcounts=0
      sendcounts(:active_blocks)=b*chunk*nlocal
      if(rank<active_blocks)recvcounts=b*chunk*counts
      do colbase=1,ng,chunk
      ncols=min(chunk,ng-colbase+1)
      rows=max(0,min(b,nfirst-base-rank*b+1))
      do r=1,rows
        a=first(base+rank*b+r-1)
        do g=1,ncols
          pair_delta(:,g,r)=op%point(:,a)-op%point(:,colbase+g-1)
        enddo
      enddo
      do component=1,4
      spin_left=1+modulo(component-1,2);spin_right=1+(component-1)/2
#ifdef HAVE_EXX_OPENBLAS_THREADS
      blas_state=0
      if(density_workers>1)then
        call exx_blas_thread_control(1,blas_state)
        if(blas_state(1)==0)density_workers=1
      endif
#endif
      !$omp parallel if(density_workers>1) num_threads(density_workers) &
      !$omp private(worker,p,j,lo,dest_rows,band,r,right_start) &
#if defined(HAVE_EXX_FUJITSU_THREADS) || defined(HAVE_EXX_NVPL_THREADS)
      !$omp private(blas_state)
      blas_state=0
      if(density_workers>1)call exx_blas_thread_control(1,blas_state)
#else
      !$omp shared(blas_state)
#endif
      worker=1
!$  worker=omp_get_thread_num()+1
      !$omp do collapse(2) schedule(static)
      do p=1,active_blocks
        do j=1,nlocal
          lo=(p-1)*b+1;dest_rows=min(b,nactive-lo+1)
          if(nsource==0)then
            sb(:,:,j,p)=zero
            cycle
          endif
          ! The GEMM writes the active rectangle directly into its MPI block.
          sb(dest_rows+1:,:,j,p)=zero
          sb(:dest_rows,ncols+1:,j,p)=zero
          do band=1,nsource
            do r=1,dest_rows
              left(r,band,worker)=source(points(lo+r-1)+(spin_left-1)*ng,band,j)
            enddo
          enddo
          right_start=colbase+(spin_right-1)*ng+(j-1)*2*ng*nsource
          call zgemm('N','C',dest_rows,ncols,nsource,one,left(1,1,worker),b, &
            source_flat(right_start:),2*ng,zero,sb(:,:,j,p),b)
        enddo
      enddo
      !$omp end do nowait
#if defined(HAVE_EXX_FUJITSU_THREADS) || defined(HAVE_EXX_NVPL_THREADS)
      call exx_blas_thread_control(0,blas_state)
#endif
      !$omp end parallel
#ifdef HAVE_EXX_OPENBLAS_THREADS
      call exx_blas_thread_control(0,blas_state)
#endif
      if(op%profile)call mark(1)
      call transpose_tiles(send,recv,sendcounts,recvcounts,displs,comm)
      if(op%profile)call mark(2)
      if(rows>0)then
        call execute_orbit_lanes(op%omp_orbit_tile,flat,fr,fs,slot,pair_delta,rows,ncols,ns, &
          op%shift,op%kernel,op%contiguous_fft,op%threads_used)
        op%pair_fft_pairs=op%pair_fft_pairs+int(rows,int64)*ncols
      endif
      ! Packing, FFTs and kernel multiplication share one OpenMP region.
      if(op%profile)call mark(3)
      call transpose_tiles(send,recv,recvcounts,sendcounts,displs,comm)
      if(op%profile)call mark(2)
      ! Linearity permits one spin block at a time: U rho U-dagger is
      ! contracted componentwise, avoiding four retained communication buffers.
#ifdef HAVE_EXX_OPENBLAS_THREADS
      blas_state=0
      if(workers>1)then
        call exx_blas_thread_control(1,blas_state)
        if(blas_state(1)==0)workers=1
      endif
#endif
      !$omp parallel if(workers>1) num_threads(workers) &
      !$omp private(worker,first_id,p,r,a,op_id,j,global_k,source_k,local_k,g,source_g,i,delta,ip, &
      !$omp         spin_out,spin_target,coefficient) &
#if defined(HAVE_EXX_FUJITSU_THREADS) || defined(HAVE_EXX_NVPL_THREADS)
      !$omp private(blas_state)
      blas_state=0
      if(workers>1)call exx_blas_thread_control(1,blas_state)
#else
      !$omp shared(blas_state)
#endif
      worker=1
!$  worker=omp_get_thread_num()+1
      !$omp do schedule(static)
      do dest=1,ng
        first_id=owner(dest)
        if(first_id<base.or.first_id>=base+nactive)cycle
        p=(first_id-base)/b+1;r=mod(first_id-base,b)+1
        a=first(first_id);op_id=operation(dest)
        do j=1,nrep
          global_k=map%first(rep_start+j-1);source_k=kperm(global_k,op_id)
          local_k=source_k-starts(rank+1)+1
          do g=1,ncols
            source_g=colbase+g-1
            i=forward(source_g,op_id)
            delta=wrap(:,source_g,op_id)-wrap(:,a,op_id)
            if(any(abs(delta)>1))error stop 'EXX spinor orbit wrap exceeds one cell'
            ip=1+delta(1)+1+3*(delta(2)+1+3*(delta(3)+1))
            row(g,worker)=rb(r,g,local_k,p)*phase(ip,j)*gauge(dest,op_id)*conjg(gauge(i,op_id))
          enddo
          local_k=global_k-starts(rank+1)+1
          do spin_out=1,2;do spin_target=1,2
            coefficient=spin_rotation(spin_out,spin_left,op_id)*conjg(spin_rotation(spin_target,spin_right,op_id))
            if(coefficient==zero)cycle
            do g=1,ncols
              i=forward(colbase+g-1,op_id)
              t(g,:,worker)=target(i+(spin_target-1)*ng,:,j)*op%phase(i,local_k)
            enddo
            call zgemm('N','N',1,nt,ncols,-coefficient/real(nk,c_double),row(1,worker),1, &
              t(1,1,worker),chunk,zero,contracted(1,worker),1)
            action(dest+(spin_out-1)*ng,:,j)=action(dest+(spin_out-1)*ng,:,j)+ &
              contracted(:,worker)*conjg(op%phase(dest,local_k))
          enddo;enddo
        enddo
      enddo
      !$omp end do nowait
#if defined(HAVE_EXX_FUJITSU_THREADS) || defined(HAVE_EXX_NVPL_THREADS)
      call exx_blas_thread_control(0,blas_state)
#endif
      !$omp end parallel
#ifdef HAVE_EXX_OPENBLAS_THREADS
      call exx_blas_thread_control(0,blas_state)
#endif
      if(op%profile)call mark(6)
      enddo
      enddo
    enddo
    ierr=0
  contains
    subroutine mark(stage)
      integer,intent(in) :: stage
      real(c_double) :: now
      now=clock_now();op%seconds(stage)=op%seconds(stage)+now-stamp;stamp=now
    end subroutine
  end subroutine
  real(c_double) function clock_now()
    use iso_fortran_env, only: int64
    integer(int64) :: tick,rate
    call system_clock(tick,rate)
    clock_now=real(tick,c_double)/real(rate,c_double)
  end function
end module
