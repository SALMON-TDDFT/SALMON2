#include "config.h"
! Occupied-space gauge for real-space HSE sources. No LCFO basis/projection.
module hse_grid_wannier
#ifdef USE_MPI
  use mpi
#endif
  use lcfo_seed, only: lcfo_seed_gamma
  use lcfo_mlwf_links, only: lcfo_initial_links
  use lcfo_dist_dense, only: lcfo_distributed_polar
  use hse_wannier_gauge, only: gauge_minimize_gamma_inplace
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: s_hse_grid_wannier,grid_wannier_source,grid_wannier_stage,grid_wannier_reset
  type s_hse_grid_wannier
    private
    logical :: initialized=.false.,enabled=.false.,in_step=.false.
    integer :: step=0,interval=1,comm=0
    real(8) :: radius=0d0,dv=0d0,lengths(3)=0d0
    complex(8),allocatable :: u(:,:),anchor(:,:),start_u(:,:),start_anchor(:,:)
    real(8),allocatable :: centers(:,:)
    logical,allocatable :: protected(:)
  end type
contains
  subroutine grid_wannier_reset(state)
    implicit none
    type(s_hse_grid_wannier),intent(out) :: state
    ! INTENT(OUT) releases all old allocatable components and applies defaults.
  end subroutine

  subroutine grid_wannier_stage(state,stage)
    implicit none
    type(s_hse_grid_wannier),intent(inout) :: state
    integer,intent(in) :: stage
    if(.not.state%initialized)error stop 'Grid MLWF: stage before source initialization'
    select case(stage)
    case(0)
      if(state%in_step)error stop 'Grid MLWF: nested physical step'
      state%step=state%step+1;state%in_step=.true.
      if(state%enabled)then
        state%start_u=state%u;state%start_anchor=state%anchor
      endif
    case(1)
      if(.not.state%in_step)error stop 'Grid MLWF: rollback outside physical step'
      if(state%enabled)then
        state%u=state%start_u;state%anchor=state%start_anchor
      endif
    case(2)
      if(.not.state%in_step)error stop 'Grid MLWF: end outside physical step'
      state%in_step=.false.
      ! Accepted source refresh occurs AFTER stage 2. Do not commit predictor U.
    case default
      error stop 'Grid MLWF: unknown Taylor stage'
    end select
  end subroutine

  subroutine grid_wannier_source(state,c,position,length,dv,comm,enabled,radius,maxiter,tolerance, &
                                  u_interval,seed_distributed,source)
    implicit none
    type(s_hse_grid_wannier),intent(inout) :: state
    complex(8),intent(in),contiguous :: c(:,:)
    real(8),intent(in) :: position(:,:),length(3),dv,radius,tolerance
    integer,intent(in) :: comm,maxiter,u_interval
    logical,intent(in) :: enabled,seed_distributed
    complex(8),allocatable,intent(out) :: source(:,:)
    complex(8),allocatable :: weighted(:,:),trial(:,:)
    real(8) :: minimum,displacement(3),r2
    integer :: i,j,a,status,bad,total_bad,ierr
    logical :: transport,first_source
    first_source=.not.state%initialized
    bad=0
    if(size(c,2)<1.or.any(shape(position)/=[3,size(c,1)]))bad=1
    if(.not.ieee_is_finite(dv))bad=1
    if(.not.ieee_is_finite(radius))bad=1
    if(.not.ieee_is_finite(tolerance))bad=1
    if(dv<=0d0.or.radius<0d0.or.tolerance<=0d0.or.maxiter<1.or.u_interval<1)bad=1
    if(.not.enabled.and.(radius>0d0.or.u_interval/=1))bad=1
    do a=1,3
      if(.not.ieee_is_finite(length(a)))bad=1
      if(length(a)<=0d0)bad=1
    enddo
    do j=1,size(c,2);do i=1,size(c,1)
      if(.not.ieee_is_finite(real(c(i,j))))bad=1
      if(.not.ieee_is_finite(aimag(c(i,j))))bad=1
    enddo;enddo
    do j=1,size(position,2);do a=1,size(position,1)
      if(.not.ieee_is_finite(position(a,j)))bad=1
    enddo;enddo
    total_bad=bad
#ifdef USE_MPI
    call MPI_Allreduce(bad,total_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
#endif
    if(total_bad/=0)error stop 'Grid MLWF: invalid source, geometry or controls'
    if(state%initialized)then
      if(state%enabled.neqv.enabled)error stop 'Grid MLWF: cannot change mode after initialization'
      if(state%radius/=radius.or.state%dv/=dv.or.any(state%lengths/=length).or. &
         state%interval/=u_interval.or.state%comm/=comm)error stop 'Grid MLWF: changed controls or geometry'
      if(enabled)then
        if(any(shape(state%anchor)/=shape(c)))error stop 'Grid MLWF: changed wavefunction shape'
      endif
    else
      state%enabled=enabled;state%radius=radius;state%dv=dv;state%lengths=length
      state%interval=u_interval;state%comm=comm
      if(enabled)call initialize_grid_gauge(state,c,position,maxiter,tolerance,seed_distributed)
      state%initialized=.true.
      if(enabled)then
        source=matmul(c,state%u)
        call report_coverage(state,source,position)
      endif
    endif
    if(.not.enabled)then
      source=c
      return
    endif
    ! The step-start anchor is immutable across all predictor/corrector calls.
    ! Held-U frames never become transport anchors (avoids secular phase drift).
    transport=state%step==0.or.mod(max(0,state%step-1),state%interval)==0
    if(transport.and..not.first_source)then
      weighted=sqrt(dv)*c;allocate(trial(size(c,2),size(c,2)))
      if(state%step>0)then
        call lcfo_distributed_polar(weighted,state%start_anchor,comm,trial,minimum,status)
      else
        call lcfo_distributed_polar(weighted,state%anchor,comm,trial,minimum,status)
      endif
      if(status/=0)error stop 'Grid MLWF: singular temporal overlap'
      call check_rotation(trial)
      state%u=trial;state%anchor=matmul(weighted,trial)
    endif
    source=matmul(c,state%u)
    if(radius==0d0.or.radius>=.5d0*sqrt(sum(length**2)))return
    r2=radius**2
    do j=1,size(c,2)
      if(state%protected(j))cycle
      do i=1,size(c,1)
        displacement=position(:,i)-state%centers(:,j)
        displacement=displacement-length*anint(displacement/length)
        if(sum(displacement**2)>r2)source(i,j)=0d0
      enddo
    enddo
  end subroutine

  subroutine initialize_grid_gauge(state,c,position,maxiter,tolerance,distributed)
    implicit none
    type(s_hse_grid_wannier),intent(inout) :: state
    complex(8),intent(in),contiguous :: c(:,:)
    real(8),intent(in) :: position(:,:),tolerance
    integer,intent(in) :: maxiter
    logical,intent(in) :: distributed
    complex(8),allocatable :: weighted(:,:),u(:,:,:),links(:,:,:,:)
    integer,allocatable :: counts(:)
    real(8) :: b(3,6),weights(6),delta,spread,gradient
    integer :: rank,np,ierr,ng,n,a,j,seed_status,status,iterations
    rank=0;np=1;ng=size(c,1);n=size(c,2)
#ifdef USE_MPI
    call MPI_Comm_rank(state%comm,rank,ierr);call MPI_Comm_size(state%comm,np,ierr)
#endif
    allocate(counts(np));counts=ng
#ifdef USE_MPI
    call MPI_Allgather(ng,1,MPI_INTEGER,counts,1,MPI_INTEGER,state%comm,ierr)
#endif
    weighted=sqrt(state%dv)*c;allocate(u(n,n,1))
    call lcfo_seed_gamma(weighted,counts,state%comm,u(:,:,1),seed_status,distributed=distributed)
    b=0d0
    do a=1,3
      delta=2d0*acos(-1d0)/state%lengths(a);b(a,a)=delta;b(a,a+3)=-delta
      weights(a)=1d0/(2d0*delta**2);weights(a+3)=weights(a)
    enddo
    call lcfo_initial_links(c,position,state%lengths,state%dv,state%comm,links)
    if(rank==0)then
      if(seed_status/=0)then
        u=0d0
        do j=1,n;u(j,j,1)=1d0;enddo
      endif
      call gauge_minimize_gamma_inplace(u,links,b,weights,maxiter,tolerance,spread,gradient,iterations,status)
      write(*,'(a,3i7,2es17.8)')'Grid MLWF initial iterations/status/seed/spread/gradient:', &
        iterations,status,seed_status,spread,gradient
    endif
    deallocate(links)
#ifdef USE_MPI
    call MPI_Bcast(status,1,MPI_INTEGER,0,state%comm,ierr)
    call MPI_Bcast(u,size(u),MPI_DOUBLE_COMPLEX,0,state%comm,ierr)
#endif
    if(status/=0.and.state%radius>0d0.and.state%radius<.5d0*sqrt(sum(state%lengths**2))) &
      error stop 'Grid MLWF: unconverged localization for finite support'
    call check_rotation(u(:,:,1))
    state%u=u(:,:,1);state%anchor=matmul(weighted,state%u)
    call initial_centers(state,matmul(c,state%u),position)
  end subroutine

  subroutine check_rotation(u)
    implicit none
    complex(8),intent(in) :: u(:,:)
    complex(8),allocatable :: overlap(:,:)
    integer :: i,j
    do j=1,size(u,2);do i=1,size(u,1)
      if(.not.ieee_is_finite(real(u(i,j))))error stop 'Grid MLWF: nonfinite rotation'
      if(.not.ieee_is_finite(aimag(u(i,j))))error stop 'Grid MLWF: nonfinite rotation'
    enddo;enddo
    overlap=matmul(conjg(transpose(u)),u)
    do j=1,size(u,2);overlap(j,j)=overlap(j,j)-1d0;enddo
    if(maxval(abs(overlap))>1d-10)error stop 'Grid MLWF: nonunitary rotation'
  end subroutine

  subroutine initial_centers(state,wf,position)
    implicit none
    type(s_hse_grid_wannier),intent(inout) :: state
    complex(8),intent(in) :: wf(:,:)
    real(8),intent(in) :: position(:,:)
    complex(8),allocatable :: local_moment(:,:),moment(:,:)
    real(8),allocatable :: local_norm(:),norm(:)
    real(8) :: pi
    integer :: n,j,a,ierr
    n=size(wf,2);pi=acos(-1d0)
    allocate(local_moment(3,n),moment(3,n),local_norm(n),norm(n))
    do j=1,n
      local_norm(j)=sum(abs(wf(:,j))**2)*state%dv
      do a=1,3
        local_moment(a,j)=sum(abs(wf(:,j))**2* &
          exp(cmplx(0d0,2d0*pi*position(a,:)/state%lengths(a),8)))*state%dv
      enddo
    enddo
    norm=local_norm;moment=local_moment
#ifdef USE_MPI
    call MPI_Allreduce(local_norm,norm,n,MPI_DOUBLE_PRECISION,MPI_SUM,state%comm,ierr)
    call MPI_Allreduce(local_moment,moment,3*n,MPI_DOUBLE_COMPLEX,MPI_SUM,state%comm,ierr)
#endif
    allocate(state%centers(3,n),state%protected(n));state%protected=.false.
    do j=1,n
      if(.not.ieee_is_finite(norm(j)))error stop 'Grid MLWF: nonfinite norm'
      if(norm(j)<=0d0)error stop 'Grid MLWF: empty occupied source'
      do a=1,3
        if(.not.ieee_is_finite(real(moment(a,j))))error stop 'Grid MLWF: nonfinite periodic moment'
        if(.not.ieee_is_finite(aimag(moment(a,j))))error stop 'Grid MLWF: nonfinite periodic moment'
        state%centers(a,j)=modulo(atan2(aimag(moment(a,j)),real(moment(a,j)))*state%lengths(a)/(2d0*pi), &
          state%lengths(a))
        if(abs(moment(a,j))/norm(j)<.1d0)state%protected(j)=.true.
      enddo
    enddo
    ! Initial centers stay fixed through propagation; source radius does not
    ! restrict propagation, and poorly defined periodic centers remain uncut.
  end subroutine

  subroutine report_coverage(state,wf,position)
    implicit none
    type(s_hse_grid_wannier),intent(in) :: state
    complex(8),intent(in) :: wf(:,:)
    real(8),intent(in) :: position(:,:)
    real(8),allocatable :: local(:,:),total(:,:)
    real(8) :: displacement(3),amount,fraction
    integer :: n,j,i,rank,ierr,iu,below
    n=size(wf,2);allocate(local(2,n),total(2,n));local=0d0;rank=0
    do j=1,n;do i=1,size(wf,1)
      amount=abs(wf(i,j))**2*state%dv;local(1,j)=local(1,j)+amount
      displacement=position(:,i)-state%centers(:,j)
      displacement=displacement-state%lengths*anint(displacement/state%lengths)
      if(state%radius==0d0.or.sum(displacement**2)<=state%radius**2)local(2,j)=local(2,j)+amount
    enddo;enddo
    total=local
#ifdef USE_MPI
    call MPI_Comm_rank(state%comm,rank,ierr)
    call MPI_Allreduce(local,total,2*n,MPI_DOUBLE_PRECISION,MPI_SUM,state%comm,ierr)
#endif
    if(rank/=0)return
    open(newunit=iu,file='grid_mlwf_radius.dat',status='replace')
    write(iu,'(a,es24.16)')'# Initial mesh WF geometric sphere coverage; radius_bohr (0=full): ',state%radius
    write(iu,'(a)')'# wf total_norm sphere_norm sphere_fraction protected_uncut'
    below=0
    do j=1,n
      fraction=total(2,j)/total(1,j)
      write(iu,'(i10,3es25.16,i4)')j,total(:,j),fraction,merge(1,0,state%protected(j))
      if(fraction<.999d0)below=below+1
    enddo
    close(iu)
    if(below>0)write(*,'(a,i8,a,i8)')'WARNING Grid MLWF radius: sphere norm below 99.9% for ',below,' of ',n
  end subroutine
end module
