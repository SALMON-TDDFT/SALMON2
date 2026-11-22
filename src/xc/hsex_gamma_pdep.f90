! Gamma insulating screening in a real, time-reversal-adapted optical basis.
! Source bands are MPI distributed. No dense G² or grid-by-G storage for PDEP.
module hsex_gamma_pdep
  use, intrinsic :: iso_c_binding
  implicit none
  private
  include 'fftw3.f03'
  public :: gamma_pdep_screen
  abstract interface
    subroutine gamma_h(x,y,status)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
    subroutine gamma_velocity(x,y,status)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:,:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine gamma_pdep_screen(n,boxes,cutoff,cell,occupied,levels,band_first,band_last,comm,rank, &
      requested,tolerance,maxiter,apply,velocity,macro,kernel,status,cache_read,cache_write)
    use hsex_screen_math, only: screen_geometry,screen_finite
    use hsex_sternheimer, only: sternheimer_solve,s_sternheimer_stats
    use hsex_pdep, only: pdep_modes
    use hsex_longwave, only: longwave_cell,longwave_lowrank_cell
    use hsex_pdep_cache, only: pdep_cache_read,pdep_cache_write
    use communication, only: comm_get_max,comm_summation,comm_is_root,comm_bcast,comm_sync_all
    implicit none
    integer,intent(in) :: n(3),boxes,band_first,band_last,comm,rank,requested,maxiter
    real(8),intent(in) :: cutoff,cell(3,3),levels(:),tolerance
    complex(8),intent(in) :: occupied(:,:)
    procedure(gamma_h) :: apply
    procedure(gamma_velocity) :: velocity
    real(8),intent(out) :: macro(3,3)
    complex(8),intent(out) :: kernel(:,:)
    integer,intent(out) :: status
    character(*),intent(in) :: cache_read,cache_write
    type(s_sternheimer_stats) :: statistics
    type(c_ptr) :: forward,backward
    complex(c_double_complex),allocatable,target :: fft(:)
    complex(8),allocatable :: optical(:,:,:),v(:,:),rhs(:),solution(:),potential(:),density(:),coeff(:),raw(:)
    complex(8),allocatable :: basis(:,:),initial(:,:),physical(:,:),matrix(:,:),chi(:,:),point(:,:),gram(:,:)
    real(8),allocatable :: weight(:),values(:),residual(:)
    integer,allocatable :: g(:,:),slots(:),reverse(:)
    real(8) :: rec(3,3),vol,r(3),pi,times(4),time_max(4),counts(3),count_sum(3),worst(2),worst_max(2)
    integer(8) :: started,prepared,built,finished,rate
    integer :: ng,nv,nb,nopt,rr,imax(3),ix,iy,iz,j,i,k,side,axis,local_code,code,calls,active,at,close_code
    logical :: cached
    call comm_sync_all(comm);call system_clock(started,rate)
    status=1;macro=0d0;kernel=0d0;statistics=s_sternheimer_stats();forward=c_null_ptr;backward=c_null_ptr
    ng=product(n);nv=size(levels);pi=acos(-1d0);local_code=0;calls=0
    if(any(n<1).or.boxes<1.or.nv<1.or.size(occupied,1)/=ng.or.size(occupied,2)/=nv)local_code=1
    if(.not.screen_finite(occupied).or..not.screen_finite(levels).or..not.screen_finite([cutoff,tolerance]))local_code=1
    if(.not.screen_finite(cell))local_code=1
    if(cutoff<=0d0.or.tolerance<=0d0.or.maxiter<1.or.requested< -1)local_code=1
    call comm_get_max(local_code,comm)
    if(local_code/=0)return
    if(any(mod(n,boxes)/=0).or.any(shape(kernel)/=[boxes**3,boxes**3]))return
    call screen_geometry(cell,rec,vol,code)
    if(code/=0)return
    gram=matmul(conjg(transpose(occupied)),occupied)
    do j=1,nv
      gram(j,j)=gram(j,j)-1d0
    enddo
    local_code=0
    if(.not.screen_finite(gram).or.maxval(abs(gram))>1d-8)local_code=1
    call comm_get_max(local_code,comm)
    if(local_code/=0)return
    deallocate(gram)
    imax=ceiling(sqrt(2d0*cutoff)*sqrt(sum(cell**2,dim=1))/(2d0*pi))+1;nb=0
    do ix=-imax(1),imax(1)
      do iy=-imax(2),imax(2)
        do iz=-imax(3),imax(3)
          if(all([ix,iy,iz]==0))cycle
          if(sum(matmul(rec,dble([ix,iy,iz]))**2)/2d0<=cutoff)nb=nb+1
        enddo
      enddo
    enddo
    nopt=nb+3;rr=requested
    if(rr==0)rr=min(nopt,2*nv)
    if(rr>nopt)return
    allocate(g(3,nb),slots(nb),reverse(nb),weight(nopt),stat=local_code)
    call comm_get_max(local_code,comm)
    if(local_code/=0)return
    i=0;weight(:3)=sqrt(4d0*pi)
    do ix=-imax(1),imax(1)
      do iy=-imax(2),imax(2)
        do iz=-imax(3),imax(3)
          if(all([ix,iy,iz]==0))cycle
          r=matmul(rec,dble([ix,iy,iz]))
          if(sum(r*r)/2d0>cutoff)cycle
          if(any(2*abs([ix,iy,iz])>=n))return
          i=i+1;g(:,i)=[ix,iy,iz];weight(i+3)=sqrt(4d0*pi/sum(r*r))
          slots(i)=1+modulo(ix,n(1))+n(1)*(modulo(iy,n(2))+n(2)*modulo(iz,n(3)))
        enddo
      enddo
    enddo
    do i=1,nb
      reverse(i)=0
      do j=1,nb
        if(all(g(:,j)==-g(:,i)))reverse(i)=j
      enddo
      if(reverse(i)==0)return
    enddo
    active=max(0,band_last-band_first+1)
    allocate(optical(ng,3,active),v(ng,3),rhs(ng),solution(ng),potential(ng),density(ng),coeff(nopt),raw(nopt),stat=local_code)
    call comm_get_max(local_code,comm)
    if(local_code/=0)return
    allocate(fft(ng),point(boxes**3,boxes**3),stat=local_code)
    call comm_get_max(local_code,comm)
    if(local_code/=0)return
    fft=0d0
    forward=fftw_plan_dft_3d(n(3),n(2),n(1),fft,fft,FFTW_FORWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    backward=fftw_plan_dft_3d(n(3),n(2),n(1),fft,fft,FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    local_code=0
    if(.not.c_associated(forward).or..not.c_associated(backward))local_code=1
    call comm_get_max(local_code,comm)
    if(local_code/=0)then
      call release_fft();return
    endif
    do j=band_first,band_last
      call apply(occupied(:,j),solution,code)
      if(code/=0.or..not.screen_finite(solution))then
        local_code=1;exit
      endif
      if(norm2(abs(solution-levels(j)*occupied(:,j)))>1d-7)then
        local_code=1;exit
      endif
      call velocity(occupied(:,j),v,code)
      if(code/=0.or..not.screen_finite(v))then
        local_code=1;exit
      endif
      do axis=1,3
        call sternheimer_solve(apply,occupied,levels(j),v(:,axis),optical(:,axis,j-band_first+1),code,statistics)
        if(code/=0)exit
      enddo
      if(code/=0)then
        local_code=1;exit
      endif
    enddo
    call comm_get_max(local_code,comm)
    if(local_code/=0)then
      call release_fft();return
    endif
    cached=cache_read/='none'
    if(cached.and.rr<0)then
      call release_fft();return
    endif
    if(cached)then
      allocate(initial(nopt,rr),stat=local_code)
      call comm_get_max(local_code,comm)
      if(local_code/=0)then
        call release_fft();return
      endif
      call pdep_cache_read(cache_read,n,g,cell,cutoff,nv,initial,comm,rank,code)
      if(code/=0)then
        call release_fft();return
      endif
      if(comm_is_root(rank))write(*,'(a)')'HSEX_GAMMA_PDEP_CACHE_READ_PASS'
    endif
    call system_clock(prepared)
    if(rr<0)then
      allocate(matrix(nopt,nopt),chi(nopt,nopt),basis(nopt,nopt),stat=local_code)
      call comm_get_max(local_code,comm)
      if(local_code/=0)then
        call release_fft();return
      endif
      matrix=0d0
      do i=1,nopt
        coeff=0d0;coeff(i)=1d0
        call polarization(coeff,matrix(:,i),code)
        if(code/=0)then
          call release_fft();return
        endif
      enddo
      calls=nopt
      ! L B_real L^H restores the original optical-head and complex-G basis.
      do i=1,nopt
        call to_physical(matrix(:,i),basis(:,i))
      enddo
      do i=1,nopt
        call to_physical(conjg(basis(i,:)),raw)
        chi(i,:)=conjg(raw)
      enddo
      do i=1,nopt
        chi(:,i)=chi(:,i)/(2d0/vol*weight*weight(i))
      enddo
      call system_clock(built)
      code=0
      if(comm_is_root(rank))call longwave_cell(n,boxes,g,rec,vol,1,chi,macro,point,kernel,code)
    else
      allocate(basis(nopt,rr),values(rr),residual(rr),physical(nopt,rr),stat=local_code)
      call comm_get_max(local_code,comm)
      if(local_code/=0)then
        call release_fft();return
      endif
      if(cached)then
        call pdep_modes(polarization,nopt,rr,tolerance,maxiter,basis,values,residual,calls,code,initial, &
          synchronize=collective_status)
      else
        call pdep_modes(polarization,nopt,rr,tolerance,maxiter,basis,values,residual,calls,code,synchronize=collective_status)
      endif
      if(code/=0)then
        if(comm_is_root(rank))write(*,'(a,i0,a,es13.5)')'HSEX_GAMMA_PDEP_REJECT status=',code, &
          ' Ritz residual=',maxval(residual)
        call release_fft();return
      endif
      do i=1,rr
        call to_physical(basis(:,i),physical(:,i))
      enddo
      call system_clock(built)
      if(comm_is_root(rank))call longwave_lowrank_cell(n,boxes,g,rec,vol,1,physical,values,macro,point,kernel,code)
    endif
    call comm_get_max(code,comm)
    if(code/=0)then
      call release_fft();return
    endif
    call comm_bcast(macro,comm,0);call comm_bcast(kernel,comm,0)
    local_code=0
    if(.not.screen_finite(kernel))local_code=1
    if(maxval(abs(aimag(kernel)))>1d-7.or.maxval(abs(kernel-transpose(kernel)))>1d-8)local_code=1
    call comm_get_max(local_code,comm)
    if(local_code/=0)then
      call release_fft();return
    endif
    if(rr>0.and.cache_write/='none')then
      call pdep_cache_write(cache_write,n,g,cell,cutoff,nv,basis,comm,rank,code)
      if(code/=0)then
        call release_fft();return
      endif
      if(comm_is_root(rank))write(*,'(a)')'HSEX_GAMMA_PDEP_CACHE_WRITE_PASS'
    endif
    call release_fft();call comm_sync_all(comm);call system_clock(finished)
    times=dble([prepared-started,built-prepared,finished-built,finished-started])/rate
    call comm_get_max(times,time_max,4,comm)
    counts=dble([statistics%rhs,statistics%iterations,statistics%h_actions])
    call comm_summation(counts,count_sum,3,comm)
    worst=[statistics%max_relative_residual,0d0]
    if(rr>0)worst(2)=maxval(residual)
    call comm_get_max(worst,worst_max,2,comm)
    if(comm_is_root(rank))then
      write(*,'(a,3i10)')'HSEX_GAMMA_PDEP dimensions G+optical/rank/actions ',nopt,rr,calls
      write(*,'(a,4es16.7)')'HSEX_GAMMA_W_SECONDS prepare/basis/box+cache/total ',time_max
      write(*,'(a,3f20.0)')'HSEX_GAMMA_W_SOLVES rhs/iterations/Hpsi ',count_sum
      write(*,'(a,2es16.7)')'HSEX_GAMMA_W_RESIDUAL CG/Ritz ',worst_max
      do i=1,3
        write(*,'(a,i0,3es19.10)')'HSEX_GAMMA_EPS row=',i,macro(i,:)
      enddo
      write(*,'(a)')'HSEX_GAMMA_SCREEN_PASS empty_band_sum=0'
    endif
    status=0
  contains
    subroutine collective_status(code)
      implicit none
      integer,intent(inout) :: code
      call comm_get_max(code,comm)
    end subroutine
    subroutine release_fft()
      implicit none
      if(c_associated(forward))call fftw_destroy_plan(forward)
      if(c_associated(backward))call fftw_destroy_plan(backward)
      forward=c_null_ptr;backward=c_null_ptr
    end subroutine
    subroutine to_physical(x,y)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer :: l,other
      y(:3)=cmplx(0d0,1d0,8)*x(:3)
      do l=1,nb
        other=reverse(l)
        if(l>other)cycle
        y(l+3)=(x(l+3)+cmplx(0d0,1d0,8)*x(other+3))/sqrt(2d0)
        y(other+3)=(x(l+3)-cmplx(0d0,1d0,8)*x(other+3))/sqrt(2d0)
      enddo
    end subroutine
    subroutine polarization(x,y,code)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: code
      complex(8) :: physical_x(nopt),total_raw(nopt)
      real(8) :: local(nopt),total(nopt)
      integer :: source,loc,l,other,a,local_error
      code=1;y=0d0
      if(.not.screen_finite(x).or.maxval(abs(aimag(x)))>1d-8)return
      call to_physical(x,physical_x)
      fft=0d0
      do l=1,nb
        fft(slots(l))=weight(l+3)*physical_x(l+3)
      enddo
      call fftw_execute_dft(backward,fft,fft);potential=fft
      local=0d0;local_error=0;density=0d0
      do source=band_first,band_last
        loc=source-band_first+1
        rhs=occupied(:,source)*potential+sqrt(4d0*pi)*matmul(optical(:,:,loc),physical_x(:3))
        call sternheimer_solve(apply,occupied,levels(source),rhs,solution,code,statistics)
        if(code/=0)then
          local_error=code;exit
        endif
        do a=1,3
          raw(a)=sqrt(4d0*pi)*dot_product(optical(:,a,loc),solution)
        enddo
        local(:3)=local(:3)+aimag(raw(:3))
        density=density+conjg(occupied(:,source))*solution
      enddo
      call comm_get_max(local_error,comm)
      if(local_error/=0)then
        code=local_error
        return
      endif
      ! Sum the occupied-band density before transforming: one forward FFT
      ! per rank/trial, rather than one per source band.
      if(active>0)then
        fft=density;call fftw_execute_dft(forward,fft,fft)
        do l=1,nb
          raw(l+3)=weight(l+3)*fft(slots(l))
        enddo
        do l=1,nb
          other=reverse(l)
          if(l>other)cycle
          local(l+3)=real(raw(l+3)+raw(other+3),8)/sqrt(2d0)
          local(other+3)=aimag(raw(l+3)-raw(other+3))/sqrt(2d0)
        enddo
      endif
      call comm_summation(local,total,nopt,comm)
      y=cmplx(4d0/vol*total,0d0,8);code=0
    end subroutine
  end subroutine gamma_pdep_screen
end module hsex_gamma_pdep
