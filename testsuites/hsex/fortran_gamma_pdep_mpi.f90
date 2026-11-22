program test_gamma_pdep_mpi
  use mpi
  use omp_lib, only: omp_get_num_threads
  use hsex_gamma_pdep, only: gamma_pdep_screen
  use hsex_longwave, only: longwave_band,longwave_cell
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  integer,parameter :: ng=64,nv=1,nb=6,nopt=9
  complex(8) :: occupied(ng,1),velocity(ng,3),body(ng,nb),chi(nopt,nopt),raw(nopt,nopt)
  complex(8) :: kernel(1,1),expected(1,1),point(1,1)
  real(8) :: cell(3,3),rec(3,3),levels(1),macro(3,3),wanted(3,3),position(3),scale,signs(nopt),tol
  integer :: g(3,nb),reverse(nopt),rank,np,provided,ierr,first,last,i,j,ix,iy,iz,slot,code,requested,limit,hcalls
  character(256) :: path,mode,read_path,write_path
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  call get_command_argument(2,path);call get_command_argument(3,mode)
!$omp parallel
!$omp master
  if(rank==0)write(*,'(a,i0)')'HSEX_W_OMP team=',omp_get_num_threads()
!$omp end master
!$omp end parallel
  cell=0d0;rec=0d0
  do i=1,3
    cell(i,i)=10d0;rec(i,i)=2d0*acos(-1d0)/10d0
  enddo
  occupied=1d0/8d0;levels=0d0;scale=1d0
  if(mode=='changed')scale=1.7d0
  slot=0
  do ix=-1,1
    do iy=-1,1
      do iz=-1,1
        if(ix*ix+iy*iy+iz*iz/=1)cycle
        slot=slot+1;g(:,slot)=[ix,iy,iz]
      enddo
    enddo
  enddo
  reverse(:3)=[1,2,3];signs=[-1d0,-1d0,-1d0,1d0,1d0,1d0,1d0,1d0,1d0]
  do j=1,nb
    do i=1,nb
      if(all(g(:,i)==-g(:,j)))reverse(j+3)=i+3
    enddo
    do i=0,ng-1
      position=dble([mod(i,4),mod(i/4,4),i/16])/4d0
      body(i+1,j)=occupied(i+1,1)*exp(cmplx(0d0,2d0*acos(-1d0)*sum(position*g(:,j)),8))
    enddo
  enddo
  hcalls=0
  if(mode=='reference'.or.mode=='cold'.or.mode=='warm'.or.mode=='changed')then
  call velocity_action(occupied(:,1),velocity,code)
  call longwave_band(action,occupied,0d0,velocity,body,raw,code)
  if(code/=0)error stop 'Gamma independent response fixture failed'
  do j=1,nopt
    do i=1,nopt
      chi(i,j)=raw(i,j)+signs(i)*signs(j)*conjg(raw(reverse(i),reverse(j)))
    enddo
  enddo
  call longwave_cell([4,4,4],1,g,rec,1000d0,1,chi,wanted,point,expected,code)
  if(code/=0)error stop 'Gamma independent optical cell failed'
  endif
  requested=nopt;limit=40;tol=1d-8;read_path='none';write_path='none'
  if(mode=='reference')requested=-1
  if(mode=='cold')write_path=path
  if(mode=='warm'.or.mode=='changed'.or.mode=='bad-cache'.or.mode=='bad-cache-rank'.or. &
    mode=='bad-cache-cell')read_path=path
  if(mode=='bad-cache-rank')requested=8
  if(mode=='bad-cache-cell')cell(1,1)=10.1d0
  if(mode=='write-exists')write_path=path
  if(mode=='maxiter')then
    requested=1;limit=1;tol=1d-14
  endif
  if(mode=='nan-state'.and.rank==0)occupied(ng,1)=ieee_value(0d0,ieee_quiet_nan)
  if(mode=='inf-state'.and.rank==0)occupied(ng,1)=ieee_value(0d0,ieee_positive_inf)
  if(mode=='bad-gram'.and.rank==0)occupied=2d0*occupied
  first=1+nv*rank/np;last=nv*(rank+1)/np;hcalls=0
  call gamma_pdep_screen([4,4,4],1,.25d0,cell,occupied,levels,first,last,MPI_COMM_WORLD,rank, &
    requested,tol,limit,action,velocity_action,macro,kernel,code,trim(read_path),trim(write_path))
  if(mode=='reference'.or.mode=='cold'.or.mode=='warm'.or.mode=='changed')then
    if(code/=0)error stop 'Gamma valid model rejected'
    if(maxval(abs(macro-wanted))>1d-10.or.maxval(abs(kernel-expected))>1d-10) &
      error stop 'Gamma FFT head/wings differ from independent real-space vertices'
  else
    if(code==0)error stop 'Gamma invalid case accepted'
  endif
  write(*,'(a,i0,a,a)')'GAMMA_PDEP_MPI_PASS rank=',rank,' case=',trim(mode)
  call MPI_Finalize(ierr)
contains
  subroutine action(x,y,status)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: status
    hcalls=hcalls+1
    y=scale*(x-sum(x)/64d0);status=0
    if(mode=='fail-cg'.and.rank==np-1.and.hcalls>15)status=1
  end subroutine
  subroutine velocity_action(x,y,status)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:,:)
    integer,intent(out) :: status
    real(8) :: f(64),pos(3)
    integer :: a,site
    status=0
    do a=1,3
      do site=0,63
        pos=dble([mod(site,4),mod(site/4,4),site/16])/4d0
        f(site+1)=sqrt(2d0)*sin(2d0*acos(-1d0)*pos(a))/8d0
      enddo
      y(:,a)=cmplx(0d0,dble(a)*scale,8)*(f*sum(x)/8d0-sum(f*x)/8d0)
    enddo
  end subroutine
end program
