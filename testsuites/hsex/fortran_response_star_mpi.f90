program test_response_star_mpi
  use mpi
  use omp_lib, only: omp_get_num_threads
  use hsex_reduced_q, only: reduced_q_test
  use hsex_screen_math, only: screen_sites
  implicit none
  integer,parameter :: ng=216,nk=64,nv=1
  complex(8) :: u(ng,ng),h(ng,ng)
  complex(8),allocatable :: states(:,:,:)
  real(8) :: levels(2,nk),cell(3,3),k(3,nk),frac(3,nk),r(3),e(ng),atoms(3,1)
  integer :: sites(3,nk),gsites(3,ng),rank,np,provided,ierr,ks,ke,i,j,code,destination,species(1)
  character(32) :: mode
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  call get_command_argument(2,mode)
!$omp parallel
!$omp master
  if(rank==0)write(*,'(a,i0)')'HSEX_W_OMP team=',omp_get_num_threads()
!$omp end master
!$omp end parallel
  cell=0d0;atoms=0d0;species=1
  if(mode=='broken-atoms')atoms(:,1)=[.13d0,.21d0,.37d0]
  do i=1,3
    cell(i,i)=10d0
  enddo
  call screen_sites([4,4,4],sites);call screen_sites([6,6,6],gsites)
  do j=1,ng
    do i=0,ng-1
      r=dble([mod(i,6),mod(i/6,6),i/36])/6d0
      u(i+1,j)=exp(cmplx(0d0,2d0*acos(-1d0)*sum(r*gsites(:,j)),8))/sqrt(dble(ng))
    enddo
    e(j)=1d0+sum(gsites(:,j)**2)
    if(mode=='broken-h')e(j)=e(j)+.1d0*gsites(1,j)**2
  enddo
  e(1)=-1d0;h=0d0
  do j=1,ng
    do i=1,ng
      h(:,i)=h(:,i)+e(j)*u(:,j)*conjg(u(i,j))
    enddo
  enddo
  do j=0,nk-1
    frac(:,j+1)=-.375d0+.25d0*[mod(j,4),mod(j/4,4),j/16]
    k(:,j+1)=2d0*acos(-1d0)/10d0*frac(:,j+1)
    levels(:,j+1)=[-1d0,2d0]+.001d0*sum(frac(:,j+1)**2)
  enddo
  ks=1+nk*rank/np;ke=nk*(rank+1)/np
  allocate(states(ng,nv,ks:ke))
  do j=ks,ke
    states(:,:,j)=u(:,:nv)
  enddo
  destination=64
  call reduced_q_test([6,6,6],[4,4,4],[4,4,4],.12d0,cell,k,levels,states,nv,ks,ke, &
    MPI_COMM_WORLD,rank,destination,6,64,1d-9,apply,code,.true.,atoms,species)
  if(mode=='fail')then
    if(code==0)error stop 'symmetry owner callback failure accepted'
  else
    if(code/=0)error stop 'symmetry shared response diagnostic failed'
  endif
  write(*,'(a,i0,a,a)')'RESPONSE_STAR_MPI_PASS rank=',rank,' case=',trim(mode)
  call MPI_Finalize(ierr)
contains
  subroutine apply(ik,x,y,status)
    implicit none
    integer,intent(in) :: ik
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: status
    y=matmul(h,x)+.001d0*sum(frac(:,ik)**2)*x;status=0
    if(ik<ks.or.ik>ke.or.mode=='fail')status=1
  end subroutine
end program
