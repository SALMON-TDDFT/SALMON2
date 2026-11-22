program test_density_gamma_full
  use mpi
  use hsex_frozen_full, only: frozen_full_w
  use hsex_screen_math, only: screen_sites
  use hsex_longwave, only: longwave_cell
  implicit none
  integer,parameter :: n(3)=[4,4,4],ng=64,mesh(3)=[8,8,8],nk=512
  complex(8),allocatable :: psi(:,:,:)
  complex(8) :: optical(3,3),point(1,1),average(1,1),form,phase
  real(8) :: cell(3,3),rec(3,3),k(3,nk),coarse(nk),expected(nk),position(3),q(3),radius
  real(8) :: volume,pi,macro(3,3),matrix(nk),total(nk),w(1,nk)
  character(16) :: mode
  integer :: rank,np,ierr,ks,ke,count,first,status,i,j,a,b,c,site,sites(3,nk),g0(3,0),calls
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  pi=acos(-1d0);cell=0d0;rec=0d0
  do i=1,3
    cell(i,i)=12d0;rec(i,i)=2d0*pi/12d0
  enddo
  volume=12d0**3
  call screen_sites(mesh,sites)
  do i=1,nk
    k(:,i)=matmul(rec,(dble(modulo(sites(:,i),mesh))+.5d0)/mesh-.5d0)
  enddo
  ks=1+nk*rank/np;ke=nk*(rank+1)/np
  allocate(psi(ng,1,max(0,ke-ks+1)));psi=0d0
  if(size(psi,3)>0)psi(1,1,:)=1d0
  first=rank+1;count=merge(1,0,rank==0);calls=0
  mode='legacy'
  if(command_argument_count()>0)call get_command_argument(1,mode)
  if(mode=='hermitian')then
    call frozen_full_w(n,mesh,mesh,1,.15d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
      first,count,ng,1d-8,.01d0,.05d0,no_action,w(:count,:),status,density_gamma=.true., &
      prepare=prepare,nproc_q=min(2,np),closure='hermitian')
  else if(mode=='pairs')then
    call frozen_full_w(n,mesh,mesh,1,.15d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
      first,count,ng,1d-8,.01d0,.05d0,no_action,w(:count,:),status,density_gamma=.true., &
      prepare=prepare,nproc_q=min(2,np))
  else
    call frozen_full_w(n,mesh,mesh,1,.15d0,cell,k,psi,1,ks,ke,MPI_COMM_WORLD,rank,np, &
      first,count,ng,1d-8,.01d0,.05d0,no_action,w(:count,:),status,density_gamma=.true.)
  endif
  if(status/=0.or.calls/=0)error stop 'complete density Gamma zero-response pipeline failed'
  matrix=0d0
  if(rank==0)matrix=w(1,:)
  call MPI_Allreduce(matrix,total,nk,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(rank==0)then
    ! Independent bare-Coulomb reference, plus analytic spherical Gamma cell.
    optical=0d0
    call longwave_cell(n,1,g0,rec,volume,nk,optical,macro,point,average,status)
    if(status/=0)error stop 'vacuum longwave reference failed'
    coarse=0d0;coarse(1)=real(average(1,1),8)
    do i=2,nk
      do a=-1,1
        do b=-1,1
          do c=-1,1
            q=[a,b,c]+dble(sites(:,i))/mesh;radius=sum(matmul(rec,q)**2)
            if(radius/2d0>.15d0)cycle
            form=0d0
            do site=0,ng-1
              position=dble([mod(site,4),mod(site/4,4),site/16])/4d0
              form=form+exp(cmplx(0d0,2d0*pi*sum(q*position),8))/ng
            enddo
            coarse(i)=coarse(i)+4d0*pi*abs(form)**2/(volume*radius)
          enddo
        enddo
      enddo
    enddo
    expected=0d0
    do j=1,nk
      do i=1,nk
        phase=exp(cmplx(0d0,2d0*pi*sum(dble(sites(:,i)*sites(:,j))/mesh),8))
        expected(j)=expected(j)+real(phase,8)*coarse(i)/nk
      enddo
    enddo
    if(maxval(abs(total-expected))>.05d0*coarse(1)/nk)error stop 'density Gamma pipeline reference mismatch'
    print '(a,2es20.10)','DENSITY_GAMMA_FULL_ERROR max/allowed ', &
      maxval(abs(total-expected)),.05d0*coarse(1)/nk
    print *,'DENSITY_GAMMA_FULL_BARE_W_NO_DOUBLE_HEAD_PASS'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine prepare(source,target,worker,source_owner,target_owner,code)
    implicit none
    integer,intent(in) :: source,target,worker,source_owner,target_owner
    integer,intent(out) :: code
    code=0
  end subroutine
  subroutine no_action(index,x,y,code)
    implicit none
    integer,intent(in) :: index
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    calls=calls+1;y=0d0;code=1
  end subroutine
end program
