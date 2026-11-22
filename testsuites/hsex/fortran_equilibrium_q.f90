program test_equilibrium_q
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  use hsex_equilibrium_q, only: equilibrium_q_response
  implicit none
  integer,parameter :: ng=64,nv=1,nk=2,nz=3
  complex(8) :: h(ng,ng,nk),vectors(ng,ng,nk),work(4*ng),z(nz),chi(nz,3),wanted(nz,3),vertex(ng)
  complex(8),allocatable :: psi(:,:,:),hocc(:,:,:)
  real(8) :: eigen(ng,nk),rw(3*ng),cell(3,3),kcart(3,nk),frac(2),residual,phase,angle,weight,delta
  integer :: rank,np,ierr,ks,ke,i,j,k,axis,site(3),other(3),index,sign,source,target,fold,trial,b,qtest
  integer :: points(3,3),qp(3),status
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  ks=1+nk*rank/np;ke=nk*(rank+1)/np
  allocate(psi(ng,1,max(0,ke-ks+1)),hocc(1,1,max(0,ke-ks+1)))
  cell=0d0
  do i=1,3
    cell(i,i)=12d0
  enddo
  frac=[-.25d0,.25d0];kcart=0d0;kcart(1,:)=frac*2d0*acos(-1d0)/12d0
  h=0d0
  do k=1,nk
    phase=kcart(1,k)*3d0
    do i=0,ng-1
      site=[mod(i,4),mod(i/4,4),i/16]
      h(i+1,i+1,k)=.6d0+.02d0*sum(cos(2d0*acos(-1d0)*(dble(site)+.17d0)/4d0))
      do axis=1,3
        other=site;other(axis)=mod(other(axis)+1,4);index=1+other(1)+4*other(2)+16*other(3)
        angle=0d0
        if(axis==1)angle=phase
        h(i+1,index,k)=h(i+1,index,k)-.1d0*exp(cmplx(0d0,angle,8))
        h(index,i+1,k)=conjg(h(i+1,index,k))
      enddo
    enddo
    vectors(:,:,k)=h(:,:,k)
    call zheev('V','L',ng,vectors(:,:,k),ng,eigen(:,k),work,size(work),rw,status)
    if(status/=0)error stop 'reference diagonalization failed'
  enddo
  do k=ks,ke
    psi(:,1,k-ks+1)=vectors(:,1,k);hocc(1,1,k-ks+1)=eigen(1,k)
  enddo
  z=[cmplx(0d0,.01d0,8),cmplx(.12d0,.01d0,8),cmplx(.2d0,.01d0,8)]
  do qtest=0,1
    qp=[qtest,0,0]
    ! A per-trial space of 16 suffices; a combined three-RHS space does not.
    call equilibrium_q_response([4,4,4],[2,1,1],qp,cell,kcart,psi,hocc,ks,ke, &
      MPI_COMM_WORLD,rank,np,16,1d-9,z,action,chi,residual,status)
    if(status/=0.or.residual>1d-9)error stop 'equilibrium q response rejected'
    points=0
    if(qtest==0)then
      points(1,1)=1;points(2,2)=1;points(3,3)=1
    else
      points(1,2)=1;points(2,3)=1
    endif
    wanted=0d0
    do sign=-1,1,2
      do source=1,nk
        target=source
        if(qtest/=0)target=3-source
        fold=nint(frac(source)+sign*qtest*.5d0-frac(target))
        do trial=1,3
          do i=0,ng-1
            site=[mod(i,4),mod(i/4,4),i/16]
            angle=2d0*acos(-1d0)*(sum(dble(site)*sign*points(:,trial))+site(1)*fold)/4d0
            vertex(i+1)=exp(cmplx(0d0,angle,8))*vectors(i+1,1,source)
          enddo
          do b=2,ng
            weight=abs(dot_product(vectors(:,b,target),vertex))**2/(12d0**3)
            delta=eigen(b,target)-eigen(1,source)
            if(sign==1)then
              wanted(:,trial)=wanted(:,trial)+weight/(z-delta)
            else
              wanted(:,trial)=wanted(:,trial)+weight*conjg(1d0/(-conjg(z)-delta))
            endif
          enddo
        enddo
      enddo
    enddo
    if(maxval(abs(chi-wanted))>1d-9)error stop 'equilibrium q spectral mismatch'
  enddo
  z(2)=cmplx(ieee_value(1d0,ieee_quiet_nan),.01d0,8)
  call equilibrium_q_response([4,4,4],[2,1,1],qp,cell,kcart,psi,hocc,ks,ke, &
    MPI_COMM_WORLD,rank,np,16,1d-9,z,action,chi,residual,status)
  if(status==0)error stop 'NaN frequency accepted'
  z(2)=cmplx(.12d0,.01d0,8)
  if(rank==np-1.and.ke>=ks)psi(1,1,1)=cmplx(ieee_value(1d0,ieee_positive_inf),0d0,8)
  call equilibrium_q_response([4,4,4],[2,1,1],qp,cell,kcart,psi,hocc,ks,ke, &
    MPI_COMM_WORLD,rank,np,16,1d-9,z,action,chi,residual,status)
  if(status==0)error stop 'remote infinite occupied state accepted'
  if(rank==0)print *, 'EQUILIBRIUM_Q_OCCUPIED_ONLY_DYNAMIC_FOLD_SPECTRAL_MPI_PASS'
  call MPI_Finalize(ierr)
contains
  subroutine action(source,target,x,y,status)
    implicit none
    integer,intent(in) :: source,target
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: status
    status=1
    if(source<ks.or.source>ke)return
    y=matmul(h(:,:,target),x);status=0
  end subroutine
end program
