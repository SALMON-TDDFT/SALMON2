program test_density_gamma
  use mpi
  use hsex_density_gamma, only: density_gamma_stencil,density_gamma_replace
  use hsex_screen_math, only: screen_sites
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  integer,parameter :: mesh(3)=[8,8,8],nq=512
  integer :: rank,np,ierr,status,na,first,count,i,j,a,b
  integer :: sites(3,nq)
  integer,allocatable :: selected(:)
  real(8) :: rec(3,3),q(3),q2,qc,pi,vol
  complex(8),allocatable :: coarse(:,:),before(:,:),expected(:,:)
  character(32) :: mode
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  mode='valid'
  if(command_argument_count()>0)call get_command_argument(1,mode)
  na=2
  if(mode=='empty')na=1
  first=1+na*na*rank/np;count=na*na*(rank+1)/np-first+1
  allocate(coarse(count,nq),before(count,nq),expected(na,na))
  pi=acos(-1d0);rec=0d0
  do i=1,3
    rec(i,i)=.4d0
  enddo
  vol=(2d0*pi)**3/.4d0**3;qc=(6d0*pi*pi/(vol*nq))**(1d0/3d0)
  call density_gamma_stencil(mesh,rec,selected,status)
  if(status/=0.or.size(selected)/=56)error stop 'nested stencil or tied-shell selection failed'
  call screen_sites(mesh,sites)
  coarse=0d0
  do i=1,nq
    q=matmul(rec,dble(sites(:,i))/mesh);q2=sum(q*q)
    do j=1,count
      a=mod(first+j-2,na)+1;b=(first+j-2)/na+1
      if(q2==0d0)then
        coarse(j,i)=1234d0
      else if(a==1.and.b==1)then
        coarse(j,i)=.2d0/q2+2d0
        if(mode=='metal')coarse(j,i)=2d0
        if(mode=='qfit')coarse(j,i)=coarse(j,i)+100d0*q2
        if(mode=='negative')coarse(j,i)=-.01d0/q2+2d0
      else if(a==2.and.b==2)then
        coarse(j,i)=3d0
      else
        coarse(j,i)=cmplx(0d0,merge(1d0,-1d0,a==1)*.01d0*q(1)/q2,8)
        if(mode=='metal')coarse(j,i)=0d0
      endif
    enddo
  enddo
  before=coarse
  if(mode=='nan'.and.count>0.and.rank==np-1)coarse(1,selected(1))=cmplx(ieee_value(1d0,ieee_quiet_nan),0d0,8)
  if(mode=='infinite'.and.count>0.and.rank==np-1) &
    coarse(1,selected(1))=cmplx(ieee_value(1d0,ieee_positive_inf),0d0,8)
  if(mode=='layout'.and.rank==np-1)first=1
  if(mode=='geometry')rec(:,3)=0d0
  if(mode=='short')then
    call density_gamma_replace(coarse(:,:8),first,na,[2,2,2],rec,vol,.05d0,MPI_COMM_WORLD,rank,status)
  else
    call density_gamma_replace(coarse,first,na,mesh,rec,vol,1d-8,MPI_COMM_WORLD,rank,status)
  endif
  if(mode=='nan'.or.mode=='infinite'.or.mode=='layout'.or.mode=='negative'.or.mode=='short'.or.mode=='qfit'.or.mode=='geometry')then
    if(status==0)error stop 'invalid Gamma candidate accepted'
    if(any(coarse(:,1)/=before(:,1)))error stop 'failed fit overwrote old Gamma cell'
    if(rank==0)print *,'DENSITY_GAMMA_REJECT_PRESERVES_INPUT_PASS'
  else
    if(status/=0)error stop 'density Gamma fit rejected'
    expected=0d0;expected(1,1)=.6d0/qc**2+2d0
    if(mode=='metal')expected(1,1)=2d0
    if(na==2)expected(2,2)=3d0
    do j=1,count
      a=mod(first+j-2,na)+1;b=(first+j-2)/na+1
      if(abs(coarse(j,1)-expected(a,b))>1d-8)error stop 'Gamma cell integral, wing or double-counting error'
    enddo
    if(any(coarse(:,2:)/=before(:,2:)))error stop 'finite q values changed'
    if(rank==0)print *,'DENSITY_GAMMA_DISTRIBUTED_CELL_PASS'
  endif
  call MPI_Finalize(ierr)
end program
