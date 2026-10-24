! Two-rank x-domain regression against an independently evaluated full-grid path.
program test_rspace
  use mpi
  use structures
  use nc_mgga23
  use sendrecv_grid, only: init_sendrecv_grid,dealloc_cache
  implicit none
  integer,parameter :: ng=16
  type(s_dft_system) :: sys,distributed
  type(s_rgrid) :: full,local
  type(s_parallel_info) :: refinfo,info
  type(s_stencil) :: st
  type(s_orbital) :: psi,hpsi,p,hp
  type(s_pp_nlcc) :: core,lcore
  type(s_sendrecv_grid) :: rs,ro,ls,lo
  real(8) :: eref,elocal,err,globalerr,theta,banderr
  integer :: rank,nproc,i,x,y,z,io,ix,ierr,neighbors(2,3)
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=2) stop 'real-space test requires two ranks'
  full%is=1;full%ie=ng;full%num=ng;full%is_array=1;full%ie_array=ng;full%nd=4
  full%is_overlap=-3;full%ie_overlap=ng+4
  allocate(full%idx(-3:ng+4),full%idy(-3:ng+4),full%idz(-3:ng+4))
  do i=-3,ng+4
    full%idx(i)=modulo(i-1,ng)+1;full%idy(i)=full%idx(i);full%idz(i)=full%idx(i)
  end do
  refinfo%im_s=1;refinfo%im_e=1;refinfo%ik_s=1;refinfo%ik_e=1;refinfo%io_s=1;refinfo%io_e=2
  refinfo%numo=2;refinfo%numk=1;refinfo%numm=1
  refinfo%if_divide_rspace=.false.;refinfo%icomm_r=MPI_COMM_SELF;refinfo%icomm_ko=MPI_COMM_SELF
  sys%nspin=2;sys%hvol=0.31d0
  allocate(sys%vec_k(3,1),sys%wtk(1),sys%rocc(2,1,2))
  sys%vec_k(:,1)=[0.17d0,-0.13d0,0.09d0];sys%wtk=1d0
  sys%rocc(:,1,1)=[0.8d0,0.57d0];sys%rocc(:,1,2)=sys%rocc(:,1,1)
  sys%rmatrix_B=reshape([1d0,0.12d0,0.04d0,0d0,1.1d0,-0.08d0,0.03d0,0d0,0.95d0],[3,3])
  st%if_orthogonal=.true.
  do i=1,3
    st%coef_nab(:,i)=[4d0/5d0,-1d0/5d0,4d0/105d0,-1d0/280d0]
    st%coef_lap(:,i)=[8d0/5d0,-1d0/5d0,8d0/315d0,-1d0/560d0]
  end do
  st%coef_lap0=sum(st%coef_lap)
  allocate(psi%zwf(ng,ng,ng,2,2,1,1),core%rho_nlcc(ng,ng,ng));core%rho_nlcc=0d0
  do z=1,ng
  do y=1,ng
  do x=1,ng
    theta=2d0*acos(-1d0)/ng
    do io=1,2
      psi%zwf(x,y,z,1,io,1,1)=cmplx(0.34d0+0.03d0*cos(theta*x*io),0.02d0*sin(theta*y),8)
      psi%zwf(x,y,z,2,io,1,1)=cmplx((-1d0)**io*0.19d0+0.02d0*sin(theta*y),0.03d0*cos(theta*z*io),8)
    end do
  end do
  end do
  end do
  hpsi=psi;hpsi%zwf=0d0
  call nc_mgga23_refresh(sys,full,refinfo,psi,st,rs,ro,core,eref)
  call nc_mgga23_add_action(sys,full,refinfo,psi,hpsi,st,ro)
  local=full
  local%is(1)=rank*ng/2+1;local%ie(1)=(rank+1)*ng/2
  local%num=local%ie-local%is+1
  local%is_array=local%is;local%ie_array=local%ie
  local%is_array(1)=local%is(1)-4;local%ie_array(1)=local%ie(1)+4
  local%is_overlap=local%is-4;local%ie_overlap=local%ie+4
  deallocate(local%idx)
  allocate(local%idx(local%is(1)-4:local%ie(1)+4))
  do i=local%is(1)-4,local%ie(1)+4
    local%idx(i)=i
  end do
  info=refinfo;info%if_divide_rspace=.true.;info%icomm_r=MPI_COMM_WORLD
  neighbors=MPI_PROC_NULL;neighbors(:,1)=1-rank
  call init_sendrecv_grid(ls,local,1,MPI_COMM_WORLD,neighbors)
  call init_sendrecv_grid(lo,local,4,MPI_COMM_WORLD,neighbors)
  allocate(p%zwf(local%is_array(1):local%ie_array(1),1:ng,1:ng,2,2,1,1))
  p%zwf=0d0
  do x=local%is(1),local%ie(1)
    p%zwf(x,:,:,:,:,:,:)=psi%zwf(x,:,:,:,:,:,:)
  end do
  hp=p;hp%zwf=0d0
  allocate(lcore%rho_nlcc(local%is(1):local%ie(1),1:ng,1:ng));lcore%rho_nlcc=0d0
  distributed=sys
  call nc_mgga23_refresh(distributed,local,info,p,st,ls,lo,lcore,elocal)
  call nc_mgga23_add_action(distributed,local,info,p,hp,st,lo)
  err=maxval(abs(hp%zwf(local%is(1):local%ie(1),:,:,:,:,:,:) &
    -hpsi%zwf(local%is(1):local%ie(1),:,:,:,:,:,:)))
  call MPI_Allreduce(err,globalerr,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  banderr=abs(distributed%nc_mgga23_payload%band_expectation-sys%nc_mgga23_payload%band_expectation)
  if(globalerr>2d-12.or.abs(elocal-eref)>2d-11.or.banderr>2d-11) then
    print *, 'real-space failure',rank,globalerr,elocal,eref,banderr
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end if
  if(rank==0) print *, 'real-space halo PASS',globalerr,abs(elocal-eref),banderr
  call dealloc_cache(ls);call dealloc_cache(lo)
  call MPI_Finalize(ierr)
end program
