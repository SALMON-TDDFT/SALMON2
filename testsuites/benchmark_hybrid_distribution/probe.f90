program probe
  use mpi
  use iso_c_binding, only: c_int64_t
  use hse_spatial
  use hse_ace, only: hse_ace_state
  use exx_orbitals
  use lcfo_scalapack
  implicit none
  interface
    function peak_rss_bytes() bind(C) result(bytes)
      import c_int64_t
      integer(c_int64_t) :: bytes
    end function
  end interface
  type(spatial_exx_state) :: exchange
  type(hse_ace_state) :: ace
  type(lcfo_dense_state) :: dense
  complex(8),allocatable :: psi(:,:,:),w(:,:,:),action(:,:,:),block(:,:)
  real(8),allocatable :: occupation(:,:)
  integer :: rank,np,err,comm_r,comm_o,comm(2),coords(2),dims(2),orb,rs,rr,ro
  integer :: kd(3),partner,ngrid,no,n(3),m(3),lo(3),first,last,nlocal,x,y,z,g,j,k,mode,nb,dim,status,code
  integer(c_int64_t) :: baseline,rss
  real(8) :: t,times(4),h(3),phase,pi,volume,energy,global_energy,error,global_error,omega,he,oe,re
  character(64) :: arg
  call MPI_Init(err)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,err);call MPI_Comm_size(MPI_COMM_WORLD,np,err)
  call get_command_argument(1,arg);read(arg,*,iostat=code)mode
  if(code/=0)call MPI_Abort(MPI_COMM_WORLD,1,err)
  baseline=peak_rss_bytes();times=0d0;global_energy=0d0;global_error=0d0
  if(mode==1)then
    call get_command_argument(2,arg);read(arg,*)ngrid
    call get_command_argument(3,arg);read(arg,*)no
    call get_command_argument(4,arg);read(arg,*)orb
    call get_command_argument(5,arg);read(arg,*)omega
    if(min(ngrid,no,orb)<1.or.mod(np,orb)/=0.or.orb>no.or.no>ngrid**3)call MPI_Abort(MPI_COMM_WORLD,2,err)
    rs=np/orb;rr=mod(rank,rs);ro=rank/rs
    dims=[int(sqrt(real(rs,8))),1]
    do while(mod(rs,dims(1))/=0)
      dims(1)=dims(1)-1
    enddo
    dims(2)=rs/dims(1)
    if(any(mod(ngrid,dims)/=0))call MPI_Abort(MPI_COMM_WORLD,3,err)
    coords=[mod(rr,dims(1)),rr/dims(1)]
    call MPI_Comm_split(MPI_COMM_WORLD,ro,rr,comm_r,err)
    call MPI_Comm_split(MPI_COMM_WORLD,rr,ro,comm_o,err)
    call MPI_Comm_split(comm_r,coords(2),coords(1),comm(1),err)
    call MPI_Comm_split(comm_r,coords(1),coords(2),comm(2),err)
    kd(1)=nint(real(no,8)**(1d0/3d0))
    do while(mod(no,kd(1))/=0)
      kd(1)=kd(1)-1
    enddo
    kd(2)=int(sqrt(real(no/kd(1),8)))
    do while(mod(no/kd(1),kd(2))/=0)
      kd(2)=kd(2)-1
    enddo
    kd(3)=no/kd(1)/kd(2)
    if(mod(no,2)/=0.or.any(kd<2).or.any(kd>ngrid))call MPI_Abort(MPI_COMM_WORLD,7,err)
    n=ngrid;h=1d0;m=[ngrid,ngrid/dims(1),ngrid/dims(2)];lo=[0,coords(1)*m(2),coords(2)*m(3)]
    first=no*ro/orb+1;last=no*(ro+1)/orb;nlocal=last-first+1
    allocate(psi(product(m),nlocal,1),w(product(m),nlocal,1),action(product(m),nlocal,1),occupation(nlocal,1))
    pi=acos(-1d0);volume=real(product(n),8);occupation=2d0
    do j=first,last
      g=0;k=j-1
      do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
        g=g+1
        partner=k+1
        if(mod(k,2)==1)partner=k-1
        psi(g,j-first+1,1)=cos(.17d0)*packet([x,y+lo(2),z+lo(3)],k) &
          +merge(1d0,-1d0,mod(k,2)==0)*sin(.17d0)*packet([x,y+lo(2),z+lo(3)],partner)
      enddo;enddo;enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,err);t=MPI_Wtime()
    call spatial_exx_refresh(exchange,n,h,dims,coords,comm,comm_r,psi,3,1d-7,status,occupation,comm_o)
    call require_success(status)
    times(1)=MPI_Wtime()-t
    call MPI_Barrier(MPI_COMM_WORLD,err);t=MPI_Wtime()
    call spatial_exx_apply(exchange,n,h,dims,coords,comm,comm_r,0d0,psi,w,status,omega,comm_o)
    call require_success(status)
    times(2)=MPI_Wtime()-t
    call MPI_Barrier(MPI_COMM_WORLD,err);t=MPI_Wtime()
    call orbital_ace_build(ace,psi,w,1d0,comm_r,comm_o,status)
    call require_success(status)
    times(3)=MPI_Wtime()-t
    call MPI_Barrier(MPI_COMM_WORLD,err);t=MPI_Wtime()
    call orbital_ace_apply(ace,psi,action,comm_r,comm_o,status)
    call require_success(status)
    times(4)=MPI_Wtime()-t
    error=maxval(abs(action-w))
    energy=real(sum(conjg(psi)*w),8)
    call MPI_Allreduce(error,global_error,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,err)
    call MPI_Allreduce(energy,global_energy,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,err)
    if(global_error>1d-10)call MPI_Abort(MPI_COMM_WORLD,4,err)
  else if(mode==2)then
    call get_command_argument(2,arg);read(arg,*)dim
    call MPI_Barrier(MPI_COMM_WORLD,err);t=MPI_Wtime()
    call lcfo_dense_init(dense,dim,MPI_COMM_WORLD,status)
    call require_success(status)
    ! No full input oracle: materialize one bounded-width block, then discard.
    nb=min(32,dim)
    allocate(block(dim,nb))
    do first=1,dim,nb
      last=min(dim,first+nb-1)
      do j=first,last;do g=1,dim
        block(g,j-first+1)=cmplx(cos(.017d0*(g+j)),sin(.013d0*(g-j)),8)/dim
        if(g==j)block(g,j-first+1)=block(g,j-first+1)+1d0+real(g,8)/dim
      enddo;enddo
      call lcfo_dense_add(dense,1,first,block(:,:last-first+1))
    enddo
    deallocate(block)
    times(1)=MPI_Wtime()-t
    call MPI_Barrier(MPI_COMM_WORLD,err);t=MPI_Wtime()
    call lcfo_dense_solve(dense,dim,he,oe,re,status)
    call require_success(status)
    times(2)=MPI_Wtime()-t
    global_error=max(he,oe,re);global_energy=sum(dense%values)
  else
    call MPI_Abort(MPI_COMM_WORLD,5,err)
  endif
  if(mode==1)write(*,'(a,3(1x,i0),2(1x,es26.16e3))') &
    'LOCALIZATION',rank,exchange%iterations,exchange%localization_status,exchange%spread,exchange%gradient
  rss=peak_rss_bytes()
  if(baseline<0.or.rss<0)call MPI_Abort(MPI_COMM_WORLD,6,err)
  write(*,'(a,1x,i0,2(1x,i0),6(1x,es26.16e3))')'MEASURE',rank,baseline,rss,times,global_energy,global_error
  if(mode==2)call lcfo_dense_free(dense)
  call MPI_Finalize(err)
contains
  complex(8) function packet(point,index) result(value)
    integer,intent(in) :: point(3),index
    integer :: center(3),axis,q
    complex(8) :: factor
    center=[mod(index,kd(1)),mod(index/kd(1),kd(2)),index/(kd(1)*kd(2))]
    value=1d0/sqrt(volume)
    do axis=1,3
      factor=0d0
      do q=0,kd(axis)-1
        factor=factor+exp(cmplx(0d0,2*pi*q*(real(point(axis),8)/ngrid-real(center(axis),8)/kd(axis)),8))
      enddo
      value=value*factor/sqrt(real(kd(axis),8))
    enddo
  end function
  subroutine require_success(value)
    integer,intent(in) :: value
    integer :: total
    call MPI_Allreduce(value,total,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)
    if(total/=0)call MPI_Abort(MPI_COMM_WORLD,10,err)
  end subroutine
end program
