program probe
  use mpi
  use iso_fortran_env,only:int64
  use ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
  use exx_gamma_pairs,only:gamma_pair_action
  implicit none
  integer::rank,np,err,orb,spatial,nsp,comm_r,comm_o,cx,ax,sy,sz,owner,first,last,nloc,g,i,j,k,hx,hy,hz,hg
  integer::n(3),m(3),lo(3),dims(3),coords(3),groups(3),a(3),b(3),q(3),stat,no,caseid
  integer(int64)::pairs,scratch,totalpairs
  complex(8),allocatable::psi(:,:,:),action(:,:,:),ref(:,:,:),kernel(:)
  real(8),allocatable::weights(:),multiplier(:)
  real(8)::pi,phase,kv,error,global_error,wt,seconds,t0
  complex(8)::v,density
  character(32)::arg
  call MPI_Init(err);call MPI_Comm_rank(MPI_COMM_WORLD,rank,err);call MPI_Comm_size(MPI_COMM_WORLD,np,err)
  call get_command_argument(1,arg);nsp=1
  if(trim(arg)=='spatial'.or.trim(arg)=='weight_mismatch')nsp=2
  spatial=mod(rank,nsp);orb=rank/nsp
  call MPI_Comm_split(MPI_COMM_WORLD,orb,spatial,comm_r,err)
  call MPI_Comm_split(MPI_COMM_WORLD,spatial,orb,comm_o,err)
  call MPI_Comm_split(comm_r,0,spatial,cx,err)
  n=[4,4,4];dims=[nsp,1,1];coords=[spatial,0,0];groups=[cx,MPI_COMM_SELF,MPI_COMM_SELF]
  m=n/dims;lo=coords*m;no=7;first=orb*no/(np/nsp)+1;last=(orb+1)*no/(np/nsp);nloc=last-first+1
  allocate(psi(product(m),nloc,1),action(product(m),nloc,1),ref(product(m),nloc,1))
  allocate(weights(nloc),multiplier(product(m)),kernel(product(n)))
  pi=acos(-1d0)
  do caseid=1,3
    do j=1,nloc
      weights(j)=occupation(first+j-1,caseid)
      g=0
      do sz=0,m(3)-1;do sy=0,m(2)-1;do ax=0,m(1)-1
        g=g+1;a=[ax,sy,sz]+lo;psi(g,j,1)=field(a,first+j-1)
      enddo;enddo;enddo
    enddo
    g=0
    do sz=0,m(3)-1;do sy=0,m(2)-1;do ax=0,m(1)-1
      g=g+1;q=[ax,sy,sz]+lo;where(q>=(n+1)/2)q=q-n
      multiplier(g)=1d0/(1d0+sum(q*q))
    enddo;enddo;enddo
    if(trim(arg)=='source_nan'.and.rank==np-1.and.caseid==1.and.nloc>0) &
      psi(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
    if(trim(arg)=='source_inf'.and.rank==np-1.and.caseid==1.and.nloc>0) &
      psi(1,1,1)=cmplx(0d0,ieee_value(0d0,ieee_positive_inf),8)
    if(trim(arg)=='nan'.and.rank==np-1.and.caseid==1.and.nloc>0)weights(1)=ieee_value(0d0,ieee_quiet_nan)
    if(trim(arg)=='inf'.and.rank==np-1.and.caseid==1.and.nloc>0) &
      weights(1)=ieee_value(0d0,ieee_positive_inf)
    if(trim(arg)=='weight_mismatch'.and.rank==np-1.and.caseid==1.and.nloc>0)weights(1)=.5d0
    if(trim(arg)=='negative'.and.rank==np-1.and.caseid==1.and.nloc>0)weights(1)=-1d0
    t0=MPI_Wtime()
    call gamma_pair_action(n,dims,coords,groups,comm_r,comm_o,psi,weights,multiplier,action,stat,pairs,scratch)
    seconds=MPI_Wtime()-t0
    if(trim(arg)=='nan'.or.trim(arg)=='negative'.or.trim(arg)=='inf'.or.trim(arg)=='weight_mismatch'.or. &
       trim(arg)=='source_nan'.or.trim(arg)=='source_inf')then
      if(stat==0)call MPI_Abort(MPI_COMM_WORLD,1,err)
      if(rank==0)print *,'GAMMA_PAIR_REJECT_PASS'
      call MPI_Finalize(err);stop
    endif
    if(stat/=0)call MPI_Abort(MPI_COMM_WORLD,2,err)
    ! Independent periodic direct convolution, including complex orbitals and fractional weights.
    do k=1,product(n)
      b=[mod(k-1,n(1)),mod((k-1)/n(1),n(2)),(k-1)/(n(1)*n(2))];v=0d0
      do hz=0,n(3)-1;do hy=0,n(2)-1;do hx=0,n(1)-1
        q=[hx,hy,hz];where(q>=(n+1)/2)q=q-n
        phase=2*pi*sum(real(q*b,8)/n);v=v+cmplx(cos(phase),sin(phase),8)/(1d0+sum(q*q))
      enddo;enddo;enddo
      kernel(k)=v/product(n)
    enddo
    ref=0d0;g=0
    do sz=0,m(3)-1;do sy=0,m(2)-1;do ax=0,m(1)-1
      g=g+1;a=[ax,sy,sz]+lo
      do j=1,nloc
        do i=1,no
          v=0d0
          do hz=0,n(3)-1;do hy=0,n(2)-1;do hx=0,n(1)-1
            b=[hx,hy,hz];q=modulo(a-b,n);hg=1+q(1)+n(1)*(q(2)+n(2)*q(3))
            v=v+kernel(hg)*conjg(field(b,i))*field(b,first+j-1)
          enddo;enddo;enddo
          ref(g,j,1)=ref(g,j,1)-occupation(i,caseid)*field(a,i)*v
        enddo
      enddo
    enddo;enddo;enddo
    error=0d0
    if(nloc>0)error=maxval(abs(action-ref))
    call MPI_Allreduce(error,global_error,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,err)
    call MPI_Allreduce(pairs,totalpairs,1,MPI_INTEGER8,MPI_SUM,comm_o,err)
    if(global_error>1d-11.or.totalpairs/=int(no*(no+1)/2,int64))call MPI_Abort(MPI_COMM_WORLD,3,err)
    if(rank==0)print *,'GAMMA_PAIR_PASS',caseid,global_error,totalpairs,scratch,seconds
  enddo
  call MPI_Finalize(err)
contains
  real(8) function occupation(i,caseid) result(w)
    integer,intent(in)::i,caseid
    w=1d0
    if(caseid==2)w=real(mod(i,4),8)/4d0
    if(caseid==3)w=0d0
  end function
  complex(8) function field(a,i) result(z)
    integer,intent(in)::a(3),i
    real(8)::x
    x=real(1+a(1)+3*a(2)+7*a(3),8)
    z=cmplx(sin(.13d0*i*x),cos(.07d0*(i+1)*x),8)
  end function
end program
