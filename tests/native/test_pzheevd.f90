program test_pzheevd
  use mpi
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer :: rank,procs,ierr,ctxt,nprow,npcol,myrow,mycol,n,nr,nc,nb=64,info,j,i,gi,gj,caseid
  integer :: desc(9),numroc,indxl2g
  real(8) :: err,orth,expected
  real(8),allocatable :: eval(:)
  complex(8),allocatable :: mat(:,:),vec(:,:),globalz(:,:),fulla(:,:),part(:,:)
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,procs,ierr)
  nprow=1;npcol=procs
  if(procs==4) then
    nprow=2;npcol=2
  end if
  call blacs_get(-1,0,ctxt)
  call blacs_gridinit(ctxt,'R',nprow,npcol)
  call blacs_gridinfo(ctxt,nprow,npcol,myrow,mycol)
  do caseid=1,4
    select case(caseid)
    case(1);n=3
    case(2);n=65
    case(3);n=129
    case(4);n=7
    end select
    nr=numroc(n,nb,myrow,0,nprow);nc=numroc(n,nb,mycol,0,npcol)
    allocate(mat(max(1,nr),max(1,nc)),vec(max(1,nr),max(1,nc)),eval(n), &
             fulla(n,n),globalz(n,n),part(n,n))
    fulla=(0d0,0d0)
    do i=1,n
      fulla(i,i)=4d0
      if(i<n.and.caseid/=4) then
        fulla(i,i+1)=(1d0,1d0);fulla(i+1,i)=(1d0,-1d0)
      end if
    end do
    mat=0;vec=0
    do j=1,nc
      gj=indxl2g(j,nb,mycol,0,npcol)
      do i=1,nr
        gi=indxl2g(i,nb,myrow,0,nprow)
        mat(i,j)=fulla(gi,gj)
      end do
    end do
    call descinit(desc,n,n,nb,nb,0,0,ctxt,max(1,nr),info)
    if(info/=0) error stop 'descriptor'
    call solve(n,mat,eval,vec,desc,info)
    if(info/=0) error stop 'solve'
    part=0
    do j=1,nc
      gj=indxl2g(j,nb,mycol,0,npcol)
      do i=1,nr
        gi=indxl2g(i,nb,myrow,0,nprow);part(gi,gj)=vec(i,j)
      end do
    end do
    call MPI_Allreduce(part,globalz,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
    part=matmul(fulla,globalz)
    do j=1,n
      part(:,j)=part(:,j)-eval(j)*globalz(:,j)
      expected=4d0
      if(caseid/=4) expected=4d0-2d0*sqrt(2d0)*cos(j*acos(-1d0)/(n+1))
      if(abs(eval(j)-expected)>1d-11) then
        print *, 'eigenvalue mismatch',rank,n,j,eval(j),expected
        error stop 'eigenvalues'
      end if
    end do
    err=maxval(abs(part))
    part=matmul(conjg(transpose(globalz)),globalz)
    do i=1,n
      part(i,i)=part(i,i)-1d0
    end do
    orth=maxval(abs(part))
    if(err>1d-11.or.orth>1d-11) error stop 'residual/orthogonality'
    if(rank==0) print *, 'PASS n, ranks, residual, orthogonality:',n,procs,err,orth
    deallocate(mat,vec,eval,fulla,globalz,part)
  end do
  call blacs_gridexit(ctxt)
  call MPI_Finalize(ierr)
contains
  subroutine solve(nactive,projected,eigenvalue,z,desc,istat)
    integer,intent(in) :: nactive,desc(9)
    complex(8),intent(inout) :: projected(:,:)
    complex(8),intent(out) :: z(:,:)
    real(8),intent(out) :: eigenvalue(:)
    integer,intent(out) :: istat
    integer :: n,nlocked,i,lwork0,lrwork0,liwork0,info0,alloc_status,iwork_query(1)
    complex(8) :: work_query(1)
    real(8) :: rwork_query(1)
    complex(8),allocatable :: work0(:)
    real(8),allocatable :: rwork0(:)
    integer,allocatable :: iwork0(:)
    n=nactive;nlocked=0
      call pzheevd('V','L',nactive,projected,1,1,desc, &
        eigenvalue(nlocked+1:n),z,1,1,desc,work_query,-1,rwork_query,-1,iwork_query,-1,info0)
      istat=merge(1,0,info0/=0)
      call sync_status(istat)
      if(istat/=0) return
      ! Check scalar query values before converting to the default integer ABI.
      if(.not.ieee_is_finite(real(work_query(1),8))) istat=1
      if(.not.ieee_is_finite(rwork_query(1))) istat=1
      call sync_status(istat)
      if(istat/=0) return
      if(real(work_query(1),8)<1d0 .or. rwork_query(1)<1d0 .or. iwork_query(1)<1) istat=1
      if(real(work_query(1),8)>real(huge(lwork0)-1,8)) istat=1
      if(rwork_query(1)>real(huge(lrwork0)-1,8)) istat=1
      call sync_status(istat)
      if(istat/=0) return
      lwork0=ceiling(real(work_query(1),8))
      lrwork0=ceiling(rwork_query(1))
      liwork0=iwork_query(1)
      call comm_get_max(lwork0,MPI_COMM_WORLD)
      call comm_get_max(lrwork0,MPI_COMM_WORLD)
      call comm_get_max(liwork0,MPI_COMM_WORLD)
      allocate(work0(lwork0),rwork0(lrwork0),iwork0(liwork0),stat=alloc_status)
      istat=merge(1,0,alloc_status/=0)
      call sync_status(istat)
      if(istat/=0) then
        if(allocated(work0)) deallocate(work0)
        if(allocated(rwork0)) deallocate(rwork0)
        if(allocated(iwork0)) deallocate(iwork0)
        return
      end if
      call pzheevd('V','L',nactive,projected,1,1,desc, &
        eigenvalue(nlocked+1:n),z,1,1,desc,work0,lwork0,rwork0,lrwork0,iwork0,liwork0,info0)
      istat=merge(1,0,info0/=0)
      do i=nlocked+1,n
        if(.not.ieee_is_finite(eigenvalue(i))) istat=1
      end do
      call sync_status(istat)
      deallocate(work0,rwork0,iwork0)
      if(istat/=0) return
    call MPI_Bcast(eigenvalue,nactive,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info0)
    if(rank==0) print *, 'workspace counts:',lwork0,lrwork0,liwork0
  end subroutine
  subroutine sync_status(value)
    integer,intent(inout) :: value
    call comm_get_max(value,MPI_COMM_WORLD)
  end subroutine
  subroutine comm_get_max(value,comm)
    integer,intent(inout) :: value
    integer,intent(in) :: comm
    integer :: tmp, ierr0
    call MPI_Allreduce(value,tmp,1,MPI_INTEGER,MPI_MAX,comm,ierr0)
    value=tmp
  end subroutine
end program
