module wannier_pair_owner
  use iso_c_binding
  use iso_fortran_env, only: int64
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: orbital_owners, pair_owner_apply
  include 'fftw3.f03'
contains
  subroutine orbital_owners(norb,peers,owners)
    integer,intent(in) :: norb,peers
    integer,intent(out) :: owners(norb)
    integer(int64) :: load(peers)
    integer :: i,r
    load=0
    do i=norb,1,-1
      r=minloc(load,dim=1)
      owners(i)=r-1
      load(r)=load(r)+int(i,int64)
    enddo
  end subroutine

  subroutine pair_owner_apply(shape,norb,ids,psi,kernel,comm,batch,action,status,pairs)
    integer,intent(in) :: shape(3),norb,ids(:),comm,batch
    complex(8),intent(in) :: psi(:,:)
    real(8),intent(in) :: kernel(:)
    complex(8),intent(out) :: action(:,:)
    integer,intent(out) :: status
    integer(int64),intent(out) :: pairs
    integer :: rank,peers,e,bad,params(5),pmin(5),pmax(5),ng,i,j,nlocal,nsource,first,nb,t,ii,jj
    integer :: x,y,z,q,neg,owner
    integer,allocatable :: owners(:),source_ids(:)
    complex(8),allocatable :: source(:,:),reverse(:,:),returned(:,:),a(:),spectrum(:),kspectrum(:)
    real(8),allocatable :: reference_kernel(:)
    integer(int64) :: local_pairs,volume
    type(c_ptr) :: forward,backward
    action=0; pairs=0
    call MPI_Comm_rank(comm,rank,e)
    call MPI_Comm_size(comm,peers,e)
    params=[shape,norb,batch]
    call MPI_Allreduce(params,pmin,5,MPI_INTEGER,MPI_MIN,comm,e)
    call MPI_Allreduce(params,pmax,5,MPI_INTEGER,MPI_MAX,comm,e)
    bad=0
    if(any(pmin/=pmax).or.any(shape<1).or.norb<1.or.batch<1)bad=1
    volume=1_int64
    do i=1,3
      if(shape(i)<1)cycle
      if(volume>int(huge(ng),int64)/int(shape(i),int64))then
        bad=1
      else
        volume=volume*int(shape(i),int64)
      endif
    enddo
    if(volume*int(max(1,batch),int64)>int(huge(ng),int64))bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    ng=product(shape);nlocal=size(ids)
    allocate(owners(norb))
    call orbital_owners(norb,peers,owners)
    if(size(psi,1)/=ng.or.size(psi,2)/=nlocal)bad=1
    if(size(action,1)/=ng.or.size(action,2)/=nlocal.or.size(kernel)/=ng)bad=1
    if(count(owners==rank)/=nlocal)bad=1
    t=0
    do i=1,norb
      if(owners(i)/=rank)cycle
      t=t+1
      if(t<=nlocal)then
        if(ids(t)/=i)bad=1
      endif
    enddo
    do j=1,size(psi,2)
      do i=1,size(psi,1)
        if(.not.ieee_is_finite(real(psi(i,j),8)))bad=1
        if(.not.ieee_is_finite(aimag(psi(i,j))))bad=1
      enddo
    enddo
    do i=1,size(kernel)
      if(.not.ieee_is_finite(kernel(i)))bad=1
    enddo
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    allocate(reference_kernel(ng));reference_kernel=kernel
    call MPI_Bcast(reference_kernel,ng,MPI_DOUBLE_PRECISION,0,comm,e)
    if(any(abs(reference_kernel-kernel)>0d0))bad=1
    do z=0,shape(3)-1
      do y=0,shape(2)-1
        do x=0,shape(1)-1
          q=1+x+shape(1)*(y+shape(2)*z)
          neg=1+modulo(-x,shape(1))+shape(1)*(modulo(-y,shape(2))+shape(2)*modulo(-z,shape(3)))
          if(abs(kernel(q)-kernel(neg))>0d0)bad=1
        enddo
      enddo
    enddo
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    allocate(a(ng),spectrum(ng),kspectrum(ng))
    forward=fftw_plan_dft_3d(shape(3),shape(2),shape(1),a,spectrum,FFTW_FORWARD,FFTW_ESTIMATE)
    backward=fftw_plan_dft_3d(shape(3),shape(2),shape(1),spectrum,a,FFTW_BACKWARD,FFTW_ESTIMATE)
    bad=0
    if(.not.c_associated(forward).or..not.c_associated(backward))bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)then
      if(c_associated(forward))call fftw_destroy_plan(forward)
      if(c_associated(backward))call fftw_destroy_plan(backward)
      return
    endif
    a=cmplx(kernel,0d0,8)
    call fftw_execute_dft(forward,a,spectrum)
    kspectrum=spectrum
    allocate(source(ng,batch),reverse(ng,batch),returned(ng,batch),source_ids(norb))
    local_pairs=0
    do owner=0,peers-1
      nsource=0
      do i=1,norb
        if(owners(i)/=owner)cycle
        nsource=nsource+1;source_ids(nsource)=i
      enddo
      do first=1,nsource,batch
        nb=min(batch,nsource-first+1);source=0;reverse=0
        if(rank==owner)source(:,:nb)=psi(:,first:first+nb-1)
        call MPI_Bcast(source,ng*nb,MPI_DOUBLE_COMPLEX,owner,comm,e)
        do j=1,nb
          jj=source_ids(first+j-1)
          do i=1,nlocal
            ii=ids(i)
            if(jj>ii)cycle
            a=conjg(source(:,j))*psi(:,i)
            call fftw_execute_dft(forward,a,spectrum)
            spectrum=spectrum*kspectrum
            call fftw_execute_dft(backward,spectrum,a)
            a=a/real(ng,8)
            action(:,i)=action(:,i)-source(:,j)*a
            if(ii/=jj)reverse(:,j)=reverse(:,j)-psi(:,i)*conjg(a)
            local_pairs=local_pairs+1
          enddo
        enddo
        call MPI_Reduce(reverse,returned,ng*nb,MPI_DOUBLE_COMPLEX,MPI_SUM,owner,comm,e)
        if(rank==owner)action(:,first:first+nb-1)=action(:,first:first+nb-1)+returned(:,:nb)
      enddo
    enddo
    call fftw_destroy_plan(forward)
    call fftw_destroy_plan(backward)
    call MPI_Allreduce(local_pairs,pairs,1,MPI_INTEGER8,MPI_SUM,comm,e)
    bad=0
    do j=1,nlocal
      do i=1,ng
        if(.not.ieee_is_finite(real(action(i,j),8)))bad=1
        if(.not.ieee_is_finite(aimag(action(i,j))))bad=1
      enddo
    enddo
    if(pairs/=int(norb,int64)*int(norb+1,int64)/2_int64)bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
  end subroutine
end module
