module wannier_patch_owner
  use iso_c_binding
  use iso_fortran_env, only: int64
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use wannier_pair_owner, only: orbital_owners
  implicit none
  private
  public :: orbital_patch,patch_metrics,patch_apply,pack_exact_patch
  include 'fftw3.f03'
  type orbital_patch
    integer,allocatable :: row(:)
    complex(8),allocatable :: value(:)
  end type
  type patch_metrics
    integer(int64) :: pairs=0,evaluated=0,zero_overlap=0,fft_points=0,max_fft_points=0
    integer(int64) :: source_payload=0,reverse_payload=0
  end type
contains
  subroutine pack_exact_patch(dense,patch)
    complex(8),intent(in) :: dense(:)
    type(orbital_patch),intent(out) :: patch
    integer :: i,k,n
    ! Remove exact zeros only. No support threshold or approximation is introduced.
    n=0
    do i=1,size(dense)
      if(keep_value(dense(i)))n=n+1
    enddo
    allocate(patch%row(n),patch%value(n));k=0
    do i=1,size(dense)
      if(.not.keep_value(dense(i)))cycle
      k=k+1;patch%row(k)=i;patch%value(k)=dense(i)
    enddo
  end subroutine

  subroutine patch_apply(shape,norb,ids,patches,kernel,comm,batch,action,status,metrics)
    integer,intent(in) :: shape(3),norb,ids(:),comm,batch
    type(orbital_patch),intent(in) :: patches(:)
    real(8),intent(in) :: kernel(:)
    complex(8),intent(out) :: action(:,:)
    integer,intent(out) :: status
    type(patch_metrics),intent(out) :: metrics
    type(orbital_patch),allocatable :: sources(:)
    type(patch_metrics) :: local
    integer :: rank,peers,e,bad,par(5),pmin(5),pmax(5),ng,i,j,k,t,owner,nsource,first,nb,ns,lo,hi
    integer :: rh(3),x,y,z,q,neg,nlocal
    integer,allocatable :: owners(:),source_ids(:),lengths(:),wire_rows(:)
    complex(8),allocatable :: wire_values(:),reverse(:,:),onto_i(:),onto_j(:)
    real(8),allocatable :: reference_kernel(:)
    integer(int64) :: volume,fp,all_counts(6),local_counts(6),sent_bytes
    action=0;metrics=patch_metrics();local=patch_metrics()
    call MPI_Comm_rank(comm,rank,e);call MPI_Comm_size(comm,peers,e)
    par=[shape,norb,batch]
    call MPI_Allreduce(par,pmin,5,MPI_INTEGER,MPI_MIN,comm,e)
    call MPI_Allreduce(par,pmax,5,MPI_INTEGER,MPI_MAX,comm,e)
    bad=0;volume=1_int64
    if(any(pmin/=pmax).or.any(shape<1).or.norb<1.or.batch<1)bad=1
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
    ng=int(volume);nlocal=size(ids)
    allocate(owners(norb));call orbital_owners(norb,peers,owners)
    if(size(patches)/=nlocal.or.count(owners==rank)/=nlocal)bad=1
    if(size(action,1)/=ng.or.size(action,2)/=nlocal.or.size(kernel)/=ng)bad=1
    t=0
    do i=1,norb
      if(owners(i)/=rank)cycle
      t=t+1
      if(t<=nlocal)then
        if(ids(t)/=i)bad=1
      endif
    enddo
    do j=1,size(patches)
      if(.not.allocated(patches(j)%row).or..not.allocated(patches(j)%value))then
        bad=1;cycle
      endif
      if(size(patches(j)%row)/=size(patches(j)%value))then
        bad=1;cycle
      endif
      do i=1,size(patches(j)%row)
        if(patches(j)%row(i)<1.or.patches(j)%row(i)>ng)bad=1
        if(i>1)then
          if(patches(j)%row(i)<=patches(j)%row(i-1))bad=1
        endif
        if(.not.ieee_is_finite(real(patches(j)%value(i),8)))bad=1
        if(.not.ieee_is_finite(aimag(patches(j)%value(i))))bad=1
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
    rh=0
    do z=0,shape(3)-1
      do y=0,shape(2)-1
        do x=0,shape(1)-1
          q=1+x+shape(1)*(y+shape(2)*z)
          neg=1+modulo(-x,shape(1))+shape(1)*(modulo(-y,shape(2))+shape(2)*modulo(-z,shape(3)))
          if(abs(kernel(q)-kernel(neg))>0d0)bad=1
          if(abs(kernel(q))>0d0)rh=max(rh,[min(x,shape(1)-x),min(y,shape(2)-y),min(z,shape(3)-z)])
        enddo
      enddo
    enddo
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    allocate(source_ids(norb),lengths(batch),sources(batch),reverse(ng,batch))
    do owner=0,peers-1
      nsource=0
      do i=1,norb
        if(owners(i)/=owner)cycle
        nsource=nsource+1;source_ids(nsource)=i
      enddo
      do first=1,nsource,batch
        nb=min(batch,nsource-first+1);lengths=0
        if(rank==owner)then
          do j=1,nb
            lengths(j)=size(patches(first+j-1)%row)
          enddo
        endif
        call MPI_Bcast(lengths,nb,MPI_INTEGER,owner,comm,e)
        ns=sum(lengths(:nb));allocate(wire_rows(max(1,ns)),wire_values(max(1,ns)))
        lo=1
        if(rank==owner)then
          do j=1,nb
            hi=lo+lengths(j)-1
            wire_rows(lo:hi)=patches(first+j-1)%row
            wire_values(lo:hi)=patches(first+j-1)%value
            lo=hi+1
          enddo
        endif
        call MPI_Bcast(wire_rows,ns,MPI_INTEGER,owner,comm,e)
        call MPI_Bcast(wire_values,ns,MPI_DOUBLE_COMPLEX,owner,comm,e)
        if(rank==owner)then
          local%source_payload=local%source_payload+(20_int64*ns+4_int64*nb)*int(peers-1,int64)
        endif
        lo=1
        do j=1,nb
          hi=lo+lengths(j)-1
          allocate(sources(j)%row(lengths(j)),sources(j)%value(lengths(j)))
          sources(j)%row=wire_rows(lo:hi);sources(j)%value=wire_values(lo:hi);lo=hi+1
        enddo
        reverse=0
        do j=1,nb
          do i=1,nlocal
            if(source_ids(first+j-1)>ids(i))cycle
            call patch_pair_fft(shape,kernel,rh,sources(j),patches(i),onto_i,onto_j,status,fp)
            if(status/=0)then
              bad=1;cycle
            endif
            local%pairs=local%pairs+1
            if(fp==0)then
              local%zero_overlap=local%zero_overlap+1
            else
              local%evaluated=local%evaluated+1
              local%fft_points=local%fft_points+fp
              local%max_fft_points=max(local%max_fft_points,fp)
            endif
            do k=1,size(sources(j)%row)
              q=sources(j)%row(k);action(q,i)=action(q,i)+onto_i(k)
            enddo
            if(source_ids(first+j-1)/=ids(i))then
              do k=1,size(patches(i)%row)
                q=patches(i)%row(k);reverse(q,j)=reverse(q,j)+onto_j(k)
              enddo
            endif
            deallocate(onto_i,onto_j)
          enddo
        enddo
        call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
        if(status/=0)return
        call return_sparse(comm,owner,reverse(:,:nb),action,first,sent_bytes)
        local%reverse_payload=local%reverse_payload+sent_bytes
        do j=1,nb
          deallocate(sources(j)%row,sources(j)%value)
        enddo
        deallocate(wire_rows,wire_values)
      enddo
    enddo
    local_counts=[local%pairs,local%evaluated,local%zero_overlap,local%fft_points, &
      local%source_payload,local%reverse_payload]
    call MPI_Allreduce(local_counts,all_counts,6,MPI_INTEGER8,MPI_SUM,comm,e)
    metrics%pairs=all_counts(1);metrics%evaluated=all_counts(2);metrics%zero_overlap=all_counts(3)
    metrics%fft_points=all_counts(4);metrics%source_payload=all_counts(5);metrics%reverse_payload=all_counts(6)
    call MPI_Allreduce(local%max_fft_points,metrics%max_fft_points,1,MPI_INTEGER8,MPI_MAX,comm,e)
    if(metrics%pairs/=int(norb,int64)*int(norb+1,int64)/2_int64)bad=1
    do j=1,nlocal
      do i=1,ng
        if(.not.ieee_is_finite(real(action(i,j),8)))bad=1
        if(.not.ieee_is_finite(aimag(action(i,j))))bad=1
      enddo
    enddo
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
  end subroutine

  subroutine patch_pair_fft(shape,kernel,rh,source,target,onto_i,onto_j,status,fft_points)
    integer,intent(in) :: shape(3),rh(3)
    real(8),intent(in) :: kernel(:)
    type(orbital_patch),intent(in) :: source,target
    complex(8),allocatable,intent(out) :: onto_i(:),onto_j(:)
    integer,intent(out) :: status
    integer(int64),intent(out) :: fft_points
    integer :: i,j,k,n,q,axis,t,prev,first,last,gap,maxgap,start,origin(3),length(3),xyz(3),d(3)
    integer :: x,y,z,idx,extent,nfft
    integer,allocatable :: overlap(:),si(:),ti(:)
    logical,allocatable :: occupied(:)
    complex(8),allocatable :: a(:),spectrum(:),kspectrum(:)
    type(c_ptr) :: forward,backward
    allocate(onto_i(size(source%row)),onto_j(size(target%row)))
    onto_i=0;onto_j=0;status=0;fft_points=0
    allocate(overlap(min(size(source%row),size(target%row))))
    allocate(si(size(overlap)),ti(size(overlap)))
    i=1;j=1;n=0
    do while(i<=size(source%row).and.j<=size(target%row))
      if(source%row(i)<target%row(j))then
        i=i+1
      else if(source%row(i)>target%row(j))then
        j=j+1
      else
        n=n+1;overlap(n)=source%row(i);si(n)=i;ti(n)=j;i=i+1;j=j+1
      endif
    enddo
    if(n==0)return
    ! Minimum cyclic interval for the pair density on each periodic axis.
    do axis=1,3
      allocate(occupied(shape(axis)));occupied=.false.
      do k=1,n
        call coordinates(overlap(k),shape,xyz);occupied(xyz(axis)+1)=.true.
      enddo
      first=-1;prev=-1;last=-1;maxgap=-1;start=0
      do t=0,shape(axis)-1
        if(.not.occupied(t+1))cycle
        if(first<0)first=t
        if(prev>=0)then
          gap=t-prev
          if(gap>maxgap)then
            maxgap=gap;start=t
          endif
        endif
        prev=t;last=t
      enddo
      gap=first+shape(axis)-last
      if(gap>maxgap)then
        maxgap=gap;start=first
      endif
      extent=shape(axis)-maxgap+1
      length(axis)=int(min(int(shape(axis),int64),int(extent,int64)+2_int64*int(rh(axis),int64)))
      origin(axis)=modulo(start-rh(axis),shape(axis))
      deallocate(occupied)
    enddo
    nfft=product(length);fft_points=int(nfft,int64)
    allocate(a(nfft),spectrum(nfft),kspectrum(nfft))
    forward=fftw_plan_dft_3d(length(3),length(2),length(1),a,spectrum,FFTW_FORWARD,FFTW_ESTIMATE)
    backward=fftw_plan_dft_3d(length(3),length(2),length(1),spectrum,a,FFTW_BACKWARD,FFTW_ESTIMATE)
    if(.not.c_associated(forward).or..not.c_associated(backward))then
      status=1
      if(c_associated(forward))call fftw_destroy_plan(forward)
      if(c_associated(backward))call fftw_destroy_plan(backward)
      return
    endif
    a=0
    do z=0,length(3)-1
      do y=0,length(2)-1
        do x=0,length(1)-1
          d=[x,y,z]
          do axis=1,3
            if(d(axis)>length(axis)/2)d(axis)=d(axis)-length(axis)
          enddo
          if(any(abs(d)>rh))cycle
          xyz=modulo(d,shape)
          q=1+xyz(1)+shape(1)*(xyz(2)+shape(2)*xyz(3))
          idx=1+x+length(1)*(y+length(2)*z);a(idx)=cmplx(kernel(q),0d0,8)
        enddo
      enddo
    enddo
    call fftw_execute_dft(forward,a,spectrum);kspectrum=spectrum;a=0
    do k=1,n
      call coordinates(overlap(k),shape,xyz);xyz=modulo(xyz-origin,shape)
      idx=1+xyz(1)+length(1)*(xyz(2)+length(2)*xyz(3))
      a(idx)=conjg(source%value(si(k)))*target%value(ti(k))
    enddo
    call fftw_execute_dft(forward,a,spectrum);spectrum=spectrum*kspectrum
    call fftw_execute_dft(backward,spectrum,a);a=a/real(nfft,8)
    do k=1,size(source%row)
      call coordinates(source%row(k),shape,xyz);xyz=modulo(xyz-origin,shape)
      if(any(xyz>=length))cycle
      idx=1+xyz(1)+length(1)*(xyz(2)+length(2)*xyz(3))
      onto_i(k)=-source%value(k)*a(idx)
    enddo
    do k=1,size(target%row)
      call coordinates(target%row(k),shape,xyz);xyz=modulo(xyz-origin,shape)
      if(any(xyz>=length))cycle
      idx=1+xyz(1)+length(1)*(xyz(2)+length(2)*xyz(3))
      onto_j(k)=-target%value(k)*conjg(a(idx))
    enddo
    call fftw_destroy_plan(forward);call fftw_destroy_plan(backward)
  end subroutine

  logical function keep_value(value)
    complex(8),intent(in) :: value
    keep_value=.true.
    if(.not.ieee_is_finite(real(value,8)))return
    if(.not.ieee_is_finite(aimag(value)))return
    keep_value=abs(value)>0d0
  end function

  subroutine return_sparse(comm,owner,reverse,action,first,sent_bytes)
    integer,intent(in) :: comm,owner,first
    complex(8),intent(in) :: reverse(:,:)
    complex(8),intent(inout) :: action(:,:)
    integer(int64),intent(out) :: sent_bytes
    integer :: rank,peers,e,ng,nb,nflat,nz,i,j,k,sender,header(2),count,wire_mode,idx
    integer,allocatable :: keys(:)
    complex(8),allocatable :: values(:),summed(:)
    call MPI_Comm_rank(comm,rank,e);call MPI_Comm_size(comm,peers,e)
    ng=size(reverse,1);nb=size(reverse,2);nflat=ng*nb;nz=0;sent_bytes=0
    do j=1,nb
      do i=1,ng
        if(keep_value(reverse(i,j)))nz=nz+1
      enddo
    enddo
    wire_mode=0
    if(20_int64*nz>16_int64*nflat)wire_mode=1
    if(rank==owner)then
      allocate(summed(nflat));summed=0
    endif
    ! Stream one sender at a time: root never holds all peers' contributions.
    do sender=0,peers-1
      if(rank==sender)then
        if(wire_mode==0)then
          count=nz;allocate(keys(max(1,count)),values(max(1,count)));k=0
          do j=1,nb
            do i=1,ng
              if(.not.keep_value(reverse(i,j)))cycle
              k=k+1;keys(k)=i+ng*(j-1);values(k)=reverse(i,j)
            enddo
          enddo
        else
          count=nflat;allocate(values(count));values=reshape(reverse,[nflat])
        endif
        header=[count,wire_mode]
        if(rank/=owner)then
          call MPI_Send(header,2,MPI_INTEGER,owner,2100,comm,e)
          if(wire_mode==0.and.count>0)call MPI_Send(keys,count,MPI_INTEGER,owner,2101,comm,e)
          if(count>0)call MPI_Send(values,count,MPI_DOUBLE_COMPLEX,owner,2102,comm,e)
          sent_bytes=8_int64+16_int64*count
          if(wire_mode==0)sent_bytes=sent_bytes+4_int64*count
        endif
      endif
      if(rank==owner)then
        if(sender/=owner)then
          call MPI_Recv(header,2,MPI_INTEGER,sender,2100,comm,MPI_STATUS_IGNORE,e)
          count=header(1)
          allocate(values(max(1,count)))
          if(header(2)==0)then
            allocate(keys(max(1,count)))
            if(count>0)call MPI_Recv(keys,count,MPI_INTEGER,sender,2101,comm,MPI_STATUS_IGNORE,e)
          endif
          if(count>0)call MPI_Recv(values,count,MPI_DOUBLE_COMPLEX,sender,2102,comm,MPI_STATUS_IGNORE,e)
        endif
        if(header(2)==0)then
          do k=1,header(1)
            idx=keys(k);summed(idx)=summed(idx)+values(k)
          enddo
        else
          summed=summed+values(:nflat)
        endif
      endif
      if(rank==sender.or.rank==owner)then
        if(allocated(keys))deallocate(keys)
        if(allocated(values))deallocate(values)
      endif
    enddo
    if(rank==owner)action(:,first:first+nb-1)=action(:,first:first+nb-1)+reshape(summed,[ng,nb])
  end subroutine

  subroutine coordinates(row,shape,xyz)
    integer,intent(in) :: row,shape(3)
    integer,intent(out) :: xyz(3)
    xyz(1)=modulo(row-1,shape(1));xyz(2)=modulo((row-1)/shape(1),shape(2))
    xyz(3)=(row-1)/(shape(1)*shape(2))
  end subroutine
end module
