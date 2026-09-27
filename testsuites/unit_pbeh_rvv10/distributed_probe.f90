program probe
  use fftw_pencils, only: pencil_transform,fftw_pencil_plans_created
  use mpi
  use rvv10, only: rvv10_evaluate,rvv10_periodic
  use rvv10_distributed, only: rvv10_evaluate_distributed
  implicit none
  integer :: ierr,rank,nproc,axis,comm(3),coords(3),dims(3),n(3),lo(3),local_n(3),i,j,k,g,l,status,color,key,nq,channel,plans_before,batch_count
  character(16) :: argument
  logical :: used
  real(8) :: err,total_err,h(3),coef(4,3),matrix(3,3)
  real(8),allocatable :: density(:,:,:),grad(:,:,:,:),ep(:,:,:),vp(:,:,:),flux(:,:,:),potential(:,:,:)
  real(8),allocatable :: local_parts(:,:),parts(:,:)
  complex(8),allocatable :: batch_input(:,:),batch_output(:,:),batch_back(:,:),fft_a(:),fft_b(:)
  complex(8),allocatable :: old_input(:),old_a(:),old_b(:),old_reference(:)
  real(8),allocatable :: rho(:),sigma(:),e(:),v(:),w(:),lr(:),ls(:),le(:),lv(:),lw(:)
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  nq=16
  call get_command_argument(1,argument)
  if(len_trim(argument)>0)read(argument,*)nq
  n=[16,12,8];h=[.7d0,.9d0,1.1d0]
  allocate(rho(product(n)),sigma(product(n)),e(product(n)),v(product(n)),w(product(n)))
  do k=0,n(3)-1;do j=0,n(2)-1;do i=0,n(1)-1
    g=1+i+n(1)*(j+n(2)*k)
    rho(g)=.01d0+.007d0*cos(2*acos(-1d0)*i/n(1))+.001d0*sin(2*acos(-1d0)*j/n(2)) &
      +.0007d0*sin(2*acos(-1d0)*k/n(3)) &
      +.0003d0*cos(2*acos(-1d0)*(real(i,8)/n(1)+real(j,8)/n(2)+real(k,8)/n(3)))
    sigma(g)=.0002d0*(1+.7d0*cos(2*acos(-1d0)*k/n(3)))
  enddo;enddo;enddo
  allocate(density(n(1),n(2),n(3)),grad(n(1),n(2),n(3),3),ep(n(1),n(2),n(3)), &
    vp(n(1),n(2),n(3)),flux(n(1),n(2),n(3)),potential(n(1),n(2),n(3)))
  allocate(local_parts(product(n),2),parts(product(n),2))
  density=reshape(rho,n);grad=0d0;matrix=0d0
  if(maxval(abs(density(:,:,2)-density(:,:,1)))<1d-8)error stop 'fixture lacks z modes'
  do i=1,3
    matrix(i,i)=1d0;coef(:,i)=[.8d0,-.2d0,4d0/105,-1d0/280]/h(i)
    do j=1,4
      grad(:,:,:,i)=grad(:,:,:,i)+coef(j,i)*(cshift(density,j,i)-cshift(density,-j,i))
    enddo
  enddo
  sigma=reshape(sum(grad**2,dim=4),[product(n)])
  call rvv10_periodic(n,h,coef,matrix,density,5.3d0,.0093d0,nq,ep,vp,status)
  if(status/=0)error stop 'serial potential failed'
  call rvv10_evaluate(n,h,rho,sigma,5.3d0,.0093d0,nq,e,v,w,status)
  if(status/=0)error stop 'serial failed'
  do axis=1,6
    if(axis>3.and.nproc/=4)cycle
    dims=1
    if(axis<=3)then
      dims(axis)=nproc
    else if(axis==4)then
      dims=[1,2,2]
    else if(axis==5)then
      dims=[2,2,1]
    else
      dims=[2,1,2]
    endif
    coords=[modulo(rank,dims(1)),modulo(rank/dims(1),dims(2)),rank/(dims(1)*dims(2))]
    do i=1,3
      color=0
      do j=1,3
        if(j/=i)color=color*nproc+coords(j)
      enddo
      key=coords(i)
      call MPI_Comm_split(MPI_COMM_WORLD,color,key,comm(i),ierr)
    enddo
    local_n=n/dims;lo=coords*local_n+1
    allocate(batch_input(n(1)*local_n(2)*local_n(3),5), &
      batch_output(n(1)*local_n(2)*local_n(3),5),batch_back(n(1)*local_n(2)*local_n(3),5), &
      fft_a(n(1)*local_n(2)*local_n(3)),fft_b(n(1)*local_n(2)*local_n(3)))
    do channel=1,5;do i=1,size(batch_input,1)
      batch_input(i,channel)=cmplx(sin(.17d0*i+channel+rank),cos(.09d0*i-channel+rank),8)
    enddo;enddo
    do batch_count=2,5
    call pencil_transform(n,dims(2:3),coords(2:3),comm(2:3),batch_input(:,:batch_count),batch_output(:,:batch_count),-1,status)
    if(status/=0)error stop 'FFTW pencils failed'
    plans_before=fftw_pencil_plans_created
    call pencil_transform(n,dims(2:3),coords(2:3),comm(2:3),batch_output(:,:batch_count),batch_back(:,:batch_count),1,status)
    if(status/=0.or.maxval(abs(batch_back(:,:batch_count)-batch_input(:,:batch_count)))>1d-11)error stop 'FFTW pencil round trip'
    if(plans_before/=fftw_pencil_plans_created)error stop 'plans were not reused'
    do channel=1,batch_count
      fft_a=batch_input(:,channel)
      call pzfft3dv_rvv10(fft_a,fft_b,n(1),n(2),n(3),dims(2),dims(3),-1,comm(2),comm(3))
      if(maxval(abs(fft_b-batch_output(:,channel)))>1d-10)error stop 'FFTW vs FFTE transform'
    enddo
    enddo
    deallocate(batch_input,batch_output,batch_back,fft_a,fft_b)

    allocate(lr(product(local_n)),ls(product(local_n)),le(product(local_n)),lv(product(local_n)),lw(product(local_n)))
    l=0
    do k=lo(3),lo(3)+local_n(3)-1;do j=lo(2),lo(2)+local_n(2)-1;do i=lo(1),lo(1)+local_n(1)-1
      l=l+1;g=i+n(1)*((j-1)+n(2)*(k-1));lr(l)=rho(g);ls(l)=sigma(g)
    enddo;enddo;enddo
    allocate(old_input(8*(8/dims(2))*(8/dims(3))),old_a(8*(8/dims(2))*(8/dims(3))), &
      old_b(8*(8/dims(2))*(8/dims(3))),old_reference(8*(8/dims(2))*(8/dims(3))))
    do i=1,size(old_input)
      old_input(i)=cmplx(sin(real(i+rank,8)),cos(real(2*i+rank,8)),8)
    enddo
    old_a=old_input
    call pzfft3dv_mod(old_a,old_b,8,8,8,dims(2),dims(3),0,comm(2),comm(3))
    call pzfft3dv_mod(old_a,old_b,8,8,8,dims(2),dims(3),-1,comm(2),comm(3))
    old_reference=old_b
    call rvv10_evaluate_distributed(n,lo,local_n,dims,coords,comm,MPI_COMM_WORLD,h,lr,ls, &
      5.3d0,.0093d0,nq,le,lv,lw,used,status,.true.)
    if(.not.used.or.status/=0)error stop 'distributed unavailable'
    err=0;l=0
    do k=lo(3),lo(3)+local_n(3)-1;do j=lo(2),lo(2)+local_n(2)-1;do i=lo(1),lo(1)+local_n(1)-1
      l=l+1;g=i+n(1)*((j-1)+n(2)*(k-1))
      err=max(err,abs(le(l)-e(g)),abs(lv(l)-v(g)),abs(lw(l)-w(g)))
    enddo;enddo;enddo
    call MPI_Allreduce(err,total_err,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(rank==0)write(*,*)axis,total_err
    if(total_err>1d-10)error stop 'distributed differs from serial'
    local_parts=0d0;l=0
    do k=lo(3),lo(3)+local_n(3)-1;do j=lo(2),lo(2)+local_n(2)-1;do i=lo(1),lo(1)+local_n(1)-1
      l=l+1;g=i+n(1)*((j-1)+n(2)*(k-1));local_parts(g,:)=[lv(l),lw(l)]
    enddo;enddo;enddo
    call MPI_Allreduce(local_parts,parts,size(parts),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    potential=reshape(parts(:,1),n)
    do i=1,3
      flux=2*reshape(parts(:,2),n)*grad(:,:,:,i)
      do j=1,4
        potential=potential-coef(j,i)*(cshift(flux,j,i)-cshift(flux,-j,i))
      enddo
    enddo
    err=maxval(abs(potential-vp))
    if(rank==0)write(*,*)'Full potential comparison:',err
    if(err>1d-10)error stop 'full potential differs from serial'
    old_a=old_input
    call pzfft3dv_mod(old_a,old_b,8,8,8,dims(2),dims(3),-1,comm(2),comm(3))
    if(maxval(abs(old_b-old_reference))>1d-12)error stop 'Poisson FFT tables overwritten'
    deallocate(old_input,old_a,old_b,old_reference)
    ! Unsupported geometry must return consistently before collective FFTs.
    call rvv10_evaluate_distributed([14,12,8],lo,local_n,dims,coords,comm,MPI_COMM_WORLD,h,lr,ls, &
      5.3d0,.0093d0,nq,le,lv,lw,used,status,.true.)
    if(used.or.status/=0)error stop 'unsupported layout did not fall back'
    ! An invalid density on one rank must fail on all ranks without deadlock.
    if(rank==0)lr(1)=-1d0
    call rvv10_evaluate_distributed(n,lo,local_n,dims,coords,comm,MPI_COMM_WORLD,h,lr,ls, &
      5.3d0,.0093d0,nq,le,lv,lw,used,status,.true.)
    if(.not.used.or.status==0)error stop 'invalid density accepted'
    do i=1,3
      call MPI_Comm_free(comm(i),ierr)
    enddo
    deallocate(lr,ls,le,lv,lw)
  enddo
  call MPI_Finalize(ierr)
end program
