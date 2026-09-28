program exchange_driver
  use mpi
  use iso_fortran_env, only: int64
  use fftw_pencils, only: fftw_pencil_transposes
  use hse_ace
  use exx_orbitals
  use hse_spatial
  use hse_wannier
  use hse_wannier_gauge, only: gauge_transport
  implicit none
  integer,parameter :: n(3)=[8,8,8],no=3,nt=5
  real(8),parameter :: h(3)=[.7d0,.8d0,.9d0]
  complex(8) :: psi(product(n),no,1),target(product(n),nt,1),reference(product(n),nt,1)
  complex(8),allocatable :: local(:,:,:),trial(:,:,:),action(:,:,:),previous_saved(:,:,:)
  type(spatial_exx_state) :: spatial,full_spatial,partitioned,masked,masked_partitioned
  type(s_hse_wannier) :: serial
  type(hse_ace_state) :: distributed_ace,reference_ace,old_ace,average_ace
  complex(8),allocatable :: local_w(:,:,:),ace_ref(:,:,:),ace_result(:,:,:),old_action(:,:,:)
  integer(int64) :: before_transposes
  integer :: np,rank,err,status,dims(2),coords(2),comm(2),m(3),lo(3),g,l,x,y,z,j,k,stage
  integer :: bad_coords(2),comm_r,comm_o,orb_rank,orb_size,spatial_rank,spatial_size,first_o,last_o,first_t,last_t
  real(8) :: error,global_error,dv,bad_dv,minimum,omega,occupation(no,1)
  character(32) :: argument
  complex(8) :: transported(no,no,1),metric(nt,nt),metric_sum(nt,nt)
  integer :: point(3),center(3)
  call get_command_argument(1,argument)
  read(argument,*)omega
  call MPI_Init(err)
  call MPI_Comm_size(MPI_COMM_WORLD,np,err)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,err)
  orb_size=1
  call get_command_argument(2,argument)
  if(len_trim(argument)>0)read(argument,*)orb_size
  if(modulo(np,orb_size)/=0)error stop 'invalid partition'
  spatial_size=np/orb_size;orb_rank=rank/spatial_size;spatial_rank=modulo(rank,spatial_size)
  call MPI_Comm_split(MPI_COMM_WORLD,orb_rank,spatial_rank,comm_r,err)
  call MPI_Comm_split(MPI_COMM_WORLD,spatial_rank,orb_rank,comm_o,err)
  dims=[spatial_size,1]
  if(spatial_size==4)dims=[2,2]
  coords=[modulo(spatial_rank,dims(1)),spatial_rank/dims(1)]
  call MPI_Comm_split(comm_r,coords(2),coords(1),comm(1),err)
  call MPI_Comm_split(comm_r,coords(1),coords(2),comm(2),err)
  m=[n(1),n(2)/dims(1),n(3)/dims(2)];lo=[0,coords(1)*m(2),coords(2)*m(3)]
  dv=product(h)
  do j=1,no
    do g=1,product(n)
      psi(g,j,1)=cmplx(sin(.07d0*g*j)+cos(.03d0*g*(j+1)),sin(.05d0*g*(j+2)),8)
    enddo
    do k=1,j-1
      psi(:,j,1)=psi(:,j,1)-psi(:,k,1)*sum(conjg(psi(:,k,1))*psi(:,j,1))*dv
    enddo
    psi(:,j,1)=psi(:,j,1)/sqrt(sum(abs(psi(:,j,1))**2)*dv)
  enddo
  do j=1,nt
    do g=1,product(n)
      target(g,j,1)=cmplx(cos(.013d0*g*j),sin(.017d0*g*(j+1)),8)
    enddo
  enddo
  allocate(local(product(m),no,1),trial(product(m),nt,1),action(product(m),nt,1))
  call wannier_init(serial,n,[1,1,1],h,reshape([0d0,0d0,0d0],[3,1]),omega,status,2.5d0)
  if(status/=0)error stop 'serial init'
  do stage=1,4
    occupation=2d0
    if(stage==2)occupation(:,1)=[2d0,.7d0,0d0]
    if(stage==3)occupation(:,1)=[1.3d0,.2d0,.1d0]
    if(stage==4)occupation=0d0
    if(stage==2)psi=psi*cmplx(cos(.13d0),sin(.13d0),8)
    call wannier_refresh_source(serial,psi,occupation,3,1d-7,status)
    if(status/=0)error stop 'serial refresh'
    call wannier_apply(serial,target,reference,status)
    if(status/=0)error stop 'serial action'
    l=0
    do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
      l=l+1;g=1+x+n(1)*(y+lo(2)+n(2)*(z+lo(3)))
      local(l,:,1)=psi(g,:,1);trial(l,:,1)=target(g,:,1)
    enddo;enddo;enddo
    call spatial_exx_refresh(spatial,n,h,dims,coords,comm,comm_r,local, &
      merge(3,0,stage==1),1d-7,status,occupation=occupation)
    if(status/=0)error stop 'spatial refresh'
    call spatial_exx_refresh(full_spatial,n,h,[1,1],[0,0],[MPI_COMM_SELF,MPI_COMM_SELF], &
      MPI_COMM_SELF,psi,merge(3,0,stage==1),1d-7,status,occupation=occupation)
    if(status/=0)error stop 'full grid localization'
    l=0
    do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
      l=l+1;g=1+x+n(1)*(y+lo(2)+n(2)*(z+lo(3)))
      if(maxval(abs(spatial%source(l,:)-full_spatial%source(g,:)))>1d-10)error stop 'localization mismatch'
    enddo;enddo;enddo
    if(stage==1)then
      previous_saved=spatial%previous
    else
      if(abs(spatial%min_singular-1d0)>1d-10)error stop 'transport overlap mismatch'
      if(maxval(abs(spatial%previous-previous_saved))>1d-10)error stop 'transport gauge mismatch'
    endif
    before_transposes=fftw_pencil_transposes
    call spatial_exx_apply(spatial,n,h,dims,coords,comm,comm_r,2.5d0,trial,action,status,omega=omega)
    if(fftw_pencil_transposes-before_transposes/=4*no*((nt+3)/4))error stop 'redundant FFT transpose'
    if(status/=0)error stop 'spatial action'
    first_o=no*orb_rank/orb_size+1;last_o=no*(orb_rank+1)/orb_size
    first_t=nt*orb_rank/orb_size+1;last_t=nt*(orb_rank+1)/orb_size
    call spatial_exx_refresh(partitioned,n,h,dims,coords,comm,comm_r, &
      local(:,first_o:last_o,:),merge(3,0,stage==1),1d-7,status, &
      occupation=occupation(first_o:last_o,:),comm_o=comm_o)
    if(status/=0)error stop 'distributed localization'
    if(any(shape(partitioned%source)/=[product(m),last_o-first_o+1]))error stop 'replicated source columns'
    if(any(abs(partitioned%previous-spatial%previous(:,first_o:last_o,:))>1d-9))error stop 'distributed transport'
    if(any(abs(partitioned%source-spatial%source(:,first_o:last_o))>1d-9))error stop 'distributed source'
    if(allocated(local_w))deallocate(local_w,ace_ref,ace_result)
    allocate(local_w(product(m),no,1),ace_ref(product(m),nt,1),ace_result(product(m),last_t-first_t+1,1))
    call spatial_exx_apply(spatial,n,h,dims,coords,comm,comm_r,2.5d0,local,local_w,status,omega=omega)
    if(status/=0)error stop 'reference source action'
    call hse_ace_build(reference_ace,local,local_w,dv,status,sum_grid)
    if(status/=0)error stop 'reference ACE build'
    call orbital_ace_build(distributed_ace,local(:,first_o:last_o,:),local_w(:,first_o:last_o,:), &
      dv,comm_r,comm_o,status)
    if(status/=0)error stop 'distributed ACE build'
    if(any(shape(distributed_ace%factors)/=[product(m),last_o-first_o+1,1]))error stop 'replicated ACE columns'
    call hse_ace_apply(reference_ace,trial,ace_ref,status,sum_grid)
    if(status/=0)error stop 'reference ACE apply'
    call orbital_ace_apply(distributed_ace,trial(:,first_t:last_t,:),ace_result,comm_r,comm_o,status)
    if(status/=0)error stop 'distributed ACE apply'
    if(any(abs(ace_result-ace_ref(:,first_t:last_t,:))>1d-10))error stop 'ACE target mismatch'
    ! ACE must reproduce exact exchange on every source column.
    call orbital_ace_apply(distributed_ace,local(:,first_o:last_o,:),action(:,:last_o-first_o+1,:), &
      comm_r,comm_o,status)
    if(status/=0)error stop 'ACE source apply'
    if(any(abs(action(:,:last_o-first_o+1,:)-local_w(:,first_o:last_o,:))>1d-10))error stop 'ACE source mismatch'
    if(stage>1)then
      call hse_ace_average(old_ace,distributed_ace,average_ace,status)
      if(status/=0)error stop 'distributed ACE average'
      call orbital_ace_apply(average_ace,trial(:,first_t:last_t,:),ace_result,comm_r,comm_o,status)
      if(status/=0)error stop 'averaged ACE apply'
      if(any(abs(ace_result-.5d0*(old_action+ace_ref(:,first_t:last_t,:)))>1d-10))error stop 'ACE average mismatch'
    endif
    old_ace=distributed_ace;old_action=ace_ref(:,first_t:last_t,:)
    ! Recompute the reference outside the distributed slice overwritten above.
    call spatial_exx_apply(spatial,n,h,dims,coords,comm,comm_r,2.5d0,trial,action,status,omega=omega)
    call spatial_exx_apply(partitioned,n,h,dims,coords,comm,comm_r,2.5d0, &
      trial(:,first_t:last_t,:),action(:,first_t:last_t,:),status,omega=omega,comm_o=comm_o)
    if(status/=0)error stop 'orbital partition action'
    error=0d0;l=0
    do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
      l=l+1;g=1+x+n(1)*(y+lo(2)+n(2)*(z+lo(3)))
      error=max(error,maxval(abs(action(l,:,1)-reference(g,:,1))))
    enddo;enddo;enddo
    call MPI_Allreduce(error,global_error,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,err)
    if(global_error>1d-10)error stop 'exchange mismatch'
    if(rank==0)print *, 'PASS spatial exchange ranks/stage/error ',np,stage,global_error,omega
  enddo
  ! The same masked sources must give the same action with global and compact
  ! convolution, including orbital streaming with uneven/empty target groups.
  allocate(masked%source(product(m),no),masked_partitioned%source(product(m),last_o-first_o+1))
  l=0
  do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
    l=l+1;point=[x,y,z]+lo
    do j=1,no
      center=[j-1,0,0]
      masked%source(l,j)=0d0
      if(all(modulo(point-center+1,n)<=1)) &
        masked%source(l,j)=cmplx(.1d0*j+.003d0*sum(point),.07d0*j,8)
    enddo
  enddo;enddo;enddo
  masked_partitioned%source=masked%source(:,first_o:last_o)
  masked%compact=.false.
  call spatial_exx_apply(masked,n,h,dims,coords,comm,comm_r,2.5d0,trial,ace_ref,status,omega=omega)
  if(status/=0)error stop 'masked global reference'
  masked%compact=.true.
  call spatial_exx_apply(masked,n,h,dims,coords,comm,comm_r,2.5d0,trial,action,status,omega=omega)
  if(status/=0)error stop 'masked compact reference'
  if(maxval(abs(action-ace_ref))>1d-11)error stop 'compact global mismatch'
  if(masked%local_pairs/=no*nt.or.masked%global_pairs/=0)error stop 'compact path not used'
  if(masked%local_points>=int(no*nt*product(n),int64))error stop 'compact FFT volume not reduced'
  metric=matmul(conjg(transpose(trial(:,:,1))),action(:,:,1))
  call MPI_Allreduce(metric,metric_sum,nt*nt,MPI_DOUBLE_COMPLEX,MPI_SUM,comm_r,err)
  if(maxval(abs(metric_sum-conjg(transpose(metric_sum))))>1d-10)error stop 'compact Hermiticity'
  masked_partitioned%compact=.false.
  call spatial_exx_apply(masked_partitioned,n,h,dims,coords,comm,comm_r,2.5d0, &
    trial(:,first_t:last_t,:),ace_result,status,omega=omega,comm_o=comm_o)
  if(status/=0)error stop 'masked orbital global'
  if(any(abs(ace_result-ace_ref(:,first_t:last_t,:))>1d-11))error stop 'masked orbital global mismatch'
  masked_partitioned%compact=.true.
  call spatial_exx_apply(masked_partitioned,n,h,dims,coords,comm,comm_r,2.5d0, &
    trial(:,first_t:last_t,:),ace_result,status,omega=omega,comm_o=comm_o)
  if(status/=0)error stop 'masked orbital compact'
  if(any(abs(ace_result-ace_ref(:,first_t:last_t,:))>1d-11))error stop 'masked orbital compact mismatch'
  if(masked_partitioned%local_pairs/=no*(last_t-first_t+1).or.masked_partitioned%global_pairs/=0) &
    error stop 'orbital compact path not used'
  if(rank==0)print *, 'PASS compact orbital exchange ranks/orbitals ',np,orb_size,omega
  if(omega>0d0)then
    ! Exercise all-dropped pairs with uneven/empty source and target owners.
    ! Huge budget is a stress test, not a recommended physical tolerance.
    masked_partitioned%screen_tolerance=1d6
    do stage=1,2
      masked_partitioned%screen_mode=stage
      call spatial_exx_apply(masked_partitioned,n,h,dims,coords,comm,comm_r,2.5d0, &
        trial(:,first_t:last_t,:),ace_result,status,omega=omega,comm_o=comm_o)
      if(status/=0.or.masked_partitioned%screen_candidates/=no*nt)error stop 'partitioned pair count'
      if(stage==1)then
        if(any(abs(ace_result-ace_ref(:,first_t:last_t,:))>1d-11))error stop 'pair diagnosis changed action'
      else
        if(masked_partitioned%screen_skipped/=no*nt)error stop 'partitioned skipped count'
        error=sum(abs(ace_result-ace_ref(:,first_t:last_t,:))**2)*dv
        call MPI_Allreduce(error,global_error,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,err)
        if(sqrt(global_error)>masked_partitioned%screen_bound+1d-10)error stop 'partitioned bound'
        if(any(ace_result/=(0d0,0d0)))error stop 'all-dropped action not zero'
      endif
    enddo
    masked_partitioned%screen_mode=0
    if(rank==0)print *, 'PASS pair screening orbital/empty layouts ',np,orb_size,omega
  endif
  call spatial_exx_apply(spatial,n,h,dims,coords,comm,comm_r,2.5d0,trial,action,status,omega=-.1d0)
  if(status==0)error stop 'negative screening accepted'
  ! Force an FFT validation failure after the first source broadcast.
  ! Last orbital group always has a target, including layouts with empty peers.
  bad_coords=coords
  if(orb_rank==orb_size-1)bad_coords(1)=dims(1)
  call spatial_exx_apply(partitioned,n,h,dims,bad_coords,comm,comm_r,2.5d0, &
    trial(:,first_t:last_t,:),action(:,first_t:last_t,:),status,omega=omega,comm_o=comm_o)
  if(status==0)error stop 'FFT failure not propagated'
  ! One bad orbital group must make every group return before source broadcasts.
  if(orb_rank==0)deallocate(partitioned%source)
  call spatial_exx_apply(partitioned,n,h,dims,coords,comm,comm_r,2.5d0, &
    trial(:,first_t:last_t,:),action(:,first_t:last_t,:),status,omega=omega,comm_o=comm_o)
  if(status==0)error stop 'invalid orbital source accepted'
  bad_dv=dv
  if(spatial_rank==0)bad_dv=-1d0
  call gauge_transport(local,spatial%previous,bad_dv,transported,minimum,status,sum_grid)
  if(status==0)error stop 'invalid local transport volume accepted'
  bad_dv=dv
  if(rank==0)bad_dv=-1d0
  call orbital_ace_build(distributed_ace,local(:,first_o:last_o,:),local_w(:,first_o:last_o,:), &
    bad_dv,comm_r,comm_o,status)
  if(status==0)error stop 'invalid orbital ACE volume accepted'
  ! Build failure clears factors collectively, so subsequent action must fail too.
  call orbital_ace_apply(distributed_ace,trial(:,first_t:last_t,:),ace_result,comm_r,comm_o,status)
  if(status==0)error stop 'uninitialized orbital ACE accepted'
  call wannier_destroy(serial)
  call MPI_Comm_free(comm(1),err);call MPI_Comm_free(comm(2),err)
  call MPI_Comm_free(comm_r,err);call MPI_Comm_free(comm_o,err)
  call MPI_Finalize(err)
contains
  subroutine sum_grid(a)
    complex(8),intent(inout) :: a(:,:)
    integer :: code
    call MPI_Allreduce(MPI_IN_PLACE,a,size(a),MPI_DOUBLE_COMPLEX,MPI_SUM,comm_r,code)
  end subroutine
end program
