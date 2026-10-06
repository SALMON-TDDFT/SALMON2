! CPU oracle for the stateful prepare/apply/release contract.
module resident_oracle
 use iso_c_binding
 use exx_batch_backend, only: s_exx_batch_backend
 use exx_local_fft, only: s_exx_local_fft,exx_local_apply,exx_local_destroy
 implicit none
 include 'fftw3.f03'
 integer :: releases=0
 type,extends(s_exx_batch_backend) :: cpu_backend
  type(s_exx_local_fft) :: fft
  complex(8),allocatable :: source(:)
  integer :: prepares=0,batches=0,capacity=0
  logical :: fail_prepare=.false.,fail_apply=.false.,fail_release=.false.
 contains
  procedure :: prepare=>prepare_cpu
  procedure :: apply=>apply_cpu
  procedure :: release=>release_cpu
 end type
contains
 subroutine prepare_cpu(self,padded,indices,filter,source,capacity,status)
  implicit none
  class(cpu_backend),target,intent(inout) :: self
  integer,intent(in) :: padded(3),indices(:),capacity
  complex(8),intent(in) :: filter(:,:,:),source(:)
  integer,intent(out) :: status
  self%prepares=self%prepares+1;status=1
  if(self%fail_prepare)return
  call exx_local_destroy(self%fft)
  self%capacity=capacity;self%source=source
  self%fft%padded=padded;self%fft%fft_points=product(padded)
  self%fft%indices=indices;self%fft%filter=filter
  allocate(self%fft%work(padded(1),padded(2),padded(3)))
  self%fft%forward=fftw_plan_dft_3d(padded(3),padded(2),padded(1), &
    self%fft%work,self%fft%work,FFTW_FORWARD,FFTW_ESTIMATE)
  self%fft%backward=fftw_plan_dft_3d(padded(3),padded(2),padded(1), &
    self%fft%work,self%fft%work,FFTW_BACKWARD,FFTW_ESTIMATE)
  if(.not.c_associated(self%fft%forward).or..not.c_associated(self%fft%backward))return
  self%fft%ready=.true.;status=0
 end subroutine
 subroutine apply_cpu(self,targets,action,status)
  implicit none
  class(cpu_backend),target,intent(inout) :: self
  complex(8),intent(in) :: targets(:,:)
  complex(8),intent(out) :: action(:,:)
  integer,intent(out) :: status
  complex(8) :: potential(size(targets,1))
  integer :: j
  self%batches=self%batches+1;status=1;action=0d0
  if(self%fail_apply)return
  if(size(targets,2)>self%capacity)error stop 'stateful capacity'
  do j=1,size(targets,2)
   call exx_local_apply(self%fft,conjg(self%source)*targets(:,j),potential,status)
   if(status/=0)return
   action(:,j)=-self%source*potential
  enddo
  status=0
 end subroutine
 subroutine release_cpu(self,status)
  implicit none
  class(cpu_backend),target,intent(inout) :: self
  integer,intent(out) :: status
  call exx_local_destroy(self%fft)
  if(allocated(self%source))deallocate(self%source)
  releases=releases+1;status=0
  if(self%fail_release)status=-1001
 end subroutine
end module
program probe
 use resident_oracle, only: cpu_backend,releases
 use mpi
 use iso_fortran_env, only: int64
 use iso_c_binding
 use exx_local_fft, only: s_exx_local_fft,exx_local_apply,exx_local_destroy
 use fftw_blocks, only: mesh_transform
 use exx_spatial_local
 implicit none
 include 'fftw3.f03'
 type(s_exx_spatial_local) :: plan
 integer :: n(3)=[16,12,8],dims(2),coords(2),comm(2),m(3),lo(3),p(3),rank,np,ierr,g,x,y,z,j,status,kind
 integer(int64) :: pairs,points,reference_points
 complex(8),allocatable :: source(:),target(:,:),action(:,:),ref(:,:),work(:,:),spec(:,:)
 real(8),allocatable :: multiplier(:)
 complex(8) :: metric(7,7),total(7,7)
 real(8) :: q2,err,err_all
 logical :: used,skip(7),oracle_fail=.false.
 integer :: batch_sizes(4)=[1,2,3,8],ibatch,callback_calls,all_calls
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 dims=[np,1];if(np==4)dims=[2,2]
 coords=[mod(rank,dims(1)),rank/dims(1)]
 call MPI_Comm_split(MPI_COMM_WORLD,coords(2),coords(1),comm(1),ierr)
 call MPI_Comm_split(MPI_COMM_WORLD,coords(1),coords(2),comm(2),ierr)
 m=[n(1),n(2)/dims(1),n(3)/dims(2)];lo=[0,coords(1)*m(2),coords(2)*m(3)]
 allocate(source(product(m)),target(product(m),7),action(product(m),7),ref(product(m),7), &
 work(product(m),7),spec(product(m),7),multiplier(product(m)))
 g=0
 do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
 g=g+1;p=[x,y,z]+lo
 where(p>=(n+1)/2)p=p-n
 q2=sum((2*acos(-1d0)*p/n)**2);multiplier(g)=1/(.7d0+q2)
 enddo;enddo;enddo
 call spatial_local_init(plan,n,dims,coords,comm,multiplier,status)
 if(status/=0)error stop 'init'
 do kind=1,3
 g=0
 do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
 g=g+1;p=[x,y,z]+lo
 source(g)=0
 if(all(modulo(p+1,n)<=2))source(g)=cmplx(.3d0+.01d0*sum(p),.2d0,8)
 if(kind==2)source(g)=source(g)*exp(cmplx(0d0,.83d0,8))
 if(kind==3)source(g)=1
 do j=1,6
 target(g,j)=cmplx(sin(real(sum(p*[2,3,5])+j,8)),cos(real(sum(p*[5,1,2])+2*j,8)),8)
 enddo
 target(g,7)=0
 enddo;enddo;enddo
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points)
 if(status/=0)error stop 'apply'
 if(kind==3)then
 if(used.or.any(action/=(0d0,0d0)))error stop 'broad fallback'
 cycle
 endif
 if(.not.used.or.pairs/=6.or.points>=6*product(n))error stop 'compact counter'
 do j=1,7
 work(:,j)=conjg(source)*target(:,j)
 enddo
 call mesh_transform(n,dims,coords,comm,work(:,1:4),spec(:,1:4),-1,status)
 call mesh_transform(n,dims,coords,comm,work(:,5:7),spec(:,5:7),-1,status)
 if(status/=0)error stop 'forward'
 do j=1,7
 spec(:,j)=spec(:,j)*multiplier
 enddo
 call mesh_transform(n,dims,coords,comm,spec(:,1:4),work(:,1:4),1,status)
 call mesh_transform(n,dims,coords,comm,spec(:,5:7),work(:,5:7),1,status)
 if(status/=0)error stop 'inverse'
 do j=1,7
 ref(:,j)=-source*work(:,j)
 enddo
 err=maxval(abs(ref-action));call MPI_Allreduce(err,err_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 if(err_all>1d-11)error stop 'global parity'
 metric=matmul(conjg(transpose(target)),action)
 call MPI_Allreduce(metric,total,49,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
 if(maxval(abs(total-conjg(transpose(total))))>1d-10)error stop 'Hermiticity'
 enddo
 source=0
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points)
 if(status/=0.or..not.used.or.pairs/=0.or.any(action/=(0d0,0d0)))error stop 'zero support'
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target(:,1:0),action(:,1:0),used,status,pairs,points)
 if(status/=0.or..not.used.or.pairs/=0)error stop 'empty targets'
 ! A changed displacement tile (single point) must replace the previous plan.
 source=0
 if(all(lo==0))source(1)=cmplx(.7d0,.4d0,8)
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points)
 if(status/=0.or..not.used.or.pairs/=6.or.points/=6)error stop 'single point'
 do j=1,7
 work(:,j)=conjg(source)*target(:,j)
 enddo
 call mesh_transform(n,dims,coords,comm,work(:,1:4),spec(:,1:4),-1,status)
 call mesh_transform(n,dims,coords,comm,work(:,5:7),spec(:,5:7),-1,status)
 do j=1,7
 spec(:,j)=spec(:,j)*multiplier
 enddo
 call mesh_transform(n,dims,coords,comm,spec(:,1:4),work(:,1:4),1,status)
 call mesh_transform(n,dims,coords,comm,spec(:,5:7),work(:,5:7),1,status)
 do j=1,7
 ref(:,j)=-source*work(:,j)
 enddo
 err=maxval(abs(ref-action));call MPI_Allreduce(err,err_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 if(err_all>1d-11)error stop 'single point parity'
 ! Exercise optional batched dispatch using the production FFTW scalar oracle.
 ! Seven columns, two explicit skips and a zero column leave partial owner batches.
 g=0
 do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
 g=g+1;p=[x,y,z]+lo
 source(g)=0d0
 if(all(modulo(p+1,n)<=2))source(g)=cmplx(.3d0+.01d0*sum(p),.2d0,8)
 enddo;enddo;enddo
 skip=.false.;skip(2)=.true.;skip(6)=.true.
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,ref,used,status,pairs,reference_points,skip)
 if(status/=0.or..not.used.or.pairs/=4)error stop 'CPU skipped-pair reference'
 plan%batch_action=>cpu_batch_action
 do ibatch=1,size(batch_sizes)
 plan%batch_size=batch_sizes(ibatch);callback_calls=0
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points,skip)
 if(status/=0.or..not.used.or.pairs/=4.or.points/=reference_points)error stop 'batched counters'
 call MPI_Allreduce(callback_calls,all_calls,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
 if(all_calls<1)error stop 'batch callback not exercised'
 if(any(action(:,2)/=(0d0,0d0)).or.any(action(:,6:7)/=(0d0,0d0)))error stop 'skip/zero columns'
 err=maxval(abs(ref-action));call MPI_Allreduce(err,err_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 if(err_all>1d-11)error stop 'batched MPI scatter parity'
 enddo
 ! A backend error on one owner must reach all ranks instead of hanging or passing.
 oracle_fail=(rank==0)
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points,skip)
 if(status==0)error stop 'batch callback error was not propagated'
 oracle_fail=.false.
 nullify(plan%batch_action)
 allocate(cpu_backend::plan%backend)
 do ibatch=1,size(batch_sizes)
 plan%batch_size=batch_sizes(ibatch)
 select type(backend=>plan%backend)
 type is(cpu_backend)
 backend%prepares=0;backend%batches=0
 end select
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points,skip)
 if(status/=0.or..not.used.or.pairs/=4.or.points/=reference_points)error stop 'resident counters'
 select type(backend=>plan%backend)
 type is(cpu_backend)
 if(backend%prepares/=1)error stop 'prepare must be once per source, outside batches'
 if(ibatch==1.and.rank==0.and.backend%batches<2)error stop 'multiple batches not exercised'
 end select
 err=maxval(abs(ref-action));call MPI_Allreduce(err,err_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 if(err_all>1d-11)error stop 'resident MPI scatter parity'
 enddo
 ! Reusing the object with a changed source must not use a stale source.
 source=2d0*source
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points,skip)
 if(status/=0)error stop 'resident new source'
 err=maxval(abs(4d0*ref-action));call MPI_Allreduce(err,err_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 if(err_all>1d-11)error stop 'resident source update'
 select type(backend=>plan%backend)
 type is(cpu_backend)
 backend%fail_prepare=(rank==0)
 end select
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points,skip)
 if(status==0)error stop 'resident prepare error propagation'
 select type(backend=>plan%backend)
 type is(cpu_backend)
 backend%fail_prepare=.false.;backend%fail_apply=(rank==0)
 end select
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points,skip)
 if(status==0)error stop 'resident apply error propagation'
 select type(backend=>plan%backend)
 type is(cpu_backend)
 backend%fail_apply=.false.;backend%prepares=0;backend%batches=0
 end select
 callback_calls=0
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target(:,1:0),action(:,1:0),used,status,pairs,points)
 if(status/=0.or..not.used.or.pairs/=0.or.callback_calls/=0)error stop 'batch empty targets'
 source=0d0
 call spatial_local_apply(plan,MPI_COMM_WORLD,source,target,action,used,status,pairs,points)
 if(status/=0.or..not.used.or.pairs/=0.or.callback_calls/=0)error stop 'batch zero support'
 select type(backend=>plan%backend)
 type is(cpu_backend)
 if(backend%prepares/=0.or.backend%batches/=0)error stop 'empty work prepared resident backend'
 end select
 call spatial_local_destroy(plan,status)
 if(status/=0.or.releases/=1.or.allocated(plan%backend))error stop 'resident release'
 call spatial_local_destroy(plan,status)
 if(status/=0.or.releases/=1)error stop 'repeat destroy'
 allocate(cpu_backend::plan%backend)
 select type(backend=>plan%backend)
 type is(cpu_backend)
 backend%fail_release=(rank==0)
 end select
 call spatial_local_init(plan,n,dims,coords,comm,multiplier,status)
 if(status==0)error stop 'negative reinit release error lost on peer'
 if(allocated(plan%backend).or.releases/=2)error stop 'reinit release lifecycle'
 call spatial_local_destroy(plan)
 if(rank==0)print *, 'PASS spatial local and batched callback (resident prepare/apply/release)',np
 call MPI_Finalize(ierr)
contains
 subroutine cpu_batch_action(padded,indices,filter,source,targets,action,status)
  implicit none
  integer,intent(in) :: padded(3),indices(:)
  complex(8),intent(in) :: filter(:,:,:),source(:),targets(:,:)
  complex(8),intent(out) :: action(:,:)
  integer,intent(out) :: status
  type(s_exx_local_fft) :: oracle
  complex(8),allocatable :: density(:),potential(:)
  integer :: column
  callback_calls=callback_calls+1
  if(size(targets,2)>plan%batch_size)error stop 'backend batch exceeds configured maximum'
  status=1;action=0d0
  if(oracle_fail)return
  oracle%padded=padded;oracle%fft_points=product(padded);oracle%indices=indices
  oracle%filter=filter
  allocate(oracle%work(padded(1),padded(2),padded(3)),density(size(source)),potential(size(source)))
  oracle%forward=fftw_plan_dft_3d(padded(3),padded(2),padded(1), &
    oracle%work,oracle%work,FFTW_FORWARD,FFTW_ESTIMATE)
  oracle%backward=fftw_plan_dft_3d(padded(3),padded(2),padded(1), &
    oracle%work,oracle%work,FFTW_BACKWARD,FFTW_ESTIMATE)
  if(.not.c_associated(oracle%forward).or..not.c_associated(oracle%backward)) &
    error stop 'callback FFTW plan'
  oracle%ready=.true.
  do column=1,size(targets,2)
    density=conjg(source)*targets(:,column)
    call exx_local_apply(oracle,density,potential,status)
    if(status/=0)error stop 'callback FFTW apply'
    action(:,column)=-source*potential
  enddo
  call exx_local_destroy(oracle)
  status=0
 end subroutine
end program
