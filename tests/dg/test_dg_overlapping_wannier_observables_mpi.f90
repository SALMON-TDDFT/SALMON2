#include "config.h"
program test_dg_overlapping_wannier_observables_mpi
  use mpi
  use dg_overlapping_wannier_observables,only:assemble_dg_overlapping_wannier_observables
  implicit none
  integer::comm,rank,nproc,ierr,p,i,nlocal,index,owned
  integer(8),allocatable::ids(:)
  real(8),allocatable::weights(:),coordinates(:,:)
  complex(8),allocatable::values(:,:),gradients(:,:,:),position(:,:,:),momentum(:,:,:),&
    derivative(:,:,:),velocity(:,:,:),nonlocal_velocity(:,:,:),rotated_values(:,:),&
    rotated_gradients(:,:,:)
  complex(8)::metric(2,2),metric_inverse(2,2),hlocal(2,2),hnonlocal(2,2),gauge(2,2),basis_map(2,2),&
    reference_position(3,2,2),reference_velocity(3,2,2),direct(2,2)
  real(8)::origin(3),cell(3),x
  logical::ok
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,4)])
  allocate(ids(nlocal),weights(nlocal),coordinates(3,nlocal),values(2,nlocal),gradients(3,2,nlocal))
  index=0
  do p=1,4
    if(mod(p-1,nproc)/=rank)cycle
    index=index+1;ids(index)=p;weights(index)=1d0;coordinates(:,index)=[0.2d0*p,0.1d0*p,0.05d0*p]
    values(:,index)=0.5d0*[cmplx(1d0,0d0,8),cmplx((-1d0)**p,0d0,8)]
    gradients(:,1,index)=(0d0,0d0);gradients(:,2,index)=(0d0,0d0)
    gradients(1,1,index)=0.4d0*values(2,index)
    gradients(1,2,index)=-0.4d0*values(1,index)
  enddo
  metric=(0d0,0d0);metric(1,1)=1d0;metric(2,2)=1d0;metric_inverse=metric
  hlocal=reshape([cmplx(1d0,0d0,8),cmplx(0.2d0,0d0,8),&
    cmplx(0.2d0,0d0,8),cmplx(2d0,0d0,8)],[2,2])
  hnonlocal=reshape([cmplx(0.1d0,0d0,8),cmplx(0d0,-0.05d0,8),&
    cmplx(0d0,0.05d0,8),cmplx(0.3d0,0d0,8)],[2,2])
  basis_map=reshape([cmplx(1d0,0d0,8),cmplx(0d0,0d0,8),&
    cmplx(0.25d0,0d0,8),cmplx(1d0,0d0,8)],[2,2])
  values=matmul(basis_map,values)
  do p=1,nlocal;do i=1,3
    gradients(i,:,p)=matmul(basis_map,gradients(i,:,p))
  enddo;enddo
  metric=matmul(conjg(basis_map),matmul(metric,transpose(basis_map)))
  metric_inverse(1,1)=metric(2,2);metric_inverse(2,2)=metric(1,1)
  metric_inverse(1,2)=-metric(1,2);metric_inverse(2,1)=-metric(2,1)
  metric_inverse=metric_inverse/(metric(1,1)*metric(2,2)-metric(1,2)*metric(2,1))
  hlocal=matmul(conjg(basis_map),matmul(hlocal,transpose(basis_map)))
  hnonlocal=matmul(conjg(basis_map),matmul(hnonlocal,transpose(basis_map)))
  origin=[0d0,0d0,0d0];cell=[2d0,2d0,2d0]
  call assemble_dg_overlapping_wannier_observables(comm,2,ids,weights,coordinates,origin,cell,&
    'cell_wrapped',values,gradients,metric,metric_inverse,hlocal,hnonlocal,4_8,1d-12,&
    position,derivative,momentum,velocity,nonlocal_velocity,owned,ok,message)
  call require(ok,trim(message));call require(owned==4,'observable core ownership')
  do i=1,3
    call require(maxval(abs(position(i,:,:)-conjg(transpose(position(i,:,:)))))<1d-13,&
      'Hermitian position')
    call require(maxval(abs(derivative(i,:,:)+conjg(transpose(derivative(i,:,:)))))<1d-13,&
      'anti-Hermitian derivative')
    call require(maxval(abs(momentum(i,:,:)-conjg(transpose(momentum(i,:,:)))))<1d-13,&
      'Hermitian canonical momentum')
    call require(maxval(abs(velocity(i,:,:)-conjg(transpose(velocity(i,:,:)))))<1d-13,&
      'Hermitian physical velocity')
    direct=cmplx(0d0,1d0,8)*(matmul(hlocal+hnonlocal,matmul(metric_inverse,position(i,:,:)))-&
      matmul(position(i,:,:),matmul(metric_inverse,hlocal+hnonlocal)))
    call require(maxval(abs(velocity(i,:,:)-direct))<1d-13,'complete commutator velocity')
    direct=cmplx(0d0,1d0,8)*(matmul(hnonlocal,matmul(metric_inverse,position(i,:,:)))-&
      matmul(position(i,:,:),matmul(metric_inverse,hnonlocal)))
    call require(maxval(abs(nonlocal_velocity(i,:,:)-direct))<1d-13,&
      'explicit nonlocal commutator contribution')
  enddo
  call require(abs(momentum(1,1,2))>1d-8,'canonical momentum diagnostic')
  reference_position=position;reference_velocity=velocity
  gauge=reshape([cmplx(sqrt(0.5d0),0d0,8),cmplx(-sqrt(0.5d0),0d0,8),&
    cmplx(sqrt(0.5d0),0d0,8),cmplx(sqrt(0.5d0),0d0,8)],[2,2])
  rotated_values=matmul(gauge,values);allocate(rotated_gradients(3,2,nlocal))
  do p=1,nlocal;do i=1,3
    rotated_gradients(i,:,p)=matmul(gauge,gradients(i,:,p))
  enddo;enddo
  metric=matmul(conjg(gauge),matmul(metric,transpose(gauge)))
  metric_inverse=matmul(conjg(gauge),matmul(metric_inverse,transpose(gauge)))
  hlocal=matmul(conjg(gauge),matmul(hlocal,transpose(gauge)))
  hnonlocal=matmul(conjg(gauge),matmul(hnonlocal,transpose(gauge)))
  call assemble_dg_overlapping_wannier_observables(comm,2,ids,weights,coordinates,origin,cell,&
    'cell_wrapped',rotated_values,rotated_gradients,metric,metric_inverse,hlocal,hnonlocal,4_8,1d-12,&
    position,derivative,momentum,velocity,nonlocal_velocity,owned,ok,message)
  call require(ok,trim(message))
  do i=1,3
    call require(maxval(abs(position(i,:,:)-matmul(conjg(gauge),&
      matmul(reference_position(i,:,:),transpose(gauge)))))<1d-12,'position gauge covariance')
    call require(maxval(abs(velocity(i,:,:)-matmul(conjg(gauge),&
      matmul(reference_velocity(i,:,:),transpose(gauge)))))<1d-12,'velocity gauge covariance')
  enddo
  if(nproc>1)then
    call assemble_dg_overlapping_wannier_observables(comm,merge(3,2,rank==0),ids,weights,coordinates,&
      origin,cell,'cell_wrapped',rotated_values,rotated_gradients,metric,metric_inverse,hlocal,&
      hnonlocal,4_8,1d-12,position,derivative,momentum,velocity,nonlocal_velocity,owned,ok,message)
    call require(.not.ok,'rank-inconsistent observable contract rejection')
  endif
  if(rank==0)then
    write(*,'(a,i0,a,4(es24.16,1x))')'OBSERVABLES ranks=',nproc,' values=',&
      real(reference_position(1,1,2)),aimag(reference_position(1,1,2)),&
      real(reference_velocity(1,1,2)),aimag(reference_velocity(1,1,2))
    write(*,'(a,i0,a)')'PASS overlapping-Wannier observables on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::lf,gf
    lf=merge(0,1,condition);call MPI_Allreduce(lf,gf,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(gf/=0)error stop label
  end subroutine
end program
