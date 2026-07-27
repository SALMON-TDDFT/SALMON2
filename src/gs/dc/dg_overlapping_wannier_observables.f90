#include "config.h"
module dg_overlapping_wannier_observables
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::assemble_dg_overlapping_wannier_observables
contains
  subroutine assemble_dg_overlapping_wannier_observables(comm,nwann,core_ids,weights,coordinates,&
      origin,cell_length,position_convention,values,gradients,metric,metric_inverse,local_hamiltonian,&
      nonlocal_hamiltonian,expected_core_count,antihermitian_tolerance,position,derivative,&
      canonical_momentum,velocity,nonlocal_velocity,ownership_count,ok,message)
    integer,intent(in)::comm,nwann
    integer(int64),intent(in)::core_ids(:),expected_core_count
    real(real64),intent(in)::weights(:),coordinates(:,:),origin(3),cell_length(3)
    character(*),intent(in)::position_convention
    complex(real64),intent(in)::values(:,:),gradients(:,:,:),metric(:,:),metric_inverse(:,:),&
      local_hamiltonian(:,:),nonlocal_hamiltonian(:,:)
    real(real64),intent(in)::antihermitian_tolerance
    complex(real64),allocatable,intent(out)::position(:,:,:),derivative(:,:,:),&
      canonical_momentum(:,:,:),velocity(:,:,:),nonlocal_velocity(:,:,:)
    integer,intent(out)::ownership_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nproc,ierr,local_bad,global_bad,total_count,i,j,p,axis,matrix_count,nwann_min,nwann_max
    integer,allocatable::counts(:),displs(:)
    integer(int64),allocatable::all_ids(:)
    integer(int64)::matrix_count64,expected_min,expected_max
    complex(real64),allocatable::local_position(:,:,:),local_derivative(:,:,:),identity(:,:),&
      total_hamiltonian(:,:),reference_matrix(:,:)
    real(real64)::wrapped_coordinate,tolerance_min,tolerance_max,position_defect,position_scale,&
      reference_geometry(6)
    character(len(position_convention))::reference_convention
    logical::shape_ok,finite_payload
    ok=.false.;message='';ownership_count=0;local_bad=0
    call MPI_Allreduce(nwann,nwann_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(nwann,nwann_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_core_count,expected_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    call MPI_Allreduce(antihermitian_tolerance,tolerance_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(antihermitian_tolerance,tolerance_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(nwann_min/=nwann_max.or.expected_min/=expected_max.or.tolerance_min/=tolerance_max)then
      message='inconsistent observable assembly contract across ranks';return
    endif
    shape_ok=size(weights)==size(core_ids).and.size(coordinates,1)==3.and.&
      size(coordinates,2)==size(core_ids).and.size(values,1)==nwann.and.&
      size(values,2)==size(core_ids).and.size(gradients,1)==3.and.&
      size(gradients,2)==nwann.and.size(gradients,3)==size(core_ids).and.&
      all(shape(metric)==[nwann,nwann]).and.all(shape(metric_inverse)==[nwann,nwann]).and.&
      all(shape(local_hamiltonian)==[nwann,nwann]).and.all(shape(nonlocal_hamiltonian)==[nwann,nwann])
    if(nwann<=0.or.expected_core_count<=0_int64.or.antihermitian_tolerance<=0d0)local_bad=1
    if(.not.shape_ok.or.trim(position_convention)/='cell_wrapped')local_bad=1
    if(any(core_ids<=0_int64).or.any(weights<=0d0).or.any(cell_length<=0d0))local_bad=1
    finite_payload=all(ieee_is_finite(weights)).and.all(ieee_is_finite(coordinates)).and.&
      all(ieee_is_finite(origin)).and.all(ieee_is_finite(cell_length))
    if(shape_ok)then
      finite_payload=finite_payload.and.finite_complex_2d(values).and.finite_complex_3d(gradients).and.&
        finite_complex_2d(metric).and.finite_complex_2d(metric_inverse).and.&
        finite_complex_2d(local_hamiltonian).and.finite_complex_2d(nonlocal_hamiltonian)
    endif
    if(.not.finite_payload)local_bad=1
    if(nwann>0.and.int(nwann,int64)<=huge(1_int64)/int(nwann,int64))then
      matrix_count64=int(nwann,int64)*int(nwann,int64)
      if(matrix_count64>int(huge(matrix_count),int64))local_bad=1
    else
      matrix_count64=0_int64;local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid observable unique-core payload or position convention';return;endif
    reference_geometry=[origin,cell_length]
    call MPI_Bcast(reference_geometry,6,MPI_DOUBLE_PRECISION,0,comm,ierr)
    reference_convention=position_convention
    call MPI_Bcast(reference_convention,len(reference_convention),MPI_CHARACTER,0,comm,ierr)
    local_bad=merge(0,1,maxval(abs(reference_geometry-[origin,cell_length]))<=antihermitian_tolerance.and.&
      reference_convention==position_convention)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='inconsistent observable cell/origin convention across ranks';return;endif
    matrix_count=int(matrix_count64)
    allocate(reference_matrix(nwann,nwann))
    reference_matrix=metric
    call MPI_Bcast(reference_matrix,matrix_count,MPI_DOUBLE_COMPLEX,0,comm,ierr)
    if(maxval(abs(reference_matrix-metric))>antihermitian_tolerance)local_bad=1
    reference_matrix=metric_inverse
    call MPI_Bcast(reference_matrix,matrix_count,MPI_DOUBLE_COMPLEX,0,comm,ierr)
    if(maxval(abs(reference_matrix-metric_inverse))>antihermitian_tolerance)local_bad=1
    reference_matrix=local_hamiltonian
    call MPI_Bcast(reference_matrix,matrix_count,MPI_DOUBLE_COMPLEX,0,comm,ierr)
    if(maxval(abs(reference_matrix-local_hamiltonian))>antihermitian_tolerance)local_bad=1
    reference_matrix=nonlocal_hamiltonian
    call MPI_Bcast(reference_matrix,matrix_count,MPI_DOUBLE_COMPLEX,0,comm,ierr)
    if(maxval(abs(reference_matrix-nonlocal_hamiltonian))>antihermitian_tolerance)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='inconsistent observable metric or Hamiltonian across ranks';return;endif
    allocate(identity(nwann,nwann));identity=matmul(metric,metric_inverse)
    do i=1,nwann;identity(i,i)=identity(i,i)-1d0;enddo
    if(maxval(abs(identity))>antihermitian_tolerance.or.&
        maxval(abs(metric-conjg(transpose(metric))))>antihermitian_tolerance.or.&
        maxval(abs(local_hamiltonian-conjg(transpose(local_hamiltonian))))>antihermitian_tolerance.or.&
        maxval(abs(nonlocal_hamiltonian-conjg(transpose(nonlocal_hamiltonian))))>antihermitian_tolerance)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid metric inverse or non-Hermitian Hamiltonian';return;endif

    call MPI_Comm_size(comm,nproc,ierr);allocate(counts(nproc),displs(nproc))
    call MPI_Allgather(size(core_ids),1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    total_count=0;displs(1)=0
    do i=1,nproc
      if(counts(i)<0.or.total_count>huge(total_count)-counts(i))then;local_bad=1;exit;endif
      if(i>1)displs(i)=total_count
      total_count=total_count+counts(i)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.int(total_count,int64)/=expected_core_count)then
      message='missing or extra observable core owner';return
    endif
    allocate(all_ids(total_count))
    call MPI_Allgatherv(core_ids,size(core_ids),MPI_INTEGER8,all_ids,counts,displs,MPI_INTEGER8,comm,ierr)
    call sort_ids(all_ids)
    do i=1,total_count
      if(all_ids(i)/=int(i,int64))then;message='duplicate or missing observable core owner';return;endif
    enddo

    allocate(local_position(3,nwann,nwann),local_derivative(3,nwann,nwann))
    local_position=(0d0,0d0);local_derivative=(0d0,0d0)
    do p=1,size(core_ids);do j=1,nwann;do i=1,nwann
      do axis=1,3
        wrapped_coordinate=origin(axis)+modulo(coordinates(axis,p)-origin(axis),cell_length(axis))
        local_position(axis,i,j)=local_position(axis,i,j)+weights(p)*wrapped_coordinate*&
          conjg(values(i,p))*values(j,p)
        local_derivative(axis,i,j)=local_derivative(axis,i,j)+weights(p)*conjg(values(i,p))*&
          gradients(axis,j,p)
      enddo
    enddo;enddo;enddo
    allocate(position(3,nwann,nwann),derivative(3,nwann,nwann),canonical_momentum(3,nwann,nwann))
    allocate(velocity(3,nwann,nwann),nonlocal_velocity(3,nwann,nwann))
    do axis=1,3
      call MPI_Allreduce(local_position(axis,:,:),position(axis,:,:),matrix_count,&
        MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      call MPI_Allreduce(local_derivative(axis,:,:),derivative(axis,:,:),matrix_count,&
        MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      position_scale=max(1d0,maxval(abs(position(axis,:,:))))
      position_defect=maxval(abs(position(axis,:,:)-conjg(transpose(position(axis,:,:)))))
      if(position_defect>antihermitian_tolerance*position_scale)local_bad=1
      position(axis,:,:)=0.5d0*(position(axis,:,:)+conjg(transpose(position(axis,:,:))))
      if(maxval(abs(derivative(axis,:,:)+conjg(transpose(derivative(axis,:,:)))))>&
          antihermitian_tolerance)then
        local_bad=1
      endif
      canonical_momentum(axis,:,:)=-cmplx(0d0,1d0,real64)*derivative(axis,:,:)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      message='position is not Hermitian or canonical derivative is not anti-Hermitian';return
    endif
    allocate(total_hamiltonian,source=local_hamiltonian+nonlocal_hamiltonian)
    do axis=1,3
      velocity(axis,:,:)=cmplx(0d0,1d0,real64)*(matmul(total_hamiltonian,&
        matmul(metric_inverse,position(axis,:,:)))-matmul(position(axis,:,:),&
        matmul(metric_inverse,total_hamiltonian)))
      nonlocal_velocity(axis,:,:)=cmplx(0d0,1d0,real64)*(matmul(nonlocal_hamiltonian,&
        matmul(metric_inverse,position(axis,:,:)))-matmul(position(axis,:,:),&
        matmul(metric_inverse,nonlocal_hamiltonian)))
      if(maxval(abs(velocity(axis,:,:)-conjg(transpose(velocity(axis,:,:)))))>&
          antihermitian_tolerance.or.maxval(abs(nonlocal_velocity(axis,:,:)-&
          conjg(transpose(nonlocal_velocity(axis,:,:)))))>antihermitian_tolerance)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='physical or nonlocal velocity is not Hermitian';return;endif
    ownership_count=total_count;ok=.true.
#else
    ok=.false.;message='overlapping-Wannier observables require MPI';ownership_count=0
#endif
  end subroutine
  logical function finite_complex_2d(values)
    complex(real64),intent(in)::values(:,:)
    finite_complex_2d=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function
  logical function finite_complex_3d(values)
    complex(real64),intent(in)::values(:,:,:)
    finite_complex_3d=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function
  subroutine sort_ids(ids)
    integer(int64),intent(inout)::ids(:)
    integer::i,j
    integer(int64)::key
    do i=2,size(ids)
      key=ids(i);j=i-1
      do while(j>=1)
        if(ids(j)<=key)exit
        ids(j+1)=ids(j);j=j-1
      enddo
      ids(j+1)=key
    enddo
  end subroutine
end module dg_overlapping_wannier_observables
