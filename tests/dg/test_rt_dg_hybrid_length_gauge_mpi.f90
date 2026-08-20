#include "config.h"
program test_rt_dg_hybrid_length_gauge_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_get_halting_mode,ieee_set_halting_mode,ieee_set_flag,&
    ieee_invalid,ieee_divide_by_zero
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
  use rt_dg_hybrid_length_gauge,only:propagate_rt_dg_hybrid_length_gauge
  implicit none
  integer,parameter::n=3
  integer::comm,rank,nproc,ierr,nowned,i,j,k,step,iterations,permutation(n)
  complex(real64)::dense_s(n,n),dense_h(n,n),dense_z(3,n,n),initial(n),reference(n),phase(n),expected
  complex(real64)::original_s(n,n),original_h(n,n),original_z(3,n,n)
  complex(real64),allocatable::coeff(:),next(:)
  type(s_dg_hybrid_sparse_metric)::metric
  type(s_dg_hybrid_sparse_operators)::operators
  real(real64)::norm_value,energy,polarization(3),previous(3),periods(3),defect,initial_norm,initial_energy,wrapped
  real(real64)::field(3)
  integer(int64)::workspace,fingerprint,reference_fingerprint
  logical::ok
  character(256)::message
  external::zgesv
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  dense_s=reshape([(1.2d0,0d0),(0.08d0,-0.03d0),(0d0,0d0),&
    (0.08d0,0.03d0),(1.05d0,0d0),(0.04d0,0.01d0),&
    (0d0,0d0),(0.04d0,-0.01d0),(0.95d0,0d0)],[n,n])
  dense_h=reshape([(-0.5d0,0d0),(0.12d0,-0.02d0),(0.03d0,0.01d0),&
    (0.12d0,0.02d0),(0.2d0,0d0),(-0.08d0,0.04d0),&
    (0.03d0,-0.01d0),(-0.08d0,-0.04d0),(0.7d0,0d0)],[n,n])
  dense_z=(0d0,0d0)
  dense_z(1,:,:)=reshape([(0.1d0,0d0),(0.2d0,-0.1d0),(0d0,0d0),&
    (0.2d0,0.1d0),(0.5d0,0d0),(0.15d0,0.02d0),&
    (0d0,0d0),(0.15d0,-0.02d0),(0.9d0,0d0)],[n,n])
  dense_z(2,:,:)=reshape([(0d0,0d0),(0.1d0,0.04d0),(0.05d0,0d0),&
    (0.1d0,-0.04d0),(0.2d0,0d0),(-0.03d0,0.02d0),&
    (0.05d0,0d0),(-0.03d0,-0.02d0),(0.4d0,0d0)],[n,n])
  dense_z(3,:,:)=reshape([(0.3d0,0d0),(0d0,0d0),(0.06d0,-0.02d0),&
    (0d0,0d0),(0.1d0,0d0),(0.07d0,0.01d0),&
    (0.06d0,0.02d0),(0.07d0,-0.01d0),(0.2d0,0d0)],[n,n])
  original_s=dense_s;original_h=dense_h;original_z=dense_z
  call distribute_system(dense_s,dense_h,dense_z,metric,operators)
  nowned=size(metric%owned_row_ids);allocate(coeff(nowned))
  initial=[(0.7d0,0.1d0),(-0.2d0,0.3d0),(0.4d0,-0.1d0)]
  initial=initial/sqrt(real(dot_product(initial,matmul(dense_s,initial))))
  do i=1,nowned;coeff(i)=initial(int(metric%owned_row_ids(i)));enddo
  field=(0d0,0d0);previous=0d0;periods=10d0
  call propagate_rt_dg_hybrid_length_gauge(comm,metric,operators,coeff,field,0.01d0,1d-12,24,previous,periods,&
    next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
  call require(ok,trim(message));reference_fingerprint=fingerprint;initial_norm=norm_value;initial_energy=energy
  call dense_step(initial,field,0.01d0,reference)
  defect=owned_defect(next,reference);call require(defect<2d-10,'one-step generalized propagation differs from dense oracle')
  call require(abs(norm_value-1d0)<2d-11.and.workspace>0_int64.and.iterations>0,'zero-field propagation receipts are invalid')

  coeff=next
  do step=2,20
    previous=polarization
    call propagate_rt_dg_hybrid_length_gauge(comm,metric,operators,coeff,field,0.01d0,1d-12,24,previous,periods,&
      next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
    call require(ok,trim(message));coeff=next
  enddo
  call require(abs(norm_value-initial_norm)<3d-10.and.abs(energy-initial_energy)<3d-10,&
    'zero-field metric norm or energy drifted')

  do i=1,nowned;coeff(i)=initial(int(metric%owned_row_ids(i)));enddo
  field=[(0.15d0,0d0),(-0.04d0,0d0),(0.02d0,0d0)];previous=[10.1d0,0d0,0d0]
  call propagate_rt_dg_hybrid_length_gauge(comm,metric,operators,coeff,field,0.03d0,1d-12,24,previous,periods,&
    next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
  call require(ok,trim(message));call dense_step(initial,field,0.03d0,reference)
  defect=owned_defect(next,reference);call require(defect<3d-10,'field-driven generalized propagation differs from dense oracle')
  wrapped=real(dot_product(reference,matmul(dense_z(1,:,:),reference)))
  wrapped=wrapped+anint((previous(1)-wrapped)/periods(1))*periods(1)
  call require(abs(polarization(1)-wrapped)<3d-10,'propagation and observable did not use the identical Z operator')
  call require(abs(polarization(1)-previous(1))<5d0,'polarization branch was not continued near previous value')
  call require(abs(polarization(1)-real(dot_product(reference,matmul(dense_z(1,:,:),reference))))>1d0,&
    'branch-continuity fixture did not cross a polarization quantum')

  permutation=[2,3,1]
  call permute_system(permutation,original_s,original_h,original_z,dense_s,dense_h,dense_z)
  call distribute_system(dense_s,dense_h,dense_z,metric,operators)
  do i=1,nowned;coeff(i)=initial(permutation(int(metric%owned_row_ids(i))));enddo
  previous=0d0
  call propagate_rt_dg_hybrid_length_gauge(comm,metric,operators,coeff,field,0.03d0,1d-12,24,previous,periods,&
    next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
  call require(ok,trim(message));defect=0d0
  do i=1,nowned
    defect=max(defect,abs(next(i)-reference(permutation(int(metric%owned_row_ids(i))))))
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(defect<3d-10,'length-gauge propagation is not permutation covariant')

  phase=[(1d0,0d0),exp(cmplx(0d0,0.31d0,real64)),exp(cmplx(0d0,-0.27d0,real64))]
  dense_s=original_s;dense_h=original_h;dense_z=original_z
  call gauge_transform_system(phase,dense_s,dense_h,dense_z)
  call distribute_system(dense_s,dense_h,dense_z,metric,operators)
  do i=1,nowned;coeff(i)=conjg(phase(int(metric%owned_row_ids(i))))*initial(int(metric%owned_row_ids(i)));enddo
  previous=0d0
  call propagate_rt_dg_hybrid_length_gauge(comm,metric,operators,coeff,field,0.03d0,1d-12,24,previous,periods,&
    next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
  call require(ok,trim(message));defect=0d0
  do i=1,nowned
    expected=conjg(phase(int(metric%owned_row_ids(i))))*reference(int(metric%owned_row_ids(i)))
    defect=max(defect,abs(next(i)-expected))
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(defect<3d-10,'length-gauge propagation is not complex-gauge covariant')
  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_LENGTH_GAUGE ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid length gauge on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine distribute_system(s,h,z,distributed_metric,distributed_operators)
    complex(real64),intent(in)::s(:,:),h(:,:),z(:,:,:)
    type(s_dg_hybrid_sparse_metric),intent(out)::distributed_metric
    type(s_dg_hybrid_sparse_operators),intent(out)::distributed_operators
    integer::row,column,position,edge
    nowned=count([(mod(row-1,nproc)==rank,row=1,n)])
    allocate(distributed_metric%owned_row_ids(nowned),distributed_metric%row_offsets(nowned+1),&
      distributed_metric%column_ids(nowned*n),distributed_metric%values(nowned*n),&
      distributed_metric%active_rows(n),distributed_metric%packet_ids(n))
    allocate(distributed_operators%owned_row_ids(nowned),distributed_operators%row_offsets(nowned+1),&
      distributed_operators%column_ids(nowned*n),distributed_operators%metric_values(nowned*n),&
      distributed_operators%hamiltonian_values(nowned*n),distributed_operators%position_values(3,nowned*n))
    position=0;edge=0;distributed_metric%row_offsets(1)=1
    do row=n,1,-1
      if(mod(row-1,nproc)/=rank)cycle
      position=position+1;distributed_metric%owned_row_ids(position)=row
      do column=1,n
        edge=edge+1;distributed_metric%column_ids(edge)=column;distributed_metric%values(edge)=s(row,column)
        distributed_operators%column_ids(edge)=column;distributed_operators%metric_values(edge)=s(row,column)
        distributed_operators%hamiltonian_values(edge)=h(row,column);distributed_operators%position_values(:,edge)=z(:,row,column)
      enddo
      distributed_metric%row_offsets(position+1)=edge+1
    enddo
    distributed_metric%valid=.true.;distributed_metric%global_count=n;distributed_metric%numerical_rank=n
    distributed_metric%condition_estimate=2d0;distributed_metric%maximum_value=maxval(abs(s))
    distributed_metric%max_row_nnz=n;distributed_metric%fingerprint=8181_int64
    distributed_metric%active_rows=.true.;distributed_metric%packet_ids=1
    distributed_operators%valid=.true.;distributed_operators%global_count=n
    distributed_operators%owned_row_ids=distributed_metric%owned_row_ids
    distributed_operators%row_offsets=distributed_metric%row_offsets
    distributed_operators%metric_fingerprint=distributed_metric%fingerprint
    distributed_operators%fingerprint=7171_int64
  end subroutine distribute_system
  subroutine dense_step(input,electric,dt,output)
    complex(real64),intent(in)::input(:);real(real64),intent(in)::electric(:),dt
    complex(real64),intent(out)::output(:)
    complex(real64)::term(n),rhs(n,1),matrix(n,n),a(n,n)
    integer::order,pivots(n),info
    logical::halt_invalid,halt_zero
    a=dense_h+electric(1)*dense_z(1,:,:)+electric(2)*dense_z(2,:,:)+electric(3)*dense_z(3,:,:)
    output=input;term=input
    do order=1,24
      rhs(:,1)=matmul(a,term);matrix=dense_s
      call ieee_get_halting_mode(ieee_invalid,halt_invalid);call ieee_get_halting_mode(ieee_divide_by_zero,halt_zero)
      call ieee_set_halting_mode(ieee_invalid,.false.);call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
      call zgesv(n,1,matrix,n,pivots,rhs,n,info)
      call ieee_set_flag(ieee_invalid,.false.);call ieee_set_flag(ieee_divide_by_zero,.false.)
      call ieee_set_halting_mode(ieee_invalid,halt_invalid);call ieee_set_halting_mode(ieee_divide_by_zero,halt_zero)
      if(info/=0)error stop 'dense length-gauge solve failed'
      term=cmplx(0d0,-dt/real(order,real64),real64)*rhs(:,1);output=output+term
    enddo
  end subroutine dense_step
  real(real64) function owned_defect(actual,global_reference)
    complex(real64),intent(in)::actual(:),global_reference(:);integer::q
    owned_defect=0d0
    do q=1,nowned;owned_defect=max(owned_defect,abs(actual(q)-global_reference(int(metric%owned_row_ids(q)))));enddo
    call MPI_Allreduce(MPI_IN_PLACE,owned_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  end function owned_defect
  subroutine gauge_transform_system(q,s,h,z)
    complex(real64),intent(in)::q(:);complex(real64),intent(inout)::s(:,:),h(:,:),z(:,:,:)
    integer::a,b,d
    do a=1,n;do b=1,n
      s(a,b)=conjg(q(a))*q(b)*s(a,b);h(a,b)=conjg(q(a))*q(b)*h(a,b)
      do d=1,3;z(d,a,b)=conjg(q(a))*q(b)*z(d,a,b);enddo
    enddo;enddo
  end subroutine gauge_transform_system
  subroutine permute_system(order,s0,h0,z0,s,h,z)
    integer,intent(in)::order(:)
    complex(real64),intent(in)::s0(:,:),h0(:,:),z0(:,:,:)
    complex(real64),intent(out)::s(:,:),h(:,:),z(:,:,:)
    integer::a,b,d
    do a=1,n;do b=1,n
      s(a,b)=s0(order(a),order(b));h(a,b)=h0(order(a),order(b))
      do d=1,3;z(d,a,b)=z0(d,order(a),order(b));enddo
    enddo;enddo
  end subroutine permute_system
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label;integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)error stop label
  end subroutine require
end program test_rt_dg_hybrid_length_gauge_mpi
