#include "config.h"
program test_rt_dg_hybrid_length_gauge_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_get_halting_mode,ieee_set_halting_mode,ieee_set_flag,&
    ieee_invalid,ieee_divide_by_zero
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
  use rt_dg_hybrid_length_gauge,only:propagate_rt_dg_hybrid_length_gauge
  use rt_dg_hybrid_sparse_exchange,only:s_rt_dg_sparse_exchange,build_rt_dg_sparse_exchange
  implicit none
  integer,parameter::n=3
  integer::comm,rank,nproc,ierr,nowned,i,j,k,a,b,step,iterations,permutation(n),claimed_rank
  complex(real64)::dense_s(n,n),dense_h(n,n),dense_z(3,n,n),initial(n),reference(n),phase(n),expected
  complex(real64)::original_s(n,n),original_h(n,n),original_z(3,n,n)
  complex(real64),allocatable::coeff(:),next(:),symmetry_next(:),related_next(:),invalid_coeff(:)
  type(s_dg_hybrid_sparse_metric)::metric,invalid_metric
  type(s_dg_hybrid_sparse_operators)::operators,invalid_operators
  type(s_rt_dg_sparse_exchange)::metric_exchange,operator_exchange
  real(real64)::norm_value,energy,polarization(3),previous(3),periods(3),defect,initial_norm,initial_energy,wrapped,&
    odd_amplitude
  real(real64)::field(3),related_field(3),representation_sign(n),cartesian_rotation(3,3)
  integer(int64)::workspace,global_workspace,fingerprint,reference_fingerprint
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
  allocate(invalid_coeff(0))
  call propagate_rt_dg_hybrid_length_gauge(comm,n,invalid_metric,invalid_operators,invalid_coeff,&
    [0d0,0d0,0d0],0.01d0,1d-12,24,[0d0,0d0,0d0],[10d0,10d0,10d0],next,&
    norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
  call require(.not.ok,'partially allocated length-gauge state was not rejected cleanly')
  invalid_metric=metric;invalid_operators=operators
  if(rank==0.and.size(invalid_metric%owned_row_ids)>0)then
    invalid_metric%owned_row_ids(1)=int(n+1,int64)
    invalid_operators%owned_row_ids(1)=int(n+1,int64)
  endif
  call propagate_rt_dg_hybrid_length_gauge(comm,n,invalid_metric,invalid_operators,&
    [(cmplx(0d0,0d0,real64),i=1,size(invalid_metric%owned_row_ids))],[0d0,0d0,0d0],&
    0.01d0,1d-12,24,[0d0,0d0,0d0],[10d0,10d0,10d0],next,norm_value,energy,&
    polarization,iterations,workspace,fingerprint,ok,message)
  call require(.not.ok,'out-of-range certified row ID was not rejected cleanly')
  nowned=size(metric%owned_row_ids);allocate(coeff(nowned))
  initial=[(0.7d0,0.1d0),(-0.2d0,0.3d0),(0.4d0,-0.1d0)]
  initial=initial/sqrt(real(dot_product(initial,matmul(dense_s,initial))))
  do i=1,nowned;coeff(i)=initial(int(metric%owned_row_ids(i)));enddo
  field=(0d0,0d0);previous=0d0;periods=10d0
  claimed_rank=n;if(rank==0)claimed_rank=n-1
  call propagate_rt_dg_hybrid_length_gauge(comm,claimed_rank,metric,operators,coeff,field,0.01d0,1d-12,24,previous,periods,&
    next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
  call require(.not.ok,'rank-disagreeing certified RT extent was accepted')
  call propagate_rt_dg_hybrid_length_gauge(comm,n,metric,operators,coeff,field,0.01d0,1d-12,24,previous,periods,&
    next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
  call require(ok,trim(message));reference_fingerprint=fingerprint;initial_norm=norm_value;initial_energy=energy
  call dense_step(initial,field,0.01d0,reference)
  defect=owned_defect(next,reference);call require(defect<2d-10,'one-step generalized propagation differs from dense oracle')
  call MPI_Allreduce(workspace,global_workspace,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
  call require(abs(norm_value-1d0)<2d-11.and.global_workspace>0_int64.and.iterations>0,&
    'zero-field propagation receipts are invalid')

  coeff=next
  do step=2,20
    previous=polarization
    call propagate_rt_dg_hybrid_length_gauge(comm,n,metric,operators,coeff,field,0.01d0,1d-12,24,previous,periods,&
      next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
    call require(ok,trim(message));coeff=next
  enddo
  call require(abs(norm_value-initial_norm)<3d-10.and.abs(energy-initial_energy)<3d-10,&
    'zero-field metric norm or energy drifted')

  ! Reject the third packet while retaining stored cross-packet S/H/Z edges.
  ! Propagation and every observable must use the identical active subspace.
  call distribute_system(original_s,original_h,original_z,metric,operators)
  metric%active_rows=[.true.,.true.,.false.];metric%packet_ids=[1,1,2];metric%numerical_rank=2
  initial=[(0.7d0,0.1d0),(-0.2d0,0.3d0),(0d0,0d0)]
  dense_s=original_s;dense_h=original_h;dense_z=original_z
  dense_s(3,:)=(0d0,0d0);dense_s(:,3)=(0d0,0d0)
  dense_s(3,3)=(1d0,0d0)
  dense_h(3,:)=(0d0,0d0);dense_h(:,3)=(0d0,0d0)
  dense_z(:,3,:)=(0d0,0d0);dense_z(:,:,3)=(0d0,0d0)
  initial=initial/sqrt(real(dot_product(initial,matmul(dense_s,initial))))
  do i=1,nowned;coeff(i)=initial(int(metric%owned_row_ids(i)));enddo
  do i=1,nowned;if(metric%owned_row_ids(i)==3_int64)coeff(i)=(5d-13,0d0);enddo
  field=[0.1d0,-0.03d0,0.02d0];previous=0d0
  call propagate_rt_dg_hybrid_length_gauge(comm,n,metric,operators,coeff,field,0.02d0,1d-12,24,previous,periods,&
    next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
  call require(ok,trim(message));call dense_step(initial,field,0.02d0,reference)
  defect=owned_defect(next,reference)
  call require(defect<3d-10,'active-subspace length-gauge propagation differs from dense oracle')
  defect=0d0
  do i=1,nowned;if(metric%owned_row_ids(i)==3_int64)defect=max(defect,abs(next(i)));enddo
  call MPI_Allreduce(MPI_IN_PLACE,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(defect==0d0,'inactive coefficient was not canonicalized')

  dense_s=original_s;dense_h=original_h;dense_z=original_z
  call distribute_system(dense_s,dense_h,dense_z,metric,operators)
  initial=[(0.7d0,0.1d0),(-0.2d0,0.3d0),(0.4d0,-0.1d0)]
  initial=initial/sqrt(real(dot_product(initial,matmul(dense_s,initial))))

  do i=1,nowned;coeff(i)=initial(int(metric%owned_row_ids(i)));enddo
  field=[(0.15d0,0d0),(-0.04d0,0d0),(0.02d0,0d0)];previous=[10.1d0,0d0,0d0]
  call propagate_rt_dg_hybrid_length_gauge(comm,n,metric,operators,coeff,field,0.03d0,1d-12,24,previous,periods,&
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
  call propagate_rt_dg_hybrid_length_gauge(comm,n,metric,operators,coeff,field,0.03d0,1d-12,24,previous,periods,&
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
  call propagate_rt_dg_hybrid_length_gauge(comm,n,metric,operators,coeff,field,0.03d0,1d-12,24,previous,periods,&
    next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
  call require(ok,trim(message));defect=0d0
  do i=1,nowned
    expected=conjg(phase(int(metric%owned_row_ids(i))))*reference(int(metric%owned_row_ids(i)))
    defect=max(defect,abs(next(i)-expected))
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(defect<3d-10,'length-gauge propagation is not complex-gauge covariant')
  call build_rt_dg_sparse_exchange(comm,n,metric%fingerprint,metric%owned_row_ids,metric%column_ids,&
    metric_exchange,ok,message);call require(ok,trim(message))
  call build_rt_dg_sparse_exchange(comm,n,operators%fingerprint,operators%owned_row_ids,operators%column_ids,&
    operator_exchange,ok,message);call require(ok,trim(message))
  do i=1,nowned
    do k=operators%row_offsets(i),operators%row_offsets(i+1)-1
      if(operators%column_ids(k)==int(operators%owned_row_ids(i)))&
        operators%hamiltonian_values(k)=operators%hamiltonian_values(k)+(0.05d0,0d0)
    enddo
  enddo
  call propagate_rt_dg_hybrid_length_gauge(comm,n,metric,operators,coeff,field,0.01d0,1d-12,24,previous,periods,&
    next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message,&
    metric_exchange,operator_exchange)
  call require(ok,trim(message));call require(operator_exchange%catalog_fingerprint==operators%fingerprint,&
    'operator value update invalidated the structure-keyed exchange schedule')

  ! A nonidentity reflection has D=diag(1,-1,1) in the certified basis and
  ! Q=diag(-1,1,1) in Cartesian space.  Z_x is odd while Z_y and Z_z are even,
  ! so D^H Z_a D = sum_b Q_ab Z_b.
  dense_s=(0d0,0d0);dense_h=(0d0,0d0);dense_z=(0d0,0d0)
  do i=1,n;dense_s(i,i)=(1d0,0d0);enddo
  dense_h(1,1)=(-0.4d0,0d0);dense_h(2,2)=(0.2d0,0d0);dense_h(3,3)=(0.6d0,0d0)
  dense_z(1,1,2)=(0.3d0,-0.04d0);dense_z(1,2,1)=conjg(dense_z(1,1,2))
  dense_z(1,2,3)=(-0.11d0,0.02d0);dense_z(1,3,2)=conjg(dense_z(1,2,3))
  dense_z(2,1,1)=(0.1d0,0d0);dense_z(2,2,2)=(0.25d0,0d0);dense_z(2,3,3)=(-0.2d0,0d0)
  dense_z(2,1,3)=(0.06d0,0.01d0);dense_z(2,3,1)=conjg(dense_z(2,1,3))
  dense_z(3,1,1)=(-0.05d0,0d0);dense_z(3,2,2)=(0.08d0,0d0);dense_z(3,3,3)=(0.12d0,0d0)
  representation_sign=[1d0,-1d0,1d0];cartesian_rotation=0d0
  cartesian_rotation(1,1)=-1d0;cartesian_rotation(2,2)=1d0;cartesian_rotation(3,3)=1d0
  defect=0d0
  do a=1,3;do i=1,n;do j=1,n
    expected=(0d0,0d0)
    do b=1,3;expected=expected+cartesian_rotation(a,b)*dense_z(b,i,j);enddo
    defect=max(defect,abs(representation_sign(i)*representation_sign(j)*dense_z(a,i,j)-expected))
  enddo;enddo;enddo
  call require(defect<1d-14.and.maxval(abs(dense_z(1,:,:)))>0.1d0,&
    'nonidentity vector-covariance fixture is invalid or vacuous')
  call distribute_system(dense_s,dense_h,dense_z,metric,operators)

  ! A field in the Q-invariant plane must not drive an initially invariant
  ! coefficient vector out of the D-invariant sector.
  initial=[(0.8d0,0d0),(0d0,0d0),(0.6d0,0d0)]
  do i=1,nowned;coeff(i)=initial(int(metric%owned_row_ids(i)));enddo
  field=[0d0,0.12d0,-0.04d0];previous=0d0
  call propagate_rt_dg_hybrid_length_gauge(comm,n,metric,operators,coeff,field,0.025d0,1d-12,24,previous,periods,&
    next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
  call require(ok,trim(message));defect=0d0
  do i=1,nowned
    j=int(metric%owned_row_ids(i));defect=max(defect,abs((representation_sign(j)-1d0)*next(i)))
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(defect<3d-10,'subgroup-preserving field escaped the invariant certified sector')

  ! A generic field may lower the equilibrium symmetry, but the solutions at
  ! E and Q^T E must remain D-related when their initial states are D-related.
  initial=[(0.8d0,0d0),(0d0,0d0),(0.6d0,0d0)]
  initial=initial/sqrt(real(dot_product(initial,initial)))
  do i=1,nowned;coeff(i)=initial(int(metric%owned_row_ids(i)));enddo
  field=[0.17d0,0.08d0,-0.03d0];related_field=matmul(transpose(cartesian_rotation),field);previous=0d0
  call propagate_rt_dg_hybrid_length_gauge(comm,n,metric,operators,coeff,field,0.025d0,1d-12,24,previous,periods,&
    next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
  call require(ok,trim(message));allocate(symmetry_next,source=next);odd_amplitude=0d0
  do i=1,nowned
    j=int(metric%owned_row_ids(i))
    if(representation_sign(j)<0d0)odd_amplitude=max(odd_amplitude,abs(symmetry_next(i)))
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,odd_amplitude,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(odd_amplitude>1d-6,'symmetry-lowering field did not generate an odd-sector amplitude')
  do i=1,nowned
    j=int(metric%owned_row_ids(i));coeff(i)=representation_sign(j)*initial(j)
  enddo
  call propagate_rt_dg_hybrid_length_gauge(comm,n,metric,operators,coeff,related_field,0.025d0,1d-12,24,previous,periods,&
    next,norm_value,energy,polarization,iterations,workspace,fingerprint,ok,message)
  call require(ok,trim(message));allocate(related_next,source=next);defect=0d0
  do i=1,nowned
    j=int(metric%owned_row_ids(i));defect=max(defect,abs(related_next(i)-representation_sign(j)*symmetry_next(i)))
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  call require(defect<3d-10,'symmetry-related fields did not produce D-related certified coefficients')
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
