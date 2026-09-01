#include "config.h"
program test_dg_hybrid_certified_rt_basis_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_quiet_nan,ieee_value
  use dg_hybrid_low_energy_symmetry,only:evaluate_dg_hybrid_low_energy_symmetry
  use dg_hybrid_certified_rt_basis,only:s_dg_hybrid_certified_rt_basis,&
    build_dg_hybrid_certified_rt_basis,validate_dg_hybrid_certified_rt_basis
  implicit none
  integer,parameter::n=4,nstress=64,ncert=2,nop=3,nscalar=2,nvector=1,ntensor=1
  integer::comm,rank,nproc,ierr,nowned,nowned_stress,i,a,b,position,localizer_mode,localizer_calls
  integer::observed_certified_rank,observed_global_count
  integer(int64)::bits_before,bits_after
  integer(int64),allocatable::row_ids(:)
  integer(int64),allocatable::stress_row_ids(:)
  integer(int64),allocatable::callback_expected_row_ids(:)
  complex(real64),allocatable::metric_rows(:,:),bad_metric_rows(:,:),c_cert_rows(:,:),&
    scaled_c_cert_rows(:,:),full_c_result(:,:),full_b_result(:,:),stress_metric(:,:),&
    stress_metric_rows(:,:),stress_c_rows(:,:),callback_expected_c(:,:)
  complex(real64)::metric(n,n),metric_variant(n,n),full_coefficients(n,n),&
    construction_representation(n,n,nop)
  complex(real64)::certified_representation(ncert,ncert,nop)
  complex(real64)::bad_representation(ncert,ncert,nop)
  complex(real64)::scalar_operators(ncert,ncert,nscalar),vector_operators(ncert,ncert,3,nvector)
  complex(real64)::bad_scalar_operators(ncert,ncert,nscalar)
  complex(real64)::bad_vector_operators(ncert,ncert,3,nvector)
  complex(real64)::tensor_operators(ncert,ncert,3,3,ntensor)
  complex(real64)::bad_tensor_operators(ncert,ncert,3,3,ntensor)
  complex(real64)::zero_scalar_operators(ncert,ncert,nscalar)
  complex(real64)::zero_vector_operators(ncert,ncert,3,nvector)
  complex(real64)::zero_tensor_operators(ncert,ncert,3,3,ntensor)
  complex(real64)::expected_u(ncert,ncert),expected_h(ncert,ncert),expected_a(ncert,1)
  complex(real64)::expected_scalar(ncert,ncert,nscalar),expected_vector(ncert,ncert,3,nvector)
  complex(real64)::expected_tensor(ncert,ncert,3,3,ntensor),identity(ncert,ncert)
  complex(real64)::projector_before(n,n),projector_after(n,n)
  real(real64)::rotations(3,3,nop),bad_rotations(3,3,nop),full_eigenvalues(n),&
    certified_eigenvalues(ncert),bad_eigenvalues(ncert)
  real(real64)::occupied_defect,target_defect,energy_defect,phase,nan_value_main
  logical::ok,observed_require_unconstrained,verify_callback_payload,observed_callback_payload
  character(256)::message
  type(s_dg_hybrid_certified_rt_basis)::result,tampered

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nowned=0
  do i=rank+1,n,nproc;nowned=nowned+1;enddo
  allocate(row_ids(nowned),metric_rows(nowned,n),bad_metric_rows(nowned,n),c_cert_rows(nowned,ncert),&
    scaled_c_cert_rows(nowned,ncert),callback_expected_row_ids(nowned),callback_expected_c(nowned,ncert))
  metric=(0d0,0d0);full_coefficients=(0d0,0d0);construction_representation=(0d0,0d0)
  do i=1,n
    metric(i,i)=real(i,real64);full_coefficients(i,i)=1d0/sqrt(real(i,real64))
    construction_representation(i,i,1)=1d0
  enddo
  construction_representation(1,1,2)=1d0;construction_representation(2,2,2)=-1d0
  phase=acos(-1d0)/4d0
  construction_representation(1,1,3)=cmplx(cos(-phase),sin(-phase),kind=real64)
  construction_representation(2,2,3)=cmplx(cos(phase),sin(phase),kind=real64)
  full_eigenvalues=[-1d0,0.5d0,2d0,3d0]
  position=0
  do i=n,1,-1
    if(mod(i-1,nproc)/=rank)cycle
    position=position+1;row_ids(position)=int(i,int64)
    metric_rows(position,:)=metric(i,:);c_cert_rows(position,:)=full_coefficients(i,1:ncert)
  enddo

  call evaluate_dg_hybrid_low_energy_symmetry(comm,metric,construction_representation,full_coefficients,&
    full_eigenvalues,1,ncert,1d-12,occupied_defect,target_defect,energy_defect,ok,message)
  call require(ok,'certified LCFO prefix fixture is not symmetry closed: '//trim(message))
  call evaluate_dg_hybrid_low_energy_symmetry(comm,metric,construction_representation,full_coefficients,&
    full_eigenvalues,1,n,1d-12,occupied_defect,target_defect,energy_defect,ok,message)
  call require(.not.ok.and.target_defect>0.5d0,&
    'construction fixture is accidentally closed outside the certified LCFO prefix')

  certified_eigenvalues=full_eigenvalues(1:ncert)
  certified_representation=(0d0,0d0)
  certified_representation(1,1,1)=1d0;certified_representation(2,2,1)=1d0
  certified_representation(1,1,2)=1d0;certified_representation(2,2,2)=-1d0
  certified_representation(1,1,3)=cmplx(cos(-phase),sin(-phase),kind=real64)
  certified_representation(2,2,3)=cmplx(cos(phase),sin(phase),kind=real64)
  rotations=0d0
  do i=1,3;rotations(i,i,1)=1d0;enddo
  rotations(1,1,2)=-1d0;rotations(2,2,2)=-1d0;rotations(3,3,2)=1d0
  rotations(1,2,3)=-1d0;rotations(2,1,3)=1d0;rotations(3,3,3)=1d0
  scalar_operators=(0d0,0d0);scalar_operators(1,1,1)=2d0;scalar_operators(2,2,1)=1d0
  scalar_operators(2,2,2)=1d-14
  vector_operators=(0d0,0d0)
  vector_operators(1,2,1,1)=1d0;vector_operators(2,1,1,1)=1d0
  vector_operators(1,2,2,1)=cmplx(0d0,-1d0,kind=real64)
  vector_operators(2,1,2,1)=cmplx(0d0,1d0,kind=real64)
  vector_operators(1,1,3,1)=1d0;vector_operators(2,2,3,1)=-1d0
  tensor_operators=(0d0,0d0)
  do b=1,3
    do a=1,3
      tensor_operators(:,:,a,b,1)=matmul(vector_operators(:,:,a,1),vector_operators(:,:,b,1))
    enddo
  enddo
  call set_expected_u(expected_u)

  localizer_mode=0;localizer_calls=0;observed_require_unconstrained=.false.
  callback_expected_row_ids=row_ids;callback_expected_c=c_cert_rows
  verify_callback_payload=.true.;observed_callback_payload=.false.
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,certified_eigenvalues,1,&
    certified_representation,rotations,scalar_operators,vector_operators,tensor_operators,1d-11,&
    fixture_localizer,result,ok,message)
  call require(ok,'certified RT basis construction failed: '//trim(message))
  call require(localizer_calls==1.and.observed_require_unconstrained.and.&
    observed_certified_rank==ncert.and.observed_global_count==n,&
    'second unconstrained localizer was not invoked exactly once at fixed certified rank')
  call require(observed_callback_payload,&
    'second unconstrained localizer did not receive the certified row IDs and coefficients')
  verify_callback_payload=.false.
  call require(result%valid.and.result%localization_converged.and.&
    .not.result%localization_symmetry_constrained.and.result%localization_iterations==11,&
    'certified RT result lost the unconstrained convergence receipt')
  call require(result%global_count==n.and.result%certified_rank==ncert.and.result%noccupied==1.and.&
    all(result%owned_row_ids==row_ids).and.maxval(abs(result%c_cert-c_cert_rows))<1d-14.and.&
    maxval(abs(result%u_rt-expected_u))<1d-14.and.&
    maxval(abs(result%b_rt-matmul(c_cert_rows,expected_u)))<1d-14,&
    'certified construction-to-RT embedding is incorrect')

  identity=(0d0,0d0);do i=1,ncert;identity(i,i)=1d0;enddo
  expected_h=matmul(conjg(transpose(expected_u)),&
    matmul(diagonal_matrix(certified_eigenvalues),expected_u))
  expected_a(:,1)=conjg(expected_u(1,:))
  call require(maxval(abs(result%metric_rt-identity))<1d-14.and.&
    maxval(abs(result%hamiltonian_rt-expected_h))<1d-13.and.&
    maxval(abs(result%initial_occupied_amplitudes-expected_a))<1d-14,&
    'S_rt, H_rt(0), or A_occ(0) violates the certified gauge formula')
  call require(result%spread_improvement>0d0.and.result%spread_improvement<1d-12.and.&
    abs(result%spread_before_total-3d0)<1d-14.and.result%spread_after_total<result%spread_before_total,&
    'finite localization spread improvement was thresholded or recorded incorrectly')
  call require(max(result%transform_unitarity_defect,result%certified_metric_defect,&
    result%rt_metric_defect,result%embedding_defect,result%projector_invariance_defect,&
    result%target_symmetry_defect_before,result%target_symmetry_defect_after,&
    result%energy_symmetry_defect_before,result%energy_symmetry_defect_after,&
    result%symmetry_defect_invariance)<1d-11,&
    'second localization changed the certified projector or symmetry defects')

  do i=1,nscalar
    expected_scalar(:,:,i)=matmul(conjg(transpose(expected_u)),matmul(scalar_operators(:,:,i),expected_u))
  enddo
  do i=1,nvector;do a=1,3
    expected_vector(:,:,a,i)=matmul(conjg(transpose(expected_u)),&
      matmul(vector_operators(:,:,a,i),expected_u))
  enddo;enddo
  do i=1,ntensor;do b=1,3;do a=1,3
    expected_tensor(:,:,a,b,i)=matmul(conjg(transpose(expected_u)),&
      matmul(tensor_operators(:,:,a,b,i),expected_u))
  enddo;enddo;enddo
  call require(maxval(abs(result%scalar_operators_rt-expected_scalar))<1d-13.and.&
    maxval(abs(result%vector_operators_rt-expected_vector))<1d-13.and.&
    maxval(abs(result%tensor_operators_rt-expected_tensor))<1d-13,&
    'exact projected operators were not transformed into the RT gauge')
  call require(maxval(abs(result%scalar_operators_rt(:,:,2)-expected_scalar(:,:,2)))<1d-28.and.&
    maxval(abs(result%scalar_operators_rt(:,:,2)))>0d0,&
    'tiny exact projected operator family was independently pruned')
  call require(max(result%scalar_covariance_defect,result%vector_covariance_defect,&
    result%tensor_covariance_defect)<1d-11,&
    'scalar, vector, or rank-two tensor transformation law is violated')
  call require(abs(scalar_operators(2,2,2))>0d0.and.result%operator_fingerprint/=0_int64,&
    'tiny exact operator element was independently discarded')

  allocate(full_c_result(n,ncert),full_b_result(n,ncert));full_c_result=(0d0,0d0);full_b_result=(0d0,0d0)
  do i=1,nowned
    full_c_result(int(row_ids(i)),:)=result%c_cert(i,:);full_b_result(int(row_ids(i)),:)=result%b_rt(i,:)
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,full_c_result,size(full_c_result),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call MPI_Allreduce(MPI_IN_PLACE,full_b_result,size(full_b_result),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  projector_before=matmul(full_c_result,matmul(conjg(transpose(full_c_result)),metric))
  projector_after=matmul(full_b_result,matmul(conjg(transpose(full_b_result)),metric))
  call require(maxval(abs(projector_after-projector_before))<1d-12,&
    'localized certified basis changed the certified construction-space projector')

  tampered=result
  tampered%scalar_operators_rt=2d0*tampered%scalar_operators_rt
  call validate_dg_hybrid_certified_rt_basis(comm,tampered,1d-11,ok,message)
  call require(.not.ok,'covariance-preserving projected-operator tamper bypassed the fingerprint')

  tampered=result;tampered%spread_before_total=tampered%spread_before_total+1d0
  call validate_dg_hybrid_certified_rt_basis(comm,tampered,1d-11,ok,message)
  call require(.not.ok,'stored certified RT receipt tamper bypassed validation')

  tampered=result
  bits_before=transfer(tampered%spread_before_total,0_int64)
  bits_after=transfer(tampered%spread_after_total,0_int64)
  tampered%spread_before_total=transfer(ieor(bits_before,1_int64),0d0)
  tampered%spread_after_total=transfer(ieor(bits_after,ishftc(1_int64,18)),0d0)
  call validate_dg_hybrid_certified_rt_basis(comm,tampered,1d-9,ok,message)
  call require(.not.ok,'two-word certified RT receipt tamper collided in the fingerprint')

  nan_value_main=ieee_value(0d0,ieee_quiet_nan)
  call validate_dg_hybrid_certified_rt_basis(comm,result,nan_value_main,ok,message)
  call require(.not.ok,'nonfinite certified RT validation tolerance was accepted')
  tampered=result;tampered%spreads_after(1)=nan_value_main
  call validate_dg_hybrid_certified_rt_basis(comm,tampered,1d-11,ok,message)
  call require(.not.ok,'nonfinite stored localization spread was accepted')

  tampered=result;tampered%transform_unitarity_defect=-huge(1d0)
  call validate_dg_hybrid_certified_rt_basis(comm,tampered,1d-11,ok,message)
  call require(.not.ok.and.index(message,'invalid certified RT result payload')>0,&
    'finite negative huge stored defect was not rejected before validator arithmetic')

  tampered=result;tampered%spread_before_total=-huge(1d0)
  call validate_dg_hybrid_certified_rt_basis(comm,tampered,1d-11,ok,message)
  call require(.not.ok.and.index(message,'invalid certified RT result payload')>0,&
    'finite negative huge spread total was not rejected before validator arithmetic')

  localizer_mode=7;localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,certified_eigenvalues,1,&
    certified_representation,rotations,scalar_operators,vector_operators,tensor_operators,1d-11,&
    fixture_localizer,tampered,ok,message)
  call require(ok.and.tampered%spread_improvement<0d0.and.localizer_calls==1,&
    'finite converged fixed-rank localization was rejected because its spread did not improve')

  bad_metric_rows=metric_rows
  do i=1,nowned
    if(row_ids(i)==3_int64)bad_metric_rows(i,3)=-3d0
  enddo
  localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,bad_metric_rows,c_cert_rows,&
    certified_eigenvalues,1,certified_representation,rotations,scalar_operators,vector_operators,&
    tensor_operators,1d-11,fixture_localizer,tampered,ok,message)
  call require(.not.ok.and.localizer_calls==0.and.index(message,'singular or indefinite')>0,&
    'singular or indefinite construction metric was accepted')

  nowned_stress=0
  do i=rank+1,nstress,nproc;nowned_stress=nowned_stress+1;enddo
  allocate(stress_row_ids(nowned_stress),stress_metric(nstress,nstress),&
    stress_metric_rows(nowned_stress,nstress),stress_c_rows(nowned_stress,ncert))
  stress_metric=(0d0,0d0);stress_c_rows=(0d0,0d0)
  stress_metric(1,1)=(2d0**(-20))**2
  do i=2,nstress-1
    stress_metric(i,i)=(2d0**(-20))**2+0.5d0**2
    stress_metric(i,i-1)=2d0**(-21);stress_metric(i-1,i)=2d0**(-21)
  enddo
  stress_metric(nstress,nstress)=1d0
  stress_metric(nstress,1)=1d0;stress_metric(1,nstress)=1d0
  position=0
  do i=nstress,1,-1
    if(mod(i-1,nproc)/=rank)cycle
    position=position+1;stress_row_ids(position)=int(i,int64)
    stress_metric_rows(position,:)=stress_metric(i,:)
  enddo
  localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,nstress,stress_row_ids,stress_metric_rows,stress_c_rows,&
    certified_eigenvalues,1,certified_representation,rotations,scalar_operators,vector_operators,&
    tensor_operators,1d-11,fixture_localizer,tampered,ok,message)
  call require(.not.ok.and.localizer_calls==0.and.index(message,'singular or indefinite')>0,&
    'finite indefinite metric overflowed distributed Cholesky before safe rejection')

  bad_metric_rows=metric_rows
  do i=1,nowned
    if(row_ids(i)==3_int64)bad_metric_rows(i,3)=1d-12
  enddo
  localizer_mode=0;localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,bad_metric_rows,c_cert_rows,&
    certified_eigenvalues,1,certified_representation,rotations,scalar_operators,vector_operators,&
    tensor_operators,1d-11,fixture_localizer,tampered,ok,message)
  call require(ok.and.tampered%valid.and.localizer_calls==1,&
    'positive construction metric was coupled to the symmetry tolerance')

  bad_metric_rows=1d-20*metric_rows;scaled_c_cert_rows=1d10*c_cert_rows
  localizer_mode=0;localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,bad_metric_rows,scaled_c_cert_rows,&
    certified_eigenvalues,1,certified_representation,rotations,scalar_operators,vector_operators,&
    tensor_operators,1d-11,fixture_localizer,tampered,ok,message)
  call require(ok.and.tampered%valid.and.localizer_calls==1,&
    'well-conditioned positive metric was rejected after a uniform unit rescaling')

  metric_variant=metric
  metric_variant(1,2)=cmplx(0d0,-0.25d0,kind=real64)
  metric_variant(2,1)=conjg(metric_variant(1,2))
  metric_variant(2,2)=1d0+0.25d0**2
  scaled_c_cert_rows=(0d0,0d0)
  do i=1,nowned
    bad_metric_rows(i,:)=metric_variant(int(row_ids(i)),:)
    if(row_ids(i)==1_int64)then
      scaled_c_cert_rows(i,1)=1d0
      scaled_c_cert_rows(i,2)=cmplx(0d0,0.25d0,kind=real64)
    else if(row_ids(i)==2_int64)then
      scaled_c_cert_rows(i,2)=1d0
    endif
  enddo
  localizer_mode=0;localizer_calls=0
  callback_expected_row_ids=row_ids;callback_expected_c=scaled_c_cert_rows
  verify_callback_payload=.true.;observed_callback_payload=.false.
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,bad_metric_rows,scaled_c_cert_rows,&
    certified_eigenvalues,1,certified_representation,rotations,scalar_operators,vector_operators,&
    tensor_operators,1d-11,fixture_localizer,tampered,ok,message)
  call require(ok.and.tampered%valid.and.localizer_calls==1.and.observed_callback_payload,&
    'complex certified-support metric, coefficients, or row permutation was rejected')
  verify_callback_payload=.false.

  metric_variant(2,1)=metric_variant(2,1)+cmplx(0d0,2d-10,kind=real64)
  do i=1,nowned
    bad_metric_rows(i,:)=metric_variant(int(row_ids(i)),:)
  enddo
  localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,bad_metric_rows,c_cert_rows,&
    certified_eigenvalues,1,certified_representation,rotations,scalar_operators,vector_operators,&
    tensor_operators,1d-11,fixture_localizer,tampered,ok,message)
  call require(.not.ok.and.localizer_calls==0.and.index(message,'not Hermitian')>0,&
    'one-sided non-Hermitian metric perturbation was accepted')

  bad_metric_rows=1d200*metric_rows;scaled_c_cert_rows=1d-100*c_cert_rows
  localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,bad_metric_rows,scaled_c_cert_rows,&
    certified_eigenvalues,1,certified_representation,rotations,scalar_operators,vector_operators,&
    tensor_operators,1d-11,fixture_localizer,tampered,ok,message)
  call require(ok.and.tampered%valid.and.localizer_calls==1,&
    'large finite metric unit scaling overflowed the distributed certification')

  bad_metric_rows=1d-200*metric_rows;scaled_c_cert_rows=1d100*c_cert_rows
  localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,bad_metric_rows,scaled_c_cert_rows,&
    certified_eigenvalues,1,certified_representation,rotations,scalar_operators,vector_operators,&
    tensor_operators,1d-11,fixture_localizer,tampered,ok,message)
  call require(ok.and.tampered%valid.and.localizer_calls==1,&
    'small finite metric unit scaling underflowed the distributed certification')

  localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,certified_eigenvalues,1,&
    certified_representation,rotations,scalar_operators,vector_operators,tensor_operators,nan_value_main,&
    fixture_localizer,tampered,ok,message)
  call require(.not.ok.and.localizer_calls==0,&
    'nonfinite certified RT construction tolerance was accepted')

  bad_eigenvalues=[certified_eigenvalues(2),certified_eigenvalues(1)];localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,bad_eigenvalues,1,&
    certified_representation,rotations,scalar_operators,vector_operators,tensor_operators,1d-11,&
    fixture_localizer,tampered,ok,message)
  call require(.not.ok.and.localizer_calls==0,'unordered certified eigenvalues were accepted')

  bad_rotations=rotations;bad_rotations(1,1,1)=2d0;localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,certified_eigenvalues,1,&
    certified_representation,bad_rotations,scalar_operators,vector_operators,tensor_operators,1d-11,&
    fixture_localizer,tampered,ok,message)
  call require(.not.ok.and.localizer_calls==0,'nonorthogonal Cartesian rotation was accepted')

  bad_representation=certified_representation
  bad_representation(1,1,1)=cmplx(1d0+8d-12,0d0,kind=real64)
  zero_scalar_operators=(0d0,0d0);zero_vector_operators=(0d0,0d0);zero_tensor_operators=(0d0,0d0)
  localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,certified_eigenvalues,1,&
    bad_representation,rotations,zero_scalar_operators,zero_vector_operators,zero_tensor_operators,1d-11,&
    fixture_localizer,tampered,ok,message)
  call require(.not.ok.and.localizer_calls==1,&
    'normalized symmetry defect above the requested tolerance was accepted')

  bad_scalar_operators=(0d0,0d0)
  bad_scalar_operators(1,1,1)=1d0;bad_scalar_operators(2,2,1)=1d0
  bad_scalar_operators(1,2,1)=6d-12;localizer_mode=9;localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,certified_eigenvalues,1,&
    certified_representation,rotations,bad_scalar_operators,zero_vector_operators,&
    zero_tensor_operators,1d-11,fixture_localizer,tampered,ok,message)
  call require(.not.ok.and.localizer_calls==1,&
    'localized operator covariance defect above tolerance was hidden by matrix-size normalization')

  bad_vector_operators=(0d0,0d0);bad_vector_operators(1,1,1,1)=6d-12
  localizer_mode=9;localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,certified_eigenvalues,1,&
    certified_representation,rotations,zero_scalar_operators,bad_vector_operators,&
    zero_tensor_operators,1d-11,fixture_localizer,tampered,ok,message)
  call require(.not.ok.and.localizer_calls==1,&
    'noncovariant localized vector operator was accepted')

  bad_tensor_operators=(0d0,0d0);bad_tensor_operators(1,1,1,3,1)=6d-12
  localizer_mode=9;localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,certified_eigenvalues,1,&
    certified_representation,rotations,zero_scalar_operators,zero_vector_operators,&
    bad_tensor_operators,1d-11,fixture_localizer,tampered,ok,message)
  call require(.not.ok.and.localizer_calls==1,&
    'noncovariant localized rank-two tensor operator was accepted')

  bad_scalar_operators=cmplx(0.75d0*huge(1d0),0d0,kind=real64)
  localizer_mode=0;localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,certified_eigenvalues,1,&
    certified_representation,rotations,bad_scalar_operators,zero_vector_operators,&
    zero_tensor_operators,1d-11,fixture_localizer,tampered,ok,message)
  call require(.not.ok.and.localizer_calls==0,&
    'finite projected operator with an unsafe gauge-transform bound was accepted')

  do localizer_mode=1,6
    localizer_calls=0
    call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,certified_eigenvalues,1,&
      certified_representation,rotations,scalar_operators,vector_operators,tensor_operators,1d-11,&
      fixture_localizer,tampered,ok,message)
    call require(.not.ok.and..not.tampered%valid.and.localizer_calls==1,&
      'invalid second-localization result was accepted or retried')
    if(localizer_mode==6)call require(index(message,'intentional localizer failure')>0,&
      'localizer failure detail was discarded')
  enddo

  localizer_mode=8;localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,certified_eigenvalues,1,&
    certified_representation,rotations,scalar_operators,vector_operators,tensor_operators,1d-11,&
    fixture_localizer,tampered,ok,message)
  call require(.not.ok.and..not.tampered%valid.and.localizer_calls==1,&
    'nonfinite localization center was accepted')

  localizer_mode=10;localizer_calls=0
  call build_dg_hybrid_certified_rt_basis(comm,n,row_ids,metric_rows,c_cert_rows,certified_eigenvalues,1,&
    certified_representation,rotations,scalar_operators,vector_operators,tensor_operators,1d-11,&
    fixture_localizer,tampered,ok,message)
  call require(.not.ok.and..not.tampered%valid.and.localizer_calls==1,&
    'finite localization spreads with a nonfinite aggregate were accepted')

  if(rank==0)then
    write(*,'(a,i0,5(a,i0))')'HYBRID_CERTIFIED_RT_BASIS ranks=',nproc,&
      ' c=',result%c_cert_fingerprint,' u=',result%localization_fingerprint,&
      ' b=',result%b_rt_fingerprint,' operators=',result%operator_fingerprint,&
      ' fingerprint=',result%fingerprint
    write(*,'(a,i0,a)')'PASS certified Hybrid RT basis on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine fixture_localizer(comm_arg,global_count_arg,certified_rank_arg,row_ids_arg,c_cert_arg,&
      require_unconstrained,transform,centers,spreads_before,spreads_after,iterations,&
      symmetry_constrained,converged,callback_ok,callback_message)
    integer,intent(in)::comm_arg,global_count_arg,certified_rank_arg
    integer(int64),intent(in)::row_ids_arg(:)
    complex(real64),intent(in)::c_cert_arg(:,:)
    logical,intent(in)::require_unconstrained
    complex(real64),allocatable,intent(out)::transform(:,:)
    real(real64),allocatable,intent(out)::centers(:,:),spreads_before(:),spreads_after(:)
    integer,intent(out)::iterations
    logical,intent(out)::symmetry_constrained,converged,callback_ok
    character(*),intent(out)::callback_message
    real(real64)::nan_value
    localizer_calls=localizer_calls+1;observed_global_count=global_count_arg
    observed_certified_rank=certified_rank_arg;observed_require_unconstrained=require_unconstrained
    if(verify_callback_payload)then
      observed_callback_payload=.false.
      if(comm_arg==comm.and.size(row_ids_arg)==size(callback_expected_row_ids).and.&
          all(shape(c_cert_arg)==shape(callback_expected_c)))then
        observed_callback_payload=all(row_ids_arg==callback_expected_row_ids).and.&
          all(abs(c_cert_arg-callback_expected_c)<1d-14)
      endif
    endif
    if(localizer_mode==2)then
      allocate(transform(ncert,ncert-1),centers(3,ncert-1),spreads_before(ncert-1),spreads_after(ncert-1))
      transform=expected_u(:,1:ncert-1);centers=0d0;spreads_before=1d0;spreads_after=0.5d0
    else
      allocate(transform(ncert,ncert),centers(3,ncert),spreads_before(ncert),spreads_after(ncert))
      transform=expected_u;centers=reshape([0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0],[3,ncert])
      spreads_before=[2d0,1d0];spreads_after=[2d0-spacing(2d0),1d0]
    endif
    if(localizer_mode==7)spreads_after=[2.25d0,1d0]
    if(localizer_mode==9)then
      transform=(0d0,0d0);transform(1,1)=1d0;transform(2,2)=1d0
    endif
    if(localizer_mode==10)then
      spreads_before=0.75d0*huge(1d0);spreads_after=0.75d0*huge(1d0)
    endif
    if(localizer_mode==1)transform(1,1)=2d0*transform(1,1)
    if(localizer_mode==5.and.rank==0)then
      nan_value=ieee_value(0d0,ieee_quiet_nan);spreads_after(1)=nan_value
    endif
    if(localizer_mode==8.and.rank==0)then
      nan_value=ieee_value(0d0,ieee_quiet_nan);centers(1,1)=nan_value
    endif
    iterations=11;symmetry_constrained=localizer_mode==4
    converged=.true.;if(localizer_mode==3.and.rank==0)converged=.false.
    callback_ok=.true.;if(localizer_mode==6.and.rank==0)callback_ok=.false.
    callback_message='';if(.not.callback_ok)callback_message='intentional localizer failure'
  end subroutine fixture_localizer

  subroutine set_expected_u(value_arg)
    complex(real64),intent(out)::value_arg(ncert,ncert)
    value_arg=(0d0,0d0)
    value_arg(1,1)=1d0/sqrt(2d0);value_arg(1,2)=cmplx(0d0,1d0,kind=real64)/sqrt(2d0)
    value_arg(2,1)=1d0/sqrt(2d0);value_arg(2,2)=cmplx(0d0,-1d0,kind=real64)/sqrt(2d0)
  end subroutine set_expected_u

  function diagonal_matrix(values)result(matrix)
    real(real64),intent(in)::values(:)
    complex(real64)::matrix(size(values),size(values))
    integer::index
    matrix=(0d0,0d0);do index=1,size(values);matrix(index,index)=values(index);enddo
  end function diagonal_matrix

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(label)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require
end program test_dg_hybrid_certified_rt_basis_mpi
