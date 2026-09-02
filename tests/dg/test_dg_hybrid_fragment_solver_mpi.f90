#include "config.h"
program test_dg_hybrid_fragment_solver_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis,build_dg_hybrid_fragment_basis
  use dg_hybrid_fragment_solver,only:solve_dg_hybrid_fragment_basis,&
    solve_dg_hybrid_fragment_spectrum,reconstruct_dg_hybrid_fragment_density
  implicit none
  integer,parameter::npoint=4,nbasis=4,nstate=2
  integer::comm,rank,nproc,ierr,nowned,i,j,k,h_calls,s_calls
  integer(int64),allocatable::wf_ids(:),pw_ids(:)
  complex(real64),allocatable::wf_values(:,:),pw_values(:,:),coefficients(:,:),&
    wrapper_coefficients(:,:),invalid_coefficients(:,:)
  complex(real64)::vectors(npoint,nbasis),eigenvectors(npoint,nbasis),hop(npoint,npoint),sop(npoint,npoint)
  real(real64)::occupations(nstate),invalid_occupations(nstate),point_weights(npoint),&
    reference_weights(npoint),eigenvalues(nstate),wrapper_eigenvalues(nstate),core_norms(nstate),&
    density(npoint),wrapper_density(npoint),residual,wrapper_residual,orthogonality,&
    wrapper_orthogonality,electron_count,wrapper_electron_count,invalid_weights(npoint)
  real(real64),allocatable::invalid_eigenvalues(:),invalid_core_norms(:)
  logical::core_mask(npoint),invalid_core_mask(npoint),ok
  type(s_dg_hybrid_fragment_basis)::basis,invalid_basis
  integer(int64)::workspace,wrapper_workspace,fingerprint,wrapper_fingerprint
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  h_calls=0;s_calls=0;reference_weights=[0.5d0,2d0,1d0,1d0]
  vectors=(0d0,0d0);do j=1,nbasis;vectors(j,j)=1d0;enddo
  eigenvectors=(0d0,0d0)
  eigenvectors(1,1)=sqrt(0.75d0);eigenvectors(3,1)=sqrt(0.25d0)
  eigenvectors(1,3)=-sqrt(0.25d0);eigenvectors(3,3)=sqrt(0.75d0)
  eigenvectors(2,2)=sqrt(0.25d0);eigenvectors(4,2)=sqrt(0.75d0)
  eigenvectors(2,4)=-sqrt(0.75d0);eigenvectors(4,4)=sqrt(0.25d0)
  hop=(0d0,0d0);sop=(0d0,0d0)
  do j=1,npoint;sop(j,j)=1d0;enddo
  do i=1,npoint;do j=1,npoint;do k=1,nbasis
    hop(i,j)=hop(i,j)+eigenvectors(i,k)*real(k,real64)*conjg(eigenvectors(j,k))*&
      sqrt(reference_weights(j)/reference_weights(i))
  enddo;enddo;enddo
  nowned=count([(mod(j-1,nproc)==rank,j=1,nbasis)])
  allocate(wf_ids(count([(mod(j-1,nproc)==rank.and.j<=2,j=1,nbasis)])),&
    pw_ids(count([(mod(j-1,nproc)==rank.and.j>2,j=1,nbasis)])))
  allocate(wf_values(npoint,size(wf_ids)),pw_values(npoint,size(pw_ids)));j=0;k=0
  if(size(wf_ids)>0)wf_values=(0d0,0d0);if(size(pw_ids)>0)pw_values=(0d0,0d0)
  do nowned=1,nbasis
    if(mod(nowned-1,nproc)/=rank)cycle
    if(nowned<=2)then;j=j+1;wf_ids(j)=100+nstated(nowned);wf_values(:,j)=vectors(:,nowned)
    else;k=k+1;pw_ids(k)=100+nstated(nowned);pw_values(:,k)=vectors(:,nowned);endif
  enddo
  call build_dg_hybrid_fragment_basis(comm,1,wf_ids,wf_values,pw_ids,pw_values,0,0,basis,ok,message)
  call require(ok,trim(message));occupations=[1.5d0,0.5d0]
  do i=1,npoint
    basis%buffer_point_ids(i)=int(mod(i+rank-1,npoint)+1,int64)
    point_weights(i)=reference_weights(int(basis%buffer_point_ids(i)))
    core_mask(i)=basis%buffer_point_ids(i)<=2_int64
    do j=1,size(basis%global_ids)
      basis%buffer_values(i,j)=vectors(int(basis%buffer_point_ids(i)),basis_ordinal(basis%global_ids(j)))
    enddo
  enddo
  call solve_dg_hybrid_fragment_spectrum(comm,basis,nstate,core_mask,point_weights,apply_h,apply_s,1d-12,&
    coefficients,eigenvalues,core_norms,residual,orthogonality,workspace,fingerprint,ok,message)
  call require(ok,trim(message));call require(maxval(abs(eigenvalues-[1d0,2d0]))<1d-12,'fragment eigenvalues mismatch')
  call require(residual<1d-12.and.orthogonality<1d-12,'fragment eigensystem receipts mismatch')
  call require(maxval(abs(core_norms-[0.75d0,0.25d0]))<1d-12,'fragment core norms mismatch')
  call reconstruct_dg_hybrid_fragment_density(comm,basis,coefficients,occupations,2d0,core_mask,&
    point_weights,density,electron_count,ok,message)
  call require(ok,trim(message))
  call require(abs(electron_count-1.25d0)<1d-12.and.&
    abs(electron_count-dot_product(occupations,core_norms))<1d-12.and.&
    abs(density(findloc(basis%buffer_point_ids,1_int64,dim=1))-&
      2.25d0)<1d-12.and.&
    abs(density(findloc(basis%buffer_point_ids,2_int64,dim=1))-&
      0.0625d0)<1d-12.and.count(abs(density)>1d-12)==2,'fragment split core density mismatch')
  call require(h_calls==1.and.s_calls==1,'density phase unexpectedly reapplied an operator')
  call require(size(coefficients,1)==size(basis%global_ids).and.size(coefficients,2)==nstate,&
    'fragment coefficients are not basis-row distributed')

  call solve_dg_hybrid_fragment_basis(comm,basis,nstate,occupations,core_mask,point_weights,apply_h,apply_s,1d-12,&
    wrapper_coefficients,wrapper_eigenvalues,wrapper_density,wrapper_electron_count,wrapper_residual,&
    wrapper_orthogonality,wrapper_workspace,wrapper_fingerprint,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(wrapper_coefficients-coefficients))<1d-13.and.&
    maxval(abs(wrapper_eigenvalues-eigenvalues))<1d-13.and.maxval(abs(wrapper_density-density))<1d-13.and.&
    abs(wrapper_electron_count-electron_count)<1d-13.and.wrapper_residual==residual.and.&
    wrapper_orthogonality==orthogonality.and.wrapper_workspace==workspace.and.workspace>0_int64.and.&
    wrapper_fingerprint==fingerprint,&
    'fragment compatibility wrapper differs from split phases')
  call require(h_calls==2.and.s_calls==2,'compatibility wrapper did not perform exactly one spectrum solve')

  invalid_occupations=occupations;if(rank==0)invalid_occupations(1)=-1d-3
  call reconstruct_dg_hybrid_fragment_density(comm,basis,coefficients,invalid_occupations,2d0,core_mask,&
    point_weights,density,electron_count,ok,message)
  call require(.not.ok,'negative fragment occupation was accepted')
  invalid_occupations=occupations;if(rank==0)invalid_occupations(1)=2d0+1d-3
  call reconstruct_dg_hybrid_fragment_density(comm,basis,coefficients,invalid_occupations,2d0,core_mask,&
    point_weights,density,electron_count,ok,message)
  call require(.not.ok,'fragment occupation above spin degeneracy was accepted')
  invalid_occupations=occupations
  if(rank==0)invalid_occupations(1)=ieee_value(0d0,ieee_quiet_nan)
  call reconstruct_dg_hybrid_fragment_density(comm,basis,coefficients,invalid_occupations,2d0,core_mask,&
    point_weights,density,electron_count,ok,message)
  call require(.not.ok,'non-finite fragment occupation was accepted')

  call reconstruct_dg_hybrid_fragment_density(comm,basis,coefficients(:,1:1),occupations,2d0,core_mask,&
    point_weights,density,electron_count,ok,message)
  call require(.not.ok,'fragment coefficient/occupation extent mismatch was accepted')
  if(nproc>1)then
    if(rank==0)then
      allocate(invalid_coefficients(size(coefficients,1),1))
    else
      allocate(invalid_coefficients(size(coefficients,1),nstate))
    endif
    invalid_coefficients=coefficients(:,1:size(invalid_coefficients,2))
    call reconstruct_dg_hybrid_fragment_density(comm,basis,invalid_coefficients,occupations,2d0,core_mask,&
      point_weights,density,electron_count,ok,message)
    call require(.not.ok,'rank-local fragment coefficient extent mismatch was accepted')
    deallocate(invalid_coefficients)
  endif
  call reconstruct_dg_hybrid_fragment_density(comm,basis,coefficients,occupations,2d0,core_mask,&
    point_weights(:npoint-1),density,electron_count,ok,message)
  call require(.not.ok,'fragment point-weight extent mismatch was accepted')

  allocate(invalid_coefficients,source=coefficients)
  if(rank==0.and.size(invalid_coefficients)>0)&
    invalid_coefficients(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,kind=real64)
  call reconstruct_dg_hybrid_fragment_density(comm,basis,invalid_coefficients,occupations,2d0,core_mask,&
    point_weights,density,electron_count,ok,message)
  call require(.not.ok,'non-finite fragment coefficient was accepted')
  deallocate(invalid_coefficients)

  invalid_basis=basis
  if(nproc==1)then
    invalid_basis%global_ids(2)=invalid_basis%global_ids(1)
  elseif(rank>0.and.size(invalid_basis%global_ids)>0)then
    invalid_basis%global_ids(1)=100_int64+nstated(1)
  endif
  call reconstruct_dg_hybrid_fragment_density(comm,invalid_basis,coefficients,occupations,2d0,core_mask,&
    point_weights,density,electron_count,ok,message)
  call require(.not.ok,'duplicate fragment basis-row ownership was accepted')

  invalid_basis=basis
  if(rank==0)invalid_basis%buffer_point_ids(2)=invalid_basis%buffer_point_ids(1)
  call reconstruct_dg_hybrid_fragment_density(comm,invalid_basis,coefficients,occupations,2d0,core_mask,&
    point_weights,density,electron_count,ok,message)
  call require(.not.ok,'duplicate fragment buffer point ID was accepted')

  invalid_basis=basis;invalid_basis%provenance_fingerprint=0_int64
  call reconstruct_dg_hybrid_fragment_density(comm,invalid_basis,coefficients,occupations,2d0,core_mask,&
    point_weights,density,electron_count,ok,message)
  call require(.not.ok,'unfinalized fragment basis provenance was accepted')

  if(nproc>1)then
    invalid_weights=point_weights
    if(rank==0)invalid_weights(findloc(basis%buffer_point_ids,2_int64,dim=1))=&
      invalid_weights(findloc(basis%buffer_point_ids,2_int64,dim=1))+1d-3
    call reconstruct_dg_hybrid_fragment_density(comm,basis,coefficients,occupations,2d0,core_mask,&
      invalid_weights,density,electron_count,ok,message)
    call require(.not.ok,'rank-inconsistent fragment point weight was accepted')
    invalid_core_mask=core_mask
    if(rank==0)invalid_core_mask(findloc(basis%buffer_point_ids,2_int64,dim=1))=.false.
    call reconstruct_dg_hybrid_fragment_density(comm,basis,coefficients,occupations,2d0,invalid_core_mask,&
      point_weights,density,electron_count,ok,message)
    call require(.not.ok,'rank-inconsistent fragment core ownership was accepted')
    call reconstruct_dg_hybrid_fragment_density(comm,basis,coefficients,occupations,&
      merge(2d0,1.5d0,rank==0),core_mask,point_weights,density,electron_count,ok,message)
    call require(.not.ok,'rank-inconsistent maximum occupation was accepted')
    allocate(invalid_eigenvalues(nstate),invalid_core_norms(nstate))
    call solve_dg_hybrid_fragment_spectrum(comm,basis,nstate,core_mask,point_weights,apply_h,apply_s,&
      merge(1d-12,1.5d-12,rank==0),invalid_coefficients,invalid_eigenvalues,invalid_core_norms,&
      residual,orthogonality,wrapper_workspace,wrapper_fingerprint,ok,message)
    call require(.not.ok,'rank-inconsistent fragment spectrum tolerance was accepted')
    call require(h_calls==2.and.s_calls==2,'invalid spectrum controls reached fragment operators')
  endif
  if(rank==0)write(*,'(a,i0,a,i0,a,i0)')'HYBRID_FRAGMENT_SOLVER ranks=',nproc,&
    ' fingerprint=',fingerprint,' workspace=',workspace
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid fragment solver on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  integer(int64) function nstated(value)
    integer,intent(in)::value
    integer(int64),parameter::noncontiguous_ids(nbasis)=[1_int64,3_int64,6_int64,10_int64]
    nstated=noncontiguous_ids(value)
  end function nstated
  integer function basis_ordinal(global_id)
    integer(int64),intent(in)::global_id
    integer::value
    basis_ordinal=0
    do value=1,nbasis
      if(global_id==100_int64+nstated(value))basis_ordinal=value
    enddo
  end function basis_ordinal
  subroutine apply_h(input,output,callback_ok)
    complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:);logical,intent(out)::callback_ok
    complex(real64)::global_input(npoint,size(input,2)),global_output(npoint,size(input,2));integer::p
    h_calls=h_calls+1;global_input=(0d0,0d0)
    do p=1,size(input,1);global_input(int(basis%buffer_point_ids(p)),:)=input(p,:);enddo
    global_output=matmul(hop,global_input)
    do p=1,size(output,1);output(p,:)=global_output(int(basis%buffer_point_ids(p)),:);enddo
    callback_ok=.true.
  end subroutine apply_h
  subroutine apply_s(input,output,callback_ok)
    complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:);logical,intent(out)::callback_ok
    complex(real64)::global_input(npoint,size(input,2)),global_output(npoint,size(input,2));integer::p
    s_calls=s_calls+1;global_input=(0d0,0d0)
    do p=1,size(input,1);global_input(int(basis%buffer_point_ids(p)),:)=input(p,:);enddo
    global_output=matmul(sop,global_input)
    do p=1,size(output,1);output(p,:)=global_output(int(basis%buffer_point_ids(p)),:);enddo
    callback_ok=.true.
  end subroutine apply_s
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text;logical::global_condition
    call MPI_Allreduce(condition,global_condition,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(.not.global_condition)then;if(rank==0)write(0,'(a)')trim(text);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
end program test_dg_hybrid_fragment_solver_mpi
