! Test-only legacy reference; excluded from the SALMON production build.
#include "config.h"
module dg_hybrid_real_space_residual
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::evaluate_dg_hybrid_real_space_residual,evaluate_dg_hybrid_face_action_residuals
contains
  subroutine evaluate_dg_hybrid_face_action_residuals(comm,global_basis_count,row_ids,component_rows,&
      hamiltonian_rows,coefficients,metric_eigen_action,interface_scale,reconstructed_actions,residuals,ok,message)
    integer,intent(in)::comm,global_basis_count
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::component_rows(:,:,:),hamiltonian_rows(:,:),coefficients(:,:),&
      metric_eigen_action(:,:),reconstructed_actions(:,:,:)
    real(real64),intent(in)::interface_scale
    real(real64),intent(out)::residuals(3)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,component,ierr,nowned,nstate,local_bad,global_bad
    integer,allocatable::row_ownership(:)
    complex(real64),allocatable::local_coefficients(:,:),global_coefficients(:,:),volume_rows(:,:),&
      volume_actions(:,:),candidate_actions(:,:),reference_component(:,:)
    real(real64)::local_norms(3,3),global_norms(3,3)
    ok=.false.;message='';residuals=huge(1d0)
    nowned=size(row_ids);nstate=size(coefficients,2);local_bad=0
    if(global_basis_count<1.or.nstate<1.or.any(shape(component_rows)/=[nowned,global_basis_count,3]).or.&
        any(shape(hamiltonian_rows)/=[nowned,global_basis_count]).or.size(coefficients,1)/=nowned.or.&
        any(shape(metric_eigen_action)/=[nowned,nstate]).or.&
        any(shape(reconstructed_actions)/=[nowned,nstate,3]).or..not.ieee_is_finite(interface_scale).or.&
        any(row_ids<1_int64).or.any(row_ids>int(global_basis_count,int64)).or.&
        .not.all(ieee_is_finite(real(component_rows))).or..not.all(ieee_is_finite(aimag(component_rows))).or.&
        .not.finite_complex(hamiltonian_rows).or..not.finite_complex(coefficients).or.&
        .not.finite_complex(metric_eigen_action).or..not.all(ieee_is_finite(real(reconstructed_actions))).or.&
        .not.all(ieee_is_finite(aimag(reconstructed_actions))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid SIPG face-action residual contract';return;endif
    allocate(row_ownership(global_basis_count));row_ownership=0
    do i=1,nowned;row_ownership(int(row_ids(i)))=row_ownership(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,row_ownership,global_basis_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(row_ownership/=1))then
      message='SIPG face-action residual ownership is incomplete';return
    endif
    allocate(local_coefficients(global_basis_count,nstate),global_coefficients(global_basis_count,nstate))
    local_coefficients=(0d0,0d0)
    do i=1,nowned;local_coefficients(int(row_ids(i)),:)=coefficients(i,:);enddo
    call MPI_Allreduce(local_coefficients,global_coefficients,size(global_coefficients),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='SIPG face-action coefficient assembly failed';return;endif
    allocate(volume_rows(nowned,global_basis_count),volume_actions(nowned,nstate),&
      candidate_actions(nowned,nstate),reference_component(nowned,nstate))
    volume_rows=hamiltonian_rows-interface_scale*sum(component_rows,dim=3)
    volume_actions=matmul(volume_rows,global_coefficients);local_norms=0d0
    do component=1,3
      candidate_actions=volume_actions
      do i=1,3
        if(i==component)then
          candidate_actions=candidate_actions+interface_scale*reconstructed_actions(:,:,i)
        else
          reference_component=matmul(component_rows(:,:,i),global_coefficients)
          candidate_actions=candidate_actions+interface_scale*reference_component
        endif
      enddo
      local_norms(1,component)=sum(abs(candidate_actions-metric_eigen_action)**2)
      local_norms(2,component)=sum(abs(candidate_actions)**2)
      local_norms(3,component)=sum(abs(metric_eigen_action)**2)
    enddo
    call MPI_Allreduce(local_norms,global_norms,9,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='SIPG face-action residual reduction failed';return;endif
    do component=1,3
      residuals(component)=sqrt(global_norms(1,component))/&
        max(1d0,sqrt(global_norms(2,component)),sqrt(global_norms(3,component)))
    enddo
    ok=all(ieee_is_finite(residuals))
    if(ok)then;message='';else;message='nonfinite SIPG face-action residual';endif
#else
    residuals=huge(1d0);ok=.false.;message='SIPG face-action residual requires MPI'
#endif
  end subroutine evaluate_dg_hybrid_face_action_residuals

  subroutine evaluate_dg_hybrid_real_space_residual(comm,global_point_count,point_ids,weights,&
      global_basis_count,row_ids,basis_values,full_action_values,coefficients,eigenvalues,residual,ok,message)
    integer,intent(in)::comm,global_point_count,global_basis_count
    integer(int64),intent(in)::point_ids(:),row_ids(:)
    real(real64),intent(in)::weights(:),eigenvalues(:)
    complex(real64),intent(in)::basis_values(:,:),full_action_values(:,:),coefficients(:,:)
    real(real64),intent(out)::residual
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,ierr,nlocal,nowned,nstate,local_bad,global_bad
    integer,allocatable::point_ownership(:),row_ownership(:)
    complex(real64),allocatable::local_coefficients(:,:),global_coefficients(:,:),psi(:,:),hpsi(:,:),eps_psi(:,:)
    real(real64)::local_norms(3),global_norms(3)
    ok=.false.;message='';residual=huge(1d0)
    nlocal=size(point_ids);nowned=size(row_ids);nstate=size(eigenvalues);local_bad=0
    if(global_point_count<1.or.global_basis_count<1.or.nstate<1.or.size(weights)/=nlocal.or.&
        any(shape(basis_values)/=[global_basis_count,nlocal]).or.&
        any(shape(full_action_values)/=shape(basis_values)).or.&
        any(shape(coefficients)/=[nowned,nstate]).or.any(point_ids<1_int64).or.&
        any(point_ids>int(global_point_count,int64)).or.any(row_ids<1_int64).or.&
        any(row_ids>int(global_basis_count,int64)).or.any(weights<=0d0).or.&
        .not.all(ieee_is_finite(weights)).or..not.all(ieee_is_finite(eigenvalues)).or.&
        .not.finite_complex(basis_values).or..not.finite_complex(full_action_values).or.&
        .not.finite_complex(coefficients))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid real-space DG residual contract';return;endif
    allocate(point_ownership(global_point_count),row_ownership(global_basis_count));point_ownership=0;row_ownership=0
    do i=1,nlocal;point_ownership(int(point_ids(i)))=point_ownership(int(point_ids(i)))+1;enddo
    do i=1,nowned;row_ownership(int(row_ids(i)))=row_ownership(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,point_ownership,global_point_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(MPI_IN_PLACE,row_ownership,global_basis_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(point_ownership/=1).or.any(row_ownership/=1))then
      message='real-space DG residual ownership is incomplete';return
    endif
    allocate(local_coefficients(global_basis_count,nstate),global_coefficients(global_basis_count,nstate))
    local_coefficients=(0d0,0d0)
    do i=1,nowned;local_coefficients(int(row_ids(i)),:)=coefficients(i,:);enddo
    call MPI_Allreduce(local_coefficients,global_coefficients,size(global_coefficients),&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='real-space DG coefficient assembly failed';return;endif
    allocate(psi(nstate,nlocal),hpsi(nstate,nlocal),eps_psi(nstate,nlocal))
    psi=matmul(transpose(global_coefficients),basis_values)
    hpsi=matmul(transpose(global_coefficients),full_action_values)
    eps_psi=psi
    do i=1,nstate;eps_psi(i,:)=eigenvalues(i)*eps_psi(i,:);enddo
    local_norms=0d0
    do i=1,nstate
      local_norms(1)=local_norms(1)+sum(weights*abs(hpsi(i,:)-eps_psi(i,:))**2)
      local_norms(2)=local_norms(2)+sum(weights*abs(hpsi(i,:))**2)
      local_norms(3)=local_norms(3)+sum(weights*abs(eps_psi(i,:))**2)
    enddo
    call MPI_Allreduce(local_norms,global_norms,3,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='real-space DG residual reduction failed';return;endif
    residual=sqrt(global_norms(1))/max(1d0,sqrt(global_norms(2)),sqrt(global_norms(3)))
    ok=ieee_is_finite(residual)
    if(ok)then;message='';else;message='nonfinite real-space DG residual';endif
#else
    residual=huge(1d0);ok=.false.;message='real-space DG residual requires MPI'
#endif
  end subroutine evaluate_dg_hybrid_real_space_residual

  logical function finite_complex(values)
    complex(real64),intent(in)::values(:,:)
    finite_complex=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex
end module dg_hybrid_real_space_residual
