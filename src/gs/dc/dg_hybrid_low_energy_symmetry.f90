#include "config.h"
module dg_hybrid_low_energy_symmetry
  use mpi
  use,intrinsic::iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::select_dg_hybrid_symmetry_target,evaluate_dg_hybrid_low_energy_symmetry
contains
  subroutine select_dg_hybrid_symmetry_target(comm,eigenvalues,requested_rank,tolerance,target_rank,ok,message)
    integer,intent(in)::comm,requested_rank
    real(real64),intent(in)::eigenvalues(:),tolerance
    integer,intent(out)::target_rank
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::ierr,rank,local_bad,global_bad,minimum_request,maximum_request
    real(real64),allocatable::reference(:)
    real(real64)::local_difference,global_difference,scale

    ok=.false.;message='';target_rank=0;local_bad=0
    if(requested_rank<1.or.requested_rank>size(eigenvalues).or.tolerance<=0d0.or.&
      .not.all(ieee_is_finite(eigenvalues)).or..not.ieee_is_finite(tolerance))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid low-energy symmetry target contract';return
    endif
    call MPI_Allreduce(requested_rank,minimum_request,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(requested_rank,maximum_request,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_request/=maximum_request)then
      message='rank-disagreeing low-energy symmetry target';return
    endif
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then;message='target rank query failed';return;endif
    allocate(reference(size(eigenvalues)));reference=eigenvalues
    call MPI_Bcast(reference,size(reference),MPI_DOUBLE_PRECISION,0,comm,ierr)
    local_difference=maxval(abs(eigenvalues-reference))
    call MPI_Allreduce(local_difference,global_difference,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_difference>tolerance)then
      message='rank-disagreeing low-energy spectrum';return
    endif
    target_rank=requested_rank
    do while(target_rank<size(eigenvalues))
      scale=max(1d0,abs(eigenvalues(target_rank)),abs(eigenvalues(target_rank+1)))
      if(abs(eigenvalues(target_rank+1)-eigenvalues(target_rank))>tolerance*scale)exit
      target_rank=target_rank+1
    enddo
    if(target_rank==size(eigenvalues).and.requested_rank<size(eigenvalues))then
      scale=max(1d0,abs(eigenvalues(target_rank)),abs(eigenvalues(target_rank-1)))
      if(abs(eigenvalues(target_rank)-eigenvalues(target_rank-1))<=tolerance*scale)then
        message='low-energy symmetry target ends in an unresolved degenerate cluster';return
      endif
    endif
    ok=.true.
  end subroutine select_dg_hybrid_symmetry_target

  subroutine evaluate_dg_hybrid_low_energy_symmetry(comm,metric,representation,coefficients,eigenvalues,&
      occupied_rank,target_rank,tolerance,occupied_defect,target_defect,energy_defect,ok,message)
    integer,intent(in)::comm,occupied_rank,target_rank
    complex(real64),intent(in)::metric(:,:),representation(:,:,:),coefficients(:,:)
    real(real64),intent(in)::eigenvalues(:),tolerance
    real(real64),intent(out)::occupied_defect,target_defect,energy_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::target_coefficients(:,:),occupied_coefficients(:,:),target_gram(:,:),&
      occupied_gram(:,:),target_action(:,:),occupied_action(:,:),target_identity(:,:),occupied_identity(:,:),&
      target_energy(:,:),transformed_energy(:,:)
    integer::n,nstate,noperation,operation,i,ierr,local_bad,global_bad
    real(real64)::orthogonality_defect,local_values(3),global_values(3),target_scale,occupied_scale,energy_scale

    occupied_defect=huge(1d0);target_defect=huge(1d0);energy_defect=huge(1d0)
    ok=.false.;message='';local_bad=0;n=size(metric,1);nstate=size(coefficients,2);noperation=size(representation,3)
    if(n<1.or.size(metric,2)/=n.or.any(shape(representation)/=[n,n,noperation]).or.noperation<1.or.&
      size(coefficients,1)/=n.or.size(eigenvalues)/=nstate.or.occupied_rank<1.or.&
      occupied_rank>target_rank.or.target_rank>nstate.or.tolerance<=0d0.or.&
      .not.finite_matrix(metric).or..not.finite_matrix(coefficients).or.&
      .not.finite_tensor(representation).or..not.all(ieee_is_finite(eigenvalues)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid low-energy symmetry evaluation contract';return
    endif
    allocate(target_coefficients(n,target_rank),occupied_coefficients(n,occupied_rank),&
      target_gram(target_rank,target_rank),occupied_gram(occupied_rank,occupied_rank),&
      target_action(target_rank,target_rank),occupied_action(occupied_rank,occupied_rank),&
      target_identity(target_rank,target_rank),occupied_identity(occupied_rank,occupied_rank),&
      target_energy(target_rank,target_rank),transformed_energy(target_rank,target_rank))
    target_coefficients=coefficients(:,1:target_rank);occupied_coefficients=coefficients(:,1:occupied_rank)
    target_identity=(0d0,0d0);occupied_identity=(0d0,0d0);target_energy=(0d0,0d0)
    do i=1,target_rank
      target_identity(i,i)=1d0;target_energy(i,i)=eigenvalues(i)
    enddo
    do i=1,occupied_rank;occupied_identity(i,i)=1d0;enddo
    target_gram=matmul(conjg(transpose(target_coefficients)),matmul(metric,target_coefficients))
    occupied_gram=matmul(conjg(transpose(occupied_coefficients)),matmul(metric,occupied_coefficients))
    orthogonality_defect=max(frobenius(target_gram-target_identity)/sqrt(real(target_rank,real64)),&
      frobenius(occupied_gram-occupied_identity)/sqrt(real(occupied_rank,real64)))
    if(.not.ieee_is_finite(orthogonality_defect).or.orthogonality_defect>tolerance)then
      message='low-energy symmetry coefficients are not metric orthonormal';return
    endif
    local_values=0d0
    target_scale=sqrt(real(target_rank,real64));occupied_scale=sqrt(real(occupied_rank,real64))
    energy_scale=max(1d0,frobenius(target_energy))
    do operation=1,noperation
      target_action=matmul(conjg(transpose(target_coefficients)),&
        matmul(metric,matmul(representation(:,:,operation),target_coefficients)))
      occupied_action=matmul(conjg(transpose(occupied_coefficients)),&
        matmul(metric,matmul(representation(:,:,operation),occupied_coefficients)))
      transformed_energy=matmul(conjg(transpose(target_action)),matmul(target_energy,target_action))
      local_values(1)=max(local_values(1),frobenius(matmul(conjg(transpose(occupied_action)),occupied_action)-&
        occupied_identity)/occupied_scale)
      local_values(2)=max(local_values(2),frobenius(matmul(conjg(transpose(target_action)),target_action)-&
        target_identity)/target_scale)
      local_values(3)=max(local_values(3),frobenius(transformed_energy-target_energy)/energy_scale)
    enddo
    call MPI_Allreduce(local_values,global_values,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(.not.ieee_is_finite(global_values)))then
      message='low-energy symmetry defect reduction failed';return
    endif
    occupied_defect=global_values(1);target_defect=global_values(2);energy_defect=global_values(3)
    ok=maxval(global_values)<=tolerance
    if(ok)then;message='';else;message='low-energy LCFO eigenspace is not symmetry closed';endif
  end subroutine evaluate_dg_hybrid_low_energy_symmetry

  pure logical function finite_matrix(values)
    complex(real64),intent(in)::values(:,:)
    finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_matrix

  pure logical function finite_tensor(values)
    complex(real64),intent(in)::values(:,:,:)
    finite_tensor=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_tensor

  pure real(real64) function frobenius(values)
    complex(real64),intent(in)::values(:,:)
    frobenius=sqrt(sum(abs(values)**2))
  end function frobenius
end module dg_hybrid_low_energy_symmetry
