#include "config.h"
module dg_hybrid_low_energy_symmetry
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_ground_state_types,only:s_dg_hybrid_spectral_certification
  implicit none
  private
  public::select_dg_hybrid_symmetry_target,evaluate_dg_hybrid_low_energy_symmetry,&
    certify_dg_hybrid_energy_window
contains
  subroutine select_dg_hybrid_symmetry_target(comm,eigenvalues,requested_rank,tolerance,target_rank,ok,message)
    integer,intent(in)::comm,requested_rank
    real(real64),intent(in)::eigenvalues(:),tolerance
    integer,intent(out)::target_rank
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::ierr,rank,local_bad,global_bad,minimum_request,maximum_request,minimum_size,maximum_size
    real(real64),allocatable::reference(:)
    real(real64)::local_difference,global_difference,scale

    ok=.false.;message='';target_rank=0;local_bad=0
    call MPI_Allreduce(size(eigenvalues),minimum_size,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='low-energy spectrum size agreement failed';return;endif
    call MPI_Allreduce(size(eigenvalues),maximum_size,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_size/=maximum_size)then
      message='rank-disagreeing low-energy spectrum size';return
    endif
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
      occupied_rank,target_rank,tolerance,occupied_defect,target_defect,energy_defect,ok,message,worst_operation)
    integer,intent(in)::comm,occupied_rank,target_rank
    complex(real64),intent(in)::metric(:,:),representation(:,:,:),coefficients(:,:)
    real(real64),intent(in)::eigenvalues(:),tolerance
    real(real64),intent(out)::occupied_defect,target_defect,energy_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,intent(out),optional::worst_operation
    complex(real64),allocatable::target_coefficients(:,:),occupied_coefficients(:,:),target_gram(:,:),&
      occupied_gram(:,:),target_action(:,:),occupied_action(:,:),target_identity(:,:),occupied_identity(:,:),&
      target_energy(:,:),transformed_energy(:,:)
    integer::n,nstate,neigen,noperation,operation,i,ierr,local_bad,global_bad,operation_maximum
    integer::minimum_integer,maximum_integer
    integer(int64)::minimum_bits,maximum_bits
    real(real64)::orthogonality_defect,global_orthogonality_defect,global_values(3),target_scale,&
      occupied_scale,energy_scale
    real(real64),allocatable::local_operation_defects(:,:),global_operation_defects(:,:)
    real(real64)::operation_score,maximum_operation_score

    occupied_defect=huge(1d0);target_defect=huge(1d0);energy_defect=huge(1d0)
    if(present(worst_operation))worst_operation=0
    ok=.false.;message='';local_bad=0;n=size(metric,1);nstate=size(coefficients,2)
    neigen=size(eigenvalues);noperation=size(representation,3)
    call agree_integer_value(comm,n,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='rank-disagreeing low-energy basis size';return
    endif
    call agree_integer_value(comm,nstate,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='rank-disagreeing low-energy state count';return
    endif
    call agree_integer_value(comm,neigen,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='rank-disagreeing low-energy eigenvalue count';return
    endif
    call agree_integer_value(comm,noperation,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='rank-disagreeing low-energy operation count';return
    endif
    call agree_integer_value(comm,occupied_rank,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='rank-disagreeing low-energy occupied rank';return
    endif
    call agree_integer_value(comm,target_rank,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='rank-disagreeing low-energy target rank';return
    endif
    call agree_real_value(comm,tolerance,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='rank-disagreeing low-energy symmetry tolerance';return
    endif
    do i=1,neigen
      call agree_real_value(comm,eigenvalues(i),minimum_bits,maximum_bits,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
        message='rank-disagreeing low-energy eigenvalues';return
      endif
    enddo
    if(n<1.or.size(metric,2)/=n.or.any(shape(representation)/=[n,n,noperation]).or.noperation<1.or.&
      size(coefficients,1)/=n.or.neigen/=nstate.or.occupied_rank<1.or.&
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
      target_energy(target_rank,target_rank),transformed_energy(target_rank,target_rank),&
      local_operation_defects(noperation,3),global_operation_defects(noperation,3))
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
    local_bad=merge(0,1,ieee_is_finite(orthogonality_defect))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='low-energy symmetry orthogonality reduction failed';return
    endif
    call MPI_Allreduce(orthogonality_defect,global_orthogonality_defect,1,&
      MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='low-energy symmetry orthogonality reduction failed';return;endif
    orthogonality_defect=global_orthogonality_defect
    if(orthogonality_defect>tolerance)then
      message='low-energy symmetry coefficients are not metric orthonormal';return
    endif
    local_operation_defects=0d0
    target_scale=sqrt(real(target_rank,real64));occupied_scale=sqrt(real(occupied_rank,real64))
    energy_scale=max(1d0,frobenius(target_energy))
    do operation=1,noperation
      target_action=matmul(conjg(transpose(target_coefficients)),&
        matmul(metric,matmul(representation(:,:,operation),target_coefficients)))
      occupied_action=matmul(conjg(transpose(occupied_coefficients)),&
        matmul(metric,matmul(representation(:,:,operation),occupied_coefficients)))
      transformed_energy=matmul(conjg(transpose(target_action)),matmul(target_energy,target_action))
      local_operation_defects(operation,1)=frobenius(matmul(conjg(transpose(occupied_action)),occupied_action)-&
        occupied_identity)/occupied_scale
      local_operation_defects(operation,2)=frobenius(matmul(conjg(transpose(target_action)),target_action)-&
        target_identity)/target_scale
      local_operation_defects(operation,3)=frobenius(transformed_energy-target_energy)/energy_scale
    enddo
    call MPI_Allreduce(local_operation_defects,global_operation_defects,3*noperation,&
      MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(.not.ieee_is_finite(global_operation_defects)))then
      message='low-energy symmetry defect reduction failed';return
    endif
    global_values=[maxval(global_operation_defects(:,1)),maxval(global_operation_defects(:,2)),&
      maxval(global_operation_defects(:,3))]
    operation_maximum=1;maximum_operation_score=-1d0
    do operation=1,noperation
      operation_score=maxval(global_operation_defects(operation,:))
      if(operation_score>maximum_operation_score)then
        maximum_operation_score=operation_score;operation_maximum=operation
      endif
    enddo
    if(present(worst_operation))worst_operation=operation_maximum
    occupied_defect=global_values(1);target_defect=global_values(2);energy_defect=global_values(3)
    ok=maxval(global_values)<=tolerance
    if(ok)then;message='';else;message='low-energy LCFO eigenspace is not symmetry closed';endif
  end subroutine evaluate_dg_hybrid_low_energy_symmetry

  subroutine precompute_dg_hybrid_symmetry_prefix_defects(comm,metric,representation,coefficients,&
      eigenvalues,noccupied,orthogonality_defects,occupied_defect,target_defects,energy_defects,&
      worst_operations,worst_defects,ok,message)
    integer,intent(in)::comm,noccupied
    complex(real64),intent(in)::metric(:,:),representation(:,:,:),coefficients(:,:)
    real(real64),intent(in)::eigenvalues(:)
    real(real64),allocatable,intent(out)::orthogonality_defects(:),target_defects(:),energy_defects(:),&
      worst_defects(:)
    real(real64),intent(out)::occupied_defect
    integer,allocatable,intent(out)::worst_operations(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::metric_coefficients(:,:),metric_gram(:,:),transformed(:,:),&
      metric_transformed(:,:),action(:,:),closure_gram(:,:),transformed_energy(:,:)
    real(real64),allocatable::local_orthogonality(:),local_operation_defects(:,:),&
      global_operation_defects(:,:)
    complex(real64)::closure_difference,energy_difference
    integer::n,noperation,operation,boundary,i,j,ierr,local_bad,global_bad
    real(real64)::orthogonality_squared,closure_squared,energy_squared,energy_norm_squared,&
      operation_score

    ok=.false.;message='';occupied_defect=huge(1d0)
    n=size(eigenvalues);noperation=size(representation,3)
    if(n<1.or.noccupied<1.or.noccupied>n.or.size(metric,1)/=n.or.size(metric,2)/=n.or.&
      size(coefficients,1)/=n.or.size(coefficients,2)/=n.or.noperation<1.or.&
      any(shape(representation)/=[n,n,noperation]))then
      message='invalid symmetry prefix precomputation contract';return
    endif
    allocate(orthogonality_defects(n),target_defects(n),energy_defects(n),worst_defects(n),&
      worst_operations(n),local_orthogonality(n),local_operation_defects(n,2),&
      global_operation_defects(n,2),metric_coefficients(n,n),metric_gram(n,n),transformed(n,n),&
      metric_transformed(n,n),action(n,n),closure_gram(n,n),transformed_energy(n,n))

    metric_coefficients=matmul(metric,coefficients)
    metric_gram=matmul(conjg(transpose(coefficients)),metric_coefficients)
    orthogonality_squared=0d0
    do boundary=1,n
      closure_difference=metric_gram(boundary,boundary)-(1d0,0d0)
      orthogonality_squared=orthogonality_squared+abs(closure_difference)**2
      do i=1,boundary-1
        orthogonality_squared=orthogonality_squared+abs(metric_gram(i,boundary))**2+&
          abs(metric_gram(boundary,i))**2
      enddo
      local_orthogonality(boundary)=sqrt(orthogonality_squared/real(boundary,real64))
    enddo
    local_bad=merge(0,1,all(ieee_is_finite(local_orthogonality)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='symmetry prefix orthogonality reduction failed';return
    endif
    call MPI_Allreduce(local_orthogonality,orthogonality_defects,n,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='symmetry prefix orthogonality reduction failed';return;endif

    occupied_defect=-1d0;target_defects=-1d0;energy_defects=-1d0
    worst_operations=0;worst_defects=-1d0
    do operation=1,noperation
      transformed=matmul(representation(:,:,operation),coefficients)
      metric_transformed=matmul(metric,transformed)
      action=matmul(conjg(transpose(coefficients)),metric_transformed)
      closure_gram=(0d0,0d0);transformed_energy=(0d0,0d0)
      energy_norm_squared=0d0
      do boundary=1,n
        if(boundary>1)then
          do j=1,boundary-1
            do i=1,boundary-1
              closure_gram(i,j)=closure_gram(i,j)+conjg(action(boundary,i))*action(boundary,j)
              transformed_energy(i,j)=transformed_energy(i,j)+&
                conjg(action(boundary,i))*eigenvalues(boundary)*action(boundary,j)
            enddo
          enddo
          do i=1,boundary-1
            closure_gram(i,boundary)=sum(conjg(action(1:boundary,i))*action(1:boundary,boundary))
            closure_gram(boundary,i)=conjg(closure_gram(i,boundary))
            transformed_energy(i,boundary)=sum(conjg(action(1:boundary,i))*&
              eigenvalues(1:boundary)*action(1:boundary,boundary))
            transformed_energy(boundary,i)=conjg(transformed_energy(i,boundary))
          enddo
        endif
        closure_gram(boundary,boundary)=sum(abs(action(1:boundary,boundary))**2)
        transformed_energy(boundary,boundary)=sum(conjg(action(1:boundary,boundary))*&
          eigenvalues(1:boundary)*action(1:boundary,boundary))
        closure_squared=0d0;energy_squared=0d0
        do j=1,boundary
          do i=1,boundary
            closure_difference=closure_gram(i,j)
            energy_difference=transformed_energy(i,j)
            if(i==j)then
              closure_difference=closure_difference-(1d0,0d0)
              energy_difference=energy_difference-cmplx(eigenvalues(i),0d0,kind=real64)
            endif
            closure_squared=closure_squared+abs(closure_difference)**2
            energy_squared=energy_squared+abs(energy_difference)**2
          enddo
        enddo
        energy_norm_squared=energy_norm_squared+eigenvalues(boundary)**2
        local_operation_defects(boundary,1)=sqrt(closure_squared/real(boundary,real64))
        local_operation_defects(boundary,2)=sqrt(energy_squared)/max(1d0,sqrt(energy_norm_squared))
      enddo
      local_bad=merge(0,1,all(ieee_is_finite(local_operation_defects)))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        message='symmetry prefix defect reduction failed';return
      endif
      call MPI_Allreduce(local_operation_defects,global_operation_defects,2*n,&
        MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='symmetry prefix defect reduction failed';return;endif
      occupied_defect=max(occupied_defect,global_operation_defects(noccupied,1))
      do boundary=1,n
        target_defects(boundary)=max(target_defects(boundary),global_operation_defects(boundary,1))
        energy_defects(boundary)=max(energy_defects(boundary),global_operation_defects(boundary,2))
        operation_score=max(global_operation_defects(noccupied,1),&
          global_operation_defects(boundary,1),global_operation_defects(boundary,2))
        if(operation_score>worst_defects(boundary))then
          worst_defects(boundary)=operation_score;worst_operations(boundary)=operation
        endif
      enddo
    enddo
    ok=.true.
  end subroutine precompute_dg_hybrid_symmetry_prefix_defects

  subroutine certify_dg_hybrid_energy_window(comm,metric,representation,coefficients,eigenvalues,&
      noccupied,e_homo,energy_window,dynamic_requested_rank,final_orbital_tolerance,symmetry_tolerance,&
      occupied_projector_defect,density_defect,certification,ok,message)
    integer,intent(in)::comm,noccupied,dynamic_requested_rank
    complex(real64),intent(in)::metric(:,:),representation(:,:,:),coefficients(:,:)
    real(real64),intent(in)::eigenvalues(:),e_homo,energy_window,final_orbital_tolerance,&
      symmetry_tolerance,occupied_projector_defect,density_defect
    type(s_dg_hybrid_spectral_certification),intent(out)::certification
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::n,nstate,noperation,i,candidate,requested,ierr,rank,local_bad,global_bad
    integer::minimum_integer,maximum_integer
    integer(int64)::minimum_bits,maximum_bits,bits,hash_value
    real(real64)::requested_cutoff,cluster_tolerance,occupied_defect
    real(real64),allocatable::orthogonality_defects(:),target_defects(:),energy_defects(:),worst_defects(:)
    integer,allocatable::worst_operations(:)
    logical::precompute_ok,selection_ok
    character(256)::precompute_message,selection_message

    certification=s_dg_hybrid_spectral_certification();ok=.false.;message=''
    n=size(metric,1);nstate=size(eigenvalues);noperation=size(representation,3)
    call agree_integer_value(comm,nstate,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='rank-disagreeing complete LCFO spectrum size';return
    endif
    call agree_integer_value(comm,n,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='rank-disagreeing construction-basis size';return
    endif
    call agree_integer_value(comm,noperation,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='rank-disagreeing symmetry operation count';return
    endif
    call agree_integer_value(comm,noccupied,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='rank-disagreeing occupied rank';return
    endif
    call agree_integer_value(comm,dynamic_requested_rank,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='rank-disagreeing dynamic requested rank';return
    endif
    call agree_real_value(comm,e_homo,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='rank-disagreeing HOMO energy';return
    endif
    call agree_real_value(comm,energy_window,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='rank-disagreeing Hybrid energy window';return
    endif
    call agree_real_value(comm,final_orbital_tolerance,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='rank-disagreeing final orbital tolerance';return
    endif
    call agree_real_value(comm,symmetry_tolerance,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='rank-disagreeing symmetry tolerance';return
    endif
    call agree_real_value(comm,occupied_projector_defect,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='rank-disagreeing occupied-projector defect';return
    endif
    call agree_real_value(comm,density_defect,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='rank-disagreeing density defect';return
    endif
    do i=1,nstate
      call agree_real_value(comm,eigenvalues(i),minimum_bits,maximum_bits,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
        message='rank-disagreeing complete LCFO spectrum';return
      endif
    enddo

    local_bad=0
    if(n<1.or.nstate/=n.or.noperation<1.or.size(metric,2)/=n.or.&
      any(shape(representation)/=[n,n,noperation]).or.size(coefficients,1)/=n.or.&
      size(coefficients,2)/=nstate)local_bad=1
    if(noccupied<1.or.noccupied>nstate.or.dynamic_requested_rank<1.or.&
      dynamic_requested_rank>nstate)local_bad=1
    if(.not.ieee_is_finite(e_homo).or..not.ieee_is_finite(energy_window).or.&
      .not.ieee_is_finite(final_orbital_tolerance).or.final_orbital_tolerance<=0d0.or.&
      .not.ieee_is_finite(symmetry_tolerance).or.symmetry_tolerance<=0d0.or.&
      .not.ieee_is_finite(occupied_projector_defect).or.occupied_projector_defect<0d0.or.&
      .not.ieee_is_finite(density_defect).or.density_defect<0d0)local_bad=1
    if(energy_window<0d0.and.energy_window/=-1d0)local_bad=1
    if(.not.all(ieee_is_finite(eigenvalues)).or..not.finite_matrix(metric).or.&
      .not.finite_matrix(coefficients).or..not.finite_tensor(representation))local_bad=1
    if(nstate>1)then;if(any(eigenvalues(2:)<eigenvalues(:nstate-1)))local_bad=1;endif
    if(noccupied>=1.and.noccupied<=nstate)then
      if(e_homo/=eigenvalues(noccupied))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid Hybrid energy-window certification contract';return
    endif

    certification%full_rank=nstate;certification%occupied_rank=noccupied
    certification%energy_window=energy_window;certification%e_homo=e_homo
    certification%occupied_projector_defect=occupied_projector_defect
    certification%density_defect=density_defect
    if(occupied_projector_defect>symmetry_tolerance)then
      message='occupied-projector symmetry failure cannot be repaired by adding empty states';return
    endif
    if(density_defect>symmetry_tolerance)then
      message='density symmetry failure cannot be repaired by adding empty states';return
    endif

    if(energy_window==-1d0)then
      certification%compatibility_dynamic_rank=.true.
      call MPI_Comm_rank(comm,rank,ierr)
      if(ierr/=MPI_SUCCESS)then;message='compatibility warning rank query failed';return;endif
      if(rank==0)write(0,'(a)')&
        '[HYBRID-SYMMETRY-COMPATIBILITY-WARNING] energy_window=-1 uses dynamic requested-rank selection'
      requested=dynamic_requested_rank
      call select_dg_hybrid_symmetry_target(comm,eigenvalues,requested,symmetry_tolerance,&
        candidate,selection_ok,selection_message)
      if(.not.selection_ok)then
        if(index(selection_message,'unresolved degenerate cluster')==0)then
          message=trim(selection_message);return
        endif
        call select_dg_hybrid_symmetry_target(comm,eigenvalues,nstate,symmetry_tolerance,&
          candidate,selection_ok,selection_message)
        if(.not.selection_ok)then;message=trim(selection_message);return;endif
      endif
      certification%requested_rank=requested
      certification%requested_cutoff=eigenvalues(requested)
    else
      requested_cutoff=e_homo+energy_window
      if(.not.ieee_is_finite(requested_cutoff))then
        message='Hybrid requested energy cutoff is not finite';return
      endif
      requested=noccupied
      do i=noccupied+1,nstate
        cluster_tolerance=degeneracy_tolerance(eigenvalues(i-1),eigenvalues(i),requested_cutoff)
        if(eigenvalues(i)>requested_cutoff+cluster_tolerance)exit
        requested=i
      enddo
      do while(requested<nstate)
        cluster_tolerance=degeneracy_tolerance(eigenvalues(requested),eigenvalues(requested+1),requested_cutoff)
        if(eigenvalues(requested+1)-eigenvalues(requested)>cluster_tolerance)exit
        requested=requested+1
      enddo
      certification%requested_rank=requested;certification%requested_cutoff=requested_cutoff
      if(requested==nstate)then
        message='Hybrid basis capacity has no proof state above the requested cluster';return
      endif
    endif
    call precompute_dg_hybrid_symmetry_prefix_defects(comm,metric,representation,coefficients,eigenvalues,&
      noccupied,orthogonality_defects,occupied_defect,target_defects,energy_defects,worst_operations,&
      worst_defects,precompute_ok,precompute_message)
    if(.not.precompute_ok)then;message=trim(precompute_message);return;endif
    certification%occupied_subspace_defect=occupied_defect
    if(occupied_defect>symmetry_tolerance)then
      message='occupied LCFO symmetry failure cannot be repaired by adding empty states';return
    endif

    if(energy_window==-1d0)then
      call record_evaluation(candidate)
      if(max(orthogonality_defects(noccupied),orthogonality_defects(candidate))>symmetry_tolerance)then
        message='low-energy symmetry coefficients are not metric orthonormal';return
      endif
      if(target_defects(candidate)>symmetry_tolerance.or.energy_defects(candidate)>symmetry_tolerance)then
        message='dynamic-rank compatibility target is not symmetry closed';return
      endif
      call publish_certification(candidate)
      return
    endif

    candidate=requested
    do
      call record_evaluation(candidate)
      if(max(orthogonality_defects(noccupied),orthogonality_defects(candidate))>symmetry_tolerance)then
        message='low-energy symmetry coefficients are not metric orthonormal';return
      endif
      if(target_defects(candidate)<=symmetry_tolerance.and.energy_defects(candidate)<=symmetry_tolerance)then
        call publish_certification(candidate);return
      endif
      candidate=candidate+1
      do while(candidate<nstate)
        cluster_tolerance=degeneracy_tolerance(eigenvalues(candidate),eigenvalues(candidate+1),requested_cutoff)
        if(eigenvalues(candidate+1)-eigenvalues(candidate)>cluster_tolerance)exit
        candidate=candidate+1
      enddo
      if(candidate>=nstate)then
        message='Hybrid basis capacity has no passing complete cluster with a proof state';return
      endif
    enddo
  contains
    real(real64) function degeneracy_tolerance(left_energy,right_energy,cutoff)
      real(real64),intent(in)::left_energy,right_energy,cutoff
      degeneracy_tolerance=max(final_orbital_tolerance,64d0*epsilon(1d0))*&
        max(1d0,abs(left_energy),abs(right_energy),abs(cutoff))
    end function degeneracy_tolerance

    subroutine record_evaluation(boundary_rank)
      integer,intent(in)::boundary_rank
      certification%boundary_cluster_rank=boundary_rank
      certification%occupied_subspace_defect=occupied_defect
      certification%target_subspace_defect=target_defects(boundary_rank)
      certification%target_energy_defect=energy_defects(boundary_rank)
      certification%worst_operation=worst_operations(boundary_rank)
      certification%worst_operation_defect=worst_defects(boundary_rank)
      certification%maximum_physical_defect=max(worst_defects(boundary_rank),occupied_projector_defect,density_defect)
    end subroutine record_evaluation

    subroutine publish_certification(certified_rank)
      integer,intent(in)::certified_rank
      certification%certified_rank=certified_rank
      certification%certified_cutoff=eigenvalues(certified_rank)
      certification%extension_states=certified_rank-certification%requested_rank
      certification%extension_energy=max(0d0,certification%certified_cutoff-certification%requested_cutoff)
      certification%proof_state_present=certified_rank<nstate
      if(certification%proof_state_present)certification%proof_energy=eigenvalues(certified_rank+1)
      hash_value=1469598103934665603_int64
      hash_value=ieor(ishftc(hash_value,7),int(certification%full_rank,int64))
      hash_value=ieor(ishftc(hash_value,7),int(certification%occupied_rank,int64))
      hash_value=ieor(ishftc(hash_value,7),int(certification%requested_rank,int64))
      hash_value=ieor(ishftc(hash_value,7),int(certification%boundary_cluster_rank,int64))
      hash_value=ieor(ishftc(hash_value,7),int(certification%certified_rank,int64))
      hash_value=ieor(ishftc(hash_value,7),int(certification%extension_states,int64))
      hash_value=ieor(ishftc(hash_value,7),int(certification%worst_operation,int64))
      bits=transfer(certification%energy_window,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      bits=transfer(certification%e_homo,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      bits=transfer(certification%requested_cutoff,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      bits=transfer(certification%certified_cutoff,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      bits=transfer(certification%extension_energy,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      bits=transfer(certification%proof_energy,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      bits=transfer(certification%occupied_subspace_defect,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      bits=transfer(certification%occupied_projector_defect,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      bits=transfer(certification%target_subspace_defect,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      bits=transfer(certification%target_energy_defect,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      bits=transfer(certification%density_defect,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      bits=transfer(certification%worst_operation_defect,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      bits=transfer(certification%maximum_physical_defect,bits);hash_value=ieor(ishftc(hash_value,7),bits)
      if(certification%compatibility_dynamic_rank)hash_value=ieor(hash_value,97_int64)
      if(certification%proof_state_present)hash_value=ieor(hash_value,193_int64)
      if(hash_value==0_int64)hash_value=389_int64
      call MPI_Allreduce(hash_value,minimum_bits,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='spectral certification fingerprint reduction failed';return;endif
      call MPI_Allreduce(hash_value,maximum_bits,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
        message='rank-disagreeing spectral certification fingerprint';return
      endif
      certification%fingerprint=hash_value;certification%valid=.true.;ok=.true.;message=''
    end subroutine publish_certification
  end subroutine certify_dg_hybrid_energy_window

  subroutine agree_integer_value(comm,value,minimum_value,maximum_value,status)
    integer,intent(in)::comm,value
    integer,intent(out)::minimum_value,maximum_value,status
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,status)
  end subroutine agree_integer_value

  subroutine agree_real_value(comm,value,minimum_value,maximum_value,status)
    integer,intent(in)::comm
    real(real64),intent(in)::value
    integer(int64),intent(out)::minimum_value,maximum_value
    integer,intent(out)::status
    integer(int64)::value_bits
    value_bits=transfer(value,value_bits)
    call MPI_Allreduce(value_bits,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS)return
    call MPI_Allreduce(value_bits,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,status)
  end subroutine agree_real_value

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
