module dg_hybrid_fragment_solver
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_get_halting_mode,ieee_set_halting_mode,ieee_set_flag,&
    ieee_invalid,ieee_divide_by_zero,ieee_overflow
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  implicit none
  private
  abstract interface
    subroutine fragment_apply(input,output,ok)
      import real64
      complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:);logical,intent(out)::ok
    end subroutine fragment_apply
  end interface
  public::solve_dg_hybrid_fragment_basis,solve_dg_hybrid_fragment_spectrum,&
    reconstruct_dg_hybrid_fragment_density,measure_dg_hybrid_fragment_core_norms
contains
  subroutine solve_dg_hybrid_fragment_basis(comm,basis,nstate,occupations,core_mask,point_weights,apply_h,apply_s,tolerance,&
      coefficients,eigenvalues,core_density,core_electron_count,maximum_residual,orthogonality_defect,&
      workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,nstate
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    real(real64),intent(in)::occupations(:),point_weights(:),tolerance
    logical,intent(in)::core_mask(:)
    procedure(fragment_apply)::apply_h,apply_s
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    real(real64),intent(out)::eigenvalues(:),core_density(:),core_electron_count,maximum_residual,orthogonality_defect
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64),allocatable::core_norms(:)
    integer::allocation_status,allocation_bad,ierr

    ok=.false.;message=''
    allocate(core_norms(max(0,nstate)),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then
      message='fragment compatibility wrapper allocation failed';return
    endif
    call solve_dg_hybrid_fragment_spectrum(comm,basis,nstate,core_mask,point_weights,apply_h,apply_s,tolerance,&
      coefficients,eigenvalues,core_norms,maximum_residual,orthogonality_defect,workspace_peak_bytes,&
      fingerprint,ok,message)
    if(.not.ok)return
    call reconstruct_dg_hybrid_fragment_density(comm,basis,coefficients,occupations,2d0,core_mask,point_weights,&
      core_density,core_electron_count,ok,message)
  end subroutine solve_dg_hybrid_fragment_basis

  subroutine solve_dg_hybrid_fragment_spectrum(comm,basis,nstate,core_mask,point_weights,apply_h,apply_s,tolerance,&
      coefficients,eigenvalues,core_norms,maximum_residual,orthogonality_defect,workspace_peak_bytes,&
      fingerprint,ok,message)
    integer,intent(in)::comm,nstate
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    logical,intent(in)::core_mask(:)
    real(real64),intent(in)::point_weights(:),tolerance
    procedure(fragment_apply)::apply_h,apply_s
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    real(real64),intent(out)::eigenvalues(:),core_norms(:),maximum_residual,orthogonality_defect
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,k,nowned,nbasis,npoint,nproc,ierr,local_bad,global_bad,info,lwork,position,rank,&
      allocation_status,allocation_bad
    integer::nstate_min,nstate_max
    integer,allocatable::point_order(:)
    integer(int64),allocatable::ordered_ids(:)
    complex(real64),allocatable::vectors(:,:),hvectors(:,:),svectors(:,:),hmat(:,:),smat(:,:),hcopy(:,:),scopy(:,:),&
      wavefunctions(:,:),work(:)
    real(real64),allocatable::all_eigenvalues(:),rwork(:),reference_norms(:)
    complex(real64)::query(1),value
    real(real64)::scale,tolerance_min,tolerance_max,local_difference,maximum_difference,agreement_limit
    integer(int64)::bits,matrix_element_count,point_basis_count,solver_workspace_bytes,layout_workspace_bytes
    logical::callback_ok,halt_invalid,halt_zero,halt_overflow,halting_disabled
    external::zhegv

    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64;halting_disabled=.false.
    maximum_residual=huge(1d0);orthogonality_defect=huge(1d0)
    if(size(eigenvalues)>0)eigenvalues=0d0
    if(size(core_norms)>0)core_norms=0d0
    local_bad=0
    if(nstate<1.or.size(eigenvalues)/=nstate.or.size(core_norms)/=nstate)local_bad=1
    if(.not.ieee_is_finite(tolerance))then
      local_bad=1
    elseif(tolerance<1d-15.or.tolerance>1d-2)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid fragment spectrum controls or output shape';return
    endif
    call MPI_Allreduce(nstate,nstate_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment spectrum state minimum reduction failed';return;endif
    call MPI_Allreduce(nstate,nstate_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment spectrum state maximum reduction failed';return;endif
    call MPI_Allreduce(tolerance,tolerance_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment spectrum tolerance minimum reduction failed';return;endif
    call MPI_Allreduce(tolerance,tolerance_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.nstate_min/=nstate_max.or.tolerance_min/=tolerance_max)then
      message='fragment spectrum controls differ between ranks';return
    endif
    call prepare_fragment_layout(comm,basis,core_mask,point_weights,ordered_ids,point_order,vectors,ok,message)
    if(.not.ok)return
    ok=.false.;nowned=size(basis%global_ids);npoint=size(vectors,1);nbasis=size(ordered_ids)
    if(nstate>nbasis)then;message='fragment spectrum state count exceeds basis extent';return;endif
    matrix_element_count=int(nbasis,int64)*int(nbasis,int64)
    point_basis_count=int(npoint,int64)*int(nbasis,int64)
    if(matrix_element_count>int(huge(0),int64))then
      message='fragment projected matrices exceed MPI extent';return
    endif

    allocate(hvectors(npoint,nbasis),svectors(npoint,nbasis),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then;message='fragment operator workspace allocation failed';return;endif
    call apply_h(vectors,hvectors,callback_ok)
    call collective_callback_status(comm,callback_ok,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='fragment Hamiltonian application failed';return
    endif
    local_bad=merge(0,1,fragment_values_are_finite(hvectors))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='fragment Hamiltonian application returned non-finite values';return
    endif
    call apply_s(vectors,svectors,callback_ok)
    call collective_callback_status(comm,callback_ok,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='fragment metric application failed';return;endif
    local_bad=merge(0,1,fragment_values_are_finite(svectors))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='fragment metric application returned non-finite values';return
    endif
    do j=1,nbasis
      hvectors(:,j)=point_weights*hvectors(:,j)
      svectors(:,j)=point_weights*svectors(:,j)
    enddo
    allocate(hmat(nbasis,nbasis),smat(nbasis,nbasis),hcopy(nbasis,nbasis),scopy(nbasis,nbasis),&
      stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then
      message='fragment projected matrix allocation failed';return
    endif
    hmat=matmul(conjg(transpose(vectors)),hvectors)
    smat=matmul(conjg(transpose(vectors)),svectors)
    local_bad=merge(0,1,fragment_values_are_finite(hmat).and.fragment_values_are_finite(smat))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='fragment projected matrices contain non-finite values';return
    endif

    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment communicator rank query failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment communicator size query failed';return;endif
    if(rank==0)then
      hcopy=hmat;scopy=smat
    endif
    call MPI_Bcast(hcopy,int(matrix_element_count),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment Hamiltonian matrix broadcast failed';return;endif
    call MPI_Bcast(scopy,int(matrix_element_count),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment projected matrix broadcast failed';return;endif
    local_difference=max(maxval(abs(hmat-hcopy)),maxval(abs(smat-scopy)))
    call MPI_Allreduce(local_difference,maximum_difference,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    scale=max(1d0,maxval(abs(hcopy)),maxval(abs(scopy)))
    agreement_limit=max(100d0*tolerance,&
      max(4096d0,64d0*real(npoint,real64))*epsilon(1d0))*scale
    if(ierr/=MPI_SUCCESS)then
      message='fragment projected matrix comparison failed';return
    elseif(.not.ieee_is_finite(maximum_difference))then
      message='fragment projected matrix comparison is non-finite';return
    elseif(maximum_difference>agreement_limit)then
      message='fragment projected matrices differ between ranks';return
    endif
    hmat=hcopy;smat=scopy
    allocate(all_eigenvalues(nbasis),rwork(max(1,3*nbasis-2)),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then;message='fragment eigensolver workspace allocation failed';return;endif
    info=0;lwork=1
    if(rank==0)then
      call ieee_get_halting_mode(ieee_invalid,halt_invalid)
      call ieee_get_halting_mode(ieee_divide_by_zero,halt_zero)
      call ieee_get_halting_mode(ieee_overflow,halt_overflow)
      call ieee_set_halting_mode(ieee_invalid,.false.)
      call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
      call ieee_set_halting_mode(ieee_overflow,.false.)
      halting_disabled=.true.;lwork=-1
      call zhegv(1,'V','U',nbasis,hmat,nbasis,smat,nbasis,all_eigenvalues,query,lwork,rwork,info)
      if(info==0)then
        if(.not.ieee_is_finite(real(query(1),real64)))then
          info=-1000
        elseif(real(query(1),real64)<1d0.or.real(query(1),real64)>real(huge(0),real64))then
          info=-1000
        else
          lwork=max(1,int(real(query(1),real64)))
        endif
      endif
      call restore_spectrum_halting_modes()
    endif
    call MPI_Bcast(info,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.info/=0)then;message='fragment ZHEGV workspace query failed';return;endif
    call MPI_Bcast(lwork,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment ZHEGV workspace broadcast failed';return;endif
    allocate(work(merge(lwork,1,rank==0)),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then;message='fragment ZHEGV work allocation failed';return;endif
    if(rank==0)then
      hmat=hcopy;smat=scopy
      call ieee_get_halting_mode(ieee_invalid,halt_invalid)
      call ieee_get_halting_mode(ieee_divide_by_zero,halt_zero)
      call ieee_get_halting_mode(ieee_overflow,halt_overflow)
      call ieee_set_halting_mode(ieee_invalid,.false.)
      call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
      call ieee_set_halting_mode(ieee_overflow,.false.)
      halting_disabled=.true.
      call zhegv(1,'V','U',nbasis,hmat,nbasis,smat,nbasis,all_eigenvalues,work,lwork,rwork,info)
      call restore_spectrum_halting_modes()
    endif
    call MPI_Bcast(info,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.info/=0)then;message='fragment generalized eigensolve failed';return;endif
    call MPI_Bcast(all_eigenvalues,nbasis,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment eigenvalue broadcast failed';return;endif
    call MPI_Bcast(hmat,int(matrix_element_count),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment eigensystem broadcast failed';return;endif
    local_bad=0
    if(any(.not.ieee_is_finite(all_eigenvalues)).or..not.fragment_values_are_finite(hmat))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='fragment eigensystem is non-finite';return;endif

    allocate(coefficients(nowned,nstate),wavefunctions(npoint,nstate),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then;message='fragment spectrum output allocation failed';return;endif
    do i=1,nowned
      position=find_sorted_fragment_id(ordered_ids,basis%global_ids(i))
      coefficients(i,:)=hmat(position,1:nstate)
    enddo
    eigenvalues=all_eigenvalues(1:nstate)
    wavefunctions=matmul(vectors,hmat(:,1:nstate))
    local_bad=merge(0,1,fragment_values_are_finite(wavefunctions).and.&
      fragment_values_are_finite(coefficients))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='fragment eigenvectors or wavefunctions are non-finite';return
    endif
    core_norms=0d0
    do i=1,npoint
      position=point_order(i)
      if(.not.core_mask(position))cycle
      do j=1,nstate
        core_norms(j)=core_norms(j)+point_weights(position)*abs(wavefunctions(position,j))**2
      enddo
    enddo
    local_bad=merge(0,1,all(ieee_is_finite(core_norms)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='fragment core norms are non-finite';return;endif
    allocate(reference_norms(nstate),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then;message='fragment core norm receipt allocation failed';return;endif
    if(rank==0)reference_norms=core_norms
    call MPI_Bcast(reference_norms,nstate,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment core norm broadcast failed';return;endif
    local_difference=maxval(abs(core_norms-reference_norms))
    call MPI_Allreduce(local_difference,maximum_difference,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    scale=max(1d0,maxval(abs(reference_norms)))
    agreement_limit=max(100d0*tolerance,&
      max(4096d0,16d0*real(npoint,real64))*epsilon(1d0))*scale
    if(ierr/=MPI_SUCCESS)then
      message='fragment core norm comparison failed';return
    elseif(any(.not.ieee_is_finite(reference_norms)).or..not.ieee_is_finite(maximum_difference))then
      message='fragment core norm comparison is non-finite';return
    elseif(maximum_difference>agreement_limit)then
      message='fragment core norms differ between ranks';return
    endif
    core_norms=reference_norms

    maximum_residual=0d0;orthogonality_defect=0d0
    if(rank==0)then
      do j=1,nstate
        maximum_residual=max(maximum_residual,maxval(abs(matmul(hcopy,hmat(:,j))-&
          eigenvalues(j)*matmul(scopy,hmat(:,j)))))
        do k=1,nstate
          value=dot_product(hmat(:,j),matmul(scopy,hmat(:,k)))
          if(j==k)value=value-(1d0,0d0)
          orthogonality_defect=max(orthogonality_defect,abs(value))
        enddo
      enddo
    endif
    call MPI_Bcast(maximum_residual,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment residual broadcast failed';return;endif
    call MPI_Bcast(orthogonality_defect,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment eigensystem receipt broadcast failed';return;endif
    scale=max(1d0,maxval(abs(eigenvalues)))
    if(.not.ieee_is_finite(maximum_residual).or..not.ieee_is_finite(orthogonality_defect))then
      message='fragment generalized eigensystem receipt is non-finite';return
    elseif(maximum_residual>100d0*tolerance*scale.or.orthogonality_defect>100d0*tolerance)then
      message='fragment generalized eigensystem receipt failed';return
    endif
    if(rank==0)then
      solver_workspace_bytes=16_int64*(3_int64*point_basis_count+4_int64*matrix_element_count+&
        int(npoint,int64)*int(nstate,int64)+int(lwork,int64))+&
        8_int64*(int(nbasis,int64)+int(size(rwork),int64)+int(npoint,int64))+&
        4_int64*int(nbasis,int64)
      layout_workspace_bytes=32_int64*point_basis_count+68_int64*int(npoint,int64)+&
        20_int64*int(nbasis,int64)+8_int64*int(nproc,int64)
      workspace_peak_bytes=max(solver_workspace_bytes,layout_workspace_bytes)
      fingerprint=ieor(basis%provenance_fingerprint,int(z'9E3779B97F4A7C15',int64))
      do j=1,nstate
        bits=transfer(eigenvalues(j),bits)
        fingerprint=ieor(ishftc(fingerprint,11),bits)
      enddo
      if(fingerprint==0_int64)fingerprint=1223_int64
    endif
    call MPI_Bcast(workspace_peak_bytes,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment workspace receipt broadcast failed';return;endif
    call MPI_Bcast(fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.workspace_peak_bytes<=0_int64.or.fingerprint==0_int64)then
      message='fragment eigensystem metadata broadcast failed';return
    endif
    ok=.true.;message=''
  contains
    subroutine restore_spectrum_halting_modes()
      if(.not.halting_disabled)return
      call ieee_set_flag(ieee_invalid,.false.)
      call ieee_set_flag(ieee_divide_by_zero,.false.)
      call ieee_set_flag(ieee_overflow,.false.)
      call ieee_set_halting_mode(ieee_invalid,halt_invalid)
      call ieee_set_halting_mode(ieee_divide_by_zero,halt_zero)
      call ieee_set_halting_mode(ieee_overflow,halt_overflow)
      halting_disabled=.false.
    end subroutine restore_spectrum_halting_modes
  end subroutine solve_dg_hybrid_fragment_spectrum

  subroutine reconstruct_dg_hybrid_fragment_density(comm,basis,coefficients,occupations,maximum_occupation,&
      core_mask,point_weights,core_density,core_electron_count,ok,message)
    integer,intent(in)::comm
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    complex(real64),intent(in)::coefficients(:,:)
    real(real64),intent(in)::occupations(:),maximum_occupation,point_weights(:)
    logical,intent(in)::core_mask(:)
    real(real64),intent(out)::core_density(:),core_electron_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    call reconstruct_fragment_density_and_norms(comm,basis,coefficients,occupations,maximum_occupation,&
      core_mask,point_weights,core_density,core_electron_count,ok,message)
  end subroutine reconstruct_dg_hybrid_fragment_density

  subroutine measure_dg_hybrid_fragment_core_norms(comm,basis,coefficients,core_mask,point_weights,core_norms,ok,message)
    ! Current coefficient-column order, including non-eigenstates from bounded
    ! updates. Shares the density reconstruction/ownership checks; no H/S or solve.
    integer,intent(in)::comm
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    complex(real64),intent(in)::coefficients(:,:)
    logical,intent(in)::core_mask(:)
    real(real64),intent(in)::point_weights(:)
    real(real64),allocatable,intent(out)::core_norms(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64),allocatable::zero_occupations(:),density(:),working_norms(:)
    real(real64)::electron_count
    integer::stat,bad,ierr
    logical::halt_invalid,halt_zero,halt_overflow
    ok=.false.;message='cannot allocate fragment core measurement'
    allocate(zero_occupations(size(coefficients,2)),density(size(core_mask)),&
      working_norms(size(coefficients,2)),stat=stat)
    call collective_allocation_status(comm,stat,bad,ierr)
    if(ierr/=MPI_SUCCESS.or.bad/=0)return
    zero_occupations=0d0
    call ieee_get_halting_mode(ieee_invalid,halt_invalid)
    call ieee_get_halting_mode(ieee_divide_by_zero,halt_zero)
    call ieee_get_halting_mode(ieee_overflow,halt_overflow)
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.)
    call reconstruct_fragment_density_and_norms(comm,basis,coefficients,zero_occupations,1d0,&
      core_mask,point_weights,density,electron_count,ok,message,working_norms)
    call ieee_set_flag(ieee_invalid,.false.)
    call ieee_set_flag(ieee_divide_by_zero,.false.)
    call ieee_set_flag(ieee_overflow,.false.)
    call ieee_set_halting_mode(ieee_invalid,halt_invalid)
    call ieee_set_halting_mode(ieee_divide_by_zero,halt_zero)
    call ieee_set_halting_mode(ieee_overflow,halt_overflow)
    if(ok)call move_alloc(working_norms,core_norms)
  end subroutine measure_dg_hybrid_fragment_core_norms

  subroutine reconstruct_fragment_density_and_norms(comm,basis,coefficients,occupations,maximum_occupation,&
      core_mask,point_weights,core_density,core_electron_count,ok,message,state_core_norms)
    integer,intent(in)::comm
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    complex(real64),intent(in)::coefficients(:,:)
    real(real64),intent(in)::occupations(:),maximum_occupation,point_weights(:)
    logical,intent(in)::core_mask(:)
    real(real64),intent(out)::core_density(:),core_electron_count
    real(real64),optional,intent(out)::state_core_norms(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,nstate,nstate_min,nstate_max,nowned,nbasis,npoint,position,rank,ierr,local_bad,global_bad,&
      allocation_status,allocation_bad
    integer::norm_controls(2),norm_min(2),norm_max(2)
    integer,allocatable::point_order(:)
    integer(int64),allocatable::ordered_ids(:)
    integer(int64)::coefficient_element_count
    complex(real64),allocatable::vectors(:,:),full_coefficients(:,:),wavefunctions(:,:)
    real(real64),allocatable::occupation_min(:),occupation_max(:),effective_occupations(:),core_norms(:),&
      reference_norms(:)
    real(real64)::maximum_min,maximum_max,effective_maximum,local_difference,maximum_difference,scale,&
      agreement_limit,reference_count,expected_count

    ok=.false.;message='';core_electron_count=0d0
    norm_controls=[merge(1,0,present(state_core_norms)),size(coefficients,2)]
    if(present(state_core_norms))then
      state_core_norms=0d0;norm_controls(2)=size(state_core_norms)
    endif
    call MPI_Allreduce(norm_controls,norm_min,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment core measurement agreement failed';return;endif
    call MPI_Allreduce(norm_controls,norm_max,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(norm_min/=norm_max))then
      message='fragment density/core measurement phase or extent differs between ranks';return
    endif
    local_bad=0
    if(.not.allocated(basis%global_ids).or..not.allocated(basis%buffer_point_ids).or.&
        .not.allocated(basis%buffer_values))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='fragment density basis is not allocated';return;endif
    nowned=size(basis%global_ids);npoint=size(basis%buffer_values,1);nstate=size(occupations)
    call MPI_Allreduce(nstate,nstate_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment density state minimum reduction failed';return;endif
    call MPI_Allreduce(nstate,nstate_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment density state maximum reduction failed';return;endif
    local_bad=0
    if(nstate<1.or.nstate_min/=nstate_max.or.size(coefficients,1)/=nowned.or.&
        size(coefficients,2)/=nstate.or.size(core_density)/=npoint)local_bad=1
    if(present(state_core_norms))then
      if(size(state_core_norms)/=nstate)local_bad=1
    endif
    if(.not.ieee_is_finite(maximum_occupation))then
      local_bad=1
    elseif(maximum_occupation<=0d0)then
      local_bad=1
    endif
    if(any(.not.ieee_is_finite(occupations)))then
      local_bad=1
    elseif(any(occupations<0d0).or.any(occupations>maximum_occupation))then
      local_bad=1
    endif
    if(.not.fragment_values_are_finite(coefficients))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid fragment density controls, occupations, or coefficient shape';return
    endif
    core_density=0d0
    call MPI_Allreduce(maximum_occupation,maximum_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment maximum occupation minimum reduction failed';return;endif
    call MPI_Allreduce(maximum_occupation,maximum_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    scale=max(1d0,abs(maximum_min),abs(maximum_max))
    agreement_limit=4096d0*epsilon(1d0)*scale
    if(ierr/=MPI_SUCCESS.or.abs(maximum_max-maximum_min)>agreement_limit)then
      message='fragment maximum occupation differs between ranks';return
    endif
    allocate(occupation_min(nstate),occupation_max(nstate),effective_occupations(nstate),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then;message='fragment occupation workspace allocation failed';return;endif
    call MPI_Allreduce(occupations,occupation_min,nstate,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation minimum reduction failed';return;endif
    call MPI_Allreduce(occupations,occupation_max,nstate,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    scale=max(1d0,maxval(abs(occupation_min)),maxval(abs(occupation_max)))
    agreement_limit=4096d0*epsilon(1d0)*scale
    if(ierr/=MPI_SUCCESS.or.maxval(abs(occupation_max-occupation_min))>agreement_limit)then
      message='fragment occupations differ between ranks';return
    endif
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment communicator rank query failed';return;endif
    effective_occupations=occupations;effective_maximum=maximum_occupation
    call MPI_Bcast(effective_occupations,nstate,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment occupation vector broadcast failed';return;endif
    call MPI_Bcast(effective_maximum,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(effective_occupations<0d0).or.&
        any(effective_occupations>effective_maximum))then
      message='fragment occupation broadcast failed';return
    endif

    call prepare_fragment_layout(comm,basis,core_mask,point_weights,ordered_ids,point_order,vectors,ok,message)
    if(.not.ok)return
    ok=.false.;nbasis=size(ordered_ids)
    if(nstate>nbasis)then;message='fragment density state count exceeds basis extent';return;endif
    coefficient_element_count=int(nbasis,int64)*int(nstate,int64)
    if(coefficient_element_count>int(huge(0),int64))then
      message='fragment coefficient collection exceeds MPI extent';return
    endif
    allocate(full_coefficients(nbasis,nstate),source=(0d0,0d0),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then
      message='fragment coefficient workspace allocation failed';return
    endif
    do i=1,nowned
      position=find_sorted_fragment_id(ordered_ids,basis%global_ids(i))
      full_coefficients(position,:)=coefficients(i,:)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,full_coefficients,int(coefficient_element_count),&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment coefficient data exchange failed';return;endif
    local_bad=merge(0,1,fragment_values_are_finite(full_coefficients))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='fragment coefficient collection failed';return;endif
    allocate(wavefunctions(npoint,nstate),core_norms(nstate),reference_norms(nstate),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then;message='fragment density workspace allocation failed';return;endif
    wavefunctions=matmul(vectors,full_coefficients)
    local_bad=merge(0,1,fragment_values_are_finite(wavefunctions))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='fragment wavefunctions are non-finite';return;endif
    do j=1,nstate
      core_density=core_density+effective_occupations(j)*abs(wavefunctions(:,j))**2
    enddo
    where(.not.core_mask)core_density=0d0
    core_norms=0d0;core_electron_count=0d0
    do i=1,npoint
      position=point_order(i)
      if(.not.core_mask(position))cycle
      core_electron_count=core_electron_count+point_weights(position)*core_density(position)
      do j=1,nstate
        core_norms(j)=core_norms(j)+point_weights(position)*abs(wavefunctions(position,j))**2
      enddo
    enddo
    local_bad=0
    if(any(.not.ieee_is_finite(core_density)).or.any(.not.ieee_is_finite(core_norms)).or.&
        .not.ieee_is_finite(core_electron_count))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='fragment density is non-finite';return;endif
    if(rank==0)then
      reference_norms=core_norms;reference_count=core_electron_count
    endif
    call MPI_Bcast(reference_norms,nstate,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment density core norm broadcast failed';return;endif
    call MPI_Bcast(reference_count,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment density electron count broadcast failed';return;endif
    local_difference=max(maxval(abs(core_norms-reference_norms)),abs(core_electron_count-reference_count))
    call MPI_Allreduce(local_difference,maximum_difference,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    scale=max(1d0,maxval(abs(reference_norms)),abs(reference_count))
    agreement_limit=max(4096d0,16d0*real(npoint,real64)*real(nstate,real64))*epsilon(1d0)*scale
    expected_count=dot_product(effective_occupations,reference_norms)
    if(ierr/=MPI_SUCCESS)then
      message='fragment core density comparison failed';return
    elseif(.not.ieee_is_finite(maximum_difference).or..not.ieee_is_finite(expected_count))then
      message='fragment core density comparison is non-finite';return
    elseif(maximum_difference>agreement_limit.or.abs(reference_count-expected_count)>agreement_limit)then
      message='fragment core density receipt differs between ranks';return
    endif
    core_electron_count=reference_count
    if(present(state_core_norms))state_core_norms=reference_norms
    ok=.true.;message=''
  end subroutine reconstruct_fragment_density_and_norms

  pure logical function fragment_values_are_finite(values)
    complex(real64),intent(in)::values(:,:)
    fragment_values_are_finite=all(ieee_is_finite(real(values,real64))).and.&
      all(ieee_is_finite(aimag(values)))
  end function fragment_values_are_finite

  subroutine sort_fragment_ids_with_order(ids,order)
    integer(int64),intent(inout)::ids(:)
    integer,intent(out)::order(:)
    integer::i,last
    integer(int64)::id_buffer
    if(size(order)/=size(ids))error stop 'fragment ID sort order extent mismatch'
    order=[(i,i=1,size(ids))]
    do i=size(ids)/2,1,-1
      call sift_down(i,size(ids))
    enddo
    do last=size(ids),2,-1
      id_buffer=ids(1);ids(1)=ids(last);ids(last)=id_buffer
      i=order(1);order(1)=order(last);order(last)=i
      call sift_down(1,last-1)
    enddo
  contains
    subroutine sift_down(first,finish)
      integer,intent(in)::first,finish
      integer::root,child,order_buffer
      integer(int64)::value_buffer
      root=first
      do while(2*root<=finish)
        child=2*root
        if(child<finish)then
          if(ids(child)<ids(child+1))child=child+1
        endif
        if(ids(root)>=ids(child))exit
        value_buffer=ids(root);ids(root)=ids(child);ids(child)=value_buffer
        order_buffer=order(root);order(root)=order(child);order(child)=order_buffer
        root=child
      enddo
    end subroutine sift_down
  end subroutine sort_fragment_ids_with_order

  pure integer function find_sorted_fragment_id(ids,target)
    integer(int64),intent(in)::ids(:),target
    integer::left,right,middle
    left=1;right=size(ids);find_sorted_fragment_id=0
    do while(left<=right)
      middle=left+(right-left)/2
      if(ids(middle)<target)then
        left=middle+1
      elseif(ids(middle)>target)then
        right=middle-1
      else
        find_sorted_fragment_id=middle;return
      endif
    enddo
  end function find_sorted_fragment_id

  subroutine collective_allocation_status(comm,allocation_status,bad,status)
    integer,intent(in)::comm,allocation_status
    integer,intent(out)::bad,status
    integer::local_bad
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,bad,1,MPI_INTEGER,MPI_MAX,comm,status)
  end subroutine collective_allocation_status

  subroutine collective_callback_status(comm,local_ok,bad,status)
    integer,intent(in)::comm
    logical,intent(in)::local_ok
    integer,intent(out)::bad,status
    integer::local_value
    local_value=merge(0,1,local_ok)
    call MPI_Allreduce(local_value,bad,1,MPI_INTEGER,MPI_MAX,comm,status)
  end subroutine collective_callback_status

  subroutine prepare_fragment_layout(comm,basis,core_mask,point_weights,ordered_ids,point_order,&
      vectors,ok,message)
    integer,intent(in)::comm
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    logical,intent(in)::core_mask(:)
    real(real64),intent(in)::point_weights(:)
    integer(int64),allocatable,intent(out)::ordered_ids(:)
    integer,allocatable,intent(out)::point_order(:)
    complex(real64),allocatable,intent(out)::vectors(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,p,q,nowned,nbasis,npoint,npoint_min,npoint_max,nproc,ierr,local_bad,global_bad,position,&
      allocation_status,allocation_bad
    integer::metadata(2),metadata_min(2),metadata_max(2)
    integer,allocatable::counts(:),displacements(:),sorted_masks(:),mask_min(:),mask_max(:),point_position(:),&
      basis_order(:)
    integer(int64),allocatable::all_ids(:),ordered_point_ids(:),point_id_min(:),point_id_max(:)
    integer(int64)::fingerprint_min,fingerprint_max,element_count
    complex(real64),allocatable::canonical_vectors(:,:)
    real(real64),allocatable::sorted_weights(:),weight_min(:),weight_max(:)
    real(real64)::scale,agreement_limit

    ok=.false.;message='';local_bad=0
    if(.not.allocated(basis%global_ids).or..not.allocated(basis%buffer_point_ids).or.&
        .not.allocated(basis%buffer_values))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='fragment solver basis is not allocated';return
    endif
    nowned=size(basis%global_ids);npoint=size(basis%buffer_values,1)
    local_bad=0
    if(npoint<1.or.size(basis%buffer_values,2)/=nowned.or.size(basis%buffer_point_ids)/=npoint.or.&
        size(core_mask)/=npoint.or.size(point_weights)/=npoint)local_bad=1
    call MPI_Allreduce(npoint,npoint_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment point extent minimum reduction failed';return;endif
    call MPI_Allreduce(npoint,npoint_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment point extent maximum reduction failed';return;endif
    metadata=[basis%fragment_id,basis%generation]
    call MPI_Allreduce(metadata,metadata_min,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment metadata minimum reduction failed';return;endif
    call MPI_Allreduce(metadata,metadata_max,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment metadata maximum reduction failed';return;endif
    call MPI_Allreduce(basis%provenance_fingerprint,fingerprint_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment provenance minimum reduction failed';return;endif
    call MPI_Allreduce(basis%provenance_fingerprint,fingerprint_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(npoint_min/=npoint_max.or.any(metadata_min/=metadata_max).or.fingerprint_min/=fingerprint_max)local_bad=1
    if(basis%fragment_id<=0.or.basis%generation<=0.or.basis%provenance_fingerprint==0_int64)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid or rank-inconsistent fragment basis shape';return
    endif
    local_bad=0
    if(any(basis%buffer_point_ids<=0_int64).or..not.fragment_values_are_finite(basis%buffer_values))local_bad=1
    if(any(.not.ieee_is_finite(point_weights)))then
      local_bad=1
    elseif(any(point_weights<=0d0))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid fragment point catalog, weights, or basis values';return
    endif

    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment communicator size query failed';return;endif
    allocate(counts(nproc),displacements(nproc),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then
      message='fragment ownership metadata allocation failed';return
    endif
    call MPI_Allgather(nowned,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment basis ownership count exchange failed';return;endif
    displacements(1)=0
    do i=2,nproc
      displacements(i)=displacements(i-1)+counts(i-1)
    enddo
    nbasis=sum(counts)
    if(nbasis<1)then;message='fragment basis is empty';return;endif
    allocate(all_ids(nbasis),ordered_ids(nbasis),basis_order(nbasis),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then;message='fragment basis ID allocation failed';return;endif
    call MPI_Allgatherv(basis%global_ids,nowned,MPI_INTEGER8,all_ids,counts,displacements,&
      MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment basis ownership exchange failed';return;endif
    ordered_ids=all_ids;call sort_fragment_ids_with_order(ordered_ids,basis_order)
    local_bad=merge(0,1,all(ordered_ids>0_int64))
    do i=2,nbasis
      if(ordered_ids(i)==ordered_ids(i-1))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid fragment basis ownership';return;endif

    allocate(ordered_point_ids(npoint),point_order(npoint),point_position(npoint),point_id_min(npoint),&
      point_id_max(npoint),sorted_weights(npoint),weight_min(npoint),weight_max(npoint),sorted_masks(npoint),&
      mask_min(npoint),mask_max(npoint),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then;message='fragment point metadata allocation failed';return;endif
    ordered_point_ids=basis%buffer_point_ids
    call sort_fragment_ids_with_order(ordered_point_ids,point_order)
    local_bad=0
    do i=2,npoint
      if(ordered_point_ids(i)==ordered_point_ids(i-1))local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='duplicate fragment physical point ID';return;endif
    call MPI_Allreduce(ordered_point_ids,point_id_min,npoint,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment point catalog minimum reduction failed';return;endif
    call MPI_Allreduce(ordered_point_ids,point_id_max,npoint,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(point_id_min/=point_id_max))then
      message='fragment physical point catalog differs between ranks';return
    endif
    ordered_point_ids=point_id_min
    do q=1,npoint
      p=point_order(q)
      point_position(p)=q
      sorted_weights(q)=point_weights(p)
      sorted_masks(q)=merge(1,0,core_mask(p))
    enddo
    call MPI_Allreduce(sorted_weights,weight_min,npoint,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment point weight minimum reduction failed';return;endif
    call MPI_Allreduce(sorted_weights,weight_max,npoint,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment point weight maximum reduction failed';return;endif
    call MPI_Allreduce(sorted_masks,mask_min,npoint,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment core mask minimum reduction failed';return;endif
    call MPI_Allreduce(sorted_masks,mask_max,npoint,MPI_INTEGER,MPI_MAX,comm,ierr)
    scale=max(1d0,maxval(abs(weight_min)),maxval(abs(weight_max)))
    agreement_limit=4096d0*epsilon(1d0)*scale
    if(ierr/=MPI_SUCCESS.or.maxval(abs(weight_max-weight_min))>agreement_limit.or.any(mask_min/=mask_max))then
      message='fragment point weights or core mask differ between ranks';return
    endif

    element_count=int(npoint,int64)*int(nbasis,int64)
    if(element_count>int(huge(0),int64))then
      message='fragment basis collection exceeds MPI extent';return
    endif
    allocate(canonical_vectors(npoint,nbasis),source=(0d0,0d0),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then;message='fragment basis collection allocation failed';return;endif
    do i=1,nowned
      position=find_sorted_fragment_id(ordered_ids,basis%global_ids(i))
      do p=1,npoint
        q=point_position(p)
        canonical_vectors(q,position)=basis%buffer_values(p,i)
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,canonical_vectors,int(element_count),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment basis column data exchange failed';return;endif
    local_bad=merge(0,1,fragment_values_are_finite(canonical_vectors))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='fragment basis column collection failed';return;endif
    allocate(vectors(npoint,nbasis),stat=allocation_status)
    call collective_allocation_status(comm,allocation_status,allocation_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.allocation_bad/=0)then;message='fragment basis output allocation failed';return;endif
    do q=1,npoint
      vectors(point_order(q),:)=canonical_vectors(q,:)
    enddo
    ok=.true.;message=''
  end subroutine prepare_fragment_layout
end module dg_hybrid_fragment_solver
