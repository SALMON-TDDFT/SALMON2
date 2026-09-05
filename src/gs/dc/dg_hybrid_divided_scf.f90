module dg_hybrid_divided_scf
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dc_scf_convergence,only:reduce_dc_density_convergence
  implicit none
  private
  abstract interface
    subroutine update_total_potential_interface(density,ok)
      import real64
      real(real64),intent(in)::density(:);logical,intent(out)::ok
    end subroutine
    subroutine solve_fragments_interface(iteration,ok)
      integer,intent(in)::iteration;logical,intent(out)::ok
    end subroutine
    subroutine assemble_core_density_interface(density,electron_count,ok)
      import real64
      real(real64),intent(out)::density(:),electron_count;logical,intent(out)::ok
    end subroutine
    subroutine mix_dc_density_interface(iteration,input_density,new_density,mixed_density,ok)
      import real64
      integer,intent(in)::iteration
      real(real64),intent(in)::input_density(:),new_density(:)
      real(real64),intent(out)::mixed_density(:);logical,intent(out)::ok
    end subroutine
  end interface
  public::run_dg_hybrid_divided_scf
contains
  subroutine run_dg_hybrid_divided_scf(comm,global_point_count,core_ids,initial_density,&
      cell_volume,convergence_mode,threshold,update_total_potential,solve_fragments,&
      assemble_core_density,mix_dc_density,maximum_iterations,core_weights,&
      expected_electron_count,electron_tolerance,converged_density,iterations,&
      convergence_value,electron_defect,ok,message)
    integer,intent(in)::comm,global_point_count,maximum_iterations
    integer(int64),intent(in)::core_ids(:)
    real(real64),intent(in)::initial_density(:),cell_volume,threshold,core_weights(:)
    real(real64),intent(in)::expected_electron_count,electron_tolerance
    character(*),intent(in)::convergence_mode
    procedure(update_total_potential_interface)::update_total_potential
    procedure(solve_fragments_interface)::solve_fragments
    procedure(assemble_core_density_interface)::assemble_core_density
    procedure(mix_dc_density_interface)::mix_dc_density
    real(real64),allocatable,intent(out)::converged_density(:)
    integer,intent(out)::iterations
    real(real64),intent(out)::convergence_value,electron_defect
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,ierr,rank,nproc,nlocal,ntotal,mode_code,local_invalid,global_invalid
    integer::integer_controls(3),minimum_integers(3),maximum_integers(3)
    integer,allocatable::counts(:),displacements(:)
    integer(int64),allocatable::all_ids(:)
    real(real64),allocatable::density(:),new_density(:),mixed_density(:)
    real(real64)::electron_count,local_absolute_sum,local_square_sum,mixed_electron_defect,input_electron_defect
    real(real64)::real_controls(4),minimum_reals(4),maximum_reals(4)
    logical::callback_ok
    character(256)::electron_message

    ok=.false.;message='';iterations=0;convergence_value=huge(1d0);electron_defect=huge(1d0)
    select case(trim(adjustl(convergence_mode)))
    case('rho_dne');mode_code=1
    case('norm_rho');mode_code=2
    case('norm_rho_dng');mode_code=3
    case default;mode_code=0
    end select
    local_invalid=0
    if(global_point_count<=0.or.maximum_iterations<=0.or.&
      size(core_ids)/=size(initial_density).or.size(core_ids)/=size(core_weights).or.&
      mode_code==0)local_invalid=1
    if(.not.ieee_is_finite(threshold).or..not.ieee_is_finite(cell_volume).or.&
      .not.ieee_is_finite(expected_electron_count).or.&
      .not.ieee_is_finite(electron_tolerance))then
      local_invalid=1
    else if(threshold<=0d0.or.cell_volume<=0d0.or.expected_electron_count<=0d0.or.&
      electron_tolerance<=0d0)then
      local_invalid=1
    endif
    if(.not.all(ieee_is_finite(initial_density)).or.&
      .not.all(ieee_is_finite(core_weights)))then
      local_invalid=1
    else if(any(core_weights<=0d0))then
      local_invalid=1
    endif
    call MPI_Allreduce(local_invalid,global_invalid,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided SCF control validation failed';return;endif
    if(global_invalid/=0)then;message='invalid divided SCF controls';return;endif
    integer_controls=[global_point_count,maximum_iterations,mode_code]
    call MPI_Allreduce(integer_controls,minimum_integers,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided SCF integer agreement failed';return;endif
    call MPI_Allreduce(integer_controls,maximum_integers,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided SCF integer agreement failed';return;endif
    real_controls=[threshold,cell_volume,expected_electron_count,electron_tolerance]
    call MPI_Allreduce(real_controls,minimum_reals,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided SCF real agreement failed';return;endif
    call MPI_Allreduce(real_controls,maximum_reals,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided SCF real agreement failed';return;endif
    if(any(minimum_integers/=maximum_integers).or.any(minimum_reals/=maximum_reals))then
      message='rank-disagreeing divided SCF controls';return
    endif

    nlocal=size(core_ids);call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided SCF rank lookup failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr)
    allocate(counts(nproc),displacements(nproc))
    call MPI_Allgather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided SCF ownership count exchange failed';return;endif
    displacements(1)=0
    do i=2,nproc;displacements(i)=displacements(i-1)+counts(i-1);enddo
    ntotal=sum(counts);allocate(all_ids(ntotal))
    call MPI_Allgatherv(core_ids,nlocal,MPI_INTEGER8,all_ids,counts,displacements,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided SCF ownership exchange failed';return;endif
    if(ntotal/=global_point_count)then;message='divided SCF core ownership is incomplete';return;endif
    do i=1,ntotal
      if(all_ids(i)<1_int64.or.all_ids(i)>int(global_point_count,int64))then
        message='divided SCF core ID is outside the global range';return
      endif
      do j=i+1,ntotal
        if(all_ids(i)==all_ids(j))then;message='divided SCF core ownership is duplicated';return;endif
      enddo
    enddo

    call validate_density_electron_count(initial_density,0d0,.false.,'initial density',&
      electron_defect,callback_ok,electron_message)
    if(.not.callback_ok)then;message=trim(electron_message);return;endif
    input_electron_defect=electron_defect
    allocate(density(nlocal),new_density(nlocal),mixed_density(nlocal));density=initial_density
    do iterations=1,maximum_iterations
      call update_total_potential(density,callback_ok)
      if(.not.collective_success(callback_ok))then;message='divided SCF potential update failed';return;endif
      call solve_fragments(iterations,callback_ok)
      if(.not.collective_success(callback_ok))then;message='divided SCF fragment solve failed';return;endif
      call assemble_core_density(new_density,electron_count,callback_ok)
      if(.not.collective_success(callback_ok))then;message='divided SCF core density assembly failed';return;endif
      call validate_density_electron_count(new_density,electron_count,.true.,&
        'callback or independent density',electron_defect,callback_ok,electron_message)
      if(.not.callback_ok)then;message=trim(electron_message);return;endif
      ! Like the DC seed convergence path, finite charge drift delays
      ! convergence; it does not abort the early/mixed-density iteration.
      electron_defect=max(electron_defect,input_electron_defect)
      local_absolute_sum=sum(abs(new_density-density))
      local_square_sum=sum((new_density-density)**2)
      call reduce_dc_density_convergence(comm,convergence_mode,local_absolute_sum,local_square_sum,&
        cell_volume,expected_electron_count,global_point_count,convergence_value,callback_ok,message)
      if(.not.callback_ok)return
      if(rank==0.and.(iterations==1.or.mod(iterations,50)==0))write(*,'(a,i0,2(a,es12.4))')&
        '[DG-HYBRID-DIVIDED-SCF] iteration=',iterations,' convergence=',convergence_value,&
        ' electron_defect=',electron_defect
      if(convergence_value<=threshold.and.electron_defect<=electron_tolerance)then
        call update_total_potential(new_density,callback_ok)
        if(.not.collective_success(callback_ok))then
          message='divided SCF terminal potential refresh failed';return
        endif
        allocate(converged_density(nlocal),source=new_density);ok=.true.;message='';return
      endif
      call mix_dc_density(iterations,density,new_density,mixed_density,callback_ok)
      if(.not.collective_success(callback_ok))then;message='divided SCF density mixing failed';return;endif
      call validate_density_electron_count(mixed_density,0d0,.false.,'mixed density',&
        mixed_electron_defect,callback_ok,electron_message)
      if(.not.callback_ok)then
        electron_defect=mixed_electron_defect;message=trim(electron_message);return
      endif
      input_electron_defect=mixed_electron_defect
      density=mixed_density
    enddo
    write(message,'(a,i0,2(a,es12.4))')'divided SCF did not converge within maximum_iterations=',&
      maximum_iterations,' convergence=',convergence_value,' electron_defect=',electron_defect
  contains
    logical function collective_success(local_ok)
      logical,intent(in)::local_ok
      call MPI_Allreduce(local_ok,collective_success,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS)collective_success=.false.
    end function collective_success
    subroutine validate_density_electron_count(candidate_density,reported_electron_count,&
        check_reported_count,density_kind,defect,valid,detail)
      real(real64),intent(in)::candidate_density(:),reported_electron_count
      logical,intent(in)::check_reported_count
      character(*),intent(in)::density_kind
      real(real64),intent(out)::defect
      logical,intent(out)::valid
      character(*),intent(out)::detail
      integer::density_invalid,global_density_invalid,reduction_error
      real(real64)::local_integral,global_integral,reported_minimum,reported_maximum,consistency_defect

      valid=.false.;detail='';defect=huge(1d0);density_invalid=0
      if(.not.all(ieee_is_finite(candidate_density)))density_invalid=1
      if(check_reported_count.and..not.ieee_is_finite(reported_electron_count))density_invalid=1
      call MPI_Allreduce(density_invalid,global_density_invalid,1,MPI_INTEGER,MPI_MAX,comm,reduction_error)
      if(reduction_error/=MPI_SUCCESS)then
        detail='divided SCF density validation reduction failed';return
      endif
      if(global_density_invalid/=0)then
        detail='non-finite divided SCF density or callback electron count';return
      endif
      local_integral=sum(core_weights*candidate_density)
      call MPI_Allreduce(local_integral,global_integral,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,reduction_error)
      if(reduction_error/=MPI_SUCCESS)then
        detail='divided SCF independent electron integration failed';return
      endif
      if(.not.ieee_is_finite(global_integral))then
        detail='non-finite divided SCF independent electron count';return
      endif
      defect=abs(global_integral-expected_electron_count)
      if(check_reported_count)then
        call MPI_Allreduce(reported_electron_count,reported_minimum,1,MPI_DOUBLE_PRECISION,&
          MPI_MIN,comm,reduction_error)
        if(reduction_error/=MPI_SUCCESS)then
          detail='divided SCF callback electron minimum failed';return
        endif
        call MPI_Allreduce(reported_electron_count,reported_maximum,1,MPI_DOUBLE_PRECISION,&
          MPI_MAX,comm,reduction_error)
        if(reduction_error/=MPI_SUCCESS)then
          detail='divided SCF callback electron maximum failed';return
        endif
        ! Reported count and independently integrated density describe the
        ! same state: disagreement here is an error, not SCF nonconvergence.
        consistency_defect=max(abs(reported_minimum-global_integral),abs(reported_maximum-global_integral))
        if(.not.ieee_is_finite(consistency_defect).or.consistency_defect>electron_tolerance)then
          defect=max(defect,consistency_defect)
          detail='divided SCF '//trim(density_kind)//' electron count mismatch';return
        endif
      endif
      if(.not.ieee_is_finite(defect))then
        detail='non-finite divided SCF electron defect';return
      endif
      valid=.true.
    end subroutine validate_density_electron_count
  end subroutine run_dg_hybrid_divided_scf
end module dg_hybrid_divided_scf
