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
      cell_volume,electron_count_target,&
      convergence_mode,threshold,update_total_potential,solve_fragments,assemble_core_density,&
      mix_dc_density,maximum_iterations,converged_density,iterations,convergence_value,ok,message)
    integer,intent(in)::comm,global_point_count,maximum_iterations
    integer(int64),intent(in)::core_ids(:)
    real(real64),intent(in)::initial_density(:),cell_volume,electron_count_target,threshold
    character(*),intent(in)::convergence_mode
    procedure(update_total_potential_interface)::update_total_potential
    procedure(solve_fragments_interface)::solve_fragments
    procedure(assemble_core_density_interface)::assemble_core_density
    procedure(mix_dc_density_interface)::mix_dc_density
    real(real64),allocatable,intent(out)::converged_density(:)
    integer,intent(out)::iterations
    real(real64),intent(out)::convergence_value
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,ierr,nproc,nlocal,ntotal,mode_code,local_invalid,global_invalid
    integer::integer_controls(3),minimum_integers(3),maximum_integers(3)
    integer,allocatable::counts(:),displacements(:)
    integer(int64),allocatable::all_ids(:)
    real(real64),allocatable::density(:),new_density(:),mixed_density(:)
    real(real64)::electron_count,local_absolute_sum,local_square_sum
    real(real64)::real_controls(3),minimum_reals(3),maximum_reals(3)
    logical::callback_ok

    ok=.false.;message='';iterations=0;convergence_value=huge(1d0)
    select case(trim(adjustl(convergence_mode)))
    case('rho_dne');mode_code=1
    case('norm_rho');mode_code=2
    case('norm_rho_dng');mode_code=3
    case default;mode_code=0
    end select
    local_invalid=0
    if(global_point_count<=0.or.maximum_iterations<=0.or.&
      size(core_ids)/=size(initial_density).or.mode_code==0)local_invalid=1
    if(.not.ieee_is_finite(threshold).or..not.ieee_is_finite(cell_volume).or.&
      .not.ieee_is_finite(electron_count_target))then
      local_invalid=1
    else if(threshold<=0d0.or.cell_volume<=0d0.or.electron_count_target<=0d0)then
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
    real_controls=[threshold,cell_volume,electron_count_target]
    call MPI_Allreduce(real_controls,minimum_reals,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided SCF real agreement failed';return;endif
    call MPI_Allreduce(real_controls,maximum_reals,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='divided SCF real agreement failed';return;endif
    if(any(minimum_integers/=maximum_integers).or.any(minimum_reals/=maximum_reals))then
      message='rank-disagreeing divided SCF controls';return
    endif

    nlocal=size(core_ids);call MPI_Comm_size(comm,nproc,ierr)
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

    allocate(density(nlocal),new_density(nlocal),mixed_density(nlocal));density=initial_density
    do iterations=1,maximum_iterations
      call update_total_potential(density,callback_ok)
      if(.not.collective_success(callback_ok))then;message='divided SCF potential update failed';return;endif
      call solve_fragments(iterations,callback_ok)
      if(.not.collective_success(callback_ok))then;message='divided SCF fragment solve failed';return;endif
      call assemble_core_density(new_density,electron_count,callback_ok)
      if(.not.collective_success(callback_ok))then;message='divided SCF core density assembly failed';return;endif
      local_absolute_sum=sum(abs(new_density-density))
      local_square_sum=sum((new_density-density)**2)
      call reduce_dc_density_convergence(comm,convergence_mode,local_absolute_sum,local_square_sum,&
        cell_volume,electron_count_target,global_point_count,convergence_value,callback_ok,message)
      if(.not.callback_ok)return
      if(convergence_value<=threshold)then
        allocate(converged_density(nlocal),source=new_density);ok=.true.;message='';return
      endif
      call mix_dc_density(iterations,density,new_density,mixed_density,callback_ok)
      if(.not.collective_success(callback_ok))then;message='divided SCF density mixing failed';return;endif
      density=mixed_density
    enddo
    message='divided SCF did not converge within maximum_iterations'
  contains
    logical function collective_success(local_ok)
      logical,intent(in)::local_ok
      call MPI_Allreduce(local_ok,collective_success,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS)collective_success=.false.
    end function collective_success
  end subroutine run_dg_hybrid_divided_scf
end module dg_hybrid_divided_scf
