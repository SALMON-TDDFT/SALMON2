module dg_hybrid_divided_scf
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
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
      convergence_mode,threshold,update_total_potential,solve_fragments,assemble_core_density,&
      mix_dc_density,maximum_iterations,converged_density,iterations,convergence_value,ok,message)
    integer,intent(in)::comm,global_point_count,maximum_iterations
    integer(int64),intent(in)::core_ids(:)
    real(real64),intent(in)::initial_density(:),threshold
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
    integer::i,j,ierr,nproc,nlocal,ntotal
    integer,allocatable::counts(:),displacements(:)
    integer(int64),allocatable::all_ids(:)
    real(real64),allocatable::density(:),new_density(:),mixed_density(:)
    real(real64)::electron_count,local_values(3),global_values(3)
    logical::callback_ok

    ok=.false.;message='';iterations=0;convergence_value=huge(1d0)
    if(global_point_count<=0.or.maximum_iterations<=0.or.threshold<=0d0)then
      message='invalid divided SCF controls';return
    endif
    if(size(core_ids)/=size(initial_density))then
      message='divided SCF core density shape mismatch';return
    endif
    select case(trim(adjustl(convergence_mode)))
    case('rho_dne','norm_rho','norm_rho_dng')
    case default
      message='unsupported divided SCF convergence quantity';return
    end select

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
      local_values=[maxval(abs(new_density-density)),sum((new_density-density)**2),sum(new_density**2)]
      call MPI_Allreduce(local_values,global_values,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      call MPI_Allreduce(MPI_IN_PLACE,global_values(2:3),2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='divided SCF convergence reduction failed';return;endif
      select case(trim(adjustl(convergence_mode)))
      case('rho_dne');convergence_value=global_values(1)
      case('norm_rho');convergence_value=sqrt(global_values(2)/real(global_point_count,real64))
      case('norm_rho_dng');convergence_value=sqrt(global_values(2)/max(global_values(3),tiny(1d0)))
      end select
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
