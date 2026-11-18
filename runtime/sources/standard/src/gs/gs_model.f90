! Model-owned input boundary for sequential standard GS calls.
! Internal Fortran API: no retained orbitals, simultaneous TDDFT or thread reentry.
module salmonpy_gs_model
  use salmonpy_gs_storage, only: s_gs_storage
  use salmonpy_gs_solver, only: gs_solve
  use inputoutput, only: s_input_parser,read_input,initialize_inputoutput_units
  use salmonpy_model_inputs, only: s_model_inputs
  use code_optimization, only: s_model_optimization,park_model_optimization,activate_model_optimization
  use salmon_xc, only: s_model_xc_workspace,park_model_xc_workspace,activate_model_xc_workspace
  implicit none
  private
  public :: s_gs_model
  type :: s_gs_model
    private
    logical :: ready=.false.
    type(s_input_parser) :: parser
    type(s_model_optimization) :: optimization
    type(s_model_xc_workspace) :: xc_workspace
  contains
    procedure :: load => gs_model_load
    procedure :: run => gs_model_run
    procedure :: clear => gs_model_clear
  end type
contains
  subroutine gs_model_load(self,path,status)
    use salmon_global, only: theory,yn_periodic,spin,yn_md,yn_opt,yn_dc,yn_symmetry,yn_spinorbit,xc
    class(s_gs_model),intent(inout) :: self
    character(*),intent(in) :: path
    integer,intent(out) :: status
    type(s_input_parser) :: candidate
    call read_input(path,candidate)
    status=2
    ! Preserve the unrestricted legacy entry; this explicit boundary has a measured scope.
    if(theory/='dft'.or.yn_periodic/='y'.or.spin/='unpolarized')return
    if(yn_md/='n'.or.yn_opt/='n'.or.yn_dc/='n'.or.yn_symmetry/='n'.or.yn_spinorbit/='n')return
    select case(trim(xc))
    case('pz','tbmbj','libxc_pz','libxc_pbe')
    case default
      return
    end select
    call self%clear()
    self%parser=candidate
    self%ready=.true.
    status=0
  end subroutine

  subroutine gs_model_run(self,status,state)
    use salmon_global, only: calc_mode
    use code_optimization, only: modx,mody,modz
    use plusU_global, only: PLUS_U_ON
    use communication, only: comm_is_attached,comm_begin_children,comm_end_children,comm_release_children
    use sym_sub, only: swap_symmetry_switch
    class(s_gs_model),intent(inout) :: self
    integer,intent(out) :: status
    type(s_gs_storage),optional,target,intent(inout) :: state
    type(s_gs_storage),target :: temporary
    type(s_gs_storage),pointer :: selected
    type(s_model_inputs) :: previous
    type(s_model_optimization) :: old_optimization
    type(s_model_xc_workspace) :: old_xc
    character(len(calc_mode)) :: old_mode
    integer :: restore_status,child_owner
    logical :: had_optimization,old_symmetry,old_directions(3)
    status=1
    if(.not.self%ready)return
    if(PLUS_U_ON.or..not.comm_is_attached())then
      status=2
      return
    endif
    selected=>temporary
    if(present(state))then
      if(state%active)return
      selected=>state
    endif
    call comm_begin_children(child_owner,status)
    if(status/=0)return
    call previous%capture()
    old_mode=calc_mode
    call self%parser%values%activate(status)
    if(status/=0)then
      call comm_end_children(child_owner)
      call comm_release_children(child_owner)
      return
    endif
    call initialize_inputoutput_units
    had_optimization=allocated(modx).or.allocated(mody).or.allocated(modz)
    if(had_optimization)call park_model_optimization(old_optimization)
    call park_model_xc_workspace(old_xc)
    call activate_model_optimization(self%optimization)
    call activate_model_xc_workspace(self%xc_workspace)
    old_symmetry=.false.;old_directions=.false.
    call swap_symmetry_switch(old_symmetry,old_directions)
    calc_mode='GS'
    call gs_solve(selected,.true.,present(state))
    call comm_end_children(child_owner)
    if(present(state))then
      state%mpi_child_owner=child_owner
    else
      call comm_release_children(child_owner)
    endif
    call swap_symmetry_switch(old_symmetry,old_directions)
    call park_model_optimization(self%optimization)
    call park_model_xc_workspace(self%xc_workspace)
    call activate_model_optimization(old_optimization)
    call activate_model_xc_workspace(old_xc)
    call previous%activate(restore_status)
    call initialize_inputoutput_units
    calc_mode=old_mode
    status=restore_status
    if(present(state))then
      if(state%report_status/=0)status=3
    endif
  end subroutine

  subroutine gs_model_clear(self)
    class(s_gs_model),intent(inout) :: self
    type(s_gs_model) :: empty
    self%ready=.false.
    self%parser=empty%parser
    self%optimization=empty%optimization
    self%xc_workspace=empty%xc_workspace
  end subroutine
end module
