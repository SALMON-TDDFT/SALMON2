module lcfo_wannier_sawf_collective
  implicit none
  private

  public :: reduce_sawf_fragment_alignment_failure
  public :: validate_sawf_density_contribution,assemble_sawf_density_unique

contains

  subroutine validate_sawf_density_contribution(local_density,active_contributor,rank,ok,message)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(8), intent(in) :: local_density(:)
    logical, intent(in) :: active_contributor
    integer, intent(in) :: rank
    logical, intent(out) :: ok
    character(*), intent(out) :: message

    ok=.true.;message=''
    if(.not.active_contributor)return
    if(.not.all(ieee_is_finite(local_density)))then
      ok=.false.
      write(message,'(a,i0,a)')'rank-local SAWF density on rank ',rank,' must be finite'
    else if(any(local_density<0d0))then
      ok=.false.
      write(message,'(a,i0,a)')'rank-local SAWF density on rank ',rank,' must be nonnegative'
    end if
    if(.not.ok)write(*,'(1x,a)')'[DC-LCFO-SAWF-DENSITY] '//trim(message)
  end subroutine validate_sawf_density_contribution

  subroutine assemble_sawf_density_unique(local_density,active_contributor,communicator,global_density)
    use communication, only: comm_summation
    real(8), intent(in) :: local_density(:)
    logical, intent(in) :: active_contributor
    integer, intent(in) :: communicator
    real(8), intent(out) :: global_density(:)
    real(8), allocatable :: contribution(:)

    allocate(contribution(size(local_density)))
    contribution=0d0
    if(active_contributor)contribution=local_density
    call comm_summation(contribution,global_density,size(contribution),communicator)
    deallocate(contribution)
  end subroutine assemble_sawf_density_unique

  subroutine reduce_sawf_fragment_alignment_failure(local_failure,communicator,rank,operation, &
      grid_map_ok,fragment_map_ok,max_targets_per_source,max_grid_residual,details,global_failure)
    use communication, only: comm_get_max
    integer, intent(in) :: local_failure,communicator,rank,operation,max_targets_per_source
    logical, intent(in) :: grid_map_ok,fragment_map_ok
    real(8), intent(in) :: max_grid_residual
    character(*), intent(in) :: details
    integer, intent(out) :: global_failure

    global_failure=merge(0,1,local_failure==0)
    if(global_failure/=0) write(*,'(1x,a,i0,a,i0,2(a,l1),a,i0,a,es13.5,2a)') &
      '[DC-LCFO-SAWF-ALIGN] rank=',rank,' operation=',operation, &
      ' grid_map_ok=',grid_map_ok,' fragment_map_ok=',fragment_map_ok, &
      ' max_targets_per_source=',max_targets_per_source, &
      ' max_grid_residual=',max_grid_residual,' details=',trim(details)
    call comm_get_max(global_failure,communicator)
  end subroutine reduce_sawf_fragment_alignment_failure

end module lcfo_wannier_sawf_collective
