module lcfo_wannier_sawf_collective
  implicit none
  private

  public :: reduce_sawf_fragment_alignment_failure

contains

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
