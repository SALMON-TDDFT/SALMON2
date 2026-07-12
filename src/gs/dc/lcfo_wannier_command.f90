module lcfo_wannier_command
  implicit none
  private
  public :: select_wannier90_command, execute_wannier90_command
  public :: is_wannier90_reuse_command, is_wannier90_export_only_command
  public :: is_wannier90_import_only_command
  public :: cache_resolved_wannier90_command, get_cached_wannier90_command
  character(1024), save :: cached_command = ''
  logical, save :: command_cached = .false.

contains

  subroutine select_wannier90_command(namelist_command, environment_command, compiled_default, resolved)
    character(*), intent(in) :: namelist_command, environment_command, compiled_default
    character(*), intent(out) :: resolved

    if(len_trim(namelist_command) > 0) then
      resolved = trim(namelist_command)
    else if(len_trim(environment_command) > 0) then
      resolved = trim(environment_command)
    else
      resolved = trim(compiled_default)
    end if
  end subroutine select_wannier90_command

  subroutine execute_wannier90_command(command_line, ok, message)
    character(*), intent(in) :: command_line
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: exit_status, command_status

    ok = .false.; message = ''
    call execute_command_line(trim(command_line), exitstat=exit_status, cmdstat=command_status)
    if(command_status /= 0) then
      write(message,'(a,i0)') 'Wannier90 execute_command_line command status=', command_status
      return
    end if
    if(exit_status /= 0) then
      write(message,'(a,i0)') 'Wannier90 external command exit status=', exit_status
      return
    end if
    ok = .true.
  end subroutine execute_wannier90_command

  subroutine cache_resolved_wannier90_command(command)
    character(*), intent(in) :: command
    cached_command = trim(command)
    command_cached = .true.
  end subroutine cache_resolved_wannier90_command

  subroutine get_cached_wannier90_command(command, found)
    character(*), intent(out) :: command
    logical, intent(out) :: found
    command = cached_command
    found = command_cached
  end subroutine get_cached_wannier90_command

  logical function is_wannier90_reuse_command(command) result(is_reuse)
    character(*), intent(in) :: command
    character(len(command)) :: value
    value = lowercase(adjustl(command))
    is_reuse = trim(value) == 'skip' .or. trim(value) == 'reuse'
  end function is_wannier90_reuse_command

  logical function is_wannier90_export_only_command(command) result(export_only)
    character(*), intent(in) :: command
    character(len(command)) :: value
    value = lowercase(adjustl(command))
    export_only = trim(value) == 'export_only' .or. trim(value) == 'export-only' .or. &
      trim(value) == 'seed_only' .or. trim(value) == 'seed-only'
  end function is_wannier90_export_only_command

  logical function is_wannier90_import_only_command(command) result(import_only)
    character(*), intent(in) :: command
    character(len(command)) :: value
    value = lowercase(adjustl(command))
    import_only = trim(value) == 'import_only' .or. trim(value) == 'import-only'
  end function is_wannier90_import_only_command

  pure function lowercase(text) result(lower)
    character(*), intent(in) :: text
    character(len(text)) :: lower
    integer :: i, code
    lower = text
    do i=1,len(text)
      code = iachar(text(i:i))
      if(code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code+32)
    end do
  end function lowercase

end module lcfo_wannier_command
