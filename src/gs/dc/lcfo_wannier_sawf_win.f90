module lcfo_wannier_sawf_win
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  private

  integer, parameter :: max_temp_attempts = 128

  type, public :: t_atomic_win_writer
    character(:), allocatable :: final_path, temp_path
    integer :: unit = 0
    logical :: active = .false., owns_temp = .false.
  end type t_atomic_win_writer

  public :: activate_sawf_win, deactivate_sawf_win
  public :: begin_atomic_win, finish_atomic_win, abort_atomic_win
  public :: write_sawf_local_preprocess_win

  interface
    integer(c_int) function c_rename(old_path, new_path) bind(C, name='rename')
      import :: c_char, c_int
      character(c_char), intent(in) :: old_path(*), new_path(*)
    end function c_rename
    integer(c_int) function c_remove(path) bind(C, name='remove')
      import :: c_char, c_int
      character(c_char), intent(in) :: path(*)
    end function c_remove
  end interface

contains

  subroutine write_sawf_local_preprocess_win(filename,num_bands,num_wann,lattice,atoms_fractional,ok,message)
    character(*),intent(in)::filename
    integer,intent(in)::num_bands,num_wann
    real(8),intent(in)::lattice(3,3),atoms_fractional(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(t_atomic_win_writer)::writer
    integer::unit,ios,axis,atom
    real(8)::determinant

    ok=.false.;message='';ios=0
    determinant=lattice(1,1)*(lattice(2,2)*lattice(3,3)-lattice(2,3)*lattice(3,2))- &
      lattice(1,2)*(lattice(2,1)*lattice(3,3)-lattice(2,3)*lattice(3,1))+ &
      lattice(1,3)*(lattice(2,1)*lattice(3,2)-lattice(2,2)*lattice(3,1))
    if(num_bands<=0.or.num_wann<=0.or.num_wann>num_bands.or.size(atoms_fractional,1)/=3.or. &
        size(atoms_fractional,2)<=0.or..not.all(ieee_is_finite(lattice)).or. &
        .not.all(ieee_is_finite(atoms_fractional)).or.abs(determinant)<=1d-14)then
      message='SAWF local preprocess WIN inputs are invalid';return
    end if
    call begin_atomic_win(writer,filename,unit,ok,message)
    if(.not.ok)return
    write(unit,'(a,i0)',iostat=ios)'num_bands = ',num_bands
    if(ios==0)write(unit,'(a,i0)',iostat=ios)'num_wann = ',num_wann
    if(ios==0)write(unit,'(a)',iostat=ios)'num_iter = 0'
    if(ios==0)write(unit,'(a)',iostat=ios)'mp_grid = 1 1 1'
    if(ios==0)write(unit,'(a)',iostat=ios)'gamma_only = true'
    if(ios==0)write(unit,'(a)',iostat=ios)'begin unit_cell_cart'
    if(ios==0)write(unit,'(a)',iostat=ios)'bohr'
    do axis=1,3;if(ios==0)write(unit,'(3es23.15)',iostat=ios)lattice(:,axis);end do
    if(ios==0)write(unit,'(a)',iostat=ios)'end unit_cell_cart'
    if(ios==0)write(unit,'(a)',iostat=ios)'begin atoms_frac'
    do atom=1,size(atoms_fractional,2)
      if(ios==0)write(unit,'(a,3(1x,es23.15))',iostat=ios)'X',modulo(atoms_fractional(:,atom),1d0)
    end do
    if(ios==0)write(unit,'(a)',iostat=ios)'end atoms_frac'
    if(ios==0)write(unit,'(a)',iostat=ios)'begin kpoints'
    if(ios==0)write(unit,'(3f12.6)',iostat=ios)0d0,0d0,0d0
    if(ios==0)write(unit,'(a)',iostat=ios)'end kpoints'
    if(ios/=0)then;call abort_atomic_win(writer);message='SAWF local preprocess WIN write failed';ok=.false.;return;end if
    call finish_atomic_win(writer,ok,message)
  end subroutine write_sawf_local_preprocess_win

  subroutine activate_sawf_win(win_path, tolerance, ok, message, temp_nonce)
    character(*), intent(in) :: win_path
    real(8), intent(in) :: tolerance
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, intent(in), optional :: temp_nonce

    if(.not. ieee_is_finite(tolerance) .or. tolerance <= 0.0d0) then
      ok = .false.
      message = 'SAWF win activation tolerance must be finite and positive'
      return
    end if
    call rewrite_sawf_keywords(win_path, .true., tolerance, ok, message, temp_nonce)
  end subroutine activate_sawf_win

  subroutine deactivate_sawf_win(win_path, ok, message, temp_nonce)
    character(*), intent(in) :: win_path
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, intent(in), optional :: temp_nonce

    call rewrite_sawf_keywords(win_path, .false., 1.0d0, ok, message, temp_nonce)
  end subroutine deactivate_sawf_win

  subroutine begin_atomic_win(writer, final_path, unit, ok, message, temp_nonce)
    type(t_atomic_win_writer), intent(inout) :: writer
    character(*), intent(in) :: final_path
    integer, intent(out) :: unit
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, intent(in), optional :: temp_nonce
    integer(int64) :: nonce
    integer :: attempt, io_status
    character(128) :: suffix
    character(512) :: io_message

    call abort_atomic_win(writer)
    ok = .false.; message = ''; unit = 0
    if(len_trim(final_path) == 0) then
      message = 'atomic win writer requires a nonempty final path'
      return
    end if
    nonce = make_nonce(temp_nonce)
    writer%final_path = trim(final_path)
    do attempt=0,max_temp_attempts-1
      write(suffix,'(a,i0,a,i0)') '.tmp.', nonce, '.', attempt
      writer%temp_path = trim(final_path)//trim(suffix)
      open(newunit=writer%unit, file=writer%temp_path, status='new', action='write', &
        form='formatted', iostat=io_status, iomsg=io_message)
      if(io_status == 0) then
        writer%active = .true.; writer%owns_temp = .true.
        unit = writer%unit; ok = .true.
        return
      end if
      writer%unit = 0
    end do
    message = 'atomic win writer could not exclusively create a temporary file: '//trim(io_message)
  end subroutine begin_atomic_win

  subroutine finish_atomic_win(writer, ok, message)
    type(t_atomic_win_writer), intent(inout) :: writer
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: io_status, close_status, rename_status
    character(512) :: io_message

    ok = .false.; message = ''
    if(.not. writer%active .or. writer%unit == 0 .or. .not. writer%owns_temp) then
      message = 'atomic win writer is not active'
      return
    end if
    flush(writer%unit, iostat=io_status, iomsg=io_message)
    close_status = 0
    close(writer%unit, iostat=close_status)
    writer%unit = 0
    if(io_status /= 0 .or. close_status /= 0) then
      message = 'atomic win writer flush/close failed: '//trim(io_message)
      call abort_atomic_win(writer)
      return
    end if
    rename_status = rename_path(writer%temp_path, writer%final_path)
    if(rename_status /= 0) then
      write(message,'(a,i0)') 'atomic win writer same-directory rename failed; status=', rename_status
      call abort_atomic_win(writer)
      return
    end if
    writer%active = .false.; writer%owns_temp = .false.
    ok = .true.
  end subroutine finish_atomic_win

  subroutine abort_atomic_win(writer)
    type(t_atomic_win_writer), intent(inout) :: writer
    integer :: io_status

    if(writer%unit /= 0) close(writer%unit, iostat=io_status)
    writer%unit = 0
    if(writer%owns_temp .and. allocated(writer%temp_path)) call remove_path(writer%temp_path)
    writer%active = .false.; writer%owns_temp = .false.
  end subroutine abort_atomic_win

  subroutine rewrite_sawf_keywords(win_path, append_keywords, tolerance, ok, message, temp_nonce)
    character(*), intent(in) :: win_path
    logical, intent(in) :: append_keywords
    real(8), intent(in) :: tolerance
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, intent(in), optional :: temp_nonce
    character(c_char), allocatable :: source(:), filtered(:), output(:)
    character(128) :: tolerance_line
    character(512) :: io_message
    character(:), allocatable :: temp_path
    integer(int64) :: file_size
    integer :: source_unit, temp_unit, io_status, close_status, rename_status
    integer :: source_count, filtered_count, output_count, removed_count
    logical :: exists, owns_temp, crlf

    ok = .false.; message = ''; source_unit = 0; temp_unit = 0; owns_temp = .false.
    if(len_trim(win_path) == 0) then
      message = 'SAWF win rewrite requires a nonempty path'
      return
    end if
    inquire(file=trim(win_path), exist=exists, size=file_size, iostat=io_status, iomsg=io_message)
    if(io_status /= 0) then
      message = 'SAWF win rewrite inquiry failed: '//trim(io_message)
      return
    end if
    if(.not. exists) then
      if(append_keywords) then
        message = 'SAWF win activation requires the current-run non-SAWF input file'
        return
      end if
      ok = .true.
      return
    end if
    if(file_size < 0_int64 .or. file_size > int(huge(0)-512,int64)) then
      message = 'SAWF win file is too large for byte-preserving rewrite'
      return
    end if
    source_count = int(file_size)
    allocate(source(source_count), filtered(source_count), stat=io_status)
    if(io_status /= 0) then
      message = 'SAWF win byte-buffer allocation failed'
      return
    end if
    open(newunit=source_unit, file=trim(win_path), status='old', action='read', &
      access='stream', form='unformatted', iostat=io_status, iomsg=io_message)
    if(io_status == 0 .and. source_count > 0) read(source_unit, iostat=io_status, iomsg=io_message) source
    close_status = 0
    if(source_unit /= 0) close(source_unit, iostat=close_status)
    source_unit = 0
    if(io_status /= 0 .or. close_status /= 0) then
      message = 'SAWF win byte-preserving read failed: '//trim(io_message)
      return
    end if

    call filter_keyword_lines(source, filtered, filtered_count, removed_count, crlf)
    if(.not. append_keywords .and. removed_count == 0) then
      ok = .true.
      return
    end if
    allocate(output(source_count+512), stat=io_status)
    if(io_status /= 0) then
      message = 'SAWF win output-buffer allocation failed'
      return
    end if
    output_count = 0
    if(filtered_count > 0) call append_bytes(output, output_count, filtered(1:filtered_count))
    if(append_keywords) then
      if(output_count > 0 .and. output(output_count) /= achar(10,kind=c_char)) &
        call append_newline(output, output_count, crlf)
      call append_text(output, output_count, 'site_symmetry = true')
      call append_newline(output, output_count, crlf)
      write(tolerance_line,'(a,1x,es23.15)') 'symmetrize_eps =', tolerance
      call append_text(output, output_count, trim(tolerance_line))
      call append_newline(output, output_count, crlf)
    end if

    call open_byte_temp(win_path, temp_nonce, temp_unit, temp_path, owns_temp, ok, message)
    if(.not. ok) return
    io_status = 0
    if(output_count > 0) write(temp_unit, iostat=io_status, iomsg=io_message) output(1:output_count)
    if(io_status == 0) flush(temp_unit, iostat=io_status, iomsg=io_message)
    close_status = 0
    close(temp_unit, iostat=close_status)
    temp_unit = 0
    if(io_status /= 0 .or. close_status /= 0) then
      message = 'SAWF win byte-preserving temporary write failed: '//trim(io_message)
      if(owns_temp) call remove_path(temp_path)
      return
    end if
    rename_status = rename_path(temp_path, trim(win_path))
    if(rename_status /= 0) then
      write(message,'(a,i0)') 'SAWF win byte-preserving atomic rename failed; status=', rename_status
      if(owns_temp) call remove_path(temp_path)
      return
    end if
    owns_temp = .false.; ok = .true.
  end subroutine rewrite_sawf_keywords

  subroutine filter_keyword_lines(source, filtered, filtered_count, removed_count, crlf)
    character(c_char), intent(in) :: source(:)
    character(c_char), intent(out) :: filtered(:)
    integer, intent(out) :: filtered_count, removed_count
    logical, intent(out) :: crlf
    integer :: first, last, line_end, lf

    filtered_count = 0; removed_count = 0; crlf = .false.
    if(size(source) >= 1) then
      if(source(1) /= achar(10,kind=c_char)) then
        do lf=2,size(source)
          if(source(lf) == achar(10,kind=c_char)) then
            crlf = source(lf-1) == achar(13,kind=c_char)
            exit
          end if
        end do
      end if
    end if
    first = 1
    do while(first <= size(source))
      lf = 0
      do last=first,size(source)
        if(source(last) == achar(10,kind=c_char)) then
          lf = last
          exit
        end if
      end do
      if(lf == 0) then
        last = size(source)
        line_end = last
      else
        last = lf
        line_end = lf-1
        if(line_end >= first .and. source(line_end) == achar(13,kind=c_char)) line_end = line_end-1
      end if
      if(is_keyword_assignment(source, first, line_end)) then
        removed_count = removed_count+1
      else
        call append_bytes(filtered, filtered_count, source(first:last))
      end if
      first = last+1
    end do
  end subroutine filter_keyword_lines

  logical function is_keyword_assignment(bytes, first, last) result(matches)
    character(c_char), intent(in) :: bytes(:)
    integer, intent(in) :: first, last
    integer :: pos

    matches = .false.; pos = first
    do while(pos <= last .and. horizontal_space(bytes(pos)))
      pos = pos+1
    end do
    if(pos > last) return
    matches = assignment_for_key(bytes, pos, last, 'site_symmetry') .or. &
      assignment_for_key(bytes, pos, last, 'symmetrize_eps')
  end function is_keyword_assignment

  logical function assignment_for_key(bytes, first, last, key) result(matches)
    character(c_char), intent(in) :: bytes(:)
    integer, intent(in) :: first, last
    character(*), intent(in) :: key
    integer :: i, pos

    matches = .false.
    if(last-first+1 < len(key)+1) return
    do i=1,len(key)
      if(lower_byte(bytes(first+i-1)) /= key(i:i)) return
    end do
    pos = first+len(key)
    do while(pos <= last .and. horizontal_space(bytes(pos)))
      pos = pos+1
    end do
    matches = pos <= last .and. bytes(pos) == '='
  end function assignment_for_key

  logical function horizontal_space(byte) result(is_space)
    character(c_char), intent(in) :: byte
    is_space = byte == ' ' .or. byte == achar(9,kind=c_char)
  end function horizontal_space

  character function lower_byte(byte) result(lower)
    character(c_char), intent(in) :: byte
    integer :: code
    code = iachar(byte)
    lower = achar(code)
    if(code >= iachar('A') .and. code <= iachar('Z')) lower = achar(code+32)
  end function lower_byte

  subroutine open_byte_temp(final_path, temp_nonce, unit, temp_path, owns_temp, ok, message)
    character(*), intent(in) :: final_path
    integer, intent(in), optional :: temp_nonce
    integer, intent(out) :: unit
    character(:), allocatable, intent(out) :: temp_path
    logical, intent(out) :: owns_temp, ok
    character(*), intent(out) :: message
    integer(int64) :: nonce
    integer :: attempt, io_status
    character(128) :: suffix
    character(512) :: io_message

    unit = 0; owns_temp = .false.; ok = .false.; message = ''
    nonce = make_nonce(temp_nonce)
    do attempt=0,max_temp_attempts-1
      write(suffix,'(a,i0,a,i0)') '.tmp.', nonce, '.', attempt
      temp_path = trim(final_path)//trim(suffix)
      open(newunit=unit, file=temp_path, status='new', action='write', access='stream', &
        form='unformatted', iostat=io_status, iomsg=io_message)
      if(io_status == 0) then
        owns_temp = .true.; ok = .true.
        return
      end if
      unit = 0
    end do
    message = 'SAWF win rewrite could not exclusively create a temporary file: '//trim(io_message)
  end subroutine open_byte_temp

  subroutine append_bytes(destination, count, bytes)
    character(c_char), intent(inout) :: destination(:)
    integer, intent(inout) :: count
    character(c_char), intent(in) :: bytes(:)
    destination(count+1:count+size(bytes)) = bytes
    count = count+size(bytes)
  end subroutine append_bytes

  subroutine append_text(destination, count, text)
    character(c_char), intent(inout) :: destination(:)
    integer, intent(inout) :: count
    character(*), intent(in) :: text
    integer :: i
    do i=1,len(text)
      count = count+1
      destination(count) = text(i:i)
    end do
  end subroutine append_text

  subroutine append_newline(destination, count, crlf)
    character(c_char), intent(inout) :: destination(:)
    integer, intent(inout) :: count
    logical, intent(in) :: crlf
    if(crlf) then
      count = count+1; destination(count) = achar(13,kind=c_char)
    end if
    count = count+1; destination(count) = achar(10,kind=c_char)
  end subroutine append_newline

  integer(int64) function make_nonce(temp_nonce) result(nonce)
    integer, intent(in), optional :: temp_nonce
    integer :: clock_count
    call system_clock(clock_count)
    nonce = modulo(int(clock_count,int64), huge(0_int64))
    if(present(temp_nonce)) nonce = modulo(int(temp_nonce,int64), huge(0_int64))
  end function make_nonce

  integer function rename_path(old_path, new_path) result(status)
    character(*), intent(in) :: old_path, new_path
    character(c_char), allocatable :: old_c(:), new_c(:)
    call make_c_string(old_path, old_c); call make_c_string(new_path, new_c)
    status = int(c_rename(old_c, new_c)); deallocate(old_c, new_c)
  end function rename_path

  subroutine remove_path(path)
    character(*), intent(in) :: path
    character(c_char), allocatable :: path_c(:)
    integer(c_int) :: ignored
    call make_c_string(path, path_c); ignored = c_remove(path_c); deallocate(path_c)
  end subroutine remove_path

  subroutine make_c_string(text, c_text)
    character(*), intent(in) :: text
    character(c_char), allocatable, intent(out) :: c_text(:)
    integer :: i, n
    n = len_trim(text); allocate(c_text(n+1))
    do i=1,n; c_text(i) = text(i:i); end do
    c_text(n+1) = c_null_char
  end subroutine make_c_string

end module lcfo_wannier_sawf_win
