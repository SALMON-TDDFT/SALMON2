module sym_sub
  implicit none
contains
  subroutine read_symmetry_file(filename, symmetry, ok, message)
    character(*), intent(in) :: filename
    real(8), allocatable, intent(out) :: symmetry(:,:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message

    allocate(symmetry(3,4,0))
    ok = .false.
    message = 'focused Spglib test does not provide file symmetry: '//trim(filename)
  end subroutine read_symmetry_file
end module sym_sub
