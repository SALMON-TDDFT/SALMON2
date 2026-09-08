module structures
  implicit none
  type :: s_dcdft
    logical :: optimized_fragment_geometry = .false.
    integer :: nxyz_domain(3) = 0
    integer, allocatable :: nxyz_domain_frag(:,:)
  end type s_dcdft
end module structures
