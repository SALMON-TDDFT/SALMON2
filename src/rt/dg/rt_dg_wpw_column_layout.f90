module rt_dg_wpw_column_layout
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  private

  character(*), parameter :: windowed_kg = 'windowed_kg'

  type, public :: s_dg_wpw_column_layout
    character(16) :: basis_kind = ''
    character(16) :: ownership_kind = ''
    integer :: n_global_columns = 0
    integer :: n_g_modes = 0
    integer :: owned_fragment_id = 0
    integer, allocatable :: pw_fragment_ids(:)
    integer, allocatable :: pw_g_ids(:)
    integer, allocatable :: pw_owner(:)
    integer, allocatable :: owned_column_ids(:)
  end type s_dg_wpw_column_layout

  public :: initialize_wpw_column_layout
  public :: initialize_wpw_fragment_root_layout
  public :: wpw_column_id
  public :: wpw_column_pair
  public :: wpw_column_owner
  public :: wpw_fragment_root_owner

contains

  subroutine initialize_wpw_column_layout(layout, n_fragments, n_g_modes, rank_id, nrank, info)
    type(s_dg_wpw_column_layout), intent(out) :: layout
    integer, intent(in) :: n_fragments, n_g_modes, rank_id, nrank
    integer, intent(out) :: info

    integer(int64) :: ncolumns64, first64, last64
    integer :: i, column_id, fragment_id, g_id, nowned

    info = 0
    if (n_fragments <= 0 .or. n_g_modes <= 0) then
      info = 1
      return
    end if
    if (nrank <= 0 .or. rank_id < 0 .or. rank_id >= nrank) then
      info = 2
      return
    end if
    ncolumns64 = int(n_fragments, int64) * int(n_g_modes, int64)
    if (ncolumns64 > int(huge(layout%n_global_columns), int64)) then
      info = 3
      return
    end if

    layout%basis_kind = windowed_kg
    layout%ownership_kind = 'arithmetic'
    layout%n_global_columns = int(ncolumns64)
    layout%n_g_modes = n_g_modes
    first64 = (int(rank_id, int64) * ncolumns64 + int(nrank, int64) - 1_int64) / int(nrank, int64) + 1_int64
    last64 = (int(rank_id + 1, int64) * ncolumns64 + int(nrank, int64) - 1_int64) / int(nrank, int64)
    nowned = max(0, int(last64 - first64 + 1_int64))
    allocate(layout%pw_fragment_ids(nowned), layout%pw_g_ids(nowned), &
             layout%pw_owner(nowned), layout%owned_column_ids(nowned))
    do i = 1, nowned
      column_id = int(first64) + i - 1
      call wpw_column_pair(column_id, n_g_modes, fragment_id, g_id, info)
      if (info /= 0) return
      layout%owned_column_ids(i) = column_id
      layout%pw_fragment_ids(i) = fragment_id
      layout%pw_g_ids(i) = g_id
      layout%pw_owner(i) = rank_id
    end do
  end subroutine initialize_wpw_column_layout

  subroutine initialize_wpw_fragment_root_layout(layout,n_fragments,n_g_modes,owned_fragment_id, &
      rank_id,nrank,info)
    type(s_dg_wpw_column_layout), intent(out) :: layout
    integer, intent(in) :: n_fragments,n_g_modes,owned_fragment_id,rank_id,nrank
    integer, intent(out) :: info
    integer(int64) :: ncolumns64,first64,last64
    integer :: i,column_id

    info=0
    if(n_fragments<=0 .or. n_g_modes<=0 .or. owned_fragment_id<1 .or. &
       owned_fragment_id>n_fragments .or. nrank/=n_fragments .or. rank_id/=owned_fragment_id-1)then
      info=1;return
    endif
    ncolumns64=int(n_fragments,int64)*int(n_g_modes,int64)
    first64=int(owned_fragment_id-1,int64)*int(n_g_modes,int64)+1_int64
    last64=first64+int(n_g_modes-1,int64)
    if(ncolumns64>int(huge(layout%n_global_columns),int64) .or. &
       last64>int(huge(layout%n_global_columns),int64))then
      info=2;return
    endif
    layout%basis_kind=windowed_kg
    layout%ownership_kind='fragment_root'
    layout%n_global_columns=int(ncolumns64)
    layout%n_g_modes=n_g_modes
    layout%owned_fragment_id=owned_fragment_id
    allocate(layout%pw_fragment_ids(n_g_modes),layout%pw_g_ids(n_g_modes), &
      layout%pw_owner(n_g_modes),layout%owned_column_ids(n_g_modes))
    do i=1,n_g_modes
      column_id=int(first64)+i-1
      layout%owned_column_ids(i)=column_id
      layout%pw_fragment_ids(i)=owned_fragment_id
      layout%pw_g_ids(i)=i
      layout%pw_owner(i)=rank_id
    enddo
  end subroutine initialize_wpw_fragment_root_layout

  integer function wpw_column_id(fragment_id, g_id, n_g_modes, info) result(column_id)
    integer, intent(in) :: fragment_id, g_id, n_g_modes
    integer, intent(out) :: info
    integer(int64) :: column64

    column_id = 0
    info = 0
    if (fragment_id <= 0 .or. n_g_modes <= 0 .or. g_id <= 0 .or. g_id > n_g_modes) then
      info = 1
      return
    end if
    column64 = int(fragment_id - 1, int64) * int(n_g_modes, int64) + int(g_id, int64)
    if (column64 > int(huge(column_id), int64)) then
      info = 2
      return
    end if
    column_id = int(column64)
  end function wpw_column_id

  subroutine wpw_column_pair(column_id, n_g_modes, fragment_id, g_id, info)
    integer, intent(in) :: column_id, n_g_modes
    integer, intent(out) :: fragment_id, g_id, info

    fragment_id = 0
    g_id = 0
    info = 0
    if (column_id <= 0 .or. n_g_modes <= 0) then
      info = 1
      return
    end if
    fragment_id = (column_id - 1) / n_g_modes + 1
    g_id = modulo(column_id - 1, n_g_modes) + 1
  end subroutine wpw_column_pair

  integer function wpw_column_owner(column_id, n_global_columns, nrank, info) result(owner)
    integer, intent(in) :: column_id, n_global_columns, nrank
    integer, intent(out) :: info

    owner = -1
    info = 0
    if (column_id <= 0 .or. column_id > n_global_columns .or. n_global_columns <= 0 .or. nrank <= 0) then
      info = 1
      return
    end if
    owner = int((int(column_id - 1, int64) * int(nrank, int64)) / int(n_global_columns, int64))
  end function wpw_column_owner

  integer function wpw_fragment_root_owner(column_id,n_g_modes,n_fragments,info) result(owner)
    integer, intent(in) :: column_id,n_g_modes,n_fragments
    integer, intent(out) :: info
    integer(int64) :: ncolumns64
    owner=-1;info=0
    if(n_g_modes<=0 .or. n_fragments<=0 .or. column_id<=0)then
      info=1;return
    endif
    ncolumns64=int(n_g_modes,int64)*int(n_fragments,int64)
    if(int(column_id,int64)>ncolumns64)then
      info=1;return
    endif
    owner=(column_id-1)/n_g_modes
  end function wpw_fragment_root_owner

end module rt_dg_wpw_column_layout
