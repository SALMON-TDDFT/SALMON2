program test_dg_wpw_fragment_support
  use dg_wpw_fragment_support, only: s_dg_wpw_support_record, build_dg_wpw_fragment_support,&
    dg_wpw_support_overlap_box
  implicit none
  type(s_dg_wpw_support_record), allocatable :: records(:)
  integer, allocatable :: unique_fragments(:)
  integer :: info, i, count_fragment_two
  integer :: overlap_lo(3),overlap_hi(3)

  call build_dg_wpw_fragment_support(14, [3,3,3], [1,1,1], records, unique_fragments, info)
  call require(info == 0, '3D support construction failed')
  call require(size(records) == 26, '3D support is not the bounded 26-neighbor stencil')
  call require(size(unique_fragments) == 27, '3D deduplicated support fragments are incomplete')
  call require(all([(all(abs(records(i)%image_shift) <= 1), i=1,size(records))]), &
    'support image shift exceeds one-neighbor stencil')
  call require(all(record_periodic_shift(records,[1,0,0])==[0,0,0]),&
    'interior neighbor was incorrectly marked as a periodic image')

  call build_dg_wpw_fragment_support(1, [2,1,1], [1,1,1], records, unique_fragments, info)
  call require(info == 0, 'two-fragment periodic support construction failed')
  call require(size(records) == 2, 'periodic minus and plus overlap records were collapsed')
  count_fragment_two = count([(records(i)%fragment_id == 2, i=1,size(records))])
  call require(count_fragment_two == 2, 'two periodic images do not map to the neighbor fragment')
  call require(any([(records(i)%image_shift(1) == -1, i=1,size(records))]), 'minus image missing')
  call require(any([(records(i)%image_shift(1) == 1, i=1,size(records))]), 'plus image missing')
  call require(all(record_periodic_shift(records,[-1,0,0])==[-1,0,0]),&
    'wrapped minus neighbor has wrong periodic shift')
  call require(all(record_periodic_shift(records,[1,0,0])==[0,0,0]),&
    'unwrapped plus neighbor has wrong periodic shift')
  call require(size(unique_fragments) == 2, 'normalization fragment list double-counts periodic images')
  call require(all(unique_fragments == [1,2]), 'normalization fragment IDs are not deterministic')
  call dg_wpw_support_overlap_box([1,1,1],[4,4,4],[1,0,0],[2,2,2],overlap_lo,overlap_hi,info)
  call require(info==0.and.all(overlap_lo==[3,1,1]).and.all(overlap_hi==[4,4,4]),&
    'plus-x neighbor buffered overlap box is wrong')
  call dg_wpw_support_overlap_box([1,1,1],[4,4,4],[-1,-1,0],[2,2,2],overlap_lo,overlap_hi,info)
  call require(info==0.and.all(overlap_lo==[1,1,1]).and.all(overlap_hi==[2,2,4]),&
    'minus edge-neighbor buffered overlap box is wrong')
  call dg_wpw_support_overlap_box([1,1,1],[4,4,4],[1,0,0],[5,2,2],overlap_lo,overlap_hi,info)
  call require(info/=0,'buffer wider than one fragment did not fail closed')

  call build_dg_wpw_fragment_support(1, [2,3,4], [1,1,1], records, unique_fragments, info)
  call require(info == 0, 'anisotropic fragment support construction failed')
  call require(record_fragment(records,[1,0,0]) == 13, 'fragment ID is not SALMON z-fast for +x')
  call require(record_fragment(records,[0,1,0]) == 5, 'fragment ID is not SALMON z-fast for +y')
  call require(record_fragment(records,[0,0,1]) == 2, 'fragment ID is not SALMON z-fast for +z')

  call build_dg_wpw_fragment_support(1, [2,1,1], [2,1,1], records, unique_fragments, info)
  call require(info /= 0, 'support beyond one neighbor did not fail closed')

  print *, 'PASS bounded fragment window-support schedule'

contains
  function record_periodic_shift(entries,shift)result(periodic)
    type(s_dg_wpw_support_record),intent(in)::entries(:)
    integer,intent(in)::shift(3)
    integer::periodic(3)
    integer::k
    periodic=99
    do k=1,size(entries)
      if(all(entries(k)%image_shift==shift))then;periodic=entries(k)%periodic_shift;return;endif
    enddo
  end function record_periodic_shift
  integer function record_fragment(entries,shift)result(id)
    type(s_dg_wpw_support_record),intent(in)::entries(:)
    integer,intent(in)::shift(3)
    integer::k
    id=0
    do k=1,size(entries)
      if(all(entries(k)%image_shift==shift))then;id=entries(k)%fragment_id;return;endif
    enddo
  end function
  subroutine require(condition,message)
    logical,intent(in)::condition
    character(*),intent(in)::message
    if(.not.condition)then
      write(*,'(a)')'FAIL: '//trim(message)
      error stop 1
    endif
  end subroutine
end program test_dg_wpw_fragment_support
