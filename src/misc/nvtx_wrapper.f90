module nvtx_wrapper
  use iso_c_binding
  implicit none

#ifdef USE_NVTX
  use nvtx, only nvtxStartRange => nvtxStartRange_org, nvtxEndRange => nvtxEndRange_org
#endif

contains

  subroutine nvtxStartRange(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id

#ifdef USE_NVTX
    call nvtxStartRange_org(name, id)
#endif
  end subroutine nvtxStartRange

  subroutine nvtxEndRange
#ifdef USE_NVTX
    call nvtxEndRange_org
#endif
  end subroutine nvtxEndRange
end module nvtx_wrapper
