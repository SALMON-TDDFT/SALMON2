module nvtx
#ifdef USE_NVTX
  use iso_c_binding
  implicit none

  integer,private:: col(13) = [ &
       Z'0000ff00', Z'00007fff', Z'00ff00ff', Z'00ff7f00', &
       Z'0000ff7f', Z'000000ff', Z'00ff007f', Z'00ffff00', &
       Z'0000ffff', Z'007f00ff', Z'00ff0000', Z'007fff00', &
       Z'00000000' ]
  character(len=256),private:: tempName

  type, bind(C):: nvtxEventAttributes
    integer(C_INT16_T):: version=1
    integer(C_INT16_T):: size=48 !
    integer(C_INT):: category=0
    integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
    integer(C_INT):: color
    integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
    integer(C_INT):: reserved0
    integer(C_INT64_T):: payload   ! union uint,int,double
    integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII = 1
    type(C_PTR):: message  ! ascii char
  end type nvtxEventAttributes

  interface nvtxRangePush
    ! push range with custom label and standard color
    subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
      use iso_c_binding
      character(kind=C_CHAR,len=*) :: name
    end subroutine nvtxRangePushA

    ! push range with custom label and custom color
    subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
      use iso_c_binding
      import:: nvtxEventAttributes
      type(nvtxEventAttributes):: event
    end subroutine nvtxRangePushEx
  end interface

  interface nvtxRangePop
    subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
    end subroutine nvtxRangePop
  end interface

contains

  subroutine nvtxStartRange(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id

    integer:: r, g, b
    type(nvtxEventAttributes):: event

    tempName=trim(name)//c_null_char

    if ( .not. present(id)) then
      call nvtxRangePush(tempName)
    else
      !event%color=col(mod(id,13)+1)
      r = mod(id * 47, 256)
      g = mod(id * 79, 256)
      b = mod(id *113, 256)
      event%color = b + 256 * (g + 256 * r)
      event%message=c_loc(tempName)
      call nvtxRangePushEx(event)
    end if
  end subroutine nvtxStartRange

  subroutine nvtxEndRange
    call nvtxRangePop
  end subroutine nvtxEndRange
#else
  implicit none

contains

  subroutine nvtxStartRange(name,id)
    character(len=*) :: name
    integer, optional:: id
  end subroutine nvtxStartRange

  subroutine nvtxEndRange
  end subroutine nvtxEndRange
#endif
end module nvtx
