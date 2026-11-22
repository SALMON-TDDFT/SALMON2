! Absolute RT step scheduling, shared by projection and screening updates.
module hsex_update_schedule
  implicit none
  private
  public :: screen_update_due
contains
  pure logical function screen_update_due(step,first,interval) result(due)
    implicit none
    integer,intent(in) :: step,first,interval
    due=.false.
    if(step<1.or.first<1.or.interval<0)return
    if(step<first)return
    if(interval==0)then
      due=step==first
    else
      due=mod(step-first,interval)==0
    endif
  end function
end module
