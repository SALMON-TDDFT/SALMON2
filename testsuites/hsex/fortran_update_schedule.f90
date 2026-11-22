program test_update_schedule
  use hsex_update_schedule, only: screen_update_due
  implicit none
  integer :: step,count
  if(.not.screen_update_due(5,5,0))error stop 'legacy event missing'
  if(screen_update_due(6,5,0))error stop 'legacy event repeated'
  if(screen_update_due(1,0,0).or.screen_update_due(1,0,1))error stop 'disabled event fired'
  if(screen_update_due(5,5,-1).or.screen_update_due(0,1,1))error stop 'invalid event fired'
  count=0
  do step=1,9
    if(screen_update_due(step,2,3))then
      if(step/=2.and.step/=5.and.step/=8)error stop 'wrong event'
      count=count+1
    endif
  enddo
  if(count/=3)error stop 'missing periodic event'
  ! Restart at5 continues at6, preserving the absolute next event8.
  count=0
  do step=6,9
    if(screen_update_due(step,2,3))then
      if(step/=8)error stop 'restart event duplicated'
      count=count+1
    endif
  enddo
  if(count/=1)error stop 'restart event missing'
  if(.not.screen_update_due(huge(step),huge(step)-3,3))error stop 'large-index event'
  print *,'HSEX_UPDATE_SCHEDULE_LEGACY_PERIODIC_RESTART_PASS'
end program
