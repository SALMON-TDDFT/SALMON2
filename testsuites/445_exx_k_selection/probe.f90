program probe
  use exx_k_selection, only: resolve_k_algorithm
  implicit none
  call check('auto',.false.,' ',' ',' ','nyn',0)
  call check('auto',.true.,' ',' ',' ','nnn',0)
  call check('reference',.false.,' ',' ',' ','nnn',0)
  call check('reference',.true.,' ',' ',' ','nnn',0)
  call check('pair',.false.,' ',' ',' ','nyn',0)
  call check('pair',.true.,' ',' ',' ','nny',0)
  call check('orbits',.false.,' ',' ',' ','ynn',0)
  call check('orbits',.true.,' ',' ',' ','ynn',0)
  call check('auto',.false.,'n','n','n','nnn',0)
  call check('auto',.true.,' ',' ','y','nny',0)
  call check('reference',.true.,' ',' ','n','nnn',0)
  call check('pair',.true.,' ',' ','y','nny',0)
  call check('reference',.true.,' ',' ','y','nnn',1)
  call check('pair',.true.,' ',' ','n','nnn',1)
  call check('orbits',.false.,' ','y',' ','nnn',1)
  call check('invalid',.false.,' ',' ',' ','nnn',1)
  call check('auto',.false.,'z',' ',' ','nnn',1)
  call check('auto',.true.,'y',' ','y','nnn',1)
  print *, 'EXX_K_SELECTION_PASS'
contains
  subroutine check(mode,spinor,old_orbits,old_pairs,old_spinor,expected,error)
    implicit none
    character(*),intent(in) :: mode,expected
    character(1),intent(in) :: old_orbits,old_pairs,old_spinor
    logical,intent(in) :: spinor
    integer,intent(in) :: error
    character(1) :: orbits,pairs,spin_pairs
    character(120) :: reason
    integer :: ierr
    orbits=old_orbits;pairs=old_pairs;spin_pairs=old_spinor
    call resolve_k_algorithm(mode,spinor,orbits,pairs,spin_pairs,reason,ierr)
    if(ierr/=error)then
      print *, mode,old_orbits,old_pairs,old_spinor,reason,ierr,error
      stop 1
    endif
    if(ierr/=0)return
    if(orbits//pairs//spin_pairs/=expected)then
      print *, mode,orbits//pairs//spin_pairs,expected
      stop 2
    endif
    if(len_trim(reason)==0)stop 3
  end subroutine check
end program probe
