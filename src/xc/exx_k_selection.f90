! Resolve native k-exchange policy without introducing numerical approximations.
module exx_k_selection
  implicit none
  private
  public :: resolve_k_algorithm
contains
  subroutine resolve_k_algorithm(mode,spinor,orbits,pairs,spin_pairs,reason,ierr)
    implicit none
    character(*),intent(in) :: mode
    logical,intent(in) :: spinor
    character(1),intent(inout) :: orbits,pairs,spin_pairs
    character(*),intent(out) :: reason
    integer,intent(out) :: ierr
    character(1) :: selected(3),legacy(3)
    integer :: i
    logical :: explicit_legacy
    ierr=1;reason='Invalid exx_k_algorithm or legacy yes/no value'
    legacy=[orbits,pairs,spin_pairs]
    do i=1,3
      if(legacy(i)/=' '.and.legacy(i)/='y'.and.legacy(i)/='n')return
    enddo
    explicit_legacy=any(legacy/=' ')
    selected=['n','n','n']
    select case(mode)
    case('auto')
      if(.not.spinor)selected(2)='y'
      reason='auto: eligible scalar pairs; reference otherwise'
      if(spinor)reason='auto: conservative spinor reference; pair performance not certified'
      do i=1,3
        if(legacy(i)/=' ')selected(i)=legacy(i)
      enddo
      if(explicit_legacy)reason='legacy switches: compatibility override of auto policy'
    case('reference')
      reason='manual reference: k-exchange reuse disabled'
    case('pair')
      if(spinor)then
        selected(3)='y'
      else
        selected(2)='y'
      endif
      reason='manual pair: spinor full pairs or scalar Hermitian pairs'
    case('orbits')
      selected(1)='y'
      reason='manual orbits: crystal-symmetry FFT reuse'
    case default
      return
    end select
    if(mode/='auto')then
      do i=1,3
        if(legacy(i)/=' '.and.legacy(i)/=selected(i))then
          reason='Conflicting exx_k_algorithm and legacy exchange switches'
          return
        endif
      enddo
    endif
    if(selected(1)=='y'.and.selected(3)=='y')then
      reason='Conflicting orbit reuse and spinor full-pair selection'
      return
    endif
    orbits=selected(1);pairs=selected(2);spin_pairs=selected(3)
    ierr=0
  end subroutine resolve_k_algorithm
end module exx_k_selection
