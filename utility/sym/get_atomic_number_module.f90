module atomic_number_module

  implicit none

  private
  public :: get_atomic_number
  public :: get_element_name

  character(2) :: name(112)
  logical :: flag_init=.true.

contains

  subroutine init
    implicit none
    name =(/ "H" ,"He", &
             "Li","Be","B" ,"C" ,"N" ,"O" ,"F" ,"Ne", &
             "Na","Mg","Al","Si","P" ,"S" ,"Cl","Ar", &
             "K" ,"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co", &
             "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", &
             "Rb","Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh", &
             "Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe", &
             "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu", &
             "Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", &
             "Hf","Ta","W" ,"Re","Os","Ir", &
             "Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", &
             "Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am", &
             "Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf", &
             "Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn" /)
    flag_init=.false.
  end subroutine init

  subroutine get_element_name( z, element_name )
    implicit none
    integer,intent(in) :: z
    character(2),intent(out) :: element_name
    if ( flag_init ) call init
    element_name=name(z)
  end subroutine get_element_name

  subroutine get_atomic_number( element_name, z )
    implicit none
    character(*),intent(in) :: element_name
    integer,intent(out) :: z
    character(2) :: a
    integer :: i
    if ( flag_init ) call init
    z=0
    do i=1,size(name)
       if ( element_name(1:2) == name(i) ) then
          z=i
          exit
       else
          a=name(i)
          if ( element_name(1:1) == a(1:len_trim(a)) ) z=i
       end if
    end do
  end subroutine get_atomic_number

end module atomic_number_module
