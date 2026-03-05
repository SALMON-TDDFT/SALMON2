module paw_variables

  implicit none

  private
  public :: paw

  type s_paw_info
     real(8),allocatable :: rad(:,:)
     real(8),allocatable :: rab(:,:)    ! (dr/dx)*dx
     real(8),allocatable :: aewf(:,:,:)
     real(8),allocatable :: pswf(:,:,:)
     real(8),allocatable :: occ(:,:)
     real(8),allocatable :: rhops(:,:)
     real(8),allocatable :: rhoae(:,:)
     real(8),allocatable :: rhocps(:,:)
     real(8),allocatable :: rhocae(:,:)
     real(8),allocatable :: vion(:,:)
     real(8),allocatable :: vzero(:,:)
     real(8),allocatable :: beta(:,:,:)
     real(8),allocatable :: vlocps(:,:)
     real(8),allocatable :: vlocae(:,:)
     integer,allocatable :: orb_2_l(:,:)
     integer,allocatable :: orb_2_n(:,:)
     integer,allocatable :: nl_2_orb(:,:,:)
     real(8),allocatable :: QIJL(:,:,:,:,:)
     real(8),allocatable :: DIJ(:,:,:)
     real(8),allocatable :: Q(:,:)
     real(8),allocatable :: MULTIPOLES(:,:)
  end type s_paw_info

  type(s_paw_info) :: paw

end module paw_variables
