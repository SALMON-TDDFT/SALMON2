module sym_vector_sub

  use sym_sub, only: Amat, Ainv, SymMatA, SymMatB, use_symmetry

  implicit none

  private
  public :: sym_vector_xyz
  public :: sym_vector_force_xyz

contains

  subroutine sym_vector_xyz( vxyz )
    implicit none
    real(8),intent(inout) :: vxyz(3)
    real(8) :: vtmp(3), Rvtmp(3)
    integer :: isym, nsym
    if ( .not.use_symmetry ) return
    nsym = size(SymMatB,3)
    vtmp = matmul( Ainv, vxyz )
    vxyz = 0.0d0
    do isym=1,nsym
       Rvtmp = matmul( SymMatB(:,1:3,isym), vtmp )
       vxyz = vxyz + Rvtmp
    end do
    vtmp = vxyz/nsym
    vxyz = matmul( Amat, vtmp )
  end subroutine sym_vector_xyz


  subroutine sym_vector_force_xyz( Fxyz, Rxyz )
    implicit none
    real(8),intent(inout) :: Fxyz(:,:)
    real(8),intent(in)  :: Rxyz(:,:)
    real(8) :: RotRaaa(3), RotFaaa(3)
    real(8) :: diff
    integer :: isym, nsym, iatm, natm, jatm
    real(8),allocatable :: Raaa(:,:), Faaa(:,:), Ftmp(:,:)

    if ( .not.use_symmetry ) return

    nsym = size(SymMatA,3)
    natm = size(Fxyz,2)
    allocate( Ftmp(3,natm) ); Ftmp=0.0d0
    allocate( Raaa(3,natm) ); Raaa=0.0d0
    allocate( Faaa(3,natm) ); Faaa=0.0d0

    Raaa(:,:) = matmul( Ainv, Rxyz )
    do iatm=1,natm
      call shift_in_01( Raaa(:,iatm) )
    end do

    Faaa(:,:) = matmul( Ainv, Fxyz )

    do iatm=1,natm

      do isym=1,nsym

        RotRaaa(:) = matmul( SymMatA(1:3,1:3,isym), Raaa(:,iatm) ) + SymMatA(1:3,4,isym)
        call shift_in_01( RotRaaa )

        RotFaaa(:) = matmul( SymMatA(1:3,1:3,isym), Faaa(:,iatm) )

        do jatm=1,natm
          diff = sum( (RotRaaa-Raaa(:,jatm))**2 )
          if ( diff < 1.d-5 ) then
            Ftmp(:,jatm) = Ftmp(:,jatm) + RotFaaa(:)
            exit
          end if
        end do !jatm
        if ( jatm > natm ) stop 'error@sym_vector_force_xyz'

      end do !isym

    end do !iatm

    Fxyz = Ftmp/nsym

    deallocate( Faaa )
    deallocate( Raaa )
    deallocate( Ftmp )

  contains

    subroutine shift_in_01( a )
      implicit none
      real(8),intent(inout) :: a(3)
      integer :: i
      do i=1,3
        if ( abs(a(i)) < 1.0d-8 ) a(i)=0.0d0
        if ( abs(a(i)-1.0d0) < 1.0d-8 ) a(i)=0.0d0
        do
          if ( a(i) < 0.0d0 ) then
            a(i) = a(i) + 1.0d0
          else if ( a(i) >= 1.0d0 ) then
            a(i) = a(i) - 1.0d0
          else
            exit
          end if
        end do
      end do !i
    end subroutine shift_in_01

  end subroutine sym_vector_force_xyz

end module sym_vector_sub
