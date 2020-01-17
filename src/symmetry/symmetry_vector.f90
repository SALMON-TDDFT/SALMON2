module sym_vector_sub

  use sym_sub, only: Amat, Ainv, SymMatB, use_symmetry

  implicit none

  private
  public :: sym_vector_xyz

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

end module sym_vector_sub
