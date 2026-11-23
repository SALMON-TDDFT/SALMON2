module sym_spin_density_sub
  use sym_sub, only: use_symmetry,SymMatA,SymMatB,Amat,Ainv
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: sym_spin_density_grid
contains
  subroutine sym_spin_density_grid(dm)
    complex(8),intent(inout)::dm(:,:,:,:,:)
    complex(8),allocatable::original(:,:,:,:,:)
    complex(8)::local(2,2),rotated(2,2)
    real(8)::rotation(3,3),axial(3,3),unit(3,3),determinant,point(3),frac(3),mag(3),n
    integer::g(3),index(3),op,i,j,k,nsym,q,isp,jsp
    if(.not.use_symmetry)return
    if(.not.allocated(SymMatA).or..not.allocated(SymMatB))error stop 'SOC density: symmetry unavailable'
    nsym=size(SymMatA,3)
    if(nsym<1.or.size(SymMatB,3)/=nsym)error stop 'SOC density: invalid symmetry count'
    if(size(dm,4)/=2.or.size(dm,5)/=2)error stop 'SOC density: spin matrix shape'
    do jsp=1,2
    do isp=1,2
    do k=1,size(dm,3)
    do j=1,size(dm,2)
    do i=1,size(dm,1)
      if(.not.ieee_is_finite(real(dm(i,j,k,isp,jsp))).or.&
         .not.ieee_is_finite(aimag(dm(i,j,k,isp,jsp)))) &
        error stop 'SOC density: nonfinite matrix'
    enddo
    enddo
    enddo
    enddo
    enddo
    if(maxval(abs(dm(:,:,:,1,2)-conjg(dm(:,:,:,2,1))))>1d-10*max(1d0,maxval(abs(dm)))) &
      error stop 'SOC density: non-Hermitian matrix'
    g=shape(dm(:,:,:,1,1));unit=0d0
    do q=1,3
      unit(q,q)=1d0
    enddo
    allocate(original,source=dm);dm=0d0
    do op=1,nsym
      rotation=matmul(Amat,matmul(SymMatA(:,1:3,op),Ainv))
      if(.not.matrix_is_finite(rotation).or..not.matrix_is_finite(SymMatB(:,:,op)).or. &
         .not.matrix_is_finite(SymMatA(:,:,op))) &
        error stop 'SOC density: nonfinite operation'
      if(maxval(abs(matmul(transpose(rotation),rotation)-unit))>1d-10) &
        error stop 'SOC density: nonorthogonal Cartesian rotation'
      if(maxval(abs(matmul(transpose(SymMatB(:,1:3,op)),SymMatA(:,1:3,op))-unit))>1d-10) &
        error stop 'SOC density: inconsistent inverse rotation'
      determinant=rotation(1,1)*(rotation(2,2)*rotation(3,3)-rotation(2,3)*rotation(3,2)) &
                 -rotation(1,2)*(rotation(2,1)*rotation(3,3)-rotation(2,3)*rotation(3,1)) &
                 +rotation(1,3)*(rotation(2,1)*rotation(3,2)-rotation(2,2)*rotation(3,1))
      if(abs(abs(determinant)-1d0)>1d-10)error stop 'SOC density: invalid rotation determinant'
      axial=sign(1d0,determinant)*rotation
      do k=1,g(3)
      do j=1,g(2)
      do i=1,g(1)
        frac=real([i-1,j-1,k-1],8)/g
        ! Reciprocal rotation transposed is the inverse DIRECT rotation.
        point=matmul(transpose(SymMatB(:,1:3,op)),frac-SymMatA(:,4,op))*g
        if(maxval(abs(point-anint(point)))>1d-8)error stop 'SOC density: noncommensurate grid'
        index=modulo(nint(point),g)+1
        local=original(index(1),index(2),index(3),:,:)
        n=real(local(1,1)+local(2,2),8)
        mag=[2d0*real(local(1,2),8),-2d0*aimag(local(1,2)),real(local(1,1)-local(2,2),8)]
        mag=matmul(axial,mag)
        rotated(1,1)=0.5d0*(n+mag(3));rotated(2,2)=0.5d0*(n-mag(3))
        rotated(1,2)=0.5d0*cmplx(mag(1),-mag(2),8);rotated(2,1)=conjg(rotated(1,2))
        dm(i,j,k,:,:)=dm(i,j,k,:,:)+rotated/nsym
      enddo
      enddo
      enddo
    enddo
  end subroutine
  logical function matrix_is_finite(a) result(ok)
    implicit none
    real(8),intent(in)::a(:,:)
    integer::i,j
    ok=.true.
    do j=1,size(a,2)
      do i=1,size(a,1)
        if(.not.ieee_is_finite(a(i,j)))ok=.false.
      enddo
    enddo
  end function
end module
