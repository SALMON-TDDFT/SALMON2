! Bounded destination-grid pieces for validated complex DC reconstruction.
module lcfo_mesh_tile
  implicit none
  private
  public :: lcfo_tile_coverage,lcfo_tile_contract,lcfo_reconstruction_tile_points
  integer,parameter :: lcfo_reconstruction_tile_points=65536
contains
  pure subroutine intersect_axes(core,jxyz,lo,m,indices,counts)
    integer,intent(in) :: core(3),jxyz(:,:),lo(3),m(3)
    integer,intent(out) :: indices(:,:),counts(3)
    integer :: axis,i
    counts=0
    do axis=1,3
      do i=1,core(axis)
        if(jxyz(i,axis)<lo(axis).or.jxyz(i,axis)>=lo(axis)+m(axis))cycle
        counts(axis)=counts(axis)+1
        indices(counts(axis),axis)=i
      enddo
    enddo
  end subroutine

  pure subroutine lcfo_tile_coverage(core,jxyz,lo,m,first,coverage)
    integer,intent(in) :: core(3),jxyz(:,:),lo(3),m(3),first
    integer,intent(inout) :: coverage(:)
    integer :: indices(maxval(core),3),counts(3),i,j,k,ix,iy,iz,g
    call intersect_axes(core,jxyz,lo,m,indices,counts)
    do k=1,counts(3);do j=1,counts(2);do i=1,counts(1)
      ix=indices(i,1);iy=indices(j,2);iz=indices(k,3)
      g=1+jxyz(ix,1)-lo(1)+m(1)*(jxyz(iy,2)-lo(2)+m(2)*(jxyz(iz,3)-lo(3)))
      g=g-first+1
      if(g>=1.and.g<=size(coverage))coverage(g)=coverage(g)+1
    enddo;enddo;enddo
  end subroutine

  pure subroutine lcfo_tile_contract(jxyz,lo,m,first,basis,coef,tile)
    integer,intent(in) :: jxyz(:,:),lo(3),m(3),first
    complex(8),intent(in) :: basis(:,:,:,:),coef(:)
    complex(8),intent(inout) :: tile(:)
    integer :: core(3),indices(max(size(basis,1),size(basis,2),size(basis,3)),3),counts(3),i,j,k,ix,iy,iz,g,b
    core=[size(basis,1),size(basis,2),size(basis,3)]
    if(size(coef)==0)return
    call intersect_axes(core,jxyz,lo,m,indices,counts)
    do k=1,counts(3);do j=1,counts(2);do i=1,counts(1)
      ix=indices(i,1);iy=indices(j,2);iz=indices(k,3)
      g=1+jxyz(ix,1)-lo(1)+m(1)*(jxyz(iy,2)-lo(2)+m(2)*(jxyz(iz,3)-lo(3)))
      g=g-first+1
      if(g<1.or.g>size(tile))cycle
      do b=1,size(coef)
        tile(g)=tile(g)+basis(ix,iy,iz,b)*coef(b)
      enddo
    enddo;enddo;enddo
  end subroutine
end module
