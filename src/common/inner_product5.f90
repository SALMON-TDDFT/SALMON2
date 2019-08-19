subroutine inner_product5(mg,info,itotmst,mst,iobnum,zmatbox1,zmatbox2,zbox,hvol)
  use structures, only: s_rgrid, s_orbital_parallel
  use salmon_communication, only: comm_summation
  use timer
  use calc_allob_sub
  !$ use omp_lib
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_orbital_parallel),intent(in) :: info
  integer,intent(in)  :: itotmst
  integer,intent(in)  :: mst(2)
  integer,intent(in)  :: iobnum
  complex(8) :: zmatbox1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:iobnum)
  complex(8) :: zmatbox2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:iobnum)
  complex(8),intent(out) :: zbox(1:itotmst)
  real(8),intent(in)  :: hvol
  integer :: iob,iob_allob
  integer :: ix,iy,iz
  complex(8) :: zbox2(1:itotmst)
  complex(8) :: sum0
  
  zbox2(:)=0.d0
  do iob=1,iobnum
    call calc_allob(iob,iob_allob,itotmst,mst,iobnum)
    sum0=0.d0
    !$OMP parallel do collapse(2) reduction(+ : sum0)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      sum0=sum0+conjg(zmatbox1(ix,iy,iz,iob))*zmatbox2(ix,iy,iz,iob)
    end do
    end do
    end do
    zbox2(iob_allob)=sum0*hvol
  end do

  call timer_begin(LOG_ALLREDUCE_INNER_PRODUCT5)
  call comm_summation(zbox2,zbox,itotmst,info%icomm_r)
  call timer_end(LOG_ALLREDUCE_INNER_PRODUCT5)

end subroutine inner_product5
