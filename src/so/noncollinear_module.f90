module noncollinear_module

  use spin_orbit_global, only: SPIN_ORBIT_ON

  implicit none

  private
  public :: SPIN_ORBIT_ON
  public :: calc_dm_noncollinear
  public :: rot_dm_noncollinear
  public :: rot_vxc_noncollinear
  public :: op_xc_noncollinear

  complex(8),allocatable :: den_mat(:,:,:,:,:)
  complex(8),allocatable :: vxc_mat(:,:,:,:,:)
  complex(8),allocatable :: old_mat(:,:,:)
  complex(8),parameter :: zero=(0.0d0,0.0d0)
  real(8),allocatable :: rot_ang(:,:,:,:)

contains

  subroutine calc_dm_noncollinear( psi, system, info, mg )
    use structures, only : s_dft_system, s_parallel_info, s_rgrid, s_orbital
    implicit none
    type(s_orbital),intent(in) :: psi
    type(s_dft_system),intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_rgrid),intent(in) :: mg
    integer :: io,ik,im,is,js,ix,iy,iz,m1,m2,m3,n1,n2,n3
    complex(8),allocatable :: ztmp(:,:,:)
    real(8) :: occ

    if ( .not.allocated(den_mat) ) then
       m1=mg%is(1); n1=mg%ie(1)
       m2=mg%is(2); n2=mg%ie(2)
       m3=mg%is(3); n3=mg%ie(3)
       allocate( den_mat(m1:n1,m2:n2,m3:n3,2,2) )
       den_mat=zero
    end if

    den_mat=zero

    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
       occ=system%rocc(io,ik,1)*system%wtk(ik)
       if ( abs(occ) < 1.0d-15 ) cycle
       do js=1,2
       do is=1,2
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
             den_mat(ix,iy,iz,is,js) = den_mat(ix,iy,iz,is,js) + occ &
                  * conjg( psi%zwf(ix,iy,iz,is,io,ik,im) )*psi%zwf(ix,iy,iz,js,io,ik,im)
          end do
          end do
          end do
       end do !is
       end do !js
    end do !io
    end do !ik
    end do !im

    !m1=mg%is(1); n1=mg%ie(1)
    !m2=mg%is(2); n2=mg%ie(2)
    !m3=mg%is(3); n3=mg%ie(3)
    !allocate( ztmp(m1:n1,m2:n2,m3:n3) ); ztmp=zero
    !ztmp(:,:,:)=ztmp(:,:,:)+den_mat(:,:,:,1,1)+den_mat(:,:,:,2,2)
    !write(*,*) sum(real(ztmp))*system%hvol,sum(aimag(ztmp))*system%hvol
    !write(*,*) minval(real(ztmp)),minval(aimag(ztmp))
    !write(*,*) maxval(real(ztmp)),maxval(aimag(ztmp))
    !deallocate(ztmp)
 
  end subroutine calc_dm_noncollinear


  subroutine rot_dm_noncollinear( rho, system, mg )
    use structures, only : s_dft_system, s_rgrid, s_scalar
    use communication, only : comm_summation
    use mpi, only: MPI_COMM_WORLD
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_scalar),intent(inout) :: rho(system%nspin)
    real(8) :: phi,theta,tmp,tmp1
    integer :: a,b,m1,m2,m3,n1,n2,n3,ix,iy,iz

    if ( .not.allocated(rot_ang) ) then
       m1=mg%is(1) ; n1=mg%ie(1)
       m2=mg%is(2) ; n2=mg%ie(2)
       m3=mg%is(3) ; n3=mg%ie(3)
       allocate( rot_ang(m1:n1,m2:n2,m3:n3,2) ) ; rot_ang=0.0d0
    end if

    rot_ang=0.0d0

    a=1
    b=2
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)

       phi = -atan( aimag(den_mat(ix,iy,iz,a,b))/real(den_mat(ix,iy,iz,a,b)) )
       theta = atan( 2.0d0*( real(den_mat(ix,iy,iz,a,b))*cos(phi) &
            -aimag(den_mat(ix,iy,iz,a,b))*sin(phi) ) &
            /real( den_mat(ix,iy,iz,a,a)-den_mat(ix,iy,iz,b,b) ) )

       rho(1)%f(ix,iy,iz) = 0.5d0*real( den_mat(ix,iy,iz,a,a)+den_mat(ix,iy,iz,b,b) ) &
            + 0.5d0*real( den_mat(ix,iy,iz,a,a)-den_mat(ix,iy,iz,b,b) )*cos(theta) &
            + (  real(den_mat(ix,iy,iz,a,b))*cos(phi) &
            -aimag(den_mat(ix,iy,iz,a,b))*sin(phi) )*sin(theta)

       rho(2)%f(ix,iy,iz) = 0.5d0*real( den_mat(ix,iy,iz,a,a)+den_mat(ix,iy,iz,b,b) ) &
            - 0.5d0*real( den_mat(ix,iy,iz,a,a)-den_mat(ix,iy,iz,b,b) )*cos(theta) &
            - (  real(den_mat(ix,iy,iz,a,b))*cos(phi) &
            -aimag(den_mat(ix,iy,iz,a,b))*sin(phi) )*sin(theta)

       rot_ang(ix,iy,iz,1) = phi
       rot_ang(ix,iy,iz,2) = theta

    end do !ix
    end do !iy
    end do !iz

    !write(*,*) "size(rho%f)",(size(rho(1)%f,ix),ix=1,3)
    !tmp=sum(rho(1)%f+rho(2)%f)*system%hvol
    !call comm_summation( tmp, tmp1, MPI_COMM_WORLD )
    !write(*,*) "sum(rho)@rot_dm_noncollinear",tmp1
    !write(*,*) minval(rho(1)%f),minval(rho(2)%f)
    !write(*,*) maxval(rho(1)%f),maxval(rho(2)%f)

  end subroutine rot_dm_noncollinear


  subroutine rot_vxc_noncollinear( Vxc, system, mg )
    use structures, only : s_dft_system, s_rgrid, s_scalar
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_scalar),intent(inout) :: Vxc(system%nspin)
    real(8) :: phi,theta,vxc_0,vxc_1
    integer :: ix,iy,iz,m1,m2,m3,n1,n2,n3

    if ( .not.allocated(vxc_mat) ) then
       m1=mg%is(1); n1=mg%ie(1)
       m2=mg%is(2); n2=mg%ie(2)
       m3=mg%is(3); n3=mg%ie(3)
       allocate( vxc_mat(m1:n1,m2:n2,m3:n3,2,2) )
    end if
    vxc_mat=zero

    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)

       phi = rot_ang(ix,iy,iz,1)
       theta = rot_ang(ix,iy,iz,2)

       vxc_0 = 0.5d0*( Vxc(1)%f(ix,iy,iz) + Vxc(2)%f(ix,iy,iz) )
       vxc_1 = 0.5d0*( Vxc(1)%f(ix,iy,iz) - Vxc(2)%f(ix,iy,iz) )

       vxc_mat(ix,iy,iz,1,1) = vxc_0 + vxc_1*cos(theta)
       vxc_mat(ix,iy,iz,2,1) = vxc_1*cmplx( cos(phi), sin(phi) )*sin(theta)
       vxc_mat(ix,iy,iz,1,2) = vxc_1*cmplx( cos(phi),-sin(phi) )*sin(theta)
       vxc_mat(ix,iy,iz,2,2) = vxc_0 - vxc_1*cos(theta)

    end do !ix
    end do !iy
    end do !iz

    Vxc(1)%f=0.0d0
    Vxc(2)%f=0.0d0

  end subroutine rot_vxc_noncollinear


  subroutine op_xc_noncollinear( tpsi, hpsi, info, mg )
    use structures, only : s_orbital, s_rgrid, s_parallel_info
    implicit none
    type(s_orbital),intent(in) :: tpsi
    type(s_orbital),intent(inout) :: hpsi
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    integer :: ix,iy,iz,im,ik,io
    if ( .not.allocated(vxc_mat) ) return
    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
       do iz=mg%is(3),mg%ie(3)
       do iy=mg%is(2),mg%ie(2)
       do ix=mg%is(1),mg%ie(1)
          hpsi%zwf(ix,iy,iz,1,io,ik,im) = hpsi%zwf(ix,iy,iz,1,io,ik,im) &
               + vxc_mat(ix,iy,iz,1,1)*tpsi%zwf(ix,iy,iz,1,io,ik,im) &
               + vxc_mat(ix,iy,iz,1,2)*tpsi%zwf(ix,iy,iz,2,io,ik,im)
          hpsi%zwf(ix,iy,iz,2,io,ik,im) = hpsi%zwf(ix,iy,iz,2,io,ik,im) &
               + vxc_mat(ix,iy,iz,2,1)*tpsi%zwf(ix,iy,iz,1,io,ik,im) &
               + vxc_mat(ix,iy,iz,2,2)*tpsi%zwf(ix,iy,iz,2,io,ik,im)
       end do
       end do
       end do
    end do
    end do
    end do
  end subroutine op_xc_noncollinear


end module noncollinear_module
