module prep_density_pp_sub

  implicit none

  private
  public :: calc_density_pp

contains

  subroutine calc_density_pp(lg,mg,system,info,pp,fg,poisson,rho)
    use salmon_global,only : nelem, kion, yn_ffte
    use communication, only: comm_summation
    use structures
    use sym_rho_sub, only: sym_rho
    implicit none
    type(s_rgrid)          ,intent(in) :: lg,mg
    type(s_dft_system)     ,intent(in) :: system
    type(s_parallel_info)  ,intent(in) :: info
    type(s_pp_info)        ,intent(in) :: pp
    type(s_reciprocal_grid),intent(in) :: fg
    type(s_poisson)                    :: poisson
    type(s_scalar)                     :: rho(:)
    !
    integer :: ia,i,ik,ix,iy,iz,kx,ky,kz,iiy,iiz,nspin
    real(8) :: g(3),gd,s,g2sq,const,sb,x
    complex(8) :: tmp_exp
    complex(8) :: vtmp1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    complex(8) :: vtmp2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    complex(8),allocatable :: zrhoG(:,:,:,:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)

    allocate( zrhoG(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),nelem) )
    zrhoG=zero

    const = 1.0d0/system%det_A
    !!$omp parallel
    !!$omp do private(ik,ix,iy,iz,g,g2sq,s,r1,dr,i,vloc_av) collapse(3)
    do ik=1,nelem
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        g(1) = fg%vec_G(1,ix,iy,iz)
        g(2) = fg%vec_G(2,ix,iy,iz)
        g(3) = fg%vec_G(3,ix,iy,iz)
        g2sq = sqrt(g(1)**2+g(2)**2+g(3)**2)
        if ( g2sq == 0.0d0 ) then
          zrhoG(ix,iy,iz,ik) = pp%zps(ik)*const
          cycle
        end if 
        s=0.0d0
        do i = 1, pp%mr(ik)
          x = g2sq*pp%rad(i,ik)
          ! --- Spherical Bessel (l=0) ---
          if ( x < 1.d-1 ) then
            sb=-(1.d0/39916800.d0*x**10-1.d0/362880.d0*x**8 &
                +1.d0/5040.d0*x**6-1.d0/120.d0*x**4+1.d0/6.d0*x**2-1.d0)
          else
            sb=sin(x)/x
          end if
          ! ---
          s = s + pp%rho_pp_tbl(i,ik)*sb*(pp%rad(i+1,ik)-pp%rad(i,ik))
        end do
        zrhoG(ix,iy,iz,ik) = s*const
      end do !ix
      end do !iy
      end do !iz
    end do !ik
    !!$omp end do
    !!$omp end parallel

    vtmp1=zero
    !!$omp parallel do collapse(2) private(ix,iy,iz,g,ia,ik,gd,tmp_exp)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      g(1) = fg%vec_G(1,ix,iy,iz)
      g(2) = fg%vec_G(2,ix,iy,iz)
      g(3) = fg%vec_G(3,ix,iy,iz)
      do ia=info%ia_s,info%ia_e
        ik=kion(ia)
        gd = g(1)*system%Rion(1,ia) + g(2)*system%Rion(2,ia) + g(3)*system%Rion(3,ia)
        tmp_exp = dcmplx( cos(gd), -sin(gd) )
        vtmp1(ix,iy,iz) = vtmp1(ix,iy,iz) + zrhoG(ix,iy,iz,ik)*tmp_exp
      end do
    end do
    end do
    end do
    !!$omp end parallel do
  
    call comm_summation( vtmp1, vtmp2, size(vtmp2), info%icomm_ko )

    if ( yn_ffte == 'n' ) then

      !$omp parallel workshare
      poisson%ff1x = zero
      poisson%ff1y = zero
      poisson%ff1z = zero
      !$omp end parallel workshare

      !$omp parallel do private(kz,ky,kx)
      do kz = mg%is(3),mg%ie(3)
      do ky = mg%is(2),mg%ie(2)
      do kx = mg%is(1),mg%ie(1)
        poisson%ff1z(kx,ky,kz) = vtmp2(kx,ky,kz)
      end do
      end do
      end do
      call comm_summation(poisson%ff1z,poisson%ff2z,mg%num(1)*mg%num(2)*lg%num(3),info%icomm_z)

      !$omp parallel do private(iz,ky,kx)
      do iz = mg%is(3),mg%ie(3)
      do ky = mg%is(2),mg%ie(2)
      do kx = mg%is(1),mg%ie(1)
        poisson%ff1y(kx,ky,iz) = sum(fg%egz(:,iz)*poisson%ff2z(kx,ky,:))
      end do
      end do
      end do
      call comm_summation(poisson%ff1y,poisson%ff2y,mg%num(1)*lg%num(2)*mg%num(3),info%icomm_y)

      !$omp parallel do private(iz,iy,kx)
      do iz = mg%is(3),mg%ie(3)
      do iy = mg%is(2),mg%ie(2)
      do kx = mg%is(1),mg%ie(1)
        poisson%ff1x(kx,iy,iz)=sum(fg%egy(:,iy)*poisson%ff2y(kx,:,iz))
      end do
      end do
      end do
      call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*mg%num(2)*mg%num(3),info%icomm_x)

      !$omp parallel do private(iz,iy,ix) collapse(2)
      do iz = mg%is(3),mg%ie(3)
      do iy = mg%is(2),mg%ie(2)
      do ix = mg%is(1),mg%ie(1)
        rho(1)%f(ix,iy,iz) = sum(fg%egx(:,ix)*poisson%ff2x(:,iy,iz))
      end do
      end do
      end do
  
    else ! yn_ffte==.true.
 
      poisson%b_ffte=zero
      !$omp parallel do private(iz,iy,ix,iiz,iiy) collapse(2)
      do iz = 1,mg%num(3)
      do iy = 1,mg%num(2)
      do ix = mg%is(1),mg%ie(1)
        iiz=iz+mg%is(3)-1
        iiy=iy+mg%is(2)-1
        poisson%b_ffte(ix,iy,iz) = vtmp2(ix,iiy,iiz)
      end do
      end do
      end do
      call comm_summation(poisson%b_ffte,poisson%a_ffte,size(poisson%a_ffte),info%icomm_x)

      call PZFFT3DV_MOD(poisson%a_ffte,poisson%b_ffte,lg%num(1),lg%num(2),lg%num(3),   &
                        info%isize_y,info%isize_z,1, &
                        info%icomm_y,info%icomm_z)

      !$omp parallel do private(iz,iy,iiz,iiy) collapse(2)
      do iz=1,mg%num(3)
      do iy=1,mg%num(2)
        iiz=iz+mg%is(3)-1
        iiy=iy+mg%is(2)-1
        rho(1)%f(mg%is(1):mg%ie(1),iiy,iiz) = poisson%b_ffte(mg%is(1):mg%ie(1),iy,iz)*system%ngrid
      end do
      end do

    end if

    call sym_rho( rho(1)%f(:,:,:) )

    nspin = size(rho)
    if ( nspin == 2 ) then
      rho(1)%f(:,:,:) = 0.5d0*rho(1)%f(:,:,:)
      rho(2)%f(:,:,:) = rho(1)%f(:,:,:)
    end if

    !write(*,*) "sum(rho)=",sum(rho(1)%f)*system%hvol
    !rewind 200
    !do ix = mg%is(1),mg%ie(1)
    !  write(200,*) (ix-1)*system%hgs(1), rho(1)%f(ix,1,1)
    !end do
    !stop

    return

  end subroutine calc_density_pp

end module prep_density_pp_sub
