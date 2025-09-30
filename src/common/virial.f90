module virial_sub
  implicit none

contains

!===================================================================================================================================

  subroutine calc_virial(system,pp,fg,info,mg,stencil,poisson,srg,ppg,tpsi,ewald,energy,virial)
    use structures
    use salmon_math
    use math_constants,only: pi,zI
    use salmon_global, only: kion,cutoff_g,yn_periodic,yn_spinorbit, yn_jm!, aEwald, cutoff_r, yn_jm, yn_fix_func, theory
    use communication, only: comm_summation!,comm_is_root
    use sendrecv_grid, only: update_overlap_complex8
    use nonlocal_potential, only: calc_uVpsi_rdivided, calc_uVpsi
    use stencil_sub, only: calc_gradient_psi
    use timer
    implicit none
    type(s_dft_system)      ,intent(in)    :: system
    type(s_pp_info)         ,intent(in)    :: pp
    type(s_reciprocal_grid) ,intent(in)    :: fg
    type(s_parallel_info)   ,intent(in)    :: info
    type(s_rgrid)           ,intent(in)    :: mg
    type(s_stencil)         ,intent(in)    :: stencil
    type(s_poisson)         ,intent(in)    :: poisson
    type(s_sendrecv_grid)   ,intent(inout) :: srg
    type(s_pp_grid)         ,intent(in)    :: ppg
    type(s_orbital)         ,intent(inout) :: tpsi
    type(s_ewald_ion_ion)   ,intent(in)    :: ewald
    type(s_dft_energy)      ,intent(in)    :: energy
    type(s_dft_virial)      ,intent(inout) :: virial

    integer :: im, Nspin, ik_s, ik_e, io_s, io_e, Norb, Nlma
    integer :: ix, iy, iz, ia, ik, io, ispin
    integer :: ilocal, j, ilma
    real(8) :: P_tmp_loc, P_tmp_nloc, rx, ry, rz, P_tmp_loc_sum, P_tmp_nloc_sum
    real(8) :: kAc(3), rtmp, g(3), r(3), G2, Gd, ptmp, P_wrk(2), P_sum(2), sysvol
    complex(8) :: rho_i, rho_e, eGd, VG
    complex(8) :: w(3),duVpsi(3)
    complex(8),allocatable :: gtpsi(:,:,:,:),uVpsibox(:,:,:,:,:),uVpsibox2(:,:,:,:,:)

    if(info%im_s/=1 .or. info%im_e/=1) stop "error: calc_virial_periodic"
    im = 1
    Nspin = system%Nspin
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e
    Norb = system%Nspin*info%numo*info%numk
    Nlma = ppg%Nlma
    sysvol = system%det_a

    ! Ewald
    virial%P_ion_ion  = energy%E_ion_ion

    ! Fourier
    P_tmp_loc   = 0d0
    select case(yn_periodic)
    case('y')
      ptmp = 0d0
      P_wrk = 0d0
      !E_wrk_local_1 =0d0
      !E_wrk_local_2 =0d0
#ifdef USE_OPENACC
#else
!$omp parallel do collapse(2) default(none) &
!$omp          reduction(+:P_wrk,ptmp) &
!$omp          private(ix,iy,iz,g,rho_i,rho_e,ia,r,Gd) &
!$omp          shared(mg,fg,system,sysvol,kion,poisson,ppg,info,yn_jm)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        g(1) = fg%vec_G(1,ix,iy,iz)
        g(2) = fg%vec_G(2,ix,iy,iz)
        g(3) = fg%vec_G(3,ix,iy,iz)

        rho_e = poisson%zrhoG_ele(ix,iy,iz)

        if (yn_jm=='n') then
	  rho_i = ppg%zrhoG_ion(ix,iy,iz)
          P_wrk(1) = P_wrk(1) + sysvol* fg%coef(ix,iy,iz) * (-rho_e*conjg(rho_i))     ! electron-ion (valence)

          do ia=info%ia_s,info%ia_e
            r = system%Rion(1:3,ia)
            Gd = g(1)*r(1) + g(2)*r(2) + g(3)*r(3)
            ptmp = ptmp + conjg(rho_e)*ppg%zrVG_ion(ix,iy,iz,Kion(ia))*exp(-zI*Gd)  ! electron-ion (core)
          end do
        end if
      end do
      end do
      end do
!$omp end parallel do
#endif
      call comm_summation(ptmp,P_wrk(2),info%icomm_ko)
      call comm_summation(P_wrk,P_sum,2,info%icomm_r)
      ! electron-ion pressure (local part)
      virial%P_ion_loc = -1d0 * ( P_sum(1) + P_sum(2) )
    end select

    !Nonlocal part
    if( yn_spinorbit=='y' )then
      !call calc_uVpsi_so(nspin,info,ppg,tpsi,uVpsibox2)
    else
    ! uVpsibox2 = < uV | exp(ikr) | psi >
      if (info%if_divide_rspace) then
        call calc_uVpsi_rdivided(nspin,info,ppg,tpsi,uVpsibox,uVpsibox2)
      else
        call calc_uVpsi(nspin,info,ppg,tpsi,uVpsibox2)
      end if
    end if

    if(info%if_divide_rspace) then
       if(tpsi%update_zwf_overlap) call update_overlap_complex8(srg,mg,tpsi%zwf)
    end if

    kAc = 0d0
    if (.not. allocated(gtpsi)) then
       allocate(gtpsi(3,mg%is_array(1):mg%ie_array(1) &
            ,mg%is_array(2):mg%ie_array(2) &
            ,mg%is_array(3):mg%ie_array(3)))
    end if

    P_tmp_nloc = 0d0
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
       ! gtpsi = (nabla) psi
       call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,io,ik,im),gtpsi,mg%is_array,mg%ie_array,mg%is,mg%ie &
            ,mg%idx,mg%idy,mg%idz,stencil%coef_nab,system%rmatrix_B)
       !call timer_begin(LOG_CALC_FORCE_NONLOCAL)
       if(yn_periodic=='y') kAc(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
       rtmp = 2d0 * system%rocc(io,ik,ispin) * system%wtk(ik) * system%Hvol

       if( yn_spinorbit=='y' )then
          !to be written
       else
#ifdef USE_OPENACC
#else
#ifndef __NVCOMPILER_LLVM__
! FIXME: NVIDIA compiler crashes with nested omp parallel clause.
!$omp parallel do private(ilocal,ilma,ia,duVpsi,j,ix,iy,iz,rx,ry,rz,w) reduction(+:P_tmp_nloc)
#endif
          do ilocal=1,ppg%ilocal_nlma
             ilma=ppg%ilocal_nlma2ilma(ilocal)
             ia  =ppg%ilocal_nlma2ia  (ilocal)
             duVpsi = 0d0
!OCL swp
!OCL swp_freg_rate(115)
             do j=1,ppg%mps(ia)
                ix = ppg%jxyz(1,j,ia)
                iy = ppg%jxyz(2,j,ia)
                iz = ppg%jxyz(3,j,ia)
                rx = ppg%rxyz(1,j,ia)
                ry = ppg%rxyz(2,j,ia)
                rz = ppg%rxyz(3,j,ia)
                w(1) = gtpsi(1,ix,iy,iz) + zI* kAc(1) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
                w(2) = gtpsi(2,ix,iy,iz) + zI* kAc(2) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
                w(3) = gtpsi(3,ix,iy,iz) + zI* kAc(3) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
                duVpsi(1) = duVpsi(1) + conjg(ppg%zekr_uV(j,ilma,ik)) * rx * w(1) ! < uV | exp(ikr) (nabla) | psi >
                duVpsi(2) = duVpsi(2) + conjg(ppg%zekr_uV(j,ilma,ik)) * ry * w(2) ! < uV | exp(ikr) (nabla) | psi >
                duVpsi(3) = duVpsi(3) + conjg(ppg%zekr_uV(j,ilma,ik)) * rz * w(3) ! < uV | exp(ikr) (nabla) | psi >
             end do
             P_tmp_nloc = P_tmp_nloc &
                  + rtmp * dble( conjg(duVpsi(1)) * uVpsibox2(ispin,io,ik,im,ilma) ) &
                  + rtmp * dble( conjg(duVpsi(2)) * uVpsibox2(ispin,io,ik,im,ilma) ) &
                  + rtmp * dble( conjg(duVpsi(3)) * uVpsibox2(ispin,io,ik,im,ilma) )
          end do
#ifndef __NVCOMPILER_LLVM__
! FIXME: NVIDIA compiler crashes with nested omp parallel clause.
!$omp end parallel do
#endif
#endif          
       end if
    end do !ispin
    end do !io
    end do !ik

    if (allocated(gtpsi)) deallocate(gtpsi)

    call comm_summation(P_tmp_nloc,P_tmp_nloc_sum,1,info%icomm_rko)

    virial%P_ion_nloc = P_tmp_nloc_sum +  3.d0 * energy%E_ion_nloc

    if(allocated(uVpsibox)) deallocate(uVpsibox)

    ! kinetic virial (correlation part)
    virial%P_kin_c   = energy%T_c

    ! XC virial
    virial%P_xc       = virial%P_xc - energy%T_c

    ! kinetic virial
    virial%P_kin      = ( energy%E_kin + energy%T_c )*2d0

    ! Hartree virial
    virial%P_h = energy%E_h

    ! total virial
    virial%P_tot = virial%P_kin + virial%P_h + virial%P_ion_loc &
         + virial%P_ion_nloc + virial%P_xc + virial%P_ion_ion

    !1/3V
    virial%P_tot = virial%P_tot*system%dvinv
    virial%P_kin = virial%P_kin*system%dvinv
    virial%P_h = virial%P_h*system%dvinv
    virial%P_ion_loc = virial%P_ion_loc*system%dvinv
    virial%P_ion_nloc = virial%P_ion_nloc*system%dvinv
    virial%P_xc = virial%P_xc*system%dvinv
    virial%P_kin_c = virial%P_kin_c*system%dvinv
    virial%P_ion_ion = virial%P_ion_ion*system%dvinv

    return    
  end subroutine calc_virial
end module virial_sub
