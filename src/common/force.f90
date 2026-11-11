!
!  Copyright 2019-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module force_sub
  implicit none

contains

!===================================================================================================================================

  subroutine calc_force(system,pp,fg,info,mg,stencil,poisson,srg,ppg,tpsi,ewald)
    use structures
    use math_constants,only : zi,pi
    use sendrecv_grid, only: s_sendrecv_grid, update_overlap_real8, update_overlap_complex8, dealloc_cache
    use communication, only: comm_summation
    use nonlocal_potential, only: calc_uVpsi_rdivided, calc_uVpsi
    use sym_vector_sub, only: sym_vector_force_xyz
    use sym_sub, only: use_symmetry
    use pseudo_pt_so_sub, only: calc_uVpsi_so
    use plusU_global, only: PLUS_U_ON, dm_mms_nla, U_eff
    use salmon_global, only: kion,cutoff_g,yn_periodic,yn_spinorbit
    use code_optimization, only: force_omp_mode
    use stencil_sub, only: calc_gradient_psi
    use timer
    implicit none
    type(s_dft_system)      ,intent(inout) :: system
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
    !
    integer :: ix,iy,iz,ia,nion,im,Nspin,ik_s,ik_e,io_s,io_e,nlma,ik,io,ispin,ilma,j
    integer :: m1,m2,jlma,n,l,Nproj_pairs,iprj,Nlma_ao
    real(8) :: kAc(3), rtmp,rtmp2(3),g(3),r(3),G2,Gd
    complex(8) :: rho_i, rho_e, egd, VG
    real(8),allocatable :: F_tmp(:,:),F_sum(:,:)
    real(8),allocatable :: dden(:,:,:,:)
    complex(8) :: w(3),duVpsi(3)
    complex(8),allocatable :: gtpsi(:,:,:,:),uVpsibox(:,:,:,:,:),uVpsibox2(:,:,:,:,:)
    complex(8),allocatable :: phipsibox(:,:,:,:),phipsibox2(:,:,:,:)
    complex(8),allocatable :: dphipsi_lma(:,:)
    complex(8) :: ddm_mms_nla(3), phipsi, dphipsi(3)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: zF_tmp(:,:)
    complex(8) :: ztmp
    integer :: Norb,iorb,ilocal

    call timer_begin(LOG_CALC_ION_FORCE)

    nion = system%nion
    if(.not.allocated(system%Force)) allocate(system%Force(3,nion))
    allocate( F_tmp(3,nion), F_sum(3,nion) )

    if(info%im_s/=1 .or. info%im_e/=1) stop "error: calc_force_periodic" !??????
    im = 1

    Nspin = system%Nspin
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e
    Norb = system%Nspin*info%numo*info%numk

    Nlma = ppg%Nlma

  ! Ewald sum
    call timer_begin(LOG_CALC_FORCE_ION_ION)
    call force_ewald_rspace(system%Force,F_tmp,system,info,ewald,pp,nion,info%icomm_r)
    call timer_end(LOG_CALC_FORCE_ION_ION)


    call timer_begin(LOG_CALC_FORCE_FOURIER)
    F_tmp   = 0d0

    select case(yn_periodic)
    case('n')
      allocate( dden(3,mg%is_array(1):mg%ie_array(1) &
                      ,mg%is_array(2):mg%ie_array(2) &
                      ,mg%is_array(3):mg%ie_array(3)))
      dden = 0d0
    case('y')

    allocate(dden(1,1,1,1)); dden=0d0 !! dummy allocation to avoid crush: not use in isolated system

    ! Fourier part (local part, etc)
      !$omp parallel do collapse(2) private(ix,iy,iz,ia,r,g,G2,Gd,rho_i,rho_e,rtmp,egd,VG) reduction(+:F_tmp)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        if(fg%if_Gzero(ix,iy,iz)) cycle
        g(1) = fg%vec_G(1,ix,iy,iz)
        g(2) = fg%vec_G(2,ix,iy,iz)
        g(3) = fg%vec_G(3,ix,iy,iz)
        G2 = sum(g(:)**2)
        if(G2 .gt. cutoff_g**2) cycle   !xxx

        rho_i = ppg%zrhoG_ion(ix,iy,iz)
        rho_e = poisson%zrhoG_ele(ix,iy,iz)

        !OCL swp
        do ia=info%ia_s,info%ia_e
          r = system%Rion(1:3,ia)
          Gd = sum(g(:)*r(:))
          egd = exp(zI*Gd)
          rtmp = pp%Zps(Kion(ia))* fg%coef(ix,iy,iz) * fg%exp_ewald(ix,iy,iz)
          VG = ppg%zVG_ion(ix,iy,iz,kion(ia)) - fg%coef(ix,iy,iz) * pp%zps(kion(ia))
          F_tmp(:,ia) = F_tmp(:,ia) + g(:)* ( rtmp * aimag(rho_i*egd) + aimag(egd*rho_e*conjg(VG)) )
        end do
      end do
      end do
      end do
      !$omp end parallel do
    end select
    call timer_end(LOG_CALC_FORCE_FOURIER)


    call timer_begin(LOG_CALC_FORCE_ELEC_ION)
    if(allocated(tpsi%rwf)) then
      allocate(tpsi%zwf(mg%is_array(1):mg%ie_array(1) &
                       ,mg%is_array(2):mg%ie_array(2) &
                       ,mg%is_array(3):mg%ie_array(3) &
                       ,nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e))
      tpsi%zwf = dcmplx(tpsi%rwf)
    end if

    if( yn_spinorbit=='y' )then
      call calc_uVpsi_so(nspin,info,ppg,tpsi,uVpsibox2)
    else
    ! uVpsibox2 = < uV | exp(ikr) | psi >
      if (info%if_divide_rspace) then
        call calc_uVpsi_rdivided(nspin,info,ppg,tpsi,uVpsibox,uVpsibox2)
      else
        call calc_uVpsi(nspin,info,ppg,tpsi,uVpsibox2)
      end if
    end if

    if( PLUS_U_ON )then
!!!!! future work: force for PLUS_U_ON
      
      Nlma_ao = size(ppg%ia_tbl_ao)
!      allocate( zF_tmp(3,nion) ); zF_tmp=zero
!      allocate( phipsibox(Nlma_ao,nspin,info%io_s:info%io_e,info%ik_s:info%ik_e)  ); phipsibox=zero
!      allocate( phipsibox2(Nlma_ao,nspin,info%io_s:info%io_e,info%ik_s:info%ik_e) ); phipsibox2=zero
!      allocate( dphipsi_lma(3,Nlma_ao) ); dphipsi_lma=zero
    
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin
        do ilma=1,Nlma_ao
          ia = ppg%ia_tbl_ao(ilma)
          phipsi = 0.0d0
          do j=1,ppg%mps_ao(ia)
            ix = ppg%jxyz_ao(1,j,ia)
            iy = ppg%jxyz_ao(2,j,ia)
            iz = ppg%jxyz_ao(3,j,ia)
            phipsi = phipsi + conjg(ppg%zekr_phi_ao(j,ilma,ik)) &
                            * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
          end do
          phipsi = phipsi * system%Hvol
!          phipsibox(ilma,ispin,io,ik) = phipsi
        end do
      end do
      end do
      end do
!      call comm_summation(phipsibox,phipsibox2,Nlma_ao*Norb,info%icomm_r)
    end if
    call timer_end(LOG_CALC_FORCE_ELEC_ION)

    if(info%if_divide_rspace) then
       if(tpsi%update_zwf_overlap) call update_overlap_complex8(srg,mg,tpsi%zwf)
    end if

    call timer_begin(LOG_CALC_FORCE_ELEC_ION)
    kAc = 0d0

#ifdef USE_OPENACC
    ! Without spin-orbit coupling, all local orbitals are done in one device pass;
    ! the spin-orbit force keeps the orbital loop below.
    if( yn_spinorbit=='n' )then
      call calc_force_nonlocal_gpu(system,info,mg,stencil,ppg,tpsi%zwf,uVpsibox2,F_tmp)
      if(yn_periodic=='n') call calc_dden_gpu(system,info,mg,stencil,tpsi%zwf,dden)
    else
#endif

#ifndef USE_OPENACC
!$omp parallel default(none) &
!$omp   firstprivate(kAc) &
!$omp   private(ik,io,ispin,rtmp,rtmp2,ilocal,ilma,ia,duVpsi,j,ix,iy,iz,w,dphipsi_lma,dphipsi,Nproj_pairs) &
!$omp   private(jlma,l,n,m1,m2,ddm_mms_nla,gtpsi) &
!$omp   shared(im,ik_s,ik_e,io_s,io_e,nspin,tpsi,mg,stencil,system,ppg,uVpsibox2,yn_periodic) &
!$omp   shared(PLUS_U_ON,Nlma_ao,phipsibox2,iorb,zF_tmp,U_eff,dm_mms_nla) &
!$omp   shared(yn_spinorbit,ztmp) &
!$omp   shared(dden) &
!$omp   reduction(+:F_tmp) &
!$omp   if(force_omp_mode)
#endif

    if (.not. allocated(gtpsi)) then
      allocate(gtpsi(3,mg%is_array(1):mg%ie_array(1) &
                      ,mg%is_array(2):mg%ie_array(2) &
                      ,mg%is_array(3):mg%ie_array(3)))
    end if

#ifndef USE_OPENACC
!$omp do collapse(2) reduction(+:dden)
#endif
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

       ! gtpsi = (nabla) psi
       !call timer_begin(LOG_CALC_FORCE_GTPSI)
       call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,io,ik,im),gtpsi,mg%is_array,mg%ie_array,mg%is,mg%ie &
            ,mg%idx,mg%idy,mg%idz,stencil%coef_nab,system%rmatrix_B)
       !call timer_end(LOG_CALC_FORCE_GTPSI)
       
       if(yn_periodic=='n') then
         rtmp = 2d0 * system%rocc(io,ik,ispin) * system%wtk(ik) * system%Hvol
!$omp parallel do collapse(2) private(iz,iy,ix,w,rtmp2)
         do iz=mg%is(3),mg%ie(3)
         do iy=mg%is(2),mg%ie(2)
         do ix=mg%is(1),mg%ie(1)
           w = conjg(gtpsi(:,ix,iy,iz)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
           rtmp2(:) = rtmp * dble(w(:))
           dden(:,ix,iy,iz) = dden(:,ix,iy,iz) + rtmp2(:)
         enddo
         enddo
         enddo
!$omp end parallel do
       end if

       ! nonlocal part
       !call timer_begin(LOG_CALC_FORCE_NONLOCAL)
       if(yn_periodic=='y') kAc(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
       rtmp = 2d0 * system%rocc(io,ik,ispin) * system%wtk(ik) * system%Hvol

       if( yn_spinorbit=='y' )then
         do ilma=1,size(ppg%ia_tbl_so)
           ia = ppg%ia_tbl_so(ilma)
           duVpsi=zero
           do j=1,ppg%mps(ia)
             ix = ppg%jxyz(1,j,ia)
             iy = ppg%jxyz(2,j,ia)
             iz = ppg%jxyz(3,j,ia)
             w(1) = gtpsi(1,ix,iy,iz) + zI*kAc(1) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
             w(2) = gtpsi(2,ix,iy,iz) + zI*kAc(2) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
             w(3) = gtpsi(3,ix,iy,iz) + zI*kAc(3) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
             duVpsi(1) = duVpsi(1) + conjg( ppg%zekr_uV_so(j,ilma,ik,ispin,1) )*w(1)
             duVpsi(2) = duVpsi(2) + conjg( ppg%zekr_uV_so(j,ilma,ik,ispin,1) )*w(2)
             duVpsi(3) = duVpsi(3) + conjg( ppg%zekr_uV_so(j,ilma,ik,ispin,1) )*w(3)
           end do
           ztmp=uVpsibox2(1,ilma,io,ik,im)+uVpsibox2(2,ilma,io,ik,im)
           F_tmp(1,ia) = F_tmp(1,ia) - rtmp*dble( conjg(duVpsi(1)) * ztmp )
           F_tmp(2,ia) = F_tmp(2,ia) - rtmp*dble( conjg(duVpsi(2)) * ztmp )
           F_tmp(3,ia) = F_tmp(3,ia) - rtmp*dble( conjg(duVpsi(3)) * ztmp )
         end do
       else
#ifndef __NVCOMPILER_LLVM__ 
! FIXME: NVIDIA compiler crashes with nested omp parallel clause.
!$omp parallel do private(ilocal,ilma,ia,duVpsi,j,ix,iy,iz,w) reduction(+:F_tmp)
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
             w(1) = gtpsi(1,ix,iy,iz) + zI* kAc(1) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
             w(2) = gtpsi(2,ix,iy,iz) + zI* kAc(2) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
             w(3) = gtpsi(3,ix,iy,iz) + zI* kAc(3) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
             duVpsi(1) = duVpsi(1) + conjg(ppg%zekr_uV(j,ilma,ik)) * w(1) ! < uV | exp(ikr) (nabla) | psi >
             duVpsi(2) = duVpsi(2) + conjg(ppg%zekr_uV(j,ilma,ik)) * w(2) ! < uV | exp(ikr) (nabla) | psi >
             duVpsi(3) = duVpsi(3) + conjg(ppg%zekr_uV(j,ilma,ik)) * w(3) ! < uV | exp(ikr) (nabla) | psi >
          end do
          F_tmp(1,ia) = F_tmp(1,ia) - rtmp * dble( conjg(duVpsi(1)) * uVpsibox2(ispin,io,ik,im,ilma) )
          F_tmp(2,ia) = F_tmp(2,ia) - rtmp * dble( conjg(duVpsi(2)) * uVpsibox2(ispin,io,ik,im,ilma) )
          F_tmp(3,ia) = F_tmp(3,ia) - rtmp * dble( conjg(duVpsi(3)) * uVpsibox2(ispin,io,ik,im,ilma) )
       end do
#ifndef __NVCOMPILER_LLVM__ 
! FIXME: NVIDIA compiler crashes with nested omp parallel clause.
!$omp end parallel do
#endif
       end if
       !call timer_end(LOG_CALC_FORCE_NONLOCAL)

       if( PLUS_U_ON )then
          do ilma=1,Nlma_ao
             ia = ppg%ia_tbl_ao(ilma)
             dphipsi = zero
             do j=1,ppg%mps_ao(ia)
                ix = ppg%jxyz_ao(1,j,ia)
                iy = ppg%jxyz_ao(2,j,ia)
                iz = ppg%jxyz_ao(3,j,ia)
                w  = gtpsi(:,ix,iy,iz) + zI* kAc(:) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
                dphipsi(:) = dphipsi(:) + conjg(ppg%zekr_phi_ao(j,ilma,ik)) * w(:)
             end do
!             dphipsi_lma(:,ilma) = dphipsi(:) * system%Hvol
          end do
          Nproj_pairs = size(ppg%proj_pairs_ao,2)
          do iprj=1,Nproj_pairs
             ilma=ppg%proj_pairs_ao(1,iprj)
             jlma=ppg%proj_pairs_ao(2,iprj)
             ia = ppg%proj_pairs_info_ao(1,iprj)
             l  = ppg%proj_pairs_info_ao(2,iprj)
             n  = ppg%proj_pairs_info_ao(3,iprj)
             m1 = ppg%proj_pairs_info_ao(4,iprj)
             m2 = ppg%proj_pairs_info_ao(5,iprj)
!             ddm_mms_nla(:)= & ! ddm_mms_nla(m1,m2,ispin,n,l,ia)=
!                  system%rocc(io,ispin,ik)*system%wtk(ik) &
!                  *( dphipsi_lma(:,ilma)*conjg(phipsibox2(jlma,ispin,io,ik)) &
!                  + phipsibox2(ilma,ispin,io,ik)*conjg(dphipsi_lma(:,jlma)) )
             if( m1 == m2 )then
!                zF_tmp(:,ia) = zF_tmp(:,ia) &
!                     - 0.5d0*U_eff(n,l,ia)*( 1.0d0 - 2.0d0*dm_mms_nla(m1,m2,ispin,n,l,ia) ) &
!                     * ddm_mms_nla(:)
             else
!                zF_tmp(:,ia) = zF_tmp(:,ia) &
!                     - 0.5d0*U_eff(n,l,ia)*( -2.0d0*dm_mms_nla(m1,m2,ispin,n,l,ia) ) &
!                     * ddm_mms_nla(:)
             end if
          end do !iprj
       end if

    end do !ispin
    end do !io
    end do !ik
#ifndef USE_OPENACC
!$omp end do
#endif

    if (allocated(gtpsi)) deallocate(gtpsi)

#ifndef USE_OPENACC
!$omp end parallel
#endif
#ifdef USE_OPENACC
    end if
#endif
    call timer_end(LOG_CALC_FORCE_ELEC_ION)

    if(yn_periodic=='n') then
    ! local part (based on density gradient)
#ifdef USE_OPENACC
!$acc parallel loop private(iz,iy,ix,ia)
#else
!$omp parallel do private(iz,iy,ix,ia)
#endif
      do ia=1,nion
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          F_tmp(:,ia) = F_tmp(:,ia) - dden(:,ix,iy,iz) * ppg%Vpsl_ion(ix,iy,iz,ia)
        end do
        end do
        end do
      end do
#ifdef USE_OPENACC
!$acc end loop
#else
!$omp end parallel do
#endif
      deallocate(dden)
    end if

    !do ia=1,nion
    !  write(*,'(1x,i4,2f20.10)') ia,real(zF_tmp(1,ia)),aimag(zF_tmp(1,ia))
    !  write(*,'(1x,4x,2f20.10)')    real(zF_tmp(2,ia)),aimag(zF_tmp(2,ia))
    !  write(*,'(1x,4x,2f20.10)')    real(zF_tmp(3,ia)),aimag(zF_tmp(3,ia))
    !end do

    call comm_summation(F_tmp,F_sum,3*nion,info%icomm_rko)

#ifdef USE_OPENACC
!$acc parallel loop private(ia)
#else
!$omp parallel do private(ia)
#endif
    do ia=1,nion
      system%Force(:,ia) = system%Force(:,ia) + F_sum(:,ia)
    end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$omp end parallel do
#endif
!
    if (use_symmetry) then
      call sym_vector_force_xyz( system%Force, system%Rion )
    end if

    if(allocated(tpsi%rwf)) deallocate(tpsi%zwf)
    if(allocated(gtpsi)) deallocate(gtpsi)
    if(allocated(uVpsibox)) deallocate(uVpsibox)
    deallocate(F_tmp,F_sum,uVpsibox2)
    if( PLUS_U_ON )then
!      deallocate( phipsibox, phipsibox2 )
!      deallocate( dphipsi_lma )
!      deallocate( zF_tmp )
    end if 

    call timer_end(LOG_CALC_ION_FORCE)
    return
  end subroutine calc_force

#ifdef USE_OPENACC
!===================================================================================================================================

  ! Nonlocal part of the electron-ion force for all local orbitals at once (OpenACC builds).
  ! One device thread per (orbital, local projector) forms B^T (nabla) psi at the projector's grid points
  ! with the stencil of calc_gradient_psi, so no gradient of a whole orbital is stored, and accumulates
  ! < uV | exp(ikr) (nabla + i(k+A/c)) | psi >. The pair contributions are summed over orbitals per projector
  ! on the device, and over the projectors of each atom on the host.
  subroutine calc_force_nonlocal_gpu(system,info,mg,stencil,ppg,zwf,uVpsibox2,F_tmp)
    use structures, only: s_dft_system, s_parallel_info, s_rgrid, s_stencil, s_pp_grid
    use math_constants, only: zi
    use salmon_global, only: yn_periodic
    implicit none
    type(s_dft_system)   ,intent(in)    :: system
    type(s_parallel_info),intent(in)    :: info
    type(s_rgrid)        ,intent(in)    :: mg
    type(s_stencil)      ,intent(in)    :: stencil
    type(s_pp_grid)      ,intent(in)    :: ppg
    complex(8)           ,intent(in)    :: zwf(mg%is_array(1):mg%ie_array(1) &
                                              ,mg%is_array(2):mg%ie_array(2) &
                                              ,mg%is_array(3):mg%ie_array(3) &
                                              ,system%nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e)
    complex(8)           ,intent(in)    :: uVpsibox2(system%nspin,info%io_s:info%io_e,info%ik_s:info%ik_e &
                                                    ,info%im_s:info%im_e,ppg%nlma)
    real(8)              ,intent(inout) :: F_tmp(3,system%nion)
    !
    integer :: ik,io,ispin,im,ik_s,ik_e,io_s,io_e,nspin,nlocal,ilocal,ilma,ia,j,ix,iy,iz,ic
    logical :: if_kAc
    real(8) :: nab(4,3),bmat(3,3),vec_Ac(3),hvol,rtmp,kAc1,kAc2,kAc3,fsum
    complex(8) :: psi,g1,g2,g3,w1,w2,w3,duVpsi1,duVpsi2,duVpsi3
    real(8),allocatable :: F_lma(:,:,:,:,:),F_ilocal(:,:)

    im = info%im_s
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e
    nspin = system%nspin
    nlocal = ppg%ilocal_nlma
    nab = stencil%coef_nab
    bmat = system%rmatrix_B
    vec_Ac = system%vec_Ac
    hvol = system%Hvol
    if_kAc = (yn_periodic=='y')

    allocate(F_lma(3,nlocal,nspin,io_s:io_e,ik_s:ik_e),F_ilocal(3,nlocal))

!$acc data create(F_lma) copyout(F_ilocal) copyin(nab,bmat,vec_Ac)

!$acc parallel loop collapse(4) gang vector &
!$acc   private(ik,io,ispin,ilocal,ilma,ia,j,ix,iy,iz,rtmp,kAc1,kAc2,kAc3,psi,g1,g2,g3,w1,w2,w3) &
!$acc   private(duVpsi1,duVpsi2,duVpsi3)
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,nspin
    do ilocal=1,nlocal
      ilma = ppg%ilocal_nlma2ilma(ilocal)
      ia   = ppg%ilocal_nlma2ia  (ilocal)
      if(if_kAc) then
        kAc1 = system%vec_k(1,ik) + vec_Ac(1)
        kAc2 = system%vec_k(2,ik) + vec_Ac(2)
        kAc3 = system%vec_k(3,ik) + vec_Ac(3)
      else
        kAc1 = 0d0
        kAc2 = 0d0
        kAc3 = 0d0
      end if
      rtmp = 2d0 * system%rocc(io,ik,ispin) * system%wtk(ik) * hvol
      duVpsi1 = 0d0
      duVpsi2 = 0d0
      duVpsi3 = 0d0
!$acc loop seq
      do j=1,ppg%mps(ia)
        ix = ppg%jxyz(1,j,ia)
        iy = ppg%jxyz(2,j,ia)
        iz = ppg%jxyz(3,j,ia)
        ! (nabla) psi at (ix,iy,iz), as in calc_gradient_psi
        g1 = nab(1,1)*(zwf(mg%idx(ix+1),iy,iz,ispin,io,ik,im) - zwf(mg%idx(ix-1),iy,iz,ispin,io,ik,im)) &
           + nab(2,1)*(zwf(mg%idx(ix+2),iy,iz,ispin,io,ik,im) - zwf(mg%idx(ix-2),iy,iz,ispin,io,ik,im)) &
           + nab(3,1)*(zwf(mg%idx(ix+3),iy,iz,ispin,io,ik,im) - zwf(mg%idx(ix-3),iy,iz,ispin,io,ik,im)) &
           + nab(4,1)*(zwf(mg%idx(ix+4),iy,iz,ispin,io,ik,im) - zwf(mg%idx(ix-4),iy,iz,ispin,io,ik,im))
        g2 = nab(1,2)*(zwf(ix,mg%idy(iy+1),iz,ispin,io,ik,im) - zwf(ix,mg%idy(iy-1),iz,ispin,io,ik,im)) &
           + nab(2,2)*(zwf(ix,mg%idy(iy+2),iz,ispin,io,ik,im) - zwf(ix,mg%idy(iy-2),iz,ispin,io,ik,im)) &
           + nab(3,2)*(zwf(ix,mg%idy(iy+3),iz,ispin,io,ik,im) - zwf(ix,mg%idy(iy-3),iz,ispin,io,ik,im)) &
           + nab(4,2)*(zwf(ix,mg%idy(iy+4),iz,ispin,io,ik,im) - zwf(ix,mg%idy(iy-4),iz,ispin,io,ik,im))
        g3 = nab(1,3)*(zwf(ix,iy,mg%idz(iz+1),ispin,io,ik,im) - zwf(ix,iy,mg%idz(iz-1),ispin,io,ik,im)) &
           + nab(2,3)*(zwf(ix,iy,mg%idz(iz+2),ispin,io,ik,im) - zwf(ix,iy,mg%idz(iz-2),ispin,io,ik,im)) &
           + nab(3,3)*(zwf(ix,iy,mg%idz(iz+3),ispin,io,ik,im) - zwf(ix,iy,mg%idz(iz-3),ispin,io,ik,im)) &
           + nab(4,3)*(zwf(ix,iy,mg%idz(iz+4),ispin,io,ik,im) - zwf(ix,iy,mg%idz(iz-4),ispin,io,ik,im))
        psi = zwf(ix,iy,iz,ispin,io,ik,im)
        ! B^T (nabla) psi + i(k+A/c) psi
        w1 = ( bmat(1,1)*g1 + bmat(2,1)*g2 + bmat(3,1)*g3 ) + zI*kAc1*psi
        w2 = ( bmat(1,2)*g1 + bmat(2,2)*g2 + bmat(3,2)*g3 ) + zI*kAc2*psi
        w3 = ( bmat(1,3)*g1 + bmat(2,3)*g2 + bmat(3,3)*g3 ) + zI*kAc3*psi
        duVpsi1 = duVpsi1 + conjg(ppg%zekr_uV(j,ilma,ik)) * w1 ! < uV | exp(ikr) (nabla) | psi >
        duVpsi2 = duVpsi2 + conjg(ppg%zekr_uV(j,ilma,ik)) * w2
        duVpsi3 = duVpsi3 + conjg(ppg%zekr_uV(j,ilma,ik)) * w3
      end do
      F_lma(1,ilocal,ispin,io,ik) = rtmp * dble( conjg(duVpsi1) * uVpsibox2(ispin,io,ik,im,ilma) )
      F_lma(2,ilocal,ispin,io,ik) = rtmp * dble( conjg(duVpsi2) * uVpsibox2(ispin,io,ik,im,ilma) )
      F_lma(3,ilocal,ispin,io,ik) = rtmp * dble( conjg(duVpsi3) * uVpsibox2(ispin,io,ik,im,ilma) )
    end do
    end do
    end do
    end do
!$acc end parallel

    ! sum over the local orbitals, one gang per projector and component
!$acc parallel loop collapse(2) gang private(ilocal,ic,fsum)
    do ilocal=1,nlocal
    do ic=1,3
      fsum = 0d0
!$acc loop vector collapse(3) reduction(+:fsum)
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,nspin
        fsum = fsum + F_lma(ic,ilocal,ispin,io,ik)
      end do
      end do
      end do
      F_ilocal(ic,ilocal) = fsum
    end do
    end do
!$acc end parallel

!$acc end data

    do ilocal=1,nlocal
      ia = ppg%ilocal_nlma2ia(ilocal)
      F_tmp(1,ia) = F_tmp(1,ia) - F_ilocal(1,ilocal)
      F_tmp(2,ia) = F_tmp(2,ia) - F_ilocal(2,ilocal)
      F_tmp(3,ia) = F_tmp(3,ia) - F_ilocal(3,ilocal)
    end do

    deallocate(F_lma,F_ilocal)
  end subroutine calc_force_nonlocal_gpu

!===================================================================================================================================

  ! Density-gradient term of the local force for isolated systems, for all local orbitals at once (OpenACC builds).
  ! One device thread per grid point sums 2 f_n w_k Re[ (B^T (nabla) psi)^* psi ] over the local orbitals in the order
  ! of the orbital loop in calc_force, so the result per grid point is that of the original accumulation.
  subroutine calc_dden_gpu(system,info,mg,stencil,zwf,dden)
    use structures, only: s_dft_system, s_parallel_info, s_rgrid, s_stencil
    implicit none
    type(s_dft_system)   ,intent(in)    :: system
    type(s_parallel_info),intent(in)    :: info
    type(s_rgrid)        ,intent(in)    :: mg
    type(s_stencil)      ,intent(in)    :: stencil
    complex(8)           ,intent(in)    :: zwf(mg%is_array(1):mg%ie_array(1) &
                                              ,mg%is_array(2):mg%ie_array(2) &
                                              ,mg%is_array(3):mg%ie_array(3) &
                                              ,system%nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e)
    real(8)              ,intent(inout) :: dden(3,mg%is_array(1):mg%ie_array(1) &
                                               ,mg%is_array(2):mg%ie_array(2) &
                                               ,mg%is_array(3):mg%ie_array(3))
    !
    integer :: ik,io,ispin,im,ik_s,ik_e,io_s,io_e,nspin,ix,iy,iz,ixs,ixe,iys,iye,izs,ize
    real(8) :: nab(4,3),bmat(3,3),hvol,rtmp,d1,d2,d3
    complex(8) :: psi,g1,g2,g3

    im = info%im_s
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e
    nspin = system%nspin
    ixs = mg%is(1)
    ixe = mg%ie(1)
    iys = mg%is(2)
    iye = mg%ie(2)
    izs = mg%is(3)
    ize = mg%ie(3)
    nab = stencil%coef_nab
    bmat = system%rmatrix_B
    hvol = system%Hvol

!$acc parallel loop collapse(3) gang vector copyin(nab,bmat) &
!$acc   private(ix,iy,iz,ik,io,ispin,rtmp,psi,g1,g2,g3,d1,d2,d3)
    do iz=izs,ize
    do iy=iys,iye
    do ix=ixs,ixe
      d1 = dden(1,ix,iy,iz)
      d2 = dden(2,ix,iy,iz)
      d3 = dden(3,ix,iy,iz)
!$acc loop seq
      do ik=ik_s,ik_e
!$acc loop seq
      do io=io_s,io_e
!$acc loop seq
      do ispin=1,nspin
        rtmp = 2d0 * system%rocc(io,ik,ispin) * system%wtk(ik) * hvol
        g1 = nab(1,1)*(zwf(mg%idx(ix+1),iy,iz,ispin,io,ik,im) - zwf(mg%idx(ix-1),iy,iz,ispin,io,ik,im)) &
           + nab(2,1)*(zwf(mg%idx(ix+2),iy,iz,ispin,io,ik,im) - zwf(mg%idx(ix-2),iy,iz,ispin,io,ik,im)) &
           + nab(3,1)*(zwf(mg%idx(ix+3),iy,iz,ispin,io,ik,im) - zwf(mg%idx(ix-3),iy,iz,ispin,io,ik,im)) &
           + nab(4,1)*(zwf(mg%idx(ix+4),iy,iz,ispin,io,ik,im) - zwf(mg%idx(ix-4),iy,iz,ispin,io,ik,im))
        g2 = nab(1,2)*(zwf(ix,mg%idy(iy+1),iz,ispin,io,ik,im) - zwf(ix,mg%idy(iy-1),iz,ispin,io,ik,im)) &
           + nab(2,2)*(zwf(ix,mg%idy(iy+2),iz,ispin,io,ik,im) - zwf(ix,mg%idy(iy-2),iz,ispin,io,ik,im)) &
           + nab(3,2)*(zwf(ix,mg%idy(iy+3),iz,ispin,io,ik,im) - zwf(ix,mg%idy(iy-3),iz,ispin,io,ik,im)) &
           + nab(4,2)*(zwf(ix,mg%idy(iy+4),iz,ispin,io,ik,im) - zwf(ix,mg%idy(iy-4),iz,ispin,io,ik,im))
        g3 = nab(1,3)*(zwf(ix,iy,mg%idz(iz+1),ispin,io,ik,im) - zwf(ix,iy,mg%idz(iz-1),ispin,io,ik,im)) &
           + nab(2,3)*(zwf(ix,iy,mg%idz(iz+2),ispin,io,ik,im) - zwf(ix,iy,mg%idz(iz-2),ispin,io,ik,im)) &
           + nab(3,3)*(zwf(ix,iy,mg%idz(iz+3),ispin,io,ik,im) - zwf(ix,iy,mg%idz(iz-3),ispin,io,ik,im)) &
           + nab(4,3)*(zwf(ix,iy,mg%idz(iz+4),ispin,io,ik,im) - zwf(ix,iy,mg%idz(iz-4),ispin,io,ik,im))
        psi = zwf(ix,iy,iz,ispin,io,ik,im)
        d1 = d1 + rtmp * dble( conjg( bmat(1,1)*g1 + bmat(2,1)*g2 + bmat(3,1)*g3 ) * psi )
        d2 = d2 + rtmp * dble( conjg( bmat(1,2)*g1 + bmat(2,2)*g2 + bmat(3,2)*g3 ) * psi )
        d3 = d3 + rtmp * dble( conjg( bmat(1,3)*g1 + bmat(2,3)*g2 + bmat(3,3)*g3 ) * psi )
      end do
      end do
      end do
      dden(1,ix,iy,iz) = d1
      dden(2,ix,iy,iz) = d2
      dden(3,ix,iy,iz) = d3
    end do
    end do
    end do
!$acc end parallel

  end subroutine calc_dden_gpu
#endif
  
!===================================================================================================================================

  subroutine force_ewald_rspace(F_sum,F_tmp,system,info,ewald,pp,nion,comm)
    use structures
    use math_constants,only : pi,zi
    use salmon_math
    use salmon_global, only: yn_periodic,kion,aEwald,cutoff_r
    use communication, only: comm_summation
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_ewald_ion_ion),intent(in) :: ewald
    type(s_pp_info)   ,intent(in) :: pp
    integer  ,intent(in) :: nion,comm
    real(8)              :: F_tmp(3,nion),F_tmp_l(3,nion),F_sum(3,nion)
    !
    integer :: ix,iy,iz,iia,ia,ib,ipair
    integer :: nion_mg_l,npair_l ! hoisted: not device-resident as dummy-arg-type members
    real(8) :: rr,rab(3),r(3)

    select case(yn_periodic)
    case('n')

      F_sum = 0d0
      do ia=1,nion
        do ib=1,nion
          if(ia==ib) cycle
          rab(:) = system%Rion(:,ia) - system%Rion(:,ib)
          rr = sum(rab(:)**2)
          F_sum(:,ia) = F_sum(:,ia) + pp%Zps(Kion(ia))*pp%Zps(Kion(ib)) * rab(:)/sqrt(rr)**3
        end do
      end do

    case('y')

      F_sum = 0d0
      F_tmp = 0d0
      if(ewald%yn_bookkeep=='y') then

         F_tmp_l = 0d0
         nion_mg_l = info%nion_mg
#ifdef USE_OPENACC
!$acc kernels loop private(iia,ia,ipair,ix,iy,iz,ib,r,rab,rr,npair_l) copyin(ewald, pp)
#else
!$omp parallel do private(iia,ia,ipair,ix,iy,iz,ib,r,rab,rr,npair_l)
#endif
         do iia= 1,nion_mg_l
        !do ia= system%nion_s, system%nion_e
        !do ia=1,nion
            ia = info%ia_mg(iia)
            npair_l = ewald%npair_bk(iia)
            do ipair = 1,npair_l
               ix = ewald%bk(1,ipair,iia)
               iy = ewald%bk(2,ipair,iia)
               iz = ewald%bk(3,ipair,iia)
               ib = ewald%bk(4,ipair,iia)
              !if (ix**2+iy**2+iz**2 == 0 .and. ia == ib) cycle
               r(1) = ix*system%primitive_a(1,1) &
                    + iy*system%primitive_a(1,2) &
                    + iz*system%primitive_a(1,3)
               r(2) = ix*system%primitive_a(2,1) &
                    + iy*system%primitive_a(2,2) &
                    + iz*system%primitive_a(2,3)
               r(3) = ix*system%primitive_a(3,1) &
                    + iy*system%primitive_a(3,2) &
                    + iz*system%primitive_a(3,3)
               rab(1) = system%Rion(1,ia)-r(1) - system%Rion(1,ib)
               rab(2) = system%Rion(2,ia)-r(2) - system%Rion(2,ib)
               rab(3) = system%Rion(3,ia)-r(3) - system%Rion(3,ib)
               rr = sum(rab(:)**2)
               if(rr .gt. cutoff_r**2) then
                 !!!
               else
                 F_tmp_l(:,ia) = F_tmp_l(:,ia)  &
                               - pp%Zps(Kion(ia))*pp%Zps(Kion(ib))*rab(:)/sqrt(rr)*(-erfc_salmon(sqrt(aEwald*rr))/rr &
                               -2*sqrt(aEwald/(rr*Pi))*exp(-aEwald*rr))
               end if

            end do  !ipair
         end do     !ia
#ifdef USE_OPENACC
!$acc end kernels
#else
!$omp end parallel do
#endif
         call comm_summation(F_tmp_l,F_tmp,3*nion,comm)

      else

         stop  "must use book-keeping method for periodic condition"

      endif
      F_sum = F_sum + F_tmp
      
    end select
  end subroutine force_ewald_rspace
end module force_sub
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
