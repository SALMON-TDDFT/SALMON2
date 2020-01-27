!
!  Copyright 2019 SALMON developers
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

  subroutine calc_force(system,pp,fg,info,mg,stencil,srg,ppg,tpsi,ewald)
    use structures
    use math_constants,only : zi
    use stencil_sub, only: calc_gradient_psi
    use sendrecv_grid, only: s_sendrecv_grid, update_overlap_real8, update_overlap_complex8, dealloc_cache
    use communication, only: comm_summation
    use nonlocal_potential, only: calc_uVpsi_rdivided
    use sym_vector_sub, only: sym_vector_xyz
    use plusU_global, only: PLUS_U_ON, dm_mms_nla, U_eff
    use timer
    implicit none
    type(s_dft_system)      ,intent(inout) :: system
    type(s_pp_info)         ,intent(in) :: pp
    type(s_reciprocal_grid) ,intent(in) :: fg
    type(s_orbital_parallel),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid)      :: srg
    type(s_pp_grid),intent(in) :: ppg
    type(s_orbital)            :: tpsi
    type(s_ewald_ion_ion),intent(in) :: ewald
    !
    integer :: ix,iy,iz,ia,nion,im,Nspin,ik_s,ik_e,io_s,io_e,nlma,ik,io,ispin,ilma,j
    integer :: m1,m2,jlma,n,l,Nproj_pairs,iprj,Nlma_ao
    real(8) :: kAc(3), rtmp, rtmp2(3)
    real(8),allocatable :: F_tmp(:,:),F_sum(:,:)
    real(8),allocatable :: dden(:,:,:,:)
    complex(8) :: w(3),duVpsi(3)
    complex(8),allocatable :: gtpsi(:,:,:,:),uVpsibox(:,:,:,:,:),uVpsibox2(:,:,:,:,:)
    complex(8),allocatable :: phipsibox(:,:),phipsibox2(:,:)
    complex(8),allocatable :: dphipsi_lma(:,:)
    complex(8) :: ddm_mms_nla(3), phipsi, dphipsi(3)
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8),allocatable :: zF_tmp(:,:)
    integer :: Norb,iorb

    call timer_begin(LOG_CALC_ION_FORCE)

    nion = system%nion
    if(.not.allocated(system%Force)) allocate(system%Force(3,nion))
    allocate( F_tmp(3,nion), F_sum(3,nion) )
    if( PLUS_U_ON ) then
      allocate( zF_tmp(3,nion) )
      zF_tmp=zero
    end if


    if(info%im_s/=1 .or. info%im_e/=1) stop "error: calc_force_periodic" !??????
    im = 1

    Nspin = system%Nspin
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e
    Norb = system%Nspin*info%numo*info%numk

    Nlma = ppg%Nlma

  ! Ewald sum of ion-ion interaction
  call timer_begin(LOG_CALC_FORCE_ION_ION)
    call force_ion_ion(F_sum,F_tmp,system,ewald,pp,fg,nion)
    system%Force = F_sum
  call timer_end(LOG_CALC_FORCE_ION_ION)

  ! electron-ion interaction
  call timer_begin(LOG_CALC_FORCE_ELEC_ION)
    if(allocated(tpsi%rwf)) then
      allocate(tpsi%zwf(mg%is_array(1):mg%ie_array(1) &
                       ,mg%is_array(2):mg%ie_array(2) &
                       ,mg%is_array(3):mg%ie_array(3) &
                       ,nspin,info%io_s:info%io_e,info%ik_s:info%ik_e,info%im_s:info%im_e))
      tpsi%zwf = cmplx(tpsi%rwf)
    end if

    allocate(gtpsi(3,mg%is_array(1):mg%ie_array(1) &
                    ,mg%is_array(2):mg%ie_array(2) &
                    ,mg%is_array(3):mg%ie_array(3)))

  ! uVpsibox2 = < uV | exp(ikr) | psi >
    call calc_uVpsi_rdivided(nspin,info,ppg,tpsi,uVpsibox,uVpsibox2)

    if( PLUS_U_ON )then
      Nlma_ao = size(ppg%ia_tbl_ao)
      allocate( phipsibox(Nlma_ao,Norb)  ); phipsibox=zero
      allocate( phipsibox2(Nlma_ao,Norb) ); phipsibox2=zero
      iorb = 0
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin
        iorb = iorb + 1
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
          phipsi = phipsi * ppg%Hvol
          phipsibox(ilma,iorb) = phipsi
        end do
      end do
      end do
      end do
      call comm_summation(phipsibox,phipsibox2,Nlma_ao*Norb,info%icomm_r)
    end if
  call timer_end(LOG_CALC_FORCE_ELEC_ION)

    if(info%if_divide_rspace) then
       call update_overlap_complex8(srg, mg, tpsi%zwf)
    end if

  call timer_begin(LOG_CALC_FORCE_ELEC_ION)
    kAc   = 0d0
    F_tmp = 0d0
    allocate( dden(3,mg%is_array(1):mg%ie_array(1) &
                    ,mg%is_array(2):mg%ie_array(2) &
                    ,mg%is_array(3):mg%ie_array(3)))

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz = mg%is_array(3),mg%ie_array(3)
    do iy = mg%is_array(2),mg%ie_array(2)
    do ix = mg%is_array(1),mg%ie_array(1)
       dden(:,ix,iy,iz) = 0d0
    enddo
    enddo
    enddo
!$omp end parallel do

    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

       ! gtpsi = (nabla) psi
       call timer_begin(LOG_CALC_FORCE_GTPSI)
       call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,io,ik,im),gtpsi,mg%is_array,mg%ie_array,mg%is,mg%ie &
            ,mg%idx,mg%idy,mg%idz,stencil%coef_nab,system%rmatrix_B)
       call timer_end(LOG_CALC_FORCE_GTPSI)


       call timer_begin(LOG_CALC_FORCE_DDEN)
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
       call timer_end(LOG_CALC_FORCE_DDEN)

       ! nonlocal part
       call timer_begin(LOG_CALC_FORCE_NONLOCAL)
       if(system%iperiodic==3) kAc(1:3) = system%vec_k(1:3,ik) + system%vec_Ac(1:3)
       rtmp = 2d0 * system%rocc(io,ik,ispin) * system%wtk(ik) * system%Hvol

!$omp parallel do private(ilma,ia,duVpsi,j,ix,iy,iz,w) reduction(+:F_tmp)
       do ilma=1,Nlma
          ia = ppg%ia_tbl(ilma)
          duVpsi = 0d0
          do j=1,ppg%mps(ia)
             ix = ppg%jxyz(1,j,ia)
             iy = ppg%jxyz(2,j,ia)
             iz = ppg%jxyz(3,j,ia)
             w = gtpsi(:,ix,iy,iz) + zI* kAc(:) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
             duVpsi = duVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * w ! < uV | exp(ikr) (nabla) | psi >
          end do
          F_tmp(:,ia) = F_tmp(:,ia)  &
                     - rtmp * dble( conjg(duVpsi(:)) * uVpsibox2(ispin,io,ik,im,ilma) ) 
       end do
!$omp end parallel do
       call timer_end(LOG_CALC_FORCE_NONLOCAL)

       if( PLUS_U_ON )then
          if( .not.allocated(dphipsi_lma) )then
             allocate( dphipsi_lma(3,Nlma_ao) ); dphipsi_lma=zero
          end if
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
             dphipsi_lma(:,ilma) = dphipsi(:) * ppg%Hvol
          end do
       end if
       if( PLUS_U_ON )then
          Nproj_pairs = size(ppg%proj_pairs_ao,2)
          do iprj=1,Nproj_pairs
             ilma=ppg%proj_pairs_ao(1,iprj)
             jlma=ppg%proj_pairs_ao(2,iprj)
             ia = ppg%proj_pairs_info_ao(1,iprj)
             l  = ppg%proj_pairs_info_ao(2,iprj)
             n  = ppg%proj_pairs_info_ao(3,iprj)
             m1 = ppg%proj_pairs_info_ao(4,iprj)
             m2 = ppg%proj_pairs_info_ao(5,iprj)
             ddm_mms_nla(:)= & ! ddm_mms_nla(m1,m2,ispin,n,l,ia)=
                  system%rocc(io,ispin,ik)*system%wtk(ik) &
                  *( dphipsi_lma(:,ilma)*conjg(phipsibox2(jlma,iorb)) &
                  + phipsibox2(ilma,iorb)*conjg(dphipsi_lma(:,jlma)) )
             if( m1 == m2 )then
                zF_tmp(:,ia) = zF_tmp(:,ia) &
                     - 0.5d0*U_eff(n,l,ia)*( 1.0d0 - 2.0d0*dm_mms_nla(m1,m2,ispin,n,l,ia) ) &
                     * ddm_mms_nla(:)
             else
                zF_tmp(:,ia) = zF_tmp(:,ia) &
                     - 0.5d0*U_eff(n,l,ia)*( -2.0d0*dm_mms_nla(m1,m2,ispin,n,l,ia) ) &
                     * ddm_mms_nla(:)
             end if
          end do !iprj
       end if

    end do !ispin
    end do !io
    end do !ik
  call timer_end(LOG_CALC_FORCE_ELEC_ION)

    ! local part (based on density gradient)
  call timer_begin(LOG_CALC_FORCE_LOCAL)
!$omp parallel do private(iz,iy,ix,ia)
    do ia=1,nion
       do iz=mg%is(3),mg%ie(3)
       do iy=mg%is(2),mg%ie(2)
       do ix=mg%is(1),mg%ie(1)
          F_tmp(:,ia) = F_tmp(:,ia) - dden(:,ix,iy,iz) * ppg%Vpsl_atom(ix,iy,iz,ia)
       end do
       end do
       end do
    end do
!$omp end parallel do
  call timer_end(LOG_CALC_FORCE_LOCAL)


    !do ia=1,nion
    !  write(*,'(1x,i4,2f20.10)') ia,real(zF_tmp(1,ia)),aimag(zF_tmp(1,ia))
    !  write(*,'(1x,4x,2f20.10)')    real(zF_tmp(2,ia)),aimag(zF_tmp(2,ia))
    !  write(*,'(1x,4x,2f20.10)')    real(zF_tmp(3,ia)),aimag(zF_tmp(3,ia))
    !end do
    call comm_summation(F_tmp,F_sum,3*nion,info%icomm_rko)
    system%Force = system%Force + F_sum
!
    do ia=1,nion
       call sym_vector_xyz( system%Force(:,ia) )
    end do

    if(allocated(tpsi%rwf)) deallocate(tpsi%zwf)
    deallocate(F_tmp,F_sum,gtpsi,uVpsibox,uVpsibox2)
    if( PLUS_U_ON )then
      deallocate( phipsibox, phipsibox2 )
      deallocate( dphipsi_lma )
      deallocate( zF_tmp )
    end if 

    call timer_end(LOG_CALC_ION_FORCE)
    return
  end subroutine calc_force

  subroutine force_ion_ion(F_sum,F_tmp,system,ewald,pp,fg,nion)
    use structures
    use math_constants,only : pi,zi
    use salmon_math
    use salmon_global, only: kion,NEwald,aEwald
    use communication, only: comm_summation
    use inputoutput, only: cutoff_r, cutoff_g
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_ewald_ion_ion),intent(in) :: ewald
    type(s_pp_info)   ,intent(in) :: pp
    type(s_reciprocal_grid),intent(in) :: fg
    integer  ,intent(in) :: nion
    real(8)              :: F_tmp(3,nion),F_tmp_l(3,nion),F_sum(3,nion)
    !
    integer :: ix,iy,iz,ia,ib,ig,ipair
    real(8) :: rr,rab(3),r(3),g(3),G2,Gd, rtmp(3)
    complex(8) :: rho_i, ctmp1(3), ctmp2

    select case(system%iperiodic)
    case(0)

      F_sum = 0d0
      do ia=1,nion
        do ib=1,nion
          if(ia==ib) cycle
          rab(:) = system%Rion(:,ia) - system%Rion(:,ib)
          rr = sum(rab(:)**2)
          F_sum(:,ia) = F_sum(:,ia) + pp%Zps(Kion(ia))*pp%Zps(Kion(ib)) * rab(:)/sqrt(rr)**3
        end do
      end do

    case(3)


      F_tmp   = 0d0
!$omp parallel do private(ig,ia,r,g,G2,Gd,rho_i,rtmp,ctmp1,ctmp2) reduction(+:F_tmp)
      do ig=fg%ig_s,fg%ig_e
         if(ig == fg%iGzero ) cycle
         g(1) = fg%Gx(ig)
         g(2) = fg%Gy(ig)
         g(3) = fg%Gz(ig)
         G2   = sum(g(:)**2)
         if(G2 .gt. cutoff_g**2) cycle   !xxx

         rho_i= fg%zrhoG_ion(ig)
         rtmp(:) = 0.5d0 * g(:) * (4*Pi/G2) * exp(-G2/(4*aEwald))
         ctmp1(:)= rtmp(:) * zI

         do ia=1,nion
            r = system%Rion(1:3,ia)
            Gd = sum(g(:)*r(:))
            ctmp2 = rho_i*exp(zI*Gd)
            F_tmp(:,ia) = F_tmp(:,ia)  &
                        + pp%Zps(Kion(ia)) * ctmp1(:) * (conjg(ctmp2)-ctmp2)
         end do
      end do
!$omp end parallel do
      call comm_summation(F_tmp,F_sum,3*nion,fg%icomm_G)


      F_tmp = 0d0
      if(ewald%yn_bookkeep=='y') then

         F_tmp_l = 0d0
!$omp parallel do private(ia,ipair,ix,iy,iz,ib,r,rab,rr)
         do ia= system%nion_s, system%nion_e
        !do ia=1,nion
            do ipair = 1,ewald%npair_bk(ia)
               ix = ewald%bk(1,ipair,ia)
               iy = ewald%bk(2,ipair,ia)
               iz = ewald%bk(3,ipair,ia)
               ib = ewald%bk(4,ipair,ia)
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
               if(rr .gt. cutoff_r**2) cycle
               F_tmp_l(:,ia) = F_tmp_l(:,ia)  &
                             - pp%Zps(Kion(ia))*pp%Zps(Kion(ib))*rab(:)/sqrt(rr)*(-erfc_salmon(sqrt(aEwald*rr))/rr &
                             -2*sqrt(aEwald/(rr*Pi))*exp(-aEwald*rr))

            end do  !ipair
         end do     !ia
!$omp end parallel do
         call comm_summation(F_tmp_l,F_tmp,3*nion,fg%icomm_G)

      else

         do ia=1,nion
            do ix=-NEwald,NEwald
            do iy=-NEwald,NEwald
            do iz=-NEwald,NEwald
            do ib=1,nion
               if (ix**2+iy**2+iz**2 == 0 .and. ia == ib) cycle
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
               F_tmp(:,ia) = F_tmp(:,ia) - pp%Zps(Kion(ia))*pp%Zps(Kion(ib))*rab(:)/sqrt(rr)*(-erfc_salmon(sqrt(aEwald*rr))/rr &
                                        -2*sqrt(aEwald/(rr*Pi))*exp(-aEwald*rr))
            end do
            end do
            end do
            end do
         end do

      endif
      F_sum = F_sum + F_tmp
      
    end select
  end subroutine force_ion_ion
end module force_sub
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
