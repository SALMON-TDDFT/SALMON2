!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!
! Macroscopic dielectric function / absorption from the full-frequency
! dielectric matrix (spec-b1).  The macroscopic (local-field-included)
! dielectric is the head of the inverse dielectric matrix:
!
!   eps_M(q->0, w) = 1 / [ eps^{-1}(q->0, w) ]_{G=0, G'=0}
!
! and Im eps_M(w) is the RPA absorption spectrum (to be compared with the
! RT-TDDFT delta-kick Im eps; the difference is the xc kernel, RPA vs ALDA).
! epsinv_w comes from calc_w_freq evaluated at the smallest-|q| proxy for
! q->0 (StageA: finite-q head; StageB will use the velocity matrix element).
!
! No proper nouns appear in this file per the project constraint.
!
module gw_absorption_sub
  implicit none
  private
  public :: calc_absorption
  public :: calc_absorption_velocity

  real(8), parameter :: au_to_eV = 27.211386d0
  real(8), parameter :: fpi      = 16.0d0*atan(1.0d0)   ! 4*pi

contains

  ! --------------------------------------------------------------------
  ! eps_M(w) = 1 / epsinv_w(ig0,ig0,iw), and (on root) a text dump of the
  ! spectrum to <sysname>_absorption.data with columns
  !   omega[eV]   Re eps_M   Im eps_M
  ! ig0 is the G=0 index in the dielectric G-list; omega_grid is in a.u.
  ! --------------------------------------------------------------------
  subroutine calc_absorption(ng, ig0, nomega, omega_grid, epsinv_w, eps_macro, &
                             is_root, sysname)
    implicit none
    integer,          intent(in)  :: ng, ig0, nomega
    real(8),          intent(in)  :: omega_grid(nomega)      ! a.u.
    complex(8),       intent(in)  :: epsinv_w(ng,ng,nomega)
    complex(8),       intent(out) :: eps_macro(nomega)
    logical,          intent(in)  :: is_root
    character(*),     intent(in)  :: sysname

    integer :: iw, fh
    complex(8) :: zhead

    do iw = 1, nomega
      zhead = epsinv_w(ig0,ig0,iw)
      if (abs(zhead) < 1.0d-30) then
        eps_macro(iw) = (0.0d0, 0.0d0)
      else
        eps_macro(iw) = (1.0d0, 0.0d0) / zhead
      end if
    end do

    if (is_root) then
      open(newunit=fh, file=trim(sysname)//'_absorption.data', status='replace')
      write(fh,'(A)') "# full-frequency G0W0+RPA macroscopic dielectric (finite-q head proxy)"
      write(fh,'(A)') "# omega[eV]            Re eps_M             Im eps_M"
      do iw = 1, nomega
        write(fh,'(3ES22.12)') omega_grid(iw)*au_to_eV, &
                               dble(eps_macro(iw)), aimag(eps_macro(iw))
      end do
      close(fh)
    end if
  end subroutine calc_absorption


  ! --------------------------------------------------------------------
  ! Velocity-head full-LFE macroscopic dielectric (spec-b1-5a).  Stable
  ! q->0 limit via the head/wing/body decomposition (the finite-q proxy
  ! blows up because eps wings ~ 1/q).  Polarization ê = x (cubic Si is
  ! isotropic).  Per k, with Delta=esp(c)-esp(v)>0, P=(f_v-f_c)*pole(w):
  !   pole(w) = 0.5[1/(w-Delta+i eta) - 1/(w+Delta-i eta)]
  !   vh(cv)  = ê . <c,k|-i nabla|v,k>      (calc_velocity_mtxel)
  !   M_G(cv) = <c,k|e^{iG.r}|v,k>          (calc_mtxel at q=0)
  ! Head  B(w)      = (1/Om) sum |vh/Delta|^2 P
  ! Wing  A_G(w)    = (1/Om) sum conjg(vh/Delta) M_G P     (G/=0)
  !       At_G(w)   = (1/Om) sum conjg(M_G) (vh/Delta) P   (G/=0)
  ! Body  chi0b(w)  = (1/Om) sum conjg(M_G) M_G' P         (G,G'/=0)
  ! Then eps00 = 1 - 4pi B ; eps_body = I - v(G) chi0b ; invert ->
  !   eps_M(w) = eps00 - (4pi)^2 sum_{GG'} A_G (eps_body^-1)_GG' At_G' / |G'|^2
  ! --------------------------------------------------------------------
  subroutine calc_absorption_velocity(system, info, mg, lg, stencil, srg, spsi, &
                                      esp, gvec, gg, ng, ig0, ispin, &
                                      nomega, omega_grid, eta, eps_macro)
    use structures,      only: s_dft_system, s_parallel_info, s_rgrid, s_stencil, &
                               s_sendrecv_grid, s_orbital
    use gw_velocity_sub, only: calc_velocity_mtxel
    use gw_mtxel_sub,    only: calc_mtxel
    use gw_epsilon_sub,  only: build_gminus
    implicit none
    type(s_dft_system),    intent(in)    :: system
    type(s_parallel_info), intent(in)    :: info
    type(s_rgrid),         intent(in)    :: mg, lg
    type(s_stencil),       intent(in)    :: stencil
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_orbital),       intent(inout) :: spsi   ! inout: halo update in velocity mtxel
    real(8),    intent(in)  :: esp(system%no, system%nk, system%nspin)
    integer,    intent(in)  :: ng
    real(8),    intent(in)  :: gvec(3,ng), gg(ng)
    integer,    intent(in)  :: ig0
    integer,    intent(in)  :: ispin
    integer,    intent(in)  :: nomega
    real(8),    intent(in)  :: omega_grid(nomega), eta
    complex(8), intent(out) :: eps_macro(nomega)

    real(8),    parameter   :: occ_thr = 1.0d-6, denom_tiny = 1.0d-8
    complex(8), parameter   :: zi = (0.0d0, 1.0d0)
    complex(8), allocatable :: vmat(:,:,:), mtxel(:,:,:), mcv(:,:)
    complex(8), allocatable :: vhp(:)
    real(8),    allocatable :: dwp(:), delp(:)
    integer,    allocatable :: iGm(:)
    complex(8), allocatable :: Bw(:), Aw(:,:), Atw(:,:), chi0b(:,:,:)
    complex(8), allocatable :: epsb(:,:), zwork(:)
    integer,    allocatable :: ipiv(:), bidx(:)
    integer :: no, nk, ik, iv, ic, ig, jg, ipair, npair, nfill, iw, im
    integer :: nb, a, b, lwork, linfo
    real(8) :: omega, inv_omega, fv, fc, de
    complex(8) :: zw, pole, hv, wgt, zlfe

    no = system%no;  nk = system%nk;  im = 1
    ! BZ k-sum factor 1/(N_k*Omega) (same as the Fock rnk in gw_sigma_x); the
    ! per-cell M and v_cv carry no 1/N_k, so the mesh factor is explicit here.
    omega = abs(system%det_a);  inv_omega = 1.0d0/(dble(nk)*omega)

    allocate(iGm(ng));  call build_gminus(ng, gvec, iGm)
    allocate(Bw(nomega), Aw(ng,nomega), Atw(ng,nomega), chi0b(ng,ng,nomega))
    Bw = (0d0,0d0); Aw = (0d0,0d0); Atw = (0d0,0d0)
    ! NUMA first-touch: zero chi0b (the 0.2GB array) with the SAME OMP-over-iw
    ! static schedule the accumulation uses, so each omega slice's pages land on
    ! the CMG of the thread that writes it (A64FX = 4 CMG x 12 cores).
!$omp parallel do default(shared) private(iw, ig, jg) schedule(static)
    do iw = 1, nomega
      do jg = 1, ng
        do ig = 1, ng
          chi0b(ig,jg,iw) = (0.0d0, 0.0d0)
        end do
      end do
    end do
!$omp end parallel do

    ! ---- accumulate head/wing/body over the BZ at q=0 (ikq=ik) ----
    do ik = 1, nk
      allocate(vmat(3,no,no))
      call calc_velocity_mtxel(system, info, mg, stencil, srg, spsi, ik, ispin, vmat)
      allocate(mtxel(ng,no,no))
      call calc_mtxel(system, info, mg, lg, spsi, gvec, ng, ik, ik, ispin, no, no, mtxel)

      npair = 0
      do iv = 1, no
        if (system%rocc(iv,ik,ispin) <= occ_thr) cycle
        do ic = 1, no
          if (system%rocc(ic,ik,ispin) > occ_thr) cycle
          npair = npair + 1
        end do
      end do

      if (npair > 0) then
        allocate(mcv(ng,npair), vhp(npair), dwp(npair), delp(npair))
        ipair = 0
        do iv = 1, no
          fv = system%rocc(iv,ik,ispin)
          if (fv <= occ_thr) cycle
          do ic = 1, no
            fc = system%rocc(ic,ik,ispin)
            if (fc > occ_thr) cycle
            de = esp(iv,ik,ispin) - esp(ic,ik,ispin)
            if (abs(de) < denom_tiny) cycle
            ipair = ipair + 1
            dwp(ipair)  = fv - fc
            delp(ipair) = -de                       ! Delta = e_c - e_v > 0
            vhp(ipair)  = vmat(1,ic,iv)             ! ê=x : <c|-i d/dx|v>
            do ig = 1, ng                            ! M_G(cv) at q=0 (G->-G map)
              if (iGm(ig) == 0) then
                mcv(ig,ipair) = (0.0d0, 0.0d0)
              else
                mcv(ig,ipair) = mtxel(iGm(ig), ic, iv)
              end if
            end do
          end do
        end do
        nfill = ipair

!$omp parallel do default(shared) schedule(static) &
!$omp   private(iw, zw, ipair, pole, hv, wgt, ig, jg)
        do iw = 1, nomega
          zw = cmplx(omega_grid(iw), 0.0d0, 8)
          do ipair = 1, nfill
            pole = 0.5d0*( 1.0d0/(zw - delp(ipair) + zi*eta) &
                        -  1.0d0/(zw + delp(ipair) - zi*eta) )
            wgt  = inv_omega * dwp(ipair) * pole
            hv   = vhp(ipair) / delp(ipair)          ! ê.v_cv / Delta
            Bw(iw) = Bw(iw) + wgt * (conjg(hv)*hv)
            do ig = 1, ng
              if (ig == ig0) cycle
              Aw(ig,iw)  = Aw(ig,iw)  + wgt * conjg(hv) * mcv(ig,ipair)
              Atw(ig,iw) = Atw(ig,iw) + wgt * conjg(mcv(ig,ipair)) * hv
            end do
            do jg = 1, ng
              if (jg == ig0) cycle
              do ig = 1, ng
                if (ig == ig0) cycle
                chi0b(ig,jg,iw) = chi0b(ig,jg,iw) &
                                + wgt * conjg(mcv(ig,ipair)) * mcv(jg,ipair)
              end do
            end do
          end do
        end do
!$omp end parallel do

        deallocate(mcv, vhp, dwp, delp)
      end if
      deallocate(vmat, mtxel)
    end do  ! ik

    ! ---- per-omega: eps00 - wing.body^-1.wing  (body = G,G' /= ig0) ----
    nb = ng - 1
    allocate(bidx(nb))
    a = 0
    do ig = 1, ng
      if (ig == ig0) cycle
      a = a + 1; bidx(a) = ig
    end do
    lwork = nb*max(nb,64)
    allocate(epsb(nb,nb), ipiv(nb), zwork(lwork))

    do iw = 1, nomega
      do b = 1, nb
        do a = 1, nb
          epsb(a,b) = - (fpi/gg(bidx(a))) * chi0b(bidx(a),bidx(b),iw)
        end do
        epsb(b,b) = epsb(b,b) + (1.0d0,0.0d0)
      end do
      call zgetrf(nb, nb, epsb, nb, ipiv, linfo)
      if (linfo /= 0) write(*,*) "[gw][velabs] zgetrf info=", linfo, " iw=", iw
      call zgetri(nb, epsb, nb, ipiv, zwork, lwork, linfo)
      if (linfo /= 0) write(*,*) "[gw][velabs] zgetri info=", linfo, " iw=", iw
      ! lfe = (4pi)^2 sum_ab A_{bidx(a)} ibody(a,b) At_{bidx(b)} / |G_{bidx(b)}|^2
      zlfe = (0.0d0,0.0d0)
      do b = 1, nb
        do a = 1, nb
          zlfe = zlfe + Aw(bidx(a),iw) * epsb(a,b) * Atw(bidx(b),iw) / gg(bidx(b))
        end do
      end do
      zlfe = (fpi*fpi) * zlfe
      eps_macro(iw) = (1.0d0,0.0d0) - fpi*Bw(iw) - zlfe
    end do

    deallocate(iGm, Bw, Aw, Atw, chi0b, epsb, ipiv, zwork, bidx)
  end subroutine calc_absorption_velocity

end module gw_absorption_sub
