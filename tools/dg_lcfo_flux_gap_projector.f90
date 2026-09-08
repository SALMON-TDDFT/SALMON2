program dg_lcfo_flux_gap_projector
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: au_to_ev = 27.211386245988_dp
  character(len=512) :: noflux_dir, flux_dir, arg
  character(len=32) :: mode
  integer :: nproj, nocc, nf1, nf2, nf3
  integer :: nfrag, nspin, nstate_frag, nstate_tot
  integer, allocatable :: nmat(:), nbasis(:,:), index_basis(:,:,:)
  real(dp), allocatable :: coef(:,:,:), hproj(:,:), eval(:), work(:)
  integer :: ifrag, info, lwork
  real(dp) :: gap_ev, herm_max

  interface
    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import :: dp
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      real(dp), intent(inout) :: a(lda,*), work(*)
      real(dp), intent(out) :: w(*)
      integer, intent(out) :: info
    end subroutine dsyev
  end interface

  call read_args(noflux_dir, flux_dir, nf1, nf2, nf3, nproj, nocc, mode)
  call read_wavefunction_header(trim(noflux_dir)//'/fragments/000001/wavefunctions.bin', &
       nfrag, nspin, nstate_frag, nstate_tot, nmat, nbasis, index_basis)
  if(nspin /= 1) stop 'spin-polarized diagnostic is not implemented'
  if(nfrag /= nf1*nf2*nf3) stop 'num_fragment product mismatch'
  if(nproj <= 0) nproj = min(nstate_tot, nmat(1))
  if(nproj > nstate_tot) stop 'nproj exceeds nstate_tot'
  if(nocc <= 0 .or. nocc >= nproj) stop 'nocc must be in 1..nproj-1'

  allocate(coef(nstate_frag,nproj,nfrag))
  do ifrag = 1, nfrag
    call read_wavefunction_coef(fragment_file(noflux_dir, ifrag, 'wavefunctions.bin'), &
         nfrag, nspin, nstate_frag, nstate_tot, nproj, coef(:,:,ifrag))
  end do

  allocate(hproj(nproj,nproj))
  hproj = 0.0_dp
  do ifrag = 1, nfrag
    call accumulate_fragment_flux_h(fragment_file(flux_dir, ifrag, hamiltonian_basename(mode)), &
         ifrag, nf1, nf2, nf3, nstate_frag, nproj, nbasis(:,1), coef, hproj, mode)
  end do

  herm_max = maxval(abs(hproj - transpose(hproj)))
  hproj = 0.5_dp * (hproj + transpose(hproj))

  allocate(eval(nproj))
  lwork = max(1, 8*nproj)
  allocate(work(lwork))
  call dsyev('V', 'U', nproj, hproj, nproj, eval, work, lwork, info)
  if(info /= 0) then
    write(*,'(a,i0)') 'dsyev failed: info=', info
    stop 1
  end if
  gap_ev = (eval(nocc+1) - eval(nocc)) * au_to_ev
  write(*,'(a)') '# DG LCFO no-Flux subspace / Flux-H projection diagnostic'
  write(*,'(a,1x,a)') 'noflux_dir =', trim(noflux_dir)
  write(*,'(a,1x,a)') 'flux_dir   =', trim(flux_dir)
  write(*,'(a,1x,a)') 'mode       =', trim(mode)
  write(*,'(a,4(1x,i0))') 'dims nfrag nstate_frag nstate_tot nproj =', &
       nfrag, nstate_frag, nstate_tot, nproj
  write(*,'(a,1x,i0,1x,a,1x,es14.6,1x,a,1x,es14.6)') 'nocc =', nocc, &
       'homo(eV)=', eval(nocc)*au_to_ev, 'lumo(eV)=', eval(nocc+1)*au_to_ev
  write(*,'(a,1x,es14.6)') 'gap(eV) =', gap_ev
  write(*,'(a,1x,es14.6)') 'max|H-H^T| before sym =', herm_max

contains

  subroutine read_args(noflux_dir, flux_dir, nf1, nf2, nf3, nproj, nocc, mode)
    character(len=*), intent(out) :: noflux_dir, flux_dir
    integer, intent(out) :: nf1, nf2, nf3, nproj, nocc
    character(len=*), intent(out) :: mode
    integer :: nargs
    character(len=512) :: arg
    nargs = command_argument_count()
    if(nargs < 7) then
      write(*,'(a)') 'usage: dg_lcfo_flux_gap_projector NOFLUX_DATA_DCDFT FLUX_DATA_DCDFT NF1 NF2 NF3 NPROJ NOCC [mode]'
      write(*,'(a)') 'modes: full block volume self volume_self cross_full component_full component_block'
      write(*,'(a)') '       volume_cross_full volume_cross_block'
      write(*,'(a)') '       weak_volume weak_volume_self weak_component_full weak_component_block'
      write(*,'(a)') '       weak_volume_cross_full weak_volume_cross_block'
      stop 2
    end if
    call get_command_argument(1, noflux_dir)
    call get_command_argument(2, flux_dir)
    call get_command_argument(3, arg); read(arg,*) nf1
    call get_command_argument(4, arg); read(arg,*) nf2
    call get_command_argument(5, arg); read(arg,*) nf3
    call get_command_argument(6, arg); read(arg,*) nproj
    call get_command_argument(7, arg); read(arg,*) nocc
    mode = 'full'
    if(nargs >= 8) call get_command_argument(8, mode)
    mode = adjustl(mode)
    select case(trim(mode))
    case('full','block','volume','self','volume_self','cross_full','component_full','component_block', &
         'volume_cross_full','volume_cross_block', &
         'weak_volume','weak_volume_self','weak_component_full','weak_component_block', &
         'weak_volume_cross_full','weak_volume_cross_block')
    case default
      stop 'unknown mode'
    end select
  end subroutine read_args

  character(len=64) function hamiltonian_basename(mode) result(name)
    character(len=*), intent(in) :: mode
    select case(trim(mode))
    case('weak_volume','weak_volume_self','weak_component_full','weak_component_block', &
         'weak_volume_cross_full','weak_volume_cross_block')
      name = 'hamiltonian_flux_weak_components.bin'
    case('volume','self','volume_self','cross_full','component_full','component_block', &
         'volume_cross_full','volume_cross_block')
      name = 'hamiltonian_flux_components.bin'
    case default
      name = 'hamiltonian_local.bin'
    end select
  end function hamiltonian_basename

  character(len=512) function fragment_file(data_dir, ifrag, basename) result(path)
    character(len=*), intent(in) :: data_dir, basename
    integer, intent(in) :: ifrag
    write(path,'(a,"/fragments/",i6.6,"/",a)') trim(data_dir), ifrag, trim(basename)
  end function fragment_file

  subroutine read_wavefunction_header(filename, nfrag, nspin, nstate_frag, nstate_tot, nmat, nbasis, index_basis)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: nfrag, nspin, nstate_frag, nstate_tot
    integer, allocatable, intent(out) :: nmat(:), nbasis(:,:), index_basis(:,:,:)
    integer :: unit
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) nfrag, nspin, nstate_frag, nstate_tot
    allocate(nmat(nspin), nbasis(nfrag,nspin), index_basis(nstate_frag,nfrag,nspin))
    read(unit) nmat(1:nspin)
    read(unit) nbasis(1:nfrag,1:nspin)
    read(unit) index_basis(1:nstate_frag,1:nfrag,1:nspin)
    close(unit)
  end subroutine read_wavefunction_header

  subroutine read_wavefunction_coef(filename, nfrag_expect, nspin_expect, nstate_frag_expect, nstate_tot_expect, nproj, coef_out)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nfrag_expect, nspin_expect, nstate_frag_expect, nstate_tot_expect, nproj
    real(dp), intent(out) :: coef_out(nstate_frag_expect,nproj)
    integer :: unit, nfrag, nspin, nstate_frag, nstate_tot
    integer, allocatable :: nmat_dummy(:), nbasis_dummy(:,:), index_dummy(:,:,:)
    real(dp), allocatable :: coef_all(:,:,:)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) nfrag, nspin, nstate_frag, nstate_tot
    if(nfrag /= nfrag_expect .or. nspin /= nspin_expect .or. &
       nstate_frag /= nstate_frag_expect .or. nstate_tot /= nstate_tot_expect) &
       stop 'wavefunctions.bin dimension mismatch'
    allocate(nmat_dummy(nspin), nbasis_dummy(nfrag,nspin), index_dummy(nstate_frag,nfrag,nspin))
    allocate(coef_all(nstate_frag,nstate_tot,nspin))
    read(unit) nmat_dummy(1:nspin)
    read(unit) nbasis_dummy(1:nfrag,1:nspin)
    read(unit) index_dummy(1:nstate_frag,1:nfrag,1:nspin)
    read(unit) coef_all(1:nstate_frag,1:nstate_tot,1:nspin)
    close(unit)
    coef_out(1:nstate_frag,1:nproj) = coef_all(1:nstate_frag,1:nproj,1)
    deallocate(nmat_dummy, nbasis_dummy, index_dummy, coef_all)
  end subroutine read_wavefunction_coef

  subroutine accumulate_fragment_flux_h(filename, ifrag, nf1, nf2, nf3, nstate_frag, nproj, nbasis, coef, hproj, mode)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: mode
    integer, intent(in) :: ifrag, nf1, nf2, nf3, nstate_frag, nproj, nbasis(:)
    real(dp), intent(in) :: coef(:,:,:)
    real(dp), intent(inout) :: hproj(nproj,nproj)
    real(dp), allocatable :: hloc(:,:,:), hself(:,:,:), hhalo(:,:,:), tmp(:,:)
    integer :: unit, n_halo, ih, jfrag, nb_i, nb_j, nb_t
    integer :: dvec(3)
    logical :: component_mode, include_cross, block_cross
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    allocate(hloc(nstate_frag,nstate_frag,1))
    read(unit) hloc
    component_mode = is_component_mode(mode)
    include_cross = (trim(mode) == 'full' .or. trim(mode) == 'block' .or. &
         trim(mode) == 'cross_full' .or. trim(mode) == 'component_full' .or. trim(mode) == 'component_block' .or. &
         trim(mode) == 'volume_cross_full' .or. trim(mode) == 'volume_cross_block' .or. &
         trim(mode) == 'weak_component_full' .or. trim(mode) == 'weak_component_block' .or. &
         trim(mode) == 'weak_volume_cross_full' .or. trim(mode) == 'weak_volume_cross_block')
    block_cross = (trim(mode) == 'block' .or. trim(mode) == 'component_block' .or. &
         trim(mode) == 'volume_cross_block' .or. trim(mode) == 'weak_component_block' .or. &
         trim(mode) == 'weak_volume_cross_block')
    if(component_mode) then
      allocate(hself(nstate_frag,nstate_frag,1))
      read(unit) hself
      select case(trim(mode))
      case('volume','weak_volume','volume_cross_full','volume_cross_block', &
           'weak_volume_cross_full','weak_volume_cross_block')
        ! hloc already contains the strong-form volume component.
      case('self')
        hloc = hself
      case('volume_self','component_full','component_block','weak_volume_self', &
           'weak_component_full','weak_component_block')
        hloc = hloc + hself
      case('cross_full')
        hloc = 0.0_dp
      end select
      deallocate(hself)
    end if
    read(unit) n_halo
    allocate(tmp(nstate_frag,nproj))
    nb_i = nbasis(ifrag)
    tmp = 0.0_dp
    if(trim(mode) == 'block' .or. trim(mode) == 'component_block') then
      tmp(1:nb_i,1:nproj) = matmul(0.5_dp * (hloc(1:nb_i,1:nb_i,1) + transpose(hloc(1:nb_i,1:nb_i,1))), &
           coef(1:nb_i,1:nproj,ifrag))
    else
      tmp(1:nb_i,1:nproj) = matmul(hloc(1:nb_i,1:nb_i,1), coef(1:nb_i,1:nproj,ifrag))
    end if
    hproj = hproj + matmul(transpose(coef(1:nb_i,1:nproj,ifrag)), tmp(1:nb_i,1:nproj))
    deallocate(hloc)

    do ih = 1, n_halo
      allocate(hhalo(nstate_frag,nstate_frag,1))
      read(unit) hhalo
      if(.not. include_cross) then
        deallocate(hhalo)
        cycle
      end if
      dvec = halo_dvec(ih, nf1, nf2, nf3, n_halo)
      jfrag = source_neighbor(ifrag, dvec, nf1, nf2, nf3)
      nb_j = nbasis(jfrag)
      if(block_cross) then
        nb_t = min(nb_i, size(hhalo,1), size(hhalo,2))
        tmp = 0.0_dp
        tmp(1:nb_t,1:nproj) = matmul(0.25_dp * (hhalo(1:nb_t,1:nb_t,1) + transpose(hhalo(1:nb_t,1:nb_t,1))), &
             coef(1:nb_t,1:nproj,ifrag))
        hproj = hproj + matmul(transpose(coef(1:nb_t,1:nproj,ifrag)), tmp(1:nb_t,1:nproj))
        nb_t = min(nb_j, size(hhalo,1), size(hhalo,2))
        tmp = 0.0_dp
        tmp(1:nb_t,1:nproj) = matmul(0.25_dp * (hhalo(1:nb_t,1:nb_t,1) + transpose(hhalo(1:nb_t,1:nb_t,1))), &
             coef(1:nb_t,1:nproj,jfrag))
        hproj = hproj + matmul(transpose(coef(1:nb_t,1:nproj,jfrag)), tmp(1:nb_t,1:nproj))
      else
        tmp = 0.0_dp
        tmp(1:nb_j,1:nproj) = matmul(hhalo(1:nb_j,1:nb_i,1), coef(1:nb_i,1:nproj,ifrag))
        hproj = hproj + 0.5_dp * matmul(transpose(coef(1:nb_j,1:nproj,jfrag)), tmp(1:nb_j,1:nproj))
        tmp = 0.0_dp
        tmp(1:nb_i,1:nproj) = matmul(transpose(hhalo(1:nb_j,1:nb_i,1)), coef(1:nb_j,1:nproj,jfrag))
        hproj = hproj + 0.5_dp * matmul(transpose(coef(1:nb_i,1:nproj,ifrag)), tmp(1:nb_i,1:nproj))
      end if
      deallocate(hhalo)
    end do
    close(unit)
    deallocate(tmp)
  end subroutine accumulate_fragment_flux_h

  logical function is_component_mode(mode) result(ok)
    character(len=*), intent(in) :: mode
    select case(trim(mode))
    case('volume','self','volume_self','cross_full','component_full','component_block', &
         'volume_cross_full','volume_cross_block', &
         'weak_volume','weak_volume_self','weak_component_full','weak_component_block', &
         'weak_volume_cross_full','weak_volume_cross_block')
      ok = .true.
    case default
      ok = .false.
    end select
  end function is_component_mode

  integer function source_neighbor(ifrag, dvec, nf1, nf2, nf3) result(jfrag)
    integer, intent(in) :: ifrag, dvec(3), nf1, nf2, nf3
    integer :: ix, iy, iz, rem
    ix = (ifrag - 1) / (nf2*nf3) + 1
    rem = modulo(ifrag - 1, nf2*nf3)
    iy = rem / nf3 + 1
    iz = modulo(rem, nf3) + 1
    ix = modulo(ix - 1 - dvec(1), nf1) + 1
    iy = modulo(iy - 1 - dvec(2), nf2) + 1
    iz = modulo(iz - 1 - dvec(3), nf3) + 1
    jfrag = ((ix - 1) * nf2 + (iy - 1)) * nf3 + iz
  end function source_neighbor

  function halo_dvec(ih_target, nf1, nf2, nf3, n_halo) result(dvec)
    integer, intent(in) :: ih_target, nf1, nf2, nf3, n_halo
    integer :: dvec(3), ih, lx, ly, lz, nh(3)
    logical :: face_only
    nh = 0
    if(nf1 > 1) nh(1) = 1
    if(nf2 > 1) nh(2) = 1
    if(nf3 > 1) nh(3) = 1
    face_only = (n_halo == count(nh(1:3) > 0) * 2)
    ih = 0
    dvec = 0
    do lx = -nh(1), nh(1)
    do ly = -nh(2), nh(2)
    do lz = -nh(3), nh(3)
      if(lx == 0 .and. ly == 0 .and. lz == 0) cycle
      if(face_only .and. count([lx,ly,lz] /= 0) /= 1) cycle
      ih = ih + 1
      if(ih == ih_target) then
        dvec = [lx, ly, lz]
        return
      end if
    end do
    end do
    end do
    if(ih_target > n_halo) stop 'halo index exceeds file n_halo'
  end function halo_dvec

end program dg_lcfo_flux_gap_projector
