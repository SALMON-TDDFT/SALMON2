program dg_wannier_flux_gap_projector
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: bpw_magic = -22022215
  integer, parameter :: bpwt_magic = -22022218
  real(dp), parameter :: au_to_ev = 27.211386245988_dp

  type face_trace
    integer :: axis = 0, side = 0, npts = 0
    real(dp) :: area_weight = 0.0_dp, alpha = 0.0_dp
    real(dp), allocatable :: u(:,:), dn(:,:)
  end type face_trace

  type frag_data
    integer :: nkeep = 0
    integer :: nxyz_domain(3) = 0
    integer :: nxyz_buffer(3) = 0
    integer :: nxyz_box(3) = 0
    real(dp), allocatable :: h(:,:)
    type(face_trace) :: face(6)
  end type frag_data

  character(len=512) :: data_dir
  character(len=32) :: mode
  integer :: nf1, nf2, nf3, nfrag, nocc, ntot
  integer :: ifrag, ioff, info, lwork
  integer, allocatable :: offset(:)
  real(dp), allocatable :: hmat(:,:), eval(:), work(:)
  type(frag_data), allocatable :: frag(:)
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

  call read_args(data_dir, nf1, nf2, nf3, nocc, mode)
  nfrag = nf1 * nf2 * nf3
  allocate(frag(nfrag), offset(nfrag+1))

  offset(1) = 1
  do ifrag=1,nfrag
    call read_bpw(fragment_file(data_dir, ifrag, 'buffer_periodic_wannier_basis.bin'), frag(ifrag))
    call read_trace(fragment_file(data_dir, ifrag, 'buffer_periodic_wannier_trace.bin'), frag(ifrag))
    offset(ifrag+1) = offset(ifrag) + frag(ifrag)%nkeep
  end do
  ntot = offset(nfrag+1) - 1
  if(nocc <= 0 .or. nocc >= ntot) stop 'nocc must be in 1..ntot-1'

  allocate(hmat(ntot,ntot))
  hmat = 0.0_dp
  if(include_local(mode)) then
    do ifrag=1,nfrag
      ioff = offset(ifrag)
      hmat(ioff:ioff+frag(ifrag)%nkeep-1,ioff:ioff+frag(ifrag)%nkeep-1) = &
        hmat(ioff:ioff+frag(ifrag)%nkeep-1,ioff:ioff+frag(ifrag)%nkeep-1) + frag(ifrag)%h
    end do
  end if
  if(include_self(mode)) call accumulate_self_surface(frag, offset, nfrag, hmat)
  if(include_cross(mode)) call accumulate_cross_surface(frag, offset, nf1, nf2, nf3, hmat)

  herm_max = maxval(abs(hmat - transpose(hmat)))
  hmat = 0.5_dp * (hmat + transpose(hmat))

  allocate(eval(ntot))
  lwork = max(1, 8*ntot)
  allocate(work(lwork))
  call dsyev('V', 'U', ntot, hmat, ntot, eval, work, lwork, info)
  if(info /= 0) then
    write(*,'(a,i0)') 'dsyev failed: info=', info
    stop 1
  end if
  gap_ev = (eval(nocc+1) - eval(nocc)) * au_to_ev
  write(*,'(a)') '# DG Wannier trace / Flux-H diagnostic'
  write(*,'(a,1x,a)') 'data_dir =', trim(data_dir)
  write(*,'(a,1x,a)') 'mode     =', trim(mode)
  write(*,'(a,3(1x,i0),1x,a,1x,i0)') 'fragments =', nf1, nf2, nf3, 'ntot =', ntot
  write(*,'(a,1x,i0,1x,a,1x,es14.6,1x,a,1x,es14.6)') 'nocc =', nocc, &
    'homo(eV)=', eval(nocc)*au_to_ev, 'lumo(eV)=', eval(nocc+1)*au_to_ev
  write(*,'(a,1x,es14.6)') 'gap(eV) =', gap_ev
  write(*,'(a,1x,es14.6)') 'max|H-H^T| before sym =', herm_max

contains

  subroutine read_args(data_dir, nf1, nf2, nf3, nocc, mode)
    character(len=*), intent(out) :: data_dir, mode
    integer, intent(out) :: nf1, nf2, nf3, nocc
    integer :: nargs
    character(len=512) :: arg
    nargs = command_argument_count()
    if(nargs < 5) then
      write(*,'(a)') 'usage: dg_wannier_flux_gap_projector DATA_DCDFT NF1 NF2 NF3 NOCC [mode]'
      write(*,'(a)') 'modes: local self cross full surface'
      stop 2
    end if
    call get_command_argument(1, data_dir)
    call get_command_argument(2, arg); read(arg,*) nf1
    call get_command_argument(3, arg); read(arg,*) nf2
    call get_command_argument(4, arg); read(arg,*) nf3
    call get_command_argument(5, arg); read(arg,*) nocc
    mode = 'full'
    if(nargs >= 6) call get_command_argument(6, mode)
    mode = adjustl(mode)
    select case(trim(mode))
    case('local','self','cross','full','surface')
    case default
      stop 'unknown mode'
    end select
  end subroutine read_args

  logical function include_local(mode)
    character(len=*), intent(in) :: mode
    include_local = (trim(mode) == 'local' .or. trim(mode) == 'full')
  end function include_local

  logical function include_self(mode)
    character(len=*), intent(in) :: mode
    include_self = (trim(mode) == 'self' .or. trim(mode) == 'full' .or. trim(mode) == 'surface')
  end function include_self

  logical function include_cross(mode)
    character(len=*), intent(in) :: mode
    include_cross = (trim(mode) == 'cross' .or. trim(mode) == 'full' .or. trim(mode) == 'surface')
  end function include_cross

  character(len=512) function fragment_file(data_dir, ifrag, basename) result(path)
    character(len=*), intent(in) :: data_dir, basename
    integer, intent(in) :: ifrag
    write(path,'(a,"/fragments/",i6.6,"/",a)') trim(data_dir), ifrag, trim(basename)
  end function fragment_file

  subroutine read_bpw(filename, frag)
    character(len=*), intent(in) :: filename
    type(frag_data), intent(inout) :: frag
    integer :: unit, magic, version, ifrag_file, nspin, nbasis, nkeep_legacy
    integer, allocatable :: ivec(:)
    real(dp), allocatable :: rvec(:), wcoef(:,:), r_wann(:,:,:), wcenter(:,:)
    real(dp), allocatable :: v_wann(:,:,:), aa_wann(:,:,:)
    integer :: aa_src
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic, version
    if(magic /= bpw_magic) stop 'invalid BPW magic'
    read(unit) ifrag_file, frag%nxyz_domain(1:3), frag%nxyz_buffer(1:3), frag%nxyz_box(1:3)
    read(unit) nspin, nbasis, frag%nkeep, nkeep_legacy
    if(nspin /= 1) stop 'spin-polarized BPW is not supported'
    allocate(ivec(frag%nkeep), rvec(frag%nkeep))
    allocate(wcoef(nbasis,frag%nkeep))
    allocate(r_wann(3,frag%nkeep,frag%nkeep), wcenter(3,frag%nkeep))
    allocate(frag%h(frag%nkeep,frag%nkeep))
    allocate(v_wann(3,frag%nkeep,frag%nkeep), aa_wann(3,frag%nkeep,frag%nkeep))
    read(unit) ivec(1:frag%nkeep), ivec(1:frag%nkeep)
    read(unit) rvec(1:frag%nkeep), ivec(1:frag%nkeep)
    read(unit) rvec(1:frag%nkeep), rvec(1:frag%nkeep)
    read(unit) wcoef(1:nbasis,1:frag%nkeep)
    read(unit) r_wann(1:3,1:frag%nkeep,1:frag%nkeep)
    read(unit) wcenter(1:3,1:frag%nkeep)
    read(unit) frag%h(1:frag%nkeep,1:frag%nkeep)
    read(unit) v_wann(1:3,1:frag%nkeep,1:frag%nkeep)
    read(unit) aa_wann(1:3,1:frag%nkeep,1:frag%nkeep)
    read(unit) aa_src
    close(unit)
    deallocate(ivec, rvec, wcoef, r_wann, wcenter, v_wann, aa_wann)
  end subroutine read_bpw

  subroutine read_trace(filename, frag)
    character(len=*), intent(in) :: filename
    type(frag_data), intent(inout) :: frag
    integer :: unit, magic, version, ifrag_file, nkeep_file, face
    integer :: nxyz_domain_file(3), nxyz_buffer_file(3), nxyz_box_file(3)
    real(dp) :: hgs(3), hvol
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic, version
    if(magic /= bpwt_magic) stop 'invalid BPW trace magic'
    read(unit) ifrag_file, nxyz_domain_file(1:3), nxyz_buffer_file(1:3), nxyz_box_file(1:3)
    if(any(nxyz_domain_file /= frag%nxyz_domain) .or. any(nxyz_buffer_file /= frag%nxyz_buffer) .or. &
       any(nxyz_box_file /= frag%nxyz_box)) stop 'BPW trace metadata mismatch'
    read(unit) hgs(1:3), hvol, nkeep_file
    if(nkeep_file /= frag%nkeep) stop 'BPW trace nkeep mismatch'
    do face=1,6
      read(unit) frag%face(face)%axis, frag%face(face)%side, frag%face(face)%npts, &
        frag%face(face)%area_weight, frag%face(face)%alpha
      allocate(frag%face(face)%u(frag%face(face)%npts,frag%nkeep))
      allocate(frag%face(face)%dn(frag%face(face)%npts,frag%nkeep))
      read(unit) frag%face(face)%u(1:frag%face(face)%npts,1:frag%nkeep)
      read(unit) frag%face(face)%dn(1:frag%face(face)%npts,1:frag%nkeep)
    end do
    close(unit)
  end subroutine read_trace

  subroutine accumulate_self_surface(frag, offset, nfrag, hmat)
    type(frag_data), intent(in) :: frag(:)
    integer, intent(in) :: offset(:), nfrag
    real(dp), intent(inout) :: hmat(:,:)
    integer :: ifrag, face, i, j, ip, ioff, ni
    real(dp) :: val, v_l, dnv_l, u_l, dnu_l
    do ifrag=1,nfrag
      ioff = offset(ifrag)
      ni = frag(ifrag)%nkeep
      do face=1,6
        do i=1,ni
        do j=1,ni
          val = 0.0_dp
          do ip=1,frag(ifrag)%face(face)%npts
            v_l = frag(ifrag)%face(face)%u(ip,i)
            dnv_l = frag(ifrag)%face(face)%dn(ip,i)
            u_l = frag(ifrag)%face(face)%u(ip,j)
            dnu_l = frag(ifrag)%face(face)%dn(ip,j)
            val = val + (-0.25_dp*v_l*dnu_l - 0.25_dp*dnv_l*u_l + &
              frag(ifrag)%face(face)%alpha*v_l*u_l) * frag(ifrag)%face(face)%area_weight
          end do
          hmat(ioff+i-1,ioff+j-1) = hmat(ioff+i-1,ioff+j-1) + val
        end do
        end do
      end do
    end do
  end subroutine accumulate_self_surface

  subroutine accumulate_cross_surface(frag, offset, nf1, nf2, nf3, hmat)
    type(frag_data), intent(in) :: frag(:)
    integer, intent(in) :: offset(:), nf1, nf2, nf3
    real(dp), intent(inout) :: hmat(:,:)
    integer :: ifrag, jfrag, face, opp, i, j, ip, ioff, joff, ni, nj
    integer :: axis, side
    real(dp) :: val, v_l, dnv_l, u_r, dnu_r
    do ifrag=1,nf1*nf2*nf3
      ioff = offset(ifrag)
      ni = frag(ifrag)%nkeep
      do face=1,6
        axis = frag(ifrag)%face(face)%axis
        side = frag(ifrag)%face(face)%side
        jfrag = neighbor_fragment(ifrag, axis, side, nf1, nf2, nf3)
        opp = face_index(axis, -side)
        joff = offset(jfrag)
        nj = frag(jfrag)%nkeep
        if(frag(ifrag)%face(face)%npts /= frag(jfrag)%face(opp)%npts) &
          stop 'cross face point mismatch'
        do i=1,ni
        do j=1,nj
          val = 0.0_dp
          do ip=1,frag(ifrag)%face(face)%npts
            v_l = frag(ifrag)%face(face)%u(ip,i)
            dnv_l = frag(ifrag)%face(face)%dn(ip,i)
            u_r = frag(jfrag)%face(opp)%u(ip,j)
            dnu_r = -frag(jfrag)%face(opp)%dn(ip,j)
            val = val + (-0.25_dp*v_l*dnu_r + 0.25_dp*dnv_l*u_r - &
              frag(ifrag)%face(face)%alpha*v_l*u_r) * frag(ifrag)%face(face)%area_weight
          end do
          hmat(ioff+i-1,joff+j-1) = hmat(ioff+i-1,joff+j-1) + val
        end do
        end do
      end do
    end do
  end subroutine accumulate_cross_surface

  integer function face_index(axis, side) result(face)
    integer, intent(in) :: axis, side
    face = 2*axis
    if(side < 0) face = 2*axis - 1
  end function face_index

  integer function neighbor_fragment(ifrag, axis, side, nf1, nf2, nf3) result(jfrag)
    integer, intent(in) :: ifrag, axis, side, nf1, nf2, nf3
    integer :: ix, iy, iz, rem
    ix = (ifrag - 1) / max(1, nf2*nf3) + 1
    rem = modulo(ifrag - 1, max(1, nf2*nf3))
    iy = rem / max(1, nf3) + 1
    iz = modulo(rem, max(1, nf3)) + 1
    select case(axis)
    case(1)
      ix = modulo(ix - 1 + side, nf1) + 1
    case(2)
      iy = modulo(iy - 1 + side, nf2) + 1
    case(3)
      iz = modulo(iz - 1 + side, nf3) + 1
    end select
    jfrag = ((ix - 1) * nf2 + (iy - 1)) * nf3 + iz
  end function neighbor_fragment
end program dg_wannier_flux_gap_projector
