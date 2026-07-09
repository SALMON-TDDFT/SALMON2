!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module gw_qp_output_sub
  implicit none

contains

! Write quasiparticle energies to a text file.
! Root process only. Output format mirrors write_eigen in src/io/write.f90.
!
! Arguments:
!   system    - DFT system (provides nspin, nk, no, rocc)
!   energy_ks - KS single-particle energies (no, nk, nspin), a.u.
!   eqp       - quasiparticle energies       (no, nk, nspin), a.u.
!   zfac      - quasiparticle renormalisation weights (no, nk, nspin)
!   sigx      - exchange self-energy          (no, nk, nspin), a.u.
!   sigc      - correlation self-energy       (no, nk, nspin), a.u.
!   vxc       - XC potential matrix element  (no, nk, nspin), a.u.

subroutine write_qp_energies(system, energy_ks, eqp, zfac, sigx, sigc, vxc)
  use structures,      only: s_dft_system
  use parallelization, only: nproc_id_global
  use communication,   only: comm_is_root
  use inputoutput,     only: uenergy_from_au, iperiodic, base_directory, sysname, unit_energy
  use filesystem,      only: open_filehandle
  implicit none
  type(s_dft_system),intent(in) :: system
  real(8),           intent(in) :: energy_ks(system%no, system%nk, system%nspin)
  real(8),           intent(in) :: eqp      (system%no, system%nk, system%nspin)
  real(8),           intent(in) :: zfac     (system%no, system%nk, system%nspin)
  real(8),           intent(in) :: sigx     (system%no, system%nk, system%nspin)
  real(8),           intent(in) :: sigc     (system%no, system%nk, system%nspin)
  real(8),           intent(in) :: vxc      (system%no, system%nk, system%nspin)

  integer          :: uid, io, iik, is
  character(256)   :: fname

  if (.not. comm_is_root(nproc_id_global)) return

  fname = trim(base_directory)//trim(sysname)//"_qp_energies.data"
  uid   = open_filehandle(trim(fname))

  write(uid,'("# quasiparticle energies (QP=KS passthrough scaffold)")')
  select case(unit_energy)
  case('au','a.u.')
    write(uid,'("# 1:io 2:Eks[a.u.] 3:Eqp[a.u.] 4:Z 5:<Sigx>[a.u.] 6:<Sigc>[a.u.] 7:<Vxc>[a.u.] 8:occ")')
  case default
    write(uid,'("# 1:io 2:Eks[eV] 3:Eqp[eV] 4:Z 5:<Sigx>[eV] 6:<Sigc>[eV] 7:<Vxc>[eV] 8:occ")')
  end select

  do is  = 1, system%nspin
  do iik = 1, system%nk
    if (iperiodic == 3) then
      write(uid,'("k=",1x,i5,",  spin=",1x,i5)') iik, is
    end if
    do io = 1, system%no
      write(uid,'(1x,i5,7(e26.16e3))') io,                                &
        energy_ks(io,iik,is)*uenergy_from_au,                              &
        eqp      (io,iik,is)*uenergy_from_au,                              &
        zfac     (io,iik,is),                                              &
        sigx     (io,iik,is)*uenergy_from_au,                              &
        sigc     (io,iik,is)*uenergy_from_au,                              &
        vxc      (io,iik,is)*uenergy_from_au,                              &
        system%rocc(io,iik,is)
    end do
  end do
  end do

  close(uid)

end subroutine write_qp_energies

! Write the KS and QP density of states (Gaussian-broadened sum over all
! bands and k-points) to <sysname>_qp_dos.data for an overlay plot.
! Root process only.  Energy axis always in eV (header states so).
!
!   DOS_X(E) = sum_{n,k,s} w_k * (2/nspin) * g_sigma(E - E_X(n,k,s)),  X = KS, QP
!   g_sigma  = normalised Gaussian, sigma = out_dos_width (a.u.->eV; default 0.15 eV)

subroutine write_qp_dos(system, energy_ks, eqp)
  use structures,      only: s_dft_system
  use parallelization, only: nproc_id_global
  use communication,   only: comm_is_root
  use inputoutput,     only: base_directory, sysname, out_dos_width
  use filesystem,      only: open_filehandle
  implicit none
  type(s_dft_system),intent(in) :: system
  real(8),           intent(in) :: energy_ks(system%no, system%nk, system%nspin)
  real(8),           intent(in) :: eqp      (system%no, system%nk, system%nspin)

  integer,parameter   :: nE = 1200
  real(8),parameter   :: ev = 27.211386d0
  real(8),allocatable :: egrid(:), dos_ks(:), dos_qp(:)
  real(8) :: emin, emax, sigma, sp, gnorm, arg
  integer :: uid, io, iik, is, ie

  if (.not. comm_is_root(nproc_id_global)) return

  emin =  1.0d30; emax = -1.0d30
  do is = 1, system%nspin
  do iik = 1, system%nk
  do io = 1, system%no
    emin = min(emin, energy_ks(io,iik,is)*ev, eqp(io,iik,is)*ev)
    emax = max(emax, energy_ks(io,iik,is)*ev, eqp(io,iik,is)*ev)
  end do
  end do
  end do
  emin = emin - 3.0d0; emax = emax + 3.0d0

  sigma = out_dos_width * ev                 ! a.u. -> eV
  if (sigma < 1.0d-3) sigma = 0.15d0         ! default broadening
  sp    = 2.0d0 / dble(system%nspin)         ! spin degeneracy (nspin=1 -> 2)
  gnorm = 1.0d0 / (sigma * sqrt(2.0d0*acos(-1.0d0)))

  allocate(egrid(nE), dos_ks(nE), dos_qp(nE))
  dos_ks(:) = 0.0d0; dos_qp(:) = 0.0d0
  do ie = 1, nE
    egrid(ie) = emin + (emax-emin)*dble(ie-1)/dble(nE-1)
  end do
  do is = 1, system%nspin
  do iik = 1, system%nk
  do io = 1, system%no
    do ie = 1, nE
      arg = (egrid(ie) - energy_ks(io,iik,is)*ev) / sigma
      dos_ks(ie) = dos_ks(ie) + system%wtk(iik)*sp*gnorm*exp(-0.5d0*arg*arg)
      arg = (egrid(ie) - eqp(io,iik,is)*ev) / sigma
      dos_qp(ie) = dos_qp(ie) + system%wtk(iik)*sp*gnorm*exp(-0.5d0*arg*arg)
    end do
  end do
  end do
  end do

  uid = open_filehandle(trim(base_directory)//trim(sysname)//"_qp_dos.data")
  write(uid,'("# QP vs KS density of states (Gaussian-broadened)")')
  write(uid,'("# sigma[eV]=",f8.4,"  spin_factor=",f4.1)') sigma, sp
  write(uid,'("# 1:E[eV] 2:DOS_KS[1/eV] 3:DOS_QP[1/eV]")')
  do ie = 1, nE
    write(uid,'(3(e22.12e3))') egrid(ie), dos_ks(ie), dos_qp(ie)
  end do
  close(uid)

  deallocate(egrid, dos_ks, dos_qp)

end subroutine write_qp_dos

! Read back the per-state QP energies written by write_qp_energies, for the
! true G0W0+RPA optics path (theory='gw_response', yn_gw_qp_inject='y'): the
! full per-(band,k,spin) QP spectrum is injected into chi0 instead of a rigid
! scissors.  Root parses the text file; result is broadcast to all ranks.
! Returns eqp_au in a.u. and ok=.false. if the file is absent/unreadable.

subroutine read_qp_energies(system, eqp_au, ok)
  use structures,      only: s_dft_system
  use parallelization, only: nproc_id_global, nproc_group_global
  use communication,   only: comm_is_root, comm_bcast
  use inputoutput,     only: uenergy_from_au, base_directory, sysname
  implicit none
  type(s_dft_system),intent(in)  :: system
  real(8),           intent(out) :: eqp_au(system%no, system%nk, system%nspin)
  logical,           intent(out) :: ok

  integer        :: uid, ios, ik, is, io, ic, icc, iok
  real(8)        :: vals(7)
  character(512) :: line
  character(256) :: fname
  logical        :: ex

  ok = .false.
  iok = 0
  eqp_au(:,:,:) = 0.0d0

  if (comm_is_root(nproc_id_global)) then
    fname = trim(base_directory)//trim(sysname)//"_qp_energies.data"
    inquire(file=trim(fname), exist=ex)
    if (ex) then
      open(newunit=uid, file=trim(fname), status='old', action='read')
      ik = 1; is = 1
      do
        read(uid,'(A)',iostat=ios) line
        if (ios /= 0) exit
        line = adjustl(line)
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#') cycle
        ic = index(line, 'k=')
        if (ic > 0) then                          ! "k=  N,  spin=  M" header
          read(line(ic+2:),*,iostat=ios) ik
          icc = index(line, 'spin=')
          if (icc > 0) read(line(icc+5:),*,iostat=ios) is
          cycle
        end if
        read(line,*,iostat=ios) io, vals(1:7)     ! io Eks Eqp Z Sigx Sigc Vxc occ
        if (ios /= 0) cycle
        if (io >= 1 .and. io <= system%no .and. ik >= 1 .and. ik <= system%nk &
            .and. is >= 1 .and. is <= system%nspin) then
          eqp_au(io,ik,is) = vals(2) / uenergy_from_au   ! Eqp column -> a.u.
        end if
      end do
      close(uid)
      iok = 1
    end if
  end if

  call comm_bcast(iok,    nproc_group_global)
  call comm_bcast(eqp_au, nproc_group_global)
  ok = (iok == 1)

end subroutine read_qp_energies

end module gw_qp_output_sub
