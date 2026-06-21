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

end module gw_qp_output_sub
