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
module salmon_pp
  implicit none

  contains

  subroutine init_pp(pp,nrmax,lmax,flag_nlcc)
    use structures, only: s_pp_info
    use salmon_global,only : nelem,lloc_ps
    use salmon_global,only : file_pseudo
    use salmon_global,only : n_Yabana_Bertsch_psformat,n_ABINIT_psformat, &
                             n_ABINITFHI_psformat,n_FHI_psformat, &
                             ps_format,izatom,nelem
    use salmon_parallel, only: nproc_group_global, nproc_id_global
    use salmon_communication, only: comm_bcast, comm_is_root
    implicit none
    type(s_pp_info) :: pp
    integer, parameter :: nrmax0=50000, lmax0=4
    integer,intent(in) :: nrmax,lmax
    logical,intent(in) :: flag_nlcc
    integer :: ik
    character(256) :: ps_file
    integer :: ips_type,nlen_psf

    allocate(pp%atom_symbol(nelem))
    allocate(pp%rmass(nelem))
    allocate(pp%mr(nelem))

    if (comm_is_root(nproc_id_global)) then
  
      do ik=1,nelem
         
        ps_file = trim(file_pseudo(ik))
        nlen_psf = len_trim(ps_file)

        if(ps_file(max(1,nlen_psf+1-8):nlen_psf) == '_rps.dat')then
          ips_type = n_Yabana_Bertsch_psformat
          ps_format(ik) = 'KY'
          call read_mr_yb(pp,ik,ps_file)
        else if(ps_file(max(1,nlen_psf+1-6):nlen_psf) == '.pspnc')then
          ips_type = n_ABINIT_psformat
          ps_format(ik) = 'ABINIT'
          call read_mr_abinit(pp,ik,ps_file)
        else if(ps_file(max(1,nlen_psf+1-5):nlen_psf) == '.psp8')then
          ips_type = n_ABINIT_psformat
          ps_format(ik) = 'ABINITPSP8'
          call read_mr_abinit(pp,ik,ps_file)
        else if(ps_file(max(1,nlen_psf+1-4):nlen_psf) == '.fhi')then
          ips_type = n_ABINITFHI_psformat
          ps_format(ik) = 'ABINITFHI'
          call read_mr_abinitfhi(pp,ik,ps_file)
        else if(ps_file(max(1,nlen_psf+1-4):nlen_psf) == '.cpi')then
          ips_type = n_FHI_psformat
          ps_format(ik) = 'FHI'
          call read_mr_fhi(pp,ik,ps_file)
        else if(ps_file(max(1,nlen_psf+1-4):nlen_psf) == '.vps')then
          ips_type = 0
          ps_format(ik) = 'ADPACK'
        else
          stop 'Unprepared ps_format is required input_pseudopotential_YS'
        end if

!! --- Making prefix ---
!      select case (ipsfileform(ik))
!      case(n_Yabana_Bertsch_psformat)   ; ps_postfix = '_rps.dat'
!      case(n_ABINIT_psformat)           ; ps_postfix = '.pspnc'
!      case(n_ABINITFHI_psformat)        ; ps_postfix = '.fhi'
!      case(n_FHI_psformat)              ; ps_postfix = '.cpi'
!!    case('ATOM')      ; ps_postfix = '.psf' !Not implemented yet
!      case default ; stop 'Unprepared ps_format is required input_pseudopotential_YS'
!      end select

! --- input pseudopotential and wave function ---
        select case (izatom(ik))
        case (1) ; pp%atom_symbol(ik) = 'H ' ; pp%rmass(ik)=1.d0
        case (2) ; pp%atom_symbol(ik) = 'He' ; pp%rmass(ik)=4.d0
        case (3) ; pp%atom_symbol(ik) = 'Li' ; pp%rmass(ik)=7.d0
        case (4) ; pp%atom_symbol(ik) = 'Be' ; pp%rmass(ik)=9.d0
        case (5) ; pp%atom_symbol(ik) = 'B ' ; pp%rmass(ik)=11.d0
        case (6) ; pp%atom_symbol(ik) = 'C ' ; pp%rmass(ik)=12.d0
        case (7) ; pp%atom_symbol(ik) = 'N ' ; pp%rmass(ik)=14.d0
        case (8) ; pp%atom_symbol(ik) = 'O ' ; pp%rmass(ik)=16.d0
        case (9) ; pp%atom_symbol(ik) = 'F ' ; pp%rmass(ik)=19.d0
        case(10) ; pp%atom_symbol(ik) = 'Ne' ; pp%rmass(ik)=20.d0
        case(11) ; pp%atom_symbol(ik) = 'Na' ; pp%rmass(ik)=23.d0
        case(12) ; pp%atom_symbol(ik) = 'Mg' ; pp%rmass(ik)=24.d0
        case(13) ; pp%atom_symbol(ik) = 'Al' ; pp%rmass(ik)=27.d0
        case(14) ; pp%atom_symbol(ik) = 'Si' ; pp%rmass(ik)=28.d0
        case(15) ; pp%atom_symbol(ik) = 'P ' ; pp%rmass(ik)=31.d0
        case(16) ; pp%atom_symbol(ik) = 'S ' ; pp%rmass(ik)=32.d0
        case(17) ; pp%atom_symbol(ik) = 'Cl' ; pp%rmass(ik)=35.d0
        case(18) ; pp%atom_symbol(ik) = 'Ar' ; pp%rmass(ik)=40.d0
        case(19) ; pp%atom_symbol(ik) = 'K ' ; pp%rmass(ik)=39.d0
        case(20) ; pp%atom_symbol(ik) = 'Ca' ; pp%rmass(ik)=40.d0
        case(21) ; pp%atom_symbol(ik) = 'Sc' ; pp%rmass(ik)=45.d0
        case(22) ; pp%atom_symbol(ik) = 'Ti' ; pp%rmass(ik)=48.d0
        case(23) ; pp%atom_symbol(ik) = 'V ' ; pp%rmass(ik)=51.d0
        case(24) ; pp%atom_symbol(ik) = 'Cr' ; pp%rmass(ik)=52.d0
        case(25) ; pp%atom_symbol(ik) = 'Mn' ; pp%rmass(ik)=55.d0
        case(26) ; pp%atom_symbol(ik) = 'Fe' ; pp%rmass(ik)=56.d0
        case(27) ; pp%atom_symbol(ik) = 'Co' ; pp%rmass(ik)=59.d0
        case(28) ; pp%atom_symbol(ik) = 'Ni' ; pp%rmass(ik)=59.d0
        case(29) ; pp%atom_symbol(ik) = 'Cu' ; pp%rmass(ik)=63.d0
        case(30) ; pp%atom_symbol(ik) = 'Zn' ; pp%rmass(ik)=65.d0
        case(31) ; pp%atom_symbol(ik) = 'Ga' ; pp%rmass(ik)=69.d0
        case(32) ; pp%atom_symbol(ik) = 'Ge' ; pp%rmass(ik)=73.d0
        case(33) ; pp%atom_symbol(ik) = 'As' ; pp%rmass(ik)=75.d0
        case(34) ; pp%atom_symbol(ik) = 'Se' ; pp%rmass(ik)=79.d0
        case(35) ; pp%atom_symbol(ik) = 'Br' ; pp%rmass(ik)=80.d0
        case(36) ; pp%atom_symbol(ik) = 'Kr' ; pp%rmass(ik)=84.d0
        case(37) ; pp%atom_symbol(ik) = 'Rb' ; pp%rmass(ik)=85.d0
        case(38) ; pp%atom_symbol(ik) = 'Sr' ; pp%rmass(ik)=88.d0
        case(39) ; pp%atom_symbol(ik) = 'Y ' ; pp%rmass(ik)=89.d0
        case(40) ; pp%atom_symbol(ik) = 'Zr' ; pp%rmass(ik)=91.d0
        case(41) ; pp%atom_symbol(ik) = 'Nb' ; pp%rmass(ik)=93.d0
        case(42) ; pp%atom_symbol(ik) = 'Mo' ; pp%rmass(ik)=96.d0
        case(43) ; pp%atom_symbol(ik) = 'Tc' ; pp%rmass(ik)=98.d0
        case(44) ; pp%atom_symbol(ik) = 'Ru' ; pp%rmass(ik)=101.d0
        case(45) ; pp%atom_symbol(ik) = 'Rh' ; pp%rmass(ik)=103.d0
        case(46) ; pp%atom_symbol(ik) = 'Pd' ; pp%rmass(ik)=106.d0
        case(47) ; pp%atom_symbol(ik) = 'Ag' ; pp%rmass(ik)=108.d0
        case(48) ; pp%atom_symbol(ik) = 'Cd' ; pp%rmass(ik)=112.d0
        case(49) ; pp%atom_symbol(ik) = 'In' ; pp%rmass(ik)=115.d0
        case(50) ; pp%atom_symbol(ik) = 'Sn' ; pp%rmass(ik)=119.d0
        case(51) ; pp%atom_symbol(ik) = 'Sb' ; pp%rmass(ik)=122.d0
        case(52) ; pp%atom_symbol(ik) = 'Te' ; pp%rmass(ik)=128.d0
        case(53) ; pp%atom_symbol(ik) = 'I ' ; pp%rmass(ik)=127.d0
        case(54) ; pp%atom_symbol(ik) = 'Xe' ; pp%rmass(ik)=131.d0
        case(55) ; pp%atom_symbol(ik) = 'Cs' ; pp%rmass(ik)=133.d0
        case(56) ; pp%atom_symbol(ik) = 'Ba' ; pp%rmass(ik)=137.d0
        case(57) ; pp%atom_symbol(ik) = 'La' ; pp%rmass(ik)=139.d0
        case(58) ; pp%atom_symbol(ik) = 'Ce' ; pp%rmass(ik)=140.d0
        case(59) ; pp%atom_symbol(ik) = 'Pr' ; pp%rmass(ik)=141.d0
        case(60) ; pp%atom_symbol(ik) = 'Nd' ; pp%rmass(ik)=144.d0
        case(61) ; pp%atom_symbol(ik) = 'Pm' ; pp%rmass(ik)=145.d0
        case(62) ; pp%atom_symbol(ik) = 'Sm' ; pp%rmass(ik)=150.d0
        case(63) ; pp%atom_symbol(ik) = 'Eu' ; pp%rmass(ik)=152.d0
        case(64) ; pp%atom_symbol(ik) = 'Gd' ; pp%rmass(ik)=157.d0
        case(65) ; pp%atom_symbol(ik) = 'Tb' ; pp%rmass(ik)=159.d0
        case(66) ; pp%atom_symbol(ik) = 'Dy' ; pp%rmass(ik)=164.d0
        case(67) ; pp%atom_symbol(ik) = 'Ho' ; pp%rmass(ik)=165.d0
        case(68) ; pp%atom_symbol(ik) = 'Er' ; pp%rmass(ik)=167.d0
        case(69) ; pp%atom_symbol(ik) = 'Tm' ; pp%rmass(ik)=169.d0
        case(70) ; pp%atom_symbol(ik) = 'Yb' ; pp%rmass(ik)=173.d0
        case(71) ; pp%atom_symbol(ik) = 'Lu' ; pp%rmass(ik)=175.d0
        case(72) ; pp%atom_symbol(ik) = 'Hf' ; pp%rmass(ik)=178.d0
        case(73) ; pp%atom_symbol(ik) = 'Ta' ; pp%rmass(ik)=181.d0
        case(74) ; pp%atom_symbol(ik) = 'W ' ; pp%rmass(ik)=184.d0
        case(75) ; pp%atom_symbol(ik) = 'Re' ; pp%rmass(ik)=186.d0
        case(76) ; pp%atom_symbol(ik) = 'Os' ; pp%rmass(ik)=190.d0
        case(77) ; pp%atom_symbol(ik) = 'Ir' ; pp%rmass(ik)=192.d0
        case(78) ; pp%atom_symbol(ik) = 'Pt' ; pp%rmass(ik)=195.d0
        case(79) ; pp%atom_symbol(ik) = 'Au' ; pp%rmass(ik)=197.d0
        case(80) ; pp%atom_symbol(ik) = 'Hg' ; pp%rmass(ik)=201.d0
        case(81) ; pp%atom_symbol(ik) = 'Tl' ; pp%rmass(ik)=204.d0
        case(82) ; pp%atom_symbol(ik) = 'Pb' ; pp%rmass(ik)=207.d0
        case(83) ; pp%atom_symbol(ik) = 'Bi' ; pp%rmass(ik)=209.d0
        case default ; stop 'Unprepared atomic data is called input_pseudopotential_YS'
        end select

      end do

    end if

    call comm_bcast(pp%atom_symbol,nproc_group_global)
    call comm_bcast(pp%rmass,nproc_group_global)
    call comm_bcast(ps_format,nproc_group_global)

    pp%lmax0=lmax0

    pp%nrmax0=nrmax0
  
    pp%nrmax=nrmax
    pp%lmax=lmax
  
    allocate(pp%lref(1:nelem))
    pp%lref(1:nelem)=lloc_ps(1:nelem)

    allocate(pp%nrps(1:nelem))
    allocate(pp%rps(1:nelem))
    allocate(pp%mlps(1:nelem))
    allocate(pp%nproj(0:lmax,1:nelem)); pp%nproj=0
    allocate(pp%zps(1:nelem))
    allocate(pp%nrloc(1:nelem))
    allocate(pp%rloc(1:nelem))

! # of angular momenta is lmax+1,
! and the maximum # of projector for each angular momentum is assumed to 2
! so that the maximum number of projectors becomes 2*(lmax+1).
    
    allocate(pp%anorm(0:2*lmax+1,nelem)); pp%anorm=0.0d0
    allocate(pp%inorm(0:2*lmax+1,nelem)); pp%inorm=0

    allocate(pp%anorm_so(0:2*lmax+1,nelem)); pp%anorm_so=0.0d0
    allocate(pp%inorm_so(0:2*lmax+1,nelem)); pp%inorm_so=0
  
    allocate(pp%rad(nrmax,nelem))
    allocate(pp%radnl(nrmax,nelem))
    
    allocate(pp%vloctbl(nrmax,nelem))
    allocate(pp%dvloctbl(nrmax,nelem))
    allocate(pp%udvtbl(nrmax,0:2*lmax+1,nelem)); pp%udvtbl=0.0d0
    allocate(pp%dudvtbl(nrmax,0:2*lmax+1,nelem)); pp%dudvtbl=0.0d0
    
    allocate(pp%rho_nlcc_tbl(nrmax,nelem))
    allocate(pp%tau_nlcc_tbl(nrmax,nelem))
  
    allocate(pp%vpp(0:nrmax0,0:2*lmax0+2),pp%upp(0:nrmax0,0:2*lmax0+1))
    allocate(pp%dvpp(0:nrmax0,0:2*lmax0+2),pp%dupp(0:nrmax0,0:2*lmax0+1))
    allocate(pp%vpp_f(0:nrmax0,0:2*lmax0+2,nelem),pp%upp_f(0:nrmax0,0:2*lmax0+1,nelem))

    allocate( pp%vpp_so(0:nrmax0,0:2*lmax0+2) ); pp%vpp_so=0.0d0
    allocate( pp%dvpp_so(0:nrmax0,0:2*lmax0+2) ); pp%dvpp_so=0.0d0
    allocate( pp%vpp_f_so(0:nrmax0,0:2*lmax0+2,nelem) ); pp%vpp_f_so=0.0d0
  
    pp%flag_nlcc=flag_nlcc
  
  end subroutine init_pp
!======================================================================
  subroutine read_mr_yb(pp,ik,ps_file)
    use structures, only: s_pp_info
    implicit none
    type(s_pp_info) :: pp
    integer :: ik
    character(256) :: ps_file
    
    open(4,file=ps_file,status='old')
    read(4,*) pp%mr(ik)
    close(4)
    return
  
  end subroutine read_mr_YB
!======================================================================
  subroutine read_mr_abinit(pp,ik,ps_file)
    use structures, only: s_pp_info
    implicit none
    type(s_pp_info) :: pp
    integer :: ik
    character(256) :: ps_file
    real(8) :: zatom, zion, pspdat,pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
    character(1) :: dummy_text
    
    open(4,file=ps_file,status='old')
    read(4,*) dummy_text
    read(4,*) zatom, zion, pspdat
    read(4,*) pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
    close(4)
    
    pp%mr(ik)=mmax
  
  end subroutine read_mr_ABINIT
!======================================================================
  subroutine read_mr_abinitfhi(pp,ik,ps_file)
    use structures, only: s_pp_info
    implicit none
    type(s_pp_info) :: pp
    integer :: ik
    integer :: i
    character(256) :: ps_file
    character(1) :: dummy_text
    
    open(4,file=ps_file,status='old')
    do i=1,18
      read(4,*) dummy_text
    end do
    read(4,*) pp%mr(ik)
    close(4)
    
  end subroutine read_mr_abinitfhi

!======================================================================
  subroutine read_mr_fhi(pp,ik,ps_file)
    use structures, only: s_pp_info
    implicit none
    type(s_pp_info) :: pp
    integer :: ik
    integer :: i
    character(256) :: ps_file
    character(1) :: dummy_text
    
    open(4,file=ps_file,status='old')
    do i=1,11
      read(4,*) dummy_text
    end do
    read(4,*) pp%mr(ik)
    close(4)
    
  end subroutine read_mr_fhi

!======================================================================
  subroutine calc_nlcc(pp, sys, rg, ppn)
    use salmon_global,only : kion
    use structures, only : s_dft_system, s_pp_info, s_pp_nlcc, s_rgrid
    implicit none
    
    type(s_dft_system), intent(in) :: sys
    type(s_rgrid), intent(in) :: rg
    type(s_pp_info), intent(in) :: pp
    type(s_pp_nlcc), intent(inout) :: ppn
  
    integer :: a, ik, ir, intr
    integer :: i, i1, i2, i3, j1, j2, j3
    integer :: irepr_min, irepr_max
    real(8) :: rion_repr(3)
    real(8) :: r, rc, r1, r2, r3
    real(8) :: ratio1, ratio2

    if(allocated(ppn%rho_nlcc)) deallocate(ppn%rho_nlcc,ppn%tau_nlcc)
 
    ! Allocate
    allocate(ppn%rho_nlcc( &
      & rg%is(1):rg%ie(1), &
      & rg%is(2):rg%ie(2), &
      & rg%is(3):rg%ie(3)))
    allocate(ppn%tau_nlcc( &
      & rg%is(1):rg%ie(1), &
      & rg%is(2):rg%ie(2), &
      & rg%is(3):rg%ie(3)))

    ppn%rho_nlcc = 0d0
    ppn%tau_nlcc = 0d0  
  
    if (sys%iperiodic == 0) then
      irepr_min = 0
      irepr_max = 0
    elseif (sys%iperiodic == 3) then
      irepr_min = -2
      irepr_max = +2
    else
      stop "Sorry, not implemented (calc_nlcc@prep_pp.f90)"
    endif
    
    if (.not. pp%flag_nlcc) return ! Do nothing
  
    do a=1, sys%nion
      ik = Kion(a)
      rc = 15d0 ! maximum
      do i=1, pp%nrmax
        if(pp%rho_nlcc_tbl(i,ik) + pp%tau_nlcc_tbl(i,ik) < 1d-6)then
          rc = pp%rad(i,ik)
          exit
        end if
        if(i == pp%nrmax) stop "no-cut-off found (calc_nlcc@prep_pp.f90)"
      end do
  
      do i1 = irepr_min, irepr_max
      do i2 = irepr_min, irepr_max
      do i3 = irepr_min, irepr_max
        rion_repr(1) = sys%rion(1, a) + i1 * sys%primitive_a(1, 1)
        rion_repr(2) = sys%rion(2, a) + i2 * sys%primitive_a(2, 2)
        rion_repr(3) = sys%rion(3, a) + i3 * sys%primitive_a(3, 3)
        do j1 = rg%is(1), rg%ie(1)
        do j2 = rg%is(2), rg%ie(2)
        do j3 = rg%is(3), rg%ie(3)
          r1 = (j1-1) * sys%hgs(1) - rion_repr(1) ! iwata
          r2 = (j2-1) * sys%hgs(2) - rion_repr(2) ! iwata
          r3 = (j3-1) * sys%hgs(3) - rion_repr(3) ! iwata
          r = sqrt(r1**2 + r2**2 + r3**2)
          if (r <= rc) then
            do ir = 1, pp%nrmax ! iwata
              if (pp%rad(ir,ik) .gt. r) exit
            end do
            intr = ir - 1
            if (intr.lt.0.or.intr.ge.pp%NRmax) stop 'bad intr at prep_ps'
            ratio1=(r-pp%rad(intr,ik))/(pp%rad(intr+1,ik)-pp%rad(intr,ik))
            ratio2=1-ratio1
            ppn%rho_nlcc(j1, j2, j3) = ppn%rho_nlcc(j1, j2, j3) & ! iwata
              +ratio1*pp%rho_nlcc_tbl(intr+1,ik)+ratio2*pp%rho_nlcc_tbl(intr,ik)
            ppn%tau_nlcc(j1, j2, j3) = ppn%tau_nlcc(j1, j2, j3) & ! iwata
              +ratio1*pp%tau_nlcc_tbl(intr+1,ik)+ratio2*pp%tau_nlcc_tbl(intr,ik)
          end if
        end do
        end do
        end do
      end do
      end do
      end do
    end do
  
    return
  end subroutine calc_nlcc
  

end module salmon_pp
