!
!  Copyright 2018-2020 SALMON developers
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
module input_pp_sub
  implicit none

  logical,private :: flag_potential_is_given=.true.  ! Pseudopotential data is given as the potentitals and the wavefunctions
  logical,private :: flag_beta_proj_is_given=.false. ! Pseudopotential data is given as projection operators
  logical,private :: flag_so=.false.

contains

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

subroutine input_pp(pp,hx,hy,hz)
  use structures,only : s_pp_info
  use salmon_global,only : file_pseudo, quiet, method_init_density
  use salmon_global,only : n_Yabana_Bertsch_psformat,n_ABINIT_psformat&
    &,n_ABINITFHI_psformat,n_FHI_psformat,ps_format,nelem,base_directory, &
    & yn_psmask
  use parallelization, only: nproc_group_global, nproc_id_global
  use communication, only: comm_bcast, comm_is_root
  use math_constants, only : pi
  use read_ps_upf_module, only: read_ps_upf
  use read_paw_upf_module, only: read_paw_upf
  implicit none
  type(s_pp_info) :: pp
  real(8),parameter :: Eps0=1d-10
  real(8),intent(in) :: hx,hy,hz
  integer :: ik,l,i,l0,ll,nprj,i1,nr
  real(8) :: rrc(0:pp%lmax0)
  real(8) :: r1
  real(8),allocatable :: rhor_nlcc(:,:)   !zero in radial index for taking derivative
  character(256) :: ps_file
  logical,allocatable :: flag_nlcc_element(:)

  allocate(rhor_nlcc(0:pp%nrmax0,0:2))
  rhor_nlcc=0d0

! Nonlinear core correction
  allocate(flag_nlcc_element(nelem)); flag_nlcc_element(:) = .false. ; pp%flag_nlcc = .false.

!      ps_file=trim(base_directory)//trim(pp%atom_symbol(ik))//trim(ps_postfix)

  rrc=0.0d0

  if (comm_is_root(nproc_id_global)) then

    do ik=1,nelem
         
      ps_file = trim(file_pseudo(ik))

      select case (ps_format(ik))
      case('KY')
        call read_ps_ky(pp,rrc,ik,ps_file)
      case('ABINIT')
        call read_ps_abinit(pp,rrc,ik,ps_file)
      case('ABINITFHI')
        call read_ps_abinitfhi(pp,rrc,rhor_nlcc,flag_nlcc_element,ik,ps_file)
      case('ABINITPSP8')
        call read_ps_abinitpsp8(pp,rrc,rhor_nlcc,flag_nlcc_element,ik,ps_file)
      case('FHI')
        call read_ps_fhi(pp,rrc,ik,ps_file)
      case('ADPACK')
        call read_ps_adpack(pp,rrc,rhor_nlcc,flag_nlcc_element,ik,ps_file)
      case('UPF')
        if( index(ps_file,'paw')>0 .or. index(ps_file,'PAW')>0 )then
          open(4,file=ps_file,status='old')
          call read_paw_upf(4,pp,ik)
          close(4)
        else
          call read_ps_upf(pp,rrc,rhor_nlcc,flag_nlcc_element,ik,ps_file)
        end if
        flag_beta_proj_is_given =.true.
!      case('ATOM')      ; call read_ps_ATOM
      case default ; stop 'Unprepared ps_format is required input_pseudopotential_YS'
      end select

      if ( flag_beta_proj_is_given ) then
        flag_potential_is_given=.false.
        if(method_init_density/='wf' .and. ps_format(ik)/='UPF') then
          stop "radial density is not available (method_init_density=pp...)"
        end if
      end if
      if ( any(pp%vpp_so/=0.0d0) ) flag_so=.true.

      if ( all(pp%nproj(:,ik)==0) ) pp%nproj(0:pp%mlps(ik),ik)=1

! outside mr (needed for isolated systems)
      if(pp%nrmax>pp%mr(ik))then
        do i=pp%mr(ik)+1,pp%nrmax
          if ( pp%rad(i,ik)>0.0d0 ) pp%vpp(i,0:pp%mlps(ik))=-pp%zps(ik)/pp%rad(i,ik) !iwata
          pp%upp(i,0:pp%mlps(ik))=0.d0
        end do
      end if

! Set meaning domain in the arrays 
      if ( any(rrc/=0.0d0) ) then
        pp%rps(ik)=maxval(rrc(0:pp%mlps(ik)))
      end if
      do i=1,pp%nrmax
        if(pp%rad(i,ik).gt.pp%rps(ik)) exit
      enddo
      pp%nrps(ik)=i
      if(pp%nrps(ik).ge.pp%nrmax) stop 'NRps>Nrmax at input_pseudopotential_YS'
      pp%nrloc(ik)=pp%nrps(ik)
      pp%rloc(ik)=pp%rps(ik)
      pp%radnl(:,ik)=pp%rad(:,ik)

! Set meaning domain in the arrays of radial wave functions
      i1=0
      r1=0.0d0
      nr=ubound(pp%upp,1)
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        do i=nr,1,-1
          if ( abs(pp%upp(i,l)) > 1.d-10 ) then
            r1=max( r1, pp%rad(i+1,ik) )
            i1=max(i+1,i1)
            exit
          end if
        end do
      end do
      l0=l
      end do
!      pp%rps_ao(ik)=r1
!      pp%nrps_ao(ik)=i1
      pp%rps_ao(ik)=pp%rps(ik)
      pp%nrps_ao(ik)=pp%nrps(ik)

      if( flag_beta_proj_is_given )then
        l0=0
        do ll=0,pp%mlps(ik)
        do l=l0,l0+pp%nproj(ll,ik)-1
          pp%inorm(l,ik) = 1
          if( pp%anorm(l,ik) < 0.0d0 )then
            pp%anorm(l,ik) = -pp%anorm(l,ik)
            pp%inorm(l,ik) = -1
          end if
          if ( abs(pp%anorm(l,ik)) < Eps0 ) pp%inorm(l,ik)=0
          pp%anorm(l,ik) = sqrt( pp%anorm(l,ik) )
        end do
        l0=l
        end do
        if ( flag_so ) then
          l0=0
          do ll=0,pp%mlps(ik)
          do l=l0,l0+pp%nproj(ll,ik)-1
            if ( pp%anorm_so(l,ik) == 0.0d0 ) cycle
            pp%inorm_so(l,ik) = 1
            if( pp%anorm_so(l,ik) < 0.0d0 )then
              pp%anorm_so(l,ik) = -pp%anorm_so(l,ik)
              pp%inorm_so(l,ik) = -1
            end if
            if ( abs(pp%anorm_so(l,ik)) < Eps0 ) pp%inorm_so(l,ik)=0
            pp%anorm_so(l,ik) = sqrt( pp%anorm_so(l,ik) )
          end do
          l0=l
          end do
        end if
      else
        do l=0,pp%mlps(ik)
          pp%anorm(l,ik) = 0.d0
          do i=1,pp%mr(ik)-1
            r1 = pp%rad(i+1,ik)-pp%rad(i,ik)
            pp%anorm(l,ik) = pp%anorm(l,ik)  &
                 + (pp%upp(i,l)**2*(pp%vpp(i,l)-pp%vpp(i,pp%lref(ik)))+pp%upp(i-1,l)**2*(pp%vpp(i-1,l)-pp%vpp(i-1,pp%lref(ik))))*r1
          end do
          pp%anorm(l,ik) = 0.5d0*pp%anorm(l,ik)
          pp%inorm(l,ik)=+1
          if(pp%anorm(l,ik).lt.0.d0) then
            pp%anorm(l,ik)=-pp%anorm(l,ik)
            pp%inorm(l,ik)=-1
          endif
          if(abs(pp%anorm(l,ik)).lt.Eps0)pp%inorm(l,ik)=0
          pp%anorm(l,ik)=sqrt(pp%anorm(l,ik))
        enddo
      end if

! number of non-local projectors
      if ( all(pp%nproj(:,ik)==0) ) pp%nproj(0:pp%mlps(ik),ik)=1

      if (.not. quiet) then
      write(*,*) '===================pseudopotential data==================='
      write(*,*) 'ik ,atom_symbol=',ik, pp%atom_symbol(ik)
      write(*,*) 'ps_format =',ps_format(ik)
      write(*,*) 'ps_file =',trim(ps_file)
      write(*,*) 'Zps(ik), Mlps(ik) =',pp%zps(ik), pp%mlps(ik)
      write(*,*) 'Rps(ik), NRps(ik) =',pp%rps(ik), pp%nrps(ik)
      write(*,*) 'Lref(ik) =',pp%lref(ik)
      write(*,*) 'nproj(ik,l) =',(pp%nproj(l,ik),l=0,pp%mlps(ik))
      write(*,*) 'anorm(ik,l) =',(real(pp%anorm(l,ik)),l=0,sum(pp%nproj(:,ik))-1)
      write(*,*) 'inorm(ik,l) =',(pp%inorm(l,ik),l=0,sum(pp%nproj(:,ik))-1)
      if( flag_so )then
        write(*,*) 'anorm_so(ik,l) =',(real(pp%anorm_so(l,ik)),l=0,sum(pp%nproj(:,ik))-1)
        write(*,*) 'inorm_so(ik,l) =',(pp%inorm_so(l,ik),l=0,sum(pp%nproj(:,ik))-1)
      end if
      write(*,*) 'Mass(ik) =',pp%rmass(ik)
      write(*,*) 'flag_nlcc_element(ik) =',flag_nlcc_element(ik)
      write(*,*) '=========================================================='
      end if

      if (yn_psmask == 'y') then
        call making_ps_with_masking(pp,hx,hy,hz,ik, &
                                    rhor_nlcc,flag_nlcc_element)
        if ( flag_so ) call making_ps_with_masking_so( pp,hx,hy,hz,ik )
        if (.not. quiet) then
        write(*,*) 'Following quantities are modified by masking procedure'
        write(*,*) 'Rps(ik), NRps(ik) =',pp%rps(ik), pp%nrps(ik)
        write(*,*) 'anorm(ik,l) =',(pp%anorm(l,ik),l=0,pp%mlps(ik))
        write(*,*) 'inorm(ik,l) =',(pp%inorm(l,ik),l=0,pp%mlps(ik))
        end if
      else if (yn_psmask == 'n') then
        call making_ps_without_masking(pp,ik,flag_nlcc_element,rhor_nlcc)
        if ( flag_so ) call making_ps_without_masking_so( pp, ik )
      else
        stop 'Wrong yn_psmask at input_pseudopotential_YS'
      end if

      pp%upp_f(:,:,ik)=pp%upp(:,:)
      pp%vpp_f(:,:,ik)=pp%vpp(:,:)

      open(4,file=trim(base_directory)//"PS_"//trim(pp%atom_symbol(ik))//"_"//trim(ps_format(ik))//"_"//trim(yn_psmask)//".dat")
      write(4,*) "# Mr=",pp%mr(ik)
      write(4,*) "# Rps(ik), NRps(ik)",pp%rps(ik), pp%nrps(ik)
      write(4,*) "# Mlps(ik), Lref(ik) =",pp%mlps(ik), pp%lref(ik)
      nprj=sum(pp%nproj(:,ik))
      do i=1,pp%nrps(ik)
        write(4,'(30e21.12)') pp%rad(i,ik),(pp%udvtbl(i,l,ik),l=0,nprj-1),(pp%dudvtbl(i,l,ik),l=0,nprj-1)
      end do
      close(4)

    enddo
  endif

  call comm_bcast(pp%mr,nproc_group_global)
  call comm_bcast(pp%zps,nproc_group_global)
  call comm_bcast(pp%mlps,nproc_group_global)
  call comm_bcast(pp%lref,nproc_group_global)
  call comm_bcast(pp%nproj,nproc_group_global)
  call comm_bcast(pp%rps,nproc_group_global)
  call comm_bcast(pp%nrps,nproc_group_global)
  call comm_bcast(pp%nrloc,nproc_group_global)
  call comm_bcast(pp%rloc,nproc_group_global)
  call comm_bcast(pp%anorm,nproc_group_global)
  call comm_bcast(pp%inorm,nproc_group_global)
  call comm_bcast(pp%rad,nproc_group_global)
  call comm_bcast(pp%radnl,nproc_group_global)
  call comm_bcast(pp%vloctbl,nproc_group_global)
  call comm_bcast(pp%dvloctbl,nproc_group_global)
  call comm_bcast(pp%udvtbl,nproc_group_global)
  call comm_bcast(pp%dudvtbl,nproc_group_global)
  call comm_bcast(pp%upp_f,nproc_group_global)
  call comm_bcast(pp%vpp_f,nproc_group_global)
  call comm_bcast(pp%rho_pp_tbl,nproc_group_global)
  call comm_bcast(pp%rho_nlcc_tbl,nproc_group_global)
  call comm_bcast(pp%tau_nlcc_tbl,nproc_group_global)
  call comm_bcast(pp%flag_nlcc,nproc_group_global)
  if ((.not. quiet) .and. comm_is_root(nproc_id_global)) &
    & write(*,*)"flag_nlcc = ",pp%flag_nlcc

  call comm_bcast(pp%anorm_so,nproc_group_global)
  call comm_bcast(pp%inorm_so,nproc_group_global)
  call comm_bcast(pp%udvtbl_so,nproc_group_global)
  call comm_bcast(pp%dudvtbl_so,nproc_group_global)
  
  call comm_bcast(pp%rps_ao,nproc_group_global)
  call comm_bcast(pp%nrps_ao,nproc_group_global)
  call comm_bcast(pp%upptbl_ao,nproc_group_global)

  return
end subroutine input_pp
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine read_ps_ky(pp,rrc,ik,ps_file)
  use structures,only : s_pp_info
  use salmon_global,only : Lmax_ps
  use inputoutput,only : au_length_aa, au_energy_ev
  implicit none
  type(s_pp_info),intent(inout) :: pp
!argument
  integer,intent(in) :: ik
  real(8),intent(out) :: rrc(0:pp%lmax0)
  character(256),intent(in) :: ps_file
!local variable
  integer :: l,i,irPC
  real(8) :: step,rPC,r,rhopp(0:pp%nrmax0),rzps

  open(4,file=ps_file,status='old')
  read(4,*) pp%mr(ik),step,pp%mlps(ik),rzps
  pp%zps(ik)=int(rzps+1d-10)
  if(pp%mr(ik) .gt.pp%nrmax0) stop 'Mr>Nrmax0 at Read_PS_KY'
  if(pp%mlps(ik).gt.pp%lmax0) stop 'Mlps(ik)>Lmax0 at Read_PS_KY'
  if(pp%mlps(ik).gt.pp%lmax)  stop 'Mlps(ik)>Lmax at Read_PS_KY'
  read(4,*) irPC,(rrc(l),l=0,pp%mlps(ik))
  rPC=real(irPC) !Radius for partial core correction: not working in this version
  do i=0,pp%mr(ik)
    read(4,*) r,rhopp(i),(pp%vpp(i,l),l=0,pp%mlps(ik))
  end do
  do i=0,pp%mr(ik)
    read(4,*) r,(pp%upp(i,l),l=0,pp%mlps(ik))
  end do
  close(4)

! change to atomic unit
  step=step/au_length_aa
  rrc(0:pp%mlps(ik))=rrc(0:pp%mlps(ik))/au_length_aa
  pp%vpp(0:pp%mr(ik),0:pp%mlps(ik))=pp%vpp(0:pp%mr(ik),0:pp%mlps(ik))/au_energy_ev
  pp%upp(0:pp%mr(ik),0:pp%mlps(ik))=pp%upp(0:pp%mr(ik),0:pp%mlps(ik))*sqrt(au_length_aa)

  do i=1,pp%nrmax
    pp%rad(i,ik)=(i-1)*step  !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
  enddo

  if(Lmax_ps(ik) >= 0)pp%mlps(ik) = Lmax_ps(ik) ! Maximum angular momentum given by input

  !do i = 0, pp%mr(ik)
  !  write(100,'(5g15.6)') pp%rad(i+1,ik), (pp%upp(i,l),l=0,pp%mlps(ik) )
  !end do

  !rewind 100
  !do i = 2, pp%mr(ik)+1
  !  r = pp%rad(i,ik)
  !  write(100,'(5g15.6)') pp%rad(i,ik), pp%rho_pp_tbl(i,ik)/(4.0d0*acos(-1.0d0)*r*r)
  !end do

  return
end subroutine read_ps_KY
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine read_ps_abinit(pp,rrc,ik,ps_file)
  use structures,only : s_pp_info
  use salmon_global,only : Lmax_ps
!See http://www.abinit.org/downloads/psp-links/psp-links/lda_tm
  implicit none
  type(s_pp_info),intent(inout) :: pp
!argument
  integer,intent(in) :: ik
  real(8),intent(out) :: rrc(0:pp%lmax0)
  character(256),intent(in) :: ps_file
!local variable
  integer :: i
  real(8) :: rzps
  integer :: ll
  real(8) :: zatom, pspdat,pspcod,pspxc,lmaxabinit,lloc,mmax,r2well,l
  real(8) :: e99_0,e99_9,nproj,rcpsp,rms,ekb1,ekb2,epsatm,rchrg,fchrg,qchrg
  character(1) :: dummy_text

  open(4,file=ps_file,status='old')
  read(4,*) dummy_text
  read(4,*) zatom, pp%zion, pspdat
  rzps = pp%zion
  pp%zps(ik)=int(rzps+1d-10)
  read(4,*) pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
  pp%mlps(ik)=lmaxabinit
  if(lloc .ne. pp%lref(ik)) write(*,*) "Warning! Lref(ik=",ik,") is different from intended one in ",ps_file
  pp%mr(ik) = mmax - 1
  if(pp%mr(ik).gt.pp%nrmax0) stop 'Mr>Nrmax0 at Read_PS_ABINIT'
  if(pp%mlps(ik).gt.pp%lmax0) stop 'Mlps(ik)>Lmax0 at Read_PS_ABINIT'
  if(pp%mlps(ik).gt.pp%lmax) stop 'Mlps(ik)>Lmax at Read_PS_ABINIT'
  do ll=0,pp%mlps(ik)
    read(4,*) l,e99_0,e99_9,nproj,rcpsp
    read(4,*) rms,ekb1,ekb2,epsatm
    rrc(ll) = rcpsp
  end do
  read(4,*) rchrg,fchrg,qchrg
  do ll=0,pp%mlps(ik)
    read(4,*) dummy_text
    do i=1,(pp%mr(ik)+1)/3
      read(4,*) pp%vpp(3*(i-1),ll),pp%vpp(3*(i-1)+1,ll),pp%vpp(3*(i-1)+2,ll)
    end do
  end do
  do ll=0,pp%mlps(ik)
    read(4,*) dummy_text
    do i=1,(pp%mr(ik)+1)/3
      read(4,*) pp%upp(3*(i-1),ll),pp%upp(3*(i-1)+1,ll),pp%upp(3*(i-1)+2,ll)
    end do
  end do
  close(4)

  do i=0,pp%nrmax-1
    pp%rad(i+1,ik) = 1.0d2*(dble(i)/dble(mmax-1)+1.0d-2)**5 - 1.0d-8
  end do

  if(Lmax_ps(ik) >= 0)pp%mlps(ik) = Lmax_ps(ik) ! Maximum angular momentum given by input

  return
end subroutine read_ps_abinit
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine read_ps_abinitfhi(pp,rrc,rhor_nlcc,flag_nlcc_element,ik,ps_file)
!This is for  FHI pseudopotential listed in abinit web page and not for original FHI98PP.
!See http://www.abinit.org/downloads/psp-links/lda_fhi
  use structures,only : s_pp_info
  use salmon_global,only : nelem, Lmax_ps, quiet
  use math_constants,only : pi
  implicit none
  type(s_pp_info),intent(inout) :: pp
!argument
  integer,intent(in) :: ik
  real(8),intent(out) :: rrc(0:pp%lmax0)
  character(256),intent(in) :: ps_file
  real(8),intent(out) :: rhor_nlcc(0:pp%nrmax0,0:2)
  logical,intent(inout) :: flag_nlcc_element(nelem)
!local variable
  character(50) :: temptext
  integer :: i,j
  real(8) :: step,rzps,dummy
  integer :: mr_l(0:pp%lmax0),l,ll
  real(8) :: step_l(0:pp%lmax),rrc_mat(0:pp%lmax,0:pp%lmax)

  rrc_mat(0:pp%lmax,0:pp%lmax)=-1.d0

  open(4,file=ps_file,status='old')
  if (.not. quiet) &
  write(*,*) '===================Header of ABINITFHI pseudo potential==================='
  do i=1,7
    read(4,'(a)') temptext
    write(*,*) temptext
  end do
  if (.not. quiet) &
  write(*,*) '===================Header of ABINITFHI pseudo potential==================='
  read(4,*) rzps,pp%mlps(ik)
  pp%zps(ik)=int(rzps+1d-10)
  pp%mlps(ik) = pp%mlps(ik)-1
  if(pp%mlps(ik).gt.pp%lmax0) stop 'Mlps(ik)>Lmax0 at Read_PS_FHI'
  if(pp%mlps(ik).gt.pp%lmax) stop 'Mlps(ik)>Lmax at Read_PS_FHI'
  do i=1,10
     read(4,*) dummy
  end do
  do l=0,pp%mlps(ik)
    read(4,*) mr_l(l),step_l(l)
    if(mr_l(l).gt.pp%nrmax0) stop 'Mr>Nrmax0 at Read_PS_FHI'
    do i=1,mr_l(l)
      read(4,*) j,pp%rad(i+1,ik),pp%upp(i,l),pp%vpp(i,l) !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
    end do
    pp%rad(1,ik)=0.d0
    pp%upp(0,l)=0.d0
    pp%vpp(0,l)=pp%vpp(1,l)-(pp%vpp(2,l)-pp%vpp(1,l))/(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik))
  end do
!Nonlinear core-correction
  do i=1,mr_l(0)
    read(4,*,end=940) pp%rad(i+1,ik),rhor_nlcc(i,0),rhor_nlcc(i,1),rhor_nlcc(i,2)
  end do
  rhor_nlcc(0,:)=rhor_nlcc(1,:)-(rhor_nlcc(2,:)-rhor_nlcc(1,:)) &
    /(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik))
  rhor_nlcc = rhor_nlcc/(4d0*pi)
  flag_nlcc_element(ik) = .true.
!Nonlinear core-correction


940 close(4)
  

  if(minval(mr_l(0:pp%mlps(ik))).ne.maxval(mr_l(0:pp%mlps(ik)))) then
    stop 'Mr are diffrent at Read_PS_FHI'
  else 
    pp%mr(ik) = minval(mr_l(0:pp%mlps(ik)))
  end if
  if((maxval(step_l(0:pp%mlps(ik)))-minval(step_l(0:pp%mlps(ik)))).ge.1.d-14) then
    stop 'step are different at Read_PS_FHI'
  else 
    step = minval(step_l(0:pp%mlps(ik)))
  end if

  do i=pp%mr(ik)+1,pp%nrmax-1
    pp%rad(i+1,ik) = pp%rad(i,ik)*step
  end do

  do l=0,pp%mlps(ik)
    do ll=0,pp%mlps(ik)
      do i=pp%mr(ik),1,-1
        if(abs(pp%vpp(i,l)-pp%vpp(i,ll)).gt.1.d-10) then
          rrc_mat(l,ll) = pp%rad(i+1+1,ik)
          exit
        end if
      end do
    end do
    rrc(l)=maxval(rrc_mat(l,:))
  end do

  if(Lmax_ps(ik) >= 0)pp%mlps(ik) = Lmax_ps(ik) ! Maximum angular momentum given by input

  return
end subroutine read_ps_abinitfhi
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine read_ps_abinitpsp8(pp,rrc,rhor_nlcc,flag_nlcc_element,ik,ps_file)
  use structures,only : s_pp_info
  use salmon_global,only : nelem
  implicit none
  type(s_pp_info),intent(inout) :: pp
  real(8),intent(out) :: rhor_nlcc(0:pp%nrmax0,0:2)
  logical,intent(inout) :: flag_nlcc_element(nelem)
!argument
  integer,intent(in) :: ik
  real(8),intent(out) :: rrc(0:pp%lmax0)
  character(256),intent(in) :: ps_file
!local variable
  integer :: i,j,j0,ll,nproj,lloc,lloc_check,lmaxabinit
  real(8) :: rzps,rcpsp,pi4
  real(8) :: zatom, pspdat,pspcod,pspxc,mmax,r2well,dr
  real(8) :: rchrg,fchrg,qchrg,sum_rho_pp,sum_rho_nlcc,rho_tmp,r_tmp
  character(1) :: dummy_text

  flag_beta_proj_is_given = .true.

  open(4,file=ps_file,status='old')
! header
  read(4,*) dummy_text
  read(4,*) zatom, pp%zion, pspdat
  rzps = pp%zion
  pp%zps(ik)=int(rzps+1d-10)
  read(4,*) pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
  pp%mlps(ik)=lmaxabinit
  if (lloc /= pp%lref(ik)) write(*,*) "Warning! Lref(ik=",ik,") is different from intended one in ",ps_file
  pp%mr(ik) = mmax - 1
  if ( pp%mr(ik)   > pp%nrmax0) stop 'Mr>Nrmax0 at Read_PS_ABINIT'
  if ( pp%mlps(ik) > pp%lmax0 ) stop 'Mlps(ik)>Lmax0 at Read_PS_ABINIT'
  if ( pp%mlps(ik) > pp%lmax  ) stop 'Mlps(ik)>Lmax at Read_PS_ABINIT'
  read(4,*) rchrg,fchrg,qchrg
  read(4,*) ( pp%nproj(ll,ik), ll=0,pp%mlps(ik) )
  if( maxval(pp%nproj(:,ik)) > 2 )then
    write(*,*) "stop!!! maxval(pp%nproj)=",maxval(pp%nproj(:,ik))
    goto 999
  end if
  read(4,*)
! non-local projectors
  j0=0
  do ll=0,pp%mlps(ik)
    nproj=pp%nproj(ll,ik)
    read(4,*) dummy_text, ( pp%anorm(j,ik), j=j0,j0+nproj-1 )
    do i=1,pp%mr(ik)+1
      read(4,*) dummy_text, pp%rad(i,ik), ( pp%vpp(i-1,j), j=j0,j0+nproj-1 )
    end do
    j0=j0+nproj
  end do
! local poetntial
  read(4,*) lloc_check
  if( lloc_check == lloc .or. lloc_check > lmaxabinit )then
    j=size(pp%vpp,2)-1
    do i=1,pp%mr(ik)+1
      read(4,*) dummy_text, pp%rad(i,ik), pp%vpp(i-1,j)
    end do
    pp%lref(ik)=j
    write(*,*) "Lref (ik=",ik,") is replaced to",j
  else
    write(*,*) "stop!!! lloc,lloc_check,lmaxabinit =",lloc,lloc_check,lmaxabinit
    goto 999
  end if
! constans for charge integration
  dr=pp%rad(2,ik)-pp%rad(1,ik)
  pi4=4.0d0*acos(-1.0d0)
! core charge for nonlinear core correction and its derivatives
  if( fchrg > 0.0d0 )then
    flag_nlcc_element(ik)=.true.
    do i=1,pp%mr(ik)+1
      read(4,*) dummy_text,r_tmp,rhor_nlcc(i-1,0),rhor_nlcc(i-1,1),rhor_nlcc(i-1,2) !,rhor_nlcc(i-1,3)
    end do
    rhor_nlcc=rhor_nlcc/pi4
    sum_rho_nlcc=0.0d0
    do i=1,pp%mr(ik)+1
      sum_rho_nlcc=sum_rho_nlcc+rhor_nlcc(i-1,0)*pp%rad(i,ik)**2
    end do
    write(*,*) "sum(rho_nlcc)=",sum_rho_nlcc*pi4*dr
  end if
! valence charge density
  sum_rho_pp=0.0d0
  do i=1,pp%mr(ik)+1
    read(4,*) dummy_text, r_tmp, rho_tmp
    sum_rho_pp=sum_rho_pp+rho_tmp*r_tmp**2
  end do
  write(*,*) "sum(rho_pp)=",sum_rho_pp*dr
  close(4)
! extend radial grid data
  do i=pp%mr(ik)+2,pp%nrmax
    pp%rad(i,ik)=dr*(i-1)
  end do
! cut-off radius of the non-local projector
  rrc=0.0d0
  j0=0
  do ll=0,pp%mlps(ik)
    rcpsp=0.0d0
    do j=j0,j0+pp%nproj(ll,ik)-1
      do i=pp%mr(ik)+1,1,-1
        if ( abs(pp%vpp(i-1,j)) > 1.d-12 ) then
          rcpsp = max( rcpsp, pp%rad(i,ik) ) 
          exit
        end if
      end do
    end do
    j0=j
    rrc(ll)=rcpsp
  end do

  !if ( Lmax_ps(ik) >= 0 ) pp%mlps(ik) = Lmax_ps(ik) ! Maximum angular momentum given by input

  return

999 stop "stop@atom/pp/input_pp.f90:read_ps_abinitpsp8"

end subroutine read_ps_abinitpsp8
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine read_ps_fhi(pp,rrc,ik,ps_file)
!This is for original FHI98PP and not for FHI pseudopotential listed in abinit web page
!See http://th.fhi-berlin.mpg.de/th/fhi98md/fhi98PP/
  use structures,only : s_pp_info
  use salmon_global,only : Lmax_ps
  implicit none
  type(s_pp_info),intent(inout) :: pp
!argument
  integer,intent(in) :: ik
  real(8),intent(out) :: rrc(0:pp%lmax0)
  character(256),intent(in) :: ps_file
!local variable
  integer :: i,j
  real(8) :: step,rzps,dummy
  integer :: mr_l(0:pp%lmax0),l,ll
  real(8) :: step_l(0:pp%lmax),rrc_mat(0:pp%lmax,0:pp%lmax)

  rrc_mat(0:pp%lmax,0:pp%lmax)=-1.d0

  open(4,file=ps_file,status='old')
  read(4,*) rzps,pp%mlps(ik)
  pp%zps(ik)=int(rzps+1d-10)
  pp%mlps(ik) = pp%mlps(ik)-1
  if(pp%mlps(ik).gt.pp%lmax0) stop 'Mlps(ik)>Lmax0 at Read_PS_FHI'
  if(pp%mlps(ik).gt.pp%lmax) stop 'Mlps(ik)>Lmax at Read_PS_FHI'
  do i=1,10
     read(4,*) dummy
  end do
  do l=0,pp%mlps(ik)
    read(4,*) mr_l(l),step_l(l)
    if(mr_l(l).gt.pp%nrmax0) stop 'Mr>Nrmax0 at Read_PS_FHI'
    do i=1,mr_l(l)
      read(4,*) j,pp%rad(i+1,ik),pp%upp(i,l),pp%vpp(i,l) !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
    end do
    pp%rad(1,ik)=0.d0
    pp%upp(0,l)=0.d0
    pp%vpp(0,l)=pp%vpp(1,l)-(pp%vpp(2,l)-pp%vpp(1,l))/(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik))
  end do
  close(4)

  if(minval(mr_l(0:pp%mlps(ik))).ne.maxval(mr_l(0:pp%mlps(ik)))) then
    stop 'Mr are diffrent at Read_PS_FHI'
  else 
    pp%mr(ik) = minval(mr_l(0:pp%mlps(ik)))
  end if
  if((maxval(step_l(0:pp%mlps(ik)))-minval(step_l(0:pp%mlps(ik)))).ge.1.d-14) then
    stop 'step are different at Read_PS_FHI'
  else 
    step = minval(step_l(0:pp%mlps(ik)))
  end if

  do i=pp%mr(ik)+1,pp%nrmax
    pp%rad(i+1,ik) = pp%rad(i,ik)*step
  end do

  do l=0,pp%mlps(ik)
    do ll=0,pp%mlps(ik)
      do i=pp%mr(ik),1,-1
        if(abs(pp%vpp(i,l)-pp%vpp(i,ll)).gt.1.d-10) then
          rrc_mat(l,ll) = pp%rad(i+1+1,ik)
          exit
        end if
      end do
    end do
    rrc(l)=maxval(rrc_mat(l,:))
  end do

  if(Lmax_ps(ik) >= 0)pp%mlps(ik) = Lmax_ps(ik) ! Maximum angular momentum given by input

  return
end subroutine read_ps_fhi

subroutine read_ps_adpack(pp,rrc,rhor_nlcc,flag_nlcc_element,ik,ps_file)
  use structures,only : s_pp_info
  use salmon_global,only : nelem
  implicit none
  type(s_pp_info),intent(inout) :: pp
  real(8),intent(out) :: rhor_nlcc(0:pp%nrmax0,0:2)
  logical,intent(inout) :: flag_nlcc_element(nelem)
  integer,intent(in) :: ik
  real(8),intent(out) :: rrc(0:pp%lmax0)
  character(256),intent(in) :: ps_file
  character(50) :: cbuf
  integer :: nprj,iprj,ll,icount,l,l0,i
  real(8) :: rdummy,r,x1,dx
  real(8),allocatable :: x(:)

  flag_beta_proj_is_given = .true.

  open(4,file=ps_file,status="old")

  pp%nproj(:,ik)=0

  do
    read(4,*) cbuf
    if ( cbuf == "valence.electron" ) then
      backspace(4)
      read(4,*) cbuf, pp%zion
      pp%zps(ik)=nint(pp%zion)
    end if
    if ( cbuf == "<project.energies" ) then
      read(4,*) nprj
      do iprj=0,nprj-1
        read(4,*) ll, pp%anorm(iprj,ik), pp%anorm_so(iprj,ik)
        pp%nproj(ll,ik)=pp%nproj(ll,ik)+1
      end do
      pp%mlps(ik)=ll
      exit
    end if
  end do

  pp%lref(ik)=ubound(pp%vpp,2)

  nprj=sum( pp%nproj(:,ik) )

  allocate( x(lbound(pp%rad,1):ubound(pp%rad,1)) ); x=0.0d0

  icount=-1
  do
    read(4,*) cbuf
    if ( cbuf == "<Pseudo.Potentials" ) then
      icount=0
      cycle
    end if
    if ( icount >= 0 ) then
      if ( cbuf == "Pseudo.Potentials>" ) then
        exit
      else
        icount=icount+1
        backspace(4)
        read(4,*) x(icount), pp%rad(icount,ik), pp%vpp(icount-1,pp%lref(ik)), &
             ( pp%vpp(icount-1,iprj), pp%vpp_so(icount-1,iprj), iprj=0,nprj-1 )
      end if
    end if
  end do

  do iprj=0,nprj-1
    do i=icount,1,-1
      pp%vpp(i,iprj) = pp%vpp(i-1,iprj)
      pp%vpp_so(i,iprj) = pp%vpp_so(i-1,iprj)
    end do
  end do
  do i=icount,1,-1
    pp%vpp(i,pp%lref(ik)) = pp%vpp(i-1,pp%lref(ik))
  end do

  do i=icount+1,2,-1
    pp%rad(i,ik) = pp%rad(i-1,ik)
    x(i) = x(i-1)
  end do
  pp%rad(1,ik)=0.0d0
  x(1)=0.0d0

  pp%mr(ik)=icount

  x1=x(2)
  dx=x(pp%mr(ik)+1)-x(pp%mr(ik))
  do i=pp%mr(ik)+2,pp%nrmax
    pp%rad(i,ik) = exp( x1 + dx*(i-1) )
  end do

!-------------------------------------------------------------
!
!  psp%rad(i) == exp( x(i) )
!
!  do i=1,pp%mr(ik)+1
!     write(*,*) i,x(i),pp%rad(i,ik),exp(x(i)),x(i+1)-x(i)
!  end do
!  do i=pp%mr(ik)+2,pp%nrmax
!     write(*,*) i,x1+dx*(i-1),pp%rad(i,ik),exp(x1+dx*(i-1))
!  end do
!-------------------------------------------------------------

  deallocate( x )

  l0=0
  do ll=0,pp%mlps(ik)
    r=0
    do l=l0,l0+pp%nproj(ll,ik)-1
      do i=pp%mr(ik),1,-1
        if ( abs(pp%vpp(i,l)) > 1.d-6 ) then
          r=max(r,pp%rad(i+1,ik))
          exit
        end if
      end do
    end do
    l0=l
    rrc(ll)=r
  end do

  do iprj=0,nprj-1
    do i=0,pp%mr(ik)
      pp%vpp(i,iprj) = pp%rad(i+1,ik)*pp%vpp(i,iprj)
      pp%vpp_so(i,iprj) = pp%rad(i+1,ik)*pp%vpp_so(i,iprj)
    end do
  end do

  icount=-1
  do
    read(4,*,END=9) cbuf
    if ( cbuf == "<density.PCC" ) then
      icount=0
      cycle
    end if
    if ( icount >= 0 ) then
      if ( cbuf == "density.PCC>" ) then
        exit
      else
        icount=icount+1
        backspace(4)
        read(4,*) rdummy, r, rhor_nlcc(icount-1,0)
      end if
    end if
  end do
  if ( icount >= 0 ) flag_nlcc_element(ik)=.true.

  r=0.0d0
  do i=1,pp%mr(ik)
     r=r+rhor_nlcc(i,0)*pp%rad(i+1,ik)**2*(pp%rad(i+1,ik)-pp%rad(i,ik))
  end do
  write(*,*) "r=",r, r*4.0d0*acos(-1.0d0)

9 continue

  close(4)

end subroutine read_ps_adpack

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!    subroutine read_ps_ATOM !.psf format created by ATOM for SIESTA
!      implicit none
!      return
!    end subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine making_ps_with_masking(pp,hx,hy,hz,ik, &
                          rhor_nlcc,flag_nlcc_element)
  use structures,only : s_pp_info
  use salmon_global, only: ps_format, nelem, alpha_mask, eta_mask
  use math_constants, only : pi
  implicit none
  type(s_pp_info),intent(inout) :: pp
  integer,intent(in) :: ik
  real(8),intent(in) :: hx,hy,hz
  real(8),intent(in) :: rhor_nlcc(0:pp%nrmax0,0:2)
  logical,intent(in) :: flag_nlcc_element(nelem)
  real(8) :: eta
  integer :: ncounter
  real(8) :: uvpp(0:pp%nrmax0,0:2*pp%lmax0+1),duvpp(0:pp%nrmax0,0:2*pp%lmax0+1)
  real(8) :: vpploc(0:pp%nrmax0),dvpploc(0:pp%nrmax0)
  real(8) :: grid_function(0:pp%nrmax0)
  integer :: i,l,ll,l0,mr
  real(8) :: r1,r2,r3,r4,const

  ncounter = 0
  do i=0,pp%mr(ik)
    if (pp%rad(i+1,ik) > dble(ncounter+1.d0)*max(Hx,Hy,Hz)) then
      ncounter = ncounter + 1
    end if
    if (ncounter/2*2 == ncounter) then 
      grid_function(i) = 1.d0
    else
      grid_function(i) = 0.d0
    end if
  end do

  if( flag_beta_proj_is_given )then
    mr=pp%mr(ik)
    vpploc(:) = pp%vpp(:,pp%lref(ik))
    do i=1,mr-1
      r1 = pp%rad(i+1,ik)-pp%rad(i,ik)
      r2 = pp%rad(i+1,ik)-pp%rad(i+2,ik)
      r3 = pp%rad(i+2,ik)-pp%rad(i,ik)
      r4 = r1/r2
      dvpploc(i)=(r4+1.d0)*(vpploc(i)-vpploc(i-1))/r1-(vpploc(i+1)-vpploc(i-1))/r3*r4
    end do
    dvpploc(0)=dvpploc(1)-(dvpploc(2)-dvpploc(1))/(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik)-pp%rad(1,ik))
    dvpploc(mr)=dvpploc(mr-1)+(dvpploc(mr-1)-dvpploc(mr-2))/(pp%rad(mr,ik)  &
                             -pp%rad(mr-1,ik))*(pp%rad(mr+1,ik)-pp%rad(mr,ik))
    l0=0
    do ll=0,pp%mlps(ik)
    do l=l0,l0+pp%nproj(ll,ik)-1
      do i=0,mr
        uvpp(i,l)=pp%vpp(i,l)
      end do
      do i=1,mr-1
        r1 = pp%rad(i+1,ik)-pp%rad(i,ik)
        r2 = pp%rad(i+1,ik)-pp%rad(i+2,ik)
        r3 = pp%rad(i+2,ik)-pp%rad(i,ik)
        r4 = r1/r2
        duvpp(i,l)=(r4+1.d0)*(uvpp(i,l)-uvpp(i-1,l))/r1-(uvpp(i+1,l)-uvpp(i-1,l))/r3*r4
        dvpploc(i)=(r4+1.d0)*(vpploc(i)-vpploc(i-1))/r1-(vpploc(i+1)-vpploc(i-1))/r3*r4
      end do
      duvpp(0,l)=2.d0*duvpp(1,l)-duvpp(2,l)
      duvpp(mr,l)=2.d0*duvpp(mr-1,l)-duvpp(mr-2,l)
    end do
    l0=l
    end do
  else
    vpploc(:) = pp%vpp(:,pp%lref(ik))
    do l=0,pp%mlps(ik)
      do i=0,pp%mr(ik)
        uvpp(i,l)=pp%upp(i,l)*(pp%vpp(i,l)-pp%vpp(i,pp%lref(ik)))
      end do
      do i=1,pp%mr(ik)-1
        r1 = pp%rad(i+1,ik)-pp%rad(i,ik)
        r2 = pp%rad(i+1,ik)-pp%rad(i+2,ik)
        r3 = pp%rad(i+2,ik)-pp%rad(i,ik)
        r4 = r1/r2
        duvpp(i,l)=(r4+1.d0)*(uvpp(i,l)-uvpp(i-1,l))/r1-(uvpp(i+1,l)-uvpp(i-1,l))/r3*r4
        dvpploc(i)=(r4+1.d0)*(vpploc(i)-vpploc(i-1))/r1-(vpploc(i+1)-vpploc(i-1))/r3*r4
      end do
      duvpp(0,l)=2.d0*duvpp(1,l)-duvpp(2,l)
      duvpp(pp%mr(ik),l)=2.d0*duvpp(pp%mr(ik)-1,l)-duvpp(pp%mr(ik)-2,l)
      dvpploc(0)=dvpploc(1)-(dvpploc(2)-dvpploc(1))/(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik)-pp%rad(1,ik))
      dvpploc(pp%mr(ik))=dvpploc(pp%mr(ik)-1)+(dvpploc(pp%mr(ik)-1)-dvpploc(pp%mr(ik)-2))/(pp%rad(pp%mr(ik),ik)  &
                        -pp%rad(pp%mr(ik)-1,ik))*(pp%rad(pp%mr(ik)+1,ik)-pp%rad(pp%mr(ik),ik))
    end do
  end if

  open(4,file="PSbeforemask_"//trim(pp%atom_symbol(ik))//"_"//trim(ps_format(ik))//".dat")
  write(4,*) "# Mr =",pp%mr(ik)
  write(4,*) "# Rps(ik), NRps(ik)",pp%rps(ik), pp%nrps(ik)
  write(4,*) "# Mlps(ik), Lref(ik) =",pp%mlps(ik), pp%lref(ik)
  do i=0,pp%mr(ik)
    write(4,'(30e21.12)') pp%rad(i+1,ik),(uvpp(i,l),l=0,pp%mlps(ik)),(duvpp(i,l),l=0,pp%mlps(ik)),  &
                          vpploc(i),dvpploc(i),grid_function(i)
  end do
  close(4)

  call ps_masking(pp,uvpp,duvpp,ik,hx,hy,hz)

  open(4,file="PSaftermask_"//trim(pp%atom_symbol(ik))//"_"//trim(ps_format(ik))//".dat")
  write(4,*) "# Mr =",pp%mr(ik)
  write(4,*) "# Rps(ik), NRps(ik)",pp%rps(ik), pp%nrps(ik)
  write(4,*) "# Mlps(ik), Lref(ik) =",pp%mlps(ik), pp%lref(ik)
  eta = alpha_mask*Pi*pp%rps(ik)/max(Hx,Hy,Hz)
  write(4,*) "# eta_mask, eta =",eta_mask,eta
  do i=0,pp%mr(ik)
    write(4,'(30e21.12)') pp%rad(i+1,ik),(uvpp(i,l),l=0,pp%mlps(ik)),(duvpp(i,l),l=0,pp%mlps(ik)),  &
                          vpploc(i),dvpploc(i),grid_function(i)
  end do
  close(4)

! multiply sqrt((2l+1)/4pi)/r**(l+1) for radial w.f.
  l0=0
  do ll=0,pp%mlps(ik)
    const=sqrt( (2.0d0*ll+1.0d0)/(4.0d0*pi) )
  do l=l0,l0+pp%nproj(ll,ik)-1
    do i=1,pp%mr(ik)
      uvpp(i,l)=uvpp(i,l)*const/(pp%rad(i+1,ik))**(ll+1)
      duvpp(i,l)=duvpp(i,l)*const/(pp%rad(i+1,ik))**(ll+1) &
                -(ll+1)*uvpp(i,l)/pp%rad(i+1,ik)
    end do
    uvpp(0,l)=2.d0*uvpp(1,l)-uvpp(2,l)
    duvpp(0,l)=2.d0*duvpp(1,l)-duvpp(2,l)
  end do
  l0=l
  end do

  if( flag_beta_proj_is_given )then
    do i=1,pp%mr(ik)
      pp%vloctbl(i,ik)=vpploc(i-1)
      pp%dvloctbl(i,ik)=dvpploc(i-1)
    end do
    l0=0
    do ll=0,pp%mlps(ik)
    do l=l0,l0+pp%nproj(ll,ik)-1
      do i=1,pp%mr(ik)+1
        pp%udvtbl(i,l,ik)=uvpp(i-1,l)
        pp%dudvtbl(i,l,ik)=duvpp(i-1,l)
      end do
      if (pp%inorm(l,ik) == 0) cycle
      pp%udvtbl(1:pp%mr(ik)+1,l,ik)=pp%udvtbl(1:pp%mr(ik)+1,l,ik)*pp%anorm(l,ik)
      pp%dudvtbl(1:pp%mr(ik)+1,l,ik)=pp%dudvtbl(1:pp%mr(ik)+1,l,ik)*pp%anorm(l,ik)
    end do
    l0=l  
    end do
  else
    do l=0,pp%mlps(ik)
      do i=1,pp%nrps(ik)
        pp%vloctbl(i,ik)=vpploc(i-1)
        pp%dvloctbl(i,ik)=dvpploc(i-1)
        pp%udvtbl(i,l,ik)=uvpp(i-1,l)
        pp%dudvtbl(i,l,ik)=duvpp(i-1,l)
      enddo
      if (pp%inorm(l,ik) == 0) cycle
      pp%udvtbl(1:pp%nrps(ik),l,ik)=pp%udvtbl(1:pp%nrps(ik),l,ik)/pp%anorm(l,ik)
      pp%dudvtbl(1:pp%nrps(ik),l,ik)=pp%dudvtbl(1:pp%nrps(ik),l,ik)/pp%anorm(l,ik)
    enddo
  end if

  pp%flag_nlcc = pp%flag_nlcc.or.flag_nlcc_element(ik)
  pp%rho_nlcc_tbl(:,ik)=0d0; pp%tau_nlcc_tbl(:,ik)=0d0
  if(.not.flag_nlcc_element(ik))return
!  do i=1,pp%nrps(ik)
  do i=1,pp%mr(ik)
    if(rhor_nlcc(i-1,0)/rhor_nlcc(0,0) < 1d-7)exit
    pp%rho_nlcc_tbl(i,ik)=rhor_nlcc(i-1,0)
    pp%tau_nlcc_tbl(i,ik)=0.25d0*rhor_nlcc(i-1,1)**2/rhor_nlcc(i-1,0)
  end do

  return
end subroutine making_ps_with_masking
!====
subroutine making_ps_with_masking_so(pp,hx,hy,hz,ik)
  use structures,only : s_pp_info
  use math_constants, only : pi
  implicit none
  type(s_pp_info),intent(inout) :: pp
  integer,intent(in) :: ik
  real(8),intent(in) :: hx,hy,hz
  integer :: ncounter
  real(8) :: uvpp(0:pp%nrmax0,0:2*pp%lmax0+1),duvpp(0:pp%nrmax0,0:2*pp%lmax0+1)
  real(8) :: grid_function(0:pp%nrmax0)
  integer :: i,l,ll,l0,mr
  real(8) :: r1,r2,r3,r4

  ncounter = 0
  do i=0,pp%mr(ik)
    if (pp%rad(i+1,ik) > dble(ncounter+1.d0)*max(Hx,Hy,Hz)) then
      ncounter = ncounter + 1
    end if
    if (ncounter/2*2 == ncounter) then 
      grid_function(i) = 1.d0
    else
      grid_function(i) = 0.d0
    end if
  end do

  if( flag_beta_proj_is_given )then
    mr=pp%mr(ik)
    l0=0
    do ll=0,pp%mlps(ik)
    do l=l0,l0+pp%nproj(ll,ik)-1
      do i=0,mr
        uvpp(i,l)=pp%vpp_so(i,l)
      end do
      do i=1,mr-1
        r1 = pp%rad(i+1,ik)-pp%rad(i,ik)
        r2 = pp%rad(i+1,ik)-pp%rad(i+2,ik)
        r3 = pp%rad(i+2,ik)-pp%rad(i,ik)
        r4 = r1/r2
        duvpp(i,l)=(r4+1.d0)*(uvpp(i,l)-uvpp(i-1,l))/r1-(uvpp(i+1,l)-uvpp(i-1,l))/r3*r4
      end do
      duvpp(0,l)=2.d0*duvpp(1,l)-duvpp(2,l)
      duvpp(mr,l)=2.d0*duvpp(mr-1,l)-duvpp(mr-2,l)
    end do
    l0=l
    end do
  else
    write(*,*) "so with upp-given pseudopotential has not been implemented"
    stop "stop@making_ps_with_masking_so"
  end if

!  open(4,file="PSbeforemask_so_"//trim(pp%atom_symbol(ik))//"_"//trim(ps_format(ik))//".dat")
!  write(4,*) "# Mr =",pp%mr(ik)
!  write(4,*) "# Rps(ik), NRps(ik)",pp%rps(ik), pp%nrps(ik)
!  write(4,*) "# Mlps(ik), Lref(ik) =",pp%mlps(ik), pp%lref(ik)
!  do i=0,pp%mr(ik)
!    write(4,'(30e21.12)') pp%rad(i+1,ik),(uvpp(i,l),l=0,pp%mlps(ik)),(duvpp(i,l),l=0,pp%mlps(ik)),  &
!                          vpploc(i),dvpploc(i),grid_function(i)
!  end do
!  close(4)

  call ps_masking(pp,uvpp,duvpp,ik,hx,hy,hz)

!  open(4,file="PSaftermask_so_"//trim(pp%atom_symbol(ik))//"_"//trim(ps_format(ik))//".dat")
!  write(4,*) "# Mr =",pp%mr(ik)
!  write(4,*) "# Rps(ik), NRps(ik)",pp%rps(ik), pp%nrps(ik)
!  write(4,*) "# Mlps(ik), Lref(ik) =",pp%mlps(ik), pp%lref(ik)
!  eta = alpha_mask*Pi*pp%rps(ik)/max(Hx,Hy,Hz)
!  write(4,*) "# eta_mask, eta =",eta_mask,eta
!  do i=0,pp%mr(ik)
!    write(4,'(30e21.12)') pp%rad(i+1,ik),(uvpp(i,l),l=0,pp%mlps(ik)),(duvpp(i,l),l=0,pp%mlps(ik)),  &
!                          vpploc(i),dvpploc(i),grid_function(i)
!  end do
!  close(4)

  if( flag_beta_proj_is_given )then
    l0=0
    do ll=0,pp%mlps(ik)
    do l=l0,l0+pp%nproj(ll,ik)-1
      do i=1,pp%mr(ik)+1
        pp%udvtbl_so(i,l,ik)=uvpp(i-1,l)
        pp%dudvtbl_so(i,l,ik)=duvpp(i-1,l)
      end do
      if (pp%inorm_so(l,ik) == 0) cycle
      pp%udvtbl_so(1:pp%mr(ik)+1,l,ik)=pp%udvtbl_so(1:pp%mr(ik)+1,l,ik)*pp%anorm_so(l,ik)
      pp%dudvtbl_so(1:pp%mr(ik)+1,l,ik)=pp%dudvtbl_so(1:pp%mr(ik)+1,l,ik)*pp%anorm_so(l,ik)
    end do
    l0=l  
    end do
  end if

  return
end subroutine making_ps_with_masking_so
!====
subroutine making_ps_without_masking(pp,ik,flag_nlcc_element,rhor_nlcc)
  use structures,only : s_pp_info
  use salmon_global, only: nelem, method_init_density
  use math_constants, only : pi
  implicit none
  type(s_pp_info),intent(inout) :: pp
  integer,intent(in) :: ik
  logical,intent(in) :: flag_nlcc_element(nelem)
  real(8),intent(in) :: rhor_nlcc(0:pp%nrmax0,0:2)
  integer :: i,l,l0,ll
  real(8) :: r1,r2,r3,r4,const,u
  
  if(method_init_density/='wf' .and. (.not. flag_beta_proj_is_given)) then
    pp%rho_pp_tbl(:,ik) = 0d0
    u = 0d0
    loop_l: do l = 0, pp%mlps(ik)
      loop_i: do i = 1, pp%mr(ik)
        pp%rho_pp_tbl(i,ik) = pp%rho_pp_tbl(i,ik) + dble(2*l+1)* pp%upp(i,l)**2
        u = u + dble(2*l+1)* pp%upp(i,l)**2 * (pp%rad(i+1,ik)-pp%rad(i,ik))
        if( u > pp%zps(ik) ) then
          exit loop_l
          exit loop_i
        end if
      end do loop_i
    end do loop_l
    u = 0d0
    do i = 1, pp%mr(ik)
      u = u + pp%rho_pp_tbl(i,ik)*(pp%rad(i+1,ik)-pp%rad(i,ik))
    end do
    write(*,*) "Int(rho)= ",u, " (for method_init_density=pp...)"
  end if

! multiply sqrt((2l+1)/4pi)/r**(l+1) for radial w.f.
  do l=0,pp%mlps(ik)
    do i=1,pp%mr(ik)
      pp%upp(i,l)=pp%upp(i,l)*sqrt((2*l+1.d0)/(4*pi))/(pp%rad(i+1,ik))**(l+1)
    enddo
    pp%upp(0,l)=pp%upp(1,l)
!    pp%upp(0,l)=2.d0*pp%upp(1,l)-pp%upp(2,l)
  enddo
  
! copy the radial wave functions for DFT+U
  l0=0
  do ll=0,pp%mlps(ik)
  do l=l0,l0+pp%nproj(ll,ik)-1
    do i=0,pp%mr(ik)-1
      pp%upptbl_ao(i+1,l,ik) = pp%upp(i,l)
    end do
  end do
  l0=l
  end do

  l0=0
  do ll=0,pp%mlps(ik)
  do l=l0,l0+pp%nproj(ll,ik)-1
    do i=1,pp%mr(ik)-1
      r1 = pp%rad(i+1,ik)-pp%rad(i,ik)
      r2 = pp%rad(i+1,ik)-pp%rad(i+2,ik)
      r3 = pp%rad(i+2,ik)-pp%rad(i,ik)
      r4 = r1/r2
      pp%dvpp(i,l)=(r4+1.d0)*(pp%vpp(i,l)-pp%vpp(i-1,l))/r1-(pp%vpp(i+1,l)-pp%vpp(i-1,l))/r3*r4
      pp%dupp(i,l)=(r4+1.d0)*(pp%upp(i,l)-pp%upp(i-1,l))/r1-(pp%upp(i+1,l)-pp%upp(i-1,l))/r3*r4
    end do
    pp%dvpp(0,l)=pp%dvpp(1,l)
    pp%dvpp(pp%mr(ik),l)=pp%dvpp(pp%mr(ik)-1,l)
    pp%dupp(0,l)=pp%dupp(1,l)
    pp%dupp(pp%mr(ik),l)=pp%dupp(pp%mr(ik)-1,l)
!    pp%dvpp(0,l)=2.d0*pp%dvpp(1,l)-pp%dvpp(2,l)
!    pp%dvpp(pp%mr(ik),l)=2.d0*pp%dvpp(pp%mr(ik)-1,l)-pp%dvpp(pp%mr(ik)-2,l)
!    pp%dupp(0,l)=2.d0*pp%dupp(1,l)-pp%dupp(2,l)
!    pp%dupp(pp%mr(ik),l)=2.d0*pp%dupp(pp%mr(ik)-1,l)-pp%dupp(pp%mr(ik)-2,l)
  end do
  l0=l
  end do

  if( flag_beta_proj_is_given )then
    l=pp%lref(ik)
    do i=1,pp%mr(ik)-1
      r1 = pp%rad(i+1,ik)-pp%rad(i,ik)
      r2 = pp%rad(i+1,ik)-pp%rad(i+2,ik)
      r3 = pp%rad(i+2,ik)-pp%rad(i,ik)
      r4 = r1/r2
      pp%dvpp(i,l)=(r4+1.d0)*(pp%vpp(i,l)-pp%vpp(i-1,l))/r1-(pp%vpp(i+1,l)-pp%vpp(i-1,l))/r3*r4
    end do
    pp%dvpp(0,l)=pp%dvpp(1,l)
    pp%dvpp(pp%mr(ik),l)=pp%dvpp(pp%mr(ik)-1,l)
  end if

  if( flag_beta_proj_is_given )then
    do i=1,pp%mr(ik)
      pp%vloctbl(i,ik)=pp%vpp(i-1,pp%lref(ik))
      pp%dvloctbl(i,ik)=pp%dvpp(i-1,pp%lref(ik))
    end do
    l0=0
    do ll=0,pp%mlps(ik)
      const=sqrt( (2.0d0*ll+1.0d0)/(4.0d0*pi) )
    do l=l0,l0+pp%nproj(ll,ik)-1
      do i=2,pp%mr(ik)
        pp%udvtbl(i,l,ik)=pp%vpp(i-1,l)/pp%rad(i,ik)**(ll+1)*const
        pp%dudvtbl(i,l,ik)=pp%dvpp(i-1,l)/pp%rad(i,ik)**(ll+1)*const &
                          +pp%vpp(i-1,l)*( -const*(ll+1)/pp%rad(i,ik)**(ll+2) )
      end do
      pp%udvtbl(1,l,ik)=pp%udvtbl(2,l,ik)
      pp%dudvtbl(1,l,ik)=pp%dudvtbl(2,l,ik)
      if (pp%inorm(l,ik) == 0) cycle
      pp%udvtbl(1:pp%mr(ik),l,ik)=pp%udvtbl(1:pp%mr(ik),l,ik)*pp%anorm(l,ik)
      pp%dudvtbl(1:pp%mr(ik),l,ik)=pp%dudvtbl(1:pp%mr(ik),l,ik)*pp%anorm(l,ik)
    end do
    l0=l
    end do
  else
    do l=0,pp%mlps(ik)
      do i=1,pp%nrps(ik)
        pp%vloctbl(i,ik)=pp%vpp(i-1,pp%lref(ik))
        pp%dvloctbl(i,ik)=pp%dvpp(i-1,pp%lref(ik))
        pp%udvtbl(i,l,ik)=(pp%vpp(i-1,l)-pp%vpp(i-1,pp%lref(ik)))*pp%upp(i-1,l)
        pp%dudvtbl(i,l,ik)=(pp%dvpp(i-1,l)-pp%dvpp(i-1,pp%lref(ik)))*pp%upp(i-1,l)   &
                            + (pp%vpp(i-1,l)-pp%vpp(i-1,pp%lref(ik)))*pp%dupp(i-1,l)
      enddo
      if (pp%inorm(l,ik) == 0) cycle
      pp%udvtbl(1:pp%nrps(ik),l,ik)=pp%udvtbl(1:pp%nrps(ik),l,ik)/pp%anorm(l,ik)
      pp%dudvtbl(1:pp%nrps(ik),l,ik)=pp%dudvtbl(1:pp%nrps(ik),l,ik)/pp%anorm(l,ik)
    enddo
  end if

  pp%flag_nlcc = pp%flag_nlcc.or.flag_nlcc_element(ik)
  pp%rho_nlcc_tbl(:,ik)=0d0; pp%tau_nlcc_tbl(:,ik)=0d0

  if ( .not.flag_nlcc_element(ik) ) return

  do i=1,pp%mr(ik)
    if ( rhor_nlcc(i-1,0)/rhor_nlcc(0,0) < 1d-7 ) exit
    pp%rho_nlcc_tbl(i,ik)=rhor_nlcc(i-1,0)
    pp%tau_nlcc_tbl(i,ik)=0.25d0*rhor_nlcc(i-1,1)**2/rhor_nlcc(i-1,0)
  end do

  return
end subroutine making_ps_without_masking

subroutine making_ps_without_masking_so( pp, ik )
  use structures,only : s_pp_info
  use math_constants, only : pi
  implicit none
  type(s_pp_info),intent(inout) :: pp
  integer,intent(in) :: ik
  integer :: i,l,l0,ll
  real(8) :: r1,r2,r3,r4,const

  l0=0
  do ll=0,pp%mlps(ik)
  do l=l0,l0+pp%nproj(ll,ik)-1
    do i=1,pp%mr(ik)-1
      r1 = pp%rad(i+1,ik)-pp%rad(i,ik)
      r2 = pp%rad(i+1,ik)-pp%rad(i+2,ik)
      r3 = pp%rad(i+2,ik)-pp%rad(i,ik)
      r4 = r1/r2
      pp%dvpp_so(i,l)=(r4+1.d0)*(pp%vpp_so(i,l)-pp%vpp_so(i-1,l))/r1-(pp%vpp_so(i+1,l)-pp%vpp_so(i-1,l))/r3*r4
    end do
    pp%dvpp_so(0,l)=pp%dvpp_so(1,l)
    pp%dvpp_so(pp%mr(ik),l)=pp%dvpp_so(pp%mr(ik)-1,l)
  end do
  l0=l
  end do

  if( flag_beta_proj_is_given )then
    l0=0
    do ll=0,pp%mlps(ik)
      const=sqrt( (2.0d0*ll+1.0d0)/(4.0d0*pi) )
    do l=l0,l0+pp%nproj(ll,ik)-1
      do i=2,pp%mr(ik)
        pp%udvtbl_so(i,l,ik) =pp%vpp_so(i-1,l)/pp%rad(i,ik)**(ll+1)*const
        pp%dudvtbl_so(i,l,ik)=pp%dvpp_so(i-1,l)/pp%rad(i,ik)**(ll+1)*const &
                             +pp%vpp_so(i-1,l)*( -const*(ll+1)/pp%rad(i,ik)**(ll+2) )
      end do
      pp%udvtbl_so(1,l,ik) =pp%udvtbl_so(2,l,ik)
      pp%dudvtbl_so(1,l,ik)=pp%dudvtbl_so(2,l,ik)
      if (pp%inorm_so(l,ik) == 0) cycle
      pp%udvtbl_so(1:pp%mr(ik),l,ik) =pp%udvtbl_so(1:pp%mr(ik),l,ik)*pp%anorm_so(l,ik)
      pp%dudvtbl_so(1:pp%mr(ik),l,ik)=pp%dudvtbl_so(1:pp%mr(ik),l,ik)*pp%anorm_so(l,ik)
    end do
    l0=l
    end do
  else
    write(*,*) "so with upp-given pseudopotential has not been implemented"
    stop "stop@making_ps_without_masking_so"
  end if

  return
end subroutine making_ps_without_masking_so

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine ps_masking(pp,uvpp,duvpp,ik,hx,hy,hz)
  use structures,only : s_pp_info
  use salmon_global,only :ps_format,alpha_mask,gamma_mask
  use salmon_math, only: xjl, dxjl
  use math_constants, only : pi
  implicit none
  type(s_pp_info),intent(inout) :: pp
!argument
  integer,intent(in) :: ik
  real(8),intent(inout) :: uvpp(0:pp%nrmax0,0:2*pp%lmax0+1)
  real(8),intent(out) :: duvpp(0:pp%nrmax0,0:2*pp%lmax0+1)
  real(8),intent(in) :: hx,hy,hz
!local variable
!Normalized mask function
  integer,parameter :: nkmax=1000
  integer :: i,j,l,l0,ll
  real(8) :: kmax,k1,k2,kr1,kr2,dr,dk
  real(8),allocatable :: radk(:),wk(:,:) !Fourier staffs
!Mask function
  real(8),allocatable :: mask(:),dmask(:)
!Functions

!Reconstruct radial coordinate Rps(ik) and NRps(ik)
  pp%rps(ik) = gamma_mask*pp%rps(ik)
  do i=1,pp%nrmax0
    if (pp%rad(i,ik) > pp%rps(ik)) exit
  end do
  pp%nrps(ik)=i
  pp%rps(ik) = pp%rad(pp%nrps(ik),ik)
  allocate(mask(pp%nrps(ik)),dmask(pp%nrps(ik)))

  call make_mask_function(pp,mask,dmask,ik)

!Make
  do i = 0,pp%nrps(ik)-1
    l0=0
    do ll = 0,pp%mlps(ik)
    do l = l0,l0+pp%nproj(ll,ik)-1
      uvpp(i,l) = uvpp(i,l)/mask(i+1)
    end do
    l0=l
    end do
  end do

  ll=sum(pp%nproj(:,ik))
  allocate(radk(nkmax),wk(nkmax,0:ll-1))
  wk(:,:)=0.d0
!  kmax = alpha_mask*Pi*sqrt(1.d0/Hx**2+1.d0/Hy**2+1.d0/Hz**2)
  kmax = alpha_mask*pi/max(Hx,Hy,Hz)
  do i = 1,nkmax
    radk(i) = kmax*(dble(i-1)/dble(nkmax-1))
  end do

  do i=1,nkmax
    do j=1,pp%mr(ik)-1
      kr1 = radk(i)*pp%rad(j,ik)
      kr2 = radk(i)*pp%rad(j+1,ik)
      dr = pp%rad(j+1,ik) - pp%rad(j,ik)
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        wk(i,l) = wk(i,l) &
             &+ 0.5d0*(xjl(kr1,ll)*uvpp(j-1,l) + xjl(kr2,ll)*uvpp(j,l))*dr
      end do
      l0=l
      end do
    end do
  end do

  open(4,file="PSFourier_"//trim(pp%atom_symbol(ik))//"_"//trim(ps_format(ik))//".dat")
  write(4,*) "# Kmax, NKmax =",kmax,nkmax
  write(4,*) "# Mlps(ik), Lref(ik) =",pp%mlps(ik), pp%lref(ik)
  write(4,*) "#  Pi/max(Hx,Hy,Hz) =", pi/max(Hx,Hy,Hz)
  write(4,*) "#  Pi*sqrt(1.d0/Hx**2+1.d0/Hy**2+1.d0/Hz**2) =", pi*sqrt(1.d0/Hx**2+1.d0/Hy**2+1.d0/Hz**2)
  ll=sum(pp%nproj(:,ik))
  do i=1,nkmax
    if(radk(i) < (Pi/max(Hx,Hy,Hz))) then
      write(4,'(8e21.12)') radk(i),(wk(i,l),l=0,ll-1),1.d0
    else 
      write(4,'(8e21.12)') radk(i),(wk(i,l),l=0,ll-1),0.d0
    end if
  end do
  close(4)

  uvpp = 0.d0; duvpp=0.d0
  do i=1,nkmax-1
    do j=1,pp%mr(ik)
      kr1 = radk(i)*pp%rad(j,ik)
      kr2 = radk(i+1)*pp%rad(j,ik)
      k1 = radk(i)
      k2 = radk(i+1)
      dk = radk(i+1) - radk(i)
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        uvpp(j-1,l) = uvpp(j-1,l) &
             &+ 0.5d0*(xjl(kr1,ll)*wk(i,l) + xjl(kr2,ll)*wk(i+1,l))*dk
        duvpp(j-1,l) = duvpp(j-1,l) &
             &+ 0.5d0*(k1*dxjl(kr1,ll)*wk(i,l) + k2*dxjl(kr2,ll)*wk(i+1,l))*dk
      end do
      l0=l
      end do
    end do
  end do

  l0=0
  do ll=0,pp%mlps(ik)
  do l=l0,l0+pp%nproj(ll,ik)-1
    uvpp(pp%mr(ik),l) = 2.d0*uvpp(pp%mr(ik)-1,l) - uvpp(pp%mr(ik)-2,l)
    duvpp(pp%mr(ik),l) = 2.d0*duvpp(pp%mr(ik)-1,l) - duvpp(pp%mr(ik)-2,l)
  end do
  l0=l
  end do
  uvpp = (2.d0/pi)*uvpp
  duvpp = (2.d0/pi)*duvpp

  do i=0,pp%nrps(ik)-1
    l0=0
    do ll = 0,pp%mlps(ik)
    do l=l0,l0+pp%nproj(ll,ik)-1
!Derivative calculation before constructing uvpp to avoid overwrite
      duvpp(i,l) = duvpp(i,l)*mask(i+1) + uvpp(i,l)*dmask(i+1)
      uvpp(i,l) = uvpp(i,l)*mask(i+1)
    end do
    l0=l
    end do
  end do

  deallocate(radk,wk)

  return

end subroutine ps_masking
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

subroutine make_mask_function(pp,rmask,dmask,ik)
!Subroutine Make_mask_function
!Name of variables are taken from ***
  use structures,only : s_pp_info
  use salmon_global, only : eta_mask
  use math_constants, only : pi
  implicit none
  type(s_pp_info),intent(inout) :: pp
!Arguments
  integer,intent(in) :: ik
  real(8),intent(inout) :: rmask(pp%nrps(ik)),dmask(pp%nrps(ik))
!local variables
  integer,parameter :: M = 200
!  real(8),parameter :: eta = 15.d0
  integer :: i,j, i3,i2,i1
  real(8) :: xp,xm,dx,nmask0,kx1,dk
  real(8) :: x(M),nmask(M),mat(M,M),k(M),nmask_k(M)
!Lapack dsyev
  integer :: INFO,LWORK
  real(8),allocatable :: WORK(:),W(:)

!Making normalized mask function in radial coordinate
  do i = 1,M
    x(i) = dble(i)/dble(M)
  end do
  do i = 1,M
    xp = 2.d0*x(i)
    mat(i,i) = sin(xp*eta_mask)/xp + dble(M)*pi - eta_mask
    do j = i+1,M
      xp = x(i) + x(j)
      xm = x(i) - x(j)
      mat(i,j) = sin(xp*eta_mask)/xp - sin(xm*eta_mask)/xm
      mat(j,i) = mat(i,j)
    end do
  end do

  allocate(W(M))
  LWORK = max(1,3*M - 1)
  allocate(WORK(LWORK))
  call dsyev('V','U',M,mat,M,W,WORK,LWORK,INFO)
  deallocate(WORK,W)
  nmask0 = 3.d0*mat(1,1)/x(1) - 3.d0*mat(2,1)/x(2) + mat(3,1)/x(3)
  do i = 1,M
    nmask(i) = mat(i,1)/x(i)/nmask0
  end do
  nmask0 = nmask0/nmask0

  open(4,file="nmask.dat")
  write(4,*) "# M =",M
  write(4,*) 0,nmask0
  do i= 1,M
    write(4,*) x(i),nmask(i)
  end do
  close(4)

!Taking Fourier transformation
  do i = 1,M
    k(i)=pi*dble(i)
  end do
  dx = x(2)-x(1)
  nmask_k(:) = 0.d0
  do i = 1,M
    do j = 1,M
      kx1 = k(i)*x(j)
      nmask_k(i) = nmask_k(i) + nmask(j)*kx1*sin(kx1) 
    end do
    nmask_k(i) = nmask_k(i)*dx/k(i)**2
  end do

  open(4,file="nmask_k.dat")
  write(4,*) 0,  3.d0*nmask_k(1) - 3.d0*nmask_k(2) + nmask_k(3)
  do i= 1,M
    write(4,*) k(i),nmask_k(i)
  end do
  close(4)

!  allocate(mask(M),dmask(M))!debug
!Making normalized mask function in radial coordinate
  rmask(:) = 0.d0; dmask(:)=0.d0
  dk = k(2) - k(1)
  do i=2,pp%nrps(ik) !Avoiding divide by zero
    do j = 1,M
      kx1 = k(j)*pp%rad(i,ik)/pp%rps(ik)
      rmask(i) =  rmask(i) + nmask_k(j)*kx1*sin(kx1)
      dmask(i)= dmask(i) + nmask_k(j)*(kx1**2*cos(kx1)-kx1*sin(kx1))
    end do
    rmask(i) = (2.d0/pi)* rmask(i)*dk*pp%rps(ik)**2/pp%rad(i,ik)**2
    dmask(i)= (2.d0/pi)*dmask(i)*dk*pp%rps(ik)**2/pp%rad(i,ik)**3 
  end do
  rmask(1) =  rmask(2)-( rmask(3)- rmask(2))/(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik)-pp%rad(1,ik))
  dmask(1)= dmask(2)-(dmask(3)-dmask(2))/(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik)-pp%rad(1,ik))
  i1=pp%nrps(ik)-2
  i2=pp%nrps(ik)-1
  i3=pp%nrps(ik)
   rmask(i3)=  rmask(i2)+( rmask(i2)- rmask(i1))/(pp%rad(i2,ik)-pp%rad(i1,ik))*(pp%rad(i3,ik)-pp%rad(i2,ik))
  dmask(i3)= dmask(i2)+(dmask(i2)-dmask(i1))/(pp%rad(i2,ik)-pp%rad(i1,ik))*(pp%rad(i3,ik)-pp%rad(i2,ik))

  open(4,file="mask.dat")
  write(4,*) "# rps(ik), nrps(ik) =",pp%rps(ik), pp%nrps(ik)
  do i= 1,pp%nrps(ik)
    write(4,'(8e22.10)') pp%rad(i,ik),rmask(i),dmask(i)
  end do
  close(4)

  return
end subroutine make_mask_function
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

end module input_pp_sub
