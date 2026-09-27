! Explicit nuclear derivative at fixed DC fragment orbitals/occupations.
! This diagnostic excludes fragment response and is NOT an MD force.
module dc_force
  use structures
  implicit none
  private
  public :: report_dc_frozen_force
contains
  subroutine report_dc_frozen_force(dc,system,info,mg,pp,ppg,psi,ewald,energy)
    use communication, only: comm_summation
    use force_sub, only: force_ewald_rspace,differentiate_projectors
    use salmon_global, only: kion,cutoff_g,temperature
    use math_constants, only: zi
    use dc_projector_force, only: core_projector_force
    type(s_dcdft),intent(in) :: dc
    type(s_dft_system),intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_rgrid),intent(in) :: mg
    type(s_pp_info),intent(in) :: pp
    type(s_pp_grid),intent(in) :: ppg
    type(s_orbital),intent(in) :: psi
    type(s_ewald_ion_ion),intent(in) :: ewald
    type(s_dft_energy),intent(in) :: energy
    real(8),allocatable :: local(:,:),electrostatic(:,:),realspace(:,:),fragment(:,:),assembled(:,:),gradient(:,:,:)
    integer,allocatable :: saved_kion(:)
    integer :: na,ia,atom,ix,iy,iz,ilocal,ilma,ik,io,spin,j
    real(8) :: g(3),g2,weight,ionic_factor,f,core_norm,occupation_local(2),occupation_frag(2),occupation_total(2)
    integer :: lower(3),upper(3)
    complex(8) :: phase,vg
    complex(8),allocatable :: values(:),derivatives(:,:)
    logical,allocatable :: in_core(:)
    real(8) :: projector_force(3)
    if(.not.allocated(dc%atom_global))error stop 'DC force diagnostic: atom map missing'
    if(.not.allocated(psi%zwf).or.info%if_divide_rspace) &
      error stop 'DC force diagnostic: complex unsplit fragment orbitals required'
    na=dc%system_tot%nion
    allocate(local(3,na),electrostatic(3,na),realspace(3,na),fragment(3,na),assembled(3,na))
    ! The inherited Ewald helper reads global kion. Restore the fragment view.
    call move_alloc(kion,saved_kion)
    allocate(kion(na));kion=dc%system_tot%kion
    call force_ewald_rspace(realspace,local,dc%system_tot,dc%info_tot,ewald,pp,na,dc%info_tot%icomm_r)
    deallocate(kion);call move_alloc(saved_kion,kion)
    local=0d0
    do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
    do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
    do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
      if(dc%fg_tot%if_Gzero(ix,iy,iz))cycle
      g=dc%fg_tot%vec_G(:,ix,iy,iz);g2=sum(g**2)
      if(g2>cutoff_g**2)cycle
      do ia=dc%info_tot%ia_s,dc%info_tot%ia_e
        atom=dc%system_tot%kion(ia)
        phase=exp(zi*sum(g*dc%system_tot%rion(:,ia)))
        ionic_factor=pp%zps(atom)*dc%fg_tot%coef(ix,iy,iz)*dc%fg_tot%exp_ewald(ix,iy,iz)
        vg=dc%ppg_tot%zVG_ion(ix,iy,iz,atom)-dc%fg_tot%coef(ix,iy,iz)*pp%zps(atom)
        local(:,ia)=local(:,ia)+g*(ionic_factor*aimag(dc%ppg_tot%zrhoG_ion(ix,iy,iz)*phase) &
          +aimag(phase*dc%poisson_tot%zrhoG_ele(ix,iy,iz)*conjg(vg)))
      enddo
    enddo;enddo;enddo
    call comm_summation(local,electrostatic,3*na,dc%info_tot%icomm_rko)
    electrostatic=electrostatic+realspace

    allocate(gradient(3,ppg%nps,ppg%nlma))
    call differentiate_projectors(pp,ppg,system%kion,gradient)
    allocate(values(ppg%nps),derivatives(3,ppg%nps),in_core(ppg%nps))
    local=0d0
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do spin=1,system%nspin
      do ilocal=1,ppg%ilocal_nlma
        ilma=ppg%ilocal_nlma2ilma(ilocal);ia=ppg%ilocal_nlma2ia(ilocal)
        do j=1,ppg%mps(ia)
          ix=ppg%jxyz(1,j,ia);iy=ppg%jxyz(2,j,ia);iz=ppg%jxyz(3,j,ia)
          values(j)=conjg(ppg%zekr_uV(j,ilma,ik))*psi%zwf(ix,iy,iz,spin,io,ik,1)
          ! Common atomic-center Bloch phase derivatives cancel in core*full.
          phase=exp(zi*sum((system%vec_k(:,ik)+system%vec_Ac)*ppg%rxyz(:,j,ia)))
          derivatives(:,j)=-gradient(:,j,ilma)*phase*psi%zwf(ix,iy,iz,spin,io,ik,1)
          in_core(j)=all([ix,iy,iz]>=mg%is).and.all([ix,iy,iz]<=min(mg%ie,dc%nxyz_domain))
        enddo
        weight=system%rocc(io,ik,spin)*system%wtk(ik)*system%hvol*ppg%rinv_uvu(ilma)
        atom=dc%atom_global(ia)
        call core_projector_force(values(:ppg%mps(ia)),derivatives(:,:ppg%mps(ia)), &
          in_core(:ppg%mps(ia)),weight,projector_force)
        local(:,atom)=local(:,atom)+projector_force
      enddo
    enddo;enddo;enddo
    call comm_summation(local,fragment,3*na,info%icomm_rko)
    local=0d0
    if(info%id_rko==0)local=fragment
    call comm_summation(local,assembled,3*na,dc%icomm_tot)
    assembled=assembled+electrostatic
    ! Core-weighted occupation entropy is a diagnostic of the thermal term.
    ! E-TS is not asserted stationary for truncated DC fragments.
    occupation_local=0d0;lower=mg%is;upper=min(mg%ie,dc%nxyz_domain)
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do spin=1,system%nspin
      core_norm=sum(abs(psi%zwf(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3),spin,io,ik,1))**2) &
        *system%hvol*system%wtk(ik)
      f=system%rocc(io,ik,spin)/2d0
      occupation_local(1)=occupation_local(1)+core_norm*system%rocc(io,ik,spin)
      if(f>0d0.and.f<1d0)occupation_local(2)=occupation_local(2) &
        -2d0*core_norm*(f*log(f)+(1d0-f)*log(1d0-f))*max(temperature,0d0)
    enddo;enddo;enddo
    call comm_summation(occupation_local,occupation_frag,2,info%icomm_rko)
    occupation_local=0d0
    if(info%id_rko==0)occupation_local=occupation_frag
    call comm_summation(occupation_local,occupation_total,2,dc%icomm_tot)
    if(dc%id_tot==0)then
      write(*,'(a)')'DC frozen-orbital force diagnostic (Ha/bohr)'
      do ia=1,na
        write(*,'(i8,3es25.16)')ia,assembled(:,ia)
      enddo
      write(*,'(a)')'DC force diagnostic excludes fragment response; not certified for MD'
      write(*,'(a,es25.16)')'DC occupation electron residual: ',occupation_total(1)-dc%elec_num_tot
      write(*,'(a,es25.16)')'DC occupation TS diagnostic (Ha): ',occupation_total(2)
      write(*,'(a,es25.16)')'DC E-minus-TS diagnostic (Ha): ',energy%E_tot-occupation_total(2)
    endif
  end subroutine
end module
