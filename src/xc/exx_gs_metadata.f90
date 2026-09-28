! Versioned physics contract for conventional global-hybrid GS -> native RT.
! The standard info.bin, occupation.bin and wavefunction formats stay unchanged.
module exx_gs_metadata
  implicit none
  private
  public :: exx_gs_metadata_write,exx_gs_metadata_check,exx_gs_occupation_check
contains
  function physics_values(system,grid) result(values)
    use structures, only: s_dft_system
    use salmon_global, only: xc,num_kgrid,pbeh_coulomb_radius,rvv10_b,rvv10_c,rvv10_nq, &
      alpha_mask,gamma_mask,eta_mask
    use exx_functional, only: exchange_fraction,exchange_screening
    implicit none
    type(s_dft_system),intent(in) :: system
    integer,intent(in) :: grid(3)
    real(8),allocatable :: values(:)
    real(8) :: radius,vdw(3)
    radius=pbeh_coulomb_radius
    if(radius<=0d0)radius=.5d0*minval(grid*num_kgrid*system%hgs)
    vdw=0d0
    if(xc=='pbeh40_rvv10')vdw=[rvv10_b,rvv10_c,real(rvv10_nq,8)]
    values=[exchange_fraction(),exchange_screening(),radius,vdw,system%hgs, &
      reshape(system%primitive_a,[9]),reshape(system%vec_k,[3*system%nk]),system%wtk, &
      reshape(system%Rion,[3*system%nion]),alpha_mask,gamma_mask,eta_mask]
  end function

  subroutine exx_gs_metadata_write(path,lg,system,info)
    use structures, only: s_rgrid,s_dft_system,s_parallel_info
    use communication, only: comm_is_root,comm_bcast
    use salmon_global, only: xc,nelem,num_kgrid,lmax_ps,lloc_ps,izatom,yn_psmask,exx_pre_scf_active
    use exx_functional, only: is_global_hybrid
    implicit none
    character(*),intent(in) :: path
    type(s_rgrid),intent(in) :: lg
    type(s_dft_system),intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    integer :: unit,status,ios
    logical :: exists
    status=0
    if(comm_is_root(info%id_rko))then
      if(is_global_hybrid(xc).and..not.exx_pre_scf_active)then
        open(newunit=unit,file=trim(path)//'hybrid_gs.bin',status='replace', &
          form='unformatted',action='write',iostat=status)
        if(status==0)then
          write(unit,iostat=status)1,[system%nk,system%no,system%nspin,system%nion,nelem]
          if(status==0)write(unit,iostat=status)xc
          if(status==0)write(unit,iostat=status)lg%num,num_kgrid,lg%Nd
          if(status==0)write(unit,iostat=status)physics_values(system,lg%num)
          if(status==0)write(unit,iostat=status)system%rocc
          if(status==0)write(unit,iostat=status)system%if_real_orbital
          if(status==0)write(unit,iostat=status)system%kion,lmax_ps(:nelem),lloc_ps(:nelem),izatom(:nelem),yn_psmask
          if(status==0)call pseudopotential_records(unit,.true.,status)
          close(unit,iostat=ios)
          if(ios/=0)status=ios
        endif
      else
        ! PBE/HSE output or unfinished PBE warmup must not retain a hybrid certificate.
        inquire(file=trim(path)//'hybrid_gs.bin',exist=exists)
        if(exists)then
          open(newunit=unit,file=trim(path)//'hybrid_gs.bin',status='old',iostat=status)
          if(status==0)close(unit,status='delete',iostat=status)
        endif
      endif
    endif
    call comm_bcast(status,info%icomm_rko)
    if(status/=0)error stop 'Conventional hybrid GS metadata write failed'
  end subroutine

  subroutine exx_gs_metadata_check(path,lg,system,info,occupation)
    use structures, only: s_rgrid,s_dft_system,s_parallel_info
    use communication, only: comm_is_root,comm_bcast
    use salmon_global, only: xc,nelem,num_kgrid,lmax_ps,lloc_ps,izatom,yn_psmask
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    character(*),intent(in) :: path
    type(s_rgrid),intent(in) :: lg
    type(s_dft_system),intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    real(8),allocatable,intent(out) :: occupation(:,:,:)
    real(8),allocatable :: current(:),saved(:)
    integer,allocatable :: species(:),lmax(:),lloc(:),atomic_number(:)
    integer :: unit,status,ios,version,dims(5),grid(3),kgrid(3),nd
    character(64) :: saved_xc
    character(1) :: mask
    logical :: saved_real
    allocate(occupation(system%no,system%nk,system%nspin))
    status=0
    if(comm_is_root(info%id_rko))then
      open(newunit=unit,file=trim(path)//'hybrid_gs.bin',status='old', &
        form='unformatted',action='read',iostat=status)
      if(status==0)then
        read(unit,iostat=status)version,dims
        if(status==0)then
          if(version/=1.or.any(dims/=[system%nk,system%no,system%nspin,system%nion,nelem]))status=1
        endif
        if(status==0)read(unit,iostat=status)saved_xc
        if(status==0)then
          if(saved_xc/=xc)status=1
        endif
        if(status==0)read(unit,iostat=status)grid,kgrid,nd
        if(status==0)then
          if(any(grid/=lg%num).or.any(kgrid/=num_kgrid).or.nd/=lg%Nd)status=1
        endif
        if(status==0)then
          current=physics_values(system,lg%num)
          allocate(saved(size(current)))
          read(unit,iostat=status)saved
          if(status==0)then
            if(.not.all(ieee_is_finite(saved)).or..not.all(ieee_is_finite(current)))then
              status=1
            else if(any(abs(saved-current)>1d-12*max(1d0,abs(current))))then
              status=1
            endif
          endif
        endif
        if(status==0)read(unit,iostat=status)occupation
        if(status==0)then
          ! This route evolves a complete, doubly occupied, one-spin subspace only.
          if(system%nspin/=1.or..not.all(ieee_is_finite(occupation)))then
            status=1
          else if(any(abs(occupation-2d0)>1d-12))then
            status=1
          else if(.not.all(ieee_is_finite(system%rocc)))then
            status=1
          else if(any(abs(occupation-system%rocc)>1d-12))then
            status=1
          endif
        endif
        if(status==0)read(unit,iostat=status)saved_real
        if(status==0)then
          allocate(species(system%nion),lmax(nelem),lloc(nelem),atomic_number(nelem))
          read(unit,iostat=status)species,lmax,lloc,atomic_number,mask
          if(status==0)then
            if(any(species/=system%kion).or.any(lmax/=lmax_ps(:nelem)).or.any(lloc/=lloc_ps(:nelem)).or. &
              any(atomic_number/=izatom(:nelem)).or.mask/=yn_psmask)status=1
          endif
        endif
        if(status==0)call pseudopotential_records(unit,.false.,status)
        close(unit,iostat=ios)
        if(ios/=0)status=ios
      endif
    endif
    call comm_bcast(status,info%icomm_rko)
    if(status/=0)error stop 'Conventional hybrid GS metadata missing, malformed, or mismatched'
    ! Check dimensions before the legacy reader allocates/indexes occupation and wavefunction arrays.
    if(comm_is_root(info%id_rko))call checkpoint_payload_check(path,system,occupation,saved_real,status)
    call comm_bcast(status,info%icomm_rko)
    if(status/=0)error stop 'Conventional hybrid GS occupation payload mismatch or invalid info.bin'
    call comm_bcast(occupation,info%icomm_rko)
  end subroutine

  subroutine checkpoint_payload_check(path,system,occupation,saved_real,status)
    use structures, only: s_dft_system
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    character(*),intent(in) :: path
    type(s_dft_system),intent(in) :: system
    real(8),intent(in) :: occupation(:,:,:)
    logical,intent(in) :: saved_real
    integer,intent(out) :: status
    integer :: unit,ios,nk,no,iteration,ranks
    logical :: real_orbitals
    real(8),allocatable :: payload(:,:,:)
    open(newunit=unit,file=trim(path)//'info.bin',status='old',form='unformatted',action='read',iostat=status)
    if(status/=0)return
    read(unit,iostat=status)nk
    if(status==0)read(unit,iostat=status)no
    if(status==0)read(unit,iostat=status)iteration
    if(status==0)read(unit,iostat=status)ranks
    if(status==0)read(unit,iostat=status)real_orbitals
    close(unit,iostat=ios)
    if(ios/=0)status=ios
    if(status/=0)return
    if(nk/=system%nk.or.no/=system%no.or.(real_orbitals.neqv.saved_real))then
      status=1;return
    endif
    allocate(payload(system%no,system%nk,system%nspin))
    open(newunit=unit,file=trim(path)//'occupation.bin',status='old',form='unformatted',action='read',iostat=status)
    if(status/=0)return
    read(unit,iostat=status)payload
    close(unit,iostat=ios)
    if(ios/=0)status=ios
    if(status/=0)return
    if(.not.all(ieee_is_finite(payload)))then
      status=1
    else if(any(abs(payload-occupation)>1d-12))then
      status=1
    endif
  end subroutine

  subroutine exx_gs_occupation_check(system,info,occupation)
    use structures, only: s_dft_system,s_parallel_info
    use communication, only: comm_get_max
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    real(8),intent(in) :: occupation(:,:,:)
    integer :: status
    status=0
    if(any(shape(system%rocc)/=shape(occupation)))then
      status=1
    else if(.not.all(ieee_is_finite(system%rocc)))then
      status=1
    else if(any(abs(system%rocc-occupation)>1d-12))then
      status=1
    endif
    call comm_get_max(status,info%icomm_rko)
    if(status/=0)error stop 'Conventional hybrid GS occupation payload mismatch'
  end subroutine

  subroutine pseudopotential_records(unit,writing,status)
    use salmon_global, only: nelem,file_pseudo
    implicit none
    integer,intent(in) :: unit
    logical,intent(in) :: writing
    integer,intent(out) :: status
    integer :: j,pu,nbytes,saved_bytes,ios
    character(:),allocatable :: bytes,saved_text
    status=0
    do j=1,nelem
      inquire(file=trim(file_pseudo(j)),size=nbytes,iostat=status)
      if(status/=0)return
      if(nbytes<1)then
        status=1;return
      endif
      allocate(character(nbytes)::bytes)
      open(newunit=pu,file=trim(file_pseudo(j)),form='unformatted',access='stream',status='old', &
        action='read',iostat=status)
      if(status==0)then
        read(pu,iostat=status)bytes
        close(pu,iostat=ios)
        if(ios/=0)status=ios
      endif
      if(status==0)then
        if(writing)then
          write(unit,iostat=status)nbytes
          if(status==0)write(unit,iostat=status)bytes
        else
          read(unit,iostat=status)saved_bytes
          if(status==0)then
            if(saved_bytes/=nbytes)status=1
          endif
          if(status==0)then
            allocate(character(nbytes)::saved_text)
            read(unit,iostat=status)saved_text
            if(status==0)then
              if(saved_text/=bytes)status=1
            endif
            deallocate(saved_text)
          endif
        endif
      endif
      deallocate(bytes)
      if(status/=0)return
    enddo
  end subroutine
end module
