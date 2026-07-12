module lcfo_wannier_sawf_orchestrator
  implicit none
  private

  type, public :: t_sawf_environment_receipt
    integer :: environment=0
    integer :: representative_fragment=0
    integer :: operation_index=0
    integer :: num_bands=0
    integer :: num_wann=0
    logical :: generated_independently=.false.
    logical :: requires_execution=.false.
    logical :: completed=.false.
    character(256) :: same_supercell_fingerprint=''
  end type t_sawf_environment_receipt

  type, public :: t_sawf_seed_bundle
    integer :: environment=0
    character(512) :: directory=''
    character(256) :: seedname=''
    character(256) :: same_supercell_fingerprint=''
  end type t_sawf_seed_bundle

  public :: build_sawf_environment_execution_plan
  public :: validate_sawf_environment_receipts
  public :: build_sawf_seed_bundles
  public :: complete_sawf_seed_bundle
  public :: propagate_sawf_representative_receipts

contains

  subroutine propagate_sawf_representative_receipts(receipts,ok,message)
    type(t_sawf_environment_receipt),intent(inout)::receipts(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::environment,representative

    ok=.false.;message=''
    do environment=1,size(receipts)
      if(receipts(environment)%requires_execution)then
        if(.not.receipts(environment)%completed.or.receipts(environment)%num_bands<=0.or. &
            receipts(environment)%num_wann<=0)then
          message='SAWF representative receipt is incomplete';return
        end if
        cycle
      end if
      representative=receipts(environment)%representative_fragment
      if(representative<1.or.representative>size(receipts).or. &
          .not.receipts(representative)%requires_execution.or. &
          .not.receipts(representative)%completed.or.receipts(environment)%operation_index<=0.or. &
          trim(receipts(environment)%same_supercell_fingerprint)/= &
          trim(receipts(representative)%same_supercell_fingerprint))then
        message='SAWF replica receipt has invalid representative provenance';return
      end if
      receipts(environment)%num_bands=receipts(representative)%num_bands
      receipts(environment)%num_wann=receipts(representative)%num_wann
      receipts(environment)%completed=.true.
    end do
    ok=.true.
  end subroutine propagate_sawf_representative_receipts

  subroutine complete_sawf_seed_bundle(bundle,receipt,ok,message)
    type(t_sawf_seed_bundle),intent(in)::bundle
    type(t_sawf_environment_receipt),intent(inout)::receipt
    logical,intent(out)::ok
    character(*),intent(out)::message
    character(16),parameter::suffix(6)=[character(16)::'win','eig','amn','mmn','dmn','chk']
    character(1024)::filename
    character(256)::stored_fingerprint,header
    integer::item,unit,ios,num_bands,num_kpoints,num_wann,file_size,iband,ikpoint,entry
    integer::mmn_bands,mmn_kpoints,mmn_neighbors
    real(8)::energy
    logical::exists

    ok=.false.;message='';receipt%completed=.false.;receipt%num_bands=0;receipt%num_wann=0
    if(.not.receipt%requires_execution.or.bundle%environment/=receipt%environment.or. &
        trim(bundle%same_supercell_fingerprint)/=trim(receipt%same_supercell_fingerprint))then
      message='SAWF seed bundle and execution receipt provenance disagree';return
    end if
    do item=1,size(suffix)
      filename=trim(bundle%directory)//'/'//trim(bundle%seedname)//'.'//trim(suffix(item))
      inquire(file=trim(filename),exist=exists,size=file_size)
      if(.not.exists.or.file_size<=0)then
        message='SAWF seed bundle artifact is missing: '//trim(filename);return
      end if
    end do
    filename=trim(bundle%directory)//'/'//trim(bundle%seedname)//'.amn'
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
    if(ios/=0)then;message='SAWF local .amn header is missing';return;end if
    read(unit,'(a)',iostat=ios)header
    if(ios==0)read(unit,*,iostat=ios)num_bands,num_kpoints,num_wann
    close(unit)
    if(ios/=0.or.num_bands<=0.or.num_kpoints/=1.or.num_wann<=0.or.num_wann>num_bands)then
      message='SAWF local .amn dimensions are invalid';return
    end if
    filename=trim(bundle%directory)//'/'//trim(bundle%seedname)//'.eig'
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
    if(ios/=0)then;message='SAWF local .eig is unreadable';return;end if
    do entry=1,num_bands
      read(unit,*,iostat=ios)iband,ikpoint,energy
      if(ios/=0.or.iband/=entry.or.ikpoint/=1)exit
    end do
    if(ios==0)read(unit,*,iostat=ios)
    close(unit)
    if(entry<=num_bands.or.ios==0)then;message='SAWF local .eig dimensions are invalid';return;end if
    filename=trim(bundle%directory)//'/'//trim(bundle%seedname)//'.mmn'
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
    if(ios/=0)then;message='SAWF local .mmn is unreadable';return;end if
    read(unit,'(a)',iostat=ios)header
    if(ios==0)read(unit,*,iostat=ios)mmn_bands,mmn_kpoints,mmn_neighbors
    close(unit)
    if(ios/=0.or.mmn_bands/=num_bands.or.mmn_kpoints/=1.or.mmn_neighbors<=0)then
      message='SAWF local .mmn dimensions are invalid';return
    end if
    filename=trim(bundle%directory)//'/'//trim(bundle%seedname)//'.sawf-fingerprint'
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
    if(ios/=0)then;message='SAWF seed bundle fingerprint is missing';return;end if
    read(unit,'(a)',iostat=ios)stored_fingerprint;close(unit)
    if(ios/=0.or.trim(stored_fingerprint)/=trim(receipt%same_supercell_fingerprint))then
      message='SAWF seed bundle belongs to a different supercell fingerprint';return
    end if
    receipt%num_bands=num_bands;receipt%num_wann=num_wann
    receipt%completed=.true.;ok=.true.
  end subroutine complete_sawf_seed_bundle

  subroutine build_sawf_seed_bundles(receipts,root_directory,base_seedname,bundles,ok,message)
    type(t_sawf_environment_receipt),intent(in)::receipts(:)
    character(*),intent(in)::root_directory,base_seedname
    type(t_sawf_seed_bundle),allocatable,intent(out)::bundles(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::environment,bundle_index,nrun
    character(16)::suffix

    ok=.false.;message=''
    nrun=count(receipts%requires_execution)
    if(size(receipts)<=0.or.nrun<=0.or.len_trim(root_directory)==0.or.len_trim(base_seedname)==0)then
      message='SAWF representative seed-bundle inputs are empty';return
    end if
    allocate(bundles(nrun));bundle_index=0
    do environment=1,size(receipts)
      if(.not.receipts(environment)%requires_execution)cycle
      if(receipts(environment)%environment/=environment.or. &
          len_trim(receipts(environment)%same_supercell_fingerprint)==0)then
        message='SAWF representative seed-bundle receipt provenance is invalid';return
      end if
      bundle_index=bundle_index+1;write(suffix,'(i6.6)')environment
      bundles(bundle_index)%environment=environment
      bundles(bundle_index)%directory=trim(root_directory)//'/environment-'//trim(suffix)
      bundles(bundle_index)%seedname=trim(base_seedname)//'-env-'//trim(suffix)
      bundles(bundle_index)%same_supercell_fingerprint= &
        trim(receipts(environment)%same_supercell_fingerprint)
    end do
    ok=.true.
  end subroutine build_sawf_seed_bundles

  subroutine build_sawf_environment_execution_plan(representative_fragment,operation_index, &
      generated_independently,supercell_fingerprint,receipts,ok,message)
    integer,intent(in)::representative_fragment(:),operation_index(:)
    logical,intent(in)::generated_independently(:)
    character(*),intent(in)::supercell_fingerprint
    type(t_sawf_environment_receipt),allocatable,intent(out)::receipts(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::environment

    ok=.false.;message=''
    if(size(representative_fragment)<=0.or.size(operation_index)/=size(representative_fragment).or. &
        size(generated_independently)/=size(representative_fragment).or. &
        len_trim(supercell_fingerprint)==0)then
      message='SAWF environment execution-plan dimensions or fingerprint are invalid';return
    end if
    if(any(representative_fragment<1).or.any(operation_index<0))then
      message='SAWF environment execution-plan provenance is invalid';return
    end if
    allocate(receipts(size(representative_fragment)))
    do environment=1,size(receipts)
      receipts(environment)%environment=environment
      receipts(environment)%representative_fragment=representative_fragment(environment)
      receipts(environment)%operation_index=operation_index(environment)
      receipts(environment)%generated_independently=generated_independently(environment)
      receipts(environment)%requires_execution=generated_independently(environment)
      receipts(environment)%same_supercell_fingerprint=trim(supercell_fingerprint)
    end do
    ok=.true.
  end subroutine build_sawf_environment_execution_plan

  subroutine validate_sawf_environment_receipts(receipts,supercell_fingerprint,ok,message)
    type(t_sawf_environment_receipt),intent(in)::receipts(:)
    character(*),intent(in)::supercell_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::environment

    ok=.false.;message=''
    if(size(receipts)<=0.or.len_trim(supercell_fingerprint)==0)then
      message='SAWF environment receipts are empty';return
    end if
    do environment=1,size(receipts)
      if(receipts(environment)%environment/=environment.or. &
          (receipts(environment)%requires_execution.and..not.receipts(environment)%completed).or. &
          trim(receipts(environment)%same_supercell_fingerprint)/=trim(supercell_fingerprint))then
        message='SAWF environment receipt is incomplete or belongs to another supercell';return
      end if
    end do
    ok=.true.
  end subroutine validate_sawf_environment_receipts

end module lcfo_wannier_sawf_orchestrator
