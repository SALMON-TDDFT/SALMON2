module lcfo_wannier_sawf_orchestrator
  implicit none
  private

  type, public :: t_sawf_environment_receipt
    integer :: environment=0
    integer :: representative_fragment=0
    integer :: operation_index=0
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

contains

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
