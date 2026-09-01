#include "config.h"
program test_dg_dc_seed_checkpoint_mpi
  use iso_fortran_env,only:int64
  use mpi
  use dg_dc_seed_checkpoint
  implicit none
  integer::comm,rank,nproc,ierr,status
  integer(int64)::publication_id,read_publication_id,committed_publication_id
  integer(int64)::first_publication_id,second_publication_id
  character(256)::directory,scenario,message
  character(512)::first_directory,second_directory,rank_directory
  type(s_dg_dc_seed_contract)::contract,expected
  type(s_dg_dc_seed_payload)::payload,restored,committed_payload
  logical::ok,exists_a,exists_b
  real(8),parameter::expected_electrons=2d0,electron_tolerance=1d-12
  real(8),parameter::density_weight=0.5d0
  real(8),parameter::write_threshold=1d-8

  call MPI_Init(ierr)
  comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr)
  call MPI_Comm_size(comm,nproc,ierr)
  call get_environment_variable('DG_DC_SEED_DIRECTORY',directory)
  call get_environment_variable('DG_DC_SEED_SCENARIO',scenario)
  call require(len_trim(directory)>0,'missing DG_DC_SEED_DIRECTORY')
  call require(len_trim(scenario)>0,'missing DG_DC_SEED_SCENARIO')
  call initialize_contract(contract)
  call initialize_payload(payload)
  publication_id=0_int64

  select case(trim(scenario))
  case('absent')
    call probe_dg_dc_seed(comm,trim(directory),contract,density_weight,expected_electrons,electron_tolerance,&
      write_threshold,status,publication_id,message)
    call require(status==DG_DC_SEED_ABSENT,trim(message))
    call read_dg_dc_seed(comm,trim(directory),contract,density_weight,expected_electrons,electron_tolerance,&
      write_threshold,restored,read_publication_id,ok,message)
    call require(.not.ok.and.read_publication_id==0_int64,'absent seed was read')

  case('roundtrip')
    call publish_valid_seed()
    call probe_valid(contract,write_threshold)
    call read_valid(contract,write_threshold)
    call require(same_payload(payload,restored),'rank-local seed payload changed in round trip')
    call require(read_publication_id==publication_id,'publication ID changed in round trip')
    ! These downstream controls are intentionally absent from the immutable seed contract.
    call downstream_only_change(3.5d0,71,0.25d0)
    call probe_valid(contract,write_threshold)
    call read_valid(contract,write_threshold)
    call require(same_payload(payload,restored),'downstream-only controls invalidated the seed')

  case('missing_manifest')
    call publish_valid_seed()
    if(rank==0)call remove_file(manifest_name(directory))
    call MPI_Barrier(comm,ierr)
    call require_invalid(contract,write_threshold,'missing manifest was accepted')

  case('missing_shard')
    call publish_valid_seed()
    if(rank==0)call remove_file(shard_name(directory,publication_id,nproc-1))
    call MPI_Barrier(comm,ierr)
    call require_invalid(contract,write_threshold,'missing rank shard was accepted')

  case('truncated_shard')
    call publish_valid_seed()
    if(rank==0)call truncate_file(shard_name(directory,publication_id,nproc-1))
    call MPI_Barrier(comm,ierr)
    call require_invalid(contract,write_threshold,'truncated rank shard was accepted')

  case('corrupt_manifest')
    call publish_valid_seed()
    if(rank==0)call corrupt_payload_byte(manifest_name(directory))
    call MPI_Barrier(comm,ierr)
    call require_invalid(contract,write_threshold,'manifest payload corruption was accepted')

  case('corrupt_shard')
    call publish_valid_seed()
    if(rank==0)call corrupt_payload_byte(shard_name(directory,publication_id,nproc-1))
    call MPI_Barrier(comm,ierr)
    call require_invalid(contract,write_threshold,'rank shard payload corruption was accepted')

  case('interrupted_publication')
    call write_dg_dc_seed(comm,trim(directory),contract,payload,density_weight,expected_electrons,&
      electron_tolerance,write_threshold,publication_id,ok,message,failure_injection_rank=0)
    call require(.not.ok.and.publication_id/=0_int64,'interrupted publication unexpectedly committed')
    call require_invalid(contract,write_threshold,'interrupted publication was accepted')

  case('interrupted_update_preserves_committed')
    call publish_valid_seed()
    committed_publication_id=publication_id
    committed_payload=payload
    call alter_payload(payload)
    call write_dg_dc_seed(comm,trim(directory),contract,payload,density_weight,expected_electrons,&
      electron_tolerance,write_threshold,publication_id,ok,message,failure_injection_rank=0)
    call require(.not.ok,'interrupted update unexpectedly committed')
    call probe_dg_dc_seed(comm,trim(directory),contract,density_weight,expected_electrons,&
      electron_tolerance,write_threshold,status,read_publication_id,message)
    call require(status==DG_DC_SEED_VALID,'interrupted update hid the committed seed: '//trim(message))
    call require(read_publication_id==committed_publication_id,&
      'interrupted update replaced the committed publication ID')
    call read_dg_dc_seed(comm,trim(directory),contract,density_weight,expected_electrons,&
      electron_tolerance,write_threshold,restored,read_publication_id,ok,message)
    call require(ok,'committed seed was unreadable after interrupted update: '//trim(message))
    call require(read_publication_id==committed_publication_id,&
      'read returned the interrupted update publication ID')
    call require(same_payload(committed_payload,restored),&
      'interrupted update changed the committed payload')

  case('mixed_publication_read')
    first_directory=trim(directory)//'/a'
    second_directory=trim(directory)//'/b'
    committed_payload=payload
    call publish_seed_at(first_directory,first_publication_id)
    call require_valid_at(first_directory,first_publication_id)
    if(nproc==1)then
      call read_seed_at(first_directory,read_publication_id,ok)
      call require(ok.and.read_publication_id==first_publication_id,&
        'single-rank publication was not readable')
      call require(same_payload(committed_payload,restored),&
        'single-rank publication payload changed')
      publication_id=first_publication_id
    else
      call alter_payload(payload)
      call publish_seed_at(second_directory,second_publication_id)
      call require_valid_at(second_directory,second_publication_id)
      if(mod(rank,2)==0)then
        rank_directory=first_directory
      else
        rank_directory=second_directory
      endif
      call probe_dg_dc_seed(comm,trim(rank_directory),contract,density_weight,expected_electrons,&
        electron_tolerance,write_threshold,status,read_publication_id,message)
      call require(status==DG_DC_SEED_INVALID,&
        'rank-mixed directories/publications were accepted by probe')
      call read_seed_at(rank_directory,read_publication_id,ok)
      call require(.not.ok,'rank-mixed directories/publications were accepted by read')
      publication_id=second_publication_id
    endif

  case('rank_directory_write_mismatch')
    first_directory=trim(directory)//'/a'
    second_directory=trim(directory)//'/b'
    if(mod(rank,2)==0)then
      rank_directory=first_directory
    else
      rank_directory=second_directory
    endif
    call write_dg_dc_seed(comm,trim(rank_directory),contract,payload,density_weight,expected_electrons,&
      electron_tolerance,write_threshold,publication_id,ok,message)
    if(nproc==1)then
      call require(ok,'single-rank directory write failed: '//trim(message))
      call require_valid_at(first_directory,publication_id)
    else
      call require(.not.ok,'rank-disagreeing directory write unexpectedly committed')
      inquire(file=trim(manifest_name(first_directory)),exist=exists_a)
      inquire(file=trim(manifest_name(second_directory)),exist=exists_b)
      call require(.not.exists_a.and..not.exists_b,&
        'rank-disagreeing directory write produced a committed manifest')
    endif

  case('rank_count_mismatch')
    call publish_valid_seed()
    expected=contract;expected%mpi_size=nproc+1
    call require_invalid(expected,write_threshold,'MPI rank-count mismatch was accepted')

  case('write_only')
    call publish_valid_seed()

  case('read_existing_rank_mismatch')
    call require_invalid(contract,write_threshold,'persisted MPI rank-count mismatch was accepted')

  case('rank_fragment_mismatch')
    call publish_valid_seed()
    expected=contract
    if(rank==0)expected%fragment_id=expected%fragment_id+1
    call require_invalid(expected,write_threshold,'rank-to-fragment mismatch was accepted')

  case('local_bound_mismatch')
    call publish_valid_seed()
    expected=contract
    if(rank==nproc-1)expected%rho_bounds(2)=expected%rho_bounds(2)+1
    call require_invalid(expected,write_threshold,'local array-bound mismatch was accepted')

  case('ownership_map_mismatch')
    call publish_valid_seed()
    expected=contract
    if(rank==0)expected%ownership_fingerprint=expected%ownership_fingerprint+1_int64
    call require_invalid(expected,write_threshold,'ownership-map mismatch was accepted')

  case('immutable_mismatch')
    call publish_valid_seed()
    expected=contract
    if(rank==nproc-1)expected%immutable_fingerprint=expected%immutable_fingerprint+1_int64
    call require_invalid(expected,write_threshold,'immutable-input mismatch was accepted')

  case('electron_count_mismatch')
    call publish_valid_seed()
    call probe_dg_dc_seed(comm,trim(directory),contract,density_weight,expected_electrons+0.5d0,&
      electron_tolerance,write_threshold,status,read_publication_id,message)
    call require(status==DG_DC_SEED_INVALID,'electron-count mismatch was accepted')

  case('residual_threshold_mismatch')
    call publish_valid_seed()
    call require_invalid(contract,1d-12,'current convergence threshold mismatch was accepted')

  case default
    call require(.false.,'unknown DG DC seed test scenario')
  end select

  if(rank==0)write(*,'(a,a,a,i0,a,i0)')'PASS DG DC seed scenario=',trim(scenario),&
    ' ranks=',nproc,' publication_id=',publication_id
  call MPI_Finalize(ierr)

contains
  subroutine initialize_contract(value)
    type(s_dg_dc_seed_contract),intent(out)::value
    value%version=1
    value%mpi_size=nproc
    value%rank=rank
    value%fragment_id=rank+1
    value%rwf_bounds=[rank,1,0,1,1,1,1, rank+1,2,0,1,2,1,1]
    value%rho_bounds=[rank,0,-1, rank+1,1,-1]
    value%vloc_bounds=[rank,0,-1, rank+1,1,-1]
    value%immutable_fingerprint=100000_int64+int(nproc,int64)
    value%ownership_fingerprint=200000_int64+int(rank,int64)
  end subroutine

  subroutine initialize_payload(value)
    type(s_dg_dc_seed_payload),intent(out)::value
    integer::i
    allocate(value%rwf(&
      contract%rwf_bounds(1):contract%rwf_bounds(8),&
      contract%rwf_bounds(2):contract%rwf_bounds(9),&
      contract%rwf_bounds(3):contract%rwf_bounds(10),&
      contract%rwf_bounds(4):contract%rwf_bounds(11),&
      contract%rwf_bounds(5):contract%rwf_bounds(12),&
      contract%rwf_bounds(6):contract%rwf_bounds(13),&
      contract%rwf_bounds(7):contract%rwf_bounds(14)))
    allocate(value%rho_tot(contract%rho_bounds(1):contract%rho_bounds(4),&
      contract%rho_bounds(2):contract%rho_bounds(5),&
      contract%rho_bounds(3):contract%rho_bounds(6)))
    allocate(value%vloc_tot(contract%vloc_bounds(1):contract%vloc_bounds(4),&
      contract%vloc_bounds(2):contract%vloc_bounds(5),&
      contract%vloc_bounds(3):contract%vloc_bounds(6)))
    allocate(value%esp(-1:1,1:1,1:1),value%rocc(-1:1,1:1,1:1))
    value%rwf=reshape([(real(1000*rank+i,8)/100d0,i=1,size(value%rwf))],shape(value%rwf))
    value%rho_tot=1d0/real(nproc,8)
    value%vloc_tot=reshape([(real(100*rank+i,8)/10d0,i=1,size(value%vloc_tot))],shape(value%vloc_tot))
    value%esp(:,1,1)=[-0.75d0,-0.25d0,0.5d0]
    value%rocc(:,1,1)=[2d0,0d0,0d0]
    value%mu=-0.5d0
    value%residual=1d-10
    value%iteration=17
  end subroutine

  subroutine publish_valid_seed()
    call write_dg_dc_seed(comm,trim(directory),contract,payload,density_weight,expected_electrons,&
      electron_tolerance,write_threshold,publication_id,ok,message)
    call require(ok,trim(message))
    call require(publication_id/=0_int64,'writer returned a zero publication ID')
  end subroutine

  subroutine publish_seed_at(path,id)
    character(*),intent(in)::path
    integer(int64),intent(out)::id
    call write_dg_dc_seed(comm,trim(path),contract,payload,density_weight,expected_electrons,&
      electron_tolerance,write_threshold,id,ok,message)
    call require(ok,'cannot publish independent seed: '//trim(message))
    call require(id/=0_int64,'independent seed writer returned a zero publication ID')
  end subroutine

  subroutine require_valid_at(path,id)
    character(*),intent(in)::path
    integer(int64),intent(in)::id
    call probe_dg_dc_seed(comm,trim(path),contract,density_weight,expected_electrons,&
      electron_tolerance,write_threshold,status,read_publication_id,message)
    call require(status==DG_DC_SEED_VALID,'independent seed is invalid: '//trim(message))
    call require(read_publication_id==id,'independent seed publication ID changed')
  end subroutine

  subroutine read_seed_at(path,id,read_ok)
    character(*),intent(in)::path
    integer(int64),intent(out)::id
    logical,intent(out)::read_ok
    call read_dg_dc_seed(comm,trim(path),contract,density_weight,expected_electrons,&
      electron_tolerance,write_threshold,restored,id,read_ok,message)
  end subroutine

  subroutine alter_payload(value)
    type(s_dg_dc_seed_payload),intent(inout)::value
    value%rwf=value%rwf+10d0
    value%vloc_tot=value%vloc_tot-3d0
    value%esp=value%esp+0.125d0
    value%mu=value%mu+0.0625d0
    value%iteration=value%iteration+1
  end subroutine

  subroutine probe_valid(value,threshold)
    type(s_dg_dc_seed_contract),intent(in)::value
    real(8),intent(in)::threshold
    call probe_dg_dc_seed(comm,trim(directory),value,density_weight,expected_electrons,electron_tolerance,&
      threshold,status,read_publication_id,message)
    call require(status==DG_DC_SEED_VALID,trim(message))
    call require(read_publication_id==publication_id,'probe returned another publication')
  end subroutine

  subroutine read_valid(value,threshold)
    type(s_dg_dc_seed_contract),intent(in)::value
    real(8),intent(in)::threshold
    call read_dg_dc_seed(comm,trim(directory),value,density_weight,expected_electrons,electron_tolerance,&
      threshold,restored,read_publication_id,ok,message)
    call require(ok,trim(message))
  end subroutine

  subroutine require_invalid(value,threshold,text)
    type(s_dg_dc_seed_contract),intent(in)::value
    real(8),intent(in)::threshold
    character(*),intent(in)::text
    call probe_dg_dc_seed(comm,trim(directory),value,density_weight,expected_electrons,electron_tolerance,&
      threshold,status,read_publication_id,message)
    call require(status==DG_DC_SEED_INVALID,text//': '//trim(message))
    call read_dg_dc_seed(comm,trim(directory),value,density_weight,expected_electrons,electron_tolerance,&
      threshold,restored,read_publication_id,ok,message)
    call require(.not.ok,text)
  end subroutine

  logical function same_payload(left,right)
    type(s_dg_dc_seed_payload),intent(in)::left,right
    same_payload=allocated(right%rwf).and.allocated(right%rho_tot).and.&
      allocated(right%vloc_tot).and.allocated(right%esp).and.allocated(right%rocc)
    if(.not.same_payload)return
    same_payload=all(lbound(left%rwf)==lbound(right%rwf)).and.&
      all(ubound(left%rwf)==ubound(right%rwf)).and.&
      all(lbound(left%rho_tot)==lbound(right%rho_tot)).and.&
      all(ubound(left%rho_tot)==ubound(right%rho_tot)).and.&
      all(lbound(left%vloc_tot)==lbound(right%vloc_tot)).and.&
      all(ubound(left%vloc_tot)==ubound(right%vloc_tot)).and.&
      all(lbound(left%esp)==lbound(right%esp)).and.&
      all(ubound(left%esp)==ubound(right%esp)).and.&
      all(lbound(left%rocc)==lbound(right%rocc)).and.&
      all(ubound(left%rocc)==ubound(right%rocc)).and.&
      all(left%rwf==right%rwf).and.all(left%rho_tot==right%rho_tot).and.&
      all(left%vloc_tot==right%vloc_tot).and.all(left%esp==right%esp).and.&
      all(left%rocc==right%rocc).and.left%mu==right%mu.and.&
      left%residual==right%residual.and.left%iteration==right%iteration
  end function

  subroutine downstream_only_change(localization_cutoff,wannier_iterations,symmetry_window)
    real(8),intent(in)::localization_cutoff,symmetry_window
    integer,intent(in)::wannier_iterations
    call require(localization_cutoff>0d0.and.wannier_iterations>0.and.symmetry_window>=0d0,&
      'invalid downstream-only test control')
  end subroutine

  function manifest_name(path)result(filename)
    character(*),intent(in)::path
    character(512)::filename
    filename=trim(path)//'/dg_dc_seed.manifest'
  end function

  function shard_name(path,id,shard_rank)result(filename)
    character(*),intent(in)::path
    integer(int64),intent(in)::id
    integer,intent(in)::shard_rank
    character(512)::filename
    write(filename,'(a,"/dg_dc_seed.",z16.16,".rank",i8.8,".shard")')&
      trim(path),id,shard_rank
  end function

  subroutine remove_file(filename)
    character(*),intent(in)::filename
    integer::unit,ios
    open(newunit=unit,file=trim(filename),status='old',iostat=ios)
    if(ios==0)close(unit,status='delete')
  end subroutine

  subroutine truncate_file(filename)
    character(*),intent(in)::filename
    integer::unit
    open(newunit=unit,file=trim(filename),status='replace',access='stream',form='unformatted')
    write(unit)'TRUNCATED'
    close(unit)
  end subroutine

  subroutine corrupt_payload_byte(filename)
    character(*),intent(in)::filename
    integer::unit,ios
    integer(int64)::file_size
    inquire(file=trim(filename),size=file_size,iostat=ios)
    call require(ios==0.and.file_size>0_int64,'cannot size file selected for payload corruption')
    open(newunit=unit,file=trim(filename),status='old',access='stream',form='unformatted',action='readwrite')
    write(unit,pos=file_size)'X'
    close(unit)
  end subroutine

  subroutine require(condition,text)
    logical,intent(in)::condition
    character(*),intent(in)::text
    if(.not.condition)then
      write(0,'(a)')trim(text)
      error stop 1
    endif
  end subroutine
end program
