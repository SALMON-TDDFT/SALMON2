#include "config.h"
program test_dg_fragment_wf_checkpoint_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int8,int64,real64
  use dg_fragment_wf_checkpoint
  implicit none
  integer::comm,rank,nproc,ierr,status
  integer(int64)::publication_id,read_publication_id,committed_publication_id
  character(512)::directory,message
  character(64)::scenario,requested_mode
  type(s_dg_fragment_wf_contract)::contract,expected
  type(s_dg_fragment_wf_payload)::payload,restored
  logical::ok,reuse,regenerate,publish,fatal

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call get_environment_variable('DG_FRAGMENT_WF_DIRECTORY',directory)
  call get_environment_variable('DG_FRAGMENT_WF_SCENARIO',scenario)
  call require(len_trim(directory)>0.and.len_trim(scenario)>0,'missing checkpoint test environment')
  call initialize_contract(contract);call initialize_payload(payload)

  select case(trim(scenario))
  case('roundtrip')
    call publish_valid()
    call read_dg_fragment_wf_checkpoint(comm,trim(directory),contract,restored,&
      read_publication_id,ok,message)
    call require(ok,trim(message));call require(read_publication_id==publication_id,&
      'publication ID changed in round trip')
    call require(same_payload(payload,restored),'fragment-WF payload changed in round trip')
  case('rank_count_mismatch')
    call publish_valid();expected=contract;expected%mpi_size=nproc+1
    call require_invalid(expected,'mpi_rank_count')
  case('rank_fragment_permutation')
    call publish_valid();expected=contract
    if(nproc>1)then
      expected%fragment_id=modulo(rank+1,nproc)+1
    else
      expected%mapping_fingerprint=expected%mapping_fingerprint+1
    endif
    call require_invalid(expected,'rank_fragment_mapping')
  case('dc_seed_mismatch')
    call publish_valid();expected=contract;expected%dc_seed_publication_id=expected%dc_seed_publication_id+1
    call require_invalid(expected,'dc_seed')
  case('dc_seed_fingerprint_mismatch')
    call publish_valid();expected=contract;expected%dc_seed_fingerprint=expected%dc_seed_fingerprint+1
    call require_invalid(expected,'dc_seed')
  case('grid_mismatch')
    call publish_valid();expected=contract;expected%grid_fingerprint=expected%grid_fingerprint+1
    call require_invalid(expected,'grid')
  case('cell_mismatch')
    call publish_valid();expected=contract;expected%cell_fingerprint=expected%cell_fingerprint+1
    call require_invalid(expected,'cell')
  case('fragment_geometry_mismatch')
    call publish_valid();expected=contract
    expected%fragment_geometry_fingerprint=expected%fragment_geometry_fingerprint+1
    call require_invalid(expected,'fragment_geometry')
  case('pseudopotential_mismatch')
    call publish_valid();expected=contract
    expected%pseudopotential_fingerprint=expected%pseudopotential_fingerprint+1
    call require_invalid(expected,'pseudopotential')
  case('boundary_mismatch')
    call publish_valid();expected=contract;expected%boundary_fingerprint=expected%boundary_fingerprint+1
    call require_invalid(expected,'boundary')
  case('inventory_mismatch')
    call publish_valid();expected=contract;expected%inventory_fingerprint=expected%inventory_fingerprint+1
    call require_invalid(expected,'inventory')
  case('ordering_mismatch')
    call publish_valid();expected=contract;expected%ordering_fingerprint=expected%ordering_fingerprint+1
    call require_invalid(expected,'inventory')
  case('selection_mismatch')
    call publish_valid();expected=contract;expected%selection_fingerprint=expected%selection_fingerprint+1
    call require_invalid(expected,'selection')
  case('local_layout_mismatch')
    call publish_valid();expected=contract
    expected%local_layout_fingerprint=expected%local_layout_fingerprint+1
    call require_invalid(expected,'grid')
  case('generation_mismatch')
    call publish_valid();expected=contract;expected%basis_generation=expected%basis_generation+1
    call require_invalid(expected,'basis_generation')
  case('gauge_version_mismatch')
    call publish_valid();expected=contract;expected%gauge_algorithm_version=expected%gauge_algorithm_version+1
    call require_invalid(expected,'gauge')
  case('gauge_mode_mismatch')
    call publish_valid();expected=contract;expected%gauge_mode='spectral'
    call require_invalid(expected,'gauge')
  case('gauge_fingerprint_mismatch')
    call publish_valid();expected=contract;expected%gauge_fingerprint=expected%gauge_fingerprint+1
    call require_invalid(expected,'gauge')
  case('missing_manifest')
    call publish_valid();if(rank==0)call remove_file(manifest_name(directory));call MPI_Barrier(comm,ierr)
    call probe_dg_fragment_wf_checkpoint(comm,trim(directory),contract,status,read_publication_id,message)
    call require(status==DG_FRAGMENT_WF_ABSENT,'incomplete publication was not an auto miss')
  case('missing_peer')
    call publish_valid();if(rank==0)call remove_file(shard_name(directory,publication_id,nproc-1))
    call MPI_Barrier(comm,ierr);call require_invalid(contract,'missing_payload')
  case('truncated_payload')
    call publish_valid();if(rank==0)call truncate_file(shard_name(directory,publication_id,nproc-1))
    call MPI_Barrier(comm,ierr);call require_invalid(contract,'payload')
  case('corrupt_payload')
    call publish_valid();if(rank==0)call corrupt_last_byte(shard_name(directory,publication_id,nproc-1))
    call MPI_Barrier(comm,ierr);call require_invalid(contract,'payload_hash')
  case('unknown_version')
    call publish_valid();if(rank==0)call overwrite_manifest_version(manifest_name(directory),99)
    call MPI_Barrier(comm,ierr);call require_invalid(contract,'version')
  case('corrupt_manifest')
    call publish_valid();if(rank==0)call corrupt_last_byte(manifest_name(directory))
    call MPI_Barrier(comm,ierr);call require_invalid(contract,'manifest')
  case('interrupted_publication')
    call write_dg_fragment_wf_checkpoint(comm,trim(directory),contract,payload,publication_id,&
      ok,message,failure_injection_rank=0)
    call require(.not.ok,'injected publication unexpectedly committed')
    call probe_dg_fragment_wf_checkpoint(comm,trim(directory),contract,status,read_publication_id,message)
    call require(status==DG_FRAGMENT_WF_ABSENT,'partial publication became visible')
  case('interrupted_update_preserves_committed')
    call publish_valid();committed_publication_id=publication_id
    payload%wannier_values=payload%wannier_values+cmplx(0.25d0,0d0,real64)
    call write_dg_fragment_wf_checkpoint(comm,trim(directory),contract,payload,publication_id,&
      ok,message,failure_injection_rank=0)
    call require(.not.ok,'injected update unexpectedly committed')
    call initialize_payload(payload)
    call read_dg_fragment_wf_checkpoint(comm,trim(directory),contract,restored,&
      read_publication_id,ok,message)
    call require(ok.and.read_publication_id==committed_publication_id,&
      'interrupted update replaced the committed generation')
    call require(same_payload(payload,restored),'interrupted update changed the committed payload')
  case('mode_disagreement')
    if(mod(rank,2)==0)then
      requested_mode='auto'
    else
      requested_mode='read'
    endif
    call decide_dg_fragment_wf_restart(comm,trim(requested_mode),DG_FRAGMENT_WF_ABSENT,reuse,&
      regenerate,publish,fatal,ok,message)
    if(nproc==1)then
      call require(ok.and.regenerate.and.publish,'single-rank auto policy failed')
    else
      call require(.not.ok.and.fatal.and.index(message,'disagree')>0,'rank-disagreeing mode was accepted')
    endif
  case('auto_invalid')
    call decide_dg_fragment_wf_restart(comm,'auto',DG_FRAGMENT_WF_INVALID,reuse,&
      regenerate,publish,fatal,ok,message)
    call require(ok.and..not.reuse.and.regenerate.and.publish.and..not.fatal,&
      'auto mode did not collectively regenerate an invalid generation')
  case('status_disagreement')
    status=merge(DG_FRAGMENT_WF_VALID,DG_FRAGMENT_WF_INVALID,mod(rank,2)==0)
    call decide_dg_fragment_wf_restart(comm,'auto',status,reuse,regenerate,publish,fatal,ok,message)
    if(nproc==1)then
      call require(ok.and.reuse,'single-rank valid status failed')
    else
      call require(.not.ok.and.fatal.and.index(message,'disagree')>0,'rank-disagreeing status was accepted')
    endif
  case('write_only')
    call publish_valid()
  case('read_existing_rank_mismatch')
    call require_invalid(contract,'mpi_rank_count')
  case default
    call require(.false.,'unknown fragment-WF checkpoint scenario')
  end select
  if(rank==0)write(*,'(3a,i0)')'PASS fragment-WF checkpoint scenario=',trim(scenario),' ranks=',nproc
  call MPI_Finalize(ierr)
contains
  subroutine initialize_contract(value)
    type(s_dg_fragment_wf_contract),intent(out)::value
    value%version=1;value%mpi_size=nproc;value%rank=rank;value%fragment_id=rank+1
    value%basis_generation=7;value%gauge_algorithm_version=3
    value%candidate_rank=3;value%retained_rank=2;value%local_row_count=3;value%seed_count=2
    value%gauge_mode='scdm'
    value%mapping_fingerprint=101_int64;value%dc_seed_publication_id=102_int64
    value%dc_seed_fingerprint=103_int64;value%grid_fingerprint=104_int64
    value%cell_fingerprint=105_int64;value%pseudopotential_fingerprint=106_int64
    value%fragment_geometry_fingerprint=107_int64;value%boundary_fingerprint=108_int64
    value%inventory_fingerprint=109_int64;value%ordering_fingerprint=110_int64
    value%selection_fingerprint=111_int64;value%gauge_fingerprint=112_int64
    value%local_layout_fingerprint=1000_int64+rank
  end subroutine initialize_contract
  subroutine initialize_payload(value)
    type(s_dg_fragment_wf_payload),intent(out)::value
    integer::i,j
    allocate(value%local_grid_ids(3),value%selected_state_ids(2),value%wannier_values(2,3),&
      value%candidate_compression(3,2),value%wannier_transform(2,2),value%centers_fractional(3,2),&
      value%dc_seed_coefficients(2,2),value%dc_seed_energies(2),value%dc_seed_occupations(2))
    value%local_grid_ids=[1_int64,2_int64,3_int64]+3_int64*rank
    value%selected_state_ids=[1_int64,3_int64];value%seed_reconstruction_defect=1d-13
    do j=1,3;do i=1,2
      value%wannier_values(i,j)=cmplx(rank+0.1d0*i,0.01d0*j,real64)
    enddo;enddo
    value%candidate_compression=reshape([(cmplx(0.1d0*i,0.02d0*i,real64),i=1,6)],[3,2])
    value%wannier_transform=reshape([cmplx(1d0,0d0,real64),cmplx(0d0,0d0,real64),&
      cmplx(0d0,0d0,real64),cmplx(1d0,0d0,real64)],[2,2])
    value%centers_fractional=reshape([(0.05d0*i+0.01d0*rank,i=1,6)],[3,2])
    value%dc_seed_coefficients=value%wannier_transform
    value%dc_seed_energies=[-0.5d0,-0.2d0];value%dc_seed_occupations=[2d0,1d0]
  end subroutine initialize_payload
  subroutine publish_valid()
    call write_dg_fragment_wf_checkpoint(comm,trim(directory),contract,payload,publication_id,ok,message)
    call require(ok,'write failed: '//trim(message));call require(publication_id/=0_int64,'zero publication ID')
  end subroutine publish_valid
  subroutine require_invalid(wanted,cause)
    type(s_dg_fragment_wf_contract),intent(in)::wanted
    character(*),intent(in)::cause
    call probe_dg_fragment_wf_checkpoint(comm,trim(directory),wanted,status,read_publication_id,message)
    call require(status==DG_FRAGMENT_WF_INVALID,'invalid checkpoint accepted: '//trim(message))
    call require(index(message,trim(cause))>0,'missing rejection cause '//trim(cause)//': '//trim(message))
    call read_dg_fragment_wf_checkpoint(comm,trim(directory),wanted,restored,read_publication_id,ok,message)
    call require(.not.ok.and..not.allocated(restored%wannier_values),'strict read partially restored invalid data')
  end subroutine require_invalid
  logical function same_payload(a,b)
    type(s_dg_fragment_wf_payload),intent(in)::a,b
    same_payload=allocated(b%wannier_values).and.allocated(b%local_grid_ids).and.&
      all(a%local_grid_ids==b%local_grid_ids).and.all(a%selected_state_ids==b%selected_state_ids).and.&
      all(a%wannier_values==b%wannier_values).and.all(a%candidate_compression==b%candidate_compression).and.&
      all(a%wannier_transform==b%wannier_transform).and.all(a%centers_fractional==b%centers_fractional).and.&
      all(a%dc_seed_coefficients==b%dc_seed_coefficients).and.&
      all(a%dc_seed_energies==b%dc_seed_energies).and.&
      all(a%dc_seed_occupations==b%dc_seed_occupations).and.&
      a%seed_reconstruction_defect==b%seed_reconstruction_defect
  end function same_payload
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    if(.not.condition)then;write(*,'(2a)')'FAIL: ',trim(label);call MPI_Abort(comm,1,ierr);endif
  end subroutine require
  character(512) function manifest_name(base)
    character(*),intent(in)::base;manifest_name=trim(base)//'/dg_fragment_wf.manifest'
  end function manifest_name
  character(512) function shard_name(base,id,owner)
    character(*),intent(in)::base;integer(int64),intent(in)::id;integer,intent(in)::owner
    write(shard_name,'(a,"/dg_fragment_wf.publication-",z16.16,".rank-",i6.6,".bin")')trim(base),id,owner
  end function shard_name
  subroutine remove_file(path)
    character(*),intent(in)::path;integer::u,ios
    open(newunit=u,file=trim(path),status='old',iostat=ios);if(ios==0)close(u,status='delete')
  end subroutine remove_file
  subroutine truncate_file(path)
    character(*),intent(in)::path;integer::u
    open(newunit=u,file=trim(path),status='replace',access='stream',form='unformatted');write(u)'TRUNCATED';close(u)
  end subroutine truncate_file
  subroutine corrupt_last_byte(path)
    character(*),intent(in)::path;integer::u,ios;integer(int64)::n;integer(int8)::byte
    open(newunit=u,file=trim(path),status='old',access='stream',form='unformatted',action='readwrite')
    inquire(unit=u,size=n);read(u,pos=n,iostat=ios)byte;if(ios==0)write(u,pos=n)ieor(byte,1_int8);close(u)
  end subroutine corrupt_last_byte
  subroutine overwrite_manifest_version(path,new_version)
    character(*),intent(in)::path;integer,intent(in)::new_version;integer::u
    open(newunit=u,file=trim(path),status='old',access='stream',form='unformatted',action='readwrite')
    write(u,pos=33)new_version;close(u)
  end subroutine overwrite_manifest_version
end program test_dg_fragment_wf_checkpoint_mpi
