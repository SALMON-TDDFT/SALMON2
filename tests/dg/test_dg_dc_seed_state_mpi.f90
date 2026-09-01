#include "config.h"
program test_dg_dc_seed_state_mpi
  use iso_fortran_env,only:int64
  use mpi
  use dg_dc_seed_checkpoint
  implicit none
  integer::comm,rank,nproc,ierr,status,iteration
  integer(int64)::publication_id,immutable_inputs(6),ownership_map(8)
  integer::rwf_bounds(14),rho_bounds(6),vloc_bounds(6)
  type(s_dg_dc_seed_contract)::contract,repeat_contract,changed_contract
  type(s_dg_dc_seed_payload)::payload,loaded
  real(8),allocatable::rwf(:,:,:,:,:,:,:),rho_owned(:,:,:),vloc_owned(:,:,:)
  real(8),allocatable::esp(:,:,:),rocc(:,:,:)
  real(8)::mu,residual
  real(8)::downstream_localization_tolerance,downstream_pw_cutoff,downstream_lcfo_window
  integer::downstream_w90_iterations
  logical::ok,run_scf,load_seed,publish_seed,fatal,scf_skipped
  character(512)::directory,message
  real(8),parameter::density_weight=0.25d0,expected_electrons=2d0
  real(8),parameter::electron_tolerance=1d-12,current_threshold=1d-8

  call MPI_Init(ierr)
  comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr)
  call MPI_Comm_size(comm,nproc,ierr)
  call get_environment_variable('DG_DC_SEED_DIRECTORY',directory)
  call require(len_trim(directory)>0,'missing DG_DC_SEED_DIRECTORY')

  rwf_bounds=[rank-1,2,-1,1,0,1,1, rank,3,-1,1,1,1,1]
  rho_bounds=[2*rank-1,0,-2, 2*rank,1,-1]
  vloc_bounds=[2*rank-1,0,-2, 2*rank,1,-1]
  immutable_inputs=[17_int64,23_int64,31_int64,37_int64,41_int64,43_int64]
  ownership_map=[int(rank,int64),int(nproc,int64),101_int64,103_int64,&
    int(rho_bounds(1),int64),int(rho_bounds(4),int64),107_int64,109_int64]

  call build_dg_dc_seed_contract(comm,rank+1,rwf_bounds,rho_bounds,vloc_bounds,&
    immutable_inputs,ownership_map,contract,ok,message)
  call require(ok,'contract build failed: '//trim(message))
  call require(contract%mpi_size==nproc.and.contract%rank==rank,&
    'contract lost exact MPI rank topology')
  call require(contract%fragment_id==rank+1,'contract lost rank-to-fragment mapping')
  call require(all(contract%rwf_bounds==rwf_bounds).and.all(contract%rho_bounds==rho_bounds).and.&
    all(contract%vloc_bounds==vloc_bounds),'contract lost exact local array bounds')

  downstream_localization_tolerance=1d-6
  downstream_w90_iterations=71
  downstream_pw_cutoff=3.5d0
  downstream_lcfo_window=0.25d0
  downstream_localization_tolerance=downstream_localization_tolerance*0.1d0
  downstream_w90_iterations=downstream_w90_iterations+100
  downstream_pw_cutoff=downstream_pw_cutoff+4d0
  downstream_lcfo_window=downstream_lcfo_window+0.5d0
  call build_dg_dc_seed_contract(comm,rank+1,rwf_bounds,rho_bounds,vloc_bounds,&
    immutable_inputs,ownership_map,repeat_contract,ok,message)
  call require(ok,'repeat contract build failed: '//trim(message))
  call require(same_contract(contract,repeat_contract),&
    'downstream localization/W90/PW/LCFO/window controls changed the DC seed contract')

  changed_contract=contract
  immutable_inputs(3)=immutable_inputs(3)+1_int64
  call build_dg_dc_seed_contract(comm,rank+1,rwf_bounds,rho_bounds,vloc_bounds,&
    immutable_inputs,ownership_map,changed_contract,ok,message)
  call require(ok,'changed immutable contract build failed: '//trim(message))
  call require(changed_contract%immutable_fingerprint/=contract%immutable_fingerprint,&
    'immutable conventional-DC input did not change the fingerprint')
  immutable_inputs(3)=immutable_inputs(3)-1_int64
  ownership_map(7)=ownership_map(7)+1_int64
  call build_dg_dc_seed_contract(comm,rank+1,rwf_bounds,rho_bounds,vloc_bounds,&
    immutable_inputs,ownership_map,changed_contract,ok,message)
  call require(ok,'changed ownership contract build failed: '//trim(message))
  call require(changed_contract%ownership_fingerprint/=contract%ownership_fingerprint,&
    'owned-grid mapping did not change the ownership fingerprint')
  ownership_map(7)=ownership_map(7)-1_int64

  call initialize_payload(payload)
  if(nproc>1)then
    payload%mu=-0.5d0+real(rank,8)/1000d0
    call write_dg_dc_seed(comm,trim(directory),contract,payload,density_weight,expected_electrons,&
      electron_tolerance,current_threshold,publication_id,ok,message)
    call require(.not.ok,'rank-inconsistent chemical potential was publishable')
    payload%mu=-0.5d0
    payload%residual=1d-10+real(rank,8)*1d-12
    call write_dg_dc_seed(comm,trim(directory),contract,payload,density_weight,expected_electrons,&
      electron_tolerance,current_threshold,publication_id,ok,message)
    call require(.not.ok,'rank-inconsistent convergence residual was publishable')
    payload%residual=1d-10
    payload%iteration=19+rank
    call write_dg_dc_seed(comm,trim(directory),contract,payload,density_weight,expected_electrons,&
      electron_tolerance,current_threshold,publication_id,ok,message)
    call require(.not.ok,'rank-inconsistent iteration provenance was publishable')
    payload%iteration=19
  endif
  call write_dg_dc_seed(comm,trim(directory),contract,payload,density_weight,expected_electrons,&
    electron_tolerance,current_threshold,publication_id,ok,message)
  call require(ok,'seed publication failed: '//trim(message))
  call require(publication_id/=0_int64,'seed publication returned a zero ID')
  call probe_dg_dc_seed(comm,trim(directory),contract,density_weight,expected_electrons,&
    electron_tolerance,current_threshold,status,publication_id,message)
  call require(status==DG_DC_SEED_VALID,'published seed did not probe valid: '//trim(message))

  changed_contract=contract
  changed_contract%mpi_size=nproc+1
  call require_invalid(changed_contract,'MPI rank-count mismatch was accepted')
  changed_contract=contract
  if(rank==0)changed_contract%rank=changed_contract%rank+1
  call require_invalid(changed_contract,'rank identity mismatch was accepted')
  changed_contract=contract
  if(rank==0)changed_contract%fragment_id=changed_contract%fragment_id+1
  call require_invalid(changed_contract,'rank-to-fragment mismatch was accepted')
  changed_contract=contract
  if(rank==nproc-1)changed_contract%rho_bounds(4)=changed_contract%rho_bounds(4)+1
  call require_invalid(changed_contract,'rank-local owned bounds mismatch was accepted')
  changed_contract=contract
  if(rank==0)changed_contract%ownership_fingerprint=changed_contract%ownership_fingerprint+1_int64
  call require_invalid(changed_contract,'ownership fingerprint mismatch was accepted')
  changed_contract=contract
  if(rank==nproc-1)changed_contract%immutable_fingerprint=changed_contract%immutable_fingerprint+1_int64
  call require_invalid(changed_contract,'immutable fingerprint mismatch was accepted')

  call read_dg_dc_seed(comm,trim(directory),contract,density_weight,expected_electrons,&
    electron_tolerance,current_threshold,loaded,publication_id,ok,message)
  call require(ok,'valid seed read failed: '//trim(message))
  allocate(rwf(1:1,1:1,1:1,1:1,1:1,1:1,1:1),rho_owned(1:1,1:1,1:1),&
    vloc_owned(1:1,1:1,1:1),esp(1:1,1:1,1:1),rocc(1:1,1:1,1:1))
  rwf=-1d0;rho_owned=-2d0;vloc_owned=-3d0;esp=-4d0;rocc=-5d0
  mu=huge(0d0);residual=huge(0d0);iteration=-1
  call restore_dg_dc_seed_payload(loaded,rwf,rho_owned,vloc_owned,esp,rocc,mu,residual,&
    iteration,ok,message)
  call require(ok,'seed state restoration failed: '//trim(message))
  call require(same_rank7(rwf,payload%rwf),'rwf values or exact bounds were not restored')
  call require(same_rank3(rho_owned,payload%rho_tot),&
    'owned density values or exact bounds were not restored')
  call require(same_rank3(vloc_owned,payload%vloc_tot),&
    'owned local-potential values or exact bounds were not restored')
  call require(same_rank3(esp,payload%esp),'orbital energies or exact bounds were not restored')
  call require(same_rank3(rocc,payload%rocc),'occupations or exact bounds were not restored')
  call require(mu==payload%mu.and.residual==payload%residual.and.iteration==payload%iteration,&
    'mu/residual/iteration were not restored exactly')

  call require_mode('off',DG_DC_SEED_VALID,.true.,.false.,.false.,.false.,.false.)
  call require_mode('write',DG_DC_SEED_VALID,.true.,.false.,.true.,.false.,.false.)
  call require_mode('read',DG_DC_SEED_VALID,.false.,.true.,.false.,.false.,.true.)
  call require_mode('read',DG_DC_SEED_ABSENT,.false.,.false.,.false.,.true.,.false.)
  call require_mode('read',DG_DC_SEED_INVALID,.false.,.false.,.false.,.true.,.false.)
  call require_mode('auto',DG_DC_SEED_VALID,.false.,.true.,.false.,.false.,.true.)
  call require_mode('auto',DG_DC_SEED_ABSENT,.true.,.false.,.true.,.false.,.false.)
  call require_mode('auto',DG_DC_SEED_INVALID,.false.,.false.,.false.,.true.,.false.)

  if(rank==0)write(*,'(a,i0)')'PASS DG DC seed state ranks=',nproc
  call MPI_Finalize(ierr)

contains
  subroutine initialize_payload(value)
    type(s_dg_dc_seed_payload),intent(out)::value
    integer::i
    allocate(value%rwf(rwf_bounds(1):rwf_bounds(8),rwf_bounds(2):rwf_bounds(9),&
      rwf_bounds(3):rwf_bounds(10),rwf_bounds(4):rwf_bounds(11),&
      rwf_bounds(5):rwf_bounds(12),rwf_bounds(6):rwf_bounds(13),&
      rwf_bounds(7):rwf_bounds(14)))
    allocate(value%rho_tot(rho_bounds(1):rho_bounds(4),rho_bounds(2):rho_bounds(5),&
      rho_bounds(3):rho_bounds(6)))
    allocate(value%vloc_tot(vloc_bounds(1):vloc_bounds(4),vloc_bounds(2):vloc_bounds(5),&
      vloc_bounds(3):vloc_bounds(6)))
    allocate(value%esp(-2:0,2:2,rank:rank),value%rocc(-2:0,2:2,rank:rank))
    value%rwf=reshape([(real(10000*rank+i,8)/1000d0,i=1,size(value%rwf))],shape(value%rwf))
    ! Eight owned density values per rank; the global integral is exactly two electrons.
    value%rho_tot=1d0/real(nproc,8)
    value%vloc_tot=reshape([(real(1000*rank+i,8)/100d0,i=1,size(value%vloc_tot))],&
      shape(value%vloc_tot))
    value%esp(:,2,rank)=[-0.875d0,-0.375d0,0.625d0]+real(rank,8)/100d0
    value%rocc(:,2,rank)=[2d0,0d0,0d0]
    value%mu=-0.5d0
    value%residual=1d-10
    value%iteration=19
  end subroutine initialize_payload

  subroutine require_invalid(value,text)
    type(s_dg_dc_seed_contract),intent(in)::value
    character(*),intent(in)::text
    integer::observed_status
    integer(int64)::observed_id
    call probe_dg_dc_seed(comm,trim(directory),value,density_weight,expected_electrons,&
      electron_tolerance,current_threshold,observed_status,observed_id,message)
    call require(observed_status==DG_DC_SEED_INVALID,text//': '//trim(message))
  end subroutine require_invalid

  subroutine require_mode(mode,seed_status,want_run,want_load,want_publish,want_fatal,want_skip)
    character(*),intent(in)::mode
    integer,intent(in)::seed_status
    logical,intent(in)::want_run,want_load,want_publish,want_fatal,want_skip
    call resolve_dg_dc_seed_mode(mode,seed_status,run_scf,load_seed,publish_seed,fatal,&
      scf_skipped,ok,message)
    call require(ok.neqv.want_fatal,'mode resolver success/fatal contract is inconsistent for '//trim(mode))
    call require(run_scf.eqv.want_run,'wrong run_scf decision for '//trim(mode))
    call require(load_seed.eqv.want_load,'wrong load_seed decision for '//trim(mode))
    call require(publish_seed.eqv.want_publish,'wrong publish decision for '//trim(mode))
    call require(fatal.eqv.want_fatal,'wrong fatal decision for '//trim(mode))
    call require(scf_skipped.eqv.want_skip,'wrong scf_skipped decision for '//trim(mode))
  end subroutine require_mode

  logical function same_contract(left,right)
    type(s_dg_dc_seed_contract),intent(in)::left,right
    same_contract=left%version==right%version.and.left%mpi_size==right%mpi_size.and.&
      left%rank==right%rank.and.left%fragment_id==right%fragment_id.and.&
      all(left%rwf_bounds==right%rwf_bounds).and.all(left%rho_bounds==right%rho_bounds).and.&
      all(left%vloc_bounds==right%vloc_bounds).and.&
      left%immutable_fingerprint==right%immutable_fingerprint.and.&
      left%ownership_fingerprint==right%ownership_fingerprint
  end function same_contract

  logical function same_rank7(left,right)
    real(8),allocatable,intent(in)::left(:,:,:,:,:,:,:),right(:,:,:,:,:,:,:)
    same_rank7=allocated(left).and.allocated(right)
    if(.not.same_rank7)return
    same_rank7=all(lbound(left)==lbound(right)).and.all(ubound(left)==ubound(right)).and.&
      all(left==right)
  end function same_rank7

  logical function same_rank3(left,right)
    real(8),allocatable,intent(in)::left(:,:,:),right(:,:,:)
    same_rank3=allocated(left).and.allocated(right)
    if(.not.same_rank3)return
    same_rank3=all(lbound(left)==lbound(right)).and.all(ubound(left)==ubound(right)).and.&
      all(left==right)
  end function same_rank3

  subroutine require(condition,text)
    logical,intent(in)::condition
    character(*),intent(in)::text
    integer::global_condition
    global_condition=merge(1,0,condition)
    call MPI_Allreduce(MPI_IN_PLACE,global_condition,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(global_condition==0)then
      if(rank==0)write(*,'(a)')'FAIL '//trim(text)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require
end program test_dg_dc_seed_state_mpi
