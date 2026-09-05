program test_dg_hybrid_schwarz_operator_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_quiet_nan,ieee_value
  use dg_hybrid_schwarz_operator,only:s_dg_hybrid_schwarz_schedule,&
    build_dg_hybrid_schwarz_schedule,apply_dg_hybrid_schwarz_rows,&
    apply_dg_hybrid_schwarz_hamiltonian
  implicit none
  integer::ierr,rank,nproc,fragment,global_count,local_count,f,p,q,g,slot,first_global,stat
  integer(int64)::accepted_fingerprint
  integer,allocatable::counts(:),offsets(:),owners(:),fragments(:),local_slots(:),generations(:),manifest(:)
  integer(int64),allocatable::row_ids(:)
  complex(real64),allocatable::interface_rows(:,:)
  complex(real64),allocatable::metric_rows(:,:),kinetic_rows(:,:),nonlocal_rows(:,:),potential_rows(:,:),&
    local_coefficients(:,:),full_coefficients(:,:),distributed_h(:,:),distributed_s(:,:),&
    expected_h(:,:),expected_s(:,:),accepted_h(:,:),h_zero(:,:),h_half(:,:),h_full(:,:)
  type(s_dg_hybrid_schwarz_schedule)::schedule
  logical::ok,requests_match,slots_match
  integer::state,exchange_count
  character(512)::message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call require(any(nproc==[2,4,8]),'test requires 2, 4, or 8 ranks')
  fragment=rank+1
  allocate(counts(nproc),offsets(nproc+1))
  counts=[(f+1,f=1,nproc)];offsets(1)=0
  do f=1,nproc;offsets(f+1)=offsets(f)+counts(f);enddo
  global_count=offsets(nproc+1);local_count=counts(fragment);first_global=offsets(fragment)+1
  allocate(owners(global_count),fragments(global_count),local_slots(global_count),generations(global_count))
  do f=1,nproc
    do slot=1,counts(f)
      p=offsets(f)+slot;owners(p)=f-1;fragments(p)=f;local_slots(p)=slot;generations(p)=3
    enddo
  enddo
  allocate(row_ids(local_count),interface_rows(local_count,global_count))
  row_ids=[(int(first_global+slot-1,int64),slot=1,local_count)]
  interface_rows=(0d0,0d0)
  if(fragment>1)interface_rows(1,offsets(fragment-1)+1)=cmplx(-0.25d0,-0.01d0*(fragment-1),real64)
  if(fragment<nproc)interface_rows(1,offsets(fragment+1)+1)=cmplx(-0.25d0,0.01d0*fragment,real64)

  call expected_manifest(fragment,nproc,manifest)
  call build_dg_hybrid_schwarz_schedule(MPI_COMM_WORLD,fragment,3,row_ids,owners,fragments,local_slots,&
    generations,interface_rows,81231_int64,93451_int64,73691_int64,schedule,ok,message,&
    declared_neighbors=manifest)
  call require(ok,'valid line schedule rejected: '//trim(message))
  call require(schedule%valid.and.schedule%fragment_id==fragment,'schedule was not published')
  call require(schedule%peer_count==size(manifest),'peer count differs from line topology')
  call require(all(schedule%peers==manifest),'peer ordering differs from fragment ordering')
  call require(size(schedule%receive_global_ids)==size(manifest),'unexpected receive request count')
  call require(size(schedule%send_local_slots)==size(manifest),'unexpected send request count')
  requests_match=.true.;slots_match=.true.
  do p=1,size(manifest)
    g=manifest(p)
    requests_match=requests_match.and.schedule%receive_global_ids(p)==int(offsets(g)+1,int64)
    slots_match=slots_match.and.schedule%send_local_slots(p)==1
  enddo
  call require(requests_match,'receive request does not name the adjacent first basis column')
  call require(slots_match,'send request does not map to local slot one')
  call require(schedule%fingerprint/=0_int64,'schedule fingerprint is zero')
  call fingerprint_agrees(schedule%fingerprint)
  accepted_fingerprint=schedule%fingerprint

  allocate(metric_rows(local_count,global_count),kinetic_rows(local_count,global_count),&
    nonlocal_rows(local_count,global_count),potential_rows(local_count,global_count),&
    local_coefficients(local_count,3),full_coefficients(global_count,3),&
    expected_h(local_count,3),expected_s(local_count,3))
  metric_rows=(0d0,0d0);kinetic_rows=(0d0,0d0);nonlocal_rows=(0d0,0d0);potential_rows=(0d0,0d0)
  do p=1,local_count
    q=int(row_ids(p));metric_rows(p,q)=cmplx(1d0+0.05d0*fragment,0d0,real64)
    kinetic_rows(p,q)=cmplx(0.7d0+0.03d0*q,0d0,real64)
    nonlocal_rows(p,q)=cmplx(0.04d0*fragment,0d0,real64)
    potential_rows(p,q)=cmplx(0.11d0+0.01d0*p,0d0,real64)
  enddo
  do q=1,global_count
    do state=1,3
      full_coefficients(q,state)=cmplx(0.1d0*fragments(q)+0.01d0*local_slots(q),&
        -0.02d0*state+0.001d0*q,real64)
    enddo
  enddo
  local_coefficients=full_coefficients(first_global:first_global+local_count-1,:)
  expected_h=matmul(kinetic_rows+nonlocal_rows+potential_rows,full_coefficients)
  expected_s=matmul(metric_rows,full_coefficients)
  call apply_dg_hybrid_schwarz_hamiltonian(MPI_COMM_WORLD,schedule,3,81231_int64,93451_int64,&
    73691_int64,row_ids,fragments,local_slots,kinetic_rows,nonlocal_rows,interface_rows,potential_rows,&
    0d0,local_coefficients,distributed_h,exchange_count,ok,message)
  call require(ok,'zero-interface Hamiltonian action failed: '//trim(message))
  call require(exchange_count==0,'zero-interface Hamiltonian unexpectedly contacted a peer')
  call require(maxval(abs(distributed_h-expected_h))<1d-12,'H(0) contains an interface contribution')
  allocate(h_zero,source=distributed_h)
  call require_hermitian_action(local_coefficients,h_zero,'H(0) is not Hermitian')

  expected_h=matmul(kinetic_rows+nonlocal_rows+potential_rows+0.5d0*interface_rows,full_coefficients)
  call apply_dg_hybrid_schwarz_hamiltonian(MPI_COMM_WORLD,schedule,3,81231_int64,93451_int64,&
    73691_int64,row_ids,fragments,local_slots,kinetic_rows,nonlocal_rows,interface_rows,potential_rows,&
    0.5d0,local_coefficients,distributed_h,exchange_count,ok,message)
  call require(ok,'half-interface Hamiltonian action failed: '//trim(message))
  call require(exchange_count==schedule%peer_count,'half-interface Hamiltonian missed a peer')
  call require(maxval(abs(distributed_h-expected_h))<1d-12,'H(0.5) does not scale the complete interface action')
  allocate(h_half,source=distributed_h)
  call require_hermitian_action(local_coefficients,h_half,'H(0.5) is not Hermitian')

  expected_h=matmul(kinetic_rows+nonlocal_rows+potential_rows+interface_rows,full_coefficients)
  call apply_dg_hybrid_schwarz_hamiltonian(MPI_COMM_WORLD,schedule,3,81231_int64,93451_int64,&
    73691_int64,row_ids,fragments,local_slots,kinetic_rows,nonlocal_rows,interface_rows,potential_rows,&
    1d0,local_coefficients,distributed_h,exchange_count,ok,message)
  call require(ok,'full-interface Hamiltonian action failed: '//trim(message))
  call require(exchange_count==schedule%peer_count,'Hamiltonian contacted a non-neighbor or missed a peer')
  call require(maxval(abs(distributed_h-expected_h))<1d-12,&
    'neighbor-only Hamiltonian differs from explicit complete DG matrix')
  allocate(h_full,source=distributed_h)
  call require(maxval(abs(h_half-(h_zero+0.5d0*(h_full-h_zero))))<1d-12,&
    'H(0.5) is not H(0) plus half the complete interface action')
  call require_hermitian_action(local_coefficients,h_full,'H(1) is not Hermitian')
  call apply_dg_hybrid_schwarz_rows(MPI_COMM_WORLD,schedule,3,81231_int64,93451_int64,73691_int64,&
    row_ids,fragments,local_slots,metric_rows,local_coefficients,distributed_s,exchange_count,ok,message)
  call require(ok,'neighbor-only metric action failed: '//trim(message))
  call require(exchange_count==0,'block-local metric unexpectedly exchanged coefficient rows')
  call require(maxval(abs(distributed_s-expected_s))<1d-12,&
    'distributed metric differs from explicit complete DG matrix')
  allocate(accepted_h,source=distributed_h)

  call reject_scale(-0.1d0,'negative interface scale was accepted')
  call reject_scale(1.1d0,'interface scale above one was accepted')
  call reject_scale(ieee_value(0d0,ieee_quiet_nan),'nonfinite interface scale was accepted')
  call reject_scale(merge(0.5d0,0.25d0,rank==0),'rank-disagreeing interface scale was accepted')

  distributed_h=h_full;accepted_h=h_full
  call apply_dg_hybrid_schwarz_hamiltonian(MPI_COMM_WORLD,schedule,3,81231_int64,93452_int64,&
    73691_int64,row_ids,fragments,local_slots,kinetic_rows,nonlocal_rows,interface_rows,potential_rows,&
    1d0,local_coefficients,distributed_h,exchange_count,ok,message)
  call require(.not.ok,'stale face fingerprint was accepted by Hamiltonian action')
  call require(all(distributed_h==accepted_h),'failed Hamiltonian action changed accepted output')

  if(size(manifest)>0)then
    call build_dg_hybrid_schwarz_schedule(MPI_COMM_WORLD,fragment,3,row_ids,owners,fragments,local_slots,&
      generations,interface_rows,81231_int64,93451_int64,73691_int64,schedule,ok,message,&
      declared_neighbors=manifest(:size(manifest)-1))
    call require(.not.ok,'missing declared neighbor was accepted')
    call require(schedule%fingerprint==accepted_fingerprint,'missing-neighbor failure changed schedule')
  endif

  deallocate(manifest);allocate(manifest(2),stat=stat)
  manifest=[merge(2,1,fragment==1),merge(2,1,fragment==1)]
  call build_dg_hybrid_schwarz_schedule(MPI_COMM_WORLD,fragment,3,row_ids,owners,fragments,local_slots,&
    generations,interface_rows,81231_int64,93451_int64,73691_int64,schedule,ok,message,&
    declared_neighbors=manifest)
  call require(.not.ok,'duplicate declared neighbor was accepted')
  call require(schedule%fingerprint==accepted_fingerprint,'duplicate-neighbor failure changed schedule')

  if(nproc>=4)then
    deallocate(manifest);allocate(manifest(1));manifest=modulo(fragment+1,nproc)+1
    if(manifest(1)==fragment)manifest(1)=modulo(fragment+2,nproc)+1
    call build_dg_hybrid_schwarz_schedule(MPI_COMM_WORLD,fragment,3,row_ids,owners,fragments,local_slots,&
      generations,interface_rows,81231_int64,93451_int64,73691_int64,schedule,ok,message,&
      declared_neighbors=manifest)
    call require(.not.ok,'non-neighbor destination was accepted')
    call require(schedule%fingerprint==accepted_fingerprint,'non-neighbor failure changed schedule')
  endif

  generations(global_count)=4
  call build_dg_hybrid_schwarz_schedule(MPI_COMM_WORLD,fragment,3,row_ids,owners,fragments,local_slots,&
    generations,interface_rows,81231_int64,93451_int64,73691_int64,schedule,ok,message)
  call require(.not.ok,'stale basis generation was accepted')
  call require(schedule%fingerprint==accepted_fingerprint,'stale-generation failure changed schedule')
  generations(global_count)=3

  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid Schwarz schedule on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine reject_scale(scale,detail)
    real(real64),intent(in)::scale
    character(*),intent(in)::detail
    distributed_h=h_full;accepted_h=h_full
    call apply_dg_hybrid_schwarz_hamiltonian(MPI_COMM_WORLD,schedule,3,81231_int64,93451_int64,&
      73691_int64,row_ids,fragments,local_slots,kinetic_rows,nonlocal_rows,interface_rows,potential_rows,&
      scale,local_coefficients,distributed_h,exchange_count,ok,message)
    call require(.not.ok,detail)
    call require(all(distributed_h==accepted_h),'rejected interface scale changed accepted output')
  end subroutine reject_scale
  subroutine require_hermitian_action(coefficients,action,detail)
    complex(real64),intent(in)::coefficients(:,:),action(:,:)
    character(*),intent(in)::detail
    complex(real64)::local_xy,local_yx,global_xy,global_yx
    integer::code
    local_xy=sum(conjg(coefficients(:,1))*action(:,2))
    local_yx=sum(conjg(action(:,1))*coefficients(:,2))
    call MPI_Allreduce(local_xy,global_xy,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_Allreduce(local_yx,global_yx,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,code)
    call require(abs(global_xy-global_yx)<1d-12,detail)
  end subroutine require_hermitian_action
  subroutine expected_manifest(id,count,neighbors)
    integer,intent(in)::id,count
    integer,allocatable,intent(out)::neighbors(:)
    integer::n
    n=merge(1,0,id>1)+merge(1,0,id<count);allocate(neighbors(n));n=0
    if(id>1)then;n=n+1;neighbors(n)=id-1;endif
    if(id<count)then;n=n+1;neighbors(n)=id+1;endif
  end subroutine expected_manifest
  subroutine fingerprint_agrees(value)
    integer(int64),intent(in)::value
    integer(int64)::minimum,maximum
    integer::code
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,code)
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,code)
    call require(minimum==maximum,'schedule fingerprint differs between ranks')
  end subroutine fingerprint_agrees
  subroutine require(condition,detail)
    logical,intent(in)::condition
    character(*),intent(in)::detail
    integer::local_failure,global_failure,code
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,code)
    if(global_failure/=0)then
      if(.not.condition)write(0,'(a,i0,2a)')'rank ',rank,': ',trim(detail)
      call MPI_Abort(MPI_COMM_WORLD,1,code)
    endif
  end subroutine require
end program test_dg_hybrid_schwarz_operator_mpi
