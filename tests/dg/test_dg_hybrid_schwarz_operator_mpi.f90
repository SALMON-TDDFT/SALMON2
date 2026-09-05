program test_dg_hybrid_schwarz_operator_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_schwarz_operator,only:s_dg_hybrid_schwarz_schedule,&
    build_dg_hybrid_schwarz_schedule
  implicit none
  integer::ierr,rank,nproc,fragment,global_count,local_count,f,p,g,slot,first_global,stat
  integer(int64)::accepted_fingerprint
  integer,allocatable::counts(:),offsets(:),owners(:),fragments(:),local_slots(:),generations(:),manifest(:)
  integer(int64),allocatable::row_ids(:)
  complex(real64),allocatable::interface_rows(:,:)
  type(s_dg_hybrid_schwarz_schedule)::schedule
  logical::ok,requests_match,slots_match
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
  if(fragment>1)interface_rows(1,offsets(fragment-1)+1)=cmplx(-0.25d0,0.01d0*fragment,real64)
  if(fragment<nproc)interface_rows(1,offsets(fragment+1)+1)=cmplx(-0.25d0,-0.01d0*fragment,real64)

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
