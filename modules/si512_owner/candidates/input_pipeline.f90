! Experimental owner/tile exchange. Orbital values never use collectives.
module wannier_neighbor_owner
  use iso_c_binding
  use iso_fortran_env, only: int64
  use mpi
!$ use omp_lib, only: omp_get_max_threads,omp_get_thread_num,omp_in_parallel
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use wannier_pair_owner, only: orbital_owners
  use wannier_patch_owner, only: orbital_patch,patch_metrics,pack_exact_patch
  implicit none
  private
  public :: orbital_patch,patch_metrics,pack_exact_patch,patch_apply,neighbor_metrics
  public :: neighbor_phase_names,prefetch_column_capacity
  character(16),parameter :: neighbor_phase_names(13)=[character(16) :: &
    'setup','tile_catalog','allocation','halo_extract','fetch_mpi','return_action', &
    'pair_compute','thread_reduce','finish','total','halo_cache','fft_plans','other']
  include 'fftw3.f03'
  type neighbor_metrics
    integer(int64) :: tiles=0,pair_tiles=0,plans=0,source_messages=0,return_messages=0
    integer(int64) :: prefetch_windows=0,prefetch_fallbacks=0,prefetch_parallel=0
    integer(int64) :: prefetch_payload_bytes=0,prefetch_metadata_bytes=0,prefetch_active_threads=0
    real(8) :: prefetch_plan_seconds=0d0,prefetch_extract_seconds=0d0
    real(8) :: prefetch_pack_seconds=0d0,prefetch_mpi_seconds=0d0
    integer(int64) :: source_values=0,return_values=0
    integer(int64) :: return_windows=0,return_apply_parallel=0,return_metadata_bytes=0
    integer(int64) :: return_plan_hits=0,return_plan_bytes=0,control_skips=0
    real(8) :: return_pack_seconds=0d0,return_mpi_seconds=0d0,return_apply_seconds=0d0
    integer :: omp_requested=1,omp_effective=1,omp_active=0
    integer(int64) :: thread_workspace_bytes=0,return_workspace_bytes=0
    real(8) :: fft_setup_seconds=0d0,pair_compute_seconds=0d0
    integer(int64) :: block_candidates=0,block_visits=0,fetch_requests=0
    integer(int64) :: receive_hits=0,extract_hits=0,extracted_values=0,halo_cache_bytes=0
    real(8) :: phase_min(13)=0d0,phase_max(13)=0d0,phase_mean(13)=0d0
  end type
  type return_queue_plan
    integer :: count=0,groups=0
    integer,allocatable :: next(:),groupnext(:),head(:),worker(:),owner(:),grouphead(:)
    integer,allocatable :: datatype(:),columns(:)
    logical,allocatable :: own_type(:)
  end type
  type message_type_cache

    integer :: block=0,datatype=MPI_DATATYPE_NULL
    integer(kind=MPI_ADDRESS_KIND),allocatable :: offsets(:)
  end type
  type source_window_plan
    integer :: count=0,groups=0,last_task=0
    integer,allocatable :: orb(:),worker(:),next(:),sender(:),dest(:),head(:),ncols(:)
    integer,allocatable :: recv_slot(:),send_slot(:),recv_type(:),send_type(:)
  end type source_window_plan
  type fft_cache

    integer :: shape(3)=0
    complex(c_double_complex),allocatable :: kernel(:)
    type(c_ptr),allocatable :: forward(:),backward(:)
  end type
contains
  pure integer function prefetch_column_capacity(budget,metadata,nfft,requests) result(columns)
    integer(int64),intent(in) :: budget,metadata
    integer,intent(in) :: nfft,requests
    columns=0
    if(nfft<1.or.requests<1.or.metadata<0.or.budget<metadata)return
    ! Traditional MPI counts are default integers; every grouped message must fit.
    columns=int(min(int(requests,int64),(budget-metadata)/(48_int64*nfft), &
      int(huge(0),int64)/int(nfft,int64)))
  end function prefetch_column_capacity

  pure subroutine active_masks(hi,ci,hj,cj,diagonal,left,right)
    integer(int64),intent(in) :: hi,ci,hj,cj
    logical,intent(in) :: diagonal
    integer(int64),intent(out) :: left,right
    left=0;right=0
    if(hi==0.or.hj==0)return
    ! Diagonal fetch and return use the union of both endpoints.
    ! Every halo orbital pairs with any core orbital in this same block.
    if(diagonal)then
      if(ci/=0)left=hi
      return
    endif
    if(ci/=0)then
      left=ci;right=hj
    endif
    if(cj/=0)then
      left=hi;right=ior(right,cj)
    endif
  end subroutine active_masks

  subroutine patch_apply(shape,norb,ids,patches,kernel,comm,batch,action,status,metrics, &
    detail,tile_shape,halo_cache_bytes,return_window_bytes,halo_prefetch_bytes)
    integer,intent(in) :: shape(3),norb,ids(:),comm,batch
    type(orbital_patch),intent(in) :: patches(:)
    real(8),intent(in) :: kernel(:)
    complex(8),intent(out) :: action(:,:)
    integer,intent(out) :: status
    type(patch_metrics),intent(out) :: metrics
    type(neighbor_metrics),intent(out),optional :: detail
    integer,intent(in),optional :: tile_shape(3)
    ! Per-cache payload cap (owner and receiver each), 0 disables reuse. Values
    ! are never cached across tiles or calls. This is not a physics threshold.
    integer,intent(in),optional :: halo_cache_bytes
    ! Two complex payload arrays together are bounded by this cap. Metadata
    ! is separate and bounded by 65536 records plus orbital heads/tails.
    integer,intent(in),optional :: return_window_bytes
    ! Total new prefetch workspace: metadata + extraction + retained receives
    ! + packing scratch. 0 disables it. Existing workspaces are separately capped.
    integer(int64),intent(in),optional :: halo_prefetch_bytes
    type(fft_cache) :: cache(8)
    type(message_type_cache) :: message_cache(1024)
    integer,allocatable :: requests(:),datatypes(:)
    logical,allocatable :: retained_types(:)
    integer(kind=MPI_ADDRESS_KIND),allocatable :: offsets(:)
    type(return_queue_plan) :: return_plans(64)
    type(source_window_plan) :: source_plans(64)
    integer,allocatable :: source_flags(:,:),source_recv_slot(:),source_send_slot(:)
    integer :: source_nfft,source_block,source_slots,source_index,source_windows
    integer(int64) :: source_bytes
    logical :: source_replay,source_record,source_ready
    integer,allocatable :: return_plan_flags(:,:)
    integer :: return_plan_core,return_plan_block,return_plan_capacity,return_plan_windows,return_plan_index
    integer(int64) :: return_plan_bytes,return_plan_hits
    logical :: return_plan_ready,return_plan_replay,return_plan_record
    integer(int64) :: message_cache_bytes
    integer,allocatable :: task_flags(:,:),task_hj(:),task_hi(:),task_bi(:),task_bj(:)
    integer(int64),allocatable :: task_left(:),task_right(:),task_prefix(:)
    integer :: task_block,task_count,task_pos,scan_hj,scan_hi,scan_core,scan_first,scan_last,jump_end,next_owned
    logical :: task_ready,task_found


    type(neighbor_metrics) :: counts
    integer :: rank,np,e,bad,par(11),mn(11),mx(11),edge(3),rh(3),ng,nlocal
    integer :: a,i,j,k,q,neg,x,y,z,v(3),p(3),nt(3),start(3),core(3),origin(3),length(3)
    integer :: worker,tile,nfft,ncore,ncache,slot,bi,bj,ni,nj,ii,jj,block,task
    integer(int64),allocatable :: halo_mask(:),core_mask(:)
    integer(int64) :: left_mask,right_mask
    integer :: ib,jb,offset
    integer :: nblocks,nhalo,ncoreblocks,hj,hi,first_core,first_i,last_i,cache_limit
    integer :: cache_slots,extract_slots,extract_next
    integer :: prefill_count,prefill_index
    integer,allocatable :: halo_blocks(:),core_blocks(:),receive_keys(:,:),receive_next(:),extract_keys(:)
    logical,allocatable :: block_core(:)
    complex(8),allocatable :: received_halo(:,:),extracted_halo(:,:)
    logical :: needed,dense_bounds,mpi_main,return_direct
    integer :: requested,nteam,tid,provided,max_fft,max_core,block_threads,thread_stats(3),all_thread_stats(3)
    integer,allocatable :: hits(:)
    integer(int64) :: core_bound,per_thread,thread_bytes,all_thread_bytes,n_pairs,n_eval,n_zero
    complex(8),allocatable :: work(:,:),spectrum(:,:),thread_i(:,:,:),thread_j(:,:,:)
    logical,allocatable :: active_i(:),active_j(:)
    integer,allocatable :: owners(:),local_index(:),global_flags(:,:),catalog(:,:),all_catalog(:,:)
    integer,allocatable :: rows(:),core_rows(:),core_slots(:)
    real(8),allocatable :: reference_kernel(:)
    complex(8),allocatable :: left(:,:),right(:,:),onto_i(:,:),onto_j(:,:)
    complex(8),allocatable :: owner_buffer(:),return_buffer(:,:)
    complex(8),allocatable,asynchronous :: queued_values(:,:),pending_values(:,:),spare_values(:,:)
    integer,allocatable :: pending_requests(:),pending_types(:),pending_next(:),pending_head(:)
    logical,allocatable :: pending_retained(:)
    integer :: pending_count,pending_nreq
    logical :: pending_active

    integer :: window_limit,queue_capacity,queue_count,group_count,hash_size
    integer,allocatable :: queue_orb(:),queue_next(:),group_next(:),orb_head(:),orb_tail(:)
    integer,allocatable :: group_worker(:),group_owner(:),group_head(:),group_tail(:),group_hash(:),group_slots(:)
    real(8) :: return_times(3),all_return_times(3)
    integer(int64) :: return_stats(2),all_return_stats(2)
    integer(int64) :: volume,lc(7),gc(7),dc(7),dg(7)
    real(8) :: clock_start,local_times(2),all_times(2)
    real(8) :: begin_time,tick,phase(13),extra_times(13),pair_end,reduce_start
    integer(int64) :: extra_counts(9),global_counts(9),cache_memory
    integer :: prefetch_slots,prefetch_until,pf_capacity,pf_count,pf_hash_size,pf_hit_capacity
    integer(int64) :: prefetch_limit,prefetch_min,prefetch_max
    integer(int64) :: pf_metadata,pf_metadata_limit,pf_base,pf_stats(6),pf_all_stats(6)
    integer,allocatable :: pf_thread_hits(:),pf_orb(:),pf_worker(:),pf_hash_slot(:),pf_next(:),pf_hash(:)
    integer,allocatable :: pf_seen(:),pf_local(:),pf_owned(:),pf_workers(:),pf_owners(:)
    integer,allocatable :: pf_sender(:),pf_dest(:),pf_head(:),pf_tail(:),pf_size(:),pf_owned_orb(:),pf_list(:)
    complex(8),allocatable,asynchronous :: prefetch_values(:,:),prefetch_owned(:,:),prefetch_buffer(:,:)
    real(8) :: pf_times(4),pf_all_times(4)
    logical :: prefetch_ready
    complex(8),allocatable,asynchronous :: future_values(:,:),held_values(:,:)
    integer,allocatable :: future_local(:),held_local(:),source_requests(:),source_types(:)
    logical,allocatable :: source_retained(:)
    logical :: future_valid,source_inflight
    integer :: future_until,source_nreq
    begin_time=MPI_Wtime();phase=0d0
    future_valid=.false.;source_inflight=.false.;source_nreq=0;future_until=0
    pending_active=.false.;pending_count=0;pending_nreq=0
    message_cache_bytes=0
    source_ready=.false.;source_replay=.false.;source_record=.false.;source_index=0;source_windows=0;source_bytes=0
    task_ready=.false.;task_count=0;task_pos=0;task_block=0
    return_plan_ready=.false.;return_plan_replay=.false.;return_plan_record=.false.
    return_plan_bytes=0;return_plan_hits=0;return_plan_windows=0;return_plan_index=0

    action=0;metrics=patch_metrics();counts=neighbor_metrics();bad=0;ncache=0
    if(present(detail))detail=counts
    call MPI_Comm_rank(comm,rank,e);call MPI_Comm_size(comm,np,e)
    requested=1
!$  requested=omp_get_max_threads()
!$  if(omp_in_parallel())bad=1
    if(requested>1)then
      call MPI_Query_thread(provided,e);call MPI_Is_thread_main(mpi_main,e)
      if(provided<MPI_THREAD_FUNNELED.or..not.mpi_main)bad=1
    endif
    counts%omp_requested=requested
    edge=[8,8,8]
    if(present(tile_shape))edge=tile_shape
    cache_limit=16777216
    if(present(halo_cache_bytes))cache_limit=halo_cache_bytes
    window_limit=268435456
    if(present(return_window_bytes))window_limit=return_window_bytes
    prefetch_limit=134217728_int64
    if(present(halo_prefetch_bytes))prefetch_limit=halo_prefetch_bytes
    par=[shape,norb,batch,edge,merge(1,0,present(tile_shape)),cache_limit,window_limit]
    call MPI_Allreduce(par,mn,11,MPI_INTEGER,MPI_MIN,comm,e)
    call MPI_Allreduce(par,mx,11,MPI_INTEGER,MPI_MAX,comm,e)
    if(cache_limit<0.or.cache_limit>67108864)bad=1
    if(window_limit<32)bad=1
    call MPI_Allreduce(prefetch_limit,prefetch_min,1,MPI_INTEGER8,MPI_MIN,comm,e)
    call MPI_Allreduce(prefetch_limit,prefetch_max,1,MPI_INTEGER8,MPI_MAX,comm,e)
    if(prefetch_limit<0.or.prefetch_min/=prefetch_max)bad=1
    volume=1_int64
    if(any(mn/=mx).or.any(shape<1).or.any(edge<1).or.norb<1.or.batch<1)bad=1
    do a=1,3
      if(shape(a)<1)cycle
      if(volume>int(huge(ng),int64)/int(shape(a),int64))then
        bad=1
      else
        volume=volume*int(shape(a),int64)
      endif
    enddo
    if(6_int64*int(norb,int64)>int(huge(ng),int64))bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    ng=int(volume);nlocal=size(ids)
    allocate(owners(norb),local_index(norb));call orbital_owners(norb,np,owners)
    local_index=0
    do a=0,np-1
      j=0
      do i=1,norb
        if(owners(i)/=a)cycle
        j=j+1;local_index(i)=j
      enddo
    enddo
    if(size(patches)/=nlocal.or.nlocal/=count(owners==rank))bad=1
    if(size(action,1)/=ng.or.size(action,2)/=nlocal.or.size(kernel)/=ng)bad=1
    j=0
    do i=1,norb
      if(owners(i)/=rank)cycle
      j=j+1
      if(j<=nlocal)then
        if(ids(j)/=i)bad=1
      endif
    enddo
    do j=1,size(patches)
      if(.not.allocated(patches(j)%row).or..not.allocated(patches(j)%value))then
        bad=1;cycle
      endif
      if(size(patches(j)%row)/=size(patches(j)%value))then
        bad=1;cycle
      endif
      do i=1,size(patches(j)%row)
        if(patches(j)%row(i)<1.or.patches(j)%row(i)>ng)bad=1
        if(i>1)then
          if(patches(j)%row(i)<=patches(j)%row(i-1))bad=1
        endif
        if(.not.ieee_is_finite(real(patches(j)%value(i),8)))bad=1
        if(.not.ieee_is_finite(aimag(patches(j)%value(i))))bad=1
      enddo
    enddo
    do i=1,size(kernel)
      if(.not.ieee_is_finite(kernel(i)))bad=1
    enddo
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    allocate(reference_kernel(ng));reference_kernel=kernel
    call MPI_Bcast(reference_kernel,ng,MPI_DOUBLE_PRECISION,0,comm,e)
    if(any(reference_kernel/=kernel))bad=1
    rh=0
    do z=0,shape(3)-1
      do y=0,shape(2)-1
        do x=0,shape(1)-1
          q=1+x+shape(1)*(y+shape(2)*z)
          neg=1+modulo(-x,shape(1))+shape(1)*(modulo(-y,shape(2))+shape(2)*modulo(-z,shape(3)))
          if(kernel(q)/=kernel(neg))bad=1
          if(kernel(q)/=0d0)rh=max(rh,[min(x,shape(1)-x),min(y,shape(2)-y),min(z,shape(3)-z)])
        enddo
      enddo
    enddo
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    ! Exchange only six integer support-bound descriptors per orbital, once per
    ! action. There is no collective in the tile/pair communication loop.
    allocate(global_flags(2,norb),catalog(6,norb),all_catalog(6,norb));catalog=0
    do i=1,nlocal
      call support_bounds(patches(i),shape,catalog(:,ids(i)))
    enddo
    call MPI_Allreduce(catalog,all_catalog,6*norb,MPI_INTEGER,MPI_SUM,comm,e)
    ! Avoid edge + shape overflow. Every axis has at most two distinct core lengths,
    ! hence at most eight cached FFT geometries, independent of the number of tiles.
    if(.not.present(tile_shape))then
      ! If the radius already fills an axis, do not repeat that periodic FFT for
      ! several small cores along the same fully covered axis.
      where(int(edge,int64)+2_int64*int(rh,int64)>=int(shape,int64))edge=shape
      ! A small cell with globally extended supports has no useful support
      ! locality. Avoid repeated almost-global transforms; still distribute its
      ! unordered orbital blocks, and never exceed 32**3 FFT points here.
      dense_bounds=.true.
      do i=1,norb
        if(all(all_catalog(4:6,i)==0))cycle
        if(any(all_catalog(4:6,i)/=shape))dense_bounds=.false.
      enddo
      if(dense_bounds.and.volume<=32768_int64)edge=shape
    endif
    edge=min(edge,shape)
    do a=1,3
      length(a)=int(min(int(shape(a),int64),int(edge(a),int64)+2_int64*int(rh(a),int64)))
    enddo
    if(32_int64*product(int(length,int64))>67108864_int64)bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    max_fft=product(length);max_core=product(edge)
    return_direct=.false.
    if(int(window_limit,int64)<32_int64*max_core)then
      if(present(return_window_bytes))then
        bad=1
      else
        ! Preserve every geometry accepted before batching. A single column
        ! that exceeds the default window uses the original bounded return.
        return_direct=.true.
      endif
    endif
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    if(status/=0)return
    core_bound=67108864_int64
    if(int(batch,int64)<=67108864_int64/(32_int64*max_core)) &
      core_bound=32_int64*max_core*batch
    per_thread=32_int64*max_fft+core_bound
    nteam=int(min(int(requested,int64),max(1_int64,134217728_int64/per_thread), &
      int(batch,int64)*int(batch,int64),int(norb,int64)*int(norb+1,int64)/2_int64))
    counts%omp_effective=nteam
    allocate(hits(nteam))
    ! Planner metadata must be rank-uniform even with heterogeneous OMP teams.
    call MPI_Allreduce(nteam,pf_hit_capacity,1,MPI_INTEGER,MPI_MAX,comm,e)
    pf_capacity=0;pf_metadata=0
    pf_metadata_limit=min(4194304_int64,int(prefetch_limit,int64)/4)
    pf_base=int(storage_size(0)/8,int64)*(3_int64*norb+2_int64*np+1+pf_hit_capacity)
    if(pf_metadata_limit>pf_base)then
      pf_capacity=int(min(65536_int64,int(norb,int64)*np, &
        (pf_metadata_limit-pf_base)/(13_int64*(storage_size(0)/8))))
    endif
    if(pf_capacity>0)then
      allocate(pf_thread_hits(pf_hit_capacity))
      pf_hash_size=2*pf_capacity+1
      allocate(pf_orb(pf_capacity),pf_worker(pf_capacity),pf_hash_slot(pf_capacity),pf_next(pf_capacity))
      allocate(pf_sender(pf_capacity),pf_dest(pf_capacity),pf_head(pf_capacity),pf_tail(pf_capacity),pf_size(pf_capacity))
      allocate(pf_owned_orb(pf_capacity),pf_list(pf_capacity),pf_hash(pf_hash_size))
      allocate(pf_seen(norb),pf_local(norb),pf_owned(norb),pf_workers(np),pf_owners(np))
      pf_metadata=pf_base+13_int64*(storage_size(0)/8)*pf_capacity
    endif
    counts%prefetch_metadata_bytes=pf_metadata
    phase(1)=MPI_Wtime()-begin_time
    nt=1+(shape-1)/edge;tile=0
    do z=0,nt(3)-1
      do y=0,nt(2)-1
        do x=0,nt(1)-1
          tick=MPI_Wtime()
          start=[x,y,z]*edge;core=min(edge,shape-start)
          length=core+2*min(rh,(shape-core)/2)
          do a=1,3
            if(rh(a)>=(shape(a)-core(a)+1)/2)length(a)=shape(a)
          enddo
          origin=modulo(start-rh,shape)
          ! Full-axis transforms use global origin to avoid unnecessary permutations.
          where(length==shape)origin=0
          nfft=product(length);ncore=product(core);worker=modulo(tile,np);tile=tile+1
          allocate(rows(nfft),core_rows(ncore),core_slots(ncore))
!$omp parallel do default(none) schedule(static) num_threads(nteam) if(nfft>=32768.and.nteam>1) &
!$omp shared(nfft,length,origin,shape,rows) private(k,v,p)
          do k=0,nfft-1
            v=[modulo(k,length(1)),modulo(k/length(1),length(2)),k/(length(1)*length(2))]
            p=modulo(origin+v,shape);rows(k+1)=1+p(1)+shape(1)*(p(2)+shape(2)*p(3))
          enddo
!$omp end parallel do
!$omp parallel do default(none) schedule(static) num_threads(nteam) if(ncore>=32768.and.nteam>1) &
!$omp shared(ncore,core,start,shape,origin,length,core_rows,core_slots) private(k,v,p) reduction(max:bad)
          do k=0,ncore-1
            v=[modulo(k,core(1)),modulo(k/core(1),core(2)),k/(core(1)*core(2))]
            p=start+v;core_rows(k+1)=1+p(1)+shape(1)*(p(2)+shape(2)*p(3))
            v=modulo(p-origin,shape)
            if(any(v>=length))bad=1
            core_slots(k+1)=1+v(1)+length(1)*(v(2)+length(2)*v(3))
          enddo
!$omp end parallel do
          global_flags=0
!$omp parallel do default(none) schedule(static) num_threads(nteam) if(norb>=128.and.nteam>1) &
!$omp shared(norb,start,core,all_catalog,shape,origin,length,global_flags) private(i)
          do i=1,norb
            if(boxes_intersect(start,core,all_catalog(1:3,i),all_catalog(4:6,i),shape))global_flags(1,i)=1
            if(boxes_intersect(origin,length,all_catalog(1:3,i),all_catalog(4:6,i),shape))global_flags(2,i)=1
          enddo
!$omp end parallel do
          if(rank==worker)counts%tiles=counts%tiles+1
          ! Two bounded orbital blocks stage each halo once for many pairs.
          ! Keep their complex input buffers within 64 MiB in total. Input batch
          ! is an upper bound; reducing it changes only work grouping.
          block=int(min(int(batch,int64),max(1_int64,67108864_int64/(32_int64*nfft))))
          nblocks=1+(norb-1)/block
          allocate(halo_blocks(nblocks),core_blocks(nblocks),block_core(nblocks))
          nhalo=0;ncoreblocks=0;block_core=.false.
          do ib=1,nblocks
            bi=1+(ib-1)*block;ni=min(block,norb-bi+1)
            if(any(global_flags(2,bi:bi+ni-1)/=0))then
              nhalo=nhalo+1;halo_blocks(nhalo)=ib
            endif
            if(any(global_flags(1,bi:bi+ni-1)/=0))then
              ncoreblocks=ncoreblocks+1;core_blocks(ncoreblocks)=ib;block_core(ib)=.true.
            endif
          enddo
          if(rank==0)counts%block_candidates=counts%block_candidates+ &
            int(nblocks,int64)*(int(nblocks,int64)+1_int64)/2_int64
          ! At most 64 halos per cache, and at most 4 MiB replicated key table.
          cache_slots=int(min(64_int64,int(norb,int64),int(cache_limit,int64)/(16_int64*nfft), &
            1048576_int64/int(np,int64)))
          extract_slots=min(nlocal,cache_slots)
          allocate(receive_keys(cache_slots,np),receive_next(np),extract_keys(extract_slots))
          receive_keys=0;receive_next=1;extract_keys=0;extract_next=1
          prefetch_slots=0;prefetch_until=0;prefetch_ready=.false.
          if(pf_capacity>0)then
            ! Retained current receives plus next receives/extraction/scratch.
            prefetch_slots=int(min(int(pf_capacity,int64), &
              max(0_int64,prefetch_limit-pf_metadata-4_int64*norb)/(64_int64*nfft), &
              int(huge(0),int64)/nfft))
          endif
          prefill_count=0
          if(ncoreblocks>0.and.extract_slots>0.and.prefetch_slots==0)then
            do i=1,nlocal
              if(global_flags(2,ids(i))/=0)prefill_count=prefill_count+1
            enddo
          endif
          if(prefill_count>0.and.prefill_count<=extract_slots)then
            k=0
            do i=1,nlocal
              if(global_flags(2,ids(i))==0)cycle
              k=k+1;extract_keys(k)=ids(i)
            enddo
            phase(2)=phase(2)+MPI_Wtime()-tick;tick=MPI_Wtime()
            allocate(extracted_halo(nfft,extract_slots))
            ! Every halo orbital has at least one core partner when a core
            ! exists. If they all fit, extract each owned halo once with OMP
            ! over whole orbitals, outside the row lookup loop and all MPI.
!$omp parallel do default(none) schedule(static) num_threads(min(nteam,prefill_count)) &
!$omp if(int(nfft,int64)*int(prefill_count,int64)>=32768_int64.and.nteam>1) &
!$omp shared(prefill_count,extract_keys,local_index,patches,rows,nfft,extracted_halo) private(i,k,prefill_index)
            do i=1,prefill_count
              prefill_index=local_index(extract_keys(i))
              do k=1,nfft
                extracted_halo(k,i)=value_at(patches(prefill_index),rows(k))
              enddo
            enddo
!$omp end parallel do
            counts%extracted_values=counts%extracted_values+int(nfft,int64)*prefill_count
            phase(4)=phase(4)+MPI_Wtime()-tick;tick=MPI_Wtime()
          endif
          if(block<=bit_size(left_mask))then
            allocate(halo_mask(1+(norb-1)/block),core_mask(1+(norb-1)/block))
            halo_mask=0;core_mask=0
            do i=1,norb
              if(global_flags(2,i)==0)cycle
              ib=1+(i-1)/block;offset=mod(i-1,block)
              halo_mask(ib)=ibset(halo_mask(ib),offset)
              if(global_flags(1,i)/=0)core_mask(ib)=ibset(core_mask(ib),offset)
            enddo
          endif
          phase(2)=phase(2)+MPI_Wtime()-tick
          queue_capacity=max(1,min(65536,window_limit/(32*ncore)));queue_count=0;group_count=0
          hash_size=2*queue_capacity+1
          allocate(queue_orb(queue_capacity),queue_next(queue_capacity),group_next(queue_capacity))
          allocate(group_worker(queue_capacity),group_owner(queue_capacity),group_head(queue_capacity), &
            group_tail(queue_capacity),group_hash(hash_size),group_slots(queue_capacity))
          allocate(orb_head(norb),orb_tail(norb));orb_head=0;orb_tail=0;group_hash=0
          counts%return_metadata_bytes=max(counts%return_metadata_bytes, &
            int(storage_size(0)/8,int64)*(10_int64*queue_capacity+1_int64+2_int64*norb))
          call begin_source_plans()
          call begin_return_plan()
          task=0
          allocate(active_i(min(block,norb)),active_j(min(block,norb)))
          tick=MPI_Wtime();call begin_task_plan()
          phase(2)=phase(2)+MPI_Wtime()-tick
          do
              tick=MPI_Wtime()
              call next_pair_task(task_found)
              if(.not.task_found)exit
              worker=modulo(modulo(tile-1,np)+modulo(task,np),np);task=task+1
              phase(2)=phase(2)+MPI_Wtime()-tick;tick=MPI_Wtime()
              if(rank==worker.and..not.allocated(left))then
                k=min(block,norb)
                allocate(left(nfft,k),right(nfft,k),onto_i(ncore,k),onto_j(ncore,k))
                allocate(thread_i(ncore,k,nteam),thread_j(ncore,k,nteam))
              endif
              phase(3)=phase(3)+MPI_Wtime()-tick
              if(task>prefetch_until)call advance_input_window(hj,hi,task-1)
              if(bi==bj)active_i=active_i.or.active_j
              ! A complete prefetch window has already transferred every input.
              ! Cached returns retain the global slot order without per-column calls.
              if(rank/=worker.and.prefetch_ready.and.task<=prefetch_until.and.return_plan_replay)then
                k=count(active_i(:ni))
                if(bi/=bj)k=k+count(active_j(:nj))
                call advance_return_queue(k)
                counts%control_skips=counts%control_skips+1
                if(task_ready)then
                  ! task is one-based; future owned task keeps the original rotation.
                  next_owned=task+1+modulo(rank-modulo(tile-1,np)-task,np)
                  jump_end=min(prefetch_until,task_count,next_owned-1)
                  if(jump_end>task)then
                    call advance_return_queue(int(task_prefix(jump_end)-task_prefix(task)))
                    counts%control_skips=counts%control_skips+jump_end-task
                    task=jump_end;task_pos=jump_end
                  endif
                endif
                cycle
              endif
              do ii=1,ni
                if(active_i(ii))call fetch_column(bi+ii-1,left,ii)
              enddo
              if(bi==bj)then
                if(rank==worker)then
                  do jj=1,nj
                    if(active_i(jj))right(:,jj)=left(:,jj)
                  enddo
                endif
              else
                do jj=1,nj
                  if(active_j(jj))call fetch_column(bj+jj-1,right,jj)
                enddo
              endif
              if(rank==worker)then
                if(.not.allocated(work))allocate(work(max_fft,nteam),spectrum(max_fft,nteam))
                slot=0
                do a=1,ncache
                  if(all(cache(a)%shape==length))slot=a
                enddo
                if(slot==0)then
                  ncache=ncache+1;slot=ncache
                  clock_start=MPI_Wtime()
                  call prepare_fft(cache(slot),length,shape,kernel,work,spectrum,bad)
                  counts%fft_setup_seconds=counts%fft_setup_seconds+MPI_Wtime()-clock_start
                  counts%plans=counts%plans+2_int64*nteam
                endif
                tick=MPI_Wtime()
!$omp parallel do default(none) schedule(static) num_threads(min(nteam,ni*nj)) &
!$omp shared(nteam,thread_i,thread_j) private(q)
                do q=1,nteam
                  thread_i(:,:,q)=0;thread_j(:,:,q)=0
                enddo
!$omp end parallel do
                hits=0;n_pairs=0;n_eval=0;n_zero=0
                thread_bytes=32_int64*max_fft*nteam+16_int64* &
                  (size(thread_i,kind=int64)+size(thread_j,kind=int64))
                counts%thread_workspace_bytes=max(counts%thread_workspace_bytes,thread_bytes)
                block_threads=int(min(int(nteam,int64),int(ni,int64)*int(nj,int64)))
                phase(3)=phase(3)+MPI_Wtime()-tick
                clock_start=MPI_Wtime()
                if(bad==0)then
!$omp parallel do collapse(2) schedule(static) num_threads(block_threads) default(none) &
!$omp shared(ni,nj,bi,bj,global_flags,left,right,nfft,ncore,core_slots,cache,slot,work,spectrum, &
!$omp        thread_i,thread_j,hits,block_threads) private(ii,jj,i,j,k,q,tid) reduction(+:n_pairs,n_eval,n_zero)
                  do jj=1,nj
                    do ii=1,ni
                      i=bi+ii-1;j=bj+jj-1;tid=1
!$                    tid=omp_get_thread_num()+1
                      if(.not.pair_needed(i,j))cycle
                      n_pairs=n_pairs+1
                      if(any_nonzero_product(left(:,ii),right(:,jj)))then
                        work(:nfft,tid)=conjg(right(:,jj))*left(:,ii)
                        call fftw_execute_dft(cache(slot)%forward(tid),work(:,tid),spectrum(:,tid))
                        spectrum(:nfft,tid)=spectrum(:nfft,tid)*cache(slot)%kernel
                        call fftw_execute_dft(cache(slot)%backward(tid),spectrum(:,tid),work(:,tid))
                        work(:nfft,tid)=work(:nfft,tid)/real(nfft,8)
                        do k=1,ncore
                          q=core_slots(k)
                          thread_i(k,ii,tid)=thread_i(k,ii,tid)-right(q,jj)*work(q,tid)
                          if(i/=j)then
                            if(bi==bj)then
                              thread_i(k,jj,tid)=thread_i(k,jj,tid)-left(q,ii)*conjg(work(q,tid))
                            else
                              thread_j(k,jj,tid)=thread_j(k,jj,tid)-left(q,ii)*conjg(work(q,tid))
                            endif
                          endif
                        enddo
                        n_eval=n_eval+1;hits(tid)=hits(tid)+1
                      else
                        n_zero=n_zero+1
                      endif
                    enddo
                  enddo
!$omp end parallel do
                endif
                pair_end=MPI_Wtime();phase(7)=phase(7)+pair_end-clock_start
                reduce_start=pair_end
                ! Each output element keeps ascending tid additions. Columns
                ! are independent; never change the pair-loop thread assignment.
                if(int(ncore,int64)*int(ni+nj,int64)*int(nteam,int64)>=262144_int64.and.block_threads>1)then
!$omp parallel do default(none) schedule(static) num_threads(min(block_threads,ni+nj)) &
!$omp shared(onto_i,onto_j,thread_i,thread_j,ni,nj) private(ii)
                  do ii=1,ni+nj
                    if(ii<=ni)then
                      call sum_column(onto_i(:,ii),thread_i(:,ii,:))
                    else
                      call sum_column(onto_j(:,ii-ni),thread_j(:,ii-ni,:))
                    endif
                  enddo
!$omp end parallel do
                else
                  do ii=1,ni
                    call sum_column(onto_i(:,ii),thread_i(:,ii,:))
                  enddo
                  do jj=1,nj
                    call sum_column(onto_j(:,jj),thread_j(:,jj,:))
                  enddo
                endif
                phase(8)=phase(8)+MPI_Wtime()-reduce_start
                counts%pair_compute_seconds=counts%pair_compute_seconds+MPI_Wtime()-clock_start
                counts%pair_tiles=counts%pair_tiles+n_pairs
                counts%omp_active=max(counts%omp_active,count(hits>0))
                metrics%evaluated=metrics%evaluated+n_eval;metrics%zero_overlap=metrics%zero_overlap+n_zero
                metrics%fft_points=metrics%fft_points+n_eval*int(nfft,int64)
                if(n_eval>0)metrics%max_fft_points=max(metrics%max_fft_points,int(nfft,int64))
              endif
              call queue_block(bi,onto_i,active_i(:ni))
              if(bi/=bj)call queue_block(bj,onto_j,active_j(:nj))
          enddo
          call finish_source_requests()
          tick=MPI_Wtime();call flush_returns();call finish_pending_returns()
          phase(6)=phase(6)+MPI_Wtime()-tick
          if(source_record)source_ready=.true.
          if(return_plan_record)then
            return_plan_windows=return_plan_index;return_plan_ready=.true.
          endif
          deallocate(active_i,active_j)
          if(allocated(thread_i))deallocate(thread_i,thread_j)
          deallocate(queue_orb,queue_next,group_next,group_worker,group_owner,group_head,group_tail, &
            group_hash,group_slots,orb_head,orb_tail)
          if(allocated(queued_values))deallocate(queued_values)
          if(allocated(pending_values))deallocate(pending_values)
          if(allocated(spare_values))deallocate(spare_values)
          if(allocated(halo_mask))deallocate(halo_mask,core_mask)
          if(allocated(left))deallocate(left,right,onto_i,onto_j)
          if(allocated(return_buffer))deallocate(return_buffer)
          cache_memory=0
          if(allocated(received_halo))cache_memory=cache_memory+16_int64*size(received_halo,kind=int64)
          if(allocated(extracted_halo))cache_memory=cache_memory+16_int64*size(extracted_halo,kind=int64)
          cache_memory=cache_memory+int(storage_size(0)/8,int64)* &
            (size(receive_keys,kind=int64)+size(receive_next,kind=int64)+size(extract_keys,kind=int64))
          counts%halo_cache_bytes=max(counts%halo_cache_bytes,cache_memory)
          if(allocated(prefetch_owned))deallocate(prefetch_owned)
          if(allocated(prefetch_buffer))deallocate(prefetch_buffer)
          if(allocated(prefetch_values))deallocate(prefetch_values)
          if(allocated(future_values))deallocate(future_values)
          if(allocated(future_local))deallocate(future_local)
          future_valid=.false.
          if(allocated(received_halo))deallocate(received_halo)
          if(allocated(extracted_halo))deallocate(extracted_halo)
          deallocate(receive_keys,receive_next,extract_keys,halo_blocks,core_blocks,block_core)
          deallocate(rows,core_rows,core_slots)
        enddo
      enddo
    enddo
    tick=MPI_Wtime()
    do a=1,ncache
      do tid=1,size(cache(a)%forward)
        if(c_associated(cache(a)%forward(tid)))call fftw_destroy_plan(cache(a)%forward(tid))
        if(c_associated(cache(a)%backward(tid)))call fftw_destroy_plan(cache(a)%backward(tid))
      enddo
    enddo
    call clear_source_types()
    call clear_return_types()
    do a=1,size(message_cache)
      if(message_cache(a)%datatype/=MPI_DATATYPE_NULL)call MPI_Type_free(message_cache(a)%datatype,e)
    enddo
    ! pairs counts pair/core tasks

    metrics%pairs=counts%pair_tiles
    lc=[metrics%pairs,metrics%evaluated,metrics%zero_overlap,metrics%fft_points, &
      metrics%source_payload,metrics%reverse_payload,metrics%max_fft_points]
    call MPI_Allreduce(lc,gc,6,MPI_INTEGER8,MPI_SUM,comm,e)
    call MPI_Allreduce(lc(7),gc(7),1,MPI_INTEGER8,MPI_MAX,comm,e)
    metrics%pairs=gc(1);metrics%evaluated=gc(2);metrics%zero_overlap=gc(3);metrics%fft_points=gc(4)
    metrics%source_payload=gc(5);metrics%reverse_payload=gc(6);metrics%max_fft_points=gc(7)
    dc=[counts%tiles,counts%pair_tiles,counts%plans,counts%source_messages,counts%return_messages, &
      counts%source_values,counts%return_values]
    call MPI_Allreduce(dc,dg,7,MPI_INTEGER8,MPI_SUM,comm,e)
    thread_stats=[requested,nteam,counts%omp_active]
    call MPI_Allreduce(thread_stats,all_thread_stats,3,MPI_INTEGER,MPI_MAX,comm,e)
    call MPI_Allreduce(counts%thread_workspace_bytes,all_thread_bytes,1,MPI_INTEGER8,MPI_MAX,comm,e)
    local_times=[counts%fft_setup_seconds,counts%pair_compute_seconds]
    call MPI_Allreduce(local_times,all_times,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,e)
    if(present(detail))then
      detail%fft_setup_seconds=all_times(1);detail%pair_compute_seconds=all_times(2)
      detail%omp_requested=all_thread_stats(1);detail%omp_effective=all_thread_stats(2)
      detail%omp_active=all_thread_stats(3);detail%thread_workspace_bytes=all_thread_bytes
      detail%tiles=dg(1);detail%pair_tiles=dg(2);detail%plans=dg(3)
      detail%source_messages=dg(4);detail%return_messages=dg(5)
      detail%source_values=dg(6);detail%return_values=dg(7)
    endif
    do j=1,nlocal
      do i=1,ng
        if(.not.ieee_is_finite(real(action(i,j),8)))bad=1
        if(.not.ieee_is_finite(aimag(action(i,j))))bad=1
      enddo
    enddo
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,comm,e)
    extra_counts=[counts%block_candidates,counts%block_visits,counts%fetch_requests, &
      counts%receive_hits,counts%extract_hits,counts%extracted_values,counts%halo_cache_bytes, &
      counts%return_workspace_bytes,counts%return_metadata_bytes]
    call MPI_Allreduce(extra_counts,global_counts,6,MPI_INTEGER8,MPI_SUM,comm,e)
    call MPI_Allreduce(extra_counts(7:9),global_counts(7:9),3,MPI_INTEGER8,MPI_MAX,comm,e)
    return_times=[counts%return_pack_seconds,counts%return_mpi_seconds,counts%return_apply_seconds]
    call MPI_Allreduce(return_times,all_return_times,3,MPI_DOUBLE_PRECISION,MPI_SUM,comm,e)
    return_stats=[counts%return_windows,counts%return_apply_parallel]
    call MPI_Allreduce(return_stats,all_return_stats,2,MPI_INTEGER8,MPI_SUM,comm,e)
    if(present(detail))then
      detail%control_skips=counts%control_skips
      detail%return_plan_hits=return_plan_hits;detail%return_plan_bytes=return_plan_bytes
      detail%return_pack_seconds=all_return_times(1)/np
      detail%return_mpi_seconds=all_return_times(2)/np
      detail%return_apply_seconds=all_return_times(3)/np
      detail%return_windows=all_return_stats(1);detail%return_apply_parallel=all_return_stats(2)
    endif
    pf_stats=[counts%prefetch_windows,counts%prefetch_fallbacks,counts%prefetch_parallel, &
      counts%prefetch_payload_bytes,counts%prefetch_metadata_bytes,counts%prefetch_active_threads]
    call MPI_Allreduce(pf_stats,pf_all_stats,3,MPI_INTEGER8,MPI_SUM,comm,e)
    call MPI_Allreduce(pf_stats(4:6),pf_all_stats(4:6),3,MPI_INTEGER8,MPI_MAX,comm,e)
    pf_times=[counts%prefetch_plan_seconds,counts%prefetch_extract_seconds, &
      counts%prefetch_pack_seconds,counts%prefetch_mpi_seconds]
    call MPI_Allreduce(pf_times,pf_all_times,4,MPI_DOUBLE_PRECISION,MPI_SUM,comm,e)
    if(present(detail))then
      detail%prefetch_windows=pf_all_stats(1);detail%prefetch_fallbacks=pf_all_stats(2)
      detail%prefetch_parallel=pf_all_stats(3);detail%prefetch_payload_bytes=pf_all_stats(4)
      detail%prefetch_metadata_bytes=pf_all_stats(5);detail%prefetch_active_threads=pf_all_stats(6)
      detail%prefetch_plan_seconds=pf_all_times(1)/np;detail%prefetch_extract_seconds=pf_all_times(2)/np
      detail%prefetch_pack_seconds=pf_all_times(3)/np;detail%prefetch_mpi_seconds=pf_all_times(4)/np
    endif
    phase(9)=MPI_Wtime()-tick
    phase(10)=MPI_Wtime()-begin_time
    phase(12)=counts%fft_setup_seconds
    ! Total ends before these profiling collectives and automatic teardown.
    ! 'other' exposes uninstrumented copies/teardown/control instead of hiding
    ! them. Rank maxima are not additive; the external caller measures full wall.
    phase(13)=max(0d0,phase(10)-sum(phase(1:9))-phase(11)-phase(12))
    call MPI_Allreduce(phase,extra_times,13,MPI_DOUBLE_PRECISION,MPI_MIN,comm,e)
    if(present(detail))detail%phase_min=extra_times
    call MPI_Allreduce(phase,extra_times,13,MPI_DOUBLE_PRECISION,MPI_MAX,comm,e)
    if(present(detail))detail%phase_max=extra_times
    call MPI_Allreduce(phase,extra_times,13,MPI_DOUBLE_PRECISION,MPI_SUM,comm,e)
    if(present(detail))then
      detail%phase_mean=extra_times/real(np,8)
      detail%block_candidates=global_counts(1);detail%block_visits=global_counts(2)
      detail%fetch_requests=global_counts(3);detail%receive_hits=global_counts(4)
      detail%extract_hits=global_counts(5);detail%extracted_values=global_counts(6)
      detail%halo_cache_bytes=global_counts(7)
      detail%return_workspace_bytes=global_counts(8)
      detail%return_metadata_bytes=global_counts(9)
    endif
  contains
    logical function pair_needed(i,j) result(yes)
      implicit none
      integer,intent(in) :: i,j
      yes=.false.
      if(i<j)return
      if(global_flags(2,i)==0.or.global_flags(2,j)==0)return
      if(global_flags(1,i)==0.and.global_flags(1,j)==0)return
      yes=.true.
    end function
    subroutine fetch_column(orb,columns,column)
      integer,intent(in) :: orb,column
      complex(8),allocatable,intent(inout) :: columns(:,:)
      integer :: cs,s
      real(8) :: cache_time
      cache_time=MPI_Wtime()
      if(prefetch_ready.and.task<=prefetch_until)then
        if(rank==worker)then
          counts%fetch_requests=counts%fetch_requests+1
          if(pf_local(orb)<=0)then
            ! A missing request is a planner bug, never silent uninitialized data.
            call MPI_Abort(comm,87,e)
          endif
          columns(:,column)=prefetch_values(:,pf_local(orb))
        endif
        phase(11)=phase(11)+MPI_Wtime()-cache_time
        return
      endif
      ! All ranks follow the same ordered requests, including nonparticipants.
      ! Thus cache hits skip the same send/receive without a new collective.
      cs=0
      do s=1,cache_slots
        if(receive_keys(s,worker+1)==orb)then
          cs=s;exit
        endif
      enddo
      if(rank==worker)counts%fetch_requests=counts%fetch_requests+1
      if(cs>0)then
        if(rank==worker)then
          columns(:,column)=received_halo(:,cs)
          if(owners(orb)/=worker)counts%receive_hits=counts%receive_hits+1
        endif
        phase(11)=phase(11)+MPI_Wtime()-cache_time
        return
      endif
      phase(11)=phase(11)+MPI_Wtime()-cache_time
      if(rank==worker)then
        call fetch(orb,columns(:,column))
      elseif(rank==owners(orb))then
        if(.not.allocated(owner_buffer))allocate(owner_buffer(max_fft))
        call fetch(orb,owner_buffer(:nfft))
      endif
      cache_time=MPI_Wtime()
      if(cache_slots>0)then
        cs=receive_next(worker+1)
        receive_keys(cs,worker+1)=orb
        receive_next(worker+1)=1+modulo(cs,cache_slots)
        if(rank==worker)then
          if(.not.allocated(received_halo))allocate(received_halo(nfft,cache_slots))
          received_halo(:,cs)=columns(:,column)
        endif
      endif
      phase(11)=phase(11)+MPI_Wtime()-cache_time
    end subroutine
    subroutine ensure_message_space(nrequests,noffsets)
      implicit none
      integer,intent(in) :: nrequests,noffsets
      if(allocated(requests))then
        if(size(requests)<nrequests)deallocate(requests,datatypes,retained_types)
      endif
      if(.not.allocated(requests))allocate(requests(nrequests),datatypes(nrequests),retained_types(nrequests))
      if(allocated(offsets))then
        if(size(offsets)<noffsets)deallocate(offsets)
      endif
      if(.not.allocated(offsets))allocate(offsets(noffsets))
    end subroutine ensure_message_space

    subroutine cached_message_type(ncols,block,datatype,retained)
      implicit none
      integer,intent(in) :: ncols,block
      integer,intent(out) :: datatype
      logical,intent(out) :: retained
      integer :: h,k,c,slot
      integer(int64) :: key,bytes
      key=modulo(int(block,int64),104729_int64)
      do c=1,ncols
        key=modulo(key*31+int(offsets(c),int64),104729_int64)
      enddo
      h=1+int(modulo(key,int(size(message_cache),int64)));slot=0;retained=.false.
      do k=1,8
        if(message_cache(h)%datatype==MPI_DATATYPE_NULL)then
          slot=h;exit
        endif
        if(message_cache(h)%block==block)then
          if(size(message_cache(h)%offsets)==ncols)then
            if(all(message_cache(h)%offsets==offsets(:ncols)))then
              datatype=message_cache(h)%datatype;retained=.true.
              return
            endif
          endif
        endif
        h=1+mod(h,size(message_cache))
      enddo
      call MPI_Type_create_hindexed_block(ncols,block,offsets,MPI_DOUBLE_COMPLEX,datatype,e)
      call MPI_Type_commit(datatype,e)
      bytes=int(ncols,int64)*int(storage_size(offsets)/8,int64)
      if(slot>0.and.message_cache_bytes+bytes<=8388608_int64)then
        allocate(message_cache(slot)%offsets(ncols))
        message_cache(slot)%offsets=offsets(:ncols)
        message_cache(slot)%block=block;message_cache(slot)%datatype=datatype
        message_cache_bytes=message_cache_bytes+bytes;retained=.true.
      endif
    end subroutine cached_message_type

    integer(int64) function current_input_bytes() result(bytes)
      bytes=0
      if(allocated(held_values))bytes=16_int64*size(held_values,kind=int64)
    end function current_input_bytes

    subroutine finish_source_requests()
      integer :: h
      real(8) :: t
      if(.not.source_inflight)return
      t=MPI_Wtime()
      if(source_nreq>0)call MPI_Waitall(source_nreq,source_requests,MPI_STATUSES_IGNORE,e)
      do h=1,source_nreq
        if(.not.source_retained(h))call MPI_Type_free(source_types(h),e)
      enddo
      source_inflight=.false.;source_nreq=0
      counts%prefetch_mpi_seconds=counts%prefetch_mpi_seconds+MPI_Wtime()-t
      phase(5)=phase(5)+MPI_Wtime()-t
    end subroutine finish_source_requests

    subroutine advance_input_window(start_hj,start_hi,before_task)
      integer,intent(in) :: start_hj,start_hi,before_task
      integer :: saved_until,next_task
      logical :: saved_ready
      if(future_valid)then
        call finish_source_requests()
        call move_alloc(future_values,prefetch_values)
        call move_alloc(future_local,pf_local)
        prefetch_until=future_until;prefetch_ready=.true.;future_valid=.false.
      else
        call prefetch_plan(start_hj,start_hi,before_task)
      endif
      if(.not.prefetch_ready.or..not.task_ready)return
      if(prefetch_until>=task_count)return
      saved_until=prefetch_until;saved_ready=prefetch_ready;next_task=saved_until+1
      call move_alloc(prefetch_values,held_values)
      call move_alloc(pf_local,held_local)
      call prefetch_plan(task_hj(next_task),task_hi(next_task),saved_until,.true.)
      future_valid=prefetch_ready;future_until=prefetch_until
      if(future_valid)then
        call move_alloc(prefetch_values,future_values)
        call move_alloc(pf_local,future_local)
      endif
      call move_alloc(held_values,prefetch_values)
      call move_alloc(held_local,pf_local)
      prefetch_until=saved_until;prefetch_ready=saved_ready
    end subroutine advance_input_window

    subroutine prefetch_plan(start_hj,start_hi,before_task,defer_wait)
      implicit none
      integer,intent(in) :: start_hj,start_hi,before_task
      logical,intent(in),optional :: defer_wait
      logical :: defer
      integer :: pj,pi,pjb,pib,pbi,pbj,pni,pnj,pfirst,plast,pcore,ptask,w,c,orb
      integer :: old_count,h,owner,g,previous,ngroups,ncols,s,needed_count,owned_count,max_group
      integer :: naccepted,idx,nreq,pass,req_index
      integer(int64) :: key
      logical :: ok,is_active,has_i,has_j,plan_hit
      integer :: plan_type,ss,local_columns,owned_columns
      integer(int64) :: plan_bytes
      real(8) :: t,tx
      defer=.false.;if(present(defer_wait))defer=defer_wait
      if(.not.allocated(pf_local))allocate(pf_local(norb))
      if(allocated(prefetch_values))deallocate(prefetch_values)
      prefetch_ready=.false.;prefetch_until=before_task
      if(prefetch_slots==0)return
      source_index=source_index+1
      plan_hit=source_replay.and.source_index<=source_windows
      t=MPI_Wtime();pf_count=0;pf_hash=0;pf_seen=0;pf_workers=0;pf_owners=0
      pf_local=0;pf_owned=0;ptask=before_task;naccepted=0
      if(plan_hit)then
        pf_count=source_plans(source_index)%count
        pf_orb(:pf_count)=source_plans(source_index)%orb
        pf_worker(:pf_count)=modulo(source_plans(source_index)%worker+modulo(tile-1,np),np)
        ptask=source_plans(source_index)%last_task;naccepted=ptask-before_task
      else
      ! Look ahead using the same sorted block lists and task progression.
      ! No orbital values, pair arithmetic or worker assignments change here.
      planning: do pj=start_hj,nhalo
        pjb=halo_blocks(pj);pbj=1+(pjb-1)*block;pnj=min(block,norb-pbj+1)
        if(block_core(pjb))then
          pfirst=pj;plast=nhalo
        else
          pfirst=1
          do while(pfirst<=ncoreblocks)
            if(core_blocks(pfirst)>=pjb)exit
            pfirst=pfirst+1
          enddo
          plast=ncoreblocks
        endif
        if(pj==start_hj)pfirst=start_hi
        do pi=pfirst,plast
          if(block_core(pjb))then
            pib=halo_blocks(pi)
          else
            pib=core_blocks(pi)
          endif
          pbi=1+(pib-1)*block;pni=min(block,norb-pbi+1)
          has_i=block_core(pib);has_j=block_core(pjb)
          w=modulo(modulo(tile-1,np)+modulo(ptask,np),np)
          old_count=pf_count;ok=.true.
          do c=1,pni+pnj
            if(c<=pni)then
              orb=pbi+c-1
              is_active=global_flags(2,orb)/=0.and.(has_j.or.global_flags(1,orb)/=0)
            else
              if(pbi==pbj)cycle
              orb=pbj+c-pni-1
              is_active=global_flags(2,orb)/=0.and.(has_i.or.global_flags(1,orb)/=0)
            endif
            if(.not.is_active)cycle
            key=int(w,int64)*int(norb,int64)+orb-1
            h=1+int(modulo(key,int(pf_hash_size,int64)))
            do
              s=pf_hash(h)
              if(s==0)exit
              if(pf_orb(s)==orb.and.pf_worker(s)==w)exit
              h=1+mod(h,pf_hash_size)
            enddo
            if(s/=0)cycle
            owner=owners(orb)
            if(pf_count==pf_capacity.or.pf_workers(w+1)==prefetch_slots)then
              ok=.false.;exit
            endif
            if(pf_seen(orb)==0.and.pf_owners(owner+1)==prefetch_slots)then
              ok=.false.;exit
            endif
            pf_count=pf_count+1;pf_orb(pf_count)=orb;pf_worker(pf_count)=w
            pf_hash(h)=pf_count;pf_hash_slot(pf_count)=h;pf_workers(w+1)=pf_workers(w+1)+1
            if(pf_seen(orb)==0)then
              pf_seen(orb)=pf_count;pf_owners(owner+1)=pf_owners(owner+1)+1
            endif
          enddo
          if(.not.ok)then
            ! Reverse insertion rollback retains valid open-address hash chains.
            do s=pf_count,old_count+1,-1
              orb=pf_orb(s);w=pf_worker(s);pf_workers(w+1)=pf_workers(w+1)-1
              pf_hash(pf_hash_slot(s))=0
              if(pf_seen(orb)==s)then
                pf_seen(orb)=0;pf_owners(owners(orb)+1)=pf_owners(owners(orb)+1)-1
              endif
            enddo
            pf_count=old_count
            exit planning
          endif
          ptask=ptask+1;naccepted=naccepted+1
          if(naccepted==4096)exit planning
        enddo
      enddo planning
      endif
      counts%prefetch_plan_seconds=counts%prefetch_plan_seconds+MPI_Wtime()-t
      phase(2)=phase(2)+MPI_Wtime()-t
      if(naccepted==0)then
        source_record=.false.;source_ready=.false.
        if(rank==0)counts%prefetch_fallbacks=counts%prefetch_fallbacks+1
        return
      endif
      prefetch_until=ptask;prefetch_ready=.true.
      if(rank==0)counts%prefetch_windows=counts%prefetch_windows+1
      t=MPI_Wtime();ngroups=0;pf_hash=0;pf_next=0;pf_thread_hits=0
      needed_count=0;owned_count=0
      if(plan_hit)then
        ngroups=source_plans(source_index)%groups
        pf_next(:pf_count)=source_plans(source_index)%next
        pf_sender(:ngroups)=source_plans(source_index)%sender
        pf_dest(:ngroups)=modulo(source_plans(source_index)%dest+modulo(tile-1,np),np)
        pf_head(:ngroups)=source_plans(source_index)%head
        pf_size(:ngroups)=source_plans(source_index)%ncols
      endif
      if(.not.allocated(source_recv_slot))allocate(source_recv_slot(pf_capacity),source_send_slot(pf_capacity))
      pf_workers=0;pf_owners=0;pf_seen=0
      do s=1,pf_count
        orb=pf_orb(s);w=pf_worker(s);owner=owners(orb)
        pf_workers(w+1)=pf_workers(w+1)+1;source_recv_slot(s)=pf_workers(w+1)
        if(pf_seen(orb)==0)then
          pf_owners(owner+1)=pf_owners(owner+1)+1;pf_seen(orb)=pf_owners(owner+1)
        endif
        source_send_slot(s)=pf_seen(orb)
        if(w==rank)then
          needed_count=needed_count+1;pf_local(orb)=needed_count
        endif
        if(owner==rank.and.pf_owned(orb)==0)then
          owned_count=owned_count+1;pf_owned(orb)=owned_count;pf_owned_orb(owned_count)=orb
        endif
        if(plan_hit)cycle
        key=int(owner,int64)*int(np,int64)+w
        h=1+int(modulo(key,int(pf_hash_size,int64)))
        do
          g=pf_hash(h)
          if(g==0)exit
          if(pf_sender(g)==owner.and.pf_dest(g)==w)exit
          h=1+mod(h,pf_hash_size)
        enddo
        if(g==0)then
          ngroups=ngroups+1;g=ngroups;pf_hash(h)=g
          pf_sender(g)=owner;pf_dest(g)=w;pf_head(g)=s;pf_tail(g)=0;pf_size(g)=0
        endif
        previous=pf_tail(g)
        if(previous/=0)pf_next(previous)=s
        pf_tail(g)=s;pf_size(g)=pf_size(g)+1
      enddo
      if(source_record.and.source_index<=size(source_plans))then
        plan_bytes=int(storage_size(0)/8,int64)*(5_int64*pf_count+6_int64*ngroups)
        if(source_bytes+plan_bytes<=33554432_int64)then
          source_plans(source_index)%count=pf_count;source_plans(source_index)%groups=ngroups
          source_plans(source_index)%last_task=ptask
          allocate(source_plans(source_index)%orb(pf_count),source_plans(source_index)%worker(pf_count), &
            source_plans(source_index)%next(pf_count),source_plans(source_index)%recv_slot(pf_count), &
            source_plans(source_index)%send_slot(pf_count),source_plans(source_index)%sender(ngroups), &
            source_plans(source_index)%dest(ngroups),source_plans(source_index)%head(ngroups), &
            source_plans(source_index)%ncols(ngroups),source_plans(source_index)%recv_type(ngroups), &
            source_plans(source_index)%send_type(ngroups))
          source_plans(source_index)%orb=pf_orb(:pf_count)
          source_plans(source_index)%worker=modulo(pf_worker(:pf_count)-modulo(tile-1,np),np)
          source_plans(source_index)%next=pf_next(:pf_count)
          source_plans(source_index)%recv_slot=source_recv_slot(:pf_count)
          source_plans(source_index)%send_slot=source_send_slot(:pf_count)
          source_plans(source_index)%sender=pf_sender(:ngroups)
          source_plans(source_index)%dest=modulo(pf_dest(:ngroups)-modulo(tile-1,np),np)
          source_plans(source_index)%head=pf_head(:ngroups);source_plans(source_index)%ncols=pf_size(:ngroups)
          source_plans(source_index)%recv_type=MPI_DATATYPE_NULL
          source_plans(source_index)%send_type=MPI_DATATYPE_NULL
          source_bytes=source_bytes+plan_bytes;source_windows=source_index
        else
          source_record=.false.;source_ready=.false.
        endif
      else if(source_record)then
        source_record=.false.;source_ready=.false.
      endif
      max_group=0
      do g=1,ngroups
        if(pf_sender(g)==pf_dest(g))cycle
        if(rank==pf_sender(g).or.rank==pf_dest(g))max_group=max(max_group,pf_size(g))
      enddo
allocate(prefetch_values(nfft,needed_count))
if(allocated(prefetch_owned))then
  if(size(prefetch_owned,2)<owned_count)deallocate(prefetch_owned)
endif
if(.not.allocated(prefetch_owned))allocate(prefetch_owned(nfft,owned_count))
if(allocated(prefetch_buffer))then
  if(size(prefetch_buffer,2)<max_group)deallocate(prefetch_buffer)
endif
if(.not.allocated(prefetch_buffer))allocate(prefetch_buffer(nfft,max_group))
      counts%prefetch_payload_bytes=max(counts%prefetch_payload_bytes, &
        16_int64*nfft*(int(needed_count,int64)+owned_count+max_group)+current_input_bytes())
      phase(3)=phase(3)+MPI_Wtime()-t;t=MPI_Wtime()
      if(int(nfft,int64)*owned_count>=32768_int64.and.nteam>1.and.owned_count>1)then
        counts%prefetch_parallel=counts%prefetch_parallel+1
!$omp parallel do default(none) schedule(static) num_threads(min(nteam,owned_count)) &
!$omp shared(owned_count) private(c)
        do c=1,owned_count
          call extract_prefetch_column(c)
        enddo
!$omp end parallel do
      else
        do c=1,owned_count
          call extract_prefetch_column(c)
        enddo
      endif
      counts%extracted_values=counts%extracted_values+int(nfft,int64)*owned_count
      counts%prefetch_extract_seconds=counts%prefetch_extract_seconds+MPI_Wtime()-t
      phase(4)=phase(4)+MPI_Wtime()-t
      ! Receive destinations are disjoint orbital columns. Owned values stay
      ! immutable until every request completes. Only the main thread uses MPI.
      call ensure_message_space(ngroups,maxval(pf_size(:ngroups)))
      nreq=0;tx=MPI_Wtime()
      do pass=1,2
        do g=1,ngroups
          owner=pf_sender(g);w=pf_dest(g);ncols=pf_size(g)
          if(owner==w)cycle
          if(pass==1.and.rank/=w)cycle
          if(pass==2.and.rank/=owner)cycle
          nreq=nreq+1
          plan_type=MPI_DATATYPE_NULL
          if(source_index<=source_windows)then
            if(allocated(source_plans(source_index)%recv_type))then
              if(pass==1)plan_type=source_plans(source_index)%recv_type(g)
              if(pass==2)plan_type=source_plans(source_index)%send_type(g)
            endif
          endif
          if(plan_type==MPI_DATATYPE_NULL)then
            s=pf_head(g)
            do c=1,ncols
              orb=pf_orb(s)
              if(pass==1)then
                idx=pf_local(orb)
              else
                idx=pf_owned(orb)
              endif
              offsets(c)=int(idx-1,MPI_ADDRESS_KIND)*int(nfft,MPI_ADDRESS_KIND)*16_MPI_ADDRESS_KIND
              s=pf_next(s)
            enddo
            if(source_index<=source_windows)then
              call MPI_Type_create_hindexed_block(ncols,nfft,offsets,MPI_DOUBLE_COMPLEX,plan_type,e)
              call MPI_Type_commit(plan_type,e)
              if(pass==1)source_plans(source_index)%recv_type(g)=plan_type
              if(pass==2)source_plans(source_index)%send_type(g)=plan_type
              datatypes(nreq)=plan_type;retained_types(nreq)=.true.
            else
              call cached_message_type(ncols,nfft,datatypes(nreq),retained_types(nreq))
            endif
          else
            datatypes(nreq)=plan_type;retained_types(nreq)=.true.
          endif
          if(pass==1)then
            call MPI_Irecv(prefetch_values,1,datatypes(nreq),owner,73,comm,requests(nreq),e)
          else
            call MPI_Isend(prefetch_owned,1,datatypes(nreq),w,73,comm,requests(nreq),e)
            counts%source_messages=counts%source_messages+1
            counts%source_values=counts%source_values+int(nfft,int64)*ncols
            metrics%source_payload=metrics%source_payload+16_int64*nfft*ncols
          endif
        enddo
      enddo
      ! Self-owned columns can be copied while remote requests are in flight.
!$omp parallel num_threads(nteam) default(shared) private(g,s,c,ncols)
      do g=1,ngroups
        if(pf_sender(g)/=rank.or.pf_dest(g)/=rank)cycle
        ncols=pf_size(g)
!$omp master
        s=pf_head(g)
        do c=1,ncols
          pf_list(c)=pf_orb(s);s=pf_next(s)
        enddo
!$omp end master
!$omp barrier
        call copy_prefetch_team(ncols,0)
      enddo
!$omp end parallel
      if(defer)then
        if(allocated(source_requests))deallocate(source_requests,source_types,source_retained)
        allocate(source_requests(nreq),source_types(nreq),source_retained(nreq))
        source_requests=requests(:nreq);source_types=datatypes(:nreq);source_retained=retained_types(:nreq)
        source_nreq=nreq;source_inflight=.true.
      else
      if(nreq>0)call MPI_Waitall(nreq,requests,MPI_STATUSES_IGNORE,e)
      do req_index=1,nreq
        if(.not.retained_types(req_index))call MPI_Type_free(datatypes(req_index),e)
      enddo
      endif
      counts%prefetch_mpi_seconds=counts%prefetch_mpi_seconds+MPI_Wtime()-tx
      phase(5)=phase(5)+MPI_Wtime()-tx
      counts%prefetch_active_threads=max(counts%prefetch_active_threads,int(count(pf_thread_hits>0),int64))
      ! Retain packing workspace for the next prefetch window of this tile.
    end subroutine prefetch_plan

    subroutine extract_prefetch_column(c)
      implicit none
      integer,intent(in) :: c
      integer :: ii,idx,thread
      thread=1
!$    thread=omp_get_thread_num()+1
      pf_thread_hits(thread)=1
      idx=local_index(pf_owned_orb(c))
      do ii=1,nfft
        prefetch_owned(ii,c)=value_at(patches(idx),rows(ii))
      enddo
    end subroutine extract_prefetch_column

    subroutine copy_prefetch_team(ncols,mode)
      implicit none
      integer,intent(in) :: ncols,mode
      integer :: c,point,orb
!$omp do collapse(2) schedule(static)
      do c=1,ncols
        do point=1,nfft
          orb=pf_list(c)
          select case(mode)
          case(0)
            prefetch_values(point,pf_local(orb))=prefetch_owned(point,pf_owned(orb))
          case(1)
            prefetch_buffer(point,c)=prefetch_owned(point,pf_owned(orb))
          case(2)
            prefetch_values(point,pf_local(orb))=prefetch_buffer(point,c)
          end select
        enddo
      enddo
!$omp end do
    end subroutine copy_prefetch_team

    subroutine copy_prefetch_columns(ncols,mode)
      implicit none
      integer,intent(in) :: ncols,mode
      integer :: c,point,orb
      if(int(nfft,int64)*ncols>=32768_int64.and.nteam>1)then
    !$omp parallel do collapse(2) default(none) schedule(static) &
    !$omp num_threads(min(nteam,nfft*ncols)) &
!$omp shared(ncols,mode,nfft,pf_list,pf_local,pf_owned,prefetch_values,prefetch_owned,prefetch_buffer) &
!$omp private(c,point,orb)
        do c=1,ncols
          do point=1,nfft
            orb=pf_list(c)
            select case(mode)
            case(0)
              prefetch_values(point,pf_local(orb))=prefetch_owned(point,pf_owned(orb))
            case(1)
              prefetch_buffer(point,c)=prefetch_owned(point,pf_owned(orb))
            case(2)
              prefetch_values(point,pf_local(orb))=prefetch_buffer(point,c)
            end select
          enddo
        enddo
    !$omp end parallel do
      else
        do c=1,ncols
          call copy_prefetch_column(c,mode)
        enddo
      endif
    end subroutine copy_prefetch_columns
    
    subroutine copy_prefetch_column(c,mode)
      implicit none
      integer,intent(in) :: c,mode
      integer :: orb,thread
      thread=1
!$    thread=omp_get_thread_num()+1
      pf_thread_hits(thread)=1
      orb=pf_list(c)
      select case(mode)
      case(0)
        prefetch_values(:,pf_local(orb))=prefetch_owned(:,pf_owned(orb))
      case(1)
        prefetch_buffer(:,c)=prefetch_owned(:,pf_owned(orb))
      case(2)
        prefetch_values(:,pf_local(orb))=prefetch_buffer(:,c)
      end select
    end subroutine copy_prefetch_column
    subroutine sum_column(out,parts)
      implicit none
      complex(8),intent(out) :: out(:)
      complex(8),intent(in) :: parts(:,:)
      integer :: t
      out=cmplx(0d0,0d0,8)
      do t=1,size(parts,2)
        out(:)=out(:)+parts(:,t)
      enddo
    end subroutine
    subroutine clear_source_types()
      integer :: p,g
      do p=1,size(source_plans)
        if(.not.allocated(source_plans(p)%recv_type))cycle
        do g=1,size(source_plans(p)%recv_type)
          if(source_plans(p)%recv_type(g)/=MPI_DATATYPE_NULL)call MPI_Type_free(source_plans(p)%recv_type(g),e)
          if(source_plans(p)%send_type(g)/=MPI_DATATYPE_NULL)call MPI_Type_free(source_plans(p)%send_type(g),e)
        enddo
      enddo
    end subroutine clear_source_types

    subroutine begin_source_plans()
      logical :: same
      same=.false.
      if(allocated(source_flags))then
        if(source_nfft==nfft.and.source_block==block.and.source_slots==prefetch_slots)then
          same=all(source_flags==global_flags)
        endif
      endif
      source_index=0;source_replay=same.and.source_ready
      source_record=.not.same.and.prefetch_slots>0
      if(.not.same)then
        call clear_source_types()
        source_plans=source_window_plan();source_bytes=0;source_windows=0;source_ready=.false.
        if(.not.allocated(source_flags))allocate(source_flags(2,norb))
        source_flags=global_flags;source_nfft=nfft;source_block=block;source_slots=prefetch_slots
      endif
    end subroutine begin_source_plans

    subroutine clear_return_types()
      integer :: p,g
      do p=1,size(return_plans)
        if(.not.allocated(return_plans(p)%datatype))cycle
        do g=1,size(return_plans(p)%datatype)
          if(return_plans(p)%own_type(g))call MPI_Type_free(return_plans(p)%datatype(g),e)
        enddo
      enddo
    end subroutine clear_return_types

    subroutine return_message_type(g,datatype,retained,ncols)
      integer,intent(in) :: g
      integer,intent(out) :: datatype,ncols
      logical,intent(out) :: retained
      integer :: s
      logical :: planned
      planned=.false.
      if(return_plan_replay.or.return_plan_record)then
        if(return_plan_index<=size(return_plans))then
          planned=allocated(return_plans(return_plan_index)%datatype)
        endif
      endif
      if(planned)then
        datatype=return_plans(return_plan_index)%datatype(g)
        if(datatype/=MPI_DATATYPE_NULL)then
          ncols=return_plans(return_plan_index)%columns(g);retained=.true.
          return
        endif
      endif
      ncols=0;s=group_head(g)
      do while(s/=0)
        ncols=ncols+1
        offsets(ncols)=int(s-1,MPI_ADDRESS_KIND)*int(ncore,MPI_ADDRESS_KIND)*16_MPI_ADDRESS_KIND
        s=group_next(s)
      enddo
      call cached_message_type(ncols,ncore,datatype,retained)
      if(planned)then
        return_plans(return_plan_index)%datatype(g)=datatype
        return_plans(return_plan_index)%columns(g)=ncols
        return_plans(return_plan_index)%own_type(g)=.not.retained
        retained=.true.
      endif
    end subroutine return_message_type

    subroutine reset_task_scan()
      scan_hj=1;scan_hi=0;scan_core=1;scan_first=1;scan_last=0
    end subroutine reset_task_scan

    subroutine scan_pair_task(found)
      logical,intent(out) :: found
      integer :: i,j,ii,jj,ib,jb,offset
      found=.false.
      do while(scan_hj<=nhalo)
        jb=halo_blocks(scan_hj);bj=1+(jb-1)*block;nj=min(block,norb-bj+1)
        if(scan_hi==0)then
          do while(scan_core<=ncoreblocks)
            if(core_blocks(scan_core)>=jb)exit
            scan_core=scan_core+1
          enddo
          if(block_core(jb))then
            scan_first=scan_hj;scan_last=nhalo
          else
            scan_first=scan_core;scan_last=ncoreblocks
          endif
          scan_hi=scan_first
        endif
        if(scan_hi>scan_last)then
          scan_hj=scan_hj+1;scan_hi=0;cycle
        endif
        hj=scan_hj;hi=scan_hi;scan_hi=scan_hi+1
        if(block_core(jb))then
          ib=halo_blocks(hi)
        else
          ib=core_blocks(hi)
        endif
        bi=1+(ib-1)*block;ni=min(block,norb-bi+1)
        if(rank==0)counts%block_visits=counts%block_visits+1
        active_i=.false.;active_j=.false.
        if(allocated(halo_mask))then
          call active_masks(halo_mask(ib),core_mask(ib),halo_mask(jb),core_mask(jb),bi==bj,left_mask,right_mask)
          if(left_mask==0)cycle
          do ii=1,ni
            active_i(ii)=btest(left_mask,ii-1)
          enddo
          do jj=1,nj
            active_j(jj)=btest(right_mask,jj-1)
          enddo
        else
          do jj=1,nj
            j=bj+jj-1
            do ii=1,ni
              i=bi+ii-1
              if(.not.pair_needed(i,j))cycle
              active_i(ii)=.true.;active_j(jj)=.true.
            enddo
          enddo
          if(.not.any(active_i))cycle
        endif
        found=.true.;return
      enddo
    end subroutine scan_pair_task

    subroutine begin_task_plan()
      integer :: capacity,n,entries
      logical :: same,found
      same=.false.
      if(allocated(task_flags))then
        if(task_block==block)same=all(task_flags==global_flags)
      endif
      task_pos=0
      call reset_task_scan()
      if(same)return
      task_ready=.false.;task_count=0
      if(allocated(task_hj))deallocate(task_hj,task_hi,task_bi,task_bj,task_left,task_right,task_prefix)
      if(.not.allocated(task_flags))allocate(task_flags(2,norb))
      task_flags=global_flags;task_block=block
      if(.not.allocated(halo_mask))return
      capacity=int(min(838860_int64,int(nhalo,int64)*(int(nhalo,int64)+1)/2))
      allocate(task_hj(capacity),task_hi(capacity),task_bi(capacity),task_bj(capacity), &
        task_left(capacity),task_right(capacity),task_prefix(0:capacity))
      task_prefix(0)=0
      do
        call scan_pair_task(found)
        if(.not.found)exit
        if(task_count==capacity)then
          deallocate(task_hj,task_hi,task_bi,task_bj,task_left,task_right,task_prefix)
          task_count=0;call reset_task_scan();return
        endif
        task_count=task_count+1;n=task_count
        task_hj(n)=hj;task_hi(n)=hi;task_bi(n)=bi;task_bj(n)=bj
        task_left(n)=left_mask;task_right(n)=right_mask
        if(bi==bj)then
          entries=count(active_i(:ni).or.active_j(:nj))
        else
          entries=count(active_i(:ni))+count(active_j(:nj))
        endif
        task_prefix(n)=task_prefix(n-1)+entries
      enddo
      task_ready=.true.
      call reset_task_scan()
    end subroutine begin_task_plan

    subroutine next_pair_task(found)
      logical,intent(out) :: found
      integer :: ii,jj
      if(.not.task_ready)then
        call scan_pair_task(found);return
      endif
      task_pos=task_pos+1;found=task_pos<=task_count
      if(.not.found)return
      hj=task_hj(task_pos);hi=task_hi(task_pos);bi=task_bi(task_pos);bj=task_bj(task_pos)
      ni=min(block,norb-bi+1);nj=min(block,norb-bj+1)
      active_i=.false.;active_j=.false.
      do ii=1,ni
        active_i(ii)=btest(task_left(task_pos),ii-1)
      enddo
      do jj=1,nj
        active_j(jj)=btest(task_right(task_pos),jj-1)
      enddo
    end subroutine next_pair_task

    subroutine begin_return_plan()
      implicit none
      logical :: same
      same=.false.
      if(allocated(return_plan_flags))then
        if(return_plan_core==ncore.and.return_plan_block==block.and.return_plan_capacity==queue_capacity)then
          same=all(return_plan_flags==global_flags)
        endif
      endif
      return_plan_replay=same.and.return_plan_ready;return_plan_record=.not.same.and..not.return_direct
      return_plan_index=0
      if(return_plan_record)then
        call clear_return_types()
        return_plans=return_queue_plan()
        return_plan_ready=.false.;return_plan_bytes=0;return_plan_windows=0
        if(.not.allocated(return_plan_flags))allocate(return_plan_flags(2,norb))
        return_plan_flags=global_flags
        return_plan_core=ncore;return_plan_block=block;return_plan_capacity=queue_capacity
      endif
    end subroutine begin_return_plan

    subroutine prepare_return_window()
      implicit none
      integer :: g,n,ng
      integer(int64) :: bytes
      return_plan_index=return_plan_index+1
      if(return_plan_replay)then
        n=return_plans(return_plan_index)%count;ng=return_plans(return_plan_index)%groups
        group_count=ng
        queue_next(:n)=return_plans(return_plan_index)%next
        group_next(:n)=return_plans(return_plan_index)%groupnext
        orb_head=return_plans(return_plan_index)%head
        group_owner(:ng)=return_plans(return_plan_index)%owner
        group_head(:ng)=return_plans(return_plan_index)%grouphead
        do g=1,ng
          group_worker(g)=modulo(return_plans(return_plan_index)%worker(g)+modulo(tile-1,np),np)
        enddo
        do g=1,ng
          if(rank/=group_worker(g).and.rank/=group_owner(g))cycle
          if(.not.allocated(queued_values))allocate(queued_values(ncore,queue_capacity))
          exit
        enddo
        counts%return_workspace_bytes=max(counts%return_workspace_bytes,16_int64*ncore*queue_capacity)
        return_plan_hits=return_plan_hits+1
      else if(return_plan_record)then
        n=queue_count;ng=group_count
        bytes=int(storage_size(0)/8,int64)*(2_int64*n+6_int64*ng+norb)
        if(return_plan_index>size(return_plans).or.return_plan_bytes+bytes>33554432_int64)then
          return_plan_record=.false.;return_plan_ready=.false.
          call clear_return_types()
        return_plans=return_queue_plan();return_plan_bytes=0
          return
        endif
        return_plans(return_plan_index)%count=n;return_plans(return_plan_index)%groups=ng
        allocate(return_plans(return_plan_index)%next(n),return_plans(return_plan_index)%groupnext(n), &
          return_plans(return_plan_index)%head(norb),return_plans(return_plan_index)%worker(ng), &
          return_plans(return_plan_index)%owner(ng),return_plans(return_plan_index)%grouphead(ng), &
          return_plans(return_plan_index)%datatype(ng),return_plans(return_plan_index)%columns(ng), &
          return_plans(return_plan_index)%own_type(ng))
        return_plans(return_plan_index)%datatype=MPI_DATATYPE_NULL
        return_plans(return_plan_index)%columns=0
        return_plans(return_plan_index)%own_type=.false.
        return_plans(return_plan_index)%next=queue_next(:n)
        return_plans(return_plan_index)%groupnext=group_next(:n)
        return_plans(return_plan_index)%head=orb_head
        return_plans(return_plan_index)%owner=group_owner(:ng)
        return_plans(return_plan_index)%grouphead=group_head(:ng)
        do g=1,ng
          return_plans(return_plan_index)%worker(g)=modulo(group_worker(g)-modulo(tile-1,np),np)
        enddo
        return_plan_bytes=return_plan_bytes+bytes
      endif
    end subroutine prepare_return_window

    subroutine advance_return_queue(entries)
      implicit none
      integer,intent(in) :: entries
      integer :: remaining,step
      real(8) :: t,paused,start
      t=MPI_Wtime();paused=0d0;remaining=entries
      do while(remaining>0)
        if(queue_count==queue_capacity)then
          start=MPI_Wtime();call flush_returns();paused=paused+MPI_Wtime()-start
        endif
        step=min(remaining,queue_capacity-queue_count)
        queue_count=queue_count+step;remaining=remaining-step
      enddo
      counts%return_pack_seconds=counts%return_pack_seconds+MPI_Wtime()-t-paused
      phase(6)=phase(6)+MPI_Wtime()-t
    end subroutine advance_return_queue

    subroutine queue_block(first,columns,active)
      implicit none
      integer,intent(in) :: first
      complex(8),allocatable,intent(in) :: columns(:,:)
      logical,intent(in) :: active(:)
      integer :: c,orb,owner,h,g,previous
      integer(int64) :: key
      real(8) :: t,paused,flush_start
      if(return_direct)then
        call return_immediate(first,columns,active)
        return
      endif
      t=MPI_Wtime();paused=0d0
      if(return_plan_replay)then
        do c=1,size(active)
          if(.not.active(c))cycle
          if(queue_count==queue_capacity)then
            flush_start=MPI_Wtime();call flush_returns();paused=paused+MPI_Wtime()-flush_start
          endif
          queue_count=queue_count+1
          if(rank==worker)then
            if(.not.allocated(queued_values))allocate(queued_values(ncore,queue_capacity))
            queued_values(:,queue_count)=columns(:,c)
          endif
        enddo
        counts%return_pack_seconds=counts%return_pack_seconds+MPI_Wtime()-t-paused
        phase(6)=phase(6)+MPI_Wtime()-t
        return
      endif
      do c=1,size(active)
        if(.not.active(c))cycle
        if(queue_count==queue_capacity)then
          flush_start=MPI_Wtime();call flush_returns();paused=paused+MPI_Wtime()-flush_start
        endif
        orb=first+c-1;owner=owners(orb)
        ! All ranks append exactly the same metadata. Include local returns:
        ! applying them early would reorder them relative to delayed messages.
        queue_count=queue_count+1;queue_orb(queue_count)=orb;queue_next(queue_count)=0
        previous=orb_tail(orb)
        if(previous==0)then
          orb_head(orb)=queue_count
        else
          queue_next(previous)=queue_count
        endif
        orb_tail(orb)=queue_count
        key=int(worker,int64)*int(np,int64)+owner
        h=1+int(modulo(key,int(hash_size,int64)))
        do
          g=group_hash(h)
          if(g==0)exit
          if(group_worker(g)==worker.and.group_owner(g)==owner)exit
          h=1+mod(h,hash_size)
        enddo
        if(g==0)then
          group_count=group_count+1;g=group_count;group_hash(h)=g
          group_worker(g)=worker;group_owner(g)=owner;group_head(g)=queue_count;group_tail(g)=0
        endif
        previous=group_tail(g)
        if(previous/=0)group_next(previous)=queue_count
        group_tail(g)=queue_count;group_next(queue_count)=0
        if(rank==worker.or.rank==owner)then
          if(.not.allocated(queued_values))allocate(queued_values(ncore,queue_capacity))
          counts%return_workspace_bytes=max(counts%return_workspace_bytes, &
            16_int64*ncore*queue_capacity)
        endif
        if(rank==worker)queued_values(:,queue_count)=columns(:,c)
      enddo
      counts%return_pack_seconds=counts%return_pack_seconds+MPI_Wtime()-t-paused
      phase(6)=phase(6)+MPI_Wtime()-t
    end subroutine queue_block

    subroutine flush_returns()
      implicit none
      integer :: g,sender,owner,ncols,h,s,idx,j,orb,point,nreq,pass,bundle,nbundles
      real(8) :: t
      if(queue_count==0)return
      call finish_pending_returns()
      t=MPI_Wtime();call prepare_return_window()
      counts%return_pack_seconds=counts%return_pack_seconds+MPI_Wtime()-t
      if(rank==0)counts%return_windows=counts%return_windows+1
      ! Queue slots belong to exactly one sender/owner pair. Receive slots
      ! never overlap outgoing slots; all transfers finish before ordered sums.
      call ensure_message_space(group_count,queue_count)
      if(allocated(pending_requests))then
        if(size(pending_requests)<group_count)deallocate(pending_requests,pending_types,pending_retained)
      endif
      if(.not.allocated(pending_requests))then
        allocate(pending_requests(group_count),pending_types(group_count),pending_retained(group_count))
      endif
      if(allocated(pending_next))then
        if(size(pending_next)<queue_count)deallocate(pending_next)
      endif
      if(.not.allocated(pending_next))allocate(pending_next(queue_capacity))
      if(.not.allocated(pending_head))allocate(pending_head(norb))
      pending_next(:queue_count)=queue_next(:queue_count);pending_head=orb_head
      call move_alloc(pending_values,spare_values)
      call move_alloc(queued_values,pending_values)
      call move_alloc(spare_values,queued_values)
      counts%return_workspace_bytes=max(counts%return_workspace_bytes,32_int64*ncore*queue_capacity)
      nreq=0;t=MPI_Wtime()
      do pass=1,2
        do g=1,group_count
          sender=group_worker(g);owner=group_owner(g)
          if(sender==owner)cycle
          if(pass==1.and.rank/=owner)cycle
          if(pass==2.and.rank/=sender)cycle
          nreq=nreq+1
          call return_message_type(g,pending_types(nreq),pending_retained(nreq),ncols)
          if(pass==1)then
            call MPI_Irecv(pending_values,1,pending_types(nreq),sender,72,comm,pending_requests(nreq),e)
          else
            call MPI_Isend(pending_values,1,pending_types(nreq),owner,72,comm,pending_requests(nreq),e)
            metrics%reverse_payload=metrics%reverse_payload+16_int64*ncore*ncols
            counts%return_messages=counts%return_messages+1
            counts%return_values=counts%return_values+int(ncore,int64)*ncols
          endif
        enddo
      enddo
      pending_nreq=nreq;pending_count=queue_count;pending_active=.true.
      counts%return_mpi_seconds=counts%return_mpi_seconds+MPI_Wtime()-t
      queue_count=0;group_count=0;group_hash=0;orb_head=0;orb_tail=0
    end subroutine flush_returns

    subroutine finish_pending_returns()
      implicit none
      integer :: h,j,bundle,nbundles
      real(8) :: t
      if(.not.pending_active)return
      t=MPI_Wtime()
      if(pending_nreq>0)call MPI_Waitall(pending_nreq,pending_requests,MPI_STATUSES_IGNORE,e)
      do h=1,pending_nreq
        if(.not.pending_retained(h))call MPI_Type_free(pending_types(h),e)
      enddo
      counts%return_mpi_seconds=counts%return_mpi_seconds+MPI_Wtime()-t
      t=MPI_Wtime()
      if(int(ncore,int64)*pending_count>=32768_int64.and.nteam>1.and.nlocal>0)then
        counts%return_apply_parallel=counts%return_apply_parallel+1
        nbundles=1+(ncore-1)/8
!$omp parallel do default(none) schedule(static) collapse(2) num_threads(min(nteam,nbundles*nlocal)) &
!$omp shared(nlocal,nbundles) private(j,bundle)
        do j=1,nlocal
          do bundle=1,nbundles
            call apply_pending_bundle(j,bundle)
          enddo
        enddo
!$omp end parallel do
      else
        do j=1,nlocal
          call apply_pending_orbital(j)
        enddo
      endif
      counts%return_apply_seconds=counts%return_apply_seconds+MPI_Wtime()-t
      pending_active=.false.
    end subroutine finish_pending_returns

    subroutine apply_bundle(idx,bundle)
      implicit none
      integer,intent(in) :: idx,bundle
      integer :: first,last,point,s,k
      complex(8) :: sums(8)
      first=1+8*(bundle-1);last=min(ncore,first+7)
      do point=first,last
        sums(point-first+1)=action(core_rows(point),idx)
      enddo
      s=orb_head(ids(idx))
      do while(s/=0)
!$omp simd private(k)
        do point=first,last
          k=point-first+1
          sums(k)=sums(k)+queued_values(point,s)
        enddo
        s=queue_next(s)
      enddo
      do point=first,last
        action(core_rows(point),idx)=sums(point-first+1)
      enddo
    end subroutine apply_bundle

    subroutine apply_point(idx,point)
      integer,intent(in) :: idx,point
      integer :: slot
      slot=orb_head(ids(idx))
      do while(slot/=0)
        action(core_rows(point),idx)=action(core_rows(point),idx)+queued_values(point,slot)
        slot=queue_next(slot)
      enddo
    end subroutine apply_point
    subroutine apply_orbital(idx)
      implicit none
      integer,intent(in) :: idx
      integer :: s,ii
      s=orb_head(ids(idx))
      do while(s/=0)
        do ii=1,ncore
          action(core_rows(ii),idx)=action(core_rows(ii),idx)+queued_values(ii,s)
        enddo
        s=queue_next(s)
      enddo
    end subroutine apply_orbital
    subroutine apply_pending_bundle(idx,bundle)
      implicit none
      integer,intent(in) :: idx,bundle
      integer :: first,last,point,s,k
      complex(8) :: sums(8)
      first=1+8*(bundle-1);last=min(ncore,first+7)
      do point=first,last
        sums(point-first+1)=action(core_rows(point),idx)
      enddo
      s=pending_head(ids(idx))
      do while(s/=0)
!$omp simd private(k)
        do point=first,last
          k=point-first+1
          sums(k)=sums(k)+pending_values(point,s)
        enddo
        s=pending_next(s)
      enddo
      do point=first,last
        action(core_rows(point),idx)=sums(point-first+1)
      enddo
    end subroutine apply_pending_bundle

    subroutine apply_pending_point(idx,point)
      integer,intent(in) :: idx,point
      integer :: slot
      slot=pending_head(ids(idx))
      do while(slot/=0)
        action(core_rows(point),idx)=action(core_rows(point),idx)+pending_values(point,slot)
        slot=pending_next(slot)
      enddo
    end subroutine apply_pending_point
    subroutine apply_pending_orbital(idx)
      implicit none
      integer,intent(in) :: idx
      integer :: s,ii
      s=pending_head(ids(idx))
      do while(s/=0)
        do ii=1,ncore
          action(core_rows(ii),idx)=action(core_rows(ii),idx)+pending_values(ii,s)
        enddo
        s=pending_next(s)
      enddo
    end subroutine apply_pending_orbital
    subroutine return_immediate(first,columns,active)
      implicit none
      integer,intent(in) :: first
      complex(8),allocatable,intent(in) :: columns(:,:)
      logical,intent(in) :: active(:)
      integer :: destinations(size(active)),indices(size(active)),ndest,ncol,d,c,h,owner,idx,ii
      logical :: seen
      real(8) :: t
      t=MPI_Wtime();ndest=0
      ! Group only this block. Distinct columns update distinct orbitals, so
      ! changing their communication order cannot change any element's sum.
      do c=1,size(active)
        if(.not.active(c))cycle
        owner=owners(first+c-1);seen=.false.
        do d=1,ndest
          if(destinations(d)==owner)seen=.true.
        enddo
        if(seen)cycle
        ndest=ndest+1;destinations(ndest)=owner
      enddo
      do d=1,ndest
        owner=destinations(d)
        if(rank/=worker.and.rank/=owner)cycle
        ncol=0
        do c=1,size(active)
          if(.not.active(c))cycle
          if(owners(first+c-1)/=owner)cycle
          ncol=ncol+1;indices(ncol)=c
        enddo
        if(owner/=worker)then
          ! block <= 64MiB/(32*nfft), ncore <= nfft: payload <= 32MiB.
          if(.not.allocated(return_buffer))then
            allocate(return_buffer(ncore,min(block,norb)))
            counts%return_workspace_bytes=max(counts%return_workspace_bytes,16_int64*size(return_buffer,kind=int64))
          endif
          if(rank==worker)then
            do h=1,ncol
              return_buffer(:,h)=columns(:,indices(h))
            enddo
            call MPI_Send(return_buffer,ncore*ncol,MPI_DOUBLE_COMPLEX,owner,72,comm,e)
            metrics%reverse_payload=metrics%reverse_payload+16_int64*ncore*ncol
            counts%return_messages=counts%return_messages+1
            counts%return_values=counts%return_values+int(ncore,int64)*ncol
          else
            call MPI_Recv(return_buffer,ncore*ncol,MPI_DOUBLE_COMPLEX,worker,72,comm,MPI_STATUS_IGNORE,e)
          endif
        endif
        if(rank==owner)then
!$omp parallel do default(none) schedule(static) num_threads(min(nteam,ncol)) &
!$omp if(int(ncore,int64)*int(ncol,int64)>=32768_int64.and.ncol>1.and.nteam>1) &
!$omp shared(ncol,indices,first,local_index,ncore,core_rows,action,columns,return_buffer,owner,worker) &
!$omp private(h,c,idx,ii)
          do h=1,ncol
            c=indices(h);idx=local_index(first+c-1)
            if(owner==worker)then
              do ii=1,ncore
                action(core_rows(ii),idx)=action(core_rows(ii),idx)+columns(ii,c)
              enddo
            else
              do ii=1,ncore
                action(core_rows(ii),idx)=action(core_rows(ii),idx)+return_buffer(ii,h)
              enddo
            endif
          enddo
!$omp end parallel do
        endif
      enddo
      phase(6)=phase(6)+MPI_Wtime()-t
    end subroutine
    subroutine fetch(orb,values)
      implicit none
      integer,intent(in) :: orb
      complex(8),intent(out) :: values(:)
      integer :: ii,owner,idx,cs,s
      real(8) :: t
      owner=owners(orb);idx=local_index(orb)
      if(rank==owner)then
        t=MPI_Wtime();cs=0
        do s=1,extract_slots
          if(extract_keys(s)==orb)then
            cs=s;exit
          endif
        enddo
        if(cs>0)then
          values=extracted_halo(:,cs)
          counts%extract_hits=counts%extract_hits+1
        else
          ! MPI is outside the region. The independent halo rows are the
          ! outermost available work in a single on-demand column extraction.
!$omp parallel do default(none) schedule(static) num_threads(nteam) if(nfft>=4096.and.nteam>1) &
!$omp shared(values,patches,idx,rows,nfft) private(ii)
          do ii=1,nfft
            values(ii)=value_at(patches(idx),rows(ii))
          enddo
!$omp end parallel do
          counts%extracted_values=counts%extracted_values+int(nfft,int64)
          if(extract_slots>0)then
            cs=extract_next;extract_keys(cs)=orb;extract_next=1+modulo(cs,extract_slots)
            if(.not.allocated(extracted_halo))allocate(extracted_halo(nfft,extract_slots))
            extracted_halo(:,cs)=values
          endif
        endif
        phase(4)=phase(4)+MPI_Wtime()-t
        if(owner/=worker)then
          t=MPI_Wtime()
          call MPI_Send(values,nfft,MPI_DOUBLE_COMPLEX,worker,71,comm,e)
          phase(5)=phase(5)+MPI_Wtime()-t
          metrics%source_payload=metrics%source_payload+16_int64*nfft
          counts%source_messages=counts%source_messages+1
          counts%source_values=counts%source_values+int(nfft,int64)
        endif
      elseif(rank==worker)then
        t=MPI_Wtime()
        call MPI_Recv(values,nfft,MPI_DOUBLE_COMPLEX,owner,71,comm,MPI_STATUS_IGNORE,e)
        phase(5)=phase(5)+MPI_Wtime()-t
      endif
    end subroutine

  end subroutine

  subroutine support_bounds(patch,shape,bounds)
    implicit none
    type(orbital_patch),intent(in) :: patch
    integer,intent(in) :: shape(3)
    integer,intent(out) :: bounds(6)
    integer :: axis,k,g,v(3),t,first,previous,last,gap,maxgap,start
    logical,allocatable :: occupied(:)
    bounds=0
    if(size(patch%row)==0)return
    do axis=1,3
      allocate(occupied(shape(axis)));occupied=.false.
      do k=1,size(patch%row)
        g=patch%row(k)-1
        v=[modulo(g,shape(1)),modulo(g/shape(1),shape(2)),g/(shape(1)*shape(2))]
        occupied(v(axis)+1)=.true.
      enddo
      first=-1;previous=-1;last=-1;maxgap=-1;start=0
      do t=0,shape(axis)-1
        if(.not.occupied(t+1))cycle
        if(first<0)first=t
        if(previous>=0)then
          gap=t-previous
          if(gap>maxgap)then
            maxgap=gap;start=t
          endif
        endif
        previous=t;last=t
      enddo
      gap=first+shape(axis)-last
      if(gap>maxgap)then
        maxgap=gap;start=first
      endif
      bounds(axis)=start;bounds(axis+3)=shape(axis)-maxgap+1
      deallocate(occupied)
    enddo
  end subroutine
  logical function boxes_intersect(a,la,b,lb,shape) result(overlap)
    implicit none
    integer,intent(in) :: a(3),la(3),b(3),lb(3),shape(3)
    integer :: axis
    overlap=.false.
    if(any(lb==0))return
    do axis=1,3
      if(modulo(a(axis)-b(axis),shape(axis))>=lb(axis).and. &
         modulo(b(axis)-a(axis),shape(axis))>=la(axis))return
    enddo
    overlap=.true.
  end function
  complex(8) function value_at(patch,row) result(value)
    implicit none
    type(orbital_patch),intent(in) :: patch
    integer,intent(in) :: row
    integer :: lo,hi,mid
    lo=1;hi=size(patch%row);value=0
    do while(lo<=hi)
      mid=lo+(hi-lo)/2
      if(patch%row(mid)==row)then
        value=patch%value(mid);return
      elseif(patch%row(mid)<row)then
        lo=mid+1
      else
        hi=mid-1
      endif
    enddo
  end function
  logical function nonzero(value)
    implicit none
    complex(8),intent(in) :: value
    nonzero=real(value,8)/=0d0.or.aimag(value)/=0d0
  end function
  logical function any_nonzero_product(a,b) result(found)
    implicit none
    complex(8),intent(in) :: a(:),b(:)
    integer :: i
    found=.false.
    do i=1,size(a)
      if(nonzero(a(i)).and.nonzero(b(i)))then
        found=.true.;return
      endif
    enddo
  end function
  subroutine prepare_fft(c,length,shape,kernel,work,spectrum,bad)
    implicit none
    type(fft_cache),intent(out) :: c
    integer,intent(in) :: length(3),shape(3)
    real(8),intent(in) :: kernel(:)
    integer,intent(inout) :: bad
    complex(8),intent(inout) :: work(:,:),spectrum(:,:)
    integer :: k,n,v(3),d(3),a,q,tid,nteam
    c%shape=length;n=product(length);nteam=size(work,2)
    allocate(c%kernel(n),c%forward(nteam),c%backward(nteam))
    c%forward=c_null_ptr;c%backward=c_null_ptr
    do k=0,n-1
      v=[modulo(k,length(1)),modulo(k/length(1),length(2)),k/(length(1)*length(2))]
      d=v
      do a=1,3
        if(v(a)>length(a)/2)d(a)=v(a)-length(a)
      enddo
      d=modulo(d,shape);q=1+d(1)+shape(1)*(d(2)+shape(2)*d(3))
      work(k+1,1)=cmplx(kernel(q),0d0,8)
    enddo
    ! Planning/destruction remain serial; each thread executes on its own arrays.
    do tid=1,nteam
      c%forward(tid)=fftw_plan_dft_3d(length(3),length(2),length(1), &
        work(:,tid),spectrum(:,tid),FFTW_FORWARD,FFTW_ESTIMATE)
      c%backward(tid)=fftw_plan_dft_3d(length(3),length(2),length(1), &
        spectrum(:,tid),work(:,tid),FFTW_BACKWARD,FFTW_ESTIMATE)
      if(.not.c_associated(c%forward(tid)).or..not.c_associated(c%backward(tid)))bad=1
    enddo
    if(bad/=0)return
    call fftw_execute_dft(c%forward(1),work(:,1),spectrum(:,1))
    c%kernel=spectrum(:n,1)
  end subroutine
end module
