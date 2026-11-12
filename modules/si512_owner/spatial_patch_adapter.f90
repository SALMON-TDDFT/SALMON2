module spatial_patch_adapter
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use iso_fortran_env, only: int64
  use orbital_owner_layout, only: owner_layout,layout_prepare,layout_forward,layout_reverse
  use wannier_patch_owner, only: orbital_patch,patch_metrics,old_patch_apply=>patch_apply,pack_exact_patch
  use wannier_neighbor_owner, only: patch_apply,neighbor_metrics,neighbor_phase_names
  implicit none
  private
  public :: spatial_owner_action
contains
  subroutine spatial_owner_action(n,h,m,lo,comm_r,comm_o,comm_pair,local,kernel_local,radius,batch,action,status,metrics)
    integer,intent(in) :: n(3),m(3),lo(3),comm_r,comm_o,comm_pair,batch
    real(8),intent(in) :: h(3),radius
    complex(8),intent(in) :: local(:,:),kernel_local(:)
    complex(8),intent(out) :: action(:,:)
    integer,intent(out) :: status
    type(patch_metrics),intent(out) :: metrics
    type(owner_layout) :: layout
    type(neighbor_metrics) :: neighbor
    character(32) :: selector
    integer :: env_status,env_length,neighbor_on,selection_min,selection_max
    type(orbital_patch),allocatable :: patches(:)
    complex(8),allocatable :: owned(:,:),owned_action(:,:)
    real(8),allocatable :: contribution(:),kernel(:)
    integer,allocatable :: rows(:),cols(:),counts(:)
    integer :: coupled,orb_rank,orb_peers,rank,e,ng,nt,norb,bad,x,y,z,g,p(3),q,neg,i,first
    integer :: par(4),pmin(4),pmax(4),spatial_peers,coupled_peers
    integer(int64) :: volume,local_volume
    real(8) :: reals(4),rmin(4),rmax(4)
    real(8) :: delta,scale,imaginary,discrepancy,roundoff,clock0,clock1,clock2,clock3,times(3),max_times(3)
    real(8) :: entry_time,adapter_times(5),adapter_min(5),adapter_max(5),adapter_mean(5)
    entry_time=MPI_Wtime()
    action=0;metrics=patch_metrics();status=1
    call MPI_Comm_dup(comm_pair,coupled,e)
    call MPI_Comm_rank(coupled,rank,e)
    call MPI_Comm_rank(comm_o,orb_rank,e);call MPI_Comm_size(comm_o,orb_peers,e)
    ! Constant-count collective preflight precedes products, allocations and variable-count MPI.
    selector='';neighbor_on=1
    call get_environment_variable('SALMON_HSE_WANNIER_NEIGHBOR',selector,length=env_length,status=env_status)
    bad=0
    if(env_status/=1)then
      if(env_status/=0.or.env_length/=1)then
        bad=1
      elseif(selector=='0')then
        neighbor_on=0
      elseif(selector/='1')then
        bad=1
      endif
    endif
    call MPI_Allreduce(neighbor_on,selection_min,1,MPI_INTEGER,MPI_MIN,coupled,e)
    call MPI_Allreduce(neighbor_on,selection_max,1,MPI_INTEGER,MPI_MAX,coupled,e)
    if(selection_min/=selection_max)bad=1
    nt=size(local,2);volume=1_int64;local_volume=1_int64
    par=[n,batch]
    call MPI_Allreduce(par,pmin,4,MPI_INTEGER,MPI_MIN,coupled,e)
    call MPI_Allreduce(par,pmax,4,MPI_INTEGER,MPI_MAX,coupled,e)
    if(any(pmin/=pmax).or.any(n<1).or.any(m<1).or.any(lo<0).or.batch<1)bad=1
    do i=1,3
      if(int(lo(i),int64)+int(m(i),int64)>int(n(i),int64))bad=1
      if(n(i)>0)then
        if(volume>int(huge(ng),int64)/int(n(i),int64))then
          bad=1
        else
          volume=volume*int(n(i),int64)
        endif
      endif
      if(m(i)>0)then
        if(local_volume>int(huge(ng),int64)/int(m(i),int64))then
          bad=1
        else
          local_volume=local_volume*int(m(i),int64)
        endif
      endif
      if(.not.ieee_is_finite(h(i)).or.h(i)<=0d0)bad=1
    enddo
    if(local_volume/=size(local,1,kind=int64).or.local_volume/=size(kernel_local,kind=int64))bad=1
    if(local_volume*int(nt,int64)>int(huge(ng),int64))bad=1
    if(volume*int(max(1,batch),int64)>int(huge(ng),int64))bad=1
    if(any(shape(action)/=shape(local)))bad=1
    if(.not.ieee_is_finite(radius).or.radius<=0d0)bad=1
    call MPI_Comm_size(comm_r,spatial_peers,e);call MPI_Comm_size(coupled,coupled_peers,e)
    if(int(spatial_peers,int64)*int(orb_peers,int64)/=int(coupled_peers,int64))bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,coupled,e)
    if(status/=0)then
      call MPI_Comm_free(coupled,e);return
    endif
    reals=[h,radius]
    call MPI_Allreduce(reals,rmin,4,MPI_DOUBLE_PRECISION,MPI_MIN,coupled,e)
    call MPI_Allreduce(reals,rmax,4,MPI_DOUBLE_PRECISION,MPI_MAX,coupled,e)
    if(any(abs(rmin-rmax)>0d0))bad=1
    if(radius>sqrt(huge(1d0)))bad=1
    do i=1,3
      if(h(i)>sqrt(huge(1d0))/real(n(i),8))bad=1
    enddo
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,coupled,e)
    if(status/=0)then
      call MPI_Comm_free(coupled,e);return
    endif
    ng=int(volume)
    allocate(counts(orb_peers),rows(product(m)),cols(nt),contribution(ng),kernel(ng))
    call MPI_Allgather(nt,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm_o,e)
    if(sum(int(counts,int64))>int(huge(norb),int64))bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,coupled,e)
    if(status/=0)then
      call MPI_Comm_free(coupled,e);return
    endif
    norb=sum(counts);first=sum(counts(:orb_rank))
    do i=1,nt
      cols(i)=first+i
    enddo
    g=0;contribution=0;imaginary=0;scale=0
    do z=0,m(3)-1
      do y=0,m(2)-1
        do x=0,m(1)-1
          g=g+1;p=[x,y,z]+lo;rows(g)=1+p(1)+n(1)*(p(2)+n(2)*p(3))
          where(p>=(n+1)/2)p=p-n
          if(.not.ieee_is_finite(real(kernel_local(g),8)))bad=1
          if(.not.ieee_is_finite(aimag(kernel_local(g))))bad=1
          ! Preserve the existing spherical real-space kernel mask and normalization.
          if(sum((p*h)**2)>radius**2)cycle
          contribution(rows(g))=real(kernel_local(g),8)
          imaginary=max(imaginary,abs(aimag(kernel_local(g))))
          scale=max(scale,abs(kernel_local(g)))
        enddo
      enddo
    enddo
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,coupled,e)
    if(status/=0)then
      call MPI_Comm_free(coupled,e);return
    endif
    call MPI_Allreduce(contribution,kernel,ng,MPI_DOUBLE_PRECISION,MPI_SUM,comm_r,e)
    call MPI_Allreduce(MPI_IN_PLACE,scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,coupled,e)
    call MPI_Allreduce(MPI_IN_PLACE,imaginary,1,MPI_DOUBLE_PRECISION,MPI_MAX,coupled,e)
    discrepancy=0
    do z=0,n(3)-1
      do y=0,n(2)-1
        do x=0,n(1)-1
          q=1+x+n(1)*(y+n(2)*z)
          neg=1+modulo(-x,n(1))+n(1)*(modulo(-y,n(2))+n(2)*modulo(-z,n(3)))
          discrepancy=max(discrepancy,abs(kernel(q)-kernel(neg)))
          contribution(q)=0.5d0*kernel(q)+0.5d0*kernel(neg)
        enddo
      enddo
    enddo
    roundoff=64d0*epsilon(1d0)*max(tiny(1d0),scale)
    if(imaginary>roundoff.or.discrepancy>roundoff)bad=1
    call MPI_Allreduce(bad,status,1,MPI_INTEGER,MPI_MAX,coupled,e)
    if(status/=0)then
      call MPI_Comm_free(coupled,e);return
    endif
    kernel=contribution
    if(rank==0)write(*,*)'OWNER_KERNEL_ROUNDOFF imaginary/symmetry/allowance ',imaginary,discrepancy,roundoff
    clock0=MPI_Wtime()
    call layout_prepare(layout,ng,norb,rows,cols,coupled,status)
    if(status/=0)then
      call MPI_Comm_free(coupled,e);return
    endif
    allocate(owned(ng,size(layout%ids)),owned_action(ng,size(layout%ids)),patches(size(layout%ids)))
    call layout_forward(layout,local,owned,status)
    if(status/=0)then
      call MPI_Comm_free(coupled,e);return
    endif
    do i=1,size(layout%ids)
      call pack_exact_patch(owned(:,i),patches(i))
    enddo
    deallocate(owned)
    clock1=MPI_Wtime()
    if(neighbor_on==1)then
      call patch_apply(n,norb,layout%ids,patches,kernel,coupled,batch,owned_action,status,metrics,detail=neighbor)
      if(rank==0.and.status==0)then
        write(*,*)'WANNIER_NEIGHBOR_OMP requested/effective/active/workspace_bytes ', &
          neighbor%omp_requested,neighbor%omp_effective,neighbor%omp_active,neighbor%thread_workspace_bytes
        write(*,*)'WANNIER_NEIGHBOR_OMP_TIME rankmax_setup/pair_compute ', &
          neighbor%fft_setup_seconds,neighbor%pair_compute_seconds
        write(*,*)'WANNIER_NEIGHBOR tiles/pairtiles/plans/source_messages/return_messages ', &
          neighbor%tiles,neighbor%pair_tiles,neighbor%plans,neighbor%source_messages,neighbor%return_messages
        write(*,*)'WANNIER_NEIGHBOR_DIRECT full_block_pairs/visited/fetch_requests/receive_hits/extract_hits/extracted_values ', &
          neighbor%block_candidates,neighbor%block_visits,neighbor%fetch_requests,neighbor%receive_hits, &
          neighbor%extract_hits,neighbor%extracted_values
        write(*,*)'WANNIER_NEIGHBOR_CACHE rankmax_bytes ',neighbor%halo_cache_bytes
        do i=1,size(neighbor_phase_names)
          write(*,*)'WANNIER_NEIGHBOR_PHASE ',trim(neighbor_phase_names(i)),' rankmin/mean/max ', &
            neighbor%phase_min(i),neighbor%phase_mean(i),neighbor%phase_max(i)
        enddo
      endif
    else
      call old_patch_apply(n,norb,layout%ids,patches,kernel,coupled,batch,owned_action,status,metrics)
    endif
    if(status/=0)then
      call MPI_Comm_free(coupled,e);return
    endif
    clock2=MPI_Wtime()
    call layout_reverse(layout,owned_action,action,status)
    clock3=MPI_Wtime();times=[clock1-clock0,clock2-clock1,clock3-clock2]
    call MPI_Allreduce(times,max_times,3,MPI_DOUBLE_PRECISION,MPI_MAX,coupled,e)
    adapter_times=[clock0-entry_time,clock1-clock0,clock2-clock1,clock3-clock2,clock3-entry_time]
    call MPI_Allreduce(adapter_times,adapter_min,5,MPI_DOUBLE_PRECISION,MPI_MIN,coupled,e)
    call MPI_Allreduce(adapter_times,adapter_max,5,MPI_DOUBLE_PRECISION,MPI_MAX,coupled,e)
    call MPI_Allreduce(adapter_times,adapter_mean,5,MPI_DOUBLE_PRECISION,MPI_SUM,coupled,e)
    adapter_mean=adapter_mean/real(coupled_peers,8)
    if(rank==0)then
      write(*,*)'WANNIER_OWNER_ACTION orbitals/pairtasks/evaluated/zero/fftpoints/maxfft/sourcebytes/reversebytes ', &
        norb,metrics%pairs,metrics%evaluated,metrics%zero_overlap,metrics%fft_points,metrics%max_fft_points, &
        metrics%source_payload,metrics%reverse_payload
      write(*,*)'WANNIER_OWNER_PHASE transpose/patch/inverse ',max_times
      write(*,*)'WANNIER_OWNER_PHASE_MIN setup/transpose/patch/inverse/total ',adapter_min
      write(*,*)'WANNIER_OWNER_PHASE_MEAN setup/transpose/patch/inverse/total ',adapter_mean
      write(*,*)'WANNIER_OWNER_PHASE_MAX setup/transpose/patch/inverse/total ',adapter_max
    endif
    call MPI_Comm_free(coupled,e)
  end subroutine
end module
