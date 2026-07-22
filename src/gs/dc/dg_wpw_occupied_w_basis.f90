module dg_wpw_occupied_w_basis
  use,intrinsic::iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  type,public::s_dg_wpw_occupied_w_basis
    logical::valid=.false.
    integer::local_fragment=0,local_count=0,global_count=0,epoch=-1
    integer(int64)::fingerprint=0_int64
    real(8)::source_condition=huge(1d0)
    integer::buffer_lo(3)=0,buffer_hi(3)=-1
    integer::gradient_lo(3)=0,gradient_hi(3)=-1
    integer,allocatable::owned_ids(:),stable_keys(:,:),core_grid_ids(:)
    complex(8),allocatable::core_values(:,:),buffer_values(:,:),buffer_gradients(:,:,:)
  end type
  public::initialize_dg_wpw_occupied_w_basis
  public::initialize_dg_wpw_occupied_w_basis_collective
  public::gather_dg_wpw_occupied_w_payload
  public::broadcast_dg_wpw_occupied_w_basis
  public::evaluate_dg_wpw_occupied_w_point
  public::dg_wpw_unwrapped_to_storage_index
  public::reorder_dg_wpw_fragment_buffer
  public::extract_dg_wpw_canonical_cell
contains
  integer function dg_wpw_unwrapped_to_storage_index(local_grid,domain,buffer)result(storage)
    integer,intent(in)::local_grid,domain,buffer
    storage=0
    if(domain<=0.or.buffer<0.or.local_grid<1-buffer.or.local_grid>domain+buffer)return
    if(local_grid>=1)then
      storage=local_grid
    else
      storage=domain+2*buffer+local_grid
    endif
  end function dg_wpw_unwrapped_to_storage_index

  subroutine reorder_dg_wpw_fragment_buffer(storage_values,domain,buffer,unwrapped_values,info)
    complex(8),intent(in)::storage_values(:,:,:,:)
    integer,intent(in)::domain(3),buffer(3)
    complex(8),intent(out)::unwrapped_values(:,:,:,:)
    integer,intent(out)::info
    integer::extent(3),ix,iy,iz,sx,sy,sz

    unwrapped_values=(0d0,0d0);info=1;extent=domain+2*buffer
    if(any(domain<=0).or.any(buffer<0).or.any(shape(storage_values)/=&
      [extent,size(storage_values,4)]).or.any(shape(unwrapped_values)/=&
      [extent,size(storage_values,4)]))return
    do iz=1,extent(3)
      sz=dg_wpw_unwrapped_to_storage_index(iz-buffer(3),domain(3),buffer(3))
      do iy=1,extent(2)
        sy=dg_wpw_unwrapped_to_storage_index(iy-buffer(2),domain(2),buffer(2))
        do ix=1,extent(1)
          sx=dg_wpw_unwrapped_to_storage_index(ix-buffer(1),domain(1),buffer(1))
          if(min(sx,min(sy,sz))<=0)return
          unwrapped_values(ix,iy,iz,:)=storage_values(sx,sy,sz,:)
        enddo
      enddo
    enddo
    if(.not.all(finite_complex(unwrapped_values)))then
      unwrapped_values=(0d0,0d0);return
    endif
    info=0
  end subroutine reorder_dg_wpw_fragment_buffer

  subroutine extract_dg_wpw_canonical_cell(unwrapped_values,domain,buffer,total_shape,&
      fragment_origin,canonical_values,info)
    ! fragment_origin is the zero-based global index of fragment core point 1.
    complex(8),intent(in)::unwrapped_values(:,:,:,:)
    integer,intent(in)::domain(3),buffer(3),total_shape(3),fragment_origin(3)
    complex(8),intent(out)::canonical_values(:,:,:,:)
    integer,intent(out)::info
    logical,allocatable::seen(:,:,:)
    integer::extent(3),ix,iy,iz,cx,cy,cz
    real(8)::scale,tolerance

    info=1;canonical_values=(0d0,0d0);extent=domain+2*buffer
    if(any(domain<=0).or.any(buffer<0).or.any(total_shape<=0).or.any(extent<total_shape).or.&
        any(fragment_origin<0).or.any(fragment_origin>=total_shape).or.&
        size(unwrapped_values,4)<=0.or.&
        any(shape(unwrapped_values)/=[extent,size(unwrapped_values,4)]).or.&
        any(shape(canonical_values)/=[total_shape,size(unwrapped_values,4)]).or.&
        .not.all(finite_complex(unwrapped_values)))return
    allocate(seen(total_shape(1),total_shape(2),total_shape(3)));seen=.false.
    do iz=1,extent(3);do iy=1,extent(2);do ix=1,extent(1)
      cx=modulo(fragment_origin(1)+ix-buffer(1)-1,total_shape(1))+1
      cy=modulo(fragment_origin(2)+iy-buffer(2)-1,total_shape(2))+1
      cz=modulo(fragment_origin(3)+iz-buffer(3)-1,total_shape(3))+1
      if(.not.seen(cx,cy,cz))then
        canonical_values(cx,cy,cz,:)=unwrapped_values(ix,iy,iz,:)
        seen(cx,cy,cz)=.true.
      else
        scale=max(1d0,maxval(abs(canonical_values(cx,cy,cz,:))),&
          maxval(abs(unwrapped_values(ix,iy,iz,:))))
        tolerance=100d0*epsilon(1d0)*scale
        if(maxval(abs(canonical_values(cx,cy,cz,:)-unwrapped_values(ix,iy,iz,:)))>tolerance)then
          canonical_values=(0d0,0d0);return
        endif
      endif
    enddo;enddo;enddo
    if(.not.all(seen))then;canonical_values=(0d0,0d0);return;endif
    info=0
  end subroutine extract_dg_wpw_canonical_cell

  subroutine evaluate_dg_wpw_occupied_w_point(basis,local_grid,values,gradients,info)
    type(s_dg_wpw_occupied_w_basis),intent(in)::basis
    integer,intent(in)::local_grid(3)
    complex(8),intent(out)::values(:),gradients(:,:)
    integer,intent(out)::info
    integer::extent(3),point

    values=(0d0,0d0);gradients=(0d0,0d0);info=1
    if(.not.basis%valid.or.size(values)/=basis%local_count.or.&
      any(shape(gradients)/=[3,basis%local_count]))return
    extent=basis%buffer_hi-basis%buffer_lo+1
    if(any(extent<=0).or.size(basis%buffer_values,1)/=product(extent).or.&
      any(shape(basis%buffer_gradients)/=[3,product(extent),basis%local_count]))return
    if(any(local_grid<basis%buffer_lo).or.any(local_grid>basis%buffer_hi))return
    point=local_grid(1)-basis%buffer_lo(1)+1+&
      (local_grid(2)-basis%buffer_lo(2))*extent(1)+&
      (local_grid(3)-basis%buffer_lo(3))*extent(1)*extent(2)
    values=basis%buffer_values(point,:)
    gradients=basis%buffer_gradients(:,point,:)
    if(.not.all(finite_complex(values)).or..not.all(finite_complex(gradients)))then
      values=(0d0,0d0);gradients=(0d0,0d0);return
    endif
    info=0
  end subroutine evaluate_dg_wpw_occupied_w_point

  subroutine broadcast_dg_wpw_occupied_w_basis(comm,root,basis,info)
    use mpi,only:MPI_Comm_rank,MPI_Bcast,MPI_Allreduce,MPI_INTEGER,MPI_INTEGER8,&
      MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_MAX,MPI_SUCCESS
    integer,intent(in)::comm,root
    type(s_dg_wpw_occupied_w_basis),intent(inout)::basis
    integer,intent(out)::info
    type(s_dg_wpw_occupied_w_basis)::candidate
    integer::rank,ierr,metadata(18),local_bad,global_bad,ncore,nbuffer
    integer(int64)::fingerprint

    info=1;metadata=0;fingerprint=0_int64;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    if(rank==root)then
      if(.not.basis%valid)then
        local_bad=1
      else
        candidate=basis
        metadata=[basis%local_fragment,basis%local_count,basis%global_count,basis%epoch,&
          basis%buffer_lo,basis%buffer_hi,size(basis%core_grid_ids),size(basis%buffer_values,1),&
          basis%gradient_lo,basis%gradient_hi]
        fingerprint=basis%fingerprint
      endif
    endif
    call MPI_Bcast(metadata,size(metadata),MPI_INTEGER,root,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    call MPI_Bcast(fingerprint,1,MPI_INTEGER8,root,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    if(any(metadata(1:4)<=0).or.any(metadata(8:10)<metadata(5:7)).or.&
      metadata(11)<=0.or.metadata(12)<=0.or.any(metadata(16:18)<metadata(13:15)).or.&
      any(metadata(13:15)<metadata(5:7)).or.any(metadata(16:18)>metadata(8:10)).or.&
      fingerprint==0_int64)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    ncore=metadata(11);nbuffer=metadata(12)
    if(rank/=root)then
      candidate%local_fragment=metadata(1);candidate%local_count=metadata(2)
      candidate%global_count=metadata(3);candidate%epoch=metadata(4)
      candidate%buffer_lo=metadata(5:7);candidate%buffer_hi=metadata(8:10)
      candidate%gradient_lo=metadata(13:15);candidate%gradient_hi=metadata(16:18)
      allocate(candidate%owned_ids(candidate%local_count),candidate%stable_keys(5,candidate%local_count),&
        candidate%core_grid_ids(ncore),candidate%core_values(ncore,candidate%local_count),&
        candidate%buffer_values(nbuffer,candidate%local_count),&
        candidate%buffer_gradients(3,nbuffer,candidate%local_count))
    endif
    call MPI_Bcast(candidate%source_condition,1,MPI_DOUBLE_PRECISION,root,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    call MPI_Bcast(candidate%owned_ids,size(candidate%owned_ids),MPI_INTEGER,root,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    call MPI_Bcast(candidate%stable_keys,size(candidate%stable_keys),MPI_INTEGER,root,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    call MPI_Bcast(candidate%core_grid_ids,ncore,MPI_INTEGER,root,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    call MPI_Bcast(candidate%core_values,size(candidate%core_values),MPI_DOUBLE_COMPLEX,root,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    call MPI_Bcast(candidate%buffer_values,size(candidate%buffer_values),MPI_DOUBLE_COMPLEX,root,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    call MPI_Bcast(candidate%buffer_gradients,size(candidate%buffer_gradients),MPI_DOUBLE_COMPLEX,root,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    candidate%fingerprint=fingerprint;candidate%valid=.true.
    local_bad=merge(0,1,basis_fingerprint(candidate)==fingerprint)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    basis=candidate;info=0
  end subroutine broadcast_dg_wpw_occupied_w_basis

  subroutine gather_dg_wpw_occupied_w_payload(comm,active,local_core_grid_ids,local_core_values,&
      local_buffer_point_ids,local_buffer_values,local_buffer_gradients,expected_core_points,&
      expected_buffer_points,core_grid_ids,core_values,buffer_values,buffer_gradients,info)
    use mpi,only:MPI_Comm_size,MPI_Allreduce,MPI_Allgather,MPI_Allgatherv,MPI_INTEGER,&
      MPI_DOUBLE_COMPLEX,MPI_MAX,MPI_MIN,MPI_SUCCESS
    integer,intent(in)::comm,local_core_grid_ids(:),local_buffer_point_ids(:),&
      expected_core_points,expected_buffer_points
    logical,intent(in)::active
    complex(8),intent(in)::local_core_values(:,:),local_buffer_values(:,:),local_buffer_gradients(:,:,:)
    integer,allocatable,intent(out)::core_grid_ids(:)
    complex(8),allocatable,intent(out)::core_values(:,:),buffer_values(:,:),buffer_gradients(:,:,:)
    integer,intent(out)::info
    integer::nrank,ierr,local_bad,global_bad,nsource,nsource_min,nsource_max,nsource_max_global,ncore,nbuffer,&
      total_core,total_buffer,i,j,k,rank,offset,source,axis,position
    integer,allocatable::core_counts(:),buffer_counts(:),core_displacements(:),buffer_displacements(:),&
      core_value_counts(:),buffer_value_counts(:),gradient_counts(:),core_value_displacements(:),&
      buffer_value_displacements(:),gradient_displacements(:),all_core_ids(:),all_buffer_ids(:),order(:)
    complex(8),allocatable::send_core(:),send_buffer(:),send_gradient(:),all_core(:),all_buffer(:),all_gradient(:)

    info=1;local_bad=0;nsource=size(local_core_values,2)
    if(active)then
      if(nsource<=0.or.size(local_buffer_values,2)/=nsource.or.size(local_buffer_gradients,3)/=nsource.or.&
        size(local_core_values,1)/=size(local_core_grid_ids).or.&
        size(local_buffer_values,1)/=size(local_buffer_point_ids).or.&
        any(shape(local_buffer_gradients)/=[3,size(local_buffer_point_ids),nsource]))local_bad=1
      ncore=size(local_core_grid_ids);nbuffer=size(local_buffer_point_ids)
      nsource_min=nsource;nsource_max=nsource
    else
      ncore=0;nbuffer=0;nsource_min=huge(1);nsource_max=0
    endif
    if(expected_core_points<=0.or.expected_buffer_points<=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    if(global_bad/=0)return
    call MPI_Allreduce(nsource_min,nsource,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    call MPI_Allreduce(nsource_max,nsource_max_global,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    if(nsource<=0.or.nsource/=nsource_max_global)return
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(core_counts(nrank),buffer_counts(nrank),core_displacements(nrank),buffer_displacements(nrank))
    call MPI_Allgather(ncore,1,MPI_INTEGER,core_counts,1,MPI_INTEGER,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    call MPI_Allgather(nbuffer,1,MPI_INTEGER,buffer_counts,1,MPI_INTEGER,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    core_displacements(1)=0;buffer_displacements(1)=0
    do i=2,nrank
      core_displacements(i)=core_displacements(i-1)+core_counts(i-1)
      buffer_displacements(i)=buffer_displacements(i-1)+buffer_counts(i-1)
    enddo
    total_core=sum(core_counts);total_buffer=sum(buffer_counts)
    local_bad=merge(0,1,total_core==expected_core_points.and.total_buffer==expected_buffer_points)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(all_core_ids(total_core),all_buffer_ids(total_buffer))
    call MPI_Allgatherv(local_core_grid_ids,ncore,MPI_INTEGER,all_core_ids,core_counts,core_displacements,&
      MPI_INTEGER,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    call MPI_Allgatherv(local_buffer_point_ids,nbuffer,MPI_INTEGER,all_buffer_ids,&
      buffer_counts,buffer_displacements,MPI_INTEGER,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    allocate(core_value_counts(nrank),buffer_value_counts(nrank),gradient_counts(nrank),&
      core_value_displacements(nrank),buffer_value_displacements(nrank),gradient_displacements(nrank))
    core_value_counts=nsource*core_counts;buffer_value_counts=nsource*buffer_counts
    gradient_counts=3*nsource*buffer_counts
    core_value_displacements=nsource*core_displacements
    buffer_value_displacements=nsource*buffer_displacements
    gradient_displacements=3*nsource*buffer_displacements
    allocate(send_core(nsource*ncore),send_buffer(nsource*nbuffer),send_gradient(3*nsource*nbuffer),&
      all_core(nsource*total_core),all_buffer(nsource*total_buffer),all_gradient(3*nsource*total_buffer))
    do i=1,ncore;do source=1,nsource
      send_core((i-1)*nsource+source)=local_core_values(i,source)
    enddo;enddo
    do i=1,nbuffer;do source=1,nsource
      send_buffer((i-1)*nsource+source)=local_buffer_values(i,source)
      do axis=1,3
        send_gradient(((i-1)*nsource+source-1)*3+axis)=local_buffer_gradients(axis,i,source)
      enddo
    enddo;enddo
    call MPI_Allgatherv(send_core,size(send_core),MPI_DOUBLE_COMPLEX,all_core,core_value_counts,&
      core_value_displacements,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    call MPI_Allgatherv(send_buffer,size(send_buffer),MPI_DOUBLE_COMPLEX,all_buffer,&
      buffer_value_counts,buffer_value_displacements,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    call MPI_Allgatherv(send_gradient,size(send_gradient),MPI_DOUBLE_COMPLEX,&
      all_gradient,gradient_counts,gradient_displacements,MPI_DOUBLE_COMPLEX,comm,ierr)
    if(.not.collective_mpi_status_ok(comm,ierr))return
    allocate(order(total_core));order=[(i,i=1,total_core)]
    do i=1,total_core-1
      k=i
      do j=i+1,total_core;if(all_core_ids(order(j))<all_core_ids(order(k)))k=j;enddo
      if(k/=i)then;position=order(i);order(i)=order(k);order(k)=position;endif
    enddo
    if(any(all_buffer_ids<1).or.any(all_buffer_ids>expected_buffer_points))local_bad=1
    do i=2,total_core;if(all_core_ids(order(i))<=all_core_ids(order(i-1)))local_bad=1;enddo
    allocate(core_grid_ids(total_core),core_values(total_core,nsource),buffer_values(total_buffer,nsource),&
      buffer_gradients(3,total_buffer,nsource))
    core_grid_ids=all_core_ids(order);core_values=(0d0,0d0);buffer_values=(0d0,0d0);buffer_gradients=(0d0,0d0)
    do i=1,total_core;do source=1,nsource
      core_values(i,source)=all_core((order(i)-1)*nsource+source)
    enddo;enddo
    do i=1,total_buffer
      position=all_buffer_ids(i)
      if(position>=1.and.position<=total_buffer)then
        if(any(abs(buffer_values(position,:))>0d0).or.any(abs(buffer_gradients(:,position,:))>0d0))local_bad=1
        do source=1,nsource
          buffer_values(position,source)=all_buffer((i-1)*nsource+source)
          do axis=1,3
            buffer_gradients(axis,position,source)=all_gradient(((i-1)*nsource+source-1)*3+axis)
          enddo
        enddo
      endif
    enddo
    do i=1,total_buffer
      if(count(all_buffer_ids==i)/=1)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    info=0
  end subroutine gather_dg_wpw_occupied_w_payload

  subroutine initialize_dg_wpw_occupied_w_basis_collective(basis,comm,local_fragment,local_keys,&
      core_grid_ids,core_values,buffer_lo,buffer_hi,buffer_values,buffer_gradients,&
      source_condition,epoch,expected_global_count,info,gradient_lo,gradient_hi)
    use mpi,only:MPI_Comm_size,MPI_Allreduce,MPI_Allgather,MPI_Allgatherv,MPI_INTEGER,MPI_MAX,&
      MPI_SUM,MPI_SUCCESS
    type(s_dg_wpw_occupied_w_basis),intent(inout)::basis
    integer,intent(in)::comm,local_fragment,local_keys(:,:),core_grid_ids(:),buffer_lo(3),buffer_hi(3)
    complex(8),intent(in)::core_values(:,:),buffer_values(:,:),buffer_gradients(:,:,:)
    real(8),intent(in)::source_condition
    integer,intent(in)::epoch,expected_global_count
    integer,intent(out)::info
    integer,intent(in),optional::gradient_lo(3),gradient_hi(3)
    type(s_dg_wpw_occupied_w_basis)::candidate
    integer::nrank,ierr,local_bad,global_bad,nlocal,nglobal,i
    integer,allocatable::counts(:),displacements(:),key_counts(:),key_displacements(:)
    integer,allocatable::flat_keys(:),global_keys(:,:),local_owners(:),global_owners(:)

    info=1;local_bad=0;nlocal=size(local_keys,2)
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(local_fragment<=0.or.nlocal<=0.or.size(local_keys,1)/=5.or.expected_global_count<=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(counts(nrank),displacements(nrank),key_counts(nrank),key_displacements(nrank))
    call MPI_Allgather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    displacements(1)=0
    do i=2,nrank;displacements(i)=displacements(i-1)+counts(i-1);enddo
    nglobal=sum(counts)
    key_counts=5*counts;key_displacements=5*displacements
    allocate(flat_keys(5*nglobal),global_keys(5,nglobal),local_owners(nlocal),global_owners(nglobal))
    call MPI_Allgatherv(local_keys,5*nlocal,MPI_INTEGER,flat_keys,key_counts,key_displacements,&
      MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    global_keys=reshape(flat_keys,[5,nglobal]);local_owners=local_fragment
    call MPI_Allgatherv(local_owners,nlocal,MPI_INTEGER,global_owners,counts,displacements,&
      MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    local_bad=merge(0,1,nglobal==expected_global_count)
    candidate=basis
    if(local_bad==0)call initialize_dg_wpw_occupied_w_basis(candidate,local_fragment,global_keys,&
      global_owners,local_keys,core_grid_ids,core_values,buffer_lo,buffer_hi,buffer_values,&
      buffer_gradients,source_condition,epoch,local_bad,gradient_lo,gradient_hi)
    local_bad=merge(0,1,local_bad==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    basis=candidate;info=0
  end subroutine initialize_dg_wpw_occupied_w_basis_collective

  subroutine initialize_dg_wpw_occupied_w_basis(basis,local_fragment,global_keys,global_owners,&
      local_keys,core_grid_ids,core_values,buffer_lo,buffer_hi,buffer_values,buffer_gradients,&
      source_condition,epoch,info,gradient_lo,gradient_hi)
    type(s_dg_wpw_occupied_w_basis),intent(inout)::basis
    integer,intent(in)::local_fragment,global_keys(:,:),global_owners(:),local_keys(:,:),core_grid_ids(:)
    integer,intent(in)::buffer_lo(3),buffer_hi(3),epoch
    complex(8),intent(in)::core_values(:,:),buffer_values(:,:),buffer_gradients(:,:,:)
    real(8),intent(in)::source_condition
    integer,intent(out)::info
    integer,intent(in),optional::gradient_lo(3),gradient_hi(3)
    type(s_dg_wpw_occupied_w_basis)::candidate
    integer,allocatable::order(:),expected_ids(:),expected_keys(:,:),local_order(:)
    integer::i,j,k,nlocal,nglobal,nbuffer,tmp,grad_lo(3),grad_hi(3)

    info=1;nlocal=size(local_keys,2);nglobal=size(global_keys,2)
    grad_lo=buffer_lo;grad_hi=buffer_hi
    if(present(gradient_lo))grad_lo=gradient_lo
    if(present(gradient_hi))grad_hi=gradient_hi
    if(local_fragment<=0.or.epoch<=0.or.nglobal<=0.or.nlocal<=0.or.&
       size(global_keys,1)/=5.or.size(local_keys,1)/=5.or.size(global_owners)/=nglobal.or.&
       any(global_owners<=0).or.count(global_owners==local_fragment)/=nlocal.or.&
       size(core_grid_ids)<=0.or.any(shape(core_values)/=[size(core_grid_ids),nlocal]).or.&
       any(buffer_hi<buffer_lo).or.any(grad_hi<grad_lo).or.any(grad_lo<buffer_lo).or.&
       any(grad_hi>buffer_hi).or..not.ieee_is_finite(source_condition).or.source_condition<=0d0.or.&
       .not.all(finite_complex(core_values)).or..not.all(finite_complex(buffer_values)).or.&
       .not.all(finite_complex(buffer_gradients)))return
    nbuffer=product(buffer_hi-buffer_lo+1)
    if(any(shape(buffer_values)/=[nbuffer,nlocal]).or.&
       any(shape(buffer_gradients)/=[3,nbuffer,nlocal]).or..not.strictly_increasing(core_grid_ids))return
    allocate(order(nglobal));order=[(i,i=1,nglobal)]
    do i=1,nglobal-1
      k=i
      do j=i+1,nglobal
        if(key_less(global_keys(:,order(j)),global_keys(:,order(k))))k=j
      enddo
      if(k/=i)then;tmp=order(i);order(i)=order(k);order(k)=tmp;endif
    enddo
    do i=2,nglobal
      if(all(global_keys(:,order(i))==global_keys(:,order(i-1))))return
    enddo
    allocate(expected_ids(nlocal),expected_keys(5,nlocal),local_order(nlocal));j=0
    do i=1,nglobal
      if(global_owners(order(i))/=local_fragment)cycle
      j=j+1;expected_ids(j)=i;expected_keys(:,j)=global_keys(:,order(i))
    enddo
    if(j/=nlocal)return
    local_order=0
    do i=1,nlocal
      do j=1,nlocal
        if(all(expected_keys(:,i)==local_keys(:,j)))then
          if(local_order(i)/=0)return
          local_order(i)=j
        endif
      enddo
      if(local_order(i)==0)return
    enddo
    candidate%local_fragment=local_fragment;candidate%local_count=nlocal
    candidate%global_count=nglobal;candidate%epoch=epoch;candidate%source_condition=source_condition
    candidate%buffer_lo=buffer_lo;candidate%buffer_hi=buffer_hi
    candidate%gradient_lo=grad_lo;candidate%gradient_hi=grad_hi
    candidate%owned_ids=expected_ids;candidate%stable_keys=local_keys(:,local_order)
    candidate%core_grid_ids=core_grid_ids;candidate%core_values=core_values(:,local_order)
    candidate%buffer_values=buffer_values(:,local_order)
    candidate%buffer_gradients=buffer_gradients(:,:,local_order)
    candidate%fingerprint=basis_fingerprint(candidate)
    if(candidate%fingerprint==0_int64)candidate%fingerprint=1_int64
    candidate%valid=.true.;basis=candidate;info=0
  end subroutine initialize_dg_wpw_occupied_w_basis

  logical function collective_mpi_status_ok(comm,local_ierr)result(ok)
    use mpi,only:MPI_Allreduce,MPI_INTEGER,MPI_MAX,MPI_SUCCESS
    integer,intent(in)::comm,local_ierr
    integer::local_bad,global_bad,ierr
    local_bad=merge(0,1,local_ierr==MPI_SUCCESS)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.global_bad==0
  end function collective_mpi_status_ok

  logical function key_less(left,right)result(less)
    integer,intent(in)::left(:),right(:)
    integer::i
    less=.false.
    do i=1,size(left)
      if(left(i)<right(i))then;less=.true.;return
      elseif(left(i)>right(i))then;return
      endif
    enddo
  end function key_less

  logical function strictly_increasing(values)result(ok)
    integer,intent(in)::values(:)
    integer::i
    ok=size(values)>0
    do i=2,size(values)
      if(values(i)<=values(i-1))then;ok=.false.;return;endif
    enddo
  end function strictly_increasing

  elemental logical function finite_complex(value)result(ok)
    complex(8),intent(in)::value
    ok=ieee_is_finite(real(value,8)).and.ieee_is_finite(aimag(value))
  end function finite_complex

  integer(int64) function basis_fingerprint(basis)result(hash)
    type(s_dg_wpw_occupied_w_basis),intent(in)::basis
    integer::i,j,k
    hash=1469598103934665603_int64
    call mix(hash,int(basis%local_fragment,int64));call mix(hash,int(basis%epoch,int64))
    call mix(hash,int(basis%local_count,int64));call mix(hash,int(basis%global_count,int64))
    call mix(hash,transfer(basis%source_condition,hash))
    do i=1,3
      call mix(hash,int(basis%buffer_lo(i),int64));call mix(hash,int(basis%buffer_hi(i),int64))
      call mix(hash,int(basis%gradient_lo(i),int64));call mix(hash,int(basis%gradient_hi(i),int64))
    enddo
    do i=1,size(basis%owned_ids);call mix(hash,int(basis%owned_ids(i),int64));enddo
    do i=1,size(basis%core_grid_ids);call mix(hash,int(basis%core_grid_ids(i),int64));enddo
    do j=1,size(basis%stable_keys,2);do i=1,5
      call mix(hash,int(basis%stable_keys(i,j),int64))
    enddo;enddo
    do j=1,size(basis%core_values,2);do i=1,size(basis%core_values,1)
      call mix(hash,transfer(real(basis%core_values(i,j),8),hash))
      call mix(hash,transfer(aimag(basis%core_values(i,j)),hash))
    enddo;enddo
    do j=1,size(basis%buffer_values,2);do i=1,size(basis%buffer_values,1)
      call mix(hash,transfer(real(basis%buffer_values(i,j),8),hash))
      call mix(hash,transfer(aimag(basis%buffer_values(i,j)),hash))
    enddo;enddo
    do k=1,size(basis%buffer_gradients,3);do j=1,size(basis%buffer_gradients,2);do i=1,3
      call mix(hash,transfer(real(basis%buffer_gradients(i,j,k),8),hash))
      call mix(hash,transfer(aimag(basis%buffer_gradients(i,j,k)),hash))
    enddo;enddo;enddo
  end function basis_fingerprint

  subroutine mix(hash,value)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::value
    hash=ieor(hash,value);hash=ieor(hash,shiftl(hash,13));hash=ieor(hash,shiftr(hash,7))
  end subroutine mix
end module dg_wpw_occupied_w_basis
