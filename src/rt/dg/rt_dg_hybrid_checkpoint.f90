#include "config.h"
module rt_dg_hybrid_checkpoint
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::iso_c_binding,only:c_char,c_int,c_null_char
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  integer,parameter::checkpoint_version=2,legacy_checkpoint_version=1
  character(16),parameter::checkpoint_magic='SALMON_DG_HYB01 '
  character(16),parameter::occupied_magic='SALMON_DG_OCC02 '
  integer,parameter::occupied_version=2
  public::write_rt_dg_hybrid_checkpoint,read_rt_dg_hybrid_checkpoint,&
    write_rt_dg_hybrid_occupied_checkpoint,read_rt_dg_hybrid_occupied_checkpoint
  interface
    function c_rename(old_path,new_path) bind(C,name='rename') result(status)
      import::c_char,c_int
      character(c_char),intent(in)::old_path(*),new_path(*)
      integer(c_int)::status
    end function c_rename
  end interface
contains
  subroutine write_rt_dg_hybrid_checkpoint(comm,path,catalog_fingerprint,metric,operators,coefficients_owned,&
      state_fingerprint,payload_fingerprint,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::path
    integer(int64),intent(in)::catalog_fingerprint,state_fingerprint
    type(s_dg_hybrid_sparse_metric),intent(in)::metric
    type(s_dg_hybrid_sparse_operators),intent(in)::operators
    complex(real64),intent(in)::coefficients_owned(:)
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,ierr,n,nowned,row,i,k,root,position,unit,io_status,close_status,local_bad,global_bad,allocation_status
    integer::metric_degree,operator_degree,max_metric_degree,max_operator_degree
    integer(int64)::metadata_hash,minimum_metadata_hash,maximum_metadata_hash
    integer,allocatable::ownership(:),owner(:),owner_position(:),metric_degrees(:),operator_degrees(:),&
      metric_columns(:),operator_columns(:)
    complex(real64),allocatable::metric_values(:),operator_hamiltonian(:),operator_position(:,:)
    complex(real64)::coefficient
    logical::file_opened
    character(16)::path_probe
    character(:),allocatable::temporary_path
    ok=.false.;message='';payload_fingerprint=0_int64;n=metric%global_count;nowned=size(metric%owned_row_ids);local_bad=0
    io_status=-1;file_opened=.false.
    temporary_path=trim(path)//'.tmp.'//trim(int64_string(catalog_fingerprint))//'.'//trim(int64_string(state_fingerprint))
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call validate_path(path,path_probe,comm,ierr);if(ierr/=MPI_SUCCESS)then;message='inconsistent hybrid checkpoint path';return;endif
    if(.not.metric%valid.or..not.operators%valid.or.n<1.or.operators%global_count/=n)local_bad=1
    if(.not.ieee_is_finite(metric%condition_estimate).or.metric%condition_estimate<1d0.or.&
      .not.ieee_is_finite(metric%maximum_value).or.metric%maximum_value<=0d0)local_bad=1
    if(catalog_fingerprint==0_int64.or.state_fingerprint==0_int64.or.metric%fingerprint==0_int64.or.&
      operators%fingerprint==0_int64.or.operators%selection_fingerprint==0_int64.or.&
      operators%window_fingerprint==0_int64.or.operators%packet_fingerprint==0_int64.or.&
      operators%complement_fingerprint==0_int64.or.operators%metric_fingerprint==0_int64.or.&
      operators%position_convention_fingerprint==0_int64)local_bad=1
    if(size(coefficients_owned)/=nowned.or.size(operators%owned_row_ids)/=nowned)local_bad=1
    if(nowned==huge(0))local_bad=1
    if(size(metric%active_rows)/=n.or.size(metric%packet_ids)/=n)local_bad=1
    if(local_bad==0)then
      if(size(metric%row_offsets)/=nowned+1.or.size(operators%row_offsets)/=nowned+1)local_bad=1
    endif
    if(local_bad==0)then
      if(any(metric%owned_row_ids<1_int64).or.any(metric%owned_row_ids>int(n,int64)))local_bad=1
      if(any(metric%owned_row_ids/=operators%owned_row_ids).or..not.finite_vector(coefficients_owned))local_bad=1
      if(metric%row_offsets(1)/=1.or.operators%row_offsets(1)/=1)local_bad=1
      if(any(metric%row_offsets(2:)<metric%row_offsets(:nowned)).or.&
        any(operators%row_offsets(2:)<operators%row_offsets(:nowned)))local_bad=1
      if(metric%row_offsets(nowned+1)-1/=size(metric%column_ids).or.&
        operators%row_offsets(nowned+1)-1/=size(operators%column_ids))local_bad=1
      if(size(metric%values)/=size(metric%column_ids).or.&
        size(operators%hamiltonian_values)/=size(operators%column_ids).or.&
        size(operators%position_values,1)/=3.or.size(operators%position_values,2)/=size(operators%column_ids))local_bad=1
      if(any(metric%column_ids<1).or.any(metric%column_ids>n).or.any(operators%column_ids<1).or.&
        any(operators%column_ids>n))local_bad=1
      if(.not.finite_vector(metric%values).or..not.finite_vector(operators%hamiltonian_values).or.&
        .not.finite_matrix(operators%position_values))local_bad=1
      if(metric%numerical_rank/=count(metric%active_rows).or.&
        .not.valid_packet_activity(metric%packet_ids,metric%active_rows))local_bad=1
      do i=1,nowned
        if(.not.strictly_increasing(metric%column_ids(metric%row_offsets(i):metric%row_offsets(i+1)-1)).or.&
          .not.strictly_increasing(operators%column_ids(operators%row_offsets(i):operators%row_offsets(i+1)-1)))local_bad=1
      enddo
    endif
    if(local_bad==0)then
      if(any(metric%packet_ids<1))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid checkpoint write contract';return;endif
    metadata_hash=catalog_fingerprint
    metadata_hash=mix_hash(metadata_hash,state_fingerprint);metadata_hash=mix_hash(metadata_hash,metric%fingerprint)
    metadata_hash=mix_hash(metadata_hash,operators%fingerprint);metadata_hash=mix_hash(metadata_hash,int(n,int64))
    metadata_hash=mix_hash(metadata_hash,int(metric%numerical_rank,int64))
    metadata_hash=mix_hash(metadata_hash,transfer(metric%condition_estimate,metadata_hash))
    metadata_hash=mix_hash(metadata_hash,transfer(metric%maximum_value,metadata_hash))
    metadata_hash=mix_hash(metadata_hash,operators%selection_fingerprint)
    metadata_hash=mix_hash(metadata_hash,operators%window_fingerprint)
    metadata_hash=mix_hash(metadata_hash,operators%packet_fingerprint)
    metadata_hash=mix_hash(metadata_hash,operators%complement_fingerprint)
    metadata_hash=mix_hash(metadata_hash,operators%metric_fingerprint)
    metadata_hash=mix_hash(metadata_hash,operators%position_convention_fingerprint)
    do row=1,n
      metadata_hash=mix_hash(metadata_hash,merge(1_int64,0_int64,metric%active_rows(row)))
      metadata_hash=mix_hash(metadata_hash,int(metric%packet_ids(row),int64))
    enddo
    call MPI_Allreduce(metadata_hash,minimum_metadata_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(metadata_hash,maximum_metadata_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_metadata_hash/=maximum_metadata_hash)then
      message='rank-disagreeing hybrid checkpoint metadata';return
    endif
    allocate(ownership(n),owner(n),owner_position(n),metric_degrees(n),operator_degrees(n),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate hybrid checkpoint metadata';return;endif
    ownership=0;owner=-1;owner_position=0;metric_degrees=0;operator_degrees=0
    do i=1,nowned
      row=int(metric%owned_row_ids(i));ownership(row)=ownership(row)+1;owner(row)=rank;owner_position(row)=i
      metric_degrees(row)=metric%row_offsets(i+1)-metric%row_offsets(i)
      operator_degrees(row)=operators%row_offsets(i+1)-operators%row_offsets(i)
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,n,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,owner,n,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,owner_position,n,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,metric_degrees,n,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,operator_degrees,n,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    if(any(ownership/=1).or.any(metric_degrees<0).or.any(operator_degrees<0))local_bad=1
    max_metric_degree=maxval(metric_degrees);max_operator_degree=maxval(operator_degrees)
    if(max_operator_degree>huge(0)/3)local_bad=1
    allocate(metric_columns(max(1,max_metric_degree)),metric_values(max(1,max_metric_degree)),&
      operator_columns(max(1,max_operator_degree)),&
      operator_hamiltonian(max(1,max_operator_degree)),operator_position(3,max(1,max_operator_degree)),stat=allocation_status)
    if(allocation_status/=0)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='invalid or unallocatable checkpoint row workspace';return;endif
    io_status=0
    if(rank==0)then
      open(newunit=unit,file=temporary_path,status='replace',access='stream',form='unformatted',action='write',iostat=io_status)
      file_opened=io_status==0
    endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)then;call cleanup();message='cannot open hybrid checkpoint for writing';return;endif
    payload_fingerprint=catalog_fingerprint;call hash_int(state_fingerprint);call hash_int(metric%fingerprint)
    call hash_int(operators%fingerprint);call hash_int(int(n,int64));call hash_int(int(metric%numerical_rank,int64))
    call hash_int(transfer(metric%condition_estimate,payload_fingerprint));call hash_int(transfer(metric%maximum_value,payload_fingerprint))
    call hash_int(operators%selection_fingerprint);call hash_int(operators%window_fingerprint)
    call hash_int(operators%packet_fingerprint);call hash_int(operators%complement_fingerprint)
    call hash_int(operators%metric_fingerprint);call hash_int(operators%position_convention_fingerprint)
    do row=1,n
      call hash_int(merge(1_int64,0_int64,metric%active_rows(row)))
      call hash_int(int(metric%packet_ids(row),int64));call hash_int(int(metric_degrees(row),int64))
      call hash_int(int(operator_degrees(row),int64))
    enddo
    if(rank==0)then
      write(unit,iostat=io_status)checkpoint_magic,checkpoint_version,n,metric%numerical_rank,catalog_fingerprint,&
        state_fingerprint,metric%fingerprint,metric%condition_estimate,metric%maximum_value,&
        operators%selection_fingerprint,operators%window_fingerprint,operators%packet_fingerprint,&
        operators%complement_fingerprint,operators%metric_fingerprint,operators%position_convention_fingerprint,&
        operators%fingerprint,metric%active_rows,metric%packet_ids,metric_degrees,operator_degrees
    endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    do row=1,n
      root=owner(row);position=owner_position(row);metric_degree=metric_degrees(row);operator_degree=operator_degrees(row)
      if(rank==root)then
        metric_columns(1:metric_degree)=metric%column_ids(metric%row_offsets(position):metric%row_offsets(position+1)-1)
        metric_values(1:metric_degree)=metric%values(metric%row_offsets(position):metric%row_offsets(position+1)-1)
        operator_columns(1:operator_degree)=operators%column_ids(operators%row_offsets(position):operators%row_offsets(position+1)-1)
        operator_hamiltonian(1:operator_degree)=&
          operators%hamiltonian_values(operators%row_offsets(position):operators%row_offsets(position+1)-1)
        operator_position(:,1:operator_degree)=&
          operators%position_values(:,operators%row_offsets(position):operators%row_offsets(position+1)-1)
        coefficient=coefficients_owned(position)
      endif
      call MPI_Bcast(metric_columns,metric_degree,MPI_INTEGER,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(metric_values,metric_degree,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(operator_columns,operator_degree,MPI_INTEGER,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(operator_hamiltonian,operator_degree,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(operator_position,3*operator_degree,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call MPI_Bcast(coefficient,1,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call hash_row(row,metric_degree,operator_degree,metric_columns,metric_values,operator_columns,&
        operator_hamiltonian,operator_position,coefficient)
      if(rank==0)write(unit,iostat=io_status)metric_columns(1:metric_degree),metric_values(1:metric_degree),&
        operator_columns(1:operator_degree),operator_hamiltonian(1:operator_degree),&
        operator_position(:,1:operator_degree),coefficient
      call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    enddo
    if(payload_fingerprint==0_int64)payload_fingerprint=1_int64
    if(rank==0)then
      write(unit,iostat=io_status)payload_fingerprint
      close_status=0;close(unit,iostat=close_status);file_opened=.false.
      if(io_status==0)io_status=close_status
    endif
    call sync_io(io_status,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)then;call cleanup();message='hybrid checkpoint final write failed';return;endif
    if(rank==0)call atomic_rename(temporary_path,trim(path),io_status)
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)then;call cleanup();message='hybrid checkpoint atomic publication failed';return;endif
    call cleanup();ok=.true.;return
900 message='hybrid checkpoint MPI stream failed';if(rank==0.and.file_opened)close(unit);call cleanup();return
910 message='hybrid checkpoint file write failed';if(rank==0.and.file_opened)close(unit);call cleanup();return
#else
    ok=.false.;message='hybrid checkpoint requires MPI';payload_fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine hash_int(value)
      integer(int64),intent(in)::value;payload_fingerprint=ieor(ishftc(payload_fingerprint,9),value)
    end subroutine hash_int
    subroutine hash_complex(value)
      complex(real64),intent(in)::value;call hash_int(transfer(real(value),payload_fingerprint));&
        call hash_int(transfer(aimag(value),payload_fingerprint))
    end subroutine hash_complex
    subroutine hash_row(global_row,md,od,mc,mv,oc,oh,op,c)
      integer,intent(in)::global_row,md,od,mc(:),oc(:);complex(real64),intent(in)::mv(:),oh(:),op(:,:),c
      integer::a,b
      call hash_int(int(global_row,int64));call hash_int(int(md,int64));call hash_int(int(od,int64))
      do a=1,md;call hash_int(int(mc(a),int64));call hash_complex(mv(a));enddo
      do a=1,od
        call hash_int(int(oc(a),int64));call hash_complex(oh(a))
        do b=1,3;call hash_complex(op(b,a));enddo
      enddo
      call hash_complex(c)
    end subroutine hash_row
    subroutine cleanup()
      if(allocated(ownership))deallocate(ownership);if(allocated(owner))deallocate(owner)
      if(allocated(owner_position))deallocate(owner_position);if(allocated(metric_degrees))deallocate(metric_degrees)
      if(allocated(operator_degrees))deallocate(operator_degrees);if(allocated(metric_columns))deallocate(metric_columns)
      if(allocated(metric_values))deallocate(metric_values);if(allocated(operator_columns))deallocate(operator_columns)
      if(allocated(operator_hamiltonian))deallocate(operator_hamiltonian)
      if(allocated(operator_position))deallocate(operator_position)
    end subroutine cleanup
#endif
  end subroutine write_rt_dg_hybrid_checkpoint

  subroutine read_rt_dg_hybrid_checkpoint(comm,path,expected_catalog,expected_state,expected_selection,expected_window,&
      expected_packet,expected_complement,expected_metric,expected_position,expected_operator,metric,operators,&
      coefficients_owned,payload_fingerprint,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::path
    integer(int64),intent(in)::expected_catalog,expected_state,expected_selection,expected_window,expected_packet,&
      expected_complement,expected_metric,expected_position,expected_operator
    type(s_dg_hybrid_sparse_metric),intent(out)::metric
    type(s_dg_hybrid_sparse_operators),intent(out)::operators
    complex(real64),allocatable,intent(out)::coefficients_owned(:)
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,unit,io_status,version,n,numerical_rank,row,i,k,nowned,metric_nnz,operator_nnz,&
      metric_degree,operator_degree,max_metric_degree,max_operator_degree,allocation_status,position
    integer::integer_header(3),allocation_bad
    integer(int64)::fingerprint_header(3),operator_header(7)
    integer(int64)::catalog,state_fp,metric_fp,selection_fp,window_fp,packet_fp,complement_fp,operator_metric_fp,&
      position_fp,operator_fp,stored_fingerprint
    integer(int64)::file_size
    integer(int64)::expected_values(9),minimum_expected(9),maximum_expected(9)
    real(real64)::condition,maximum_value
    logical,allocatable::active_rows(:)
    integer,allocatable::packet_ids(:),metric_degrees(:),operator_degrees(:),metric_columns(:),operator_columns(:)
    complex(real64),allocatable::metric_values(:),operator_metric(:),operator_hamiltonian(:),operator_position(:,:)
    complex(real64)::coefficient
    logical::file_opened
    character(16)::magic,path_probe
    ok=.false.;message='';payload_fingerprint=0_int64;file_opened=.false.
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    expected_values=[expected_catalog,expected_state,expected_selection,expected_window,expected_packet,&
      expected_complement,expected_metric,expected_position,expected_operator]
    call MPI_Allreduce(expected_values,minimum_expected,9,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(expected_values,maximum_expected,9,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(any(minimum_expected/=maximum_expected).or.any(expected_values==0_int64))then
      message='rank-disagreeing expected checkpoint provenance';return
    endif
    call validate_path(path,path_probe,comm,ierr);if(ierr/=MPI_SUCCESS)then;message='inconsistent hybrid restart path';return;endif
    io_status=0;file_size=0_int64
    if(rank==0)then
      open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',action='read',iostat=io_status)
      file_opened=io_status==0
      if(io_status==0)inquire(unit=unit,size=file_size,iostat=io_status)
      if(io_status==0)read(unit,iostat=io_status)magic,version,n,numerical_rank,catalog,state_fp,metric_fp,condition,maximum_value,&
        selection_fp,window_fp,packet_fp,complement_fp,operator_metric_fp,position_fp,operator_fp
    endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0)then
      if(rank==0.and.file_opened)close(unit);message='cannot read hybrid checkpoint header';return
    endif
    call MPI_Bcast(file_size,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(magic,len(magic),MPI_CHARACTER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    if(rank==0)then
      integer_header=[version,n,numerical_rank]
      fingerprint_header=[catalog,state_fp,metric_fp]
      operator_header=[selection_fp,window_fp,packet_fp,complement_fp,operator_metric_fp,position_fp,operator_fp]
    endif
    call MPI_Bcast(integer_header,3,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(fingerprint_header,3,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(condition,1,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(maximum_value,1,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Bcast(operator_header,7,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)return
    version=integer_header(1);n=integer_header(2);numerical_rank=integer_header(3)
    catalog=fingerprint_header(1);state_fp=fingerprint_header(2);metric_fp=fingerprint_header(3)
    selection_fp=operator_header(1);window_fp=operator_header(2);packet_fp=operator_header(3)
    complement_fp=operator_header(4);operator_metric_fp=operator_header(5);position_fp=operator_header(6);operator_fp=operator_header(7)
    if(magic/=checkpoint_magic.or.(version/=checkpoint_version.and.version/=legacy_checkpoint_version).or.&
      n<1.or.numerical_rank<1.or.numerical_rank>n.or.&
      n>huge(0)/3.or.int(n,int64)>file_size/4_int64)then
      if(rank==0.and.file_opened)then;close(unit);file_opened=.false.;endif;message='incompatible hybrid checkpoint version';return
    endif
    if(catalog/=expected_catalog.or.state_fp/=expected_state.or.selection_fp/=expected_selection.or.&
      window_fp/=expected_window.or.packet_fp/=expected_packet.or.complement_fp/=expected_complement.or.&
      metric_fp/=expected_metric.or.position_fp/=expected_position.or.operator_fp/=expected_operator.or.&
      operator_metric_fp/=metric_fp)then
      if(rank==0.and.file_opened)then;close(unit);file_opened=.false.;endif;message='stale hybrid checkpoint provenance';return
    endif
    allocate(active_rows(n),packet_ids(n),metric_degrees(n),operator_degrees(n),stat=allocation_status)
    if(rank==0.and.allocation_status==0)read(unit,iostat=io_status)active_rows,packet_ids,metric_degrees,operator_degrees
    allocation_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(allocation_bad,k,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.k/=0)then
      if(rank==0.and.file_opened)then;close(unit);file_opened=.false.;endif
      call cleanup_buffers();message='cannot allocate checkpoint metadata';return
    endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 920
    call MPI_Bcast(active_rows,n,MPI_LOGICAL,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    call MPI_Bcast(packet_ids,n,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    call MPI_Bcast(metric_degrees,n,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    call MPI_Bcast(operator_degrees,n,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    if(any(packet_ids<1).or.any(packet_ids>n).or.any(metric_degrees<0).or.any(operator_degrees<0).or.n>huge(0)/3.or.&
      numerical_rank/=count(active_rows).or..not.valid_packet_activity(packet_ids,active_rows).or.&
      .not.ieee_is_finite(condition).or.condition<1d0.or..not.ieee_is_finite(maximum_value).or.maximum_value<=0d0)goto 920
    if(rank>=n)then;nowned=0;else;nowned=(n-1-rank)/nproc+1;endif
    metric_nnz=0;operator_nnz=0
    do row=1,n
      if(metric_degrees(row)>n.or.operator_degrees(row)>n)goto 920
      if(mod(row-1,nproc)==rank)then
        if(metric_degrees(row)>huge(0)-metric_nnz.or.operator_degrees(row)>huge(0)-operator_nnz)goto 920
        metric_nnz=metric_nnz+metric_degrees(row);operator_nnz=operator_nnz+operator_degrees(row)
      endif
    enddo
    max_metric_degree=maxval(metric_degrees);max_operator_degree=maxval(operator_degrees)
    allocate(metric_columns(max(1,max_metric_degree)),metric_values(max(1,max_metric_degree)),&
      operator_columns(max(1,max_operator_degree)),operator_metric(max(1,max_operator_degree)),&
      operator_hamiltonian(max(1,max_operator_degree)),operator_position(3,max(1,max_operator_degree)),&
      metric%owned_row_ids(nowned),metric%row_offsets(nowned+1),metric%column_ids(metric_nnz),metric%values(metric_nnz),&
      metric%active_rows(n),metric%packet_ids(n),operators%owned_row_ids(nowned),operators%row_offsets(nowned+1),&
      operators%column_ids(operator_nnz),operators%metric_values(operator_nnz),operators%hamiltonian_values(operator_nnz),&
      operators%position_values(3,operator_nnz),coefficients_owned(nowned),stat=allocation_status)
    allocation_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(allocation_bad,k,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.k/=0)goto 920
    metric%row_offsets(1)=1;operators%row_offsets(1)=1;position=0;metric_nnz=0;operator_nnz=0
    payload_fingerprint=catalog;call hash_int_read(state_fp);call hash_int_read(metric_fp);call hash_int_read(operator_fp)
    call hash_int_read(int(n,int64));call hash_int_read(int(numerical_rank,int64))
    call hash_int_read(transfer(condition,payload_fingerprint));call hash_int_read(transfer(maximum_value,payload_fingerprint))
    call hash_int_read(selection_fp);call hash_int_read(window_fp);call hash_int_read(packet_fp);call hash_int_read(complement_fp)
    call hash_int_read(operator_metric_fp);call hash_int_read(position_fp)
    do row=1,n
      call hash_int_read(merge(1_int64,0_int64,active_rows(row)))
      call hash_int_read(int(packet_ids(row),int64));call hash_int_read(int(metric_degrees(row),int64))
      call hash_int_read(int(operator_degrees(row),int64))
    enddo
    do row=1,n
      metric_degree=metric_degrees(row);operator_degree=operator_degrees(row)
      if(rank==0)then
        if(version==legacy_checkpoint_version)then
          read(unit,iostat=io_status)metric_columns(1:metric_degree),metric_values(1:metric_degree),&
            operator_columns(1:operator_degree),operator_metric(1:operator_degree),operator_hamiltonian(1:operator_degree),&
            operator_position(:,1:operator_degree),coefficient
        else
          read(unit,iostat=io_status)metric_columns(1:metric_degree),metric_values(1:metric_degree),&
            operator_columns(1:operator_degree),operator_hamiltonian(1:operator_degree),&
            operator_position(:,1:operator_degree),coefficient
        endif
      endif
      call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 920
      call MPI_Bcast(metric_columns,metric_degree,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      call MPI_Bcast(metric_values,metric_degree,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      call MPI_Bcast(operator_columns,operator_degree,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      if(version==legacy_checkpoint_version)then
        call MPI_Bcast(operator_metric,operator_degree,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      endif
      call MPI_Bcast(operator_hamiltonian,operator_degree,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      call MPI_Bcast(operator_position,3*operator_degree,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      call MPI_Bcast(coefficient,1,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
      if(any(metric_columns(1:metric_degree)<1).or.any(metric_columns(1:metric_degree)>n).or.&
        any(operator_columns(1:operator_degree)<1).or.any(operator_columns(1:operator_degree)>n))goto 920
      if(.not.strictly_increasing(metric_columns(1:metric_degree)).or.&
        .not.strictly_increasing(operator_columns(1:operator_degree)))goto 920
      if(.not.finite_vector(metric_values(1:metric_degree)).or.&
        .not.finite_vector(operator_hamiltonian(1:operator_degree)).or.&
        .not.finite_matrix(operator_position(:,1:operator_degree)).or..not.finite_vector([coefficient]))goto 920
      if(version==legacy_checkpoint_version)then
        if(.not.finite_vector(operator_metric(1:operator_degree)))goto 920
      else
        operator_metric(1:operator_degree)=(0d0,0d0)
        do i=1,operator_degree
          k=findloc(metric_columns(1:metric_degree),operator_columns(i),dim=1)
          if(k>0)operator_metric(i)=metric_values(k)
        enddo
      endif
      call hash_row_read(row,metric_degree,operator_degree,metric_columns,metric_values,operator_columns,operator_metric,&
        operator_hamiltonian,operator_position,coefficient)
      if(mod(row-1,nproc)==rank)then
        position=position+1;metric%owned_row_ids(position)=row;operators%owned_row_ids(position)=row
        metric%column_ids(metric_nnz+1:metric_nnz+metric_degree)=metric_columns(1:metric_degree)
        metric%values(metric_nnz+1:metric_nnz+metric_degree)=metric_values(1:metric_degree);metric_nnz=metric_nnz+metric_degree
        metric%row_offsets(position+1)=metric_nnz+1
        operators%column_ids(operator_nnz+1:operator_nnz+operator_degree)=operator_columns(1:operator_degree)
        operators%metric_values(operator_nnz+1:operator_nnz+operator_degree)=operator_metric(1:operator_degree)
        operators%hamiltonian_values(operator_nnz+1:operator_nnz+operator_degree)=operator_hamiltonian(1:operator_degree)
        operators%position_values(:,operator_nnz+1:operator_nnz+operator_degree)=operator_position(:,1:operator_degree)
        operator_nnz=operator_nnz+operator_degree;operators%row_offsets(position+1)=operator_nnz+1;coefficients_owned(position)=coefficient
      endif
    enddo
    if(rank==0)then;read(unit,iostat=io_status)stored_fingerprint;close(unit);file_opened=.false.;endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 920
    call MPI_Bcast(stored_fingerprint,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 920
    if(payload_fingerprint==0_int64)payload_fingerprint=1_int64
    if(stored_fingerprint/=payload_fingerprint)then;call cleanup_read();message='corrupt hybrid checkpoint payload';return;endif
    metric%valid=.true.;metric%global_count=n;metric%numerical_rank=numerical_rank;metric%condition_estimate=condition
    metric%maximum_value=maximum_value;metric%fingerprint=metric_fp;metric%active_rows=active_rows;metric%packet_ids=packet_ids
    metric%max_row_nnz=max_metric_degree
    operators%valid=.true.;operators%global_count=n;operators%selection_fingerprint=selection_fp
    operators%window_fingerprint=window_fp;operators%packet_fingerprint=packet_fp
    operators%complement_fingerprint=complement_fp
    operators%metric_fingerprint=operator_metric_fp;operators%position_convention_fingerprint=position_fp
    operators%fingerprint=operator_fp;ok=.true.;call cleanup_buffers();return
920 if(rank==0.and.file_opened)close(unit);call cleanup_read();message='hybrid checkpoint read failed';return
#else
    ok=.false.;message='hybrid checkpoint requires MPI';payload_fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine hash_int_read(value)
      integer(int64),intent(in)::value;payload_fingerprint=ieor(ishftc(payload_fingerprint,9),value)
    end subroutine hash_int_read
    subroutine hash_complex_read(value)
      complex(real64),intent(in)::value;call hash_int_read(transfer(real(value),payload_fingerprint));&
        call hash_int_read(transfer(aimag(value),payload_fingerprint))
    end subroutine hash_complex_read
    subroutine hash_row_read(global_row,md,od,mc,mv,oc,om,oh,op,c)
      integer,intent(in)::global_row,md,od,mc(:),oc(:);complex(real64),intent(in)::mv(:),om(:),oh(:),op(:,:),c
      integer::a,b
      call hash_int_read(int(global_row,int64));call hash_int_read(int(md,int64));call hash_int_read(int(od,int64))
      do a=1,md;call hash_int_read(int(mc(a),int64));call hash_complex_read(mv(a));enddo
      do a=1,od
        call hash_int_read(int(oc(a),int64))
        if(version==legacy_checkpoint_version)call hash_complex_read(om(a))
        call hash_complex_read(oh(a))
        do b=1,3;call hash_complex_read(op(b,a));enddo
      enddo
      call hash_complex_read(c)
    end subroutine hash_row_read
    subroutine cleanup_buffers()
      if(allocated(active_rows))deallocate(active_rows);if(allocated(packet_ids))deallocate(packet_ids)
      if(allocated(metric_degrees))deallocate(metric_degrees);if(allocated(operator_degrees))deallocate(operator_degrees)
      if(allocated(metric_columns))deallocate(metric_columns);if(allocated(metric_values))deallocate(metric_values)
      if(allocated(operator_columns))deallocate(operator_columns);if(allocated(operator_metric))deallocate(operator_metric)
      if(allocated(operator_hamiltonian))deallocate(operator_hamiltonian);if(allocated(operator_position))deallocate(operator_position)
    end subroutine cleanup_buffers
    subroutine cleanup_read()
      call cleanup_buffers()
      if(allocated(coefficients_owned))deallocate(coefficients_owned)
      if(allocated(metric%owned_row_ids))deallocate(metric%owned_row_ids)
      if(allocated(metric%row_offsets))deallocate(metric%row_offsets)
      if(allocated(metric%column_ids))deallocate(metric%column_ids)
      if(allocated(metric%values))deallocate(metric%values)
      if(allocated(metric%active_rows))deallocate(metric%active_rows)
      if(allocated(metric%packet_ids))deallocate(metric%packet_ids)
      if(allocated(operators%owned_row_ids))deallocate(operators%owned_row_ids)
      if(allocated(operators%row_offsets))deallocate(operators%row_offsets)
      if(allocated(operators%column_ids))deallocate(operators%column_ids)
      if(allocated(operators%metric_values))deallocate(operators%metric_values)
      if(allocated(operators%hamiltonian_values))deallocate(operators%hamiltonian_values)
      if(allocated(operators%position_values))deallocate(operators%position_values)
      metric%valid=.false.;operators%valid=.false.
    end subroutine cleanup_read
#endif
  end subroutine read_rt_dg_hybrid_checkpoint

  subroutine write_rt_dg_hybrid_occupied_checkpoint(comm,path,global_count,row_ids,coefficients,occupations,eigenvalues,&
      catalog_fingerprint,basis_fingerprint,provenance_fingerprints,operator_fingerprint,state_fingerprint,scf_receipts,&
      maximum_scf_residual,payload_fingerprint,ok,message)
    integer,intent(in)::comm,global_count
    character(*),intent(in)::path
    integer(int64),intent(in)::row_ids(:),catalog_fingerprint,basis_fingerprint,provenance_fingerprints(6),&
      operator_fingerprint,state_fingerprint
    complex(real64),intent(in)::coefficients(:,:)
    real(real64),intent(in)::occupations(:),eigenvalues(:),scf_receipts(5),maximum_scf_residual
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,noccupied,nlocal,i,j,row,root,position,unit,io_status,close_status,bad,global_bad,status
    integer,allocatable::counts(:),owners(:),positions(:)
    complex(real64),allocatable::row_values(:)
    integer(int64)::minimum_i,maximum_i,bits
    real(real64)::minimum_r,maximum_r
    character(16)::probe
    character(:),allocatable::temporary_path
    logical::opened
    ok=.false.;message='';payload_fingerprint=0_int64;opened=.false.;io_status=-1
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    call validate_path(path,probe,comm,ierr);if(ierr/=MPI_SUCCESS)then;message='inconsistent occupied checkpoint path';return;endif
    nlocal=size(row_ids);noccupied=size(coefficients,2);bad=0
    if(global_count<1.or.noccupied<1.or.size(coefficients,1)/=nlocal.or.size(occupations)/=noccupied.or.&
      size(eigenvalues)/=noccupied)bad=1
    call agree_int(global_count);call agree_int(noccupied)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid occupied checkpoint dimensions';return;endif
    if(any(row_ids<1_int64).or.any(row_ids>int(global_count,int64)).or..not.finite_matrix(coefficients).or.&
      any(.not.ieee_is_finite(occupations)).or.any(.not.ieee_is_finite(eigenvalues)).or.&
      any(.not.ieee_is_finite(scf_receipts)).or.any(occupations<0d0))bad=1
    if(catalog_fingerprint==0_int64.or.basis_fingerprint==0_int64.or.any(provenance_fingerprints==0_int64).or.&
      operator_fingerprint==0_int64.or.state_fingerprint==0_int64)bad=1
    if(.not.ieee_is_finite(maximum_scf_residual).or.maximum_scf_residual<=0d0.or.any(scf_receipts<0d0).or.&
      any(scf_receipts>maximum_scf_residual))bad=1
    call agree_i64(catalog_fingerprint);call agree_i64(basis_fingerprint);call agree_i64(operator_fingerprint);call agree_i64(state_fingerprint)
    do i=1,6;call agree_i64(provenance_fingerprints(i));enddo
    call agree_real(maximum_scf_residual)
    do i=1,noccupied;call agree_real(occupations(i));call agree_real(eigenvalues(i));enddo
    do i=1,5;call agree_real(scf_receipts(i));enddo
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid occupied checkpoint write contract';return;endif
    allocate(counts(global_count),owners(global_count),positions(global_count),row_values(noccupied),stat=status)
    bad=merge(0,1,status==0);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call clean;message='cannot allocate occupied checkpoint workspace';return;endif
    counts=0;owners=-1;positions=0
    do i=1,nlocal;row=int(row_ids(i));counts(row)=counts(row)+1;owners(row)=rank;positions(row)=i;enddo
    call MPI_Allreduce(MPI_IN_PLACE,counts,global_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,owners,global_count,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Allreduce(MPI_IN_PLACE,positions,global_count,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    bad=merge(0,1,all(counts==1));call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call clean;message='occupied checkpoint rows are not exactly once';return;endif
    temporary_path=trim(path)//'.occupied.tmp.'//trim(int64_string(state_fingerprint))
    if(rank==0)then
      open(newunit=unit,file=temporary_path,status='replace',access='stream',form='unformatted',action='write',iostat=io_status);opened=io_status==0
    endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    payload_fingerprint=catalog_fingerprint
    call occupied_hash(basis_fingerprint);do i=1,6;call occupied_hash(provenance_fingerprints(i));enddo
    call occupied_hash(operator_fingerprint);call occupied_hash(state_fingerprint)
    call occupied_hash(int(global_count,int64));call occupied_hash(int(noccupied,int64))
    do i=1,noccupied;call occupied_hash(transfer(occupations(i),bits));call occupied_hash(transfer(eigenvalues(i),bits));enddo
    do i=1,5;call occupied_hash(transfer(scf_receipts(i),bits));enddo
    if(rank==0)write(unit,iostat=io_status)occupied_magic,occupied_version,global_count,noccupied,catalog_fingerprint,&
      basis_fingerprint,provenance_fingerprints,operator_fingerprint,state_fingerprint,occupations,eigenvalues,scf_receipts
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    do row=1,global_count
      root=owners(row);position=positions(row);if(rank==root)row_values=coefficients(position,:)
      call MPI_Bcast(row_values,noccupied,MPI_DOUBLE_COMPLEX,root,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      call occupied_hash(int(row,int64));do j=1,noccupied;call occupied_hash(transfer(real(row_values(j)),bits));call occupied_hash(transfer(aimag(row_values(j)),bits));enddo
      if(rank==0)write(unit,iostat=io_status)row_values
      call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    enddo
    if(payload_fingerprint==0_int64)payload_fingerprint=1_int64
    if(rank==0)then
      write(unit,iostat=io_status)payload_fingerprint;close_status=0;close(unit,iostat=close_status);opened=.false.;if(io_status==0)io_status=close_status
    endif
    call sync_io(io_status,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    if(rank==0)call atomic_rename(temporary_path,trim(path),io_status)
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    call clean;ok=.true.;return
900 message='occupied checkpoint MPI stream failed';if(rank==0.and.opened)close(unit);call clean;return
910 message='occupied checkpoint publication failed';if(rank==0.and.opened)close(unit);call clean;return
#else
    ok=.false.;message='occupied checkpoint requires MPI';payload_fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine agree_int(value)
      integer,intent(in)::value;integer::lo,hi
      call MPI_Allreduce(value,lo,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,hi,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.lo/=hi)bad=1
    end subroutine
    subroutine agree_i64(value)
      integer(int64),intent(in)::value
      call MPI_Allreduce(value,minimum_i,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,maximum_i,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.minimum_i/=maximum_i)bad=1
    end subroutine
    subroutine agree_real(value)
      real(real64),intent(in)::value
      call MPI_Allreduce(value,minimum_r,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,maximum_r,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.transfer(minimum_r,bits)/=transfer(maximum_r,bits))bad=1
    end subroutine
    subroutine occupied_hash(value)
      integer(int64),intent(in)::value;payload_fingerprint=mix_hash(payload_fingerprint,value)
    end subroutine
    subroutine clean
      if(allocated(counts))deallocate(counts);if(allocated(owners))deallocate(owners)
      if(allocated(positions))deallocate(positions);if(allocated(row_values))deallocate(row_values)
    end subroutine
#endif
  end subroutine write_rt_dg_hybrid_occupied_checkpoint

  subroutine read_rt_dg_hybrid_occupied_checkpoint(comm,path,expected_catalog,expected_basis,expected_operator,&
      expected_state,expected_provenance,expected_occupations,maximum_scf_residual,global_count,row_ids,coefficients,&
      occupations,eigenvalues,scf_receipts,payload_fingerprint,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::path
    integer(int64),intent(in)::expected_catalog,expected_basis,expected_operator,expected_state,expected_provenance(6)
    real(real64),intent(in)::expected_occupations(:),maximum_scf_residual
    integer,intent(out)::global_count
    integer(int64),allocatable,intent(out)::row_ids(:)
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    real(real64),allocatable,intent(out)::occupations(:),eigenvalues(:)
    real(real64),intent(out)::scf_receipts(5)
    integer(int64),intent(out)::payload_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,unit,io_status,bad,global_bad,noccupied,version,nlocal,row,position,j,status
    integer(int64)::catalog,basis,provenance(6),operator_receipt,state,stored_fingerprint,bits,lo_i,hi_i,file_size
    real(real64)::lo_r,hi_r
    complex(real64),allocatable::row_values(:)
    character(16)::magic,probe
    logical::opened
    ok=.false.;message='';payload_fingerprint=0_int64;global_count=0;scf_receipts=0d0;opened=.false.;io_status=-1
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    call validate_path(path,probe,comm,ierr);if(ierr/=MPI_SUCCESS)then;message='inconsistent occupied checkpoint path';return;endif
    bad=0
    call agree_expected(expected_catalog);call agree_expected(expected_basis);call agree_expected(expected_operator);call agree_expected(expected_state)
    do j=1,6;call agree_expected(expected_provenance(j));enddo
    call agree_expected_count(size(expected_occupations))
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='rank-disagreeing occupied checkpoint expectation';return;endif
    call agree_expected_real(maximum_scf_residual)
    if(.not.ieee_is_finite(maximum_scf_residual).or.maximum_scf_residual<=0d0)bad=1
    do j=1,size(expected_occupations);call agree_expected_real(expected_occupations(j));enddo
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='rank-disagreeing occupied checkpoint expectation';return;endif
    if(rank==0)then
      open(newunit=unit,file=trim(path),status='old',access='stream',form='unformatted',action='read',iostat=io_status);opened=io_status==0
      if(io_status==0)inquire(unit=unit,size=file_size,iostat=io_status)
      if(io_status==0)read(unit,iostat=io_status)magic,version,global_count,noccupied,catalog,basis,provenance,operator_receipt,state
    endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    call MPI_Bcast(magic,len(magic),MPI_CHARACTER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(version,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(global_count,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(noccupied,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(catalog,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(basis,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(operator_receipt,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(state,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(provenance,6,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(file_size,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    bad=0
    if(magic/=occupied_magic.or.version/=occupied_version.or.global_count<1.or.noccupied<1)bad=1
    if(catalog/=expected_catalog.or.basis/=expected_basis.or.any(provenance/=expected_provenance).or.&
      operator_receipt/=expected_operator.or.state/=expected_state)bad=1
    if(global_count>huge(0)-nproc.or.noccupied>huge(0)/max(1,global_count))bad=1
    if(file_size<1_int64.or.int(noccupied,int64)>file_size/16_int64.or.&
      int(global_count,int64)>file_size/max(16_int64,16_int64*int(noccupied,int64)))bad=1
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='stale or incompatible occupied checkpoint';goto 920;endif
    nlocal=(global_count+nproc-1-rank)/nproc
    allocate(row_ids(nlocal),coefficients(nlocal,noccupied),occupations(noccupied),eigenvalues(noccupied),row_values(noccupied),stat=status)
    bad=merge(0,1,status==0);call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate occupied checkpoint output';goto 920;endif
    if(rank==0)read(unit,iostat=io_status)occupations,eigenvalues,scf_receipts
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    call MPI_Bcast(occupations,noccupied,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(eigenvalues,noccupied,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    call MPI_Bcast(scf_receipts,5,MPI_DOUBLE_PRECISION,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    bad=merge(0,1,all(ieee_is_finite(occupations)).and.all(occupations>=0d0).and.&
      all(ieee_is_finite(eigenvalues)).and.all(ieee_is_finite(scf_receipts)).and.all(scf_receipts>=0d0).and.&
      all(scf_receipts<=maximum_scf_residual).and.size(expected_occupations)==noccupied)
    if(bad==0)then;if(any(occupations/=expected_occupations))bad=1;endif
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.global_bad/=0)goto 930
    payload_fingerprint=catalog;call read_hash(basis);do j=1,6;call read_hash(provenance(j));enddo
    call read_hash(operator_receipt);call read_hash(state)
    call read_hash(int(global_count,int64));call read_hash(int(noccupied,int64))
    do j=1,noccupied;call read_hash(transfer(occupations(j),bits));call read_hash(transfer(eigenvalues(j),bits));enddo
    do j=1,5;call read_hash(transfer(scf_receipts(j),bits));enddo
    position=0
    do row=1,global_count
      if(rank==0)read(unit,iostat=io_status)row_values
      call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
      call MPI_Bcast(row_values,noccupied,MPI_DOUBLE_COMPLEX,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
      if(.not.finite_vector(row_values))goto 930
      call read_hash(int(row,int64));do j=1,noccupied;call read_hash(transfer(real(row_values(j)),bits));call read_hash(transfer(aimag(row_values(j)),bits));enddo
      if(mod(row-1,nproc)==rank)then;position=position+1;row_ids(position)=row;coefficients(position,:)=row_values;endif
    enddo
    if(payload_fingerprint==0_int64)payload_fingerprint=1_int64
    if(rank==0)read(unit,iostat=io_status)stored_fingerprint
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS.or.io_status/=0)goto 910
    call MPI_Bcast(stored_fingerprint,1,MPI_INTEGER8,0,comm,ierr);if(ierr/=MPI_SUCCESS)goto 900
    if(rank==0)then;close(unit,iostat=io_status);opened=.false.;endif
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.io_status/=0.or.stored_fingerprint/=payload_fingerprint)goto 930
    if(allocated(row_values))deallocate(row_values);ok=.true.;return
900 message='occupied checkpoint MPI read failed';goto 920
910 message='occupied checkpoint file read failed';goto 920
930 message='corrupt occupied checkpoint payload'
920 if(rank==0.and.opened)close(unit);call cleanup_output;return
#else
    ok=.false.;message='occupied checkpoint requires MPI';payload_fingerprint=0_int64;global_count=0;scf_receipts=0d0
#endif
  contains
#ifdef USE_MPI
    subroutine agree_expected(value)
      integer(int64),intent(in)::value
      call MPI_Allreduce(value,lo_i,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,hi_i,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.lo_i/=hi_i.or.value==0_int64)bad=1
    end subroutine
    subroutine agree_expected_count(value)
      integer,intent(in)::value;integer::lo,hi
      call MPI_Allreduce(value,lo,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,hi,1,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.lo/=hi)bad=1
    end subroutine
    subroutine agree_expected_real(value)
      real(real64),intent(in)::value
      call MPI_Allreduce(value,lo_r,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)then;bad=1;return;endif
      call MPI_Allreduce(value,hi_r,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.transfer(lo_r,bits)/=transfer(hi_r,bits))bad=1
    end subroutine
    subroutine read_hash(value)
      integer(int64),intent(in)::value;payload_fingerprint=mix_hash(payload_fingerprint,value)
    end subroutine
    subroutine cleanup_output
      if(allocated(row_ids))deallocate(row_ids);if(allocated(coefficients))deallocate(coefficients)
      if(allocated(occupations))deallocate(occupations);if(allocated(eigenvalues))deallocate(eigenvalues)
      if(allocated(row_values))deallocate(row_values);global_count=0;payload_fingerprint=0_int64
    end subroutine
#endif
  end subroutine read_rt_dg_hybrid_occupied_checkpoint

  subroutine validate_path(path,probe,comm,ierr)
    character(*),intent(in)::path;character(16),intent(out)::probe;integer,intent(in)::comm;integer,intent(out)::ierr
#ifdef USE_MPI
    integer::i,status
    integer(int64)::local_hash,minimum_hash,maximum_hash
    local_hash=int(len_trim(path),int64)
    do i=1,len_trim(path);local_hash=ieor(ishftc(local_hash,7),int(iachar(path(i:i)),int64));enddo
    call MPI_Allreduce(local_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS.and.minimum_hash/=maximum_hash)ierr=MPI_ERR_OTHER
    probe=''
#else
    probe='';ierr=1
#endif
  end subroutine validate_path
  subroutine sync_io(io_status,comm,ierr)
    integer,intent(inout)::io_status;integer,intent(in)::comm;integer,intent(out)::ierr
#ifdef USE_MPI
    call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,ierr)
#else
    ierr=1
#endif
  end subroutine sync_io
  subroutine atomic_rename(old_path,new_path,status)
    character(*),intent(in)::old_path,new_path;integer,intent(out)::status
    character(c_char),allocatable::old_c(:),new_c(:)
    integer::i
    allocate(old_c(len_trim(old_path)+1),new_c(len_trim(new_path)+1))
    do i=1,len_trim(old_path);old_c(i)=old_path(i:i);enddo;old_c(size(old_c))=c_null_char
    do i=1,len_trim(new_path);new_c(i)=new_path(i:i);enddo;new_c(size(new_c))=c_null_char
    status=int(c_rename(old_c,new_c))
  end subroutine atomic_rename
  logical function finite_vector(values)
    complex(real64),intent(in)::values(:)
    finite_vector=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_vector
  logical function finite_matrix(values)
    complex(real64),intent(in)::values(:,:)
    finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_matrix
  logical function strictly_increasing(values)
    integer,intent(in)::values(:);integer::i
    strictly_increasing=.true.
    do i=2,size(values);if(values(i)<=values(i-1))then;strictly_increasing=.false.;return;endif;enddo
  end function strictly_increasing
  logical function valid_packet_activity(packet_ids,active_rows)
    integer,intent(in)::packet_ids(:);logical,intent(in)::active_rows(:)
    integer,allocatable::states(:);integer::i,status
    valid_packet_activity=.false.
    if(size(packet_ids)/=size(active_rows).or.any(packet_ids<1).or.any(packet_ids>size(packet_ids)))return
    allocate(states(size(packet_ids)),stat=status);if(status/=0)return;states=-1
    do i=1,size(packet_ids)
      if(states(packet_ids(i))<0)then
        states(packet_ids(i))=merge(1,0,active_rows(i))
      else if(states(packet_ids(i))/=merge(1,0,active_rows(i)))then
        deallocate(states);return
      endif
    enddo
    deallocate(states);valid_packet_activity=.true.
  end function valid_packet_activity
  pure integer(int64) function mix_hash(seed,value)
    integer(int64),intent(in)::seed,value
    mix_hash=ieor(ishftc(seed,9),value)
  end function mix_hash
  function int64_string(value) result(text)
    integer(int64),intent(in)::value;character(32)::text
    write(text,'(i0)')value
  end function int64_string
end module rt_dg_hybrid_checkpoint
