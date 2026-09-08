#include "config.h"
module dg_hybrid_full_cell_operator_adapter
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  ! Diagnostic full-cell oracle used before the scalable neighbor assembler is
  ! enabled. It bounds memory by materializing basis tiles, but intentionally
  ! performs global validation and is not a production sparse-assembly route.
  abstract interface
    subroutine dg_hybrid_basis_provider(first_column,column_count,tile_values,ok)
      import real64
      integer,intent(in)::first_column,column_count
      complex(real64),intent(out)::tile_values(:,:)
      logical,intent(out)::ok
    end subroutine dg_hybrid_basis_provider
    subroutine dg_hybrid_tile_operator(tile_in,tile_out,ok)
      import real64
      complex(real64),intent(in)::tile_in(:,:)
      complex(real64),intent(out)::tile_out(:,:)
      logical,intent(out)::ok
    end subroutine dg_hybrid_tile_operator
  end interface
  public::project_dg_hybrid_full_cell_sparse_operators
contains
  subroutine project_dg_hybrid_full_cell_sparse_operators(comm,global_spatial_count,global_basis_count,spatial_ids,weights,&
      coordinates,materialize_basis,row_ids,row_offsets,column_ids,expected_metric_values,tile_width,apply_tile,selection_fingerprint,&
      window_fingerprint,packet_fingerprint,complement_fingerprint,metric_fingerprint,position_convention_fingerprint,&
      tolerance,operators,&
      persistent_bytes,transient_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_spatial_count,global_basis_count,tile_width
    integer(int64),intent(in)::spatial_ids(:),row_ids(:)
    real(real64),intent(in)::weights(:),coordinates(:,:),tolerance
    integer,intent(in)::row_offsets(:),column_ids(:)
    complex(real64),intent(in)::expected_metric_values(:)
    procedure(dg_hybrid_basis_provider)::materialize_basis
    procedure(dg_hybrid_tile_operator)::apply_tile
    integer(int64),intent(in)::selection_fingerprint,window_fingerprint,packet_fingerprint,&
      complement_fingerprint,metric_fingerprint
    ! coordinates must already use the caller's chosen periodic branch and
    ! local-origin convention; this opaque receipt binds that convention.
    integer(int64),intent(in)::position_convention_fingerprint
    type(s_dg_hybrid_sparse_operators),intent(out)::operators
    integer(int64),intent(out)::persistent_bytes,transient_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,ierr,nlocal,nbasis,nowned,nnz,i,j,k,row,column,j0,j1,width,root,position
    integer::local_bad,global_bad,minimum_integer,maximum_integer,allocation_status,max_degree,degree
    integer,allocatable::spatial_count(:),row_count(:),row_owner(:),&
      row_position(:),column_buffer(:)
    integer(int64)::bits,minimum_bits,maximum_bits,persistent_complex,persistent_integer,&
      persistent_int64,transient_complex,transient_integer,quantized
    complex(real64),allocatable::tile_in(:,:),tile_out(:,:),bra_tile(:,:),local_receipt(:),reduced_receipt(:),&
      remote_metric(:),remote_hamiltonian(:),remote_position(:,:)
    real(real64)::local_scale,global_scale,safe_scale,hermitian_defect,operator_scale,quantization_limit
    logical::callback_ok
    ok=.false.;message='';persistent_bytes=0_int64;transient_peak_bytes=0_int64;fingerprint=0_int64
    local_bad=0;nlocal=size(spatial_ids);nbasis=global_basis_count;nowned=size(row_ids);nnz=size(column_ids)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call agree_integer(global_spatial_count,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent hybrid full-cell grid extent';return;endif
    call agree_integer(nbasis,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent hybrid basis extent';return;endif
    call agree_integer(tile_width,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent hybrid operator tile width';return;endif
    bits=transfer(tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='inconsistent hybrid operator tolerance';return;endif
    call agree_receipt(selection_fingerprint,'selection',local_bad)
    call agree_receipt(window_fingerprint,'window',local_bad)
    call agree_receipt(packet_fingerprint,'packet',local_bad)
    call agree_receipt(complement_fingerprint,'complement',local_bad)
    call agree_receipt(metric_fingerprint,'metric',local_bad)
    call agree_receipt(position_convention_fingerprint,'position convention',local_bad)
    if(local_bad/=0)return
    if(global_spatial_count<1.or.nbasis<1.or.tile_width<1)local_bad=1
    if(size(weights)/=nlocal.or.any(shape(coordinates)/=[3,nlocal]))local_bad=1
    if(nowned==huge(0).or.nnz==huge(0))local_bad=1
    if(local_bad==0)then
      if(size(row_offsets)/=nowned+1.or.size(row_offsets)<1)local_bad=1
    endif
    if(size(expected_metric_values)/=nnz)local_bad=1
    if(any(spatial_ids<1_int64).or.any(spatial_ids>int(global_spatial_count,int64)))local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(nbasis,int64)))local_bad=1
    if(.not.ieee_is_finite(tolerance).or..not.all(ieee_is_finite(weights)).or.&
      .not.all(ieee_is_finite(coordinates)).or.&
      .not.finite_vector(expected_metric_values))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid sparse operator shape or finite contract';return;endif
    if(tolerance<1d-15.or.tolerance>1d-2.or.any(weights<=0d0))local_bad=1
    if(row_offsets(1)/=1.or.row_offsets(nowned+1)/=nnz+1.or.any(row_offsets(2:)<row_offsets(:nowned)))local_bad=1
    max_degree=0
    do i=1,nowned
      max_degree=max(max_degree,row_offsets(i+1)-row_offsets(i))
      do k=row_offsets(i),row_offsets(i+1)-1
        if(column_ids(k)<1.or.column_ids(k)>nbasis)local_bad=1
        if(k>row_offsets(i))then;if(column_ids(k)<=column_ids(k-1))local_bad=1;endif
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,max_degree,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    if(int(tile_width,int64)*int(nlocal,int64)>int(huge(0),int64))local_bad=1
    if(3_int64*int(nnz,int64)>int(huge(0),int64))local_bad=1
    if(3_int64*int(nbasis,int64)>int(huge(0),int64))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid sparse operator graph';return;endif
    ! Receipts cover storage owned by this adapter. Caller-owned grid metadata and
    ! callback-internal hpsi workspace are deliberately outside this boundary.
    persistent_complex=0_int64;persistent_integer=0_int64;persistent_int64=int(nowned,int64)
    transient_complex=0_int64;transient_integer=0_int64
    call add_count_product(persistent_complex,5_int64,int(nnz,int64),local_bad)
    call add_count(persistent_integer,int(nowned,int64),local_bad)
    call add_count(persistent_integer,1_int64,local_bad)
    call add_count(persistent_integer,int(nnz,int64),local_bad)
    call add_count_product(transient_complex,2_int64*int(tile_width,int64),int(nlocal,int64),local_bad)
    call add_count(transient_complex,int(nlocal,int64),local_bad)
    call add_count(transient_complex,10_int64,local_bad)
    call add_count_product(transient_complex,5_int64,int(nbasis,int64),local_bad)
    call add_count(transient_integer,int(global_spatial_count,int64),local_bad)
    call add_count_product(transient_integer,3_int64,int(nbasis,int64),local_bad)
    call add_count(transient_integer,int(max_degree,int64),local_bad)
    if(.not.bytes_from_counts(persistent_complex,persistent_integer,persistent_int64,persistent_bytes))local_bad=1
    if(.not.bytes_from_counts(transient_complex,transient_integer,0_int64,transient_peak_bytes))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='hybrid sparse operator byte receipt overflow';return;endif
    allocate(spatial_count(global_spatial_count),row_count(nbasis),row_owner(nbasis),row_position(nbasis),&
      column_buffer(max(1,max_degree)),tile_in(tile_width,nlocal),tile_out(tile_width,nlocal),&
      bra_tile(1,nlocal),&
      local_receipt(5),reduced_receipt(5),remote_metric(nbasis),remote_hamiltonian(nbasis),&
      remote_position(3,nbasis),operators%owned_row_ids(nowned),operators%row_offsets(nowned+1),&
      operators%column_ids(nnz),operators%metric_values(nnz),operators%hamiltonian_values(nnz),&
      operators%position_values(3,nnz),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate hybrid sparse operator workspace';return;endif
    spatial_count=0
    do i=1,nlocal
      spatial_count(int(spatial_ids(i)))=spatial_count(int(spatial_ids(i)))+1
    enddo
    row_count=0;row_owner=-1;row_position=0
    do i=1,nowned
      row=int(row_ids(i));row_count(row)=row_count(row)+1;row_owner(row)=rank;row_position(row)=i
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,spatial_count,global_spatial_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid spatial ownership reduction failed';return;endif
    call reduce_ownership(row_count,row_owner,row_position,nbasis)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid basis ownership reduction failed';return;endif
    if(any(spatial_count/=1).or.any(row_count/=1))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='duplicate or missing hybrid operator owner';return;endif
    local_scale=0d0
    if(nlocal>0)local_scale=max(maxval(abs(coordinates)),maxval(weights))
    call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    safe_scale=sqrt(sqrt(huge(1d0)))/(16d0*real(max(global_spatial_count,nbasis),real64))
    if(ierr/=MPI_SUCCESS.or.global_scale>safe_scale)then;call cleanup();message='unsafe hybrid operator input magnitude';return;endif
    operators%metric_values=(0d0,0d0);operators%hamiltonian_values=(0d0,0d0);operators%position_values=(0d0,0d0)
    j0=1
    do
      width=min(tile_width,nbasis-j0+1);j1=j0+(width-1)
      call materialize_basis(j0,width,tile_in(1:width,:),callback_ok)
      local_bad=merge(0,1,callback_ok.and.finite_matrix(tile_in(1:width,:)))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid basis tile materialization failed';return;endif
      local_scale=0d0;if(nlocal>0)local_scale=maxval(abs(tile_in(1:width,:)))
      call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_scale>safe_scale)then;call cleanup();message='unsafe hybrid basis tile magnitude';return;endif
      call apply_tile(tile_in(1:width,:),tile_out(1:width,:),callback_ok)
      local_bad=merge(0,1,callback_ok.and.finite_matrix(tile_out(1:width,:)))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid Hamiltonian tile callback failed';return;endif
      local_scale=0d0;if(nlocal>0)local_scale=maxval(abs(tile_out(1:width,:)))
      call MPI_Allreduce(local_scale,operator_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.operator_scale>safe_scale)then;call cleanup();message='unsafe hybrid Hamiltonian tile magnitude';return;endif
      do row=1,nbasis
        call materialize_basis(row,1,bra_tile,callback_ok)
        local_bad=merge(0,1,callback_ok.and.finite_matrix(bra_tile))
        call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid bra tile materialization failed';return;endif
        local_scale=0d0;if(nlocal>0)local_scale=maxval(abs(bra_tile))
        call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.global_scale>safe_scale)then;call cleanup();message='unsafe hybrid bra tile magnitude';return;endif
        root=row_owner(row);degree=0;column_buffer=0
        if(rank==root)then
          position=row_position(row);degree=row_offsets(position+1)-row_offsets(position)
          column_buffer(1:degree)=column_ids(row_offsets(position):row_offsets(position+1)-1)
        endif
        call MPI_Bcast(degree,1,MPI_INTEGER,root,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid graph degree broadcast failed';return;endif
        call MPI_Bcast(column_buffer,degree,MPI_INTEGER,root,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid graph row broadcast failed';return;endif
        do k=1,degree
          column=column_buffer(k);if(column<j0.or.column>j1)cycle
          local_receipt(1)=sum(weights*conjg(bra_tile(1,:))*tile_in(column-j0+1,:))
          local_receipt(2)=sum(weights*conjg(bra_tile(1,:))*tile_out(column-j0+1,:))
          do i=1,3
            local_receipt(2+i)=sum(weights*conjg(bra_tile(1,:))*coordinates(i,:)*tile_in(column-j0+1,:))
          enddo
          call MPI_Reduce(local_receipt,reduced_receipt,5,MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
          if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid sparse edge reduction failed';return;endif
          if(rank==root)then
            position=row_offsets(row_position(row))+k-1
            operators%metric_values(position)=reduced_receipt(1)
            operators%hamiltonian_values(position)=reduced_receipt(2)
            operators%position_values(:,position)=reduced_receipt(3:5)
          endif
        enddo
      enddo
      if(j1==nbasis)exit
      j0=j1+1
    enddo
    local_scale=0d0
    if(nnz>0)local_scale=maxval(abs(operators%metric_values-expected_metric_values))
    call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_scale>100d0*tolerance)then
      call cleanup();message='recomputed hybrid metric does not match its provenance payload';return
    endif
    hermitian_defect=0d0;operator_scale=1d0
    do row=1,nbasis
      call broadcast_operator_row(row,remote_metric,remote_hamiltonian,remote_position,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid operator validation broadcast failed';return;endif
      operator_scale=max(operator_scale,maxval(abs(remote_metric)),maxval(abs(remote_hamiltonian)),&
        maxval(abs(remote_position)))
      do i=1,nowned
        k=find_column(row,i,row_offsets,column_ids)
        if(k>0)then
          hermitian_defect=max(hermitian_defect,abs(operators%metric_values(k)-conjg(remote_metric(int(row_ids(i))))))
          hermitian_defect=max(hermitian_defect,abs(operators%hamiltonian_values(k)-conjg(remote_hamiltonian(int(row_ids(i))))))
          do j=1,3
            hermitian_defect=max(hermitian_defect,&
              abs(operators%position_values(j,k)-conjg(remote_position(j,int(row_ids(i))))))
          enddo
        else
          hermitian_defect=max(hermitian_defect,abs(remote_metric(int(row_ids(i)))),&
            abs(remote_hamiltonian(int(row_ids(i)))),maxval(abs(remote_position(:,int(row_ids(i))))))
        endif
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,hermitian_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.hermitian_defect>100d0*tolerance*operator_scale)then
      call cleanup();message='hybrid sparse S/H/Z is not Hermitian or graph-complete';return
    endif
    quantization_limit=0.25d0*real(huge(0_int64),real64)*100d0*tolerance
    if(operator_scale>quantization_limit)then;call cleanup();message='hybrid operator fingerprint range is unsafe';return;endif
    fingerprint=selection_fingerprint
    fingerprint=ieor(ishftc(fingerprint,9),window_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,9),packet_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,9),complement_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,9),metric_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,9),position_convention_fingerprint)
    do row=1,nbasis
      call broadcast_operator_row(row,remote_metric,remote_hamiltonian,remote_position,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid operator fingerprint broadcast failed';return;endif
      fingerprint=ieor(ishftc(fingerprint,9),int(row,int64))
      do j=1,nbasis
        if(abs(remote_metric(j))+abs(remote_hamiltonian(j))+sum(abs(remote_position(:,j)))==0d0)cycle
        fingerprint=ieor(ishftc(fingerprint,9),int(j,int64))
        call hash_complex(remote_metric(j));call hash_complex(remote_hamiltonian(j))
        do i=1,3;call hash_complex(remote_position(i,j));enddo
      enddo
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    operators%valid=.true.;operators%global_count=nbasis;operators%owned_row_ids=row_ids
    operators%row_offsets=row_offsets;operators%column_ids=column_ids
    operators%selection_fingerprint=selection_fingerprint;operators%window_fingerprint=window_fingerprint
    operators%packet_fingerprint=packet_fingerprint;operators%complement_fingerprint=complement_fingerprint
    operators%metric_fingerprint=metric_fingerprint;operators%fingerprint=fingerprint
    operators%position_convention_fingerprint=position_convention_fingerprint
    operators%persistent_bytes=persistent_bytes;operators%transient_peak_bytes=transient_peak_bytes;ok=.true.
#else
    ok=.false.;message='hybrid full-cell sparse operator adapter requires MPI'
    persistent_bytes=0_int64;transient_peak_bytes=0_int64;fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine agree_receipt(value,label,bad)
      integer(int64),intent(in)::value;character(*),intent(in)::label;integer,intent(inout)::bad
      call agree_int64(value,minimum_bits,maximum_bits,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.value==0_int64)then
        message='invalid or inconsistent hybrid '//trim(label)//' provenance';bad=1
      endif
    end subroutine agree_receipt
    subroutine reduce_ownership(counts,owners,positions,count)
      integer,intent(inout)::counts(:),owners(:),positions(:);integer,intent(in)::count
      call MPI_Allreduce(MPI_IN_PLACE,counts,count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(MPI_IN_PLACE,owners,count,MPI_INTEGER,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(MPI_IN_PLACE,positions,count,MPI_INTEGER,MPI_MAX,comm,ierr)
    end subroutine reduce_ownership
    subroutine broadcast_operator_row(global_row,metric_row,hamiltonian_row,position_row,status)
      integer,intent(in)::global_row;integer,intent(out)::status
      complex(real64),intent(out)::metric_row(:),hamiltonian_row(:),position_row(:,:)
      integer::local_position_index,edge
      metric_row=(0d0,0d0);hamiltonian_row=(0d0,0d0);position_row=(0d0,0d0);root=row_owner(global_row)
      if(rank==root)then
        local_position_index=row_position(global_row)
        do edge=row_offsets(local_position_index),row_offsets(local_position_index+1)-1
          metric_row(column_ids(edge))=operators%metric_values(edge)
          hamiltonian_row(column_ids(edge))=operators%hamiltonian_values(edge)
          position_row(:,column_ids(edge))=operators%position_values(:,edge)
        enddo
      endif
      call MPI_Bcast(metric_row,nbasis,MPI_DOUBLE_COMPLEX,root,comm,status);if(status/=MPI_SUCCESS)return
      call MPI_Bcast(hamiltonian_row,nbasis,MPI_DOUBLE_COMPLEX,root,comm,status);if(status/=MPI_SUCCESS)return
      call MPI_Bcast(position_row,3*nbasis,MPI_DOUBLE_COMPLEX,root,comm,status)
    end subroutine broadcast_operator_row
    subroutine hash_complex(value)
      complex(real64),intent(in)::value
      quantized=nint(real(value)/(100d0*tolerance),int64);fingerprint=ieor(ishftc(fingerprint,9),quantized)
      quantized=nint(aimag(value)/(100d0*tolerance),int64);fingerprint=ieor(ishftc(fingerprint,9),quantized)
    end subroutine hash_complex
    subroutine cleanup()
      if(allocated(spatial_count))deallocate(spatial_count)
      if(allocated(row_count))deallocate(row_count)
      if(allocated(row_owner))deallocate(row_owner)
      if(allocated(row_position))deallocate(row_position)
      if(allocated(column_buffer))deallocate(column_buffer)
      if(allocated(tile_in))deallocate(tile_in)
      if(allocated(tile_out))deallocate(tile_out)
      if(allocated(bra_tile))deallocate(bra_tile)
      if(allocated(local_receipt))deallocate(local_receipt)
      if(allocated(reduced_receipt))deallocate(reduced_receipt)
      if(allocated(remote_metric))deallocate(remote_metric)
      if(allocated(remote_hamiltonian))deallocate(remote_hamiltonian)
      if(allocated(remote_position))deallocate(remote_position)
      if(allocated(operators%owned_row_ids))deallocate(operators%owned_row_ids)
      if(allocated(operators%row_offsets))deallocate(operators%row_offsets)
      if(allocated(operators%column_ids))deallocate(operators%column_ids)
      if(allocated(operators%metric_values))deallocate(operators%metric_values)
      if(allocated(operators%hamiltonian_values))deallocate(operators%hamiltonian_values)
      if(allocated(operators%position_values))deallocate(operators%position_values)
      operators%valid=.false.
    end subroutine cleanup
#endif
  end subroutine project_dg_hybrid_full_cell_sparse_operators

  integer function find_column(column,row_position,row_offsets,column_ids)
    integer,intent(in)::column,row_position,row_offsets(:),column_ids(:);integer::edge
    find_column=0
    do edge=row_offsets(row_position),row_offsets(row_position+1)-1
      if(column_ids(edge)==column)then;find_column=edge;return;endif
    enddo
  end function find_column
  logical function finite_matrix(values)
    complex(real64),intent(in)::values(:,:)
    finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_matrix
  logical function finite_vector(values)
    complex(real64),intent(in)::values(:)
    finite_vector=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_vector
  logical function bytes_from_counts(complex_count,integer_count,int64_count,bytes)
    integer(int64),intent(in)::complex_count,integer_count,int64_count;integer(int64),intent(out)::bytes
    bytes=0_int64;bytes_from_counts=.false.
    if(complex_count<0_int64.or.integer_count<0_int64.or.int64_count<0_int64)return
    if(complex_count>huge(bytes)/16_int64.or.integer_count>huge(bytes)/4_int64.or.&
      int64_count>huge(bytes)/8_int64)return
    bytes=16_int64*complex_count
    if(4_int64*integer_count>huge(bytes)-bytes)return
    bytes=bytes+4_int64*integer_count
    if(8_int64*int64_count>huge(bytes)-bytes)return
    bytes=bytes+8_int64*int64_count;bytes_from_counts=.true.
  end function bytes_from_counts
  subroutine add_count(total,value,bad)
    integer(int64),intent(inout)::total
    integer(int64),intent(in)::value
    integer,intent(inout)::bad
    if(bad/=0)return
    if(value<0_int64.or.total>huge(total)-value)then
      bad=1
    else
      total=total+value
    endif
  end subroutine add_count
  subroutine add_count_product(total,left,right,bad)
    integer(int64),intent(inout)::total
    integer(int64),intent(in)::left,right
    integer,intent(inout)::bad
    if(bad/=0)return
    if(left<0_int64.or.right<0_int64)then
      bad=1;return
    endif
    if(left/=0_int64)then
      if(right>huge(total)/left)then;bad=1;return;endif
    endif
    call add_count(total,left*right,bad)
  end subroutine add_count_product
#ifdef USE_MPI
  subroutine agree_integer(value,minimum_value,maximum_value,comm,ierr)
    integer,intent(in)::value,comm;integer,intent(out)::minimum_value,maximum_value,ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine agree_integer
  subroutine agree_int64(value,minimum_value,maximum_value,comm,ierr)
    integer(int64),intent(in)::value;integer,intent(in)::comm
    integer(int64),intent(out)::minimum_value,maximum_value;integer,intent(out)::ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
  end subroutine agree_int64
#endif
end module dg_hybrid_full_cell_operator_adapter
