program test_dg_wpw_nonlocal_projector_mpi
  use mpi
  use dg_wpw_nonlocal_projector,only:s_dg_wpw_projector_overlap,&
    exchange_dg_wpw_projector_overlaps,assemble_dg_wpw_nonlocal_blocks,&
    reduce_dg_wpw_p_projector_partials,reduce_dg_wpw_projector_dense,&
    validate_dg_wpw_projector_provenance,reduce_dg_wpw_projector_records,&
    validate_dg_wpw_record_counts
  implicit none
  type(s_dg_wpw_projector_overlap),allocatable::owned(:),support(:),partial(:),reduced(:)
  type(s_dg_wpw_projector_overlap),allocatable::empty_owned(:),zero_support(:),bad_partial(:)
  type(s_dg_wpw_projector_overlap),allocatable::fragment_local(:),fragment_records(:)
  integer,allocatable::asymmetric_peers(:)
  integer::rank,nrank,ierr,info,peers(1),projector_ids(3),rows(2),columns(2),i,j
  real(8)::weights(3)
  complex(8)::values(2),oracle(2),dense_overlap(4,3)
  complex(8)::dense_local(2,2),dense_sum(2,2)
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)
  if(nrank/=2)error stop 'projector exchange fixture requires two ranks'
  call validate_dg_wpw_record_counts([huge(0)/2+1],info)
  if(info==0)error stop 'overflowing packed MPI count was accepted'
  call validate_dg_wpw_record_counts([0,1],info)
  if(info/=0)error stop 'valid MPI record counts were rejected'
  if(rank==0)then;allocate(asymmetric_peers(1));asymmetric_peers=[1]
  else;allocate(asymmetric_peers(0));endif
  allocate(empty_owned(0))
  call exchange_dg_wpw_projector_overlaps(MPI_COMM_WORLD,asymmetric_peers,empty_owned,zero_support,info)
  if(info==0)error stop 'asymmetric peer schedule was accepted'
  deallocate(empty_owned)
  projector_ids=[11,13,17];weights=[2d0,-0.5d0,1.25d0]
  call validate_dg_wpw_projector_provenance(MPI_COMM_WORLD,[1,2],[4,8],&
    [1d0,merge(2d0,3d0,rank==0)],info)
  if(info==0)error stop 'mismatched projector provenance was accepted'
  dense_local=rank+1;dense_sum=0
  call reduce_dg_wpw_projector_dense(MPI_COMM_WORLD,0,dense_local,dense_sum,merge(0,9,rank==0),info)
  if(info==0)error stop 'rank-local dense overlap failure did not fail collectively'
  allocate(fragment_local(1));fragment_local=s_dg_wpw_projector_overlap(1,11,cmplx(rank+1d0,0d0,8))
  call reduce_dg_wpw_projector_records(MPI_COMM_WORLD,0,fragment_local,fragment_records,0,info)
  if(info/=0)error stop 'sparse fragment projector reduction failed'
  if(rank==0)then
    if(size(fragment_records)/=1.or.abs(fragment_records(1)%value-(3d0,0d0))>1d-13)&
      error stop 'sparse fragment projector sum is wrong'
  else if(size(fragment_records)/=0)then
    error stop 'non-root sparse fragment result is not empty'
  endif
  allocate(partial(2))
  if(rank==0)then
    partial=[s_dg_wpw_projector_overlap(1,11,(1d0,2d0)),&
      s_dg_wpw_projector_overlap(3,13,(3d0,4d0))]
  else
    partial=[s_dg_wpw_projector_overlap(1,11,(5d0,6d0)),&
      s_dg_wpw_projector_overlap(3,13,(7d0,8d0))]
  endif
  peers=[1-rank]
  call reduce_dg_wpw_p_projector_partials(MPI_COMM_WORLD,2,peers,partial,reduced,info)
  if(info/=0.or.size(reduced)/=1)error stop 'P projector partial reduction failed'
  if(rank==0)then
    if(reduced(1)%basis_id/=1.or.reduced(1)%projector_id/=11.or.&
       abs(reduced(1)%value-(6d0,8d0))>1d-13)error stop 'rank-zero P projector sum is wrong'
  else
    if(reduced(1)%basis_id/=3.or.reduced(1)%projector_id/=13.or.&
       abs(reduced(1)%value-(10d0,12d0))>1d-13)error stop 'rank-one P projector sum is wrong'
  endif
  allocate(bad_partial(1));bad_partial=s_dg_wpw_projector_overlap(2,11,(1d0,0d0))
  call reduce_dg_wpw_p_projector_partials(MPI_COMM_WORLD,2,peers,bad_partial,reduced,info,basis_offset=2)
  if(info==0)error stop 'basis ID at namespace offset was accepted'
  if(rank==0)then
    allocate(empty_owned(0))
  else
    allocate(empty_owned(1));empty_owned=s_dg_wpw_projector_overlap(301,11,(1d0,0d0))
  endif
  call exchange_dg_wpw_projector_overlaps(MPI_COMM_WORLD,peers,empty_owned,zero_support,info)
  if(info/=0.or.size(zero_support)/=1)error stop 'zero-count peer exchange failed'
  dense_overlap=reshape([&
    (1d0,0.5d0),(0.2d0,-0.1d0),(0d0,0d0),(0.7d0,0.3d0),&
    (0d0,0d0),(0.4d0,0.2d0),(0.8d0,-0.6d0),(0d0,0d0),&
    (-0.3d0,0.1d0),(0d0,0d0),(0.5d0,0.25d0),(0.9d0,-0.2d0)],shape(dense_overlap))
  if(rank==0)then
    allocate(owned(4));owned=[&
      s_dg_wpw_projector_overlap(101,11,dense_overlap(1,1)),&
      s_dg_wpw_projector_overlap(101,17,dense_overlap(1,3)),&
      s_dg_wpw_projector_overlap(102,11,dense_overlap(2,1)),&
      s_dg_wpw_projector_overlap(102,13,dense_overlap(2,2))]
    rows=[101,102];columns=[201,202];peers=[1]
  else
    allocate(owned(4));owned=[&
      s_dg_wpw_projector_overlap(201,13,dense_overlap(3,2)),&
      s_dg_wpw_projector_overlap(201,17,dense_overlap(3,3)),&
      s_dg_wpw_projector_overlap(202,11,dense_overlap(4,1)),&
      s_dg_wpw_projector_overlap(202,17,dense_overlap(4,3))]
    rows=[201,202];columns=[101,102];peers=[0]
  endif
  call exchange_dg_wpw_projector_overlaps(MPI_COMM_WORLD,peers,owned,support,info)
  if(info/=0.or.size(support)/=8)error stop 'bounded projector overlap exchange failed'
  call assemble_dg_wpw_nonlocal_blocks(support,projector_ids,weights,rows,columns,values,info)
  if(info/=0)error stop 'exchanged projector block assembly failed'
  do i=1,2
    oracle(i)=0
    do j=1,3
      oracle(i)=oracle(i)+conjg(dense_overlap(basis_pos(rows(i)),j))*weights(j)*&
        dense_overlap(basis_pos(columns(i)),j)
    enddo
  enddo
  if(maxval(abs(values-oracle))>1d-13)error stop 'exchanged projector blocks differ from dense oracle'
  if(rank==0)print '(a)','PASS bounded two-rank WPW projector overlap exchange matches dense oracle'
  call MPI_Finalize(ierr)
contains
  integer function basis_pos(id)result(pos)
    integer,intent(in)::id
    select case(id)
    case(101);pos=1
    case(102);pos=2
    case(201);pos=3
    case(202);pos=4
    case default;error stop 'unknown fixture basis ID'
    end select
  end function
end program
