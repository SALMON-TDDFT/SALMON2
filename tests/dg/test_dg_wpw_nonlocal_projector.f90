program test_dg_wpw_nonlocal_projector
  use dg_wpw_nonlocal_projector,only:s_dg_wpw_projector_overlap,&
    assemble_dg_wpw_nonlocal_blocks,build_dg_wpw_projector_overlap_partials
  implicit none
  type(s_dg_wpw_projector_overlap)::records(8),duplicate(9),bad(8)
  type(s_dg_wpw_projector_overlap),allocatable::partials(:)
  integer::projector_ids(3),rows(6),columns(6),info,i,j
  real(8)::weights(3)
  complex(8)::values(6),oracle(6),dense_overlap(4,3)
  complex(8)::w_values(1,2),p_values(1,2),projector_values(2)

  w_values=reshape([(2d0,1d0),(4d0,-1d0)],[1,2])
  p_values=reshape([(3d0,0.5d0),(5d0,-0.5d0)],[1,2])
  projector_values=[(0.5d0,0.25d0),(-1d0,0.5d0)]
  call build_dg_wpw_projector_overlap_partials([7,9],[0.2d0,0.3d0],[2],w_values,[4],p_values,&
    [7,9],[11,11],projector_values,100,partials,info)
  if(info/=0.or.size(partials)/=2)error stop 'projector quadrature partial construction failed'
  if(partials(1)%basis_id/=2.or.partials(1)%projector_id/=11.or.&
     abs(partials(1)%value-sum(conjg(projector_values)*w_values(1,:)*[0.2d0,0.3d0]))>1d-13)&
    error stop 'W projector quadrature uses the wrong convention'
  if(partials(2)%basis_id/=104.or.partials(2)%projector_id/=11.or.&
     abs(partials(2)%value-sum(conjg(projector_values)*p_values(1,:)*[0.2d0,0.3d0]))>1d-13)&
    error stop 'P projector namespace or quadrature is wrong'

  projector_ids=[11,13,17];weights=[2d0,-0.5d0,1.25d0]
  dense_overlap=reshape([&
    (1d0,0.5d0),(0.2d0,-0.1d0),(0d0,0d0),(0.7d0,0.3d0),&
    (0d0,0d0),(0.4d0,0.2d0),(0.8d0,-0.6d0),(0d0,0d0),&
    (-0.3d0,0.1d0),(0d0,0d0),(0.5d0,0.25d0),(0.9d0,-0.2d0)],shape(dense_overlap))
  records=[&
    s_dg_wpw_projector_overlap(101,11,dense_overlap(1,1)),&
    s_dg_wpw_projector_overlap(101,17,dense_overlap(1,3)),&
    s_dg_wpw_projector_overlap(102,11,dense_overlap(2,1)),&
    s_dg_wpw_projector_overlap(102,13,dense_overlap(2,2)),&
    s_dg_wpw_projector_overlap(201,13,dense_overlap(3,2)),&
    s_dg_wpw_projector_overlap(201,17,dense_overlap(3,3)),&
    s_dg_wpw_projector_overlap(202,11,dense_overlap(4,1)),&
    s_dg_wpw_projector_overlap(202,17,dense_overlap(4,3))]
  rows=[101,102,101,201,202,201];columns=[102,101,201,101,201,202]
  do i=1,size(rows)
    oracle(i)=(0d0,0d0)
    do j=1,3
      oracle(i)=oracle(i)+conjg(dense_overlap(basis_pos(rows(i)),j))*weights(j)*&
        dense_overlap(basis_pos(columns(i)),j)
    enddo
  enddo
  call assemble_dg_wpw_nonlocal_blocks(records,projector_ids,weights,rows,columns,values,info)
  if(info/=0.or.maxval(abs(values-oracle))>1d-13)error stop 'projector outer-product oracle mismatch'
  call assemble_dg_wpw_nonlocal_blocks(records,projector_ids,weights,[999],[101],values(1:1),info)
  if(info/=0.or.abs(values(1))>1d-13)error stop 'missing basis did not produce a zero block'

  duplicate(1:8)=records;duplicate(9)=records(8)
  call assemble_dg_wpw_nonlocal_blocks(duplicate,projector_ids,weights,rows,columns,values,info)
  if(info==0)error stop 'duplicate overlap record was accepted'
  bad=records;bad(1)%projector_id=19
  call assemble_dg_wpw_nonlocal_blocks(bad,projector_ids,weights,rows,columns,values,info)
  if(info==0)error stop 'unknown projector ID was accepted'
  print '(a)','PASS bounded WPW nonlocal projector outer products match dense oracle'
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
