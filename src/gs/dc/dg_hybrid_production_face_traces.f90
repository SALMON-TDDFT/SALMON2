#include "config.h"
module dg_hybrid_production_face_traces
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_sipg_operator,only:s_dg_hybrid_sipg_face_operator
#ifdef USE_MPI
  use mpi,only:MPI_Allreduce,MPI_Comm_rank,MPI_Comm_size,MPI_DOUBLE_COMPLEX,MPI_IN_PLACE,MPI_INTEGER,MPI_INTEGER8,&
    MPI_MAX,MPI_MIN,MPI_STATUS_IGNORE,MPI_SUCCESS,MPI_SUM,MPI_Sendrecv
#endif
  implicit none
  private

  type,public::s_dg_hybrid_production_face_trace
    logical::frozen=.false.
    integer::global_face_id=0,minus_fragment=0,plus_fragment=0
    integer::periodic_shift(3)=0
    real(real64)::canonical_normal(3)=0d0,h_normal=0d0
    integer(int64),allocatable::point_ids_minus(:),point_ids_plus(:)
    real(real64),allocatable::weights(:)
    integer,allocatable::basis_ids_minus(:),basis_ids_plus(:)
    complex(real64),allocatable::value_minus(:,:),value_plus(:,:)
    complex(real64),allocatable::derivative_minus(:,:),derivative_plus(:,:)
    integer(int64)::fingerprint=0_int64
  end type s_dg_hybrid_production_face_trace

  public::build_dg_hybrid_production_face_trace,assemble_dg_hybrid_production_face,&
    validate_dg_hybrid_production_face_collection,materialize_dg_hybrid_production_face_collection,&
    assemble_dg_hybrid_production_interface_rows,freeze_dg_hybrid_basis_directory
contains
  subroutine freeze_dg_hybrid_basis_directory(icomm,bases,effective_ids,basis_owner,basis_fragment,ok,message)
    integer,intent(in)::icomm,effective_ids(:)
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    integer,allocatable,intent(out)::basis_owner(:),basis_fragment(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,allocatable::ownership(:)
    integer::id_rank,fragment,i,position,local_bad,global_bad,ierr
    ok=.false.;message='';local_bad=0
    call MPI_Comm_rank(icomm,id_rank,ierr);if(ierr/=MPI_SUCCESS)return
    if(size(effective_ids)<1.or.size(bases)<1.or.any(effective_ids<=0))local_bad=1
    allocate(basis_owner(size(effective_ids)),basis_fragment(size(effective_ids)),ownership(size(effective_ids)))
    basis_owner=-1;basis_fragment=0;ownership=0
    do fragment=1,size(bases)
      if(.not.allocated(bases(fragment)%global_ids))then;local_bad=1;cycle;endif
      do i=1,size(bases(fragment)%global_ids)
        position=findloc(effective_ids,int(bases(fragment)%global_ids(i)),dim=1)
        if(position<=0.or.ownership(position)/=0)then;local_bad=1;cycle;endif
        basis_owner(position)=id_rank;basis_fragment(position)=fragment;ownership(position)=1
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid local fixed-basis directory input';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,basis_owner,size(basis_owner),MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,basis_fragment,size(basis_fragment),MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,ownership,size(ownership),MPI_INTEGER,MPI_SUM,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1).or.any(basis_owner<0).or.any(basis_fragment<1))then
      message='fixed basis rows are not owned exactly once';return
    endif
    ok=.true.
#else
    ok=.false.;message='fixed-basis directory construction requires MPI'
#endif
  end subroutine freeze_dg_hybrid_basis_directory

  subroutine assemble_dg_hybrid_production_interface_rows(icomm,global_count,row_ids,traces,penalty_factor,&
      interface_rows,ok,message)
    integer,intent(in)::icomm,global_count
    integer(int64),intent(in)::row_ids(:)
    type(s_dg_hybrid_production_face_trace),intent(in)::traces(:)
    real(real64),intent(in)::penalty_factor
    complex(real64),allocatable,intent(out)::interface_rows(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    type(s_dg_hybrid_sipg_face_operator)::face
    integer::i,j,k,row,local_bad
    logical::face_ok
    character(256)::face_message
    ok=.false.;message='';local_bad=0
    if(global_count<1.or.size(traces)<1.or..not.ieee_is_finite(penalty_factor).or.penalty_factor<=0d0)&
      local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(global_count,int64)))local_bad=1
    do i=1,size(row_ids)
      if(count(row_ids==row_ids(i))/=1)local_bad=1
    enddo
    if(local_bad/=0)then;message='invalid production interface row layout';return;endif
    allocate(interface_rows(size(row_ids),global_count));interface_rows=(0d0,0d0)
    do i=1,size(traces)
      if(.not.traces(i)%frozen)cycle
      call assemble_dg_hybrid_production_face(icomm,traces(i),penalty_factor,face,face_ok,face_message)
      if(.not.face_ok)then;message=trim(face_message);return;endif
      if(any(face%global_basis_ids<1).or.any(face%global_basis_ids>global_count))then
        message='production face basis ID lies outside the fixed catalog';return
      endif
      do j=1,size(face%global_basis_ids)
        row=findloc(row_ids,int(face%global_basis_ids(j),int64),dim=1)
        if(row==0)cycle
        do k=1,size(face%global_basis_ids)
          interface_rows(row,face%global_basis_ids(k))=&
            interface_rows(row,face%global_basis_ids(k))+face%total(j,k)
        enddo
      enddo
    enddo
    ok=.true.
#else
    ok=.false.;message='production interface row assembly requires MPI'
#endif
  end subroutine assemble_dg_hybrid_production_interface_rows

  subroutine materialize_dg_hybrid_production_face_collection(icomm,origins,sizes,global_size,hgs,coef_nab,bases,&
      basis_owner,basis_fragment,effective_ids,faces,ok,message)
    integer,intent(in)::icomm,origins(:,:),sizes(:,:),global_size(3),basis_owner(:),basis_fragment(:),&
      effective_ids(:)
    real(real64),intent(in)::hgs(3),coef_nab(:,:)
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    type(s_dg_hybrid_production_face_trace),allocatable,intent(out)::faces(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    type(s_dg_hybrid_production_face_trace),allocatable::candidate(:)
    complex(real64),allocatable::value_minus(:,:),value_plus(:,:),derivative_minus(:,:),derivative_plus(:,:)
    integer::fragment,axis,tangent(2),t1,t2,position(3),neighbor_position(3),neighbor,face_count,cell_count,&
      minus_fragment,plus_fragment,minus_point(3),plus_point(3),normal_sign,periodic_shift(3),&
      i,j,g,npoint,column,id_rank,nproc,ierr,local_bad,minimum_integer
    integer,allocatable::ids_minus(:),ids_plus(:),cell_group(:),&
      cell_axis(:),cell_minus_fragment(:),cell_plus_fragment(:),cell_normal_sign(:),cell_shift(:,:),&
      cell_minus_position(:,:),cell_plus_position(:,:),group_axis(:),group_minus_fragment(:),&
      group_plus_fragment(:),group_normal_sign(:),group_shift(:,:)
    integer(int64),allocatable::minus_ids(:),plus_ids(:)
    real(real64)::normal(3),weight
    logical::face_ok,active
    ok=.false.;message='';local_bad=0
    call MPI_Comm_rank(icomm,id_rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(icomm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    if(size(origins,1)/=3.or.any(shape(origins)/=shape(sizes)).or.size(origins,2)<2.or.&
        size(bases)/=size(origins,2).or.size(basis_owner)/=size(effective_ids).or.&
        size(basis_fragment)/=size(effective_ids).or.any(basis_owner<0).or.any(basis_owner>=nproc).or.&
        any(basis_fragment<1).or.any(basis_fragment>size(bases)).or.any(global_size<=0).or.any(sizes<=0).or.any(origins<0).or.&
        any(.not.ieee_is_finite(hgs)).or.any(hgs<=0d0).or.size(coef_nab,1)<1.or.size(coef_nab,2)/=3.or.&
        any(.not.ieee_is_finite(coef_nab)))local_bad=1
    if(local_bad==0)then
      do fragment=1,size(origins,2)
        if(any(origins(:,fragment)+sizes(:,fragment)>global_size))local_bad=1
      enddo
    endif
    if(local_bad/=0)then;message='invalid production fragment topology input';return;endif
    call validate_basis_materialization(bases,effective_ids,local_bad)
    if(local_bad/=0)then;message='invalid production fragment basis materialization';return;endif
    local_bad=0
    do fragment=1,size(bases);do i=1,size(bases(fragment)%global_ids)
      minimum_integer=findloc(effective_ids,int(bases(fragment)%global_ids(i)),dim=1)
      if(minimum_integer<=0.or.basis_owner(minimum_integer)/=id_rank.or.basis_fragment(minimum_integer)/=fragment)local_bad=1
    enddo;enddo
    if(local_bad/=0)then;message='local basis disagrees with the frozen owner directory';return;endif
    if(.not.partition_is_complete(origins,sizes,global_size))then
      message='production fragment geometry is overlapping or incomplete';return
    endif

    cell_count=count_cross_fragment_face_cells(origins,sizes,global_size)
    if(cell_count<1)then;message='production topology contains no cross-fragment face';return;endif
    allocate(cell_group(cell_count),cell_axis(cell_count),cell_minus_fragment(cell_count),&
      cell_plus_fragment(cell_count),cell_normal_sign(cell_count),cell_shift(3,cell_count),&
      cell_minus_position(3,cell_count),cell_plus_position(3,cell_count),group_axis(cell_count),&
      group_minus_fragment(cell_count),group_plus_fragment(cell_count),group_normal_sign(cell_count),&
      group_shift(3,cell_count))
    cell_count=0;face_count=0
    do fragment=1,size(bases)
      do axis=1,3
        tangent=pack([1,2,3],[1,2,3]/=axis)
        do t2=0,sizes(tangent(2),fragment)-1
          do t1=0,sizes(tangent(1),fragment)-1
            position=origins(:,fragment)
            position(axis)=origins(axis,fragment)+sizes(axis,fragment)-1
            position(tangent(1))=origins(tangent(1),fragment)+t1
            position(tangent(2))=origins(tangent(2),fragment)+t2
            neighbor_position=position;neighbor_position(axis)=neighbor_position(axis)+1
            periodic_shift=0
            if(neighbor_position(axis)>=global_size(axis))periodic_shift(axis)=1
            neighbor_position=modulo(neighbor_position,global_size)
            neighbor=find_fragment_owner(neighbor_position,origins,sizes)
            if(neighbor<=0)then;message='production face has no neighboring fragment';return;endif
            if(neighbor==fragment)cycle
            cell_count=cell_count+1
            if(fragment<neighbor)then
              minus_fragment=fragment;plus_fragment=neighbor;minus_point=position;plus_point=neighbor_position
              normal_sign=1
            else
              minus_fragment=neighbor;plus_fragment=fragment;minus_point=neighbor_position;plus_point=position
              periodic_shift=-periodic_shift;normal_sign=-1
            endif
            g=0
            do j=1,face_count
              if(group_axis(j)==axis.and.group_minus_fragment(j)==minus_fragment.and.&
                  group_plus_fragment(j)==plus_fragment.and.group_normal_sign(j)==normal_sign.and.&
                  all(group_shift(:,j)==periodic_shift))then;g=j;exit;endif
            enddo
            if(g==0)then
              face_count=face_count+1;g=face_count;group_axis(g)=axis
              group_minus_fragment(g)=minus_fragment;group_plus_fragment(g)=plus_fragment
              group_normal_sign(g)=normal_sign;group_shift(:,g)=periodic_shift
            endif
            cell_group(cell_count)=g;cell_axis(cell_count)=axis
            cell_minus_fragment(cell_count)=minus_fragment;cell_plus_fragment(cell_count)=plus_fragment
            cell_normal_sign(cell_count)=normal_sign;cell_shift(:,cell_count)=periodic_shift
            cell_minus_position(:,cell_count)=minus_point;cell_plus_position(:,cell_count)=plus_point
          enddo
        enddo
      enddo
    enddo
    allocate(candidate(face_count))
    do g=1,face_count
      axis=group_axis(g);tangent=pack([1,2,3],[1,2,3]/=axis);npoint=count(cell_group==g)
      minus_fragment=group_minus_fragment(g);plus_fragment=group_plus_fragment(g)
      normal_sign=group_normal_sign(g);periodic_shift=group_shift(:,g)
      allocate(minus_ids(npoint),plus_ids(npoint),ids_minus(count(basis_fragment==minus_fragment)),&
        ids_plus(count(basis_fragment==plus_fragment)))
      ids_minus=pack(effective_ids,basis_fragment==minus_fragment)
      ids_plus=pack(effective_ids,basis_fragment==plus_fragment)
      allocate(value_minus(npoint,size(ids_minus)),derivative_minus(npoint,size(ids_minus)),&
        value_plus(npoint,size(ids_plus)),derivative_plus(npoint,size(ids_plus)))
      value_minus=(0d0,0d0);derivative_minus=(0d0,0d0)
      value_plus=(0d0,0d0);derivative_plus=(0d0,0d0);j=0;local_bad=0
      do i=1,cell_count
        if(cell_group(i)/=g)cycle
        j=j+1;minus_ids(j)=physical_grid_id(cell_minus_position(:,i),global_size)
        plus_ids(j)=physical_grid_id(cell_plus_position(:,i),global_size)
        do column=1,size(ids_minus)
          minimum_integer=findloc(bases(minus_fragment)%global_ids,int(ids_minus(column),int64),dim=1)
          if(minimum_integer<=0)cycle
          call sample_basis_value(bases(minus_fragment),minus_ids(j),minimum_integer,value_minus(j,column),face_ok)
          if(face_ok)call sample_normal_derivative(bases(minus_fragment),cell_minus_position(:,i),global_size,axis,&
            normal_sign,coef_nab(:,axis),minimum_integer,derivative_minus(j,column),face_ok)
          if(.not.face_ok)local_bad=1
        enddo
        do column=1,size(ids_plus)
          minimum_integer=findloc(bases(plus_fragment)%global_ids,int(ids_plus(column),int64),dim=1)
          if(minimum_integer<=0)cycle
          call sample_basis_value(bases(plus_fragment),plus_ids(j),minimum_integer,value_plus(j,column),face_ok)
          if(face_ok)call sample_normal_derivative(bases(plus_fragment),cell_plus_position(:,i),global_size,axis,&
            normal_sign,coef_nab(:,axis),minimum_integer,derivative_plus(j,column),face_ok)
          if(.not.face_ok)local_bad=1
        enddo
      enddo
      if(local_bad/=0)then;message='production interface lacks stencil support';return;endif
      call exchange_face_columns(value_minus,ids_minus,[ids_minus,ids_plus],effective_ids,basis_owner,icomm,110,ierr)
      if(ierr/=MPI_SUCCESS)return
      call exchange_face_columns(derivative_minus,ids_minus,[ids_minus,ids_plus],effective_ids,basis_owner,icomm,111,ierr)
      if(ierr/=MPI_SUCCESS)return
      call exchange_face_columns(value_plus,ids_plus,[ids_minus,ids_plus],effective_ids,basis_owner,icomm,112,ierr)
      if(ierr/=MPI_SUCCESS)return
      call exchange_face_columns(derivative_plus,ids_plus,[ids_minus,ids_plus],effective_ids,basis_owner,icomm,113,ierr)
      if(ierr/=MPI_SUCCESS)return
      normal=0d0;normal(axis)=real(normal_sign,real64);weight=hgs(tangent(1))*hgs(tangent(2))
      active=any([(basis_owner(findloc(effective_ids,ids_minus(i),dim=1))==id_rank,i=1,size(ids_minus))]).or.&
        any([(basis_owner(findloc(effective_ids,ids_plus(i),dim=1))==id_rank,i=1,size(ids_plus))])
      if(active)call store_local_face_trace(g,minus_fragment,plus_fragment,periodic_shift,normal,hgs(axis),&
        minus_ids,plus_ids,[(weight,i=1,npoint)],ids_minus,ids_plus,value_minus,derivative_minus,value_plus,&
        -derivative_plus,candidate(g))
      deallocate(minus_ids,plus_ids,ids_minus,ids_plus,value_minus,derivative_minus,value_plus,derivative_plus)
    enddo
    allocate(faces(face_count));faces=candidate;ok=.true.;message=''
#else
    ok=.false.;message='production face materialization requires MPI'
#endif
  end subroutine materialize_dg_hybrid_production_face_collection

#ifdef USE_MPI
  subroutine validate_basis_materialization(bases,effective_ids,bad)
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    integer,intent(in)::effective_ids(:)
    integer,intent(out)::bad
    integer::fragment,i
    bad=0
    do fragment=1,size(bases)
      if((bases(fragment)%fragment_id/=0.and.bases(fragment)%fragment_id/=fragment).or.&
          .not.allocated(bases(fragment)%global_ids).or.&
          .not.allocated(bases(fragment)%buffer_point_ids).or..not.allocated(bases(fragment)%buffer_values))then
        bad=1;cycle
      endif
      if(any(shape(bases(fragment)%buffer_values)/=&
          [size(bases(fragment)%buffer_point_ids),size(bases(fragment)%global_ids)]))then;bad=1;cycle;endif
      do i=1,size(bases(fragment)%global_ids)
        if(bases(fragment)%global_ids(i)>huge(0).or.&
            count(effective_ids==int(bases(fragment)%global_ids(i)))/=1)bad=1
      enddo
    enddo
  end subroutine validate_basis_materialization

  logical function partition_is_complete(origins,sizes,global_size) result(complete)
    integer,intent(in)::origins(:,:),sizes(:,:),global_size(3)
    integer::position(3),count_owner,i,ix,iy,iz
    complete=.true.
    do iz=0,global_size(3)-1;do iy=0,global_size(2)-1;do ix=0,global_size(1)-1
      position=[ix,iy,iz]
      count_owner=count([(all(position>=origins(:,i)).and.all(position<origins(:,i)+sizes(:,i)),i=1,size(origins,2))])
      if(count_owner/=1)then;complete=.false.;return;endif
    enddo;enddo;enddo
  end function partition_is_complete

  integer function find_fragment_owner(position,origins,sizes) result(owner)
    integer,intent(in)::position(3),origins(:,:),sizes(:,:)
    integer::i
    owner=0
    do i=1,size(origins,2)
      if(all(position>=origins(:,i)).and.all(position<origins(:,i)+sizes(:,i)))then;owner=i;return;endif
    enddo
  end function find_fragment_owner

  integer function count_cross_fragment_face_cells(origins,sizes,global_size) result(face_count)
    integer,intent(in)::origins(:,:),sizes(:,:),global_size(3)
    integer::fragment,axis,tangent(2),t1,t2,position(3),neighbor_position(3),neighbor
    face_count=0
    do fragment=1,size(origins,2);do axis=1,3
      tangent=pack([1,2,3],[1,2,3]/=axis)
      do t2=0,sizes(tangent(2),fragment)-1;do t1=0,sizes(tangent(1),fragment)-1
        position=origins(:,fragment)
        position(axis)=origins(axis,fragment)+sizes(axis,fragment)-1
        position(tangent(1))=origins(tangent(1),fragment)+t1
        position(tangent(2))=origins(tangent(2),fragment)+t2
        neighbor_position=position;neighbor_position(axis)=modulo(neighbor_position(axis)+1,global_size(axis))
        neighbor=find_fragment_owner(neighbor_position,origins,sizes)
        if(neighbor/=fragment)face_count=face_count+1
      enddo;enddo
    enddo;enddo
  end function count_cross_fragment_face_cells

  integer(int64) function physical_grid_id(position,global_size) result(identifier)
    integer,intent(in)::position(3),global_size(3)
    identifier=1_int64+int(position(1),int64)+int(global_size(1),int64)*&
      (int(position(2),int64)+int(global_size(2),int64)*int(position(3),int64))
  end function physical_grid_id

  subroutine sample_basis_value(basis,identifier,column,value,ok)
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    integer(int64),intent(in)::identifier
    integer,intent(in)::column
    complex(real64),intent(out)::value
    logical,intent(out)::ok
    integer::position
    position=findloc(basis%buffer_point_ids,identifier,dim=1);ok=position>0
    if(ok)value=basis%buffer_values(position,column)
  end subroutine sample_basis_value

  subroutine sample_normal_derivative(basis,position,global_size,axis,normal_sign,coef_nab,column,derivative,ok)
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    integer,intent(in)::position(3),global_size(3),axis,normal_sign,column
    real(real64),intent(in)::coef_nab(:)
    complex(real64),intent(out)::derivative
    logical,intent(out)::ok
    complex(real64)::plus_value,minus_value
    integer::offset,plus_position(3),minus_position(3)
    derivative=(0d0,0d0);ok=.true.
    do offset=1,size(coef_nab)
      plus_position=position;minus_position=position
      plus_position(axis)=modulo(position(axis)+offset,global_size(axis))
      minus_position(axis)=modulo(position(axis)-offset,global_size(axis))
      call sample_basis_value(basis,physical_grid_id(plus_position,global_size),column,plus_value,ok)
      if(.not.ok)return
      call sample_basis_value(basis,physical_grid_id(minus_position,global_size),column,minus_value,ok)
      if(.not.ok)return
      derivative=derivative+real(normal_sign,real64)*coef_nab(offset)*(plus_value-minus_value)
    enddo
  end subroutine sample_normal_derivative

  subroutine exchange_face_columns(values,ids,participant_ids,effective_ids,basis_owner,icomm,tag,ierr)
    complex(real64),intent(inout)::values(:,:)
    integer,intent(in)::ids(:),participant_ids(:),effective_ids(:),basis_owner(:),icomm,tag
    integer,intent(out)::ierr
    complex(real64),allocatable::send_buffer(:),receive_buffer(:)
    integer::id_rank,nproc,peer,i,j,local_count,peer_count,position

    call MPI_Comm_rank(icomm,id_rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(icomm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    if(.not.any([(basis_owner(findloc(effective_ids,participant_ids(i),dim=1))==id_rank,&
        i=1,size(participant_ids))]))then;ierr=MPI_SUCCESS;return;endif
    local_count=count([(basis_owner(findloc(effective_ids,ids(i),dim=1))==id_rank,i=1,size(ids))])
    do peer=0,nproc-1
      if(peer==id_rank)cycle
      if(.not.any([(basis_owner(findloc(effective_ids,participant_ids(i),dim=1))==peer,&
          i=1,size(participant_ids))]))cycle
      peer_count=count([(basis_owner(findloc(effective_ids,ids(i),dim=1))==peer,i=1,size(ids))])
      allocate(send_buffer(size(values,1)*local_count),receive_buffer(size(values,1)*peer_count))
      position=0
      do j=1,size(ids)
        if(basis_owner(findloc(effective_ids,ids(j),dim=1))/=id_rank)cycle
        do i=1,size(values,1);position=position+1;send_buffer(position)=values(i,j);enddo
      enddo
      call MPI_Sendrecv(send_buffer,size(send_buffer),MPI_DOUBLE_COMPLEX,peer,tag,receive_buffer,&
        size(receive_buffer),MPI_DOUBLE_COMPLEX,peer,tag,icomm,MPI_STATUS_IGNORE,ierr)
      if(ierr/=MPI_SUCCESS)return
      position=0
      do j=1,size(ids)
        if(basis_owner(findloc(effective_ids,ids(j),dim=1))/=peer)cycle
        do i=1,size(values,1);position=position+1;values(i,j)=receive_buffer(position);enddo
      enddo
      deallocate(send_buffer,receive_buffer)
    enddo
  end subroutine exchange_face_columns

  subroutine store_local_face_trace(face_id,fragment_minus,fragment_plus,periodic_shift,normal,h_normal,&
      point_ids_minus,point_ids_plus,weights,basis_ids_minus,basis_ids_plus,value_minus,outward_minus,&
      value_plus,outward_plus,face)
    integer,intent(in)::face_id,fragment_minus,fragment_plus,periodic_shift(3)
    integer(int64),intent(in)::point_ids_minus(:),point_ids_plus(:)
    integer,intent(in)::basis_ids_minus(:),basis_ids_plus(:)
    real(real64),intent(in)::normal(3),h_normal,weights(:)
    complex(real64),intent(in)::value_minus(:,:),outward_minus(:,:),value_plus(:,:),outward_plus(:,:)
    type(s_dg_hybrid_production_face_trace),intent(out)::face
    face%global_face_id=face_id;face%minus_fragment=fragment_minus;face%plus_fragment=fragment_plus
    face%periodic_shift=periodic_shift;face%canonical_normal=normal;face%h_normal=h_normal
    allocate(face%point_ids_minus(size(point_ids_minus)),face%point_ids_plus(size(point_ids_plus)),&
      face%weights(size(weights)),face%basis_ids_minus(size(basis_ids_minus)),face%basis_ids_plus(size(basis_ids_plus)),&
      face%value_minus(size(value_minus,1),size(value_minus,2)),&
      face%derivative_minus(size(outward_minus,1),size(outward_minus,2)),&
      face%value_plus(size(value_plus,1),size(value_plus,2)),&
      face%derivative_plus(size(outward_plus,1),size(outward_plus,2)))
    face%point_ids_minus=point_ids_minus;face%point_ids_plus=point_ids_plus;face%weights=weights
    face%basis_ids_minus=basis_ids_minus;face%basis_ids_plus=basis_ids_plus
    face%value_minus=value_minus;face%derivative_minus=outward_minus
    face%value_plus=value_plus;face%derivative_plus=-outward_plus
    face%fingerprint=trace_fingerprint(face_id,fragment_minus,fragment_plus,periodic_shift,normal,h_normal,&
      point_ids_minus,point_ids_plus,weights,basis_ids_minus,basis_ids_plus,value_minus,outward_minus,value_plus,&
      outward_plus)
    if(face%fingerprint==0_int64)face%fingerprint=1_int64
    face%frozen=.true.
  end subroutine store_local_face_trace

#endif

  subroutine build_dg_hybrid_production_face_trace(icomm,face_id,fragment_minus,fragment_plus,periodic_shift,&
      normal,h_normal,point_ids_minus,point_ids_plus,weights,basis_ids_minus,basis_ids_plus,value_minus,&
      outward_minus,value_plus,outward_plus,face,ok,message)
    integer,intent(in)::icomm,face_id,fragment_minus,fragment_plus,periodic_shift(3)
    integer(int64),intent(in)::point_ids_minus(:),point_ids_plus(:)
    integer,intent(in)::basis_ids_minus(:),basis_ids_plus(:)
    real(real64),intent(in)::normal(3),h_normal,weights(:)
    complex(real64),intent(in)::value_minus(:,:),outward_minus(:,:),value_plus(:,:),outward_plus(:,:)
    type(s_dg_hybrid_production_face_trace),intent(out)::face
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,j,ierr,local_bad,global_bad
    integer(int64)::local_hash,minimum_hash,maximum_hash

    ok=.false.;message=''
    local_bad=0
    if(face_id<=0.or.fragment_minus<=0.or.fragment_plus<=fragment_minus)local_bad=1
    if(size(point_ids_minus)<1.or.size(point_ids_minus)/=size(point_ids_plus))local_bad=1
    if(size(weights)/=size(point_ids_minus))local_bad=1
    if(size(value_minus,1)/=size(weights).or.size(value_minus,2)/=size(basis_ids_minus))local_bad=1
    if(size(basis_ids_minus)<1.or.size(basis_ids_plus)<1)local_bad=1
    if(any(shape(outward_minus)/=shape(value_minus)))local_bad=1
    if(size(value_plus,1)/=size(weights).or.size(value_plus,2)/=size(basis_ids_plus))local_bad=1
    if(any(shape(outward_plus)/=shape(value_plus)))local_bad=1
    if(local_bad==0)then
      if(h_normal<=0d0.or.&
          abs(sqrt(sum(normal**2))-1d0)>=1d-12.or..not.all(ieee_is_finite(normal)).or.&
          .not.ieee_is_finite(h_normal).or..not.all(ieee_is_finite(weights)).or.any(weights<=0d0))local_bad=1
      if(.not.all(ieee_is_finite(real(value_minus))).or..not.all(ieee_is_finite(aimag(value_minus))).or.&
          .not.all(ieee_is_finite(real(outward_minus))).or..not.all(ieee_is_finite(aimag(outward_minus))).or.&
          .not.all(ieee_is_finite(real(value_plus))).or..not.all(ieee_is_finite(aimag(value_plus))).or.&
          .not.all(ieee_is_finite(real(outward_plus))).or..not.all(ieee_is_finite(aimag(outward_plus))))local_bad=1
      do i=1,size(basis_ids_minus)
        if(count(basis_ids_minus==basis_ids_minus(i))/=1.or.any(basis_ids_plus==basis_ids_minus(i)))local_bad=1
      enddo
      do i=1,size(basis_ids_plus)
        if(count(basis_ids_plus==basis_ids_plus(i))/=1)local_bad=1
      enddo
      do i=1,size(point_ids_minus)
        if(point_ids_minus(i)<=0_int64.or.point_ids_plus(i)<=0_int64.or.&
            count(point_ids_minus==point_ids_minus(i))/=1.or.count(point_ids_plus==point_ids_plus(i))/=1)local_bad=1
      enddo
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid or incomplete production face trace';return;endif
    local_hash=trace_fingerprint(face_id,fragment_minus,fragment_plus,periodic_shift,normal,h_normal,&
      point_ids_minus,point_ids_plus,weights,basis_ids_minus,basis_ids_plus,value_minus,outward_minus,value_plus,outward_plus)
    call MPI_Allreduce(local_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,icomm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='production face fingerprint minimum reduction failed';return;endif
    call MPI_Allreduce(local_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='rank-disagreeing production face topology or trace';return
    endif
    face%global_face_id=face_id
    face%minus_fragment=fragment_minus;face%plus_fragment=fragment_plus
    face%periodic_shift=periodic_shift;face%canonical_normal=normal;face%h_normal=h_normal
    allocate(face%point_ids_minus(size(point_ids_minus)),face%point_ids_plus(size(point_ids_plus)),&
      face%weights(size(weights)),&
      face%basis_ids_minus(size(basis_ids_minus)),face%basis_ids_plus(size(basis_ids_plus)),&
      face%value_minus(size(value_minus,1),size(value_minus,2)),&
      face%derivative_minus(size(outward_minus,1),size(outward_minus,2)),&
      face%value_plus(size(value_plus,1),size(value_plus,2)),&
      face%derivative_plus(size(outward_plus,1),size(outward_plus,2)))
    face%point_ids_minus=point_ids_minus;face%point_ids_plus=point_ids_plus;face%weights=weights
    face%basis_ids_minus=basis_ids_minus;face%basis_ids_plus=basis_ids_plus
    face%value_minus=value_minus;face%derivative_minus=outward_minus
    face%value_plus=value_plus;face%derivative_plus=-outward_plus
    face%fingerprint=local_hash;if(face%fingerprint==0_int64)face%fingerprint=1_int64
    face%frozen=.true.;ok=.true.
#else
    ok=.false.;message='production face traces require MPI'
#endif
  end subroutine build_dg_hybrid_production_face_trace

  subroutine assemble_dg_hybrid_production_face(icomm,trace,penalty_factor,face,ok,message)
    integer,intent(in)::icomm
    type(s_dg_hybrid_production_face_trace),intent(in)::trace
    real(real64),intent(in)::penalty_factor
    type(s_dg_hybrid_sipg_face_operator),intent(out)::face
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::jump(:),average_derivative(:)
    integer::point,i,j,n,local_bad
    integer(int64)::recomputed
    ok=.false.;message=''
    call validate_stored_face(trace,local_bad,recomputed)
    if(.not.ieee_is_finite(penalty_factor).or.penalty_factor<=0d0)local_bad=1
    if(local_bad/=0)then;message='invalid mutable production face payload';return;endif
#ifdef USE_MPI
    n=size(trace%basis_ids_minus)+size(trace%basis_ids_plus)
    face%global_face_id=trace%global_face_id;face%periodic_shift=trace%periodic_shift;face%basis_count=n
    allocate(face%global_basis_ids(n),face%consistency(n,n),face%adjoint_consistency(n,n),&
      face%raw_penalty(n,n),face%physical_penalty(n,n),face%total(n,n),jump(n),average_derivative(n))
    face%global_basis_ids=[trace%basis_ids_minus,trace%basis_ids_plus]
    face%consistency=(0d0,0d0);face%adjoint_consistency=(0d0,0d0);face%raw_penalty=(0d0,0d0)
    do point=1,size(trace%weights)
      jump=[trace%value_minus(point,:),-trace%value_plus(point,:)]
      average_derivative=0.5d0*[trace%derivative_minus(point,:),trace%derivative_plus(point,:)]
      do j=1,n;do i=1,n
        face%consistency(i,j)=face%consistency(i,j)-&
          0.5d0*trace%weights(point)*conjg(jump(i))*average_derivative(j)
        face%adjoint_consistency(i,j)=face%adjoint_consistency(i,j)-&
          0.5d0*trace%weights(point)*conjg(average_derivative(i))*jump(j)
        face%raw_penalty(i,j)=face%raw_penalty(i,j)+trace%weights(point)*(penalty_factor/trace%h_normal)*&
          conjg(jump(i))*jump(j)
      enddo;enddo
    enddo
    face%physical_penalty=0.5d0*face%raw_penalty
    face%total=face%consistency+face%adjoint_consistency+face%physical_penalty
    ok=.true.
#else
    message='production grouped face assembly requires MPI'
#endif
  end subroutine assemble_dg_hybrid_production_face

  subroutine validate_dg_hybrid_production_face_collection(icomm,faces,ok,message)
    integer,intent(in)::icomm
    type(s_dg_hybrid_production_face_trace),intent(in)::faces(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,local_bad,global_bad,ierr,minimum_count,maximum_count
    integer(int64)::recomputed,collection_hash,minimum_hash,maximum_hash
    local_bad=merge(0,1,size(faces)>0)
    collection_hash=int(z'1F83D9ABFB41BD6B',int64)
    do i=1,size(faces)
      call validate_stored_face(faces(i),global_bad,recomputed);local_bad=max(local_bad,global_bad)
      collection_hash=ieor(ishftc(collection_hash,11),recomputed)
      do j=i+1,size(faces)
        if(same_physical_face(faces(i),faces(j)))local_bad=1
      enddo
    enddo
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='production face collection validation reduction failed';return;endif
    call MPI_Allreduce(size(faces),minimum_count,1,MPI_INTEGER,MPI_MIN,icomm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='production face-count minimum reduction failed';return;endif
    call MPI_Allreduce(size(faces),maximum_count,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='production face-count maximum reduction failed';return;endif
    call MPI_Allreduce(collection_hash,minimum_hash,1,MPI_INTEGER8,MPI_MIN,icomm,ierr)
    if(ierr/=MPI_SUCCESS)then;ok=.false.;message='production collection hash minimum reduction failed';return;endif
    call MPI_Allreduce(collection_hash,maximum_hash,1,MPI_INTEGER8,MPI_MAX,icomm,ierr)
    ok=ierr==MPI_SUCCESS.and.global_bad==0.and.minimum_count==maximum_count.and.minimum_hash==maximum_hash
#else
    ok=local_bad==0
#endif
    if(ok)then;message='';else;message='duplicate or invalid physical production face collection';endif
  end subroutine validate_dg_hybrid_production_face_collection

  logical function same_physical_face(first,second) result(same)
    type(s_dg_hybrid_production_face_trace),intent(in)::first,second
    same=.false.
    if(first%minus_fragment/=second%minus_fragment.or.first%plus_fragment/=second%plus_fragment)return
    if(any(first%periodic_shift/=second%periodic_shift))return
    if(maxval(abs(first%canonical_normal-second%canonical_normal))>=1d-12)return
    if(size(first%point_ids_minus)/=size(second%point_ids_minus).or.&
        size(first%point_ids_plus)/=size(second%point_ids_plus))return
    same=all(first%point_ids_minus==second%point_ids_minus).and.&
      all(first%point_ids_plus==second%point_ids_plus)
  end function same_physical_face

  subroutine validate_stored_face(face,bad,recomputed)
    type(s_dg_hybrid_production_face_trace),intent(in)::face
    integer,intent(out)::bad
    integer(int64),intent(out)::recomputed
    bad=0;recomputed=0_int64
    if(.not.face%frozen.or.face%fingerprint==0_int64.or.&
        .not.allocated(face%weights).or..not.allocated(face%basis_ids_minus).or.&
        .not.allocated(face%basis_ids_plus).or..not.allocated(face%point_ids_minus).or.&
        .not.allocated(face%point_ids_plus).or..not.allocated(face%value_minus).or.&
        .not.allocated(face%derivative_minus).or..not.allocated(face%value_plus).or.&
        .not.allocated(face%derivative_plus))then;bad=1;return;endif
    if(size(face%weights)<1.or.size(face%point_ids_minus)/=size(face%weights).or.&
        size(face%point_ids_plus)/=size(face%weights).or.&
        any(shape(face%value_minus)/=[size(face%weights),size(face%basis_ids_minus)]).or.&
        any(shape(face%derivative_minus)/=shape(face%value_minus)).or.&
        any(shape(face%value_plus)/=[size(face%weights),size(face%basis_ids_plus)]).or.&
        any(shape(face%derivative_plus)/=shape(face%value_plus)))then;bad=1;return;endif
    recomputed=trace_fingerprint(face%global_face_id,face%minus_fragment,face%plus_fragment,&
      face%periodic_shift,face%canonical_normal,face%h_normal,face%point_ids_minus,face%point_ids_plus,&
      face%weights,face%basis_ids_minus,&
      face%basis_ids_plus,face%value_minus,face%derivative_minus,face%value_plus,-face%derivative_plus)
    if(recomputed==0_int64)recomputed=1_int64
    if(recomputed/=face%fingerprint)bad=1
  end subroutine validate_stored_face

  integer(int64) function trace_fingerprint(face_id,fragment_minus,fragment_plus,shift,normal,h_normal,&
      point_ids_minus,point_ids_plus,weights,ids_minus,ids_plus,value_minus,outward_minus,value_plus,outward_plus)&
      result(hash)
    integer,intent(in)::face_id,fragment_minus,fragment_plus,shift(3),ids_minus(:),ids_plus(:)
    integer(int64),intent(in)::point_ids_minus(:),point_ids_plus(:)
    real(real64),intent(in)::normal(3),h_normal,weights(:)
    complex(real64),intent(in)::value_minus(:,:),outward_minus(:,:),value_plus(:,:),outward_plus(:,:)
    integer::i,j
    hash=int(z'510E527FADE682D1',int64)
    call mix(int(face_id,int64));call mix(int(fragment_minus,int64))
    call mix(int(fragment_plus,int64))
    do i=1,3;call mix(int(shift(i),int64));call mix(transfer(normal(i),hash));enddo
    call mix(transfer(h_normal,hash))
    do i=1,size(point_ids_minus)
      call mix(point_ids_minus(i));call mix(point_ids_plus(i));call mix(transfer(weights(i),hash))
    enddo
    do i=1,size(ids_minus);call mix(int(ids_minus(i),int64));enddo
    do i=1,size(ids_plus);call mix(int(ids_plus(i),int64));enddo
    do j=1,size(value_minus,2);do i=1,size(value_minus,1)
      call mix_complex(value_minus(i,j));call mix_complex(outward_minus(i,j))
    enddo;enddo
    do j=1,size(value_plus,2);do i=1,size(value_plus,1)
      call mix_complex(value_plus(i,j));call mix_complex(outward_plus(i,j))
    enddo;enddo
  contains
    subroutine mix(value)
      integer(int64),intent(in)::value
      hash=ieor(ishftc(hash,11),value);hash=ieor(hash,ishftc(hash,17))
    end subroutine mix
    subroutine mix_complex(value)
      complex(real64),intent(in)::value
      call mix(transfer(real(value,real64),hash));call mix(transfer(aimag(value),hash))
    end subroutine mix_complex
  end function trace_fingerprint
end module dg_hybrid_production_face_traces
