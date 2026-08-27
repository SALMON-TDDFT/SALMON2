module dg_hybrid_sparse_operators
  use,intrinsic::iso_fortran_env,only:int64,real64
  implicit none
  private
  type,public::s_dg_hybrid_sparse_operators
    logical::valid=.false.
    integer::global_count=0
    integer(int64)::selection_fingerprint=0_int64,window_fingerprint=0_int64,&
      packet_fingerprint=0_int64,complement_fingerprint=0_int64,metric_fingerprint=0_int64,&
      position_convention_fingerprint=0_int64,&
      fingerprint=0_int64,persistent_bytes=0_int64,transient_peak_bytes=0_int64
    integer(int64),allocatable::owned_row_ids(:)
    integer,allocatable::row_offsets(:),column_ids(:)
    complex(real64),allocatable::metric_values(:),hamiltonian_values(:),position_values(:,:)
  end type s_dg_hybrid_sparse_operators

  type,public::s_dg_hybrid_coupling_envelope
    logical::valid=.false.
    integer::global_count=0
    integer(int64)::structure_fingerprint=0_int64
    integer(int64)::schedule_fingerprint=0_int64
    integer,allocatable::row_offsets(:),column_ids(:)
    complex(real64),allocatable::component_values(:,:)
  contains
    procedure::find_edge=>find_envelope_edge
  end type s_dg_hybrid_coupling_envelope

  public::build_dg_hybrid_operator_envelope,set_dg_hybrid_operator_edge
contains
  subroutine build_dg_hybrid_operator_envelope(global_count,basis_edges,nonlocal_edges,face_edges,&
      component_count,envelope,ok,message)
    integer,intent(in)::global_count,basis_edges(:,:),nonlocal_edges(:,:),face_edges(:,:),component_count
    type(s_dg_hybrid_coupling_envelope),intent(out)::envelope
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,allocatable::present(:,:)
    integer::row,column,i,j,k,nnz
    ok=.false.;message=''
    if(global_count<1.or.component_count<1.or.size(basis_edges,1)/=2.or.&
        size(nonlocal_edges,1)/=2.or.size(face_edges,1)/=2)then
      message='invalid operator coupling-envelope extent';return
    endif
    allocate(present(global_count,global_count));present=.false.
    do i=1,global_count;present(i,i)=.true.;enddo
    call add_edges(basis_edges);call add_edges(nonlocal_edges);call add_edges(face_edges)
    if(.not.ok.and.len_trim(message)>0)return
    nnz=count(present)
    allocate(envelope%row_offsets(global_count+1),envelope%column_ids(nnz),&
      envelope%component_values(component_count,nnz))
    envelope%row_offsets(1)=1;k=0
    do row=1,global_count
      do column=1,global_count
        if(.not.present(row,column))cycle
        k=k+1;envelope%column_ids(k)=column
      enddo
      envelope%row_offsets(row+1)=k+1
    enddo
    envelope%component_values=(0d0,0d0)
    envelope%structure_fingerprint=int(z'082EFA98EC4E6C89',int64)
    do row=1,global_count
      envelope%structure_fingerprint=ieor(ishftc(envelope%structure_fingerprint,7),int(row,int64))
      do j=envelope%row_offsets(row),envelope%row_offsets(row+1)-1
        envelope%structure_fingerprint=ieor(ishftc(envelope%structure_fingerprint,7),&
          int(envelope%column_ids(j),int64))
      enddo
    enddo
    if(envelope%structure_fingerprint==0_int64)envelope%structure_fingerprint=1_int64
    envelope%schedule_fingerprint=ieor(ishftc(envelope%structure_fingerprint,19),&
      int(component_count,int64))
    if(envelope%schedule_fingerprint==0_int64)envelope%schedule_fingerprint=1_int64
    envelope%global_count=global_count;envelope%valid=.true.;ok=.true.
  contains
    subroutine add_edges(edges)
      integer,intent(in)::edges(:,:)
      integer::edge,a,b
      do edge=1,size(edges,2)
        a=edges(1,edge);b=edges(2,edge)
        if(a<1.or.a>global_count.or.b<1.or.b>global_count)then
          message='operator coupling edge is outside the basis extent';ok=.false.;return
        endif
        present(a,b)=.true.;present(b,a)=.true.
      enddo
      ok=.true.
    end subroutine add_edges
  end subroutine build_dg_hybrid_operator_envelope

  subroutine set_dg_hybrid_operator_edge(envelope,component,row,column,value,ok,message)
    type(s_dg_hybrid_coupling_envelope),intent(inout)::envelope
    integer,intent(in)::component,row,column
    complex(real64),intent(in)::value
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::position
    ok=.false.;message=''
    if(.not.envelope%valid.or.component<1.or.component>size(envelope%component_values,1))then
      message='invalid operator component';return
    endif
    position=envelope%find_edge(row,column)
    if(position==0)then;message='operator value is outside the coupling envelope';return;endif
    envelope%component_values(component,position)=value;ok=.true.
  end subroutine set_dg_hybrid_operator_edge

  integer function find_envelope_edge(envelope,row,column) result(position)
    class(s_dg_hybrid_coupling_envelope),intent(in)::envelope
    integer,intent(in)::row,column
    integer::edge
    position=0
    if(.not.envelope%valid.or.row<1.or.row>envelope%global_count)return
    do edge=envelope%row_offsets(row),envelope%row_offsets(row+1)-1
      if(envelope%column_ids(edge)==column)then;position=edge;return;endif
    enddo
  end function find_envelope_edge
end module dg_hybrid_sparse_operators
