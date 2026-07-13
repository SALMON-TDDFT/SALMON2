module lcfo_wannier_sawf_templates
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter :: sawf_template_schema_version=2
  character(24),parameter :: sawf_template_magic='SALMON_SAWF_TEMPLATE_V2'

  type, public :: t_sawf_template_fingerprint
    character(256) :: geometry='', pseudopotential='', grid='', band_window=''
    character(256) :: complete_projection_shell='', symmetry='', buffer='', generator=''
    integer :: schema_version=sawf_template_schema_version
  end type

  type,public::t_sawf_ragged_local_basis
    integer::representative_fragment=0,operation_index=0
    logical::generated_independently=.false.
    complex(8),allocatable::core(:,:),buffer(:,:)
    complex(8),allocatable::gauge_unitary(:,:)
    real(8)::gauge_residual=huge(0d0)
  end type t_sawf_ragged_local_basis

  type, public :: t_sawf_template_checkpoint
    type(t_sawf_template_fingerprint) :: fingerprint
    real(8), allocatable :: centers(:,:),spreads(:)
    real(8), allocatable :: basis(:,:),buffer_basis(:,:)
    complex(8), allocatable :: orbitals(:,:),d_band(:,:,:),d_wann(:,:,:)
    complex(8),allocatable :: buffer_orbitals(:,:),band_to_wannier(:,:)
    complex(8), allocatable :: gauge_unitary(:,:)
    real(8) :: gauge_residual=huge(0d0)
  end type

  type,public::t_sawf_acceptance_checkpoint
    character(256)::supercell_fingerprint=''
    integer,allocatable::buffer_size(:)
    real(8),allocatable::center_residual(:),projector_residual(:),overlap_residual(:)
    real(8),allocatable::ww_residual(:),wp_residual(:),face_residual(:,:)
    real(8)::global_local_projector_residual=huge(0d0)
    real(8)::global_local_overlap_residual=huge(0d0)
    real(8)::global_local_fixed_h_residual=huge(0d0)
    real(8),allocatable::global_local_face_residual(:)
  end type

  public :: build_sawf_environment_orbits, validate_sawf_template_fingerprint
  public :: validate_sawf_actual_group_operation, replicate_sawf_operator_block
  public :: stitch_sawf_neighbor_gauge, validate_sawf_buffer_convergence
  public :: validate_sawf_global_local_equivalence
  public :: write_sawf_template_checkpoint, read_sawf_template_checkpoint
  public :: materialize_sawf_local_bases
  public :: apply_sawf_gauge_connection
  public :: stitch_and_apply_sawf_neighbor_pair
  public :: stitch_and_apply_sawf_neighbor_pair_real
  public :: build_sawf_neighbor_trace_overlap
  public :: whiten_sawf_buffered_overlap
  public :: stitch_sawf_buffered_neighbor_gauge
  public :: build_sawf_supercell_fingerprint, build_sawf_local_environment_fingerprint
  public :: select_sawf_environment_materialization
  public :: validate_sawf_structure_class
  public :: build_sawf_file_content_digest
  public :: measure_sawf_vacuum_occupancy
  public :: sawf_closest_periodic_cartesian
  public :: materialize_sawf_ragged_local_basis
  public :: write_sawf_materialized_basis_checkpoint
  public :: read_sawf_materialized_basis_checkpoint
  public :: stitch_sawf_materialized_neighbor_pair
  public :: validate_sawf_materialized_neighbor_closure
  public :: build_sawf_shared_buffer_point_maps
  public :: build_sawf_fragment_gauge_tree
  public :: sawf_fragments_share_face
  public :: write_sawf_acceptance_checkpoint,read_sawf_acceptance_checkpoint
  public :: validate_sawf_acceptance_checkpoint
  public :: admit_sawf_hierarchical_basis

contains
  subroutine validate_sawf_materialized_neighbor_closure(left_basis,right_basis,left_shared_points, &
      right_shared_points,grid_weight,relative_cutoff,tolerance,closure_residual, &
      alignment_residual,ok,message)
    type(t_sawf_ragged_local_basis),intent(in)::left_basis,right_basis
    integer,intent(in)::left_shared_points(:),right_shared_points(:)
    real(8),intent(in)::grid_weight,relative_cutoff,tolerance
    real(8),intent(out)::closure_residual,alignment_residual
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::cross_overlap(:,:),left_gram(:,:),right_gram(:,:),gauge(:,:)
    integer::n,i
    ok=.false.;message='';closure_residual=huge(0d0);alignment_residual=huge(0d0)
    if(.not.allocated(left_basis%buffer).or..not.allocated(right_basis%buffer).or. &
        grid_weight<=0d0.or.size(left_shared_points)<=0.or. &
        size(left_shared_points)/=size(right_shared_points).or. &
        size(left_basis%buffer,2)/=size(right_basis%buffer,2))then
      message='SAWF neighbor closure dimensions or grid weight are ambiguous';return
    end if
    n=size(left_basis%buffer,2)
    if(n<=0.or.any(left_shared_points<1).or.any(left_shared_points>size(left_basis%buffer,1)).or. &
        any(right_shared_points<1).or.any(right_shared_points>size(right_basis%buffer,1)))then
      message='SAWF neighbor closure shared-point map is invalid';return
    end if
    allocate(cross_overlap(n,n),left_gram(n,n),right_gram(n,n),gauge(n,n))
    left_gram=grid_weight*matmul(conjg(transpose(left_basis%buffer(left_shared_points,:))), &
      left_basis%buffer(left_shared_points,:))
    right_gram=grid_weight*matmul(conjg(transpose(right_basis%buffer(right_shared_points,:))), &
      right_basis%buffer(right_shared_points,:))
    cross_overlap=grid_weight*matmul(conjg(transpose(left_basis%buffer(left_shared_points,:))), &
      right_basis%buffer(right_shared_points,:))
    call stitch_sawf_buffered_neighbor_gauge(cross_overlap,left_gram,right_gram,relative_cutoff, &
      tolerance,gauge,alignment_residual,ok,message)
    if(.not.ok)return
    do i=1,n;gauge(i,i)=gauge(i,i)-(1d0,0d0);end do
    closure_residual=sqrt(sum(abs(gauge)**2)/dble(n))
    ok=closure_residual<=tolerance
    if(.not.ok)message='SAWF neighbor gauge-loop closure residual exceeds tolerance'
  end subroutine

  subroutine admit_sawf_hierarchical_basis(filename,expected_fingerprint,tolerance,checkpoint,ok,message)
    character(*),intent(in)::filename,expected_fingerprint
    real(8),intent(in)::tolerance
    type(t_sawf_acceptance_checkpoint),intent(inout)::checkpoint
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical::reusable
    call read_sawf_acceptance_checkpoint(filename,expected_fingerprint,checkpoint,reusable,ok,message)
    if(.not.ok)return
    if(.not.reusable)then
      ok=.false.;message='hierarchical SAWF acceptance provenance does not match this supercell';return
    end if
    call validate_sawf_acceptance_checkpoint(checkpoint,tolerance,ok,message)
    if(.not.ok)message='hierarchical SAWF admission rejected: '//trim(message)
  end subroutine

  subroutine validate_sawf_acceptance_checkpoint(checkpoint,tolerance,ok,message)
    type(t_sawf_acceptance_checkpoint),intent(in)::checkpoint
    real(8),intent(in)::tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::npair,nface
    ok=.false.;message=''
    if(tolerance<=0d0.or.len_trim(checkpoint%supercell_fingerprint)==0)then
      message='SAWF acceptance tolerance or fingerprint is invalid';return
    end if
    if(.not.allocated(checkpoint%buffer_size).or..not.allocated(checkpoint%center_residual).or. &
        .not.allocated(checkpoint%projector_residual).or..not.allocated(checkpoint%overlap_residual).or. &
        .not.allocated(checkpoint%ww_residual).or..not.allocated(checkpoint%wp_residual).or. &
        .not.allocated(checkpoint%face_residual).or. &
        .not.allocated(checkpoint%global_local_face_residual))then
      message='SAWF acceptance checkpoint is incomplete';return
    end if
    if(size(checkpoint%buffer_size)<3.or.any(checkpoint%buffer_size<=0).or. &
        any(checkpoint%buffer_size(2:)<=checkpoint%buffer_size(:size(checkpoint%buffer_size)-1)))then
      message='SAWF acceptance requires at least three increasing buffers';return
    end if
    npair=size(checkpoint%buffer_size)-1;nface=size(checkpoint%global_local_face_residual)
    if(nface<=0.or.size(checkpoint%face_residual,1)/=nface.or. &
        size(checkpoint%center_residual)/=npair.or. &
        size(checkpoint%projector_residual)/=npair.or.size(checkpoint%overlap_residual)/=npair.or. &
        size(checkpoint%ww_residual)/=npair.or.size(checkpoint%wp_residual)/=npair.or. &
        any(shape(checkpoint%face_residual)/=[nface,npair]))then
      message='SAWF acceptance residual series dimensions are inconsistent';return
    end if
    if(.not.all(ieee_is_finite(checkpoint%center_residual)).or. &
        .not.all(ieee_is_finite(checkpoint%projector_residual)).or. &
        .not.all(ieee_is_finite(checkpoint%overlap_residual)).or. &
        .not.all(ieee_is_finite(checkpoint%ww_residual)).or. &
        .not.all(ieee_is_finite(checkpoint%wp_residual)).or. &
        .not.all(ieee_is_finite(checkpoint%face_residual)).or. &
        .not.ieee_is_finite(checkpoint%global_local_projector_residual).or. &
        .not.ieee_is_finite(checkpoint%global_local_overlap_residual).or. &
        .not.ieee_is_finite(checkpoint%global_local_fixed_h_residual).or. &
        .not.all(ieee_is_finite(checkpoint%global_local_face_residual)))then
      message='SAWF acceptance checkpoint contains non-finite residuals';return
    end if
    if(any(checkpoint%center_residual<0d0).or.any(checkpoint%projector_residual<0d0).or. &
        any(checkpoint%overlap_residual<0d0).or.any(checkpoint%ww_residual<0d0).or. &
        any(checkpoint%wp_residual<0d0).or.any(checkpoint%face_residual<0d0).or. &
        checkpoint%global_local_projector_residual<0d0.or. &
        checkpoint%global_local_overlap_residual<0d0.or. &
        checkpoint%global_local_fixed_h_residual<0d0.or. &
        any(checkpoint%global_local_face_residual<0d0))then
      message='SAWF acceptance residuals must be nonnegative';return
    end if
    call validate_sawf_buffer_convergence(checkpoint%center_residual(npair), &
      checkpoint%projector_residual(npair),checkpoint%overlap_residual(npair), &
      checkpoint%ww_residual(npair),checkpoint%wp_residual(npair), &
      checkpoint%face_residual(:,npair),tolerance,ok,message)
    if(.not.ok)return
    call validate_sawf_global_local_equivalence(checkpoint%global_local_projector_residual, &
      checkpoint%global_local_overlap_residual,checkpoint%global_local_fixed_h_residual, &
      checkpoint%global_local_face_residual,tolerance,ok,message)
  end subroutine

  subroutine write_sawf_acceptance_checkpoint(filename,checkpoint,ok,message)
    character(*),intent(in)::filename
    type(t_sawf_acceptance_checkpoint),intent(in)::checkpoint
    logical,intent(out)::ok
    character(*),intent(out)::message
    character(32),parameter::magic='SALMON_SAWF_ACCEPTANCE_V1'
    integer::unit,ios,dims(3)
    call validate_sawf_acceptance_checkpoint(checkpoint,huge(0d0),ok,message)
    if(.not.ok)return
    dims=[size(checkpoint%buffer_size),size(checkpoint%face_residual,1), &
      size(checkpoint%global_local_face_residual)]
    open(newunit=unit,file=filename,status='replace',access='stream',form='unformatted',iostat=ios)
    if(ios/=0)then;ok=.false.;message='cannot open SAWF acceptance checkpoint';return;end if
    write(unit,iostat=ios)magic,checkpoint%supercell_fingerprint,dims,checkpoint%buffer_size, &
      checkpoint%center_residual,checkpoint%projector_residual,checkpoint%overlap_residual, &
      checkpoint%ww_residual,checkpoint%wp_residual,checkpoint%face_residual, &
      checkpoint%global_local_projector_residual,checkpoint%global_local_overlap_residual, &
      checkpoint%global_local_fixed_h_residual,checkpoint%global_local_face_residual
    close(unit);ok=ios==0
    if(.not.ok)message='cannot write SAWF acceptance checkpoint'
  end subroutine

  subroutine read_sawf_acceptance_checkpoint(filename,expected_fingerprint,checkpoint,reusable,ok,message)
    character(*),intent(in)::filename,expected_fingerprint
    type(t_sawf_acceptance_checkpoint),intent(inout)::checkpoint
    logical,intent(out)::reusable,ok
    character(*),intent(out)::message
    character(32),parameter::expected_magic='SALMON_SAWF_ACCEPTANCE_V1'
    character(32)::magic
    integer::unit,ios,dims(3),npair
    ok=.false.;reusable=.false.;message=''
    open(newunit=unit,file=filename,status='old',access='stream',form='unformatted',iostat=ios)
    if(ios/=0)then;message='cannot open SAWF acceptance checkpoint';return;end if
    read(unit,iostat=ios)magic,checkpoint%supercell_fingerprint,dims
    if(ios/=0.or.magic/=expected_magic.or.dims(1)<3.or.dims(2)<=0.or.dims(3)<=0)then
      close(unit);message='SAWF acceptance checkpoint header is invalid';return
    end if
    npair=dims(1)-1
    if(allocated(checkpoint%buffer_size))deallocate(checkpoint%buffer_size)
    if(allocated(checkpoint%center_residual))deallocate(checkpoint%center_residual)
    if(allocated(checkpoint%projector_residual))deallocate(checkpoint%projector_residual)
    if(allocated(checkpoint%overlap_residual))deallocate(checkpoint%overlap_residual)
    if(allocated(checkpoint%ww_residual))deallocate(checkpoint%ww_residual)
    if(allocated(checkpoint%wp_residual))deallocate(checkpoint%wp_residual)
    if(allocated(checkpoint%face_residual))deallocate(checkpoint%face_residual)
    if(allocated(checkpoint%global_local_face_residual))deallocate(checkpoint%global_local_face_residual)
    allocate(checkpoint%buffer_size(dims(1)),checkpoint%center_residual(npair), &
      checkpoint%projector_residual(npair),checkpoint%overlap_residual(npair), &
      checkpoint%ww_residual(npair),checkpoint%wp_residual(npair), &
      checkpoint%face_residual(dims(2),npair),checkpoint%global_local_face_residual(dims(3)))
    read(unit,iostat=ios)checkpoint%buffer_size,checkpoint%center_residual, &
      checkpoint%projector_residual,checkpoint%overlap_residual,checkpoint%ww_residual, &
      checkpoint%wp_residual,checkpoint%face_residual,checkpoint%global_local_projector_residual, &
      checkpoint%global_local_overlap_residual,checkpoint%global_local_fixed_h_residual, &
      checkpoint%global_local_face_residual
    close(unit)
    if(ios/=0)then;message='SAWF acceptance checkpoint payload is invalid';return;end if
    ok=.true.;reusable=trim(checkpoint%supercell_fingerprint)==trim(expected_fingerprint)
    if(.not.reusable)message='SAWF acceptance checkpoint supercell fingerprint mismatch'
  end subroutine

  subroutine build_sawf_fragment_gauge_tree(global_mesh,fragment_origin,fragment_shape,parent,ok,message)
    integer,intent(in)::global_mesh(3),fragment_origin(:,:),fragment_shape(:,:)
    integer,intent(out)::parent(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::fragment,candidate
    ok=.false.;message='';parent=0
    if(any(global_mesh<=0).or.size(fragment_origin,1)/=3.or. &
        any(shape(fragment_shape)/=shape(fragment_origin)).or. &
        size(parent)/=size(fragment_origin,2).or.size(parent)<=0.or. &
        any(fragment_origin<0).or.any(fragment_shape<=0))then
      message='SAWF fragment gauge tree geometry is invalid';return
    end if
    do fragment=1,size(parent)
      if(any(fragment_origin(:,fragment)+fragment_shape(:,fragment)>global_mesh))then
        message='SAWF fragment gauge tree extends outside the supercell mesh';return
      end if
    end do
    do fragment=2,size(parent)
      do candidate=1,fragment-1
        if(sawf_fragments_share_face(global_mesh,fragment_origin(:,candidate), &
            fragment_shape(:,candidate),fragment_origin(:,fragment),fragment_shape(:,fragment))) &
          parent(fragment)=candidate
      end do
      if(parent(fragment)==0)then
        write(message,'(a,i0)')'SAWF fragment gauge graph is disconnected at fragment ',fragment
        return
      end if
    end do
    ok=.true.
  end subroutine

  logical function sawf_fragments_share_face(global_mesh,left_origin,left_shape,right_origin,right_shape)
    integer,intent(in)::global_mesh(3),left_origin(3),left_shape(3),right_origin(3),right_shape(3)
    integer::axis,other,touch_count
    logical::touch,overlap
    touch_count=0
    do axis=1,3
      touch=left_origin(axis)+left_shape(axis)==right_origin(axis).or. &
        right_origin(axis)+right_shape(axis)==left_origin(axis).or. &
        (left_origin(axis)+left_shape(axis)==global_mesh(axis).and.right_origin(axis)==0).or. &
        (right_origin(axis)+right_shape(axis)==global_mesh(axis).and.left_origin(axis)==0)
      if(.not.touch)cycle
      overlap=.true.
      do other=1,3
        if(other==axis)cycle
        overlap=overlap.and.max(left_origin(other),right_origin(other))< &
          min(left_origin(other)+left_shape(other),right_origin(other)+right_shape(other))
      end do
      if(overlap)touch_count=touch_count+1
    end do
    sawf_fragments_share_face=touch_count==1
  end function

  subroutine build_sawf_shared_buffer_point_maps(global_mesh,left_origin,left_shape,right_origin, &
      right_shape,buffer_width,left_map,right_map,ok,message)
    integer,intent(in)::global_mesh(3),left_origin(3),left_shape(3),right_origin(3), &
      right_shape(3),buffer_width(3)
    integer,allocatable,intent(out)::left_map(:),right_map(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::left_id(:),right_id(:),left_order(:),right_order(:),tmp_left(:),tmp_right(:)
    integer::nl,nr,i,j,k,row,count,left_box(3),right_box(3),gid
    ok=.false.;message=''
    if(any(global_mesh<=0).or.any(left_shape<=0).or.any(right_shape<=0).or.any(buffer_width<0))then
      message='SAWF shared buffer geometry is invalid';return
    end if
    left_box=left_shape+2*buffer_width;right_box=right_shape+2*buffer_width
    nl=product(left_box);nr=product(right_box)
    allocate(left_id(nl),left_order(nl),right_id(nr),right_order(nr))
    row=0
    do k=0,left_box(3)-1;do j=0,left_box(2)-1;do i=0,left_box(1)-1
      row=row+1
      gid=modulo(left_origin(1)-buffer_width(1)+i,global_mesh(1))+ &
        global_mesh(1)*(modulo(left_origin(2)-buffer_width(2)+j,global_mesh(2))+ &
        global_mesh(2)*modulo(left_origin(3)-buffer_width(3)+k,global_mesh(3)))
      left_id(row)=gid;left_order(row)=row
    end do;end do;end do
    row=0
    do k=0,right_box(3)-1;do j=0,right_box(2)-1;do i=0,right_box(1)-1
      row=row+1
      gid=modulo(right_origin(1)-buffer_width(1)+i,global_mesh(1))+ &
        global_mesh(1)*(modulo(right_origin(2)-buffer_width(2)+j,global_mesh(2))+ &
        global_mesh(2)*modulo(right_origin(3)-buffer_width(3)+k,global_mesh(3)))
      right_id(row)=gid;right_order(row)=row
    end do;end do;end do
    call sawf_sort_grid_pairs(left_id,left_order,1,nl)
    call sawf_sort_grid_pairs(right_id,right_order,1,nr)
    if((nl>1.and.any(left_id(2:)==left_id(:nl-1))).or. &
        (nr>1.and.any(right_id(2:)==right_id(:nr-1))))then
      message='SAWF buffer contains duplicate periodic grid images';return
    end if
    allocate(tmp_left(min(nl,nr)),tmp_right(min(nl,nr)))
    i=1;j=1;count=0
    do while(i<=nl.and.j<=nr)
      if(left_id(i)<right_id(j))then;i=i+1
      else if(left_id(i)>right_id(j))then;j=j+1
      else
        count=count+1;tmp_left(count)=left_order(i);tmp_right(count)=right_order(j);i=i+1;j=j+1
      end if
    end do
    if(count==0)then;message='SAWF neighboring buffers have no shared supercell grid points';return;end if
    allocate(left_map(count),right_map(count));left_map=tmp_left(:count);right_map=tmp_right(:count)
    ok=.true.
  end subroutine

  recursive subroutine sawf_sort_grid_pairs(ids,order,lo,hi)
    integer,intent(inout)::ids(:),order(:)
    integer,intent(in)::lo,hi
    integer::i,j,pivot,tmp
    if(lo>=hi)return
    i=lo;j=hi;pivot=ids((lo+hi)/2)
    do
      do while(ids(i)<pivot);i=i+1;end do
      do while(ids(j)>pivot);j=j-1;end do
      if(i<=j)then
        tmp=ids(i);ids(i)=ids(j);ids(j)=tmp
        tmp=order(i);order(i)=order(j);order(j)=tmp
        i=i+1;j=j-1
      end if
      if(i>j)exit
    end do
    if(lo<j)call sawf_sort_grid_pairs(ids,order,lo,j)
    if(i<hi)call sawf_sort_grid_pairs(ids,order,i,hi)
  end subroutine

  subroutine stitch_sawf_materialized_neighbor_pair(left_basis,right_basis,left_shared_points, &
      right_shared_points,grid_weight, &
      relative_cutoff,alignment_tolerance,ok,message)
    type(t_sawf_ragged_local_basis),intent(in)::left_basis
    type(t_sawf_ragged_local_basis),intent(inout)::right_basis
    integer,intent(in)::left_shared_points(:),right_shared_points(:)
    real(8),intent(in)::grid_weight,relative_cutoff,alignment_tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::cross_overlap(:,:),left_gram(:,:),right_gram(:,:),gauge(:,:)
    real(8)::residual
    integer::n
    ok=.false.;message=''
    if(.not.allocated(left_basis%buffer).or..not.allocated(right_basis%buffer).or. &
        .not.allocated(right_basis%core).or..not.allocated(right_basis%gauge_unitary).or. &
        grid_weight<=0d0.or.size(left_shared_points)<=0.or. &
        size(left_shared_points)/=size(right_shared_points).or. &
        size(left_basis%buffer,2)/=size(right_basis%buffer,2).or. &
        size(right_basis%core,2)/=size(right_basis%buffer,2))then
      message='SAWF materialized neighbor basis dimensions or grid weight are ambiguous';return
    end if
    n=size(left_basis%buffer,2)
    if(n<=0.or.any(shape(right_basis%gauge_unitary)/=[n,n]).or. &
        any(left_shared_points<1).or.any(left_shared_points>size(left_basis%buffer,1)).or. &
        any(right_shared_points<1).or.any(right_shared_points>size(right_basis%buffer,1)))then
      message='SAWF materialized neighbor gauge dimensions are inconsistent';return
    end if
    allocate(cross_overlap(n,n),left_gram(n,n),right_gram(n,n),gauge(n,n))
    left_gram=grid_weight*matmul(conjg(transpose(left_basis%buffer(left_shared_points,:))), &
      left_basis%buffer(left_shared_points,:))
    right_gram=grid_weight*matmul(conjg(transpose(right_basis%buffer(right_shared_points,:))), &
      right_basis%buffer(right_shared_points,:))
    cross_overlap=grid_weight*matmul(conjg(transpose(left_basis%buffer(left_shared_points,:))), &
      right_basis%buffer(right_shared_points,:))
    call stitch_sawf_buffered_neighbor_gauge(cross_overlap,left_gram,right_gram,relative_cutoff, &
      alignment_tolerance,gauge,residual,ok,message)
    if(.not.ok)return
    right_basis%core=matmul(right_basis%core,gauge)
    right_basis%buffer=matmul(right_basis%buffer,gauge)
    right_basis%gauge_unitary=matmul(right_basis%gauge_unitary,gauge)
    right_basis%gauge_residual=max(right_basis%gauge_residual,residual)
  end subroutine

  subroutine write_sawf_materialized_basis_checkpoint(filename,supercell_fingerprint,fragment, &
      basis,ok,message)
    character(*),intent(in)::filename,supercell_fingerprint
    integer,intent(in)::fragment
    type(t_sawf_ragged_local_basis),intent(in)::basis
    logical,intent(out)::ok
    character(*),intent(out)::message
    character(32),parameter::magic='SALMON_SAWF_LOCAL_BASIS_V1'
    character(256)::stored_fingerprint
    integer::unit,ios,dims(4)
    ok=.false.;message=''
    if(fragment<=0.or.len_trim(supercell_fingerprint)==0.or..not.allocated(basis%core).or. &
        .not.allocated(basis%buffer).or..not.allocated(basis%gauge_unitary))then
      message='SAWF materialized basis checkpoint input is incomplete';return
    end if
    dims=[size(basis%core,1),size(basis%core,2),size(basis%buffer,1),size(basis%buffer,2)]
    if(dims(2)<=0.or.dims(2)/=dims(4).or.any(shape(basis%gauge_unitary)/=[dims(2),dims(2)]))then
      message='SAWF materialized basis checkpoint dimensions are inconsistent';return
    end if
    open(newunit=unit,file=filename,status='replace',access='stream',form='unformatted',iostat=ios)
    if(ios/=0)then;message='cannot open SAWF materialized basis checkpoint';return;end if
    stored_fingerprint=supercell_fingerprint
    write(unit,iostat=ios)magic,stored_fingerprint,fragment,basis%representative_fragment, &
      basis%operation_index,basis%generated_independently,dims,basis%gauge_residual, &
      basis%core,basis%buffer,basis%gauge_unitary
    close(unit)
    if(ios/=0)then;message='cannot write SAWF materialized basis checkpoint';return;end if
    ok=.true.
  end subroutine

  subroutine read_sawf_materialized_basis_checkpoint(filename,expected_supercell_fingerprint, &
      expected_fragment,basis,reusable,ok,message)
    character(*),intent(in)::filename,expected_supercell_fingerprint
    integer,intent(in)::expected_fragment
    type(t_sawf_ragged_local_basis),intent(inout)::basis
    logical,intent(out)::reusable,ok
    character(*),intent(out)::message
    character(32),parameter::expected_magic='SALMON_SAWF_LOCAL_BASIS_V1'
    character(32)::magic
    character(256)::stored_fingerprint
    integer::unit,ios,fragment,dims(4)
    ok=.false.;reusable=.false.;message=''
    open(newunit=unit,file=filename,status='old',access='stream',form='unformatted',iostat=ios)
    if(ios/=0)then;message='cannot open SAWF materialized basis checkpoint';return;end if
    read(unit,iostat=ios)magic,stored_fingerprint,fragment,basis%representative_fragment, &
      basis%operation_index,basis%generated_independently,dims,basis%gauge_residual
    if(ios/=0.or.magic/=expected_magic.or.dims(1)<=0.or.dims(2)<=0.or.dims(3)<=0.or. &
        dims(2)/=dims(4))then
      close(unit);message='SAWF materialized basis checkpoint header is invalid';return
    end if
    if(allocated(basis%core))deallocate(basis%core)
    if(allocated(basis%buffer))deallocate(basis%buffer)
    if(allocated(basis%gauge_unitary))deallocate(basis%gauge_unitary)
    allocate(basis%core(dims(1),dims(2)),basis%buffer(dims(3),dims(4)), &
      basis%gauge_unitary(dims(2),dims(2)))
    read(unit,iostat=ios)basis%core,basis%buffer,basis%gauge_unitary
    close(unit)
    if(ios/=0)then;message='SAWF materialized basis checkpoint payload is invalid';return;end if
    ok=.true.
    reusable=trim(stored_fingerprint)==trim(expected_supercell_fingerprint).and. &
      fragment==expected_fragment
    if(.not.reusable)message='SAWF materialized basis checkpoint provenance mismatch'
  end subroutine

  subroutine materialize_sawf_ragged_local_basis(representative_core,representative_buffer, &
      core_map,buffer_map,d_wann,representative_fragment,operation_index,generated_independently, &
      local_basis,ok,message)
    complex(8),intent(in)::representative_core(:,:),representative_buffer(:,:),d_wann(:,:)
    integer,intent(in)::core_map(:),buffer_map(:),representative_fragment,operation_index
    logical,intent(in)::generated_independently
    type(t_sawf_ragged_local_basis),intent(inout)::local_basis
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,allocatable::seen(:)
    integer::source,target
    if(allocated(local_basis%core))deallocate(local_basis%core)
    if(allocated(local_basis%buffer))deallocate(local_basis%buffer)
    if(allocated(local_basis%gauge_unitary))deallocate(local_basis%gauge_unitary)
    ok=.false.;message=''
    if(size(representative_core,1)<=0.or.size(representative_core,2)<=0.or. &
        size(representative_buffer,1)<=0.or.size(representative_buffer,2)/=size(representative_core,2).or. &
        size(core_map)/=size(representative_core,1).or.size(buffer_map)/=size(representative_buffer,1).or. &
        size(d_wann,1)/=size(representative_core,2).or.size(d_wann,2)/=size(representative_core,2).or. &
        representative_fragment<=0.or.operation_index<=0)then
      message='SAWF ragged local materialization dimensions or provenance are invalid';return
    end if
    allocate(local_basis%core(size(representative_core,1),size(representative_core,2)), &
      local_basis%buffer(size(representative_buffer,1),size(representative_buffer,2)))
    if(generated_independently)then
      local_basis%core=representative_core;local_basis%buffer=representative_buffer
    else
      allocate(seen(size(core_map)));seen=.false.;local_basis%core=(0d0,0d0)
      do source=1,size(core_map)
        target=core_map(source)
        if(target<1.or.target>size(core_map).or.seen(target))then
          message='SAWF ragged core grid map is not a permutation';return
        end if
        seen(target)=.true.;local_basis%core(target,:)=matmul(representative_core(source,:),d_wann)
      end do
      deallocate(seen);allocate(seen(size(buffer_map)));seen=.false.;local_basis%buffer=(0d0,0d0)
      do source=1,size(buffer_map)
        target=buffer_map(source)
        if(target<1.or.target>size(buffer_map).or.seen(target))then
          message='SAWF ragged buffer grid map is not a permutation';return
        end if
        seen(target)=.true.;local_basis%buffer(target,:)=matmul(representative_buffer(source,:),d_wann)
      end do
    end if
    local_basis%representative_fragment=representative_fragment
    local_basis%operation_index=operation_index
    local_basis%generated_independently=generated_independently
    allocate(local_basis%gauge_unitary(size(d_wann,1),size(d_wann,2)))
    local_basis%gauge_unitary=(0d0,0d0)
    do source=1,size(d_wann,1);local_basis%gauge_unitary(source,source)=(1d0,0d0);end do
    local_basis%gauge_residual=0d0
    ok=.true.
  end subroutine materialize_sawf_ragged_local_basis

  subroutine measure_sawf_vacuum_occupancy(density,threshold,fraction,ok,message)
    real(8),intent(in)::density(:),threshold;real(8),intent(out)::fraction
    logical,intent(out)::ok;character(*),intent(out)::message
    ok=.false.;message='';fraction=0d0
    if(size(density)<=0.or.threshold<=0d0.or..not.all(ieee_is_finite(density)).or. &
       .not.ieee_is_finite(threshold).or.any(density<0d0))then
      message='SAWF vacuum density input or threshold is invalid';return
    end if
    fraction=dble(count(density<threshold))/dble(size(density));ok=.true.
  end subroutine

  subroutine build_sawf_file_content_digest(filename,digest,ok,message)
    character(*),intent(in)::filename;character(*),intent(out)::digest
    logical,intent(out)::ok;character(*),intent(out)::message
    integer::unit,ios;integer(8)::h1,h2,count;character(1)::byte;character(512)::iomsg
    ok=.false.;message='';digest='';count=0
    open(newunit=unit,file=filename,access='stream',form='unformatted',status='old',action='read', &
      iostat=ios,iomsg=iomsg)
    if(ios/=0)then;message='SAWF pseudopotential content open failed: '//trim(iomsg);return;end if
    call sawf_hash_begin(h1,h2)
    do
      read(unit,iostat=ios)byte
      if(ios/=0)exit
      call sawf_hash_byte(iachar(byte),h1,h2);count=count+1
    end do
    close(unit)
    if(count<=0)then;message='SAWF pseudopotential content is empty or unreadable';return;end if
    call sawf_hash_integer8([count],h1,h2);write(digest,'(a,z16.16,z16.16)')'FILE-',h1,h2;ok=.true.
  end subroutine

  subroutine validate_sawf_structure_class(structure_class,environment_key,vacuum_fraction,orbit,ok,message)
    character(*),intent(in)::structure_class,environment_key(:)
    real(8),intent(in)::vacuum_fraction(:);integer,intent(in)::orbit(:)
    logical,intent(out)::ok;character(*),intent(out)::message
    integer::i,j,unique_count,norbit,largest_count,second_count
    integer,allocatable::orbit_count(:)
    ok=.false.;message='';unique_count=0
    if(size(environment_key)<=0.or.size(vacuum_fraction)/=size(environment_key).or. &
       size(orbit)/=size(environment_key).or.any(vacuum_fraction<0d0).or. &
       any(vacuum_fraction>1d0).or.any(orbit<=0))then
      message='SAWF structure-class measured geometry arrays are invalid';return
    end if
    do i=1,size(environment_key)
      if(len_trim(environment_key(i))==0)then;message='SAWF environment fingerprint is blank';return;end if
      do j=1,i-1;if(environment_key(j)==environment_key(i))exit;end do
      if(j==i)unique_count=unique_count+1
    end do
    norbit=maxval(orbit);allocate(orbit_count(norbit));orbit_count=0
    do i=1,size(orbit);orbit_count(orbit(i))=orbit_count(orbit(i))+1;end do
    largest_count=maxval(orbit_count);second_count=0
    do i=1,norbit
      if(orbit_count(i)<largest_count)second_count=max(second_count,orbit_count(i))
    end do
    select case(trim(structure_class))
    case('auto','amorphous')
      ok=.true.
    case('crystal')
      ok=maxval(orbit)==1
      if(.not.ok)message='SAWF crystal class has symmetry-inequivalent measured environments'
    case('defect')
      ok=.true.
    case('interface')
      ok=.true.
    case('surface')
      ok=maxval(vacuum_fraction)>1d-12
      if(.not.ok)message='SAWF surface class requires measured vacuum occupancy'
    case default
      message='SAWF structure class is invalid'
    end select
  end subroutine

  subroutine select_sawf_environment_materialization(orbit,fragment_maps,representative_fragment, &
      operation_index,generate_independently,ok,message)
    integer,intent(in)::orbit(:),fragment_maps(:,:)
    integer,intent(out)::representative_fragment(size(orbit)),operation_index(size(orbit))
    logical,intent(out)::generate_independently(size(orbit)),ok
    character(*),intent(out)::message
    integer::environment,candidate,representative,operation,norbit
    ok=.false.;message='';representative_fragment=0;operation_index=0
    generate_independently=.false.
    if(size(orbit)<=0.or.size(fragment_maps,1)/=size(orbit).or.any(orbit<=0))then
      message='SAWF materialization orbit/map dimensions are invalid';return
    end if
    norbit=maxval(orbit)
    do environment=1,size(orbit)
      representative=0
      do candidate=1,size(orbit)
        if(orbit(candidate)==orbit(environment))then;representative=candidate;exit;end if
      end do
      if(representative==0)then;message='SAWF environment orbit has no representative';return;end if
      representative_fragment(environment)=representative
      generate_independently(environment)=environment==representative
      do operation=1,size(fragment_maps,2)
        if(fragment_maps(representative,operation)==environment)then
          operation_index(environment)=operation;exit
        end if
      end do
      if(operation_index(environment)==0)then
        message='SAWF environment has no actual-group materialization operation';return
      end if
    end do
    ok=.true.
  end subroutine

  subroutine build_sawf_supercell_fingerprint(lattice,pbc,species,coordinates,grid,buffer, &
      pseudopotential_digest,band_window,projection_shell,xc,generator,key,ok,message)
    real(8),intent(in)::lattice(3,3),coordinates(:,:)
    logical,intent(in)::pbc(3);integer,intent(in)::species(:),grid(3),buffer(3)
    character(*),intent(in)::pseudopotential_digest,band_window,projection_shell,xc,generator
    character(*),intent(out)::key;logical,intent(out)::ok;character(*),intent(out)::message
    integer(8)::h1,h2;integer::i
    ok=.false.;message='';key=''
    if(size(coordinates,1)/=3.or.size(coordinates,2)/=size(species).or.size(species)<=0.or. &
       any(grid<=0).or.any(buffer<0).or..not.all(ieee_is_finite(lattice)).or. &
       .not.all(ieee_is_finite(coordinates)))then
      message='SAWF supercell fingerprint input is invalid';return
    end if
    call sawf_hash_begin(h1,h2);call sawf_hash_character('SAWF-SUPERCELL-V1',h1,h2)
    call sawf_hash_real(reshape(lattice,[9]),h1,h2);call sawf_hash_logical(pbc,h1,h2)
    call sawf_hash_integer(species,h1,h2);call sawf_hash_real(reshape(coordinates,[size(coordinates)]),h1,h2)
    call sawf_hash_integer(grid,h1,h2);call sawf_hash_integer(buffer,h1,h2)
    call sawf_hash_character(pseudopotential_digest,h1,h2);call sawf_hash_character(band_window,h1,h2)
    call sawf_hash_character(projection_shell,h1,h2);call sawf_hash_character(xc,h1,h2)
    call sawf_hash_character(generator,h1,h2)
    write(key,'(a,z16.16,z16.16)')'SC-',h1,h2;ok=.true.
  end subroutine

  subroutine build_sawf_local_environment_fingerprint(supercell_key,lattice,tolerance,species,relative_coordinates, &
      vacuum_fraction,key,ok,message)
    character(*),intent(in)::supercell_key
    real(8),intent(in)::lattice(3,3),tolerance
    integer,intent(in)::species(:);real(8),intent(in)::relative_coordinates(:,:),vacuum_fraction
    character(*),intent(out)::key;logical,intent(out)::ok;character(*),intent(out)::message
    integer(8)::h1,h2,ah1,ah2,tmp_h,qdist
    integer(8),allocatable::atom_h1(:),atom_h2(:),neighbor_q(:)
    integer,allocatable::neighbor_z(:)
    integer::i,j,k,tmp_i,n
    real(8)::delta(3),cartesian(3),scale,distance
    logical::image_ok
    ok=.false.;message='';key=''
    if(len_trim(supercell_key)==0.or.size(relative_coordinates,1)/=3.or. &
       size(relative_coordinates,2)/=size(species).or.vacuum_fraction<0d0.or.vacuum_fraction>1d0.or. &
       .not.all(ieee_is_finite(relative_coordinates)).or..not.ieee_is_finite(vacuum_fraction).or. &
       .not.all(ieee_is_finite(lattice)).or.tolerance<=0d0)then
      message='SAWF local-environment fingerprint input is invalid';return
    end if
    call sawf_hash_begin(h1,h2);call sawf_hash_character('SAWF-LOCAL-V2',h1,h2)
    call sawf_hash_character(trim(supercell_key),h1,h2)
    n=size(species);scale=max(tolerance,epsilon(1d0));allocate(atom_h1(n),atom_h2(n))
    allocate(neighbor_z(max(0,n-1)),neighbor_q(max(0,n-1)))
    do i=1,n
      k=0
      do j=1,n
        if(j==i)cycle;k=k+1;delta=relative_coordinates(:,j)-relative_coordinates(:,i)
        call sawf_closest_periodic_cartesian(lattice,delta,cartesian,image_ok)
        if(.not.image_ok)then;message='SAWF periodic closest-image search failed';return;end if
        distance=sqrt(sum(cartesian**2))
        neighbor_z(k)=species(j);neighbor_q(k)=nint(distance/scale,8)
      end do
      do j=1,max(0,n-2);do k=j+1,n-1
        if(neighbor_z(k)<neighbor_z(j).or.(neighbor_z(k)==neighbor_z(j).and.neighbor_q(k)<neighbor_q(j)))then
          tmp_i=neighbor_z(j);neighbor_z(j)=neighbor_z(k);neighbor_z(k)=tmp_i
          tmp_h=neighbor_q(j);neighbor_q(j)=neighbor_q(k);neighbor_q(k)=tmp_h
        end if
      end do;end do
      call sawf_hash_begin(ah1,ah2);call sawf_hash_integer([species(i)],ah1,ah2)
      call sawf_hash_integer(neighbor_z,ah1,ah2);call sawf_hash_integer8(neighbor_q,ah1,ah2)
      atom_h1(i)=ah1;atom_h2(i)=ah2
    end do
    do i=1,n-1;do j=i+1,n
      if(atom_h1(j)<atom_h1(i).or.(atom_h1(j)==atom_h1(i).and.atom_h2(j)<atom_h2(i)))then
        tmp_h=atom_h1(i);atom_h1(i)=atom_h1(j);atom_h1(j)=tmp_h
        tmp_h=atom_h2(i);atom_h2(i)=atom_h2(j);atom_h2(j)=tmp_h
      end if
    end do;end do
    call sawf_hash_integer8(atom_h1,h1,h2);call sawf_hash_integer8(atom_h2,h1,h2)
    qdist=nint(vacuum_fraction/scale,8);call sawf_hash_integer8([qdist],h1,h2)
    write(key,'(a,z16.16,z16.16)')'ENV-',h1,h2;ok=.true.
  end subroutine

  subroutine sawf_closest_periodic_cartesian(lattice,fractional_delta,cartesian,ok)
    real(8),intent(in)::lattice(3,3),fractional_delta(3)
    real(8),intent(out)::cartesian(3);logical,intent(out)::ok
    real(8)::inverse(3,3),det,candidate(3),candidate_cart(3),best2,bound
    integer::lo(3),hi(3),n1,n2,n3
    det=lattice(1,1)*(lattice(2,2)*lattice(3,3)-lattice(2,3)*lattice(3,2)) &
      -lattice(1,2)*(lattice(2,1)*lattice(3,3)-lattice(2,3)*lattice(3,1)) &
      +lattice(1,3)*(lattice(2,1)*lattice(3,2)-lattice(2,2)*lattice(3,1))
    ok=.false.;cartesian=0d0
    if(abs(det)<=tiny(1d0).or..not.ieee_is_finite(det))return
    inverse(1,:)=[lattice(2,2)*lattice(3,3)-lattice(2,3)*lattice(3,2), &
      lattice(1,3)*lattice(3,2)-lattice(1,2)*lattice(3,3), &
      lattice(1,2)*lattice(2,3)-lattice(1,3)*lattice(2,2)]/det
    inverse(2,:)=[lattice(2,3)*lattice(3,1)-lattice(2,1)*lattice(3,3), &
      lattice(1,1)*lattice(3,3)-lattice(1,3)*lattice(3,1), &
      lattice(1,3)*lattice(2,1)-lattice(1,1)*lattice(2,3)]/det
    inverse(3,:)=[lattice(2,1)*lattice(3,2)-lattice(2,2)*lattice(3,1), &
      lattice(1,2)*lattice(3,1)-lattice(1,1)*lattice(3,2), &
      lattice(1,1)*lattice(2,2)-lattice(1,2)*lattice(2,1)]/det
    candidate=fractional_delta-dnint(fractional_delta);cartesian=matmul(lattice,candidate)
    best2=sum(cartesian**2);bound=sqrt(sum(inverse**2))*sqrt(best2)+1d-12
    lo=floor(fractional_delta-bound)-1;hi=ceiling(fractional_delta+bound)+1
    if(any(hi-lo>200))return
    do n3=lo(3),hi(3);do n2=lo(2),hi(2);do n1=lo(1),hi(1)
      candidate=fractional_delta-real([n1,n2,n3],8);candidate_cart=matmul(lattice,candidate)
      if(sum(candidate_cart**2)<best2)then;best2=sum(candidate_cart**2);cartesian=candidate_cart;end if
    end do;end do;end do
    ok=.true.
  end subroutine

  subroutine sawf_hash_integer8(values,h1,h2)
    integer(8),intent(in)::values(:);integer(8),intent(inout)::h1,h2
    integer::i,b
    do i=1,size(values);do b=0,7
      call sawf_hash_byte(int(iand(shiftr(values(i),8*b),255_8)),h1,h2)
    end do;end do
  end subroutine

  subroutine sawf_hash_begin(h1,h2)
    integer(8),intent(out)::h1,h2;h1=1469598103_8;h2=1099511627_8
  end subroutine
  subroutine sawf_hash_byte(value,h1,h2)
    integer,intent(in)::value;integer(8),intent(inout)::h1,h2
    integer(8),parameter::prime=2147483647_8
    h1=modulo(h1*257_8+int(value,8)+1_8,prime)
    h2=modulo(h2*263_8+int(value,8)+17_8,prime)
  end subroutine
  subroutine sawf_hash_character(value,h1,h2)
    character(*),intent(in)::value;integer(8),intent(inout)::h1,h2;integer::i
    call sawf_hash_integer([len_trim(value)],h1,h2)
    do i=1,len_trim(value);call sawf_hash_byte(iachar(value(i:i)),h1,h2);end do
  end subroutine
  subroutine sawf_hash_integer(values,h1,h2)
    integer,intent(in)::values(:);integer(8),intent(inout)::h1,h2
    integer::i,b;integer(8)::word
    do i=1,size(values);word=int(values(i),8)
      do b=0,7;call sawf_hash_byte(int(iand(shiftr(word,8*b),255_8)),h1,h2);end do
    end do
  end subroutine
  subroutine sawf_hash_logical(values,h1,h2)
    logical,intent(in)::values(:);integer(8),intent(inout)::h1,h2;integer::i
    do i=1,size(values);call sawf_hash_byte(merge(1,0,values(i)),h1,h2);end do
  end subroutine
  subroutine sawf_hash_real(values,h1,h2)
    real(8),intent(in)::values(:);integer(8),intent(inout)::h1,h2
    integer::i,b;integer(8)::word
    do i=1,size(values);word=transfer(values(i),word)
      do b=0,7;call sawf_hash_byte(int(iand(shiftr(word,8*b),255_8)),h1,h2);end do
    end do
  end subroutine

  subroutine stitch_sawf_buffered_neighbor_gauge(cross_overlap,left_gram,right_gram, &
      relative_cutoff,alignment_tolerance,gauge_unitary,residual,ok,message)
    complex(8),intent(in)::cross_overlap(:,:),left_gram(:,:),right_gram(:,:)
    real(8),intent(in)::relative_cutoff,alignment_tolerance
    complex(8),intent(out)::gauge_unitary(size(cross_overlap,2),size(cross_overlap,1))
    real(8),intent(out)::residual;logical,intent(out)::ok;character(*),intent(out)::message
    complex(8),allocatable::whitened(:,:)
    allocate(whitened(size(cross_overlap,1),size(cross_overlap,2)))
    call whiten_sawf_buffered_overlap(cross_overlap,left_gram,right_gram,relative_cutoff, &
      whitened,ok,message)
    if(.not.ok)return
    call stitch_sawf_neighbor_gauge(whitened,alignment_tolerance,gauge_unitary,residual,ok,message)
  end subroutine

  subroutine whiten_sawf_buffered_overlap(cross_overlap,left_gram,right_gram,relative_cutoff, &
      whitened_overlap,ok,message)
    complex(8),intent(in)::cross_overlap(:,:),left_gram(:,:),right_gram(:,:)
    real(8),intent(in)::relative_cutoff
    complex(8),intent(out)::whitened_overlap(size(cross_overlap,1),size(cross_overlap,2))
    logical,intent(out)::ok;character(*),intent(out)::message
    complex(8),allocatable::left_inverse_sqrt(:,:),right_inverse_sqrt(:,:)
    logical::left_ok,right_ok;character(256)::detail
    ok=.false.;message='';whitened_overlap=(0d0,0d0)
    if(size(cross_overlap,1)/=size(cross_overlap,2).or. &
       any(shape(left_gram)/=shape(cross_overlap)).or.any(shape(right_gram)/=shape(cross_overlap)))then
      message='SAWF buffered overlap Gram dimensions are inconsistent';return
    end if
    call sawf_hermitian_inverse_sqrt(left_gram,relative_cutoff,left_inverse_sqrt,left_ok,detail)
    if(.not.left_ok)then;message='left '//trim(detail);return;end if
    call sawf_hermitian_inverse_sqrt(right_gram,relative_cutoff,right_inverse_sqrt,right_ok,detail)
    if(.not.right_ok)then;message='right '//trim(detail);return;end if
    whitened_overlap=matmul(left_inverse_sqrt,matmul(cross_overlap,right_inverse_sqrt))
    ok=.true.
  end subroutine

  subroutine sawf_hermitian_inverse_sqrt(gram,relative_cutoff,inverse_sqrt,ok,message)
    complex(8),intent(in)::gram(:,:);real(8),intent(in)::relative_cutoff
    complex(8),allocatable,intent(out)::inverse_sqrt(:,:)
    logical,intent(out)::ok;character(*),intent(out)::message
    complex(8),allocatable::eigenvectors(:,:),work(:),scaled(:,:)
    complex(8)::query(1);real(8),allocatable::eigenvalues(:),rwork(:)
    integer::n,info,lwork,i;external::zheev
    ok=.false.;message=''
    if(size(gram,1)<=0.or.size(gram,1)/=size(gram,2).or.relative_cutoff<=0d0)then
      message='SAWF buffered Gram dimensions or cutoff are invalid';return
    end if
    n=size(gram,1);allocate(eigenvectors(n,n),eigenvalues(n),rwork(max(1,3*n-2)))
    eigenvectors=gram;lwork=-1
    call zheev('V','U',n,eigenvectors,n,eigenvalues,query,lwork,rwork,info)
    if(info/=0)then;message='SAWF buffered Gram eigensolver query failed';return;end if
    lwork=max(1,int(real(query(1),8)));allocate(work(lwork));eigenvectors=gram
    call zheev('V','U',n,eigenvectors,n,eigenvalues,work,lwork,rwork,info)
    if(info/=0.or.maxval(eigenvalues)<=0d0)then
      message='SAWF buffered Gram eigensolver failed';return
    end if
    if(minval(eigenvalues)<=relative_cutoff*maxval(eigenvalues))then
      message='SAWF buffered Gram is rank deficient';return
    end if
    allocate(scaled(n,n),inverse_sqrt(n,n));scaled=eigenvectors
    do i=1,n;scaled(:,i)=scaled(:,i)/sqrt(eigenvalues(i));end do
    inverse_sqrt=matmul(scaled,conjg(transpose(eigenvectors)))
    ok=.true.
  end subroutine

  subroutine build_sawf_neighbor_trace_overlap(left_trace,right_trace,area_weight,overlap,ok,message)
    real(8),intent(in)::left_trace(:,:),right_trace(:,:),area_weight
    real(8),intent(out)::overlap(size(left_trace,2),size(right_trace,2))
    logical,intent(out)::ok; character(*),intent(out)::message
    ok=.false.;message='';overlap=0d0
    if(size(left_trace,1)<=0.or.size(left_trace,1)/=size(right_trace,1).or. &
       size(left_trace,2)/=size(right_trace,2).or.area_weight<=0d0)then
      message='SAWF neighbor face traces have inconsistent dimensions or weight';return
    end if
    if(.not.all(ieee_is_finite(left_trace)).or..not.all(ieee_is_finite(right_trace)).or. &
       .not.ieee_is_finite(area_weight))then
      message='SAWF neighbor face trace overlap input is non-finite';return
    end if
    overlap=area_weight*matmul(transpose(left_trace),right_trace)
    ok=.true.
  end subroutine

  subroutine stitch_and_apply_sawf_neighbor_pair_real(overlap,tolerance,neighbor_gauge, &
      basis,ww_block,wp_block,face_self_block,face_neighbor_block,gauge,residual,ok,message)
    real(8),intent(in)::overlap(:,:),neighbor_gauge(:,:),tolerance
    real(8),intent(inout)::basis(:,:),ww_block(:,:),wp_block(:,:),face_self_block(:,:), &
      face_neighbor_block(:,:)
    real(8),intent(out)::gauge(size(overlap,2),size(overlap,1)),residual
    logical,intent(out)::ok; character(*),intent(out)::message
    complex(8),allocatable::overlap_c(:,:),neighbor_c(:,:),basis_c(:,:),ww_c(:,:),wp_c(:,:), &
      self_c(:,:),face_c(:,:),gauge_c(:,:)
    real(8)::imaginary_residual
    allocate(overlap_c(size(overlap,1),size(overlap,2)),neighbor_c(size(neighbor_gauge,1), &
      size(neighbor_gauge,2)),basis_c(size(basis,1),size(basis,2)), &
      ww_c(size(ww_block,1),size(ww_block,2)),wp_c(size(wp_block,1),size(wp_block,2)), &
      self_c(size(face_self_block,1),size(face_self_block,2)), &
      face_c(size(face_neighbor_block,1),size(face_neighbor_block,2)), &
      gauge_c(size(overlap,2),size(overlap,1)))
    overlap_c=cmplx(overlap,0d0,8);neighbor_c=cmplx(neighbor_gauge,0d0,8)
    basis_c=cmplx(basis,0d0,8);ww_c=cmplx(ww_block,0d0,8);wp_c=cmplx(wp_block,0d0,8)
    self_c=cmplx(face_self_block,0d0,8);face_c=cmplx(face_neighbor_block,0d0,8)
    call stitch_and_apply_sawf_neighbor_pair(overlap_c,tolerance,neighbor_c,basis_c,ww_c,wp_c, &
      self_c,face_c,gauge_c,residual,ok,message)
    if(.not.ok)return
    imaginary_residual=maxval(abs(aimag(gauge_c)))
    if(imaginary_residual>tolerance)then
      ok=.false.; write(message,'(a,es13.5)') &
        'SAWF real GS gauge has non-negligible imaginary residual=',imaginary_residual
      return
    end if
    gauge=real(gauge_c,8);basis=real(basis_c,8);ww_block=real(ww_c,8);wp_block=real(wp_c,8)
    face_self_block=real(self_c,8);face_neighbor_block=real(face_c,8)
  end subroutine

  subroutine stitch_and_apply_sawf_neighbor_pair(overlap,tolerance,neighbor_gauge_unitary, &
      basis,ww_block,wp_block,face_self_block,face_neighbor_block,gauge_unitary,residual,ok,message)
    complex(8),intent(in)::overlap(:,:),neighbor_gauge_unitary(:,:)
    real(8),intent(in)::tolerance
    complex(8),intent(inout)::basis(:,:),ww_block(:,:),wp_block(:,:), &
      face_self_block(:,:),face_neighbor_block(:,:)
    complex(8),intent(out)::gauge_unitary(size(overlap,2),size(overlap,1))
    real(8),intent(out)::residual; logical,intent(out)::ok; character(*),intent(out)::message
    call stitch_sawf_neighbor_gauge(overlap,tolerance,gauge_unitary,residual,ok,message)
    if(.not.ok)return
    call apply_sawf_gauge_connection(gauge_unitary,neighbor_gauge_unitary,basis,ww_block, &
      wp_block,face_self_block,face_neighbor_block)
  end subroutine

  subroutine apply_sawf_gauge_connection(gauge_unitary,neighbor_gauge_unitary,basis,ww_block,wp_block, &
      face_self_block,face_neighbor_block)
    complex(8),intent(in)::gauge_unitary(:,:),neighbor_gauge_unitary(:,:)
    complex(8),intent(inout)::basis(:,:),ww_block(:,:),wp_block(:,:), &
      face_self_block(:,:),face_neighbor_block(:,:)
    complex(8),allocatable::left(:,:)
    basis=matmul(basis,gauge_unitary)
    left=conjg(transpose(gauge_unitary))
    ww_block=matmul(left,matmul(ww_block,gauge_unitary))
    wp_block=matmul(left,wp_block)
    face_self_block=matmul(left,matmul(face_self_block,gauge_unitary))
    face_neighbor_block=matmul(left,matmul(face_neighbor_block,neighbor_gauge_unitary))
  end subroutine

  subroutine materialize_sawf_local_bases(representative_basis,representative_index, &
      operation_index,independent_local,point_map,d_wann,local_basis,ok,message)
    complex(8),intent(in)::representative_basis(:,:,:),d_wann(:,:,:)
    integer,intent(in)::representative_index(:),operation_index(:),point_map(:,:)
    logical,intent(in)::independent_local(:)
    complex(8),intent(out)::local_basis(:,:,:)
    logical,intent(out)::ok; character(*),intent(out)::message
    integer::environment,source_point,target_point,representative,operation,npoint,nwann
    logical,allocatable::seen(:)
    ok=.false.; message=''; local_basis=(0d0,0d0)
    npoint=size(representative_basis,1); nwann=size(representative_basis,2)
    if(size(local_basis,1)/=npoint.or.size(local_basis,2)/=nwann.or. &
       size(local_basis,3)/=size(representative_index).or. &
       size(operation_index)/=size(representative_index).or. &
       size(independent_local)/=size(representative_index).or. &
       size(point_map,1)/=npoint.or.size(d_wann,1)/=nwann.or.size(d_wann,2)/=nwann)then
      message='SAWF local materialization dimensions are inconsistent'; return
    end if
    allocate(seen(npoint))
    do environment=1,size(representative_index)
      representative=representative_index(environment)
      if(representative<1.or.representative>size(representative_basis,3))then
        message='SAWF representative environment index is out of range'; return
      end if
      if(independent_local(environment))then
        local_basis(:,:,environment)=representative_basis(:,:,representative)
        cycle
      end if
      operation=operation_index(environment)
      if(operation<1.or.operation>size(point_map,2).or.operation>size(d_wann,3))then
        message='SAWF replicated environment operation index is invalid'; return
      end if
      seen=.false.
      do source_point=1,npoint
        target_point=point_map(source_point,operation)
        if(target_point<1.or.target_point>npoint.or.seen(target_point))then
          message='SAWF symmetry point map is not a permutation'; return
        end if
        seen(target_point)=.true.
        local_basis(target_point,:,environment)=matmul( &
          representative_basis(source_point,:,representative),d_wann(:,:,operation))
      end do
      if(.not.all(seen))then; message='SAWF symmetry point map is incomplete'; return; end if
    end do
    ok=.true.
  end subroutine

  subroutine clear_sawf_template_checkpoint(template)
    type(t_sawf_template_checkpoint),intent(inout)::template
    if(allocated(template%centers))deallocate(template%centers)
    if(allocated(template%spreads))deallocate(template%spreads)
    if(allocated(template%basis))deallocate(template%basis)
    if(allocated(template%buffer_basis))deallocate(template%buffer_basis)
    if(allocated(template%orbitals))deallocate(template%orbitals)
    if(allocated(template%buffer_orbitals))deallocate(template%buffer_orbitals)
    if(allocated(template%band_to_wannier))deallocate(template%band_to_wannier)
    if(allocated(template%d_band))deallocate(template%d_band)
    if(allocated(template%d_wann))deallocate(template%d_wann)
    if(allocated(template%gauge_unitary))deallocate(template%gauge_unitary)
    template%gauge_residual=huge(0d0)
  end subroutine

  subroutine write_sawf_template_checkpoint(filename,template,ok,message)
    character(*),intent(in)::filename
    type(t_sawf_template_checkpoint),intent(in)::template
    logical,intent(out)::ok; character(*),intent(out)::message
    integer::unit,ios,dims(19); character(512)::iomsg
    ok=.false.; message=''; dims=0
    if(.not.allocated(template%centers).or..not.allocated(template%spreads).or. &
       .not.allocated(template%basis).or..not.allocated(template%buffer_basis).or. &
       .not.allocated(template%buffer_orbitals).or..not.allocated(template%band_to_wannier).or. &
       .not.allocated(template%d_band).or..not.allocated(template%d_wann).or. &
       .not.allocated(template%gauge_unitary))then
      message='SAWF template checkpoint is missing required basis metadata'; return
    end if
    dims=[size(template%centers,1),size(template%centers,2),size(template%spreads), &
      size(template%basis,1),size(template%basis,2),size(template%buffer_basis,1), &
      size(template%buffer_basis,2),0,0, &
      size(template%d_band,1),size(template%d_band,2),size(template%d_band,3), &
      size(template%d_wann,1),size(template%d_wann,2),size(template%d_wann,3), &
      size(template%buffer_orbitals,1),size(template%buffer_orbitals,2), &
      size(template%band_to_wannier,1),size(template%band_to_wannier,2)]
    if(allocated(template%orbitals))dims(8:9)=shape(template%orbitals)
    open(newunit=unit,file=filename,access='stream',form='unformatted',status='replace', &
      action='write',iostat=ios,iomsg=iomsg)
    if(ios/=0)then; message='SAWF template checkpoint open failed: '//trim(iomsg); return; end if
    write(unit,iostat=ios,iomsg=iomsg)sawf_template_magic,template%fingerprint%schema_version, &
      template%fingerprint%geometry,template%fingerprint%pseudopotential, &
      template%fingerprint%grid,template%fingerprint%band_window, &
      template%fingerprint%complete_projection_shell,template%fingerprint%symmetry, &
      template%fingerprint%buffer,template%fingerprint%generator,dims, &
      size(template%gauge_unitary,1),size(template%gauge_unitary,2),template%gauge_residual, &
      template%centers,template%spreads,template%basis,template%buffer_basis
    if(ios==0.and.allocated(template%orbitals))write(unit,iostat=ios,iomsg=iomsg)template%orbitals
    if(ios==0)write(unit,iostat=ios,iomsg=iomsg)template%buffer_orbitals,template%band_to_wannier, &
      template%d_band,template%d_wann,template%gauge_unitary
    close(unit)
    if(ios/=0)then; message='SAWF template checkpoint write failed: '//trim(iomsg); return; end if
    ok=.true.
  end subroutine

  subroutine read_sawf_template_checkpoint(filename,expected,template,reuse,ok,message)
    character(*),intent(in)::filename
    type(t_sawf_template_fingerprint),intent(in)::expected
    type(t_sawf_template_checkpoint),intent(inout)::template
    logical,intent(out)::reuse,ok; character(*),intent(out)::message
    type(t_sawf_template_fingerprint)::stored
    integer::unit,ios,dims(19),gauge_dims(2); character(512)::iomsg
    character(24)::magic; logical::exists
    call clear_sawf_template_checkpoint(template)
    ok=.false.; reuse=.false.; message=''; inquire(file=filename,exist=exists)
    if(.not.exists)then; ok=.true.; message='SAWF cache absent; regeneration required'; return; end if
    open(newunit=unit,file=filename,access='stream',form='unformatted',status='old', &
      action='read',iostat=ios,iomsg=iomsg)
    if(ios/=0)then; message='SAWF template checkpoint open failed: '//trim(iomsg); return; end if
    read(unit,iostat=ios,iomsg=iomsg)magic,stored%schema_version,stored%geometry, &
      stored%pseudopotential,stored%grid,stored%band_window,stored%complete_projection_shell, &
      stored%symmetry,stored%buffer,stored%generator,dims,gauge_dims,template%gauge_residual
    if(ios/=0.or.magic/=sawf_template_magic)then
      close(unit); message='SAWF template checkpoint header is invalid'; return
    end if
    if(any(dims<0).or.any(gauge_dims<=0).or.dims(1)/=3.or.dims(3)/=dims(2).or. &
       any(dims(4:7)<=0).or.any(dims(16:19)<=0))then
      close(unit); message='SAWF template checkpoint dimensions are invalid'; return
    end if
    template%fingerprint=stored
    call validate_sawf_template_fingerprint(stored,expected,reuse)
    if(.not.reuse)then
      close(unit); ok=.true.; message='SAWF fingerprint mismatch; regeneration required'; return
    end if
    allocate(template%centers(dims(1),dims(2)),template%spreads(dims(3)), &
      template%basis(dims(4),dims(5)),template%buffer_basis(dims(6),dims(7)), &
      template%d_band(dims(10),dims(11),dims(12)),template%d_wann(dims(13),dims(14),dims(15)), &
      template%buffer_orbitals(dims(16),dims(17)),template%band_to_wannier(dims(18),dims(19)), &
      template%gauge_unitary(gauge_dims(1),gauge_dims(2)),stat=ios)
    if(ios/=0)then; close(unit); message='SAWF template checkpoint allocation failed'; return; end if
    if(dims(8)>0.and.dims(9)>0)allocate(template%orbitals(dims(8),dims(9)),stat=ios)
    if(ios==0)read(unit,iostat=ios,iomsg=iomsg)template%centers,template%spreads, &
      template%basis,template%buffer_basis
    if(ios==0.and.allocated(template%orbitals))read(unit,iostat=ios,iomsg=iomsg)template%orbitals
    if(ios==0)read(unit,iostat=ios,iomsg=iomsg)template%buffer_orbitals,template%band_to_wannier, &
      template%d_band,template%d_wann,template%gauge_unitary
    close(unit)
    if(ios/=0)then
      call clear_sawf_template_checkpoint(template)
      reuse=.false.; message='SAWF template checkpoint payload read failed: '//trim(iomsg); return
    end if
    ok=.true.
  end subroutine

  subroutine build_sawf_environment_orbits(equivalent,defect_intersects,orbit,regenerate,ok,message)
    logical,intent(in)::equivalent(:,:),defect_intersects(:)
    integer,intent(out)::orbit(size(defect_intersects))
    logical,intent(out)::regenerate(size(defect_intersects)),ok
    character(*),intent(out)::message
    integer::i,j,norb,head,tail,current
    integer,allocatable::queue(:)
    ok=.false.; message=''; orbit=0; regenerate=.false.
    if(size(equivalent,1)/=size(defect_intersects).or.size(equivalent,2)/=size(defect_intersects))then
      message='SAWF environment equivalence matrix has invalid dimensions'; return
    end if
    if(any(equivalent.neqv.transpose(equivalent)).or.any(.not.[(equivalent(i,i),i=1,size(orbit))]))then
      message='SAWF environment equivalence relation is not symmetric/reflexive'; return
    end if
    allocate(queue(size(orbit)));norb=0
    do i=1,size(orbit)
      if(orbit(i)==0)then
        norb=norb+1; orbit(i)=norb; regenerate(i)=.true.;head=1;tail=1;queue(1)=i
        do while(head<=tail)
          current=queue(head);head=head+1
          do j=1,size(orbit)
            if(equivalent(current,j).and.orbit(j)==0)then
              orbit(j)=norb;tail=tail+1;queue(tail)=j
            end if
          end do
        end do
      end if
    end do
    regenerate=regenerate.or.defect_intersects
    ok=.true.
  end subroutine

  subroutine validate_sawf_template_fingerprint(cached,current,reuse)
    type(t_sawf_template_fingerprint),intent(in)::cached,current
    logical,intent(out)::reuse
    reuse=cached%schema_version==current%schema_version .and. cached%geometry==current%geometry .and. &
      cached%pseudopotential==current%pseudopotential .and. cached%grid==current%grid .and. &
      cached%band_window==current%band_window .and. &
      cached%complete_projection_shell==current%complete_projection_shell .and. &
      cached%symmetry==current%symmetry .and. cached%buffer==current%buffer .and. &
      cached%generator==current%generator
  end subroutine

  subroutine validate_sawf_actual_group_operation(in_actual_group,is_parent_only,ok,message)
    logical,intent(in)::in_actual_group,is_parent_only
    logical,intent(out)::ok; character(*),intent(out)::message
    ok=in_actual_group.and..not.is_parent_only; message=''
    if(.not.ok)message='parent-crystal operation is not in the actual defective-supercell group'
  end subroutine

  subroutine replicate_sawf_operator_block(block,d_left,d_right,replica)
    complex(8),intent(in)::block(:,:),d_left(:,:),d_right(:,:)
    complex(8),intent(out)::replica(size(d_left,2),size(d_right,2))
    replica=matmul(conjg(transpose(d_left)),matmul(block,d_right))
  end subroutine

  subroutine stitch_sawf_neighbor_gauge(overlap,tolerance,gauge_unitary,residual,ok,message)
    complex(8),intent(in)::overlap(:,:); real(8),intent(in)::tolerance
    complex(8),intent(out)::gauge_unitary(size(overlap,2),size(overlap,1))
    real(8),intent(out)::residual; logical,intent(out)::ok; character(*),intent(out)::message
    complex(8),allocatable::a(:,:),u(:,:),vt(:,:),work(:),aligned(:,:)
    real(8),allocatable::singular(:),rwork(:); integer::n,info,lwork,i
    ok=.false.; message=''; residual=huge(0d0); gauge_unitary=(0d0,0d0)
    if(size(overlap,1)/=size(overlap,2).or.tolerance<=0d0)then
      message='SAWF neighbor overlap has ambiguous orbital counts'; return
    end if
    n=size(overlap,1); lwork=max(1,8*n)
    allocate(a(n,n),u(n,n),vt(n,n),work(lwork),singular(n),rwork(max(1,5*n)),aligned(n,n))
    a=overlap
    call zgesvd('A','A',n,n,a,n,singular,u,n,vt,n,work,lwork,rwork,info)
    if(info/=0)then; message='SAWF neighbor overlap SVD failed'; return; end if
    if(minval(singular)<=tolerance)then; message='SAWF neighbor overlap is rank deficient'; return; end if
    gauge_unitary=matmul(conjg(transpose(vt)),conjg(transpose(u)))
    aligned=matmul(overlap,gauge_unitary)
    do i=1,n; aligned(i,i)=aligned(i,i)-(1d0,0d0); end do
    residual=sqrt(sum(abs(aligned)**2))/sqrt(dble(n))
    ok=residual<=tolerance
    if(.not.ok)message='SAWF post-alignment gauge residual exceeds tolerance'
  end subroutine

  subroutine validate_sawf_buffer_convergence(center_residual,projector_residual,overlap_residual, &
      ww_residual,wp_residual,face_residual,tolerance,ok,message)
    real(8),intent(in)::center_residual,projector_residual,overlap_residual,ww_residual,wp_residual
    real(8),intent(in)::face_residual(:),tolerance; logical,intent(out)::ok; character(*),intent(out)::message
    ok=max(center_residual,projector_residual,overlap_residual,ww_residual,wp_residual, &
      maxval(face_residual))<=tolerance; message=''
    if(.not.ok)message='SAWF buffer convergence failed for orbital or operator/DG face data'
  end subroutine

  subroutine validate_sawf_global_local_equivalence(occupied_projector_residual,overlap_residual, &
      fixed_hamiltonian_residual,face_residual,tolerance,ok,message)
    real(8),intent(in)::occupied_projector_residual,overlap_residual,fixed_hamiltonian_residual
    real(8),intent(in)::face_residual(:),tolerance; logical,intent(out)::ok; character(*),intent(out)::message
    ok=max(occupied_projector_residual,overlap_residual,fixed_hamiltonian_residual, &
      maxval(face_residual))<=tolerance; message=''
    if(.not.ok)message='hierarchical SAWF is not operator-equivalent to monolithic global SAWF'
  end subroutine
end module lcfo_wannier_sawf_templates
