#include "config.h"
module lcfo_wannier_sawf
  ! Symmetry discovery/representation remains here; scalable environment
  ! ownership, template provenance and operator gates live in
  ! lcfo_wannier_sawf_templates so site_symmetry is never the sole gate.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_null_char
  use, intrinsic :: iso_fortran_env, only: int64
  use sym_sub, only: read_symmetry_file
  use dg_overlapping_wannier_projection, only: &
    t_sawf_projection_channel=>t_dg_projection_channel,&
    dg_complete_shell_target_count,build_dg_complete_shell_manifest,&
    validate_dg_complete_shell_manifest,&
    shared_real_harmonic_value=>dg_real_harmonic_value
  implicit none
  private

  integer(int64), parameter :: max_periodic_candidates = 1000000_int64
  integer(c_int), parameter :: sawf_spglib_success = 0_c_int
  integer(c_int), parameter :: sawf_spglib_count_failure = 2_c_int
  integer(c_int), parameter :: sawf_spglib_capacity_too_small = 3_c_int
  integer(c_int), parameter :: sawf_spglib_api_failure = 4_c_int
  integer(c_int), parameter :: sawf_spglib_allocation_failure = 5_c_int
  ! C64 has at most 48*32=1536 point/primitive-translation operations.
  integer, parameter :: max_sawf_symmetry_operations = 4096

  type, public :: t_sawf_symop
    integer :: W(3,3) = 0
    real(8) :: tau(3) = 0.0d0
    real(8) :: R(3,3) = 0.0d0
    integer, allocatable :: atom_map(:)
    real(8) :: integer_residual = 0.0d0
    real(8) :: metric_residual = 0.0d0
    real(8) :: orthogonality_residual = 0.0d0
    real(8) :: determinant_residual = 0.0d0
    real(8) :: atom_map_residual = 0.0d0
  end type t_sawf_symop

  type, public :: t_sawf_crystallographic_catalog
    integer :: space_group_number = 0
    integer :: hall_number = 0
    character(6) :: point_group_symbol = ''
    integer, allocatable :: integer_rotation(:,:,:)
    real(8), allocatable :: fractional_translation(:,:)
    type(t_sawf_symop), allocatable :: operations(:)
  end type t_sawf_crystallographic_catalog

  type :: t_sawf_operation_index
    integer, allocatable :: slots(:)
  end type t_sawf_operation_index

  public :: load_sawf_symmetry_auto, load_sawf_symmetry_file
  public :: load_sawf_crystallographic_catalog_auto
  public :: build_sawf_spd_projection_map, validate_sawf_spd_projection_map
  public :: t_sawf_projection_channel
  public :: build_sawf_wannier_representation, sawf_real_harmonic_value
  public :: sawf_spd_projection_count, write_sawf_projection_block
  public :: sawf_projection_shell_lmax
  public :: build_sawf_operation_product_table

#ifdef HAVE_SPGLIB
  interface
    integer(c_int) function salmon_sawf_spglib_get_symmetry(lattice_columns, &
        fractional_positions, species, num_atoms, tolerance, output_capacity, &
        rotations, translations, num_operations) &
        bind(C, name='salmon_sawf_spglib_get_symmetry')
      import :: c_double, c_int
      real(c_double), intent(in) :: lattice_columns(*)
      real(c_double), intent(in) :: fractional_positions(*)
      integer(c_int), intent(in) :: species(*)
      integer(c_int), value :: num_atoms
      real(c_double), value :: tolerance
      integer(c_int), value :: output_capacity
      integer(c_int), intent(out) :: rotations(*)
      real(c_double), intent(out) :: translations(*)
      integer(c_int), intent(out) :: num_operations
    end function salmon_sawf_spglib_get_symmetry
    integer(c_int) function salmon_sawf_spglib_get_dataset_metadata( &
        lattice_columns, fractional_positions, species, num_atoms, tolerance, &
        spacegroup_number, hall_number, pointgroup_symbol, pointgroup_capacity) &
        bind(C, name='salmon_sawf_spglib_get_dataset_metadata')
      import :: c_char, c_double, c_int
      real(c_double), intent(in) :: lattice_columns(*)
      real(c_double), intent(in) :: fractional_positions(*)
      integer(c_int), intent(in) :: species(*)
      integer(c_int), value :: num_atoms
      real(c_double), value :: tolerance
      integer(c_int), intent(out) :: spacegroup_number, hall_number
      character(c_char), intent(out) :: pointgroup_symbol(*)
      integer(c_int), value :: pointgroup_capacity
    end function salmon_sawf_spglib_get_dataset_metadata
  end interface
#endif

contains
  subroutine load_sawf_crystallographic_catalog_auto(lattice,frac_pos,species,tolerance,catalog,ok,message)
    real(8),intent(in)::lattice(3,3),frac_pos(:,:),tolerance
    integer,intent(in)::species(:)
    type(t_sawf_crystallographic_catalog),intent(out)::catalog
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef HAVE_SPGLIB
    integer(c_int),allocatable::c_species(:)
    integer(c_int)::status,spacegroup_number,hall_number
    character(c_char)::pointgroup_symbol(7)
    integer::i,nop
#endif
    catalog%space_group_number=0; catalog%hall_number=0; catalog%point_group_symbol=''
    call load_sawf_symmetry_auto(lattice,frac_pos,species,tolerance,catalog%operations,ok,message)
    if(.not.ok)return
#ifndef HAVE_SPGLIB
    message='SAWF automatic crystallographic catalog requires SALMON built with USE_SPGLIB=ON'
    ok=.false.; return
#else
    allocate(c_species(size(species))); c_species=int(species,c_int)
    pointgroup_symbol=c_null_char
    status=salmon_sawf_spglib_get_dataset_metadata(lattice,frac_pos,c_species, &
      int(size(species),c_int),real(tolerance,c_double),spacegroup_number,hall_number, &
      pointgroup_symbol,int(size(pointgroup_symbol),c_int))
    if(status/=sawf_spglib_success)then
      call set_spglib_error(status,'dataset metadata',message); ok=.false.; return
    end if
    catalog%space_group_number=int(spacegroup_number)
    catalog%hall_number=int(hall_number)
    do i=1,min(len(catalog%point_group_symbol),size(pointgroup_symbol))
      if(pointgroup_symbol(i)==c_null_char)exit
      catalog%point_group_symbol(i:i)=pointgroup_symbol(i)
    end do
    nop=size(catalog%operations)
    allocate(catalog%integer_rotation(3,3,nop),catalog%fractional_translation(3,nop))
    do i=1,nop
      catalog%integer_rotation(:,:,i)=catalog%operations(i)%W
      catalog%fractional_translation(:,i)=catalog%operations(i)%tau
    end do
    ok=.true.; message=''
#endif
  end subroutine load_sawf_crystallographic_catalog_auto

  subroutine build_sawf_operation_product_table(operations,lattice,lattice_inv,tolerance, &
      left_index,right_index,product_index,ok,message)
    type(t_sawf_symop),intent(in)::operations(:)
    real(8),intent(in)::lattice(3,3),lattice_inv(3,3),tolerance
    integer,allocatable,intent(out)::left_index(:),right_index(:),product_index(:)
    logical,intent(out)::ok; character(*),intent(out)::message
    type(t_sawf_operation_index)::operation_index
    integer::iop,jop,kop,n,mapped_atom,relation
    integer(int64)::product_w(3,3)
    real(8)::product_tau(3),residual
    logical::product_ok,residual_ok
    ok=.false.; message=''; n=size(operations)
    call build_operation_index(operations,operation_index,ok,message); if(.not.ok)return
    allocate(left_index(n*n),right_index(n*n),product_index(n*n)); relation=0
    do iop=1,n; do jop=1,n
      relation=relation+1
      call multiply_integer_matrices(int(operations(iop)%W,int64), &
        int(operations(jop)%W,int64),product_w,product_ok)
      if(.not.product_ok)then; message='SAWF product-table integer overflow'; ok=.false.; return; end if
      product_tau=operations(iop)%tau+matmul(real(operations(iop)%W,8),operations(jop)%tau)
      mapped_atom=operations(iop)%atom_map(operations(jop)%atom_map(1))
      call find_indexed_operation(operations,operation_index,product_w,mapped_atom,kop)
      if(kop==0)then; message='SAWF product-table operation is absent'; ok=.false.; return; end if
      call translation_residual(product_tau-operations(kop)%tau,lattice,lattice_inv,residual,residual_ok)
      if(.not.residual_ok.or.residual>tolerance)then
        message='SAWF product-table translation residual exceeds tolerance'; ok=.false.; return
      end if
      left_index(relation)=iop; right_index(relation)=jop; product_index(relation)=kop
    end do; end do
    ok=.true.
  end subroutine


  subroutine sawf_projection_shell_lmax(num_atoms, num_wann, max_l, ok, message)
    integer, intent(in) :: num_atoms, num_wann
    integer, intent(out) :: max_l
    logical, intent(out) :: ok
    character(*), intent(out) :: message

    ok = .false.; message = ''; max_l = -1
    if (num_atoms <= 0) then
      message = 'SAWF shell selection requires a positive atom count'
    else if (num_wann == 4*num_atoms) then
      max_l = 1; ok = .true.
    else if (num_wann == 9*num_atoms) then
      max_l = 2; ok = .true.
    else
      write(message,'(a,i0,a,i0,a,i0)') 'SAWF requires complete s+p+d shells or complete s+p shells: num_wann=', &
        num_wann, ' sp_required=', 4*num_atoms, ' spd_required=', 9*num_atoms
    end if
  end subroutine sawf_projection_shell_lmax

  subroutine sawf_spd_projection_count(num_atoms, channel_count, ok, message, max_l)
    integer, intent(in) :: num_atoms
    integer, intent(out) :: channel_count
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, intent(in), optional :: max_l
    integer :: shell_l

    if (num_atoms < 0) then
      ok=.false.;channel_count=0
      message = 'SAWF projection map atom count is negative'
      return
    end if
    if(num_atoms==0)then
      channel_count=0;ok=.true.;message='';return
    endif
    shell_l=2;if(present(max_l))shell_l=max_l
    call dg_complete_shell_target_count(num_atoms,shell_l,channel_count,ok,message)
  end subroutine sawf_spd_projection_count


  subroutine build_sawf_spd_projection_map(num_atoms, channels, ok, message, max_l)
    integer, intent(in) :: num_atoms
    type(t_sawf_projection_channel), allocatable, intent(out) :: channels(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, intent(in), optional :: max_l
    integer :: shell_l,atom
    integer,allocatable::atom_ids(:)

    shell_l = 2
    if (present(max_l)) shell_l = max_l
    if(num_atoms==0)then
      allocate(channels(0));ok=.true.;message='';return
    endif
    allocate(atom_ids(num_atoms));atom_ids=[(atom,atom=1,num_atoms)]
    call build_dg_complete_shell_manifest(atom_ids,shell_l,channels,ok,message)
  end subroutine build_sawf_spd_projection_map


  subroutine validate_sawf_spd_projection_map(channels, num_atoms, ok, message)
    type(t_sawf_projection_channel), intent(in) :: channels(:)
    integer, intent(in) :: num_atoms
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: max_l,atom
    integer,allocatable::atom_ids(:)

    ok = .false.
    message = ''
    call sawf_projection_shell_lmax(num_atoms, size(channels), max_l, ok, message)
    if (.not. ok) return
    allocate(atom_ids(num_atoms));atom_ids=[(atom,atom=1,num_atoms)]
    call validate_dg_complete_shell_manifest(channels,atom_ids,max_l,ok,message)
    if(.not.ok)message='SAWF pseudo_channels requires complete s+p+d shells or complete s+p shells: '//trim(message)
  end subroutine validate_sawf_spd_projection_map


  subroutine write_sawf_projection_block(iunit, site_symmetry, fractional_positions, ok, message, max_l)
    integer, intent(in) :: iunit
    character(*), intent(in) :: site_symmetry
    real(8), intent(in) :: fractional_positions(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, intent(in), optional :: max_l
    integer :: atom, io_status, shell_l
    character(256) :: io_message

    ok = .false.
    message = ''
    shell_l = 2
    if (present(max_l)) shell_l = max_l
    if (trim(site_symmetry) == 'off') then
      write(iunit,'(a)',iostat=io_status,iomsg=io_message) 'random'
      if (io_status /= 0) then
        message = 'SAWF projection block write failed: '//trim(io_message)
        return
      end if
      ok = .true.
      return
    end if
    if (trim(site_symmetry) /= 'auto' .and. trim(site_symmetry) /= 'file') then
      message = 'SAWF projection block mode must be off, auto, or file'
      return
    end if
    if (size(fractional_positions,1) /= 3 .or. size(fractional_positions,2) <= 0) then
      message = 'SAWF projection block requires 3-by-N fractional atom positions'
      return
    end if
    if (.not. all(ieee_is_finite(fractional_positions))) then
      message = 'SAWF projection block contains non-finite fractional atom positions'
      return
    end if
    do atom = 1, size(fractional_positions,2)
      if (shell_l == 1) then
        write(iunit,'(a,f18.12,a,f18.12,a,f18.12,a)',iostat=io_status,iomsg=io_message) &
          'f=', fractional_positions(1,atom), ',', fractional_positions(2,atom), ',', &
          fractional_positions(3,atom), ':s;p'
      else
        write(iunit,'(a,f18.12,a,f18.12,a,f18.12,a)',iostat=io_status,iomsg=io_message) &
          'f=', fractional_positions(1,atom), ',', fractional_positions(2,atom), ',', &
          fractional_positions(3,atom), ':s;p;d'
      end if
      if (io_status /= 0) then
        message = 'SAWF projection block write failed: '//trim(io_message)
        return
      end if
    end do
    ok = .true.
  end subroutine write_sawf_projection_block


  real(8) function sawf_real_harmonic_value(l, m, r) result(value)
    integer, intent(in) :: l, m
    real(8), intent(in) :: r(3)
    value=shared_real_harmonic_value(l,m,r)
  end function sawf_real_harmonic_value


  subroutine build_sawf_wannier_representation(operations, channels, d_wann, ok, message)
    type(t_sawf_symop), intent(in) :: operations(:)
    type(t_sawf_projection_channel), intent(in) :: channels(:)
    real(8), allocatable, intent(out) :: d_wann(:,:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8) :: p_rep(3,3), d_rep(5,5), orthogonality
    integer :: allocation_status, num_atoms, iop, source_atom, target_atom, row0, col0, channels_per_atom
    character(256) :: allocation_message

    ok = .false.
    message = ''
    if (size(channels) == 0) then
      message = 'SAWF pseudo_channels requires complete s+p+d shells'
      return
    end if
    num_atoms = maxval(channels%atom)
    call validate_sawf_spd_projection_map(channels, num_atoms, ok, message)
    if (.not. ok) return
    channels_per_atom = size(channels)/num_atoms
    if (size(operations) == 0) then
      message = 'SAWF D_wann requires at least one symmetry operation'
      ok = .false.
      return
    end if

    allocate(d_wann(size(channels),size(channels),size(operations)), &
      stat=allocation_status, errmsg=allocation_message)
    if (allocation_status /= 0) then
      message = 'SAWF D_wann allocation failed: '//trim(allocation_message)
      ok = .false.
      return
    end if
    d_wann = 0.0d0
    do iop = 1, size(operations)
      if (.not. allocated(operations(iop)%atom_map) .or. &
          size(operations(iop)%atom_map) /= num_atoms) then
        write(message,'(a,i0,a)') 'SAWF operation ', iop, ' has an invalid atom map for D_wann'
        ok = .false.
        return
      end if
      if (any(operations(iop)%atom_map < 1) .or. &
          any(operations(iop)%atom_map > num_atoms)) then
        write(message,'(a,i0,a)') 'SAWF operation ', iop, ' atom map is out of range'
        ok = .false.
        return
      end if
      if (.not. all(ieee_is_finite(operations(iop)%R))) then
        write(message,'(a,i0,a)') 'SAWF operation ', iop, &
          ' has a non-finite Cartesian rotation'
        ok = .false.
        return
      end if
      orthogonality = maxval(abs(matmul(transpose(operations(iop)%R), &
        operations(iop)%R)-identity_real()))
      if (.not. ieee_is_finite(orthogonality)) then
        write(message,'(a,i0,a)') 'SAWF operation ', iop, &
          ' produced a non-finite orthogonality residual'
        ok = .false.
        return
      end if
      if (orthogonality > 1.0d-10) then
        write(message,'(a,i0,a,es12.4)') 'SAWF operation ', iop, &
          ' is not orthogonal while building D_wann; residual=', orthogonality
        ok = .false.
        return
      end if
      call real_p_representation(operations(iop)%R, p_rep)
      call real_d_representation(operations(iop)%R, d_rep)
      do source_atom = 1, num_atoms
        target_atom = operations(iop)%atom_map(source_atom)
        col0 = channels_per_atom*(source_atom-1)
        row0 = channels_per_atom*(target_atom-1)
        d_wann(row0+1,col0+1,iop) = 1.0d0
        d_wann(row0+2:row0+4,col0+2:col0+4,iop) = p_rep
        if (channels_per_atom == 9) d_wann(row0+5:row0+9,col0+5:col0+9,iop) = d_rep
      end do
      orthogonality = maxval(abs(matmul(transpose(d_wann(:,:,iop)), &
        d_wann(:,:,iop))-identity_matrix(size(channels))))
      if (.not. all(ieee_is_finite(d_wann(:,:,iop))) .or. &
          .not. ieee_is_finite(orthogonality)) then
        write(message,'(a,i0,a)') 'SAWF D_wann operation ', iop, &
          ' produced a non-finite matrix or residual'
        ok = .false.
        return
      end if
      if (orthogonality > 5.0d-10) then
        write(message,'(a,i0,a,es12.4)') 'SAWF D_wann operation ', iop, &
          ' is not unitary; residual=', orthogonality
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine build_sawf_wannier_representation


  subroutine real_p_representation(rotation, representation)
    real(8), intent(in) :: rotation(3,3)
    real(8), intent(out) :: representation(3,3)
    real(8) :: basis(3,3)

    ! Columns are the Wannier90 mr=1,2,3 directions (z,x,y).
    basis = 0.0d0
    basis(3,1) = 1.0d0
    basis(1,2) = 1.0d0
    basis(2,3) = 1.0d0
    representation = matmul(transpose(basis),matmul(rotation,basis))
  end subroutine real_p_representation


  subroutine real_d_representation(rotation, representation)
    real(8), intent(in) :: rotation(3,3)
    real(8), intent(out) :: representation(5,5)
    real(8) :: basis(3,3,5), transformed(3,3)
    integer :: i, j

    call real_d_tensor_basis(basis)
    do j = 1, 5
      transformed = matmul(transpose(rotation),matmul(basis(:,:,j),rotation))
      do i = 1, 5
        ! Wannier90 convention: row j, column i is
        ! <Rotated Y(j)|Y(i)>, as written by pw2wannier90 compute_dmn.
        representation(j,i) = sum(basis(:,:,i)*transformed)
      end do
    end do
  end subroutine real_d_representation


  subroutine real_d_tensor_basis(basis)
    real(8), intent(out) :: basis(3,3,5)
    real(8), parameter :: inv_sqrt2 = 1.0d0/sqrt(2.0d0)
    real(8), parameter :: inv_sqrt6 = 1.0d0/sqrt(6.0d0)

    basis = 0.0d0
    basis(1,1,1)=-inv_sqrt6; basis(2,2,1)=-inv_sqrt6
    basis(3,3,1)=2.0d0*inv_sqrt6
    basis(3,1,2)=inv_sqrt2; basis(1,3,2)=inv_sqrt2
    basis(2,3,3)=inv_sqrt2; basis(3,2,3)=inv_sqrt2
    basis(1,1,4)=inv_sqrt2; basis(2,2,4)=-inv_sqrt2
    basis(1,2,5)=inv_sqrt2; basis(2,1,5)=inv_sqrt2
  end subroutine real_d_tensor_basis


  function identity_matrix(n) result(identity)
    integer, intent(in) :: n
    real(8) :: identity(n,n)
    integer :: i

    identity = 0.0d0
    do i = 1, n
      identity(i,i) = 1.0d0
    end do
  end function identity_matrix

  subroutine load_sawf_symmetry_auto(lattice, frac_pos, species, tolerance, &
      operations, ok, message)
    real(8), intent(in) :: lattice(3,3)
    real(8), intent(in) :: frac_pos(:,:)
    integer, intent(in) :: species(:)
    real(8), intent(in) :: tolerance
    type(t_sawf_symop), allocatable, intent(out) :: operations(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
#ifdef HAVE_SPGLIB
    real(8), allocatable :: raw(:,:,:)
    real(8) :: lattice_inv(3,3), metric(3,3), metric_scale
    integer(c_int), allocatable :: c_species(:), c_rotations(:,:,:)
    real(c_double), allocatable :: c_translations(:,:)
    integer(c_int) :: dummy_rotation(1), operation_count, status
    real(c_double) :: dummy_translation(1)
    integer :: allocation_status, iop
#endif

    ok = .false.
    message = ''
#ifndef HAVE_SPGLIB
    message = 'SAWF automatic symmetry loading requires SALMON built with USE_SPGLIB=ON'
    return
#else
    call prepare_sawf_validation(lattice, frac_pos, species, tolerance, &
      lattice_inv, metric, metric_scale, ok, message)
    if (.not. ok) return
    if (any(species <= 0)) then
      message = 'SAWF species identifiers must be positive for Spglib'
      ok = .false.
      return
    end if
    if (size(species,kind=int64) > int(huge(0_c_int),int64)) then
      message = 'SAWF atom count exceeds the Spglib integer range'
      ok = .false.
      return
    end if
    if (any(int(species,int64) > int(huge(0_c_int),int64))) then
      message = 'SAWF species identifiers exceed the Spglib integer range'
      ok = .false.
      return
    end if

    allocate(c_species(size(species)), stat=allocation_status)
    if (allocation_status /= 0) then
      message = 'SAWF Spglib species allocation failed'
      ok = .false.
      return
    end if
    c_species = int(species,c_int)
    status = salmon_sawf_spglib_get_symmetry(lattice, frac_pos, c_species, &
      int(size(species),c_int), real(tolerance,c_double), 0_c_int, &
      dummy_rotation, dummy_translation, operation_count)
    if (status /= sawf_spglib_success) then
      call set_spglib_error(status, 'operation-count query', message)
      ok = .false.
      return
    end if
    if (operation_count <= 0_c_int) then
      message = 'Spglib returned an invalid SAWF symmetry operation count'
      ok = .false.
      return
    end if
    if (int(operation_count,int64) > int(max_sawf_symmetry_operations,int64)) then
      write(message,'(a,i0,a,i0)') 'Spglib operation count ', operation_count, &
        ' exceeds SAWF operation limit ', max_sawf_symmetry_operations
      ok = .false.
      return
    end if

    allocate(c_rotations(3,3,operation_count), &
      c_translations(3,operation_count), raw(3,4,operation_count), &
      stat=allocation_status)
    if (allocation_status /= 0) then
      message = 'SAWF Spglib operation allocation failed'
      ok = .false.
      return
    end if
    status = salmon_sawf_spglib_get_symmetry(lattice, frac_pos, c_species, &
      int(size(species),c_int), real(tolerance,c_double), operation_count, &
      c_rotations, c_translations, operation_count)
    if (status /= sawf_spglib_success) then
      call set_spglib_error(status, 'operation generation', message)
      ok = .false.
      return
    end if
    if (size(raw,3) /= operation_count) then
      message = 'Spglib changed the SAWF symmetry operation count'
      ok = .false.
      return
    end if
    do iop = 1, size(raw,3)
      raw(:,1:3,iop) = real(c_rotations(:,:,iop),8)
      raw(:,4,iop) = real(c_translations(:,iop),8)
    end do
    call normalize_sawf_operations(raw, lattice, lattice_inv, metric, &
      metric_scale, frac_pos, species, tolerance, operations, ok, message)
#endif
  end subroutine load_sawf_symmetry_auto

  subroutine load_sawf_symmetry_file(filename, lattice, frac_pos, species, &
      tolerance, operations, ok, message)
    character(*), intent(in) :: filename
    real(8), intent(in) :: lattice(3,3)
    real(8), intent(in) :: frac_pos(:,:)
    integer, intent(in) :: species(:)
    real(8), intent(in) :: tolerance
    type(t_sawf_symop), allocatable, intent(out) :: operations(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8), allocatable :: raw(:,:,:)
    real(8) :: lattice_inv(3,3), metric(3,3), metric_scale
    logical :: read_ok
    character(512) :: read_message

    ok = .false.
    message = ''
    call prepare_sawf_validation(lattice, frac_pos, species, tolerance, &
      lattice_inv, metric, metric_scale, ok, message)
    if (.not. ok) return

    call read_symmetry_file(filename, raw, read_ok, read_message)
    if (.not. read_ok) then
      message = trim(read_message)
      ok = .false.
      return
    end if
    call normalize_sawf_operations(raw, lattice, lattice_inv, metric, &
      metric_scale, frac_pos, species, tolerance, operations, ok, message)
  end subroutine load_sawf_symmetry_file


  subroutine prepare_sawf_validation(lattice, frac_pos, species, tolerance, &
      lattice_inv, metric, metric_scale, ok, message)
    real(8), intent(in) :: lattice(3,3), frac_pos(:,:), tolerance
    integer, intent(in) :: species(:)
    real(8), intent(out) :: lattice_inv(3,3), metric(3,3), metric_scale
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8) :: determinant

    ok = .false.
    if (size(frac_pos,1) /= 3 .or. size(frac_pos,2) /= size(species)) then
      message = 'SAWF atom arrays have inconsistent dimensions'
      return
    end if
    if (size(species) == 0) then
      message = 'SAWF symmetry validation requires at least one atom'
      return
    end if
    if (.not. ieee_is_finite(tolerance)) then
      message = 'SAWF symmetry tolerance must be finite'
      return
    end if
    if (.not. all(ieee_is_finite(lattice))) then
      message = 'SAWF lattice entries must be finite'
      return
    end if
    if (.not. all(ieee_is_finite(frac_pos))) then
      message = 'SAWF fractional atom positions must be finite'
      return
    end if
    if (tolerance <= 0.0d0) then
      message = 'SAWF symmetry tolerance must be positive'
      return
    end if
    call inverse_3x3(lattice, lattice_inv, determinant, ok)
    if (.not. ok) then
      write(message,'(a,es12.4)') 'SAWF lattice is singular; determinant=', determinant
      return
    end if
    metric = matmul(transpose(lattice), lattice)
    metric_scale = max(1.0d0, maxval(abs(metric)))
  end subroutine prepare_sawf_validation


  subroutine normalize_sawf_operations(raw, lattice, lattice_inv, metric, &
      metric_scale, frac_pos, species, tolerance, operations, ok, message)
    real(8), intent(in) :: raw(:,:,:), lattice(3,3), lattice_inv(3,3)
    real(8), intent(in) :: metric(3,3), metric_scale, frac_pos(:,:), tolerance
    integer, intent(in) :: species(:)
    type(t_sawf_symop), allocatable, intent(out) :: operations(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: identity_index, iop
    type(t_sawf_operation_index) :: operation_index

    ok = .false.
    if (size(raw,1) /= 3 .or. size(raw,2) /= 4 .or. size(raw,3) == 0) then
      message = 'SAWF symmetry source returned invalid operation dimensions'
      return
    end if
    if (size(raw,3) > max_sawf_symmetry_operations) then
      write(message,'(a,i0,a,i0)') 'SAWF symmetry operation count ', size(raw,3), &
        ' exceeds operation limit ', max_sawf_symmetry_operations
      return
    end if
    allocate(operations(size(raw,3)))
    do iop = 1, size(operations)
      call normalize_operation(iop, raw(:,:,iop), lattice, lattice_inv, metric, &
        metric_scale, frac_pos, species, tolerance, operations(iop), ok, message)
      if (.not. ok) return
    end do

    identity_index = 0
    do iop = 1, size(operations)
      if (is_identity(operations(iop), lattice, lattice_inv, tolerance)) then
        identity_index = iop
        exit
      end if
    end do
    if (identity_index == 0) then
      message = 'SAWF symmetry set is missing the identity operation'
      return
    end if
    call validate_inverses(operations, lattice, lattice_inv, tolerance, ok, message)
    if (.not. ok) return
    call build_operation_index(operations, operation_index, ok, message)
    if (.not. ok) return
    call validate_closure(operations, operation_index, lattice, lattice_inv, &
      tolerance, ok, message)
  end subroutine normalize_sawf_operations


#ifdef HAVE_SPGLIB
  subroutine set_spglib_error(status, stage, message)
    integer(c_int), intent(in) :: status
    character(*), intent(in) :: stage
    character(*), intent(out) :: message

    select case (status)
    case (sawf_spglib_count_failure)
      message = 'Spglib C API returned an invalid symmetry operation count'
    case (sawf_spglib_capacity_too_small)
      message = 'SAWF Spglib output capacity is too small'
    case (sawf_spglib_api_failure)
      message = 'Spglib C API failed during '//trim(stage)
    case (sawf_spglib_allocation_failure)
      message = 'SAWF Spglib C wrapper allocation failed during '//trim(stage)
    case default
      message = 'SAWF Spglib C wrapper rejected inputs during '//trim(stage)
    end select
  end subroutine set_spglib_error
#endif


  subroutine normalize_operation(iop, raw, lattice, lattice_inv, metric, metric_scale, &
      frac_pos, species, tolerance, operation, ok, message)
    integer, intent(in) :: iop
    real(8), intent(in) :: raw(3,4), lattice(3,3), lattice_inv(3,3)
    real(8), intent(in) :: metric(3,3), metric_scale
    real(8), intent(in) :: frac_pos(:,:)
    integer, intent(in) :: species(:)
    real(8), intent(in) :: tolerance
    type(t_sawf_symop), intent(out) :: operation
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8) :: atom_residual
    integer(int64) :: W_work(3,3), determinant_integer
    logical :: map_ok, determinant_ok, periodic_search_ok

    ok = .false.
    if (.not. all(ieee_is_finite(raw))) then
      write(message,'(a,i0,a)') 'SAWF operation ', iop, &
        ' contains a non-finite rotation or translation value'
      return
    end if
    if (any(abs(raw(:,1:3)) > real(huge(0),8))) then
      write(message,'(a,i0,a)') 'SAWF operation ', iop, &
        ' rotation is not representable as default integer'
      return
    end if
    W_work = nint(raw(:,1:3), kind=int64)
    operation%integer_residual = maxval(abs(raw(:,1:3) - real(W_work,8)))
    if (operation%integer_residual > tolerance) then
      write(message,'(a,i0,a,es12.4,a,es12.4)') 'SAWF operation ', iop, &
        ' rotation is not near integer; residual=', operation%integer_residual, &
        ' tolerance=', tolerance
      return
    end if
    operation%W = int(W_work)

    call determinant_integer_3x3(W_work, determinant_integer, determinant_ok)
    if (.not. determinant_ok) then
      write(message,'(a,i0,a)') 'SAWF operation ', iop, &
        ' determinant intermediate exceeds int64 range'
      return
    end if
    operation%determinant_residual = abs(abs(real(determinant_integer,8)) - 1.0d0)
    if (determinant_integer /= 1_int64 .and. determinant_integer /= -1_int64) then
      write(message,'(a,i0,a,i0,a,es12.4)') 'SAWF operation ', iop, &
        ' determinant is not +/-1; determinant=', determinant_integer, &
        ' residual=', operation%determinant_residual
      return
    end if

    operation%tau = modulo(raw(:,4), 1.0d0)
    where (abs(operation%tau - 1.0d0) <= tolerance)
      operation%tau = 0.0d0
    end where
    operation%R = matmul(lattice, matmul(real(operation%W,8), lattice_inv))
    operation%metric_residual = maxval(abs( &
      matmul(transpose(real(operation%W,8)), matmul(metric, real(operation%W,8))) - metric))
    if (operation%metric_residual > tolerance * metric_scale) then
      write(message,'(a,i0,a,es12.4,a,es12.4)') 'SAWF operation ', iop, &
        ' violates lattice metric; residual=', operation%metric_residual, &
        ' tolerance=', tolerance * metric_scale
      return
    end if

    operation%orthogonality_residual = maxval(abs( &
      matmul(transpose(operation%R), operation%R) - identity_real()))
    if (operation%orthogonality_residual > tolerance) then
      write(message,'(a,i0,a,es12.4,a,es12.4)') 'SAWF operation ', iop, &
        ' Cartesian rotation is not orthogonal; residual=', &
        operation%orthogonality_residual, ' tolerance=', tolerance
      return
    end if

    call build_atom_map(operation%W, operation%tau, lattice, lattice_inv, &
      frac_pos, species, tolerance, operation%atom_map, atom_residual, &
      map_ok, periodic_search_ok)
    operation%atom_map_residual = atom_residual
    if (.not. periodic_search_ok) then
      write(message,'(a,i0,a)') 'SAWF operation ', iop, &
        ' lattice too ill-conditioned for exact periodic search'
      return
    end if
    if (.not. map_ok) then
      write(message,'(a,i0,a,es12.4,a,es12.4)') 'SAWF operation ', iop, &
        ' has no one-to-one same-species periodic atom map; residual=', &
        operation%atom_map_residual, ' tolerance=', tolerance
      return
    end if
    ok = .true.
  end subroutine normalize_operation


  subroutine build_atom_map(W, tau, lattice, lattice_inv, frac_pos, species, tolerance, &
      atom_map, max_residual, ok, periodic_search_ok)
    integer, intent(in) :: W(3,3)
    real(8), intent(in) :: tau(3), lattice(3,3), lattice_inv(3,3), frac_pos(:,:)
    integer, intent(in) :: species(:)
    real(8), intent(in) :: tolerance
    integer, allocatable, intent(out) :: atom_map(:)
    real(8), intent(out) :: max_residual
    logical, intent(out) :: ok, periodic_search_ok
    logical, allocatable :: visited_target(:)
    real(8), allocatable :: distance(:,:)
    real(8) :: transformed(3), delta(3), residual, rejected_residual
    integer, allocatable :: target_owner(:)
    integer :: iatom, jatom, ispecies
    logical :: distance_ok

    allocate(atom_map(size(species)), target_owner(size(species)))
    allocate(visited_target(size(species)), distance(size(species),size(species)))
    atom_map = 0
    target_owner = 0
    distance = huge(1.0d0)
    max_residual = 0.0d0
    ok = .false.
    periodic_search_ok = .true.
    do iatom = 1, size(species)
      transformed = matmul(real(W,8), frac_pos(:,iatom)) + tau
      do jatom = 1, size(species)
        if (species(jatom) /= species(iatom)) cycle
        delta = transformed - frac_pos(:,jatom)
        call closest_periodic_residual(delta, lattice, lattice_inv, &
          distance(iatom,jatom), distance_ok)
        if (.not. distance_ok) then
          periodic_search_ok = .false.
          return
        end if
      end do
    end do

    do ispecies = 1, size(species)
      if (any(species(1:ispecies-1) == species(ispecies))) cycle
      do iatom = 1, size(species)
        if (species(iatom) /= species(ispecies)) cycle
        visited_target = .false.
        rejected_residual = huge(1.0d0)
        if (.not. augment_matching(iatom)) then
          max_residual = rejected_residual
          return
        end if
      end do
    end do

    do jatom = 1, size(species)
      if (target_owner(jatom) > 0) atom_map(target_owner(jatom)) = jatom
    end do
    do iatom = 1, size(species)
      max_residual = max(max_residual, distance(iatom,atom_map(iatom)))
    end do
    ok = .true.

  contains

    recursive logical function augment_matching(source_atom) result(found)
      integer, intent(in) :: source_atom
      integer :: candidate, owner, target
      real(8) :: candidate_distance

      found = .false.
      do target = 1, size(species)
        if (species(target) /= species(source_atom)) cycle
        if (distance(source_atom,target) > tolerance) &
          rejected_residual = min(rejected_residual, distance(source_atom,target))
      end do
      do
        candidate = 0
        candidate_distance = huge(1.0d0)
        do target = 1, size(species)
          if (visited_target(target)) cycle
          if (species(target) /= species(source_atom)) cycle
          if (distance(source_atom,target) > tolerance) cycle
          if (distance(source_atom,target) < candidate_distance) then
            candidate = target
            candidate_distance = distance(source_atom,target)
          end if
        end do
        if (candidate == 0) return

        visited_target(candidate) = .true.
        owner = target_owner(candidate)
        if (owner == 0) then
          target_owner(candidate) = source_atom
          found = .true.
          return
        end if
        if (augment_matching(owner)) then
          target_owner(candidate) = source_atom
          found = .true.
          return
        end if
      end do
    end function augment_matching
  end subroutine build_atom_map


  subroutine validate_inverses(operations, lattice, lattice_inv, tolerance, ok, message)
    type(t_sawf_symop), intent(in) :: operations(:)
    real(8), intent(in) :: lattice(3,3), lattice_inv(3,3), tolerance
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: iop, jop
    integer(int64) :: product_W(3,3)
    real(8) :: residual, best
    logical :: product_ok, residual_ok

    ok = .false.
    do iop = 1, size(operations)
      best = huge(1.0d0)
      do jop = 1, size(operations)
        call multiply_integer_matrices(int(operations(iop)%W,int64), &
          int(operations(jop)%W,int64), product_W, product_ok)
        if (.not. product_ok) then
          write(message,'(a,i0,a,i0)') 'SAWF inverse integer product overflows for operations ', &
            iop, ' and ', jop
          return
        end if
        if (all(product_W == identity_integer64())) then
          call translation_residual(operations(iop)%tau + &
            matmul(real(operations(iop)%W,8), operations(jop)%tau), &
            lattice, lattice_inv, residual, residual_ok)
          if (.not. residual_ok) then
            message = 'SAWF lattice too ill-conditioned for exact periodic search'
            return
          end if
          best = min(best, residual)
          if (residual <= tolerance) exit
        end if
      end do
      if (jop > size(operations)) then
        write(message,'(a,i0,a,es12.4,a,es12.4)') 'SAWF operation ', iop, &
          ' has no inverse; residual=', best, ' tolerance=', tolerance
        return
      end if
    end do
    ok = .true.
  end subroutine validate_inverses


  subroutine build_operation_index(operations, operation_index, ok, message)
    type(t_sawf_symop), intent(in) :: operations(:)
    type(t_sawf_operation_index), intent(out) :: operation_index
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer(int64) :: minimum_size, table_size, doubled_size
    integer :: allocation_status, existing, iop, probe, slot
    logical :: arithmetic_ok

    ok = .false.
    call checked_multiply_int64(int(size(operations),int64), 2_int64, &
      minimum_size, arithmetic_ok)
    if (.not. arithmetic_ok) then
      message = 'SAWF operation-index size calculation overflowed'
      return
    end if
    table_size = 1_int64
    do while (table_size < minimum_size)
      call checked_multiply_int64(table_size, 2_int64, doubled_size, arithmetic_ok)
      if (.not. arithmetic_ok .or. doubled_size > int(huge(0),int64)) then
        message = 'SAWF operation-index allocation size is not representable'
        return
      end if
      table_size = doubled_size
    end do
    allocate(operation_index%slots(int(table_size)), stat=allocation_status)
    if (allocation_status /= 0) then
      message = 'SAWF operation-index allocation failed'
      return
    end if
    operation_index%slots = 0

    do iop = 1, size(operations)
      slot = operation_hash_slot(int(operations(iop)%W,int64), &
        operations(iop)%atom_map(1), size(operation_index%slots))
      do probe = 1, size(operation_index%slots)
        existing = operation_index%slots(slot)
        if (existing == 0) then
          operation_index%slots(slot) = iop
          exit
        end if
        if (operations(existing)%atom_map(1) == operations(iop)%atom_map(1) .and. &
            all(operations(existing)%W == operations(iop)%W)) then
          write(message,'(a,i0,a,i0)') 'SAWF operations ', existing, ' and ', iop, &
            ' have a duplicate rotation/atom-map index key'
          return
        end if
        slot = slot + 1
        if (slot > size(operation_index%slots)) slot = 1
      end do
      if (probe > size(operation_index%slots)) then
        message = 'SAWF operation index is unexpectedly full'
        return
      end if
    end do
    ok = .true.
  end subroutine build_operation_index


  integer function operation_hash_slot(W, mapped_atom, table_size) result(slot)
    integer(int64), intent(in) :: W(3,3)
    integer, intent(in) :: mapped_atom, table_size
    integer(int64) :: hash_value
    integer :: i, j

    hash_value = 7640891576956012809_int64
    do j = 1, 3
      do i = 1, 3
        hash_value = ieor(hash_value, W(i,j))
        hash_value = ior(ishft(hash_value, 11), shiftr(hash_value, 53))
      end do
    end do
    hash_value = ieor(hash_value, int(mapped_atom,int64))
    hash_value = ior(ishft(hash_value, 11), shiftr(hash_value, 53))
    slot = int(modulo(hash_value, int(table_size,int64))) + 1
  end function operation_hash_slot


  subroutine find_indexed_operation(operations, operation_index, W, mapped_atom, &
      operation_number)
    type(t_sawf_symop), intent(in) :: operations(:)
    type(t_sawf_operation_index), intent(in) :: operation_index
    integer(int64), intent(in) :: W(3,3)
    integer, intent(in) :: mapped_atom
    integer, intent(out) :: operation_number
    integer :: candidate, probe, slot

    operation_number = 0
    slot = operation_hash_slot(W, mapped_atom, size(operation_index%slots))
    do probe = 1, size(operation_index%slots)
      candidate = operation_index%slots(slot)
      if (candidate == 0) return
      if (operations(candidate)%atom_map(1) == mapped_atom .and. &
          all(int(operations(candidate)%W,int64) == W)) then
        operation_number = candidate
        return
      end if
      slot = slot + 1
      if (slot > size(operation_index%slots)) slot = 1
    end do
  end subroutine find_indexed_operation


  subroutine validate_closure(operations, operation_index, lattice, lattice_inv, &
      tolerance, ok, message)
    type(t_sawf_symop), intent(in) :: operations(:)
    type(t_sawf_operation_index), intent(in) :: operation_index
    real(8), intent(in) :: lattice(3,3), lattice_inv(3,3), tolerance
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: iop, jop, kop, mapped_atom
    integer(int64) :: product_W(3,3)
    real(8) :: product_tau(3), residual
    logical :: product_ok, residual_ok

    ok = .false.
    do iop = 1, size(operations)
      do jop = 1, size(operations)
        call multiply_integer_matrices(int(operations(iop)%W,int64), &
          int(operations(jop)%W,int64), product_W, product_ok)
        if (.not. product_ok) then
          write(message,'(a,i0,a,i0)') 'SAWF closure integer product overflows for operations ', &
            iop, ' and ', jop
          return
        end if
        product_tau = operations(iop)%tau + &
          matmul(real(operations(iop)%W,8), operations(jop)%tau)
        mapped_atom = operations(iop)%atom_map(operations(jop)%atom_map(1))
        call find_indexed_operation(operations, operation_index, product_W, &
          mapped_atom, kop)
        if (kop == 0) then
          write(message,'(a,i0,a,i0,a,es12.4,a,es12.4)') &
            'SAWF group closure fails for operations ', iop, ' and ', jop, &
            '; residual=', huge(1.0d0), ' tolerance=', tolerance
          return
        end if
        call translation_residual(product_tau - operations(kop)%tau, &
          lattice, lattice_inv, residual, residual_ok)
        if (.not. residual_ok) then
          message = 'SAWF lattice too ill-conditioned for exact periodic search'
          return
        end if
        if (residual > tolerance) then
          write(message,'(a,i0,a,i0,a,es12.4,a,es12.4)') &
            'SAWF group closure fails for operations ', iop, ' and ', jop, &
            '; residual=', residual, ' tolerance=', tolerance
          return
        end if
      end do
    end do
    ok = .true.
  end subroutine validate_closure


  logical function is_identity(operation, lattice, lattice_inv, tolerance)
    type(t_sawf_symop), intent(in) :: operation
    real(8), intent(in) :: lattice(3,3), lattice_inv(3,3), tolerance
    real(8) :: residual
    logical :: residual_ok

    is_identity = .false.
    if (any(operation%W /= identity_integer())) return
    call translation_residual(operation%tau, lattice, lattice_inv, residual, residual_ok)
    if (.not. residual_ok) return
    is_identity = residual <= tolerance
  end function is_identity


  subroutine translation_residual(tau, lattice, lattice_inv, residual, ok)
    real(8), intent(in) :: tau(3), lattice(3,3), lattice_inv(3,3)
    real(8), intent(out) :: residual
    logical, intent(out) :: ok
    call closest_periodic_residual(tau, lattice, lattice_inv, residual, ok)
  end subroutine translation_residual


  subroutine closest_periodic_residual(delta, lattice, lattice_inv, residual, ok)
    real(8), intent(in) :: delta(3), lattice(3,3), lattice_inv(3,3)
    real(8), intent(out) :: residual
    logical, intent(out) :: ok
    real(8) :: reduced(3), candidate_delta(3), cartesian(3)
    real(8) :: best_squared, candidate_squared, radius, bound, tie_tolerance
    real(8) :: lower_real, upper_real, safe_integer_limit
    integer(int64) :: best_translation(3), candidate_translation(3)
    integer(int64) :: lower(3), upper(3), n1, n2, n3
    integer(int64) :: width, width_minus_one, candidate_count, updated_count
    integer :: idir
    logical :: arithmetic_ok

    ok = .false.
    residual = huge(1.0d0)
    if (.not. all(ieee_is_finite(delta))) return
    reduced = modulo(delta, 1.0d0)
    where (reduced > 0.5d0) reduced = reduced - 1.0d0
    best_translation = nint(reduced, kind=int64)
    candidate_delta = reduced - real(best_translation,8)
    cartesian = matmul(lattice, candidate_delta)
    best_squared = sum(cartesian**2)
    if (.not. ieee_is_finite(best_squared)) return
    radius = sqrt(max(0.0d0, best_squared))
    safe_integer_limit = real(huge(0_int64),8) * &
      (1.0d0 - 4.0d0 * epsilon(1.0d0))
    do idir = 1, 3
      bound = sqrt(sum(lattice_inv(idir,:)**2)) * radius
      bound = bound + 32.0d0 * epsilon(1.0d0) * max(1.0d0, bound)
      if (.not. ieee_is_finite(bound)) return
      lower_real = reduced(idir) - bound
      upper_real = reduced(idir) + bound
      if (.not. ieee_is_finite(lower_real) .or. .not. ieee_is_finite(upper_real)) return
      if (lower_real < -safe_integer_limit .or. upper_real > safe_integer_limit) then
        return
      end if
      lower(idir) = ceiling(lower_real, kind=int64)
      upper(idir) = floor(upper_real, kind=int64)
    end do

    candidate_count = 1_int64
    do idir = 1, 3
      call checked_add_int64(upper(idir), -lower(idir), width_minus_one, arithmetic_ok)
      if (.not. arithmetic_ok) return
      call checked_add_int64(width_minus_one, 1_int64, width, arithmetic_ok)
      if (.not. arithmetic_ok .or. width <= 0_int64) return
      call checked_multiply_int64(candidate_count, width, updated_count, arithmetic_ok)
      if (.not. arithmetic_ok) return
      if (updated_count > max_periodic_candidates) return
      candidate_count = updated_count
    end do

    tie_tolerance = 64.0d0 * epsilon(1.0d0) * max(1.0d0, best_squared)
    do n3 = lower(3), upper(3)
      do n2 = lower(2), upper(2)
        do n1 = lower(1), upper(1)
          candidate_translation = [n1, n2, n3]
          candidate_delta = reduced - real(candidate_translation,8)
          cartesian = matmul(lattice, candidate_delta)
          candidate_squared = sum(cartesian**2)
          if (candidate_squared < best_squared - tie_tolerance .or. &
              (abs(candidate_squared - best_squared) <= tie_tolerance .and. &
               lexicographically_less(candidate_translation, best_translation))) then
            best_squared = candidate_squared
            best_translation = candidate_translation
          end if
        end do
      end do
    end do
    residual = sqrt(max(0.0d0, best_squared))
    ok = .true.
  end subroutine closest_periodic_residual


  pure logical function lexicographically_less(left, right)
    integer(int64), intent(in) :: left(3), right(3)
    integer :: idir
    lexicographically_less = .false.
    do idir = 1, 3
      if (left(idir) < right(idir)) then
        lexicographically_less = .true.
        return
      else if (left(idir) > right(idir)) then
        return
      end if
    end do
  end function lexicographically_less


  pure function identity_integer() result(identity)
    integer :: identity(3,3)
    identity = 0
    identity(1,1) = 1
    identity(2,2) = 1
    identity(3,3) = 1
  end function identity_integer


  pure function identity_integer64() result(identity)
    integer(int64) :: identity(3,3)
    identity = 0_int64
    identity(1,1) = 1_int64
    identity(2,2) = 1_int64
    identity(3,3) = 1_int64
  end function identity_integer64


  subroutine determinant_integer_3x3(matrix, determinant, ok)
    integer(int64), intent(in) :: matrix(3,3)
    integer(int64), intent(out) :: determinant
    logical, intent(out) :: ok

    determinant = 0_int64
    ok = .true.
    call add_determinant_term(matrix(1,1), matrix(2,2), matrix(3,3), 1, determinant, ok)
    call add_determinant_term(matrix(1,2), matrix(2,3), matrix(3,1), 1, determinant, ok)
    call add_determinant_term(matrix(1,3), matrix(2,1), matrix(3,2), 1, determinant, ok)
    call add_determinant_term(matrix(1,3), matrix(2,2), matrix(3,1), -1, determinant, ok)
    call add_determinant_term(matrix(1,2), matrix(2,1), matrix(3,3), -1, determinant, ok)
    call add_determinant_term(matrix(1,1), matrix(2,3), matrix(3,2), -1, determinant, ok)
  end subroutine determinant_integer_3x3


  subroutine add_determinant_term(a, b, c, sign, accumulator, ok)
    integer(int64), intent(in) :: a, b, c
    integer, intent(in) :: sign
    integer(int64), intent(inout) :: accumulator
    logical, intent(inout) :: ok
    integer(int64) :: product, term, updated

    if (.not. ok) return
    call checked_multiply_int64(a, b, product, ok)
    if (.not. ok) return
    call checked_multiply_int64(product, c, term, ok)
    if (.not. ok) return
    if (sign < 0) then
      term = -term
    end if
    call checked_add_int64(accumulator, term, updated, ok)
    if (ok) accumulator = updated
  end subroutine add_determinant_term


  subroutine multiply_integer_matrices(left, right, product, ok)
    integer(int64), intent(in) :: left(3,3), right(3,3)
    integer(int64), intent(out) :: product(3,3)
    logical, intent(out) :: ok
    integer(int64) :: term, updated
    integer :: i, j, k

    product = 0_int64
    ok = .true.
    do j = 1, 3
      do i = 1, 3
        do k = 1, 3
          call checked_multiply_int64(left(i,k), right(k,j), term, ok)
          if (.not. ok) return
          call checked_add_int64(product(i,j), term, updated, ok)
          if (.not. ok) return
          product(i,j) = updated
        end do
      end do
    end do
  end subroutine multiply_integer_matrices


  subroutine checked_multiply_int64(left, right, product, ok)
    integer(int64), intent(in) :: left, right
    integer(int64), intent(out) :: product
    logical, intent(out) :: ok
    integer(int64), parameter :: max_int = huge(0_int64)
    integer(int64), parameter :: min_int = -max_int

    ok = .false.
    product = 0_int64
    if (left == 0_int64 .or. right == 0_int64) then
      ok = .true.
      return
    end if
    if (left > 0_int64) then
      if (right > 0_int64) then
        if (left > max_int/right) return
      else
        if (right < min_int/left) return
      end if
    else
      if (right > 0_int64) then
        if (left < min_int/right) return
      else
        if (left < max_int/right) return
      end if
    end if
    product = left * right
    ok = .true.
  end subroutine checked_multiply_int64


  subroutine checked_add_int64(left, right, sum_value, ok)
    integer(int64), intent(in) :: left, right
    integer(int64), intent(out) :: sum_value
    logical, intent(out) :: ok
    integer(int64), parameter :: max_int = huge(0_int64)
    integer(int64), parameter :: min_int = -max_int

    ok = .false.
    sum_value = 0_int64
    if (right > 0_int64) then
      if (left > max_int - right) return
    else if (right < 0_int64) then
      if (left < min_int - right) return
    end if
    sum_value = left + right
    ok = .true.
  end subroutine checked_add_int64


  pure function identity_real() result(identity)
    real(8) :: identity(3,3)
    identity = 0.0d0
    identity(1,1) = 1.0d0
    identity(2,2) = 1.0d0
    identity(3,3) = 1.0d0
  end function identity_real


  pure real(8) function determinant_3x3(matrix)
    real(8), intent(in) :: matrix(3,3)
    determinant_3x3 = matrix(1,1) * (matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)) &
      - matrix(1,2) * (matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)) &
      + matrix(1,3) * (matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1))
  end function determinant_3x3


  subroutine inverse_3x3(matrix, inverse, determinant, ok)
    real(8), intent(in) :: matrix(3,3)
    real(8), intent(out) :: inverse(3,3), determinant
    logical, intent(out) :: ok

    determinant = determinant_3x3(matrix)
    ok = abs(determinant) > tiny(1.0d0)
    if (.not. ok) then
      inverse = 0.0d0
      return
    end if
    inverse(1,1) = matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)
    inverse(1,2) = matrix(1,3)*matrix(3,2) - matrix(1,2)*matrix(3,3)
    inverse(1,3) = matrix(1,2)*matrix(2,3) - matrix(1,3)*matrix(2,2)
    inverse(2,1) = matrix(2,3)*matrix(3,1) - matrix(2,1)*matrix(3,3)
    inverse(2,2) = matrix(1,1)*matrix(3,3) - matrix(1,3)*matrix(3,1)
    inverse(2,3) = matrix(1,3)*matrix(2,1) - matrix(1,1)*matrix(2,3)
    inverse(3,1) = matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1)
    inverse(3,2) = matrix(1,2)*matrix(3,1) - matrix(1,1)*matrix(3,2)
    inverse(3,3) = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)
    inverse = inverse / determinant
  end subroutine inverse_3x3

end module lcfo_wannier_sawf
