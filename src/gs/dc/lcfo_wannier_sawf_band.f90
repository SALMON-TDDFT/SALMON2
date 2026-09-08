module lcfo_wannier_sawf_band
  ! D_band construction in this module is consumed by the scalable-template
  ! gate only after the full operation set passes representation closure.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use lcfo_wannier_sawf, only: t_sawf_symop
  implicit none
  private

  public :: map_sawf_periodic_grid_point, build_sawf_periodic_grid_map
  public :: locate_sawf_fragment_point, validate_sawf_fragment_tiling
  public :: validate_sawf_fragment_symmetry_map
  public :: build_sawf_fragment_buffer_point_map
  public :: accumulate_sawf_dmn_band, accumulate_sawf_dmn_band_blocks
  public :: diagnose_sawf_fragment_basis_closure
  public :: validate_sawf_dmn_band, validate_sawf_representation_product
  public :: validate_sawf_operation_set_products
  public :: validate_sawf_seed_header, validate_sawf_seed_basis_metadata

contains

  subroutine validate_sawf_seed_header(expected_fragments, requested_fragment, requested_bands, &
      file_fragments, file_spins, file_fragment_states, file_total_states, ok, message)
    integer, intent(in) :: expected_fragments,requested_fragment,requested_bands
    integer, intent(in) :: file_fragments,file_spins,file_fragment_states,file_total_states
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer(int64) :: count64,coefficient_count64

    ok=.false.; message=''
    if(expected_fragments <= 0 .or. requested_fragment < 1 .or. &
        requested_fragment > expected_fragments .or. requested_bands <= 0) then
      message='SAWF seed request metadata is invalid'
      return
    end if
    if(file_fragments /= expected_fragments) then
      message='SAWF seed fragment count does not match the current DC decomposition'
      return
    end if
    if(file_spins <= 0 .or. file_fragment_states <= 0 .or. file_total_states <= 0 .or. &
        file_total_states < requested_bands) then
      message='SAWF seed state dimensions are invalid or insufficient'
      return
    end if
    count64=int(file_fragment_states,int64)*int(file_fragments,int64)
    if(count64 > huge(0_int64)/int(file_spins,int64)) then
      message='SAWF seed index metadata size overflows int64'
      return
    end if
    count64=count64*int(file_spins,int64)
    if(count64 > int(huge(0),int64)) then
      message='SAWF seed index metadata size exceeds default integer indexing'
      return
    end if
    if(int(file_fragment_states,int64) > huge(0_int64)/int(file_total_states,int64)) then
      message='SAWF seed coefficient metadata size overflows int64'
      return
    end if
    coefficient_count64=int(file_fragment_states,int64)*int(file_total_states,int64)
    if(coefficient_count64 > huge(0_int64)/int(file_spins,int64)) then
      message='SAWF seed coefficient metadata size overflows int64'
      return
    end if
    coefficient_count64=coefficient_count64*int(file_spins,int64)
    if(coefficient_count64 > int(huge(0),int64) .or. &
        int(file_total_states,int64)*int(file_spins,int64) > int(huge(0),int64)) then
      message='SAWF seed coefficient metadata exceeds default integer indexing'
      return
    end if
    ok=.true.
  end subroutine validate_sawf_seed_header


  subroutine validate_sawf_seed_basis_metadata(file_fragment_states,file_total_states, &
      n_mat,n_basis,ok,message,index_basis)
    integer, intent(in) :: file_fragment_states,file_total_states
    integer, intent(in) :: n_mat(:),n_basis(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, intent(in), optional :: index_basis(:,:,:)
    integer(int64) :: basis_sum64
    integer :: ifrag,ispin

    ok=.false.; message=''
    if(size(n_mat) <= 0 .or. size(n_basis,1) <= 0 .or. size(n_basis,2) /= size(n_mat)) then
      message='SAWF seed basis metadata dimensions are inconsistent'
      return
    end if
    if(file_total_states <= 0 .or. any(n_mat <= 0)) then
      message='SAWF seed global matrix dimensions are invalid'
      return
    end if
    if(any(n_basis < 0) .or. any(n_basis > file_fragment_states)) then
      message='SAWF seed basis count exceeds the fragment-state dimension'
      return
    end if
    do ispin=1,size(n_mat)
      basis_sum64=0_int64
      do ifrag=1,size(n_basis,1)
        if(basis_sum64 > huge(0_int64)-int(n_basis(ifrag,ispin),int64)) then
          message='SAWF seed global basis count overflows int64'
          return
        end if
        basis_sum64=basis_sum64+int(n_basis(ifrag,ispin),int64)
      end do
      if(basis_sum64 > int(huge(0),int64) .or. int(n_mat(ispin),int64) /= basis_sum64) then
        message='SAWF seed n_mat does not equal the global basis count for its spin'
        return
      end if
    end do
    if(present(index_basis)) then
      if(size(index_basis,1) /= file_fragment_states .or. &
          size(index_basis,2) /= size(n_basis,1) .or. &
          size(index_basis,3) /= size(n_mat)) then
        message='SAWF seed global basis-index dimensions are inconsistent'
        return
      end if
      do ispin=1,size(n_mat)
        if(any(index_basis(:,:,ispin) < 0) .or. &
            any(index_basis(:,:,ispin) > n_mat(ispin))) then
          message='SAWF seed global basis index exceeds n_mat for its spin'
          return
        end if
      end do
    end if
    ok=.true.
  end subroutine validate_sawf_seed_basis_metadata

  subroutine prepare_grid_operation(operation, mesh, tolerance, coefficient, shift, ok, message)
    type(t_sawf_symop), intent(in) :: operation
    integer, intent(in) :: mesh(3)
    real(8), intent(in) :: tolerance
    integer(int64), intent(out) :: coefficient(3,3), shift(3)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer(int64) :: numerator
    integer :: axis, source_axis
    real(8) :: scaled_shift, wrapped_tau

    ok = .false.
    message = ''
    coefficient = 0_int64
    shift = 0_int64
    if (any(mesh <= 0)) then
      message = 'SAWF symmetry grid mesh dimensions must be positive'
      return
    end if
    if (.not. ieee_is_finite(tolerance) .or. tolerance <= 0.0d0) then
      message = 'SAWF symmetry grid tolerance must be finite and positive'
      return
    end if
    if (.not. all(ieee_is_finite(operation%tau))) then
      message = 'SAWF symmetry grid translation contains a non-finite value'
      return
    end if
    do axis = 1, 3
      do source_axis = 1, 3
        if (operation%W(axis,source_axis) /= 0 .and. &
            abs(int(operation%W(axis,source_axis),int64)) > &
            huge(0_int64)/int(mesh(axis),int64)) then
          message = 'SAWF symmetry grid integer rotation overflows mesh arithmetic'
          return
        end if
        numerator = int(operation%W(axis,source_axis),int64)*int(mesh(axis),int64)
        if (modulo(numerator,int(mesh(source_axis),int64)) /= 0_int64) then
          write(message,'(a,i0,a,i0,a)') 'SAWF symmetry operation axis ',axis, &
            ' from axis ',source_axis,' is incompatible with the mesh'
          return
        end if
        coefficient(axis,source_axis) = numerator/int(mesh(source_axis),int64)
      end do
      wrapped_tau = modulo(operation%tau(axis),1.0d0)
      scaled_shift = wrapped_tau*dble(mesh(axis))
      shift(axis) = nint(scaled_shift,kind=int64)
      if (abs(scaled_shift-dble(shift(axis))) > tolerance*dble(mesh(axis))) then
        write(message,'(a,i0,a,es13.5)') 'SAWF symmetry translation on axis ',axis, &
          ' is incompatible with the mesh; fractional grid shift=',scaled_shift
        return
      end if
      shift(axis) = modulo(shift(axis),int(mesh(axis),int64))
    end do
    ok = .true.
  end subroutine prepare_grid_operation


  subroutine map_sawf_periodic_grid_point(operation, mesh, tolerance, source_index, &
      target_index, ok, message)
    type(t_sawf_symop), intent(in) :: operation
    integer, intent(in) :: mesh(3), source_index(3)
    real(8), intent(in) :: tolerance
    integer, intent(out) :: target_index(3)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer(int64) :: coefficient(3,3), shift(3), source0(3), target0
    integer :: axis

    target_index = 0
    call prepare_grid_operation(operation,mesh,tolerance,coefficient,shift,ok,message)
    if (.not. ok) return
    if (any(source_index < 1) .or. any(source_index > mesh)) then
      message = 'SAWF symmetry source grid index is out of range'
      ok = .false.
      return
    end if
    source0 = int(source_index-1,int64)
    do axis = 1, 3
      target0 = shift(axis)+sum(coefficient(axis,:)*source0)
      target_index(axis) = 1+int(modulo(target0,int(mesh(axis),int64)))
    end do
    ok = .true.
  end subroutine map_sawf_periodic_grid_point


  subroutine checked_grid_size(mesh, npoints, ok, message)
    integer, intent(in) :: mesh(3)
    integer, intent(out) :: npoints
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer(int64) :: product64
    integer :: axis

    ok = .false.
    message = ''
    npoints = 0
    if (any(mesh <= 0)) then
      message = 'SAWF grid dimensions must be positive'
      return
    end if
    product64 = 1_int64
    do axis = 1, 3
      if (product64 > int(huge(0),int64)/int(mesh(axis),int64)) then
        message = 'SAWF grid point count overflows default integer'
        return
      end if
      product64 = product64*int(mesh(axis),int64)
    end do
    npoints = int(product64)
    ok = .true.
  end subroutine checked_grid_size


  subroutine build_sawf_periodic_grid_map(operation, mesh, tolerance, point_map, ok, message)
    type(t_sawf_symop), intent(in) :: operation
    integer, intent(in) :: mesh(3)
    real(8), intent(in) :: tolerance
    integer, allocatable, intent(out) :: point_map(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, allocatable :: multiplicity(:)
    integer :: allocation_status, ix, iy, iz, p, npoints
    integer :: source_index(3), target_index(3), target
    character(256) :: allocation_message

    call checked_grid_size(mesh,npoints,ok,message)
    if (.not. ok) return
    allocate(point_map(npoints),multiplicity(npoints),stat=allocation_status,errmsg=allocation_message)
    if (allocation_status /= 0) then
      message = 'SAWF symmetry grid-map allocation failed: '//trim(allocation_message)
      ok = .false.
      return
    end if
    multiplicity = 0
    p = 0
    do iz=1,mesh(3)
      do iy=1,mesh(2)
        do ix=1,mesh(1)
          source_index=[ix,iy,iz]
          call map_sawf_periodic_grid_point(operation,mesh,tolerance,source_index, &
            target_index,ok,message)
          if (.not. ok) then
            deallocate(point_map,multiplicity)
            return
          end if
          p=p+1
          target=target_index(1)+(target_index(2)-1)*mesh(1) + &
            (target_index(3)-1)*mesh(1)*mesh(2)
          point_map(p)=target
          multiplicity(target)=multiplicity(target)+1
        end do
      end do
    end do
    if (any(multiplicity /= 1)) then
      message = 'SAWF symmetry grid operation is not a one-to-one periodic permutation'
      deallocate(point_map,multiplicity)
      ok = .false.
      return
    end if
    deallocate(multiplicity)
    ok = .true.
  end subroutine build_sawf_periodic_grid_map


  subroutine locate_sawf_fragment_point(global_index, mesh, fragment_origin, fragment_shape, &
      owner_fragment, local_index, ok, message)
    integer, intent(in) :: global_index(3), mesh(3)
    integer, intent(in) :: fragment_origin(:,:), fragment_shape(:,:)
    integer, intent(out) :: owner_fragment, local_index(3)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: axis, candidate, count_owner, local0(3)

    ok = .false.
    message = ''
    owner_fragment = 0
    local_index = 0
    if (size(fragment_origin,1) /= 3 .or. size(fragment_shape,1) /= 3 .or. &
        size(fragment_origin,2) /= size(fragment_shape,2) .or. size(fragment_origin,2) <= 0) then
      message = 'SAWF fragment ownership arrays have inconsistent dimensions'
      return
    end if
    if (any(mesh <= 0) .or. any(global_index < 1) .or. any(global_index > mesh) .or. &
        any(fragment_shape <= 0) .or. any(fragment_shape > spread(mesh,2,size(fragment_shape,2)))) then
      message = 'SAWF fragment ownership geometry is invalid'
      return
    end if
    count_owner = 0
    do candidate=1,size(fragment_origin,2)
      do axis=1,3
        local0(axis)=modulo(global_index(axis)-1-fragment_origin(axis,candidate),mesh(axis))
      end do
      if(all(local0 < fragment_shape(:,candidate))) then
        count_owner=count_owner+1
        owner_fragment=candidate
        local_index=local0+1
      end if
    end do
    if(count_owner == 0) then
      message = 'SAWF fragment tiling has a missing global grid point'
      return
    else if(count_owner > 1) then
      message = 'SAWF fragment tiling has duplicate ownership of a global grid point'
      return
    end if
    ok = .true.
  end subroutine locate_sawf_fragment_point


  subroutine validate_sawf_fragment_tiling(mesh, fragment_origin, fragment_shape, ok, message)
    integer, intent(in) :: mesh(3), fragment_origin(:,:), fragment_shape(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: ix, iy, iz, owner, local_index(3), global_index(3)

    do iz=1,mesh(3)
      do iy=1,mesh(2)
        do ix=1,mesh(1)
          global_index=[ix,iy,iz]
          call locate_sawf_fragment_point(global_index,mesh,fragment_origin,fragment_shape, &
            owner,local_index,ok,message)
          if(.not.ok) return
        end do
      end do
    end do
    ok=.true.
    message=''
  end subroutine validate_sawf_fragment_tiling


  subroutine build_sawf_fragment_buffer_point_map(operation,mesh,source_origin,source_shape, &
      target_origin,target_shape,buffer,tolerance,point_map,ok,message)
    type(t_sawf_symop), intent(in) :: operation
    integer, intent(in) :: mesh(3),source_origin(3),source_shape(3)
    integer, intent(in) :: target_origin(3),target_shape(3),buffer(3)
    real(8), intent(in) :: tolerance
    integer, allocatable, intent(out) :: point_map(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, allocatable :: multiplicity(:)
    integer :: source_box(3),target_box(3),source_points,target_points
    integer :: source_index(3),target_index(3),target_local(3)
    integer :: ix,iy,iz,axis,candidate,matches,p,target,allocation_status
    character(256) :: detail,allocation_message

    ok=.false.; message=''
    if(allocated(point_map)) deallocate(point_map)
    if(any(mesh<=0) .or. any(source_shape<=0) .or. any(target_shape<=0) .or. any(buffer<0)) then
      message='SAWF fragment buffer-map geometry is invalid'
      return
    end if
    source_box=source_shape+2*buffer
    target_box=target_shape+2*buffer
    call checked_grid_size(source_box,source_points,ok,detail)
    if(.not.ok) then; message=trim(detail); return; end if
    call checked_grid_size(target_box,target_points,ok,detail)
    if(.not.ok) then; message=trim(detail); return; end if
    if(source_points/=target_points) then
      message='SAWF symmetry maps fragment buffers with unequal point counts'
      ok=.false.
      return
    end if
    allocate(point_map(source_points),multiplicity(target_points),stat=allocation_status, &
      errmsg=allocation_message)
    if(allocation_status/=0) then
      message='SAWF fragment buffer-map allocation failed: '//trim(allocation_message)
      if(allocated(point_map)) deallocate(point_map)
      if(allocated(multiplicity)) deallocate(multiplicity)
      ok=.false.
      return
    end if
    multiplicity=0; p=0
    do iz=1,source_box(3)
      do iy=1,source_box(2)
        do ix=1,source_box(1)
          source_index=1+modulo(source_origin+[ix,iy,iz]-1-buffer,mesh)
          call map_sawf_periodic_grid_point(operation,mesh,tolerance,source_index, &
            target_index,ok,detail)
          if(.not.ok) then
            message=trim(detail)
            deallocate(point_map,multiplicity)
            return
          end if
          do axis=1,3
            matches=0; target_local(axis)=0
            do candidate=1,target_box(axis)
              if(1+modulo(target_origin(axis)+candidate-1-buffer(axis),mesh(axis)) == &
                  target_index(axis)) then
                matches=matches+1
                target_local(axis)=candidate
              end if
            end do
            if(matches/=1) then
              write(message,'(a,i0,a,i0)') 'SAWF mapped buffer point has non-unique target on axis ', &
                axis,' matches=',matches
              deallocate(point_map,multiplicity)
              ok=.false.
              return
            end if
          end do
          p=p+1
          target=target_local(1)+(target_local(2)-1)*target_box(1) + &
            (target_local(3)-1)*target_box(1)*target_box(2)
          point_map(p)=target
          multiplicity(target)=multiplicity(target)+1
        end do
      end do
    end do
    if(any(multiplicity/=1)) then
      message='SAWF fragment buffer map is not a one-to-one permutation'
      deallocate(point_map,multiplicity)
      ok=.false.
      return
    end if
    deallocate(multiplicity)
    ok=.true.
  end subroutine build_sawf_fragment_buffer_point_map


  subroutine validate_sawf_fragment_symmetry_map(operation,mesh,fragment_origin, &
      fragment_shape,buffer,tolerance,grid_map_ok,fragment_map_ok,max_targets_per_source, &
      source_to_target,max_grid_residual,center_available,center_grid,message)
    type(t_sawf_symop), intent(in) :: operation
    integer, intent(in) :: mesh(3),fragment_origin(:,:),fragment_shape(:,:),buffer(3)
    real(8), intent(in) :: tolerance
    logical, intent(out) :: grid_map_ok,fragment_map_ok,center_available
    integer, intent(out) :: max_targets_per_source
    integer, allocatable, intent(out) :: source_to_target(:)
    real(8), intent(out) :: max_grid_residual,center_grid(3)
    character(*), intent(out) :: message
    integer(int64) :: coefficient(3,3),shift(3),source_count64,target_count64
    integer(int64) :: source0(3),target0(3)
    integer, allocatable :: target_multiplicity(:)
    logical, allocatable :: target_seen(:)
    integer :: allocation_status,axis,source_axis,nonzero_count,nfrag,npoints
    integer :: ifrag,jfrag,ix,iy,iz,distinct_targets,owner,local_index(3),target_index(3)
    logical :: ok
    character(256) :: detail
    real(8) :: scaled_value

    grid_map_ok=.false.
    fragment_map_ok=.false.
    max_targets_per_source=0
    max_grid_residual=0.0d0
    center_available=.false.
    center_grid=0.0d0
    message=''
    nfrag=size(fragment_origin,2)
    if(size(fragment_origin,1)/=3 .or. size(fragment_shape,1)/=3 .or. &
        size(fragment_shape,2)/=nfrag .or. nfrag<=0) then
      message='SAWF fragment symmetry map geometry dimensions are inconsistent'
      return
    end if
    if(any(buffer<0)) then
      message='SAWF fragment symmetry map buffer widths must be nonnegative'
      return
    end if
    call checked_grid_size(mesh,npoints,ok,detail)
    if(.not.ok) then
      message=trim(detail)
      return
    end if

    call prepare_grid_operation(operation,mesh,tolerance,coefficient,shift,ok,detail)
    if(.not.ok) then
      message=trim(detail)
      return
    end if

    ! Residual is measured in zero-based grid-index units.  The validated
    ! int64 affine coefficients avoid an unsafe default-integer NINT here.
    do axis=1,3
      do source_axis=1,3
        scaled_value=dble(operation%W(axis,source_axis))*dble(mesh(axis))/dble(mesh(source_axis))
        max_grid_residual=max(max_grid_residual, &
          abs(scaled_value-dble(coefficient(axis,source_axis))))
      end do
      scaled_value=modulo(operation%tau(axis),1.0d0)*dble(mesh(axis))
      max_grid_residual=max(max_grid_residual,min(abs(scaled_value-dble(shift(axis))), &
        abs(dble(mesh(axis))-abs(scaled_value-dble(shift(axis))))))
    end do

    ! The fragment halo is axis based, so only signed axis permutations are supported.
    do axis=1,3
      nonzero_count=count(coefficient(axis,:)/=0_int64)
      if(nonzero_count/=1 .or. maxval(abs(coefficient(axis,:)))/=1_int64) then
        message='SAWF fragment symmetry discrete map is not a signed-axis permutation'
        return
      end if
    end do
    do source_axis=1,3
      if(count(coefficient(:,source_axis)/=0_int64)/=1) then
        message='SAWF fragment symmetry discrete map is not a signed-axis permutation'
        return
      end if
    end do
    grid_map_ok=.true.

    if(all(operation%W==-reshape([1,0,0,0,1,0,0,0,1],[3,3]))) then
      center_available=.true.
      center_grid=0.5d0*dble(shift)
    end if

    ! Correctness-first O(nops*Ngrid*nfrag) owner scan in the runtime caller.
    ! The current 32^3/eight-fragment gate is small; cache this map later.
    call validate_sawf_fragment_tiling(mesh,fragment_origin,fragment_shape,ok,detail)
    if(.not.ok) then
      message=trim(detail)
      return
    end if
    do axis=1,3
      source_axis=maxloc(abs(coefficient(axis,:)),dim=1)
      if(buffer(axis)/=buffer(source_axis)) then
        write(message,'(a,i0,a,i0)') 'SAWF fragment symmetry buffer mismatch between axes ', &
          axis,' and ',source_axis
        return
      end if
    end do

    allocate(source_to_target(nfrag),stat=allocation_status)
    if(allocation_status/=0) then
      message='SAWF fragment symmetry map work allocation failed'
      return
    end if
    allocate(target_multiplicity(nfrag),stat=allocation_status)
    if(allocation_status/=0) then
      message='SAWF fragment symmetry map work allocation failed'
      deallocate(source_to_target)
      return
    end if
    allocate(target_seen(nfrag),stat=allocation_status)
    if(allocation_status/=0) then
      message='SAWF fragment symmetry map work allocation failed'
      deallocate(source_to_target,target_multiplicity)
      return
    end if
    source_to_target=0
    target_multiplicity=0
    do ifrag=1,nfrag
      source_count64=1_int64
      do axis=1,3
        if(fragment_shape(axis,ifrag)<=0 .or. &
            source_count64>huge(0_int64)/int(fragment_shape(axis,ifrag),int64)) then
          message='SAWF source fragment core point count is invalid or overflows int64'
          deallocate(source_to_target,target_multiplicity,target_seen)
          return
        end if
        source_count64=source_count64*int(fragment_shape(axis,ifrag),int64)
      end do
      target_seen=.false.
      do iz=0,fragment_shape(3,ifrag)-1
        do iy=0,fragment_shape(2,ifrag)-1
          do ix=0,fragment_shape(1,ifrag)-1
            source0(1)=modulo(int(fragment_origin(1,ifrag),int64)+int(ix,int64),int(mesh(1),int64))
            source0(2)=modulo(int(fragment_origin(2,ifrag),int64)+int(iy,int64),int(mesh(2),int64))
            source0(3)=modulo(int(fragment_origin(3,ifrag),int64)+int(iz,int64),int(mesh(3),int64))
            do axis=1,3
              target0(axis)=shift(axis)+sum(coefficient(axis,:)*source0)
              target_index(axis)=1+int(modulo(target0(axis),int(mesh(axis),int64)))
            end do
            call locate_sawf_fragment_point(target_index,mesh,fragment_origin,fragment_shape, &
              owner,local_index,ok,detail)
            if(.not.ok) then
              message=trim(detail)
              deallocate(source_to_target,target_multiplicity,target_seen)
              return
            end if
            target_seen(owner)=.true.
          end do
        end do
      end do
      distinct_targets=count(target_seen)
      max_targets_per_source=max(max_targets_per_source,distinct_targets)
      if(distinct_targets==1) then
        source_to_target(ifrag)=findloc(target_seen,.true.,dim=1)
        jfrag=source_to_target(ifrag)
        target_count64=1_int64
        do axis=1,3
          if(target_count64>huge(0_int64)/int(fragment_shape(axis,jfrag),int64)) then
            message='SAWF target fragment core point count overflows int64'
            deallocate(source_to_target,target_multiplicity,target_seen)
            return
          end if
          target_count64=target_count64*int(fragment_shape(axis,jfrag),int64)
        end do
        if(target_count64/=source_count64) then
          source_to_target(ifrag)=0
          message='SAWF symmetry maps between fragments with unequal core point counts'
        end if
      end if
    end do
    do ifrag=1,nfrag
      if(source_to_target(ifrag)>0) &
        target_multiplicity(source_to_target(ifrag))=target_multiplicity(source_to_target(ifrag))+1
    end do
    fragment_map_ok=all(source_to_target>0) .and. all(target_multiplicity==1)
    if(.not.fragment_map_ok .and. len_trim(message)==0) then
      write(message,'(a,i0)') 'SAWF symmetry splits a source core across target fragments; max_targets=', &
        max_targets_per_source
    end if
    deallocate(target_multiplicity,target_seen)
  end subroutine validate_sawf_fragment_symmetry_map


  subroutine accumulate_sawf_dmn_band(states, point_map, hvol, first_point, last_point, &
      d_band_contribution, ok, message)
    complex(8), intent(in) :: states(:,:)
    integer, intent(in) :: point_map(:)
    real(8), intent(in) :: hvol
    integer, intent(in) :: first_point, last_point
    complex(8), intent(out) :: d_band_contribution(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, allocatable :: source_points(:), target_points(:)
    integer :: allocation_status, count_points, i

    d_band_contribution=(0.0d0,0.0d0)
    if(size(states,1) /= size(point_map) .or. first_point < 1 .or. &
        last_point > size(point_map) .or. last_point < first_point-1) then
      message='SAWF D_band local point range or state dimensions are invalid'
      ok=.false.
      return
    end if
    count_points=max(0,last_point-first_point+1)
    allocate(source_points(count_points),target_points(count_points),stat=allocation_status)
    if(allocation_status /= 0) then
      message='SAWF D_band point-list allocation failed'
      ok=.false.
      return
    end if
    do i=1,count_points
      source_points(i)=first_point+i-1
      target_points(i)=point_map(source_points(i))
    end do
    call accumulate_sawf_dmn_band_blocks(states,states,source_points,target_points,hvol, &
      d_band_contribution,ok,message)
    deallocate(source_points,target_points)
  end subroutine accumulate_sawf_dmn_band


  subroutine accumulate_sawf_dmn_band_blocks(source_states, target_states, source_points, &
      target_points, hvol, d_band_contribution, ok, message)
    complex(8), intent(in) :: source_states(:,:), target_states(:,:)
    integer, intent(in) :: source_points(:), target_points(:)
    real(8), intent(in) :: hvol
    complex(8), intent(out) :: d_band_contribution(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: source_selected(:,:), target_selected(:,:)
    integer :: allocation_status, nband, npair

    ok=.false.
    message=''
    d_band_contribution=(0.0d0,0.0d0)
    nband=size(source_states,2)
    npair=size(source_points)
    if(.not.ieee_is_finite(hvol) .or. hvol <= 0.0d0 .or. nband <= 0 .or. &
        size(target_states,2) /= nband .or. size(target_points) /= npair .or. &
        size(d_band_contribution,1) /= nband .or. size(d_band_contribution,2) /= nband) then
      message='SAWF fragment-block D_band dimensions or grid volume are invalid'
      return
    end if
    if(any(source_points < 1) .or. any(source_points > size(source_states,1)) .or. &
        any(target_points < 1) .or. any(target_points > size(target_states,1))) then
      message='SAWF fragment-block D_band point index is out of range'
      return
    end if
    if(.not.all(ieee_is_finite(real(source_states,8))) .or. &
        .not.all(ieee_is_finite(aimag(source_states))) .or. &
        .not.all(ieee_is_finite(real(target_states,8))) .or. &
        .not.all(ieee_is_finite(aimag(target_states)))) then
      message='SAWF fragment-block D_band states contain a non-finite value'
      return
    end if
    allocate(source_selected(npair,nband),target_selected(npair,nband),stat=allocation_status)
    if(allocation_status /= 0) then
      message='SAWF fragment-block D_band work allocation failed'
      return
    end if
    if(npair > 0) then
      source_selected=source_states(source_points,:)
      target_selected=target_states(target_points,:)
      ! Active source->target map: this is directly M(m,n)=<psi_m|g psi_n>.
      ! It is the forward representation used by sitesym U(Rk)=d U D^dagger.
      d_band_contribution=hvol*matmul(conjg(transpose(target_selected)),source_selected)
    end if
    deallocate(source_selected,target_selected)
    ok=.true.
  end subroutine accumulate_sawf_dmn_band_blocks


  subroutine diagnose_sawf_fragment_basis_closure(source_basis,target_basis, &
      source_points,target_points,hvol,leakage,residual_norm2,transformed_norm2,ok,message)
    real(8), intent(in) :: source_basis(:,:),target_basis(:,:),hvol
    integer, intent(in) :: source_points(:),target_points(:)
    real(8), intent(out) :: leakage,residual_norm2,transformed_norm2
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8), allocatable :: transformed_basis(:,:),projection_coefficient(:,:)
    integer :: allocation_status,p

    ok=.false.; message=''; leakage=0.0d0; residual_norm2=0.0d0; transformed_norm2=0.0d0
    if(size(source_points)/=size(target_points) .or. size(source_points)==0 .or. &
        hvol<=0.0d0 .or. size(source_basis,2)==0 .or. size(target_basis,2)==0) then
      message='invalid SAWF fragment basis closure diagnostic dimensions'
      return
    end if
    if(any(source_points<1) .or. any(source_points>size(source_basis,1)) .or. &
        any(target_points<1) .or. any(target_points>size(target_basis,1))) then
      message='SAWF fragment basis closure diagnostic index is out of range'
      return
    end if
    allocate(transformed_basis(size(target_basis,1),size(source_basis,2)), &
      projection_coefficient(size(target_basis,2),size(source_basis,2)), &
      stat=allocation_status)
    if(allocation_status/=0) then
      message='SAWF fragment basis closure diagnostic allocation failed'
      return
    end if
    transformed_basis=0.0d0
    do p=1,size(source_points)
      transformed_basis(target_points(p),:)=source_basis(source_points(p),:)
    end do
    transformed_norm2=hvol*sum(transformed_basis*transformed_basis)
    if(transformed_norm2<=tiny(1.0d0)) then
      message='SAWF transformed fragment basis block has zero norm'
      deallocate(transformed_basis,projection_coefficient)
      return
    end if
    projection_coefficient=hvol*matmul(transpose(target_basis),transformed_basis)
    transformed_basis=transformed_basis-matmul(target_basis,projection_coefficient)
    residual_norm2=hvol*sum(transformed_basis*transformed_basis)
    leakage=max(0.0d0,residual_norm2/transformed_norm2)
    deallocate(transformed_basis,projection_coefficient)
    ok=ieee_is_finite(leakage) .and. ieee_is_finite(residual_norm2) .and. &
      ieee_is_finite(transformed_norm2)
    if(.not.ok) message='SAWF fragment basis closure diagnostic is non-finite'
  end subroutine diagnose_sawf_fragment_basis_closure


  subroutine validate_sawf_dmn_band(d_band, closure_tolerance, singular_min, singular_max, &
      closure_residual, ok, message)
    complex(8), intent(in) :: d_band(:,:)
    real(8), intent(in) :: closure_tolerance
    real(8), intent(out) :: singular_min, singular_max, closure_residual
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: matrix_work(:,:), work(:), metric(:,:), dummy_u(:,:), dummy_vt(:,:)
    complex(8) :: query(1)
    real(8), allocatable :: singular_values(:), rwork(:)
    integer :: allocation_status, info, i, lwork, n
    external :: zgesvd

    ok=.false.; message=''; singular_min=0.0d0; singular_max=0.0d0
    closure_residual=huge(0.0d0)
    if(size(d_band,1) <= 0 .or. size(d_band,2) /= size(d_band,1) .or. &
        .not.ieee_is_finite(closure_tolerance) .or. closure_tolerance <= 0.0d0) then
      message='SAWF D_band closure validation dimensions or tolerance are invalid'
      return
    end if
    if(.not.all(ieee_is_finite(real(d_band,8))) .or. .not.all(ieee_is_finite(aimag(d_band)))) then
      message='SAWF D_band contains a non-finite value'
      return
    end if
    n=size(d_band,1)
    allocate(matrix_work(n,n),metric(n,n),singular_values(n),rwork(max(1,5*n)), &
      dummy_u(1,1),dummy_vt(1,1),stat=allocation_status)
    if(allocation_status /= 0) then
      message='SAWF D_band validation allocation failed'
      return
    end if
    metric=matmul(conjg(transpose(d_band)),d_band)
    do i=1,n
      metric(i,i)=metric(i,i)-(1.0d0,0.0d0)
    end do
    closure_residual=maxval(abs(metric))
    matrix_work=d_band
    lwork=-1
    call zgesvd('N','N',n,n,matrix_work,n,singular_values,dummy_u,1,dummy_vt,1, &
      query,lwork,rwork,info)
    if(info /= 0 .or. .not.ieee_is_finite(real(query(1),8))) then
      message='SAWF D_band singular-value workspace query failed'
      deallocate(matrix_work,metric,singular_values,rwork,dummy_u,dummy_vt)
      return
    end if
    lwork=max(1,int(real(query(1),8)))
    allocate(work(lwork),stat=allocation_status)
    if(allocation_status /= 0) then
      message='SAWF D_band SVD workspace allocation failed'
      deallocate(matrix_work,metric,singular_values,rwork,dummy_u,dummy_vt)
      return
    end if
    matrix_work=d_band
    call zgesvd('N','N',n,n,matrix_work,n,singular_values,dummy_u,1,dummy_vt,1, &
      work,lwork,rwork,info)
    if(info /= 0 .or. .not.all(ieee_is_finite(singular_values))) then
      message='SAWF D_band singular-value decomposition failed'
      deallocate(matrix_work,metric,singular_values,rwork,dummy_u,dummy_vt,work)
      return
    end if
    singular_min=minval(singular_values)
    singular_max=maxval(singular_values)
    if(max(abs(singular_min-1.0d0),abs(singular_max-1.0d0)) > closure_tolerance .or. &
        closure_residual > 2.0d0*closure_tolerance) then
      write(message,'(a,3(a,es13.5))') 'SAWF selected band subspace is not symmetry-closed:', &
        ' singular_min=',singular_min,' singular_max=',singular_max, &
        ' closure_residual=',closure_residual
      deallocate(matrix_work,metric,singular_values,rwork,dummy_u,dummy_vt,work)
      return
    end if
    deallocate(matrix_work,metric,singular_values,rwork,dummy_u,dummy_vt,work)
    ok=.true.
  end subroutine validate_sawf_dmn_band


  subroutine validate_sawf_representation_product(d_g2, d_g1, d_g2g1, tolerance, &
      residual, ok, message)
    complex(8), intent(in) :: d_g2(:,:), d_g1(:,:), d_g2g1(:,:)
    real(8), intent(in) :: tolerance
    real(8), intent(out) :: residual
    logical, intent(out) :: ok
    character(*), intent(out) :: message

    ok=.false.; message=''; residual=huge(0.0d0)
    if(size(d_g2,1) <= 0 .or. size(d_g2,2) /= size(d_g2,1) .or. &
        any(shape(d_g1) /= shape(d_g2)) .or. any(shape(d_g2g1) /= shape(d_g2)) .or. &
        .not.ieee_is_finite(tolerance) .or. tolerance <= 0.0d0) then
      message='SAWF representation product dimensions or tolerance are invalid'
      return
    end if
    residual=maxval(abs(d_g2g1-matmul(d_g2,d_g1)))
    if(.not.ieee_is_finite(residual) .or. residual > tolerance) then
      write(message,'(a,es13.5)') 'SAWF D_band operation product residual=',residual
      return
    end if
    ok=.true.
  end subroutine validate_sawf_representation_product


  subroutine validate_sawf_operation_set_products(d_set, left_index, right_index, &
      product_index, tolerance, max_residual, ok, message)
    complex(8), intent(in) :: d_set(:,:,:)
    integer, intent(in) :: left_index(:), right_index(:), product_index(:)
    real(8), intent(in) :: tolerance
    real(8), intent(out) :: max_residual
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8) :: residual
    integer :: irelation, noperation
    logical :: relation_ok
    character(256) :: relation_message

    ok=.false.; message=''; max_residual=0.0d0
    noperation=size(d_set,3)
    if(size(d_set,1) <= 0 .or. size(d_set,2) /= size(d_set,1) .or. &
        size(right_index) /= size(left_index) .or. &
        size(product_index) /= size(left_index) .or. size(left_index) <= 0) then
      message='SAWF operation-set representation dimensions are invalid'
      return
    end if
    if(any(left_index < 1) .or. any(left_index > noperation) .or. &
        any(right_index < 1) .or. any(right_index > noperation) .or. &
        any(product_index < 1) .or. any(product_index > noperation)) then
      message='SAWF operation-set product index is out of range'
      return
    end if
    do irelation=1,size(left_index)
      call validate_sawf_representation_product(d_set(:,:,left_index(irelation)), &
        d_set(:,:,right_index(irelation)),d_set(:,:,product_index(irelation)), &
        tolerance,residual,relation_ok,relation_message)
      max_residual=max(max_residual,residual)
      if(.not.relation_ok) then
        write(message,'(a,i0,2a)') 'SAWF operation-set product relation ',irelation, &
          ' failed: ',trim(relation_message)
        return
      end if
    end do
    ok=.true.
  end subroutine validate_sawf_operation_set_products

end module lcfo_wannier_sawf_band
