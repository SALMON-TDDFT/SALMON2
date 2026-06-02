module dc_fragment_geometry
  use structures, only: s_dcdft
  implicit none

  private
  public :: get_fragment_domain
  public :: optimize_fragment_geometry

contains

  subroutine get_fragment_domain(dc, ifrag, nxyz_domain)
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ifrag
    integer, intent(out) :: nxyz_domain(3)

    if (dc%optimized_fragment_geometry .and. allocated(dc%nxyz_domain_frag)) then
      nxyz_domain(1:3) = dc%nxyz_domain_frag(1:3, ifrag)
    else
      nxyz_domain(1:3) = dc%nxyz_domain(1:3)
    end if
  end subroutine get_fragment_domain

  subroutine optimize_fragment_geometry(dc, num_fragment, num_rgrid, al, natom, rion)
    type(s_dcdft), intent(inout) :: dc
    integer, intent(in) :: num_fragment(3), num_rgrid(3), natom
    real(8), intent(in) :: al(3), rion(3, natom)
    integer :: axis, iboundary, delta, nfrag_axis_max
    integer, allocatable :: widths(:,:), candidate_widths(:,:)
    integer :: widths_out(3, product(num_fragment))
    real(8) :: objective, candidate_objective, initial_objective
    real(8) :: load_frag(product(num_fragment))
    logical :: valid, improved

    nfrag_axis_max = maxval(num_fragment)
    allocate(widths(3, nfrag_axis_max), candidate_widths(3, nfrag_axis_max))
    widths = 0
    candidate_widths = 0

    do axis = 1, 3
      call initialize_uniform_widths(num_rgrid(axis), num_fragment(axis), widths(axis, 1:num_fragment(axis)))
    end do

    call evaluate_fragment_partition(widths, num_fragment, num_rgrid, al, natom, rion, objective, valid, load_frag)
    if (.not. valid) then
      stop "DC fragment optimization: uniform initial partition contains vacuum fragment(s)"
    end if
    initial_objective = objective
    write(*,'(1x,a)') "DC fragment optimization enabled"
    write(*,'(1x,a,es12.5)') "initial objective=", initial_objective

    improved = .true.
    do while (improved)
      improved = .false.
      do axis = 1, 3
        do iboundary = 1, num_fragment(axis) - 1
          do delta = -1, 1, 2
            candidate_widths(:, :) = widths(:, :)
            candidate_widths(axis, iboundary) = candidate_widths(axis, iboundary) + delta
            candidate_widths(axis, iboundary + 1) = candidate_widths(axis, iboundary + 1) - delta
            if (minval(candidate_widths(axis, 1:num_fragment(axis))) < 1) cycle
            call evaluate_fragment_partition(candidate_widths, num_fragment, num_rgrid, al, natom, rion, &
                                             candidate_objective, valid)
            if (valid .and. candidate_objective < objective) then
              widths(:, :) = candidate_widths(:, :)
              objective = candidate_objective
              improved = .true.
            end if
          end do
        end do
      end do
    end do

    call flatten_axis_widths(widths, num_fragment, widths_out)
    if (.not. allocated(dc%nxyz_domain_frag)) allocate(dc%nxyz_domain_frag(3, product(num_fragment)))
    dc%nxyz_domain_frag(:, :) = widths_out(:, :)
    dc%optimized_fragment_geometry = .true.
    call evaluate_fragment_partition(widths, num_fragment, num_rgrid, al, natom, rion, objective, valid, load_frag)
    call print_fragment_optimization_report(widths, num_fragment, initial_objective, objective, load_frag)

    deallocate(candidate_widths, widths)
  end subroutine optimize_fragment_geometry

  subroutine initialize_uniform_widths(num_grid, nfrag_axis, widths_axis)
    integer, intent(in) :: num_grid, nfrag_axis
    integer, intent(out) :: widths_axis(nfrag_axis)
    integer :: i

    widths_axis(:) = num_grid / nfrag_axis
    do i = 1, mod(num_grid, nfrag_axis)
      widths_axis(i) = widths_axis(i) + 1
    end do
  end subroutine initialize_uniform_widths

  subroutine flatten_axis_widths(widths, num_fragment, widths_out)
    integer, intent(in) :: widths(:, :), num_fragment(3)
    integer, intent(out) :: widths_out(3, product(num_fragment))
    integer :: ifrag, ix, iy, iz

    ifrag = 0
    do ix = 1, num_fragment(1)
      do iy = 1, num_fragment(2)
        do iz = 1, num_fragment(3)
          ifrag = ifrag + 1
          widths_out(1, ifrag) = widths(1, ix)
          widths_out(2, ifrag) = widths(2, iy)
          widths_out(3, ifrag) = widths(3, iz)
        end do
      end do
    end do
  end subroutine flatten_axis_widths

  subroutine evaluate_fragment_partition(widths, num_fragment, num_rgrid, al, natom, rion, objective, valid, load_frag_out)
    integer, intent(in) :: widths(:, :), num_fragment(3), num_rgrid(3), natom
    real(8), intent(in) :: al(3), rion(3, natom)
    real(8), intent(out) :: objective
    logical, intent(out) :: valid
    real(8), intent(out), optional :: load_frag_out(product(num_fragment))
    integer :: starts(3, maxval(num_fragment))
    real(8) :: bounds(3, maxval(num_fragment) + 1)
    real(8) :: load_frag(product(num_fragment)), target_load
    integer :: n_frag, axis, iseg

    n_frag = product(num_fragment)
    load_frag(:) = 0.0d0

    do axis = 1, 3
      if (sum(widths(axis, 1:num_fragment(axis))) /= num_rgrid(axis)) then
        valid = .false.
        objective = huge(0.0d0)
        return
      end if
      starts(axis, 1) = 0
      bounds(axis, 1) = 0.0d0
      do iseg = 1, num_fragment(axis)
        if (iseg > 1) starts(axis, iseg) = starts(axis, iseg - 1) + widths(axis, iseg - 1)
        bounds(axis, iseg + 1) = al(axis) * dble(sum(widths(axis, 1:iseg))) / dble(num_rgrid(axis))
      end do
    end do

    call accumulate_fragment_loads(bounds, num_fragment, natom, rion, load_frag)

    if (any(load_frag <= 0.0d0)) then
      valid = .false.
      objective = huge(0.0d0)
      return
    end if

    target_load = sum(load_frag(1:n_frag)) / dble(n_frag)
    objective = sum((load_frag(1:n_frag) - target_load)**2)
    if (present(load_frag_out)) load_frag_out(1:n_frag) = load_frag(1:n_frag)
    valid = .true.
  end subroutine evaluate_fragment_partition

  subroutine accumulate_fragment_loads(bounds, num_fragment, natom, rion, load_frag)
    real(8), intent(in) :: bounds(:, :), rion(3, natom)
    integer, intent(in) :: num_fragment(3), natom
    real(8), intent(inout) :: load_frag(product(num_fragment))
    integer :: iatom, ix, iy, iz, ifrag

    do iatom = 1, natom
      ix = find_axis_segment(rion(1, iatom), bounds(1, 1:num_fragment(1)+1), num_fragment(1))
      iy = find_axis_segment(rion(2, iatom), bounds(2, 1:num_fragment(2)+1), num_fragment(2))
      iz = find_axis_segment(rion(3, iatom), bounds(3, 1:num_fragment(3)+1), num_fragment(3))
      if (ix < 1 .or. iy < 1 .or. iz < 1) cycle
      ifrag = ((ix - 1) * num_fragment(2) + (iy - 1)) * num_fragment(3) + iz
      load_frag(ifrag) = load_frag(ifrag) + 1.0d0
    end do
  end subroutine accumulate_fragment_loads

  integer function find_axis_segment(coord, bounds, nseg) result(iseg)
    real(8), intent(in) :: coord
    real(8), intent(in) :: bounds(:)
    integer, intent(in) :: nseg
    integer :: i

    iseg = -1
    do i = 1, nseg
      if (coord >= bounds(i) .and. (coord < bounds(i + 1) .or. (i == nseg .and. coord <= bounds(i + 1)))) then
        iseg = i
        return
      end if
    end do
  end function find_axis_segment

  subroutine print_fragment_optimization_report(widths, num_fragment, initial_objective, final_objective, load_frag)
    integer, intent(in) :: widths(:, :), num_fragment(3)
    real(8), intent(in) :: initial_objective, final_objective, load_frag(product(num_fragment))
    integer :: axis

    write(*,'(1x,a,es12.5)') "final objective=", final_objective
    write(*,'(1x,a)') "axis widths:"
    do axis = 1, 3
      write(*,'(1x,a,i0,a,100(1x,i0))') "axis ", axis, ":", widths(axis, 1:num_fragment(axis))
    end do
    write(*,'(1x,a,1x,es12.5,1x,es12.5,1x,es12.5)') "fragment load stats:", &
      minval(load_frag), maxval(load_frag), sum(load_frag)/dble(size(load_frag))
  end subroutine print_fragment_optimization_report

end module dc_fragment_geometry
