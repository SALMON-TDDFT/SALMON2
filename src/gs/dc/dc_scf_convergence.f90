module dc_scf_convergence
  use mpi
  use,intrinsic::iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::reduce_dc_density_convergence
contains
  subroutine reduce_dc_density_convergence(comm,mode,local_absolute_sum,local_square_sum,&
      cell_volume,electron_count,global_point_count,value,ok,message)
    integer,intent(in)::comm,global_point_count
    character(*),intent(in)::mode
    real(real64),intent(in)::local_absolute_sum,local_square_sum,cell_volume,electron_count
    real(real64),intent(out)::value
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::mode_code,integer_controls(2),minimum_integers(2),maximum_integers(2)
    integer::local_bad,global_bad,ierr
    real(real64)::real_controls(2),minimum_reals(2),maximum_reals(2)
    real(real64)::local_sums(2),global_sums(2)

    ok=.false.;message='';value=huge(1d0)
    select case(trim(adjustl(mode)))
    case('rho_dne');mode_code=1
    case('norm_rho');mode_code=2
    case('norm_rho_dng');mode_code=3
    case default;mode_code=0
    end select
    local_bad=0
    if(mode_code<=0.or.global_point_count<=0)local_bad=1
    if(.not.ieee_is_finite(local_absolute_sum).or.&
      .not.ieee_is_finite(local_square_sum).or.&
      .not.ieee_is_finite(cell_volume).or.&
      .not.ieee_is_finite(electron_count))then
      local_bad=1
    else if(local_absolute_sum<0d0.or.local_square_sum<0d0.or.&
      cell_volume<=0d0.or.electron_count<=0d0)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DC density convergence validation reduction failed';return;endif
    if(global_bad/=0)then;message='invalid DC density convergence controls or accumulators';return;endif

    integer_controls=[mode_code,global_point_count]
    call MPI_Allreduce(integer_controls,minimum_integers,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DC density convergence integer agreement failed';return;endif
    call MPI_Allreduce(integer_controls,maximum_integers,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DC density convergence integer agreement failed';return;endif
    real_controls=[cell_volume,electron_count]
    call MPI_Allreduce(real_controls,minimum_reals,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DC density convergence real agreement failed';return;endif
    call MPI_Allreduce(real_controls,maximum_reals,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='DC density convergence real agreement failed';return;endif
    if(any(minimum_integers/=maximum_integers).or.any(minimum_reals/=maximum_reals))then
      message='rank-disagreeing DC density convergence controls';return
    endif

    local_sums=[local_absolute_sum,local_square_sum]
    call MPI_Allreduce(local_sums,global_sums,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      message='DC density convergence sum reduction failed';return
    endif
    if(any(.not.ieee_is_finite(global_sums)))then
      message='DC density convergence sum reduction failed';return
    endif
    select case(mode_code)
    case(1);value=global_sums(1)*cell_volume/electron_count
    case(2);value=global_sums(2)
    case(3);value=global_sums(2)/real(global_point_count,real64)
    end select
    if(.not.ieee_is_finite(value))then
      message='invalid reduced DC density convergence value';return
    endif
    if(value<0d0)then
      message='invalid reduced DC density convergence value';return
    endif
    ok=.true.;message=''
  end subroutine reduce_dc_density_convergence
end module dc_scf_convergence
