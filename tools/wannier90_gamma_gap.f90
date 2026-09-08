program wannier90_gamma_gap
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: ev_to_au = 1.0_dp / 27.211386245988_dp

  character(len=512) :: hr_file
  character(len=512) :: arg
  integer :: nocc
  integer :: nwann, nrpts, ir, i, j, info, lwork
  integer :: rx, ry, rz
  real(dp) :: hre, him, gap_ev
  real(dp), allocatable :: h(:,:), eval(:), work(:)
  character(len=1024) :: header
  integer :: unit

  interface
    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import :: dp
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      real(dp), intent(inout) :: a(lda,*), work(*)
      real(dp), intent(out) :: w(*)
      integer, intent(out) :: info
    end subroutine dsyev
  end interface

  if(command_argument_count() < 2) then
    write(*,'(a)') 'usage: wannier90_gamma_gap SEED_hr.dat NOCC'
    stop 2
  end if

  call get_command_argument(1, hr_file)
  call get_command_argument(2, arg)
  read(arg,*) nocc

  open(newunit=unit, file=trim(hr_file), status='old', action='read')
  read(unit,'(a)') header
  read(unit,*) nwann
  read(unit,*) nrpts
  allocate(h(nwann,nwann), eval(nwann))
  h = 0.0_dp

  ! Degeneracies occupy enough integer records to cover nrpts values.
  call skip_degeneracy_records(unit, nrpts)

  do ir=1,nrpts*nwann*nwann
    read(unit,*,end=100) rx, ry, rz, i, j, hre, him
    if(rx /= 0 .or. ry /= 0 .or. rz /= 0) cycle
    if(i < 1 .or. i > nwann .or. j < 1 .or. j > nwann) cycle
    h(i,j) = h(i,j) + hre * ev_to_au
  end do
100 continue
  close(unit)

  h = 0.5_dp * (h + transpose(h))
  lwork = max(1, 8*nwann)
  allocate(work(lwork))
  call dsyev('V', 'U', nwann, h, nwann, eval, work, lwork, info)
  if(info /= 0) then
    write(*,'(a,i0)') 'dsyev failed: info=', info
    stop 1
  end if
  if(nocc < 1 .or. nocc >= nwann) stop 'nocc must be in 1..nwann-1'

  gap_ev = (eval(nocc+1) - eval(nocc)) * 27.211386245988_dp
  write(*,'(a)') '# Wannier90 gamma HR gap diagnostic'
  write(*,'(a,1x,a)') 'hr_file =', trim(hr_file)
  write(*,'(a,1x,i0)') 'nwann =', nwann
  write(*,'(a,1x,i0)') 'nocc =', nocc
  write(*,'(a,1x,es16.8)') 'homo(eV) =', eval(nocc) * 27.211386245988_dp
  write(*,'(a,1x,es16.8)') 'lumo(eV) =', eval(nocc+1) * 27.211386245988_dp
  write(*,'(a,1x,es16.8)') 'gap(eV) =', gap_ev

contains

  subroutine skip_degeneracy_records(unit, nrpts)
    integer, intent(in) :: unit, nrpts
    integer :: skipped, ios
    character(len=1024) :: line

    skipped = 0
    do while(skipped < nrpts)
      read(unit,'(a)',iostat=ios) line
      if(ios /= 0) stop 'failed to read Wannier90 HR degeneracy block'
      skipped = skipped + count_tokens(line)
    end do
  end subroutine skip_degeneracy_records

  integer function count_tokens(line) result(n)
    character(len=*), intent(in) :: line
    integer :: k
    logical :: in_token

    n = 0
    in_token = .false.
    do k=1,len_trim(line)
      if(line(k:k) == ' ' .or. line(k:k) == char(9)) then
        in_token = .false.
      else if(.not. in_token) then
        n = n + 1
        in_token = .true.
      end if
    end do
  end function count_tokens

end program wannier90_gamma_gap
