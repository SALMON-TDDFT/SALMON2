! Compare actual extracted before/after routines. No electron solver is called.
program test_scf_convergence
  use iso_fortran_env, only: int64
  use mpi
  use baseline_convergence, only: original=>measure
  use owned_convergence, only: owned=>measure
  implicit none
  real(8) :: rho(3,2,2),potential(3,2,2),next_rho(3,2,2),next_potential(3,2,2)
  real(8) :: expected,actual,delta(3)
  integer :: ierr,rank,i,j,spin,spin_electrons(2),spin_case
  character(12) :: mode(5)
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  mode=[character(12) :: 'rho_dne','norm_rho','norm_rho_dng','norm_pot','norm_pot_dng']
  delta=[0d0,1d-12,.125d0]
  rho=.25d0*(rank+1);potential=-.5d0*(rank+1)
  do spin=1,2
    do spin_case=1,2
      spin_electrons=0
      if(spin_case==2)spin_electrons=[2,4]
      do i=1,size(mode)
        do j=1,size(delta)
          next_rho=rho+delta(j);next_potential=potential-delta(j)
          call original(rho,next_rho,potential,next_potential,trim(mode(i)),spin,4,spin_electrons,expected)
          call owned(rho,next_rho,potential,next_potential,trim(mode(i)),spin,4,spin_electrons,actual)
          if(transfer(actual,0_int64)/=transfer(expected,0_int64))stop 91
        enddo
      enddo
    enddo
  enddo
  if(rank==0)print *, 'SCF convergence parity PASS'
  call MPI_Finalize(ierr)
end program
