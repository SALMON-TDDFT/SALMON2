! Synthetic ownership checks; no SCF, orbitals or GS files are generated.
program test_scf_progress
  use structures, only: s_mixing,s_scf_iteration_progress
  implicit none
  type(s_mixing) :: first,second,saved
  integer :: i
  first%flag_mix_zero=.false.;first%num_rho_stock=1
  first%mixrate=.3d0;first%alpha_mb=.4d0;first%beta_p=.5d0
  first%convergence_value_prev=.6d0
  allocate(first%rho_in(2))
  allocate(first%rho_in(1)%f(-1:1,2:3,0:1))
  first%rho_in(1)%f=17d0
  allocate(first%scf_progress%rho_old%f(-1:1,2:3,0:1))
  allocate(first%scf_progress%Vlocal_old%f(-1:1,2:3,0:1))
  first%scf_progress%rho_old%f=2d0
  first%scf_progress%Vlocal_old%f=-3d0
  first%scf_progress%electron_normalization=4d0
  first%scf_progress%iterations=7
  first%scf_progress%residual=.125d0
  first%scf_progress%converged=.true.
  first%scf_progress%no_mix_converged_count=3
  first%scf_progress%rion_update=.true.
  saved=first
  second=first
  second%scf_progress%rho_old%f=99d0
  second%scf_progress%Vlocal_old%f=88d0
  second%scf_progress%iterations=100
  if(any(first%scf_progress%rho_old%f/=2d0))stop 1
  if(any(first%scf_progress%Vlocal_old%f/=-3d0))stop 2
  if(first%scf_progress%iterations/=7)stop 3
  if(any(lbound(saved%scf_progress%rho_old%f)/=[-1,2,0]))stop 4
  if(any(ubound(saved%scf_progress%rho_old%f)/=[1,3,1]))stop 5
  do i=1,3
    second%scf_progress=s_scf_iteration_progress()
    if(second%mixrate/=.3d0)stop 13
    if(any(second%rho_in(1)%f/=17d0))stop 14
    if(allocated(second%scf_progress%rho_old%f))stop 6
    if(allocated(second%scf_progress%Vlocal_old%f))stop 7
    if(second%scf_progress%iterations/=0.or.second%scf_progress%converged)stop 8
    second=saved
    if(second%scf_progress%residual/=.125d0)stop 9
    if(second%scf_progress%electron_normalization/=4d0)stop 10
    if(second%scf_progress%no_mix_converged_count/=3)stop 11
    if(.not.second%scf_progress%rion_update)stop 12
  enddo
  print *, 'SCF progress ownership PASS'
end program
