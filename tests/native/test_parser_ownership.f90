program test_parser_ownership
  use inputoutput, only: s_input_parser,read_input,read_stdin,read_input_common, &
    read_atomic_coordinates,check_bad_input,close_parser_workspace,initialize_inputoutput_units
  use salmon_global
  use math_constants, only: set_math_constants
  use parallelization, only: setup_parallel,end_parallel,nproc_id_global
  use communication, only: comm_is_root
  implicit none
  type(s_input_parser) :: a,b,again,left,right,blank
  character(32) :: mode
  integer :: status
  call set_math_constants
  call setup_parallel
  call get_command_argument(1,mode)
  if(mode=='legacy-a'.or.mode=='legacy-b')then
    if(mode=='legacy-a')call read_input('a.inp')
    if(mode=='legacy-b')call read_input('b.inp')
    if(comm_is_root(nproc_id_global))print *,'legacy parse PASS'
    call end_parallel
    stop
  endif
  a%log_path='a.variables.log';b%log_path='b.variables.log';again%log_path='again.variables.log'
  call read_input('a.inp',a)
  if(xc/='libxc_pbe'.or.natom/=2.or.mixrate/=.23d0)error stop 'first parsed input'
  if(any([a%fh_namelist,a%fh_atomic_coor,a%fh_atomic_red_coor,a%fh_reentrance]/=0)) &
    error stop 'first scratch units leaked'
  call read_input('b.inp',b)
  if(xc/='none'.or.natom/=1.or.mixrate/=.5d0)error stop 'omitted defaults inherited'
  if(any(num_kgrid/=1).or.E_amplitude1/=0d0)error stop 'k/pulse input inherited'
  if(allocated(Rion).and.size(Rion,2)/=1)error stop 'atom allocation inherited'
  if(a%inml_functional/=0.or.b%inml_functional>=0)error stop 'diagnostic ownership'
  if(.not.a%if_nml_coor.or.b%if_nml_coor.or.b%if_nml_red_coor)error stop 'coordinate flags inherited'
  if(unit_system/='au'.or.any(abs(Rion_red(:,1)-[.125d0,.25d0,.375d0])>1d-14)) &
    error stop 'external reduced coordinates or units inherited'
  call a%values%activate(status)
  call initialize_inputoutput_units
  if(status/=0.or.xc/='libxc_pbe'.or.natom/=2.or.mixrate/=.23d0)error stop 'first owner changed'
  if(unit_system/='A_eV_fs')error stop 'unit selector ownership'
  call b%values%activate(status)
  if(status/=0.or.xc/='none'.or.natom/=1)error stop 'second owner changed'
  call a%defaults%activate(status)
  if(status/=0.or.natom/=0.or.allocated(Rion))error stop 'default record inherited atoms'
  if(mixrate/=.5d0.or.xc/='none')error stop 'default record contains explicit values'
  call read_input('a.inp',again)
  if(xc/='libxc_pbe'.or.natom/=2.or.mixrate/=.23d0)error stop 'A B A mismatch'
  ! Split both inputs before consuming either: their scratch files cannot collide.
  call read_stdin('a.inp',left)
  call read_stdin('b.inp',right)
  if(comm_is_root(nproc_id_global))then
    if(left%fh_namelist==right%fh_namelist)error stop 'shared scratch connection'
  endif
  call read_input_common(left)
  call read_atomic_coordinates(left)
  call check_bad_input
  call left%values%capture()
  call read_input_common(right)
  call read_atomic_coordinates(right)
  call check_bad_input
  call right%values%capture()
  call left%values%activate(status)
  if(status/=0.or.xc/='libxc_pbe'.or.natom/=2)error stop 'interleaved scratch source'
  call right%values%activate(status)
  if(status/=0.or.xc/='none'.or.natom/=1)error stop 'second scratch source'
  call close_parser_workspace(left)
  call close_parser_workspace(right)
  blank%log_path='blank.variables.log'
  call read_input('blank.inp',blank)
  if(natom/=0.or.blank%if_nml_coor.or.blank%if_nml_red_coor)error stop 'absent atom source inherited'
  if(size(Rion,2)/=0)error stop 'absent atoms inherited'
  again%log_path='reuse.variables.log'
  call read_input('b.inp',again)
  if(xc/='none'.or.natom/=1.or.mixrate/=.5d0)error stop 'workspace reuse inherited input'
  call read_input('a.inp',again)
  if(xc/='libxc_pbe'.or.natom/=2.or.mixrate/=.23d0)error stop 'workspace reuse A B A'
  call b%values%activate(status)
  if(status/=0.or.xc/='none'.or.natom/=1)error stop 'workspace reuse changed sibling'
  if(comm_is_root(nproc_id_global))print *,'model parser ownership PASS'  
  call end_parallel
end program
