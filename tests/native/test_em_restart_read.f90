! Explicit GS-free MPI test: error status, caller survival and repaired retry.
program em_restart_read_probe
  use,intrinsic :: iso_c_binding,only:c_int,c_char,c_null_char
  use mpi
  use parallelization,only:setup_parallel,end_parallel,nproc_id_global,nproc_size_global,nproc_group_global
  use salmon_global,only:directory_read_data
  use inputoutput,only:read_input,s_input_parser
  use structures,only:s_fdtd_system
  use lorentz_drude,only:s_ld_state
  use fdtd_eh,only:ls_fdtd_eh,eh_init,eh_calc,eh_finalize
  use math_constants,only:set_math_constants
  use timer,only:timer_initialize
  use common_maxwell,only:input_r_txt_em,input_r_bin_em,output_r_txt_em,output_r_bin_em
  implicit none
  type(s_input_parser) :: parser
  type(s_fdtd_system) :: fs
  type(s_ld_state) :: ld
  type(ls_fdtd_eh) :: fe
  integer :: ios=0,lo(5)=1,hi(5)=2,u,err,low,high,value,nold
  real(8) :: x(2)=[2d0,3d0]
  complex(8) :: z(2)=[(2d0,4d0),(3d0,5d0)]
  character(32) :: mode,id,number
  character(256) :: target
  logical :: exists,opened
  interface
    integer(c_int) function rename_file(old,new) bind(C,name='rename')
      import c_int,c_char
      character(kind=c_char),intent(in) :: old(*),new(*)
    end function
  end interface
  call set_math_constants
  call setup_parallel
  directory_read_data='restart'
  call get_command_argument(1,mode)
  if(mode=='batch-full'.or.mode=='batch-half')then
    call read_input('inputfile',parser,ios,native_host=.true.)
    if(ios/=0)error stop 'batch input failed'
    call timer_initialize
    call eh_init(fs,fe,ld)
    ! Fixture-only early finish: physical nt_em and window settings stay fixed.
    if(mode=='batch-half')fe%iter_end=fe%iter_end/2
    call eh_calc(fs,fe,ld)
    call eh_finalize(fs,fe,ld)
    call end_parallel
    stop
  endif
  call get_command_argument(2,target)
  call get_command_argument(3,number)
  nold=nproc_size_global
  if(len_trim(number)>0)read(number,*)nold
  write(id,'(I0)')nproc_id_global
  open(unit=600,file='host-'//trim(id),status='new')
  write(600,'(A)')'before'
  if(mode=='prepare')then
    ! Writers use their fixed legacy path; the harness copies this set to restart/.
    call output_r_txt_em(600,1,1,'single','value',i0d=23,io_status=ios)
    call output_r_bin_em(600,lo,hi,'single','real',r1d=x,io_status=ios)
    call output_r_bin_em(600,lo,hi,'single','complex',c1d=z,io_status=ios)
    call output_r_bin_em(600,lo,hi,'all','real',r1d=x,io_status=ios)
    call output_r_bin_em(600,lo,hi,'all','complex',c1d=z,io_status=ios)
    if(ios/=0)error stop 'fixture save failed'
  else
    x=-999d0;z=(-999d0,17d0);value=-999
    call probe
    call MPI_Allreduce(ios,low,1,MPI_INTEGER,MPI_MIN,nproc_group_global,err)
    call MPI_Allreduce(ios,high,1,MPI_INTEGER,MPI_MAX,nproc_group_global,err)
    if(low/=high)error stop 'reader status differs across ranks'
    if(len_trim(target)>0)then
      if(low==0)error stop 'reader did not report a failure'
      open(newunit=u,file='failure-status-'//trim(id),status='new')
      write(u,*)ios
      close(u)
      if(nproc_id_global==0)then
        inquire(file=trim(target),exist=exists)
        if(exists)then
          err=rename_file(trim(target)//c_null_char,trim(target)//'.failed'//c_null_char)
          if(err/=0)error stop 'failed fixture could not be retained'
        endif
        err=rename_file(trim(target)//'.repair'//c_null_char,trim(target)//c_null_char)
        if(err/=0)error stop 'fixture repair failed'
      endif
      ! Fixture-only rendezvous after repairing shared test files.
      call MPI_Barrier(nproc_group_global,err)
      ios=0
      x=-998d0;z=(-998d0,19d0);value=-998
      call probe
    endif
    if(ios/=0)error stop 'normal/repaired read failed'
    select case(trim(mode))
    case('text')
      if(value/=23)error stop 'wrong scalar after read'
    case('single-real','local-real','remap-real')
      if(any(x/=[2d0,3d0]))error stop 'wrong real payload after read'
    case('single-complex','local-complex','remap-complex')
      if(any(z/=[(2d0,4d0),(3d0,5d0)]))error stop 'wrong complex payload after read'
    end select
  endif
  inquire(unit=600,opened=opened)
  if(.not.opened)error stop 'reader closed caller unit'
  write(600,'(A)')'after'
  close(600)
  call MPI_Barrier(nproc_group_global,err)
  if(nproc_id_global==0)print *,'EM restart read caller survival and retry PASS'
  call end_parallel
contains
  subroutine probe
    select case(trim(mode))
    case('text')
      call input_r_txt_em(600,1,1,'value',i0d=value,io_status=ios)
    case('single-real')
      call input_r_bin_em(600,lo,hi,nold,.true.,'single','real',r1d=x,io_status=ios)
    case('single-complex')
      call input_r_bin_em(600,lo,hi,nold,.true.,'single','complex',c1d=z,io_status=ios)
    case('local-real')
      call input_r_bin_em(600,lo,hi,nold,.true.,'all','real',r1d=x,io_status=ios)
    case('local-complex')
      call input_r_bin_em(600,lo,hi,nold,.true.,'all','complex',c1d=z,io_status=ios)
    case('remap-real')
      call input_r_bin_em(600,lo,hi,nold,.false.,'all','real',r1d=x,io_status=ios)
    case('remap-complex')
      call input_r_bin_em(600,lo,hi,nold,.false.,'all','complex',c1d=z,io_status=ios)
    case default
      error stop 'unknown probe mode'
    end select
  end subroutine probe
end program em_restart_read_probe
