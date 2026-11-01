program probe
 use mpi
 use parallelization
 use math_constants
 use timer
 use perflog
 use salmon_global,only:theory,iperiodic,yn_dc_lcfo_diag
 implicit none
 integer::ierr,fd
 call MPI_Init(ierr)
 nproc_group_global=MPI_COMM_WORLD
 call MPI_Comm_rank(MPI_COMM_WORLD,nproc_id_global,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,nproc_size_global,ierr)
 call set_math_constants()
 call timer_initialize()
 theory='dft';iperiodic=3;yn_dc_lcfo_diag='n'
 call timer_set(LOG_TOTAL,10d0)
 call timer_set(LOG_INIT_GS,1d0+real(nproc_id_global,8)*1d-9)
 call timer_set(LOG_INIT_GS_ITERATION,merge(0d0,1d-6,nproc_id_global==0))
 call timer_set(LOG_GS_ITERATION,2d0)
 call timer_set(LOG_WRITE_GS_RESULTS,real(nproc_id_global+1,8))
 call timer_set(LOG_WRITE_GS_DATA,merge(1d-320,1d-6,nproc_id_global==0))
 open(newunit=fd,file='readable.rank'//rank_string()//'.txt',status='replace')
 call write_performance(fd,write_mode_readable)
 close(fd)
 open(newunit=fd,file='csv.rank'//rank_string()//'.txt',status='replace')
 call write_performance(fd,write_mode_csv)
 close(fd)
 call MPI_Finalize(ierr)
contains
 function rank_string() result(s)
 implicit none
 character(8)::s
 write(s,'(i8.8)')nproc_id_global
 end function
end program
