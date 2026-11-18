! Copyright 2026 SALMONpy contributors. Apache-2.0.
! Run in a copied standard RT fixture with host MPI; no new solver parameters.
program output_close
 use mpi,only:MPI_Init,MPI_Finalize,MPI_Comm_rank,MPI_Comm_size,MPI_COMM_WORLD
 use iso_c_binding
 use communication,only:comm_attach,comm_detach
 use filesystem,only:open_filehandle
 use salmonpy_tddft_context,only:sp_tddft_init,sp_tddft_close
 implicit none
 interface
   integer(c_int) function transaction(mode) bind(C,name='sp_tddft_step_transaction')
     import c_int
     integer(c_int),value::mode
   end function
 end interface
 real(c_double)::a(3)=0d0,j(3)
 integer(c_int)::decomposition(5)
 integer::ierr,rank,n,status,foreign
 logical::connected
 character(32)::mode
 call get_command_argument(1,mode)
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,n,ierr)
 call comm_attach(MPI_COMM_WORLD,ierr);if(ierr/=0)stop 1
 decomposition=[1,1,1,1,n]
 status=sp_tddft_init('rt.inp'//c_null_char,.001d0,2,a,decomposition,j)
 if(status/=0)stop 2
 if(mode=='protected')then
   status=transaction(0);if(status/=0)stop 3
   status=transaction(1);if(status/=0)stop 4
   if(rank==0)foreign=open_filehandle('foreign.data','new')
 endif
 status=sp_tddft_close();if(status/=0)stop 5
 if(rank==0)then
   inquire(file='rt_rt.data',opened=connected);if(connected)stop 6
   inquire(file='rt_rt_energy.data',opened=connected);if(connected)stop 7
   if(mode=='protected')then
     inquire(unit=foreign,opened=connected);if(.not.connected)stop 8
     write(foreign,*)'Host still owns its output'
     close(foreign)
   endif
 endif
 call comm_detach(ierr);call MPI_Finalize(ierr)
 print *,'TDDFT owned output close PASS ',trim(mode)
end program
