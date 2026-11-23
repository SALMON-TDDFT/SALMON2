! Explicit file scopes never change the host's working directory.
program test_directory_context
 use iso_c_binding,only:c_null_char
 use filesystem
 use parallelization,only:setup_parallel,end_parallel,nproc_id_global
 use communication,only:comm_attach,comm_detach,comm_begin_children,comm_end_children, &
   comm_resume_children,comm_release_children,comm_create_group,comm_child_counts
 use mpi,only:MPI_COMM_WORLD,MPI_Barrier,MPI_SUCCESS
 implicit none
 character(4096)::root,a,b,host,filename
 integer::code,unit,rank,ios,owner_a,owner_b,child,live,slots
 logical::exists
 call setup_parallel
 call get_command_argument(1,root)
 rank=nproc_id_global
 write(a,'(a,a,i0)')trim(root),'/owner-a-',rank
 write(b,'(a,a,i0)')trim(root),'/owner-b-',rank
 write(host,'(a,i0)')'host-',rank
 call create_directory(trim(a))
 call create_directory(trim(b))
 open(newunit=unit,file=trim(host),status='replace')
 write(unit,'(a)')'host';close(unit)
 if(sp_io_directory_protocol()<1.or.sp_io_directory_protocol()>3)stop 1
 code=sp_io_directory_enter(trim(root)//'/missing'//c_null_char)
 if(code==0)stop 2
 code=sp_io_directory_enter(trim(a)//c_null_char)
 if(code/=0)stop 3
 code=sp_io_directory_enter(trim(b)//c_null_char)
 if(code==0)stop 4
 inquire(file=trim(host),exist=exists)
 if(.not.exists)stop 5
 if(file_exists(trim(host)))stop 6
 call create_directory('child')
 filename=resolve_path('child/value')
 open(newunit=unit,file=trim(filename),status='replace')
 write(unit,'(a)')'a';close(unit)
 if(.not.file_exists('child/value'))stop 7
 if(resolve_path(trim(b)) /= trim(b))stop 8
 call sp_io_directory_leave()
 if(resolve_path('relative')/='relative')stop 9
 code=sp_io_directory_enter(trim(b)//c_null_char)
 if(code/=0)stop 10
 if(file_exists('child/value'))stop 11
 open(newunit=unit,file=resolve_path('value'),status='replace')
 write(unit,'(a)')'b';close(unit)
 open(newunit=unit,file=resolve_path('missing'),status='old',iostat=ios)
 if(ios==0)stop 12
 call sp_io_directory_leave()
 code=sp_io_directory_enter(trim(a)//c_null_char)
 if(code/=0)stop 13
 open(newunit=unit,file=resolve_path('child/value'),status='old')
 read(unit,'(a)')filename;close(unit)
 if(trim(filename)/='a')stop 14
 call remove_file('child/value')
 if(file_exists('child/value'))stop 15
 call sp_io_directory_leave()
 inquire(file=trim(host),exist=exists)
 if(.not.exists)stop 16
 ! Attach without taking ownership of host MPI, then resume each context's scope.
 call comm_attach(MPI_COMM_WORLD,ios)
 if(ios/=MPI_SUCCESS)stop 20
 call comm_begin_children(owner_a,ios)
 if(ios/=0)stop 21
 child=comm_create_group(MPI_COMM_WORLD,0,rank)
 call comm_end_children(owner_a)
 call comm_begin_children(owner_b,ios)
 if(ios/=0)stop 22
 child=comm_create_group(MPI_COMM_WORLD,0,rank)
 call comm_end_children(owner_b)
 call comm_resume_children(owner_a,ios)
 if(ios/=0)stop 23
 call comm_resume_children(owner_b,ios)
 if(ios==0)stop 24
 child=comm_create_group(MPI_COMM_WORLD,0,rank)
 call comm_end_children(owner_a)
 call comm_resume_children(owner_b,ios)
 if(ios/=0)stop 25
 child=comm_create_group(MPI_COMM_WORLD,0,rank)
 call comm_end_children(owner_b)
 call comm_child_counts(live,slots)
 if(live/=4)stop 26
 call comm_release_children(owner_a)
 call comm_child_counts(live,slots)
 if(live/=2)stop 27
 call comm_resume_children(huge(0),ios)
 if(ios==0)stop 28
 child=comm_create_group(MPI_COMM_WORLD,0,rank)
 call comm_release_children(owner_a)
 call comm_release_children(owner_b)
 call comm_child_counts(live,slots)
 if(live/=1)stop 29
 call comm_detach(ios)
 if(ios/=MPI_SUCCESS)stop 30
 call MPI_Barrier(MPI_COMM_WORLD,ios)
 if(ios/=MPI_SUCCESS)stop 31
 call end_parallel
 print *, 'MPI resumed owners/late children/sibling/host PASS'
 print *, 'directory scope/host/absolute/relative/failure/sibling PASS'
end program
