! Explicit file scopes never change the host's working directory.
program test_directory_context
 use iso_c_binding,only:c_null_char
 use filesystem
 use parallelization,only:setup_parallel,end_parallel,nproc_id_global
 implicit none
 character(4096)::root,a,b,host,filename
 integer::code,unit,rank,ios
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
 call end_parallel
 print *, 'directory scope/host/absolute/relative/failure/sibling PASS'
end program
