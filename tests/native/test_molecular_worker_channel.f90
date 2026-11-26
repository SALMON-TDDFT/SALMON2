program test_worker_channel
 use mpi
 use communication,only:comm_attach,comm_detach
 use parallelization,only:setup_parallel,nproc_id_global,nproc_size_global
 use tddft_worker,only:worker_read_line,worker_state_reply
 implicit none
 integer::ierr,rank,size_world,group,group_size,color,u,ios,total,one
 character(128)::line,expected
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,size_world,ierr)
 if(size_world/=2.and.size_world/=4)call MPI_Abort(MPI_COMM_WORLD,1,ierr)
 group_size=size_world/2;color=rank/group_size
 call MPI_Comm_split(MPI_COMM_WORLD,color,rank,group,ierr)
 call comm_attach(group,ierr)
 if(ierr/=MPI_SUCCESS)call MPI_Abort(MPI_COMM_WORLD,2,ierr)
 call setup_parallel
 if(nproc_size_global/=group_size)call MPI_Abort(MPI_COMM_WORLD,3,ierr)
 u=-1
 write(expected,'(a,i0)')'owned-group-',color
 if(nproc_id_global==0)then
  open(newunit=u,status='scratch',action='readwrite')
  write(u,'(a)')trim(expected)
  rewind(u)
 endif
 call worker_read_line(line,ios,u)
 if(ios/=0.or.trim(line)/=trim(expected))call MPI_Abort(MPI_COMM_WORLD,4,ierr)
 call worker_read_line(line,ios,u)
 if(ios>=0.or.len_trim(line)/=0)call MPI_Abort(MPI_COMM_WORLD,5,ierr)
 call worker_state_reply(5,color,3,4)
 if(nproc_id_global==0)close(u)
 call comm_detach(ierr)
 call MPI_Comm_free(group,ierr)
 one=1
 call MPI_Allreduce(one,total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
 if(total/=size_world)call MPI_Abort(MPI_COMM_WORLD,6,ierr)
 if(rank==0)print *,'molecular subgroup control and parent survival PASS'
 call MPI_Finalize(ierr)
end program
