program conditioning_driver
 use mpi
 use hse_ace
 use exx_orbitals
 implicit none
 type(hse_ace_state) :: dense,sparse
 complex(8) :: u(2,2,1),w(2,2,1),target(2,2,1),a(2,2,1),b(2,2,1)
 integer :: status,err
 call MPI_Init(err)
 u(:,:,1)=reshape([1d0,1.1d-6,1d0,-1.1d-6],[2,2])/sqrt(2d0);w=-u
 target=0d0;target(1,1,1)=1d0;target(2,2,1)=1d0
 call orbital_ace_build(dense,u,w,1d0,MPI_COMM_SELF,MPI_COMM_SELF,status)
 if(status/=0)error stop 'accepted near-limit dense metric rejected'
 call orbital_ace_build(sparse,u,w,1d0,MPI_COMM_SELF,MPI_COMM_SELF,status,packed=.true.)
 if(status/=0)error stop 'accepted near-limit packed metric rejected'
 call orbital_ace_apply(dense,target,a,MPI_COMM_SELF,MPI_COMM_SELF,status)
 if(status/=0)error stop 'dense action'
 call orbital_ace_apply(sparse,target,b,MPI_COMM_SELF,MPI_COMM_SELF,status)
 if(status/=0)error stop 'packed action'
 print *, 'CONDITION / DENSE-SPARSE ERROR ',dense%condition,maxval(abs(a-b))
 if(maxval(abs(a-b))>1d-9)error stop 'packed metric lost dense-action accuracy'
 call MPI_Finalize(err)
end program
