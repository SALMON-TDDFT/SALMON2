program ptcn_spinor_probe
 use mpi
 use exx_ptcn_core,only:exx_ptcn_solve
 use ptcn_native_matrix
 implicit none
 type(ptcn_owner_matrix)::p
 complex(8),allocatable::u(:,:,:),x(:,:,:)
 real(8)::dt,error,gram,diff,global_diff,occup(2),total,local
 integer::s,rank,np,first,last,j,err,builds,apps,mode
 character(10)::arg
 call MPI_Init(s);call MPI_Comm_rank(MPI_COMM_WORLD,rank,s);call MPI_Comm_size(MPI_COMM_WORLD,np,s)
 mode=0;call get_command_argument(1,arg);if(len_trim(arg)>0)read(arg,*)mode
 first=rank*2/np+1;last=(rank+1)*2/np
 call ptcn_owner_create(p,4,2,first,last,MPI_COMM_WORLD,err)
 if(err/=0)call MPI_Abort(MPI_COMM_WORLD,1,s)
 allocate(u(4,last-first+1,1),x(4,last-first+1,1));u=0
 do j=first,last;u(j,j-first+1,1)=1;enddo
 occup=[0.8d0,0.2d0];dt=0.01d0
 call exx_ptcn_solve(u,dt,1d0,action,prec,norm,x,error,builds,apps,err, &
  noccupied=2,use_anderson=.true.,projection=project,consensus=agree,mixing_reduce=reduce_gram)
 if(mode==2)then
  if(err==0)call MPI_Abort(MPI_COMM_WORLD,8,s)
  if(rank==0)print *,'FINAL_RESIDUAL_REJECTION_PASS'
  call ptcn_owner_destroy(p);call MPI_Finalize(s);stop
 endif
 if(err/=0.or.error>1d-10)call MPI_Abort(MPI_COMM_WORLD,2,s)
 call ptcn_owner_gram(p,x(:,:,1),1d0,gram,err)
 if(err/=0.or.gram>1d-10)call MPI_Abort(MPI_COMM_WORLD,3,s)
 diff=0
 do j=first,last
  diff=max(diff,abs(x(j,j-first+1,1)-(1-dt*dt/4)/(1+dt*dt/4)))
  diff=max(diff,abs(x(3-j,j-first+1,1)+(0d0,1d0)*dt/(1+dt*dt/4)))
 enddo
 call MPI_Allreduce(diff,global_diff,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,s)
 if(global_diff>1d-9)call MPI_Abort(MPI_COMM_WORLD,4,s)
 if(rank==0)print *,'SPINOR_CN_SOLVE_PASS',np,error,gram,global_diff,builds,apps
 call ptcn_owner_destroy(p);call MPI_Finalize(s)
contains
 subroutine reduce_gram(matrix)
 complex(8),intent(inout)::matrix(:,:)
 complex(8)::total(size(matrix,1),size(matrix,2))
 integer::status
 call MPI_Allreduce(matrix,total,size(matrix),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,status)
 matrix=total
 end subroutine
 subroutine action(a,b,refresh)
 complex(8),intent(in)::a(:,:,:)
 complex(8),intent(out)::b(:,:,:)
 logical,intent(in)::refresh
 b=0;b(1,:,:)=a(2,:,:);b(2,:,:)=a(1,:,:)
 end subroutine
 subroutine prec(a,b)
 complex(8),intent(in)::a(:,:,:)
 complex(8),intent(out)::b(:,:,:)
 b=a
 end subroutine
 real(8) function norm(a)
 complex(8),intent(in)::a(:,:,:)
 real(8)::l,t
 integer::status
 l=sum(abs(a)**2);call MPI_Allreduce(l,t,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
 norm=sqrt(t)
 if(mode==1)norm=max(norm,6d-11*sqrt(2d0))
 if(mode==2)norm=max(norm,2d-10*sqrt(2d0))
 end function
 subroutine agree(code)
 integer,intent(inout)::code
 integer::v,status
 call MPI_Allreduce(code,v,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,status);code=v
 end subroutine
 subroutine project(a,b,dv,r,noccupied)
 complex(8),intent(in)::a(:,:,:),b(:,:,:)
 complex(8),intent(out)::r(:,:,:)
 real(8),intent(in)::dv
 integer,intent(in)::noccupied
 integer::code
 call ptcn_owner_project(p,a(:,:,1),b(:,:,1),dv,noccupied,r(:,:,1),code,occup)
 if(code/=0)call MPI_Abort(MPI_COMM_WORLD,5,s)
 end subroutine
end program
