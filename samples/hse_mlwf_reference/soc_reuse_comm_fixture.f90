module communication
use mpi
implicit none
interface comm_summation
 module procedure comm_summation_integer,comm_summation_array1d_integer
end interface
interface comm_bcast
module procedure comm_bcast_array1d_integer
end interface
contains
subroutine comm_get_groupinfo(c,r,n)
integer,intent(in)::c
integer,intent(out)::r,n
integer::e
call MPI_Comm_rank(c,r,e);call MPI_Comm_size(c,n,e)
end subroutine
subroutine comm_summation_integer(x,y,c)
integer,intent(in)::x,c
integer,intent(out)::y
integer::e
call MPI_Allreduce(x,y,1,MPI_INTEGER,MPI_SUM,c,e)
end subroutine
subroutine comm_summation_array1d_integer(x,y,n,c)
integer,intent(in)::n,x(n),c
integer,intent(out)::y(n)
integer::e
call MPI_Allreduce(x,y,n,MPI_INTEGER,MPI_SUM,c,e)
end subroutine
subroutine comm_bcast_array1d_integer(x,c,r)
integer,intent(inout)::x(:)
integer,intent(in)::c,r
integer::e
call MPI_Bcast(x,size(x),MPI_INTEGER,r,c,e)
end subroutine
end module
