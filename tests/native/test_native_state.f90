! Functional ISO_C_BINDING host; reads existing GS, never exports checkpoints.
#ifndef SALMONPY_STATE_SYMMETRY
#define SALMONPY_STATE_SYMMETRY 0
#endif
#ifdef SALMONPY_FAST_HSE
#define INFO_SYMBOL "sp_fast_hse_handle_info"
#define CLOSE_SYMBOL "sp_fast_hse_handle_close"
#else
#define INFO_SYMBOL "sp_tddft_handle_info"
#define CLOSE_SYMBOL "sp_tddft_handle_close"
#endif
program native_state_host
 use iso_c_binding
 use mpi
 implicit none
 interface
  integer(c_int) function attach(comm) bind(C,name="sp_attach")
   import c_int,c_double,c_char,c_int64_t
   integer(c_int64_t),value::comm
  end function
  integer(c_int) function detach() bind(C,name="sp_detach")
   import c_int,c_double,c_char,c_int64_t
  end function
  integer(c_int) function directory_enter(path) bind(C,name="sp_io_directory_enter")
   import c_int,c_double,c_char,c_int64_t
   character(c_char),intent(in)::path(*)
  end function
#ifdef SALMONPY_FAST_HSE
  integer(c_int) function hse_initialize(path,dt,capacity,symmetry,handle,j) bind(C,name="sp_fast_hse_handle_init")
   import c_int,c_double,c_char,c_int64_t
   character(c_char),intent(in)::path(*)
   real(c_double),value::dt
   integer(c_int),value::capacity,symmetry
   integer(c_int64_t),intent(out)::handle
   real(c_double),intent(out)::j(3)
  end function
  integer(c_int) function hse_new_state(dt,capacity,handle,j) bind(C,name="sp_fast_hse_model_new")
   import c_int,c_double,c_char,c_int64_t
   real(c_double),value::dt
   integer(c_int),value::capacity
   integer(c_int64_t),intent(out)::handle
   real(c_double),intent(out)::j(3)
  end function
  integer(c_int) function hse_advance(handle,old,a,j,energy) bind(C,name="sp_fast_hse_handle_step")
   import c_int,c_double,c_char,c_int64_t
   integer(c_int64_t),value::handle
   real(c_double),intent(in)::old(3),a(3)
   real(c_double),intent(out)::j(3),energy
  end function
#else
  integer(c_int) function initialize(path,dt,capacity,a,decomposition,handle,j) bind(C,name="sp_tddft_handle_init")
   import c_int,c_double,c_char,c_int64_t
   character(c_char),intent(in)::path(*)
   real(c_double),value::dt
   integer(c_int),value::capacity
   real(c_double),intent(in)::a(3)
   integer(c_int),intent(in)::decomposition(5)
   integer(c_int64_t),intent(out)::handle
   real(c_double),intent(out)::j(3)
  end function
  integer(c_int) function new_state(handle,j) bind(C,name="sp_tddft_model_new")
   import c_int,c_double,c_char,c_int64_t
   integer(c_int64_t),intent(out)::handle
   real(c_double),intent(out)::j(3)
  end function
  integer(c_int) function advance(handle,old,a,j) bind(C,name="sp_tddft_handle_step")
   import c_int,c_double,c_char,c_int64_t
   integer(c_int64_t),value::handle
   real(c_double),intent(in)::old(3),a(3)
   real(c_double),intent(out)::j(3)
  end function
#endif
  integer(c_int) function info(handle,step,a,j) bind(C,name=INFO_SYMBOL)
   import c_int,c_double,c_char,c_int64_t
   integer(c_int64_t),value::handle
   integer(c_int),intent(out)::step
   real(c_double),intent(out)::a(3),j(3)
  end function
  integer(c_int) function retire(handle) bind(C,name=CLOSE_SYMBOL)
   import c_int,c_double,c_char,c_int64_t
   integer(c_int64_t),value::handle
  end function
  subroutine directory_leave() bind(C,name="sp_io_directory_leave")
  end subroutine
 end interface
#ifdef SALMONPY_FAST_HSE
 real(c_double),parameter::host_dt=.001_c_double
 integer,parameter::rejected_old_a=1
#else
 real(c_double),parameter::host_dt=.003_c_double
 integer,parameter::rejected_old_a=10
#endif
 integer::ierr,provided,rank,n,point,first,last,i,k,sibling,u,d
 integer(c_int)::code,step,before_step,after_step,dec(5)
 integer(c_int64_t)::h(2),fresh
 real(c_double)::zero(3),a(3),j(3,2),old(3,2),initial(3),before_a(3),before_j(3)
 real(c_double)::after_a(3),after_j(3),bad(3),rejected(3),history(3,4),previous(3)
 logical::finalized
 character(4096)::work
 character(20)::arg
 call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
 if(ierr/=MPI_SUCCESS.or.provided<MPI_THREAD_FUNNELED)call MPI_Abort(MPI_COMM_WORLD,90,ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,n,ierr)
 call get_command_argument(1,work);call get_command_argument(2,arg);read(arg,*)point
 first=1;last=2
 if(point/=0)then
  first=point;last=point
 endif
 dec=[1,n,1,1,1];zero=0d0;old=0d0;h=0
 call check(attach(int(MPI_COMM_WORLD,c_int64_t)),0)
 if(rank==0)open(newunit=u,file=trim(work)//'/trace.tsv',status='new')
 call enter();code=initialize('rt.inp'//c_null_char,host_dt,16_c_int,zero,dec,h(first),j(:,first))
 call leave_check(code,0);initial=j(:,first)
 if(point==0)then
  call enter();code=new_state(h(2),j(:,2));call leave_check(code,0);call same(j(:,2),initial)
 endif
 do i=1,12
  do k=first,last
   sibling=3-k
   if(point==0)then
    call enter();code=info(h(sibling),before_step,before_a,before_j);call leave_check(code,0)
   endif
   call drive(i,k,a)
   if(i==3.and.k==first)then
    bad=old(:,k);bad(1)=bad(1)+1d0
    call enter();code=advance(h(k),bad,a,rejected);call leave_check(code,rejected_old_a)
    call enter();code=info(h(k),step,after_a,after_j);call leave_check(code,0)
    call check(step,i-1);call same(after_a,old(:,k));call same(after_j,j(:,k))
   endif
   call enter();code=advance(h(k),old(:,k),a,j(:,k));call leave_check(code,0);old(:,k)=a
   if(point==0)then
    call enter();code=info(h(sibling),after_step,after_a,after_j);call leave_check(code,0)
    call check(after_step,before_step);call same(after_a,before_a);call same(after_j,before_j)
    if(k==1.and.i<=4)history(:,i)=j(:,k)
   endif
   if(rank==0)write(u,'(i0,1x,i0,1x,6(es26.17e3,1x))')i,k,a,j(:,k)
  enddo
 enddo
 if(point==0)then
  call enter();code=info(h(2),before_step,before_a,before_j);call leave_check(code,0)
  call enter();code=retire(h(1));call leave_check(code,0)
  call enter();code=info(h(1),step,after_a,after_j);call leave_check(code,12)
  call enter();code=new_state(fresh,after_j);call leave_check(code,0);call same(after_j,initial)
  if(fresh==h(1))call MPI_Abort(MPI_COMM_WORLD,96,ierr)
  call enter();code=info(h(1),step,after_a,after_j);call leave_check(code,12)
  previous=0d0
  do i=1,4
   call drive(i,1,a);call enter();code=advance(fresh,previous,a,after_j);call leave_check(code,0)
   if(any(abs(after_j-history(:,i))>1d-12))call MPI_Abort(MPI_COMM_WORLD,97,ierr)
   previous=a
  enddo
  call enter();code=retire(fresh);call leave_check(code,0)
  call enter();code=info(h(2),after_step,after_a,after_j);call leave_check(code,0)
  call check(after_step,before_step);call same(after_a,before_a);call same(after_j,before_j)
  call enter();code=retire(h(2));call leave_check(code,0)
 else
  call enter();code=retire(h(first));call leave_check(code,0)
  call enter();code=info(h(first),step,after_a,after_j);call leave_check(code,12)
 endif
 call check(detach(),0)
 call MPI_Finalized(finalized,ierr)
 if(finalized)call MPI_Abort(MPI_COMM_WORLD,98,ierr)
 call MPI_Barrier(MPI_COMM_WORLD,ierr)
 if(rank==0)close(u)
 call MPI_Finalize(ierr)
contains
#ifdef SALMONPY_FAST_HSE
 integer(c_int) function initialize(path,dt,capacity,a,decomposition,handle,current) result(status)
  character(c_char),intent(in)::path(*)
  real(c_double),intent(in)::dt,a(3)
  integer(c_int),intent(in)::capacity,decomposition(5)
  integer(c_int64_t),intent(out)::handle
  real(c_double),intent(out)::current(3)
  status=hse_initialize(path,dt,capacity,int(SALMONPY_STATE_SYMMETRY,c_int),handle,current)
 end function
 integer(c_int) function new_state(handle,current) result(status)
  integer(c_int64_t),intent(out)::handle
  real(c_double),intent(out)::current(3)
  status=hse_new_state(host_dt,16_c_int,handle,current)
 end function
 integer(c_int) function advance(handle,previous,field,current) result(status)
  integer(c_int64_t),intent(in)::handle
  real(c_double),intent(in)::previous(3),field(3)
  real(c_double),intent(out)::current(3)
  real(c_double)::energy
  status=hse_advance(handle,previous,field,current,energy)
 end function
#endif
 subroutine check(got,wanted)
  integer,intent(in)::got,wanted
  if(got/=wanted)then
   write(*,*)'native status',got,'expected',wanted
   call MPI_Abort(MPI_COMM_WORLD,90,ierr)
  endif
 end subroutine
 subroutine enter()
  call check(directory_enter(trim(work)//c_null_char),0)
 end subroutine
 subroutine leave_check(status,wanted)
  integer,intent(in)::status,wanted
  call directory_leave();call check(status,wanted)
 end subroutine
 subroutine same(x,y)
  real(c_double),intent(in)::x(3),y(3)
  if(any(x/=y))call MPI_Abort(MPI_COMM_WORLD,91,ierr)
 end subroutine
 subroutine drive(iter,which,x)
  integer,intent(in)::iter,which
  real(c_double),intent(out)::x(3)
  real(c_double)::factor
  factor=1d0
  if(which==2)factor=-.5d0
  x=[factor*mod(iter,5)/1048576d0,-factor*mod(iter,7)/2097152d0,factor*mod(iter,3)/4194304d0]
#if SALMONPY_STATE_SYMMETRY
  ! This fixture's retained symmetry operations preserve the z direction.
  x(1:2)=0d0
#endif
 end subroutine
end program
