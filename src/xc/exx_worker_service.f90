! Explicit exchange service; SCF owners keep their original communicator and ACE state.
module exx_worker_service
 use mpi,only:MPI_Ibcast,MPI_Wait,MPI_Test,MPI_Comm_rank,MPI_Comm_size,MPI_Abort,MPI_Finalize, &
 MPI_Allreduce,MPI_Bcast,MPI_Wtime,MPI_STATUS_SIZE,MPI_COMM_WORLD,MPI_INTEGER,MPI_MIN,MPI_MAX, &
 MPI_DOUBLE_PRECISION
 use iso_c_binding,only:c_int
 use iso_fortran_env,only:int64,error_unit
 use exx_k_exchange,only:exx_k_kernel
 use exx_symmetry,only:s_exx_symmetry_map
 use exx_spinor_wave_workers,only:apply_wave_workers
 implicit none
 private
 public::service_enter,service_stop,service_enabled,service_exchange,service_fail
 integer,save::worldrank=0,worldsize=1,nowners=1,calls=0
 logical,save::active=.false.
 integer,save::service_divisor=16,service_column_groups=1
 interface
 integer(c_int) function sleep_microseconds(delay) bind(C,name='usleep')
 import c_int
 implicit none
 integer(c_int),value::delay
 end function
 end interface
contains
 subroutine wait_command(command)
 implicit none
 integer,intent(inout)::command
 integer::request,e,status(MPI_STATUS_SIZE),rc
 logical::ready
 ! All ranks match the same nonblocking collective. Only service helpers sleep;
 ! owners retain normal MPI progress and launch exchange without a fixed delay.
 call MPI_Ibcast(command,1,MPI_INTEGER,0,MPI_COMM_WORLD,request,e)
 if(worldrank<nowners)then
 call MPI_Wait(request,status,e)
 else
 do
 call MPI_Test(request,ready,status,e)
 if(ready)exit
 rc=sleep_microseconds(1000_c_int)
 enddo
 endif
 end subroutine
 logical function service_enabled()
 implicit none
 service_enabled=active
 end function
 subroutine service_enter(owner_ranks,tile_divisor,rprime_group_count)
 implicit none
 integer,intent(in)::owner_ranks,tile_divisor,rprime_group_count
 integer::e,command,bad,totalbad,settings(3),lowest(3),highest(3)
 call MPI_Comm_rank(MPI_COMM_WORLD,worldrank,e)
 call MPI_Comm_size(MPI_COMM_WORLD,worldsize,e)
 settings=[owner_ranks,tile_divisor,rprime_group_count]
 bad=0
 if(owner_ranks<1.or.owner_ranks>worldsize.or.tile_divisor<1.or.rprime_group_count<1)bad=1
 if(rprime_group_count>0)then
   if(modulo(worldsize,rprime_group_count)/=0)bad=1
 endif
 call MPI_Allreduce(bad,totalbad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,e)
 call MPI_Allreduce(settings,lowest,3,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,e)
 call MPI_Allreduce(settings,highest,3,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,e)
 if(totalbad/=0.or.any(lowest/=highest))call service_fail('Invalid or inconsistent EXX worker layout')
 nowners=owner_ranks
 service_divisor=tile_divisor
 service_column_groups=rprime_group_count
 active=nowners<worldsize
 if(.not.active)return
 print *,'EXX_SERVICE_LAYOUT',worldrank,nowners,worldsize
 if(worldrank<nowners)return
 do
 call wait_command(command)
 if(command==0)exit
 if(command/=1)call MPI_Abort(MPI_COMM_WORLD,83,e)
 call helper_exchange
 enddo
 print *,'EXX_SERVICE_STOP',worldrank,calls
 call MPI_Finalize(e)
 stop
 end subroutine service_enter
 subroutine service_fail(message)
 implicit none
 character(*),intent(in)::message
 integer::e
 write(error_unit,*)trim(message)
 flush(error_unit)
 call MPI_Abort(MPI_COMM_WORLD,86,e)
 error stop 'MPI_Abort returned unexpectedly'
 end subroutine
 subroutine service_stop
 implicit none
 integer::command
 if(.not.active)return
 command=0
 call wait_command(command)
 end subroutine
 subroutine service_exchange(op,map,es,target,action,ownerstarts,ownercounts,rep_start,owner_rank,owner_comm,ierr)
 implicit none
 type(exx_k_kernel),intent(inout)::op
 type(s_exx_symmetry_map),intent(inout)::map
 complex(8),intent(in)::es(:,:,:),target(:,:,:)
 complex(8),intent(out)::action(:,:,:)
 integer,intent(in)::ownerstarts(:),ownercounts(:),rep_start,owner_rank,owner_comm
 integer,intent(out)::ierr
 integer::command,e,np,rank,bad,totalbad,width,nt
 integer,allocatable::starts(:),counts(:)
 integer(int64)::payload,pairs,density,contract,scratch
 real(8)::t,meta_start,metadata_seconds
 meta_start=MPI_Wtime()
 call MPI_Comm_size(owner_comm,np,e);call MPI_Comm_rank(owner_comm,rank,e)
 bad=0
 if(np/=nowners.or.rank/=worldrank.or.owner_rank/=worldrank)bad=1
 call MPI_Allreduce(bad,totalbad,1,MPI_INTEGER,MPI_MAX,owner_comm,e)
 if(totalbad/=0)call MPI_Abort(MPI_COMM_WORLD,84,e)
 command=1;call wait_command(command)
 width=size(es,2);nt=size(target,2)
 call share_metadata(op,map,width,nt)
 allocate(starts(worldsize),counts(worldsize));starts=map%nfull+1;counts=0
 starts(:nowners)=ownerstarts;counts(:nowners)=ownercounts
 call MPI_Bcast(starts,worldsize,MPI_INTEGER,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(counts,worldsize,MPI_INTEGER,0,MPI_COMM_WORLD,e)
 metadata_seconds=MPI_Wtime()-meta_start
 t=MPI_Wtime()
 call apply_wave_workers(op,map,es,target,action,starts,counts,rep_start,worldrank,MPI_COMM_WORLD,ierr, &
  payload,pairs,density,contract,scratch,service_divisor,service_column_groups)
 calls=calls+1
 print *,'EXX_SERVICE_COMPUTE',worldrank,calls,density,contract,scratch,payload,MPI_Wtime()-t,metadata_seconds
 end subroutine
 subroutine helper_exchange
 implicit none
 type(exx_k_kernel)::op
 type(s_exx_symmetry_map)::map
 complex(8),allocatable::es(:,:,:),target(:,:,:),action(:,:,:)
 integer,allocatable::starts(:),counts(:)
 integer::width,nt,e,ierr
 integer(int64)::payload,pairs,density,contract,scratch
 real(8)::t
 width=0;nt=0
 call share_metadata(op,map,width,nt)
 allocate(starts(worldsize),counts(worldsize))
 call MPI_Bcast(starts,worldsize,MPI_INTEGER,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(counts,worldsize,MPI_INTEGER,0,MPI_COMM_WORLD,e)
 op%phase_start=starts(worldrank+1)
 allocate(op%phase(op%ng,0),es(2*op%ng,width,0),target(2*op%ng,nt,0),action(2*op%ng,nt,0))
 t=MPI_Wtime()
 call apply_wave_workers(op,map,es,target,action,starts,counts,map%nrep+1,worldrank,MPI_COMM_WORLD,ierr, &
  payload,pairs,density,contract,scratch,service_divisor,service_column_groups)
 calls=calls+1
 print *,'EXX_SERVICE_COMPUTE',worldrank,calls,density,contract,scratch,payload,MPI_Wtime()-t
 if(ierr/=0)call MPI_Abort(MPI_COMM_WORLD,85,e)
 end subroutine
 subroutine share_metadata(op,map,width,nt)
 implicit none
 type(exx_k_kernel),intent(inout)::op
 type(s_exx_symmetry_map),intent(inout)::map
 integer,intent(inout)::width,nt
 integer::meta(18),e,t
 logical::helper
 helper=worldrank>=nowners
 if(worldrank==0)meta=[op%n,op%mesh,op%ng,op%nk,op%nq,op%block,map%nfull,map%nrep,map%nsym, &
 map%max_little,shape(op%kernel),width]
 call MPI_Bcast(meta,size(meta),MPI_INTEGER,0,MPI_COMM_WORLD,e)
 t=nt;call MPI_Bcast(t,1,MPI_INTEGER,0,MPI_COMM_WORLD,e)
 if(helper)then
 op%n=meta(1:3);op%mesh=meta(4:6);op%ng=meta(7);op%nk=meta(8);op%nq=meta(9);op%block=meta(10)
 map%nfull=meta(11);map%nrep=meta(12);map%nsym=meta(13);map%max_little=meta(14);map%n=op%n
 allocate(op%point(3,op%ng),op%order(op%nk),op%shift(3,op%nk))
 allocate(op%kernel(0:meta(15)-1,0:meta(16)-1,0:meta(17)-1))
 allocate(map%full_k(3,map%nfull),map%rep_k(3,map%nrep),map%rotation(3,3,map%nsym), &
 map%translation(3,map%nsym),map%first(map%nrep+1),map%owner(map%nfull), &
 map%operations(map%nsym,map%nfull),map%multiplicity(map%nfull),map%grid(op%ng,map%nsym))
 endif
 call MPI_Bcast(op%point,size(op%point),MPI_INTEGER,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(op%order,size(op%order),MPI_INTEGER,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(op%shift,size(op%shift),MPI_INTEGER,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(op%kernel,size(op%kernel),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(map%h,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(map%full_k,size(map%full_k),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(map%rep_k,size(map%rep_k),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(map%rotation,size(map%rotation),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(map%translation,size(map%translation),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(map%first,size(map%first),MPI_INTEGER,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(map%owner,size(map%owner),MPI_INTEGER,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(map%operations,size(map%operations),MPI_INTEGER,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(map%multiplicity,size(map%multiplicity),MPI_INTEGER,0,MPI_COMM_WORLD,e)
 call MPI_Bcast(map%grid,size(map%grid),MPI_INTEGER,0,MPI_COMM_WORLD,e)
 ! Helper shape values are returned through saved scalar storage, avoiding root-only arrays.
 if(helper)then
 width=meta(18);nt=t
 endif
 end subroutine
end module
