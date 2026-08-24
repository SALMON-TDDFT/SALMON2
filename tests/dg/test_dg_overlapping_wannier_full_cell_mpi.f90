#include "config.h"
program test_dg_overlapping_wannier_full_cell_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_overlapping_wannier_full_cell,only:project_dg_full_cell_hamiltonian_tiles,&
    s_dg_full_cell_redistribution_schedule,initialize_dg_full_cell_redistribution,&
    apply_dg_full_cell_redistribution_forward,apply_dg_full_cell_redistribution_reverse,&
    clear_dg_full_cell_redistribution
  implicit none
  integer,parameter::ngrid=11,nbasis=5
  integer::comm,rank,nproc,ierr,p,i,j,k,nlocal,nowned,tile,index,component_mode,ndestination,&
    destination_first,destination_last
  integer(int64),allocatable::spatial_ids(:),row_ids(:),destination_ids(:),saved_ids(:)
  real(real64),allocatable::weights(:)
  complex(real64),allocatable::basis(:,:),rows(:,:),reference(:,:),global_basis(:,:),global_h(:,:),&
    source_values(:,:),destination_values(:,:),roundtrip_values(:,:)
  real(real64)::x,pi,defect
  integer(int64)::workspace
  logical::ok
  character(256)::message
  type(s_dg_full_cell_redistribution_schedule)::schedule
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  pi=acos(-1d0)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,ngrid)])
  nowned=count([(mod(i-1,nproc)==rank,i=1,nbasis)])
  allocate(spatial_ids(nlocal),weights(nlocal),basis(nbasis,nlocal),row_ids(nowned))
  index=0
  do p=1,ngrid
    if(mod(p-1,nproc)/=rank)cycle
    index=index+1;spatial_ids(index)=p;weights(index)=1d0/ngrid;x=2d0*pi*(p-1)/ngrid
    do i=1,nbasis
      basis(i,index)=exp(cmplx(0d0,real(i-3,real64)*x,real64))/sqrt(real(ngrid,real64))
    enddo
  enddo
  index=0
  do i=1,nbasis
    if(mod(i-1,nproc)/=rank)cycle
    index=index+1;row_ids(index)=i
  enddo

  ! The source uses cyclic ownership while the hpsi destination uses contiguous
  ! real-space ownership.  Their local orders deliberately differ on MPI>1.
  destination_first=rank*ngrid/nproc+1
  destination_last=(rank+1)*ngrid/nproc
  ndestination=max(0,destination_last-destination_first+1)
  allocate(destination_ids(ndestination),source_values(3,nlocal),&
    destination_values(3,ndestination),roundtrip_values(3,nlocal))
  do p=1,ndestination;destination_ids(p)=destination_first+p-1;enddo
  do p=1,nlocal;do i=1,3
    source_values(i,p)=cmplx(real(100*i+int(spatial_ids(p)),real64),&
      real(-10*i-int(spatial_ids(p)),real64),real64)
  enddo;enddo
  call initialize_dg_full_cell_redistribution(comm,ngrid,spatial_ids,destination_ids,&
    schedule,workspace,ok,message)
  call require(ok,trim(message))
  call apply_dg_full_cell_redistribution_forward(schedule,source_values,destination_values,ok,message)
  call require(ok,trim(message))
  defect=0d0
  do p=1,ndestination;do i=1,3
    defect=max(defect,abs(destination_values(i,p)-cmplx(real(100*i+int(destination_ids(p)),real64),&
      real(-10*i-int(destination_ids(p)),real64),real64)))
  enddo;enddo
  call require(defect<1d-14,'forward redistribution does not preserve physical-ID values')
  call apply_dg_full_cell_redistribution_reverse(schedule,destination_values,roundtrip_values,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(roundtrip_values-source_values))<1d-14,&
    'reverse redistribution does not restore source ordering')
  source_values=2d0*source_values
  call apply_dg_full_cell_redistribution_forward(schedule,source_values,destination_values,ok,message)
  call require(ok,trim(message))
  call apply_dg_full_cell_redistribution_reverse(schedule,destination_values,roundtrip_values,ok,message)
  call require(ok.and.maxval(abs(roundtrip_values-source_values))<1d-14,&
    'cached redistribution schedule cannot be reused')
  call clear_dg_full_cell_redistribution(schedule)

  allocate(saved_ids(size(spatial_ids)),source=spatial_ids)
  if(rank==0)spatial_ids(1)=2_int64
  call initialize_dg_full_cell_redistribution(comm,ngrid,spatial_ids,destination_ids,&
    schedule,workspace,ok,message)
  call require(.not.ok,'duplicate/missing source physical IDs were accepted')
  spatial_ids=saved_ids
  call clear_dg_full_cell_redistribution(schedule)
  if(rank==0.and.ndestination>0)destination_ids(1)=2_int64
  call initialize_dg_full_cell_redistribution(comm,ngrid,spatial_ids,destination_ids,&
    schedule,workspace,ok,message)
  call require(.not.ok,'duplicate/missing destination physical IDs were accepted')
  do p=1,ndestination;destination_ids(p)=destination_first+p-1;enddo
  call clear_dg_full_cell_redistribution(schedule)
  if(rank==0)spatial_ids(1)=int(ngrid+1,int64)
  call initialize_dg_full_cell_redistribution(comm,ngrid,spatial_ids,destination_ids,&
    schedule,workspace,ok,message)
  call require(.not.ok,'out-of-range source physical ID was accepted')
  spatial_ids=saved_ids
  call clear_dg_full_cell_redistribution(schedule)
  allocate(global_basis(nbasis,ngrid),global_h(nbasis,ngrid),reference(nbasis,nbasis))
  do p=1,ngrid
    x=2d0*pi*(p-1)/ngrid
    do i=1,nbasis
      global_basis(i,p)=exp(cmplx(0d0,real(i-3,real64)*x,real64))/sqrt(real(ngrid,real64))
    enddo
  enddo
  do component_mode=1,3
    call apply_reference(global_basis,global_h)
    do j=1,nbasis;do i=1,nbasis
      reference(i,j)=sum(conjg(global_basis(i,:))*global_h(j,:))/ngrid
    enddo;enddo
    do tile=1,3
      call project_dg_full_cell_hamiltonian_tiles(comm,ngrid,spatial_ids,weights,basis,row_ids,tile,&
        apply_tile,rows,workspace,ok,message)
      call require(ok,trim(message))
      defect=0d0
      do i=1,nowned
        defect=max(defect,maxval(abs(rows(i,:)-reference(int(row_ids(i)),:))))
      enddo
      call require(defect<2d-13,'tiled component rows differ from dense reference')
    enddo
  enddo
  component_mode=3;call apply_reference(global_basis,global_h)
  do j=1,nbasis;do i=1,nbasis
    reference(i,j)=sum(conjg(global_basis(i,:))*global_h(j,:))/ngrid
  enddo;enddo
  if(rank==0)then
    write(*,'(a,i0,a,4(es24.16,1x))')'FULL_CELL ranks=',nproc,' values=',&
      real(reference(1,2)),aimag(reference(1,2)),real(reference(4,5)),aimag(reference(4,5))
    write(*,'(a,i0,a)')'PASS full-cell tiled projection on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine apply_tile(tile_in,tile_out,callback_ok)
    complex(real64),intent(in)::tile_in(:,:)
    complex(real64),intent(out)::tile_out(:,:)
    logical,intent(out)::callback_ok
    complex(real64),allocatable::local_global(:,:),all_global(:,:)
    integer::q,t,loc
    allocate(local_global(size(tile_in,1),ngrid),all_global(size(tile_in,1),ngrid))
    local_global=(0d0,0d0)
    do loc=1,nlocal
      local_global(:,int(spatial_ids(loc)))=tile_in(:,loc)
    enddo
    call MPI_Allreduce(local_global,all_global,size(all_global),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do loc=1,nlocal
      q=int(spatial_ids(loc))
      do t=1,size(tile_in,1)
        select case(component_mode)
        case(1)
          tile_out(t,loc)=all_global(t,q)-0.25d0*(all_global(t,1+mod(q,ngrid))+&
            all_global(t,1+mod(q-2+ngrid,ngrid)))
        case(2)
          tile_out(t,loc)=(0.3d0+0.2d0*cos(2d0*pi*(q-1)/ngrid))*all_global(t,q)
        case default
          tile_out(t,loc)=all_global(t,q)-0.25d0*(all_global(t,1+mod(q,ngrid))+&
            all_global(t,1+mod(q-2+ngrid,ngrid)))+&
            (0.3d0+0.2d0*cos(2d0*pi*(q-1)/ngrid))*all_global(t,q)
        end select
      enddo
    enddo
    callback_ok=ierr==MPI_SUCCESS
  end subroutine
  subroutine apply_reference(input,output)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    integer::q,t
    do q=1,ngrid;do t=1,size(input,1)
      select case(component_mode)
      case(1)
        output(t,q)=input(t,q)-0.25d0*(input(t,1+mod(q,ngrid))+input(t,1+mod(q-2+ngrid,ngrid)))
      case(2)
        output(t,q)=(0.3d0+0.2d0*cos(2d0*pi*(q-1)/ngrid))*input(t,q)
      case default
        output(t,q)=input(t,q)-0.25d0*(input(t,1+mod(q,ngrid))+input(t,1+mod(q-2+ngrid,ngrid)))+&
          (0.3d0+0.2d0*cos(2d0*pi*(q-1)/ngrid))*input(t,q)
      end select
    enddo;enddo
  end subroutine
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::lf,gf
    lf=merge(0,1,condition);call MPI_Allreduce(lf,gf,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(gf/=0)error stop label
  end subroutine
end program
