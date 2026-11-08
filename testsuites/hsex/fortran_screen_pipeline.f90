program test_pipeline
  use, intrinsic :: iso_fortran_env, only: output_unit
  use mpi
  use hsex_screen_driver, only: screen_compute
  implicit none
  integer::rank,np,err,n(3),mesh(3),qmesh(3),boxes,nv,nb,ng,nk,u,first,count,rows,lo,hi,status,threads
  real(8)::cutoff,cell(3,3),eta
  real(8),allocatable::k(:,:),levels(:,:),pop(:,:),w(:,:),full(:,:),part(:,:)
  complex(8),allocatable::states(:,:,:),h(:,:,:)
  integer :: ncomponent
  character(32)::arg,mode
  logical :: distributed
  call MPI_Init(err);call MPI_Comm_rank(MPI_COMM_WORLD,rank,err);call MPI_Comm_size(MPI_COMM_WORLD,np,err)
  call get_command_argument(1,arg);read(arg,*)threads
  ncomponent=1
  if(command_argument_count()>1)then
    call get_command_argument(2,arg)
    read(arg,*)ncomponent
  endif
  open(newunit=u,file='pipeline.bin',access='stream',form='unformatted',status='old')
  read(u)n,mesh,qmesh,boxes,nv,nb,cutoff,eta,cell
  ng=ncomponent*product(n);nk=product(mesh)
  allocate(k(3,nk),levels(nb,nk),pop(nb,nk),states(ng,nb,nk),h(ng,ng,nk))
  read(u)k,levels,pop,states,h;close(u)
  rows=(boxes**6+np-1)/np;first=rank*rows+1;count=max(0,min(rows,boxes**6-first+1))
  lo=nk*rank/np+1;hi=nk*(rank+1)/np
  allocate(w(count,nk),full(boxes**6,nk),part(boxes**6,nk))
  call get_command_argument(3,mode)
  distributed=mode=='distributed'
  if(mode=='dense-distributed')then
    call screen_compute(n,mesh,qmesh,boxes,cutoff,cell,k,levels,pop,states(:,:,lo:hi),nv,eta,.05d0,.02d0,threads, &
      lo,hi,MPI_COMM_WORLD,rank,np,load_h,first,count,w,status,distributed_states=.true.)
  else if(distributed)then
    call screen_compute(n,mesh,qmesh,boxes,cutoff,cell,k,levels,pop,states(:,:,lo:hi),nv,eta,.05d0,.02d0,threads, &
      lo,hi,MPI_COMM_WORLD,rank,np,load_h,first,count,w,status,apply_h,0.002d0,distributed_states=.true.)
  else if(mode=='replicated')then
    call screen_compute(n,mesh,qmesh,boxes,cutoff,cell,k,levels,pop,states,nv,eta,.05d0,.02d0,threads, &
      lo,hi,MPI_COMM_WORLD,rank,np,load_h,first,count,w,status,apply_h,0.002d0)
  else
    call screen_compute(n,mesh,qmesh,boxes,cutoff,cell,k,levels,pop,states,nv,eta,.05d0,.02d0,threads, &
      lo,hi,MPI_COMM_WORLD,rank,np,load_h,first,count,w,status)
  endif
  if(status/=0)then
    print *, 'SCREEN_FAILURE',rank,status
    flush(output_unit)
    call MPI_Abort(MPI_COMM_WORLD,1,err)
  endif
  part=0d0;part(first:first+count-1,:)=w
  call MPI_Reduce(part,full,size(part),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,err)
  if(rank==0)then
    open(newunit=u,file='pipeline-output.bin',access='stream',form='unformatted',status='replace')
    write(u)full;close(u)
  endif
  call MPI_Finalize(err)
contains
  subroutine apply_h(j,x,y,status)
    implicit none
    integer,intent(in)::j
    complex(8),intent(in)::x(:)
    complex(8),intent(out)::y(:)
    integer,intent(out)::status
    y=matmul(h(:,:,j),x)
    status=0
  end subroutine apply_h

  subroutine load_h(j,a,status)
    implicit none
    integer,intent(in)::j
    complex(8),intent(out)::a(:,:)
    integer,intent(out)::status
    a=h(:,:,j);status=0
  end subroutine
end program
