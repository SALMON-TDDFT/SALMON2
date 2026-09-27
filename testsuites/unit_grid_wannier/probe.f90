#include "config.h"
program probe
#ifdef USE_MPI
 use mpi
#endif
 use hse_grid_wannier
 implicit none
 type(s_hse_grid_wannier) :: state,reference,full_state,masked_state,plain
 integer,parameter :: ng=512,n=2
 complex(8) :: all_c(ng,n),all_extra(ng),rotation(n,n)
 real(8) :: all_position(3,ng),length(3),dv,q,radius
 complex(8),allocatable :: c(:,:),extra(:),source(:,:),initial(:,:),other(:,:),candidate(:,:),final(:,:)
 real(8),allocatable :: position(:,:)
 integer :: rank,np,comm,ierr,lo,hi,i,j,x,y,z,g,step,mode,irad
 character(32) :: arg
 rank=0;np=1;comm=0
#ifdef USE_MPI
 call MPI_Init(ierr);comm=MPI_COMM_WORLD
 call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,np,ierr)
#endif
 mode=0;call get_command_argument(1,arg);if(len_trim(arg)>0)read(arg,*)mode
 dv=.125d0;length=4d0;g=0
 do z=0,7;do y=0,7;do x=0,7
  g=g+1;all_position(:,g)=.5d0*[x,y,z]
  all_c(g,:)=0d0;all_extra(g)=0d0
  if(x<=2)all_c(g,1)=exp(-sum((all_position(:,g)-[.5d0,1d0,1d0])**2))
  if(x>=5)all_c(g,2)=exp(-sum((all_position(:,g)-[3d0,1d0,1d0])**2))
  if(x==4)all_extra(g)=exp(-sum((all_position(:,g)-[2d0,1d0,1d0])**2))
 enddo;enddo;enddo
 do j=1,n;all_c(:,j)=all_c(:,j)/sqrt(sum(abs(all_c(:,j))**2)*dv);enddo
 all_extra=all_extra/sqrt(sum(abs(all_extra)**2)*dv)
 lo=ng*rank/np+1;hi=ng*(rank+1)/np
 c=all_c(lo:hi,:);extra=all_extra(lo:hi);position=all_position(:,lo:hi)
 if(mode==1)then
  call grid_wannier_source(state,c,position,length,dv,comm,.false.,1d0,100,1d-8,1,.false.,source)
  error stop 'Expected disabled-radius rejection'
 endif
 call grid_wannier_source(plain,c,position,length,dv,comm,.false.,0d0,100,1d-8,1,.false.,source)
 call assert_close(source,c,'disabled source changes C')
 call grid_wannier_source(state,c,position,length,dv,comm,.true.,0d0,100,1d-8,1,.false.,source)
 initial=source
 call assert_close(matmul(source,conjg(transpose(source))),matmul(c,conjg(transpose(c))), &
   'full-source density differs')
 reference=state
 call grid_wannier_stage(state,0)
 candidate=c;candidate(:,1)=cos(.23d0)*c(:,1)+sin(.23d0)*extra
 call grid_wannier_source(state,candidate,position,length,dv,comm,.true.,0d0,100,1d-8,1,.false.,source)
 call grid_wannier_stage(state,1);call grid_wannier_stage(state,2)
 final=c;final(:,2)=cos(.17d0)*c(:,2)+sin(.17d0)*extra
 call grid_wannier_source(state,final,position,length,dv,comm,.true.,0d0,100,1d-8,1,.false.,source)
 call grid_wannier_stage(reference,0);call grid_wannier_stage(reference,1);call grid_wannier_stage(reference,2)
 call grid_wannier_source(reference,final,position,length,dv,comm,.true.,0d0,100,1d-8,1,.false.,other)
 call assert_close(source,other,'predictor gauge leaked into accepted source')
 ! Long stationary occupied-space rotation: only transported frames anchor U.
 do irad=0,1
 radius=.8d0*irad
 call grid_wannier_reset(state)
 call grid_wannier_source(state,c,position,length,dv,comm,.true.,radius,100,1d-8,4,.false.,initial)
 do step=1,101
  q=.037d0*step
  rotation(1,1)=cos(q);rotation(2,2)=cos(q)
  rotation(1,2)=cmplx(0d0,sin(q),8);rotation(2,1)=rotation(1,2)
  candidate=matmul(c,rotation)
  call grid_wannier_stage(state,0)
  call grid_wannier_source(state,candidate,position,length,dv,comm,.true.,radius,100,1d-8,4,.false.,source)
  call grid_wannier_stage(state,1);call grid_wannier_stage(state,2)
  call grid_wannier_source(state,candidate,position,length,dv,comm,.true.,radius,100,1d-8,4,.false.,source)
  if(mod(step-1,4)==0)call assert_close(source,initial,'transport anchor drifts at refresh')
 enddo
 enddo
 ! Radius does not alter initialization U; masked sources retain original amplitudes.
 call grid_wannier_source(full_state,c,position,length,dv,comm,.true.,0d0,100,1d-8,1,.false.,other)
 call grid_wannier_source(masked_state,c,position,length,dv,comm,.true.,.8d0,100,1d-8,1,.false.,source)
 do j=1,n;do i=1,size(c,1)
  if(abs(source(i,j))>0d0.and.abs(source(i,j)-other(i,j))>1d-10)error stop 'Mask renormalized WF'
 enddo;enddo
 call sphere_cases
 if(rank==0)print *,'PASS mesh MLWF source, predictor rollback, U cadence, 3D periodic mask and protection'
#ifdef USE_MPI
 call MPI_Finalize(ierr)
#endif
contains
 subroutine assert_close(a,b,message)
  implicit none
  complex(8),intent(in) :: a(:,:),b(:,:)
  character(*),intent(in) :: message
  real(8) :: error,total
  error=0d0;if(size(a)>0)error=maxval(abs(a-b));total=error
#ifdef USE_MPI
  call MPI_Allreduce(error,total,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
#endif
  if(total>1d-10)then
   print *,message,total
   error stop 'Grid MLWF regression'
  endif
 end subroutine
 subroutine sphere_cases
  implicit none
  type(s_hse_grid_wannier) :: one,uncut
  complex(8) :: global(ng,1)
  complex(8),allocatable :: local(:,:),cut(:,:),whole(:,:)
  real(8) :: pos(3,ng),ll(3)
  integer :: a,ix,jj,iu,id,protected
  real(8) :: norm,sphere,fraction
  character(256) :: header
  ll=4d0
  do a=1,3
   call grid_wannier_reset(one);call grid_wannier_reset(uncut)
   global=0d0;pos=all_position
   ! Periodic neighbor at 3.5, far tail at 1.5 along each tested axis.
   pos(:,1)=0d0;pos(:,2)=0d0;pos(a,2)=3.5d0;pos(:,3)=0d0;pos(a,3)=1.5d0
   global(1,1)=sqrt(.99d0/dv);global(2,1)=sqrt(.005d0/dv);global(3,1)=sqrt(.005d0/dv)
   local=global(lo:hi,:)
   call grid_wannier_source(uncut,local,pos(:,lo:hi),ll,dv,comm,.true.,0d0,20,1d-8,1,.false.,whole)
   call grid_wannier_source(one,local,pos(:,lo:hi),ll,dv,comm,.true.,.6d0,20,1d-8,1,.false.,cut)
   do jj=lo,hi
    ix=jj-lo+1
    if(jj==3)then
     if(abs(cut(ix,1))>1d-12)error stop '3D mask missed far tail'
    else
     if(abs(cut(ix,1)-whole(ix,1))>1d-10)error stop 'Periodic nearest neighbor was cut'
    endif
   enddo
   if(rank==0)then
    open(newunit=iu,file='grid_mlwf_radius.dat',status='old')
    read(iu,'(a)')header;read(iu,'(a)')header
    read(iu,*)id,norm,sphere,fraction,protected;close(iu)
    if(abs(fraction-.995d0)>1d-10.or.protected/=0)error stop 'Wrong radius diagnostic'
   endif
  enddo
  call grid_wannier_reset(one)
  do jj=1,ng
   global(jj,1)=sqrt((1d0+.05d0*sum(cos(2d0*acos(-1d0)*all_position(:,jj)/ll)))/(ng*dv))
  enddo
  local=global(lo:hi,:)
  call grid_wannier_source(one,local,position,ll,dv,comm,.true.,.01d0,20,1d-8,1,.false.,cut)
  ! Near-uniform state has unreliable circular centers; all grid points are protected.
  if(maxval(abs(abs(cut)-abs(local)))>1d-10)error stop 'Protected state truncated'
 end subroutine
end program
