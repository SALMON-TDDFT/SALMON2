#include "config.h"
! Independent native-core candidate. Each spatial rank owns its FFT box:
! native core plus radius halo, with periodic unique points and native-core return.
module sr_group_fft
use communication,only:comm_get_max,comm_bcast,comm_get_groupinfo,comm_create_group,comm_free_group
! MPI remains for operations/buffer ranks not covered by communication.
#ifdef USE_MPI
use mpi,only: &
 MPI_Allgather, &
 MPI_DOUBLE_COMPLEX, &
 MPI_INTEGER, &
 MPI_Sendrecv, &
 MPI_STATUS_IGNORE
#endif
use structures,only:s_rgrid
use sendrecv_grid_box,only:s_box_halo,init_box_halo,update_box_halo,destroy_box_halo
use,intrinsic::ieee_arithmetic,only:ieee_is_finite
use iso_fortran_env,only:int64
use fftw_blocks,only:mesh_transform,fftw_block_seconds
implicit none
private
public::group_plan,group_prepare,group_apply,group_destroy,ensure_halo_batch,sr_all_finite
interface sr_all_finite
 module procedure finite_real_vector,finite_complex_vector,finite_complex_matrix
end interface
type point_map
 integer,allocatable::idx(:)
end type
type group_plan
 type(point_map),allocatable::wave_send(:),wave_recv(:)
 logical::wave_maps_ready=.false.
 integer(int64)::wave_map_bytes=0
 type(s_box_halo)::halo_plan
 type(s_rgrid)::halo_grid
 logical::use_halo=.false.
 integer::halo_width(3),halo_shape(3),halo_neig(2,3)
 integer::world,rank,peers,group,grank,gpeers,axis(3)
 integer::n(3),m(3),lo(3),gd(3),gc(3),length(3),fm(3),flo(3),origin(3)
 integer,allocatable::los(:,:),origins(:,:),flos(:,:),groups(:),gcoords(:,:)
 logical,allocatable::exchange(:)
 complex(8),allocatable::filter(:,:)
 real(8)::phase(5)=0d0 ! halo, FFT, multiply, core return, total
 real(8)::fft_detail(4)=0d0
 integer(int64)::messages=0,bytes=0
 logical::ready=.false.
end type
contains
real(8) function clock()
  use misc_routines,only:get_wtime
  implicit none
clock=get_wtime()
end function
subroutine group_prepare(p,n,dims,coords,comm,kernel,h,radius,status,halo_enabled)
  implicit none
type(group_plan),intent(out)::p
integer,intent(in)::n(3),dims(3),coords(3),comm
complex(8),intent(in)::kernel(:)
real(8),intent(in)::h(3),radius
integer,intent(out)::status
logical,intent(in),optional::halo_enabled
#ifdef USE_MPI
integer::e,a,j,x,y,z,g,q(3),d(3),gc0(3),core(3),halo(3),color,key,other(2),qlo(3),bad,neig(2,3),side,cc(3),jj
complex(8),allocatable::buf(:),local(:,:),spec(:,:)
p%world=comm;call comm_get_groupinfo(comm,p%rank,p%peers)
bad=0
if(any(n<1).or.any(dims<1).or.product(dims)/=p%peers)bad=1
status=bad;call comm_get_max(status,comm)
if(status/=0)return
if(any(modulo(n,dims)/=0).or.any(coords<0).or.any(coords>=dims))bad=1
if(any(h<=0).or..not.sr_all_finite(h).or..not.ieee_is_finite(radius))bad=1
if(radius<=0)bad=1
p%n=n;p%m=n/dims;p%lo=coords*p%m
if(size(kernel)/=product(p%m))bad=1
if(.not.sr_all_finite(kernel))bad=1
status=bad;call comm_get_max(status,comm)
if(status/=0)return
! One native spatial rank is the center, for every radius.
! Reaching a full periodic axis must not enlarge the FFT communicator.
p%gd=1
p%gc=0;gc0=coords
core=p%m;halo=ceiling(min(radius/h,real(n,8)/2))
p%length=min(n,core+2*halo)
p%origin=gc0*core-halo
where(p%length==n)p%origin=0
p%fm=p%length/p%gd;p%flo=p%gc*p%fm
color=gc0(1)+(dims(1)/p%gd(1))*(gc0(2)+(dims(2)/p%gd(2))*gc0(3))
key=p%gc(1)+p%gd(1)*(p%gc(2)+p%gd(2)*p%gc(3))
p%group=comm_create_group(comm,color,key)
call comm_get_groupinfo(p%group,p%grank,p%gpeers)
do a=1,3
 other=pack([1,2,3],[1,2,3]/=a)
 color=p%gc(other(1))+p%gd(other(1))*p%gc(other(2))
 p%axis(a)=comm_create_group(p%group,color,p%gc(a))
enddo
allocate(p%los(3,p%peers),p%origins(3,p%peers),p%flos(3,p%peers),p%groups(p%peers),p%gcoords(3,p%gpeers))
call MPI_Allgather(p%lo,3,MPI_INTEGER,p%los,3,MPI_INTEGER,comm,e)
call MPI_Allgather(p%origin,3,MPI_INTEGER,p%origins,3,MPI_INTEGER,comm,e)
call MPI_Allgather(p%flo,3,MPI_INTEGER,p%flos,3,MPI_INTEGER,comm,e)
call MPI_Allgather(key,1,MPI_INTEGER,p%groups,1,MPI_INTEGER,comm,e)
call MPI_Allgather(p%gc,3,MPI_INTEGER,p%gcoords,3,MPI_INTEGER,p%group,e)
allocate(p%exchange(p%peers))
do j=1,p%peers
 p%exchange(j)=intersects(p,p%los(:,j),p%origin,p%flo).or.intersects(p,p%lo,p%origins(:,j),p%flos(:,j))
enddo
! Use SALMON face exchanges only while native halo storage is smaller than
! an FFT box and does not duplicate a periodic axis. Otherwise retain original
! distributed native-block delivery, preserving the same masked kernel.
p%halo_width=halo;p%halo_shape=p%m+2*halo
p%use_halo=all(p%halo_shape<=n).and.product(p%halo_shape)<=product(p%length).and..not.all(p%length==n)
if(present(halo_enabled))p%use_halo=p%use_halo.and.halo_enabled
if(p%use_halo)then
 p%halo_grid%is=p%lo;p%halo_grid%ie=p%lo+p%m-1;p%halo_grid%nd=1
 p%halo_grid%is_array=p%lo-halo;p%halo_grid%ie_array=p%lo+p%m-1+halo
 do a=1,3;do side=1,2
  cc=coords;cc(a)=modulo(cc(a)+merge(1,-1,side==1),dims(a))
  neig(side,a)=-1
  do jj=1,p%peers
   if(all(p%los(:,jj)==cc*p%m))neig(side,a)=jj-1
  enddo
 enddo;enddo
 p%halo_neig=neig
 call init_box_halo(p%halo_plan,p%halo_grid,halo,8,comm,neig,status)
 if(status/=0)return
endif
allocate(buf(size(kernel)),local(product(p%fm),1),spec(product(p%fm),1))
local=0
! The kernel origin is displacement zero, independent of the halo origin.
do j=0,p%peers-1
 buf=kernel
 call comm_bcast(buf,comm,j)
 g=0
 do z=0,p%m(3)-1;do y=0,p%m(2)-1;do x=0,p%m(1)-1
  g=g+1;d=modulo(p%los(:,j+1)+[x,y,z]+n/2,n)-n/2
  if(sum((d*h)**2)>radius**2)cycle
  q=modulo(d,p%length)-p%flo
  if(any(q<0).or.any(q>=p%fm))cycle
  local(1+q(1)+p%fm(1)*(q(2)+p%fm(2)*q(3)),1)=buf(g)
 enddo;enddo;enddo
enddo
call mesh_transform(p%length,p%gd,p%gc,p%axis,local,spec,-1,bad)
status=bad;call comm_get_max(status,comm)
if(status/=0)return
p%filter=spec;p%ready=.true.
#else
! Native-core group FFT communication requires MPI; ordinary HSE remains available.
status=1
#endif
end subroutine
logical function intersects(p,lo,origin,flo) result(ok)
  implicit none
type(group_plan),intent(in)::p
integer,intent(in)::lo(3),origin(3),flo(3)
integer::a,k,q
logical::hit
ok=.true.
do a=1,3
 hit=.false.
 do k=0,p%m(a)-1
  q=modulo(lo(a)+k-origin(a),p%n(a))-flo(a)
  if(q>=0.and.q<p%fm(a))hit=.true.
 enddo
 ok=ok.and.hit
enddo
end function
subroutine ensure_halo_batch(p,batch,status)
  implicit none
type(group_plan),intent(inout)::p
integer,intent(in)::batch
integer,intent(out)::status
status=0
if(.not.p%use_halo)return
if(p%halo_plan%nb==batch)return
call destroy_box_halo(p%halo_plan)
call init_box_halo(p%halo_plan,p%halo_grid,p%halo_width,batch,p%world,p%halo_neig,status)
end subroutine

subroutine group_apply(p,source,targets,action,status)
  implicit none
type(group_plan),intent(inout)::p
complex(8),intent(in)::source(:),targets(:,:)
complex(8),intent(out)::action(:,:)
integer,intent(out)::status
#ifdef USE_MPI
integer::e,first,nb,k,g,j,x,y,z,q(3),round,partner,rounds,bad
integer,parameter::batch=8
complex(8),allocatable::halo_data(:,:,:,:),halo_chunk(:,:,:,:)
logical,allocatable::filled(:)
integer::hc(3),ii,gr
complex(8),allocatable::send(:,:),recv(:,:),tile(:,:),spec(:,:),back(:,:),chunk(:,:)
real(8)::beg,t,detail(4)
action=0;status=1
if(.not.p%ready)return
bad=0
if(size(source)/=product(p%m).or.size(targets,1)/=size(source))bad=1
if(any(shape(action)/=shape(targets)))bad=1
if(.not.sr_all_finite(source))bad=1
if(.not.sr_all_finite(targets))bad=1
status=bad;call comm_get_max(status,p%world)
if(status/=0)return
! Same target count contract is checked before any target-dependent collective.
nb=size(targets,2);k=nb;call comm_get_max(k,p%world)
bad=merge(0,1,nb==k);status=bad;call comm_get_max(status,p%world)
if(status/=0)return
beg=clock();detail=fftw_block_seconds
allocate(send(size(source),batch),recv(size(source),batch),tile(product(p%fm),batch), &
 spec(product(p%fm),batch),back(product(p%fm),batch),chunk(product(p%fm),batch))
if(p%use_halo)then
 allocate(halo_data(0:p%halo_shape(1)-1,0:p%halo_shape(2)-1,0:p%halo_shape(3)-1,batch), &
 halo_chunk(0:p%halo_shape(1)-1,0:p%halo_shape(2)-1,0:p%halo_shape(3)-1,batch),filled(product(p%fm)))
endif
rounds=1
do while(rounds<p%peers);rounds=rounds*2;enddo
do first=1,size(targets,2),batch
 nb=min(batch,size(targets,2)-first+1);send=0;tile=0
 t=clock()
 do k=1,nb;send(:,k)=conjg(source)*targets(:,first+k-1);enddo
 if(p%use_halo)then
  halo_data=0;g=0
  do z=0,p%m(3)-1;do y=0,p%m(2)-1;do x=0,p%m(1)-1
   g=g+1;hc=[x,y,z]+p%halo_width
   halo_data(hc(1),hc(2),hc(3),:nb)=send(g,:nb)
  enddo;enddo;enddo
  call ensure_halo_batch(p,batch,status)
  if(status/=0)return
  call update_box_halo(p%halo_plan,p%halo_grid,halo_data)
  do ii=1,size(p%halo_plan%step);do gr=1,3;do j=1,2
   partner=p%halo_plan%step(ii)%neig(j,gr)
   if(partner<0.or.partner==p%rank)cycle
   p%messages=p%messages+1
   p%bytes=p%bytes+16_int64*size(p%halo_plan%step(ii)%cache(1,j,gr)%zbuf,kind=int64)
  enddo;enddo;enddo
  filled=.false.
  ! Each fixed subgroup member contributes only points in its native halo.
  ! Choose the first available copy for every FFT point; halos are not summed.
  do round=0,p%gpeers-1
   partner=ieor(p%grank,round)
   if(partner==p%grank)then
    halo_chunk=halo_data
   else
    call MPI_Sendrecv(halo_data,size(halo_data),MPI_DOUBLE_COMPLEX,partner,2, &
     halo_chunk,size(halo_chunk),MPI_DOUBLE_COMPLEX,partner,2,p%group,MPI_STATUS_IGNORE,e)
    p%messages=p%messages+1;p%bytes=p%bytes+16_int64*size(halo_data,kind=int64)
   endif
   ! Recover native origin of this subgroup member from its group coordinates.
   hc=p%lo-p%gc*p%m+p%gcoords(:,partner+1)*p%m-p%halo_width
   ii=0
   do z=0,p%fm(3)-1;do y=0,p%fm(2)-1;do x=0,p%fm(1)-1
    ii=ii+1
    if(filled(ii))cycle
    q=modulo(p%origin+p%flo+[x,y,z]-hc,p%n)
    if(any(q>=p%halo_shape))cycle
    tile(ii,:nb)=halo_chunk(q(1),q(2),q(3),:nb);filled(ii)=.true.
   enddo;enddo;enddo
  enddo
  bad=merge(0,1,all(filled));status=bad;call comm_get_max(status,p%world)
  if(status/=0)return
 else
 do round=0,rounds-1
  partner=ieor(p%rank,round)
  if(partner>=p%peers)cycle
  if(.not.p%exchange(partner+1))cycle
  if(partner==p%rank)then
   recv=send
  else
   call MPI_Sendrecv(send,size(send),MPI_DOUBLE_COMPLEX,partner,0,recv,size(recv),MPI_DOUBLE_COMPLEX,partner, &
    0,p%world,MPI_STATUS_IGNORE,e)
   p%messages=p%messages+1;p%bytes=p%bytes+16_int64*size(send,kind=int64)
  endif
  g=0
  do z=0,p%m(3)-1;do y=0,p%m(2)-1;do x=0,p%m(1)-1
   g=g+1;q=modulo(p%los(:,partner+1)+[x,y,z]-p%origin,p%n)-p%flo
   if(any(q<0).or.any(q>=p%fm))cycle
   j=1+q(1)+p%fm(1)*(q(2)+p%fm(2)*q(3));tile(j,:nb)=recv(g,:nb)
  enddo;enddo;enddo
 enddo
 endif
 p%phase(1)=p%phase(1)+clock()-t;t=clock()
 call mesh_transform(p%length,p%gd,p%gc,p%axis,tile,spec,-1,bad)
 status=bad;call comm_get_max(status,p%world)
 if(status/=0)return
 p%phase(2)=p%phase(2)+clock()-t;t=clock()
 do k=1,nb;spec(:,k)=spec(:,k)*p%filter(:,1);enddo
 p%phase(3)=p%phase(3)+clock()-t;t=clock()
 call mesh_transform(p%length,p%gd,p%gc,p%axis,spec,back,1,bad)
 status=bad;call comm_get_max(status,p%world)
 if(status/=0)return
 p%phase(2)=p%phase(2)+clock()-t;t=clock()
 ! Only one remote FFT chunk is held at a time. Each native core point belongs
 ! to exactly one FFT chunk, even when the halo crosses the global boundary.
 if(all(p%fm==p%m).and.all(p%flo==p%lo).and.all(p%origin==0))then
  do k=1,nb;action(:,first+k-1)=-source*back(:,k);enddo
 else
 do round=0,p%gpeers-1
  partner=ieor(p%grank,round)
  if(partner==p%grank)then
   chunk=back
  else
   call MPI_Sendrecv(back,size(back),MPI_DOUBLE_COMPLEX,partner,1,chunk,size(chunk),MPI_DOUBLE_COMPLEX, &
    partner,1,p%group,MPI_STATUS_IGNORE,e)
   p%messages=p%messages+1;p%bytes=p%bytes+16_int64*size(back,kind=int64)
  endif
  g=0
  do z=0,p%m(3)-1;do y=0,p%m(2)-1;do x=0,p%m(1)-1
   g=g+1;q=modulo(p%lo+[x,y,z]-p%origin,p%n)-p%gcoords(:,partner+1)*p%fm
   if(any(q<0).or.any(q>=p%fm))cycle
   j=1+q(1)+p%fm(1)*(q(2)+p%fm(2)*q(3))
   action(g,first:first+nb-1)=-source(g)*chunk(j,:nb)
  enddo;enddo;enddo
 enddo
 endif
 p%phase(4)=p%phase(4)+clock()-t
enddo
p%phase(5)=p%phase(5)+clock()-beg;p%fft_detail=p%fft_detail+fftw_block_seconds-detail
status=0
#else
action=0d0;status=1
#endif
end subroutine
subroutine group_destroy(p)
  implicit none
type(group_plan),intent(inout)::p
integer::a,e
if(.not.p%ready)return
if(p%use_halo)call destroy_box_halo(p%halo_plan)
do a=1,3;call comm_free_group(p%axis(a));enddo
call comm_free_group(p%group);p%ready=.false.
if(allocated(p%wave_send))deallocate(p%wave_send,p%wave_recv)
p%wave_maps_ready=.false.
if(allocated(p%filter))deallocate(p%filter)
if(allocated(p%los))deallocate(p%los,p%origins,p%flos,p%groups,p%gcoords,p%exchange)
end subroutine

! Scalar IEEE inquiries avoid compiler-dependent array-expression lowering.
logical function finite_real_vector(values) result(ok)
  implicit none
  real(8),intent(in)::values(:)
  integer::i
  ok=.false.
  do i=1,size(values)
    if(.not.ieee_is_finite(values(i)))return
  enddo
  ok=.true.
end function
logical function finite_complex_vector(values) result(ok)
  implicit none
  complex(8),intent(in)::values(:)
  integer::i
  ok=.false.
  do i=1,size(values)
    if(.not.ieee_is_finite(real(values(i),8)))return
    if(.not.ieee_is_finite(aimag(values(i))))return
  enddo
  ok=.true.
end function
logical function finite_complex_matrix(values) result(ok)
  implicit none
  complex(8),intent(in)::values(:,:)
  integer::i,j
  ok=.false.
  do j=1,size(values,2)
    do i=1,size(values,1)
      if(.not.ieee_is_finite(real(values(i,j),8)))return
      if(.not.ieee_is_finite(aimag(values(i,j))))return
    enddo
  enddo
  ok=.true.
end function
end module
