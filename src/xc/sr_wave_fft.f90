#include "config.h"
module sr_wave_fft
use misc_routines,only:sr_walltime=>get_wtime
use communication,only:comm_get_max
! MPI remains for operations/buffer ranks not covered by communication.
#ifdef USE_MPI
use mpi,only: &
 MPI_Allreduce, &
 MPI_DOUBLE_COMPLEX, &
 MPI_INTEGER, &
 MPI_MAX, &
 MPI_MIN, &
 MPI_Sendrecv, &
 MPI_STATUS_IGNORE
#endif
use sr_group_fft,only:group_plan,ensure_halo_batch,sr_all_finite
use sendrecv_grid_box,only:update_box_halo
use fftw_blocks,only:mesh_transform
use iso_fortran_env,only:int64
implicit none
private
public::wave_apply,wave_apply_pairs,wave_screen_supported,prepare_wave_maps
contains
pure logical function wave_screen_supported(mode,tolerance)
  implicit none
integer,intent(in)::mode
real(8),intent(in)::tolerance
! Zero-tolerance exact support can evaluate additional zero pairs without
! changing the exchange operator. Positive omission tolerances remain rejected.
wave_screen_supported=mode==0.or.(mode==2.and.tolerance==0d0)
end function
! Fixed geometry maps persist for all source blocks of this group plan.
subroutine point_intersection(p,sender,dest,source_index,tile_index)
  implicit none
type(group_plan),intent(in)::p
integer,intent(in)::sender,dest
integer,allocatable,intent(out)::source_index(:),tile_index(:)
integer::x,y,z,j,q(3),global_point(3),hc(3),base(3),provider,round,nfound,source_shape(3)
integer,allocatable::si(:),ti(:)
allocate(si(product(p%fm)),ti(product(p%fm)));nfound=0;j=0
base=p%lo-p%gc*p%m
if(p%use_halo)then
 source_shape=p%halo_shape
else
 source_shape=p%m
endif
do z=0,p%fm(3)-1;do y=0,p%fm(2)-1;do x=0,p%fm(1)-1
 j=j+1
 if(p%use_halo)then
  global_point=modulo(p%origin+p%gcoords(:,dest+1)*p%fm+[x,y,z],p%n)
  ! Preserve original first-copy selection in destination XOR partner order.
  do round=0,p%gpeers-1
   provider=ieor(dest,round)
   hc=base+p%gcoords(:,provider+1)*p%m-p%halo_width
   q=modulo(global_point-hc,p%n)
   if(all(q<p%halo_shape))exit
  enddo
  if(round==p%gpeers.or.provider/=sender)cycle
 else
  global_point=modulo(p%origins(:,dest+1)+p%flos(:,dest+1)+[x,y,z],p%n)
  q=modulo(global_point-p%los(:,sender+1),p%n)
  if(any(q>=p%m))cycle
 endif
 nfound=nfound+1;si(nfound)=1+q(1)+source_shape(1)*(q(2)+source_shape(2)*q(3));ti(nfound)=j
enddo;enddo;enddo
allocate(source_index(nfound),tile_index(nfound));source_index=si(:nfound);tile_index=ti(:nfound)
end subroutine
subroutine prepare_wave_maps(p)
  implicit none
type(group_plan),intent(inout)::p
integer::peers,rank,partner
integer,allocatable::si(:),ti(:)
if(p%wave_maps_ready)return
peers=p%peers;rank=p%rank
if(p%use_halo)then
 peers=p%gpeers;rank=p%grank
endif
allocate(p%wave_send(peers),p%wave_recv(peers));p%wave_map_bytes=0
do partner=0,peers-1
 call point_intersection(p,rank,partner,si,ti)
 call move_alloc(si,p%wave_send(partner+1)%idx);deallocate(ti)
 call point_intersection(p,partner,rank,si,ti)
 call move_alloc(ti,p%wave_recv(partner+1)%idx);deallocate(si)
 p%wave_map_bytes=p%wave_map_bytes+int(storage_size(rank)/8,int64)* &
  (size(p%wave_send(partner+1)%idx,kind=int64)+size(p%wave_recv(partner+1)%idx,kind=int64))
enddo
p%wave_maps_ready=.true.
end subroutine
subroutine packed_delivery(p,values,tile,status,halo)
  implicit none
type(group_plan),intent(inout)::p
complex(8),intent(in)::values(:,:)
complex(8),intent(out)::tile(:,:)
integer,intent(out)::status
complex(8),intent(in),optional::halo(0:,0:,0:,:)
#ifdef USE_MPI
complex(8),allocatable::sb(:,:),rb(:,:)
logical::filled(product(p%fm))
integer::peers,rank,comm,rounds,round,partner,ns,nr,nb,k,j,idx,x,y,z,e,bad
peers=p%peers;rank=p%rank;comm=p%world
if(p%use_halo)then
 peers=p%gpeers;rank=p%grank;comm=p%group
endif
nb=size(values,2);tile=0;filled=.false.;rounds=1
 do while(rounds<peers);rounds=rounds*2;enddo
do round=0,rounds-1
 partner=ieor(rank,round)
 if(partner>=peers)cycle
 ns=size(p%wave_send(partner+1)%idx);nr=size(p%wave_recv(partner+1)%idx)
 if(ns==0.and.nr==0)cycle
 allocate(sb(ns,nb),rb(nr,nb))
 do j=1,ns
  idx=p%wave_send(partner+1)%idx(j)
  if(p%use_halo)then
   x=modulo(idx-1,p%halo_shape(1));y=modulo((idx-1)/p%halo_shape(1),p%halo_shape(2));z=(idx-1)/(p%halo_shape(1)*p%halo_shape(2))
   sb(j,:)=halo(x,y,z,:nb)
  else
   sb(j,:)=values(idx,:nb)
  endif
 enddo
 if(partner==rank)then
  rb=sb
 else
  call MPI_Sendrecv(sb,size(sb),MPI_DOUBLE_COMPLEX,partner,8,rb,size(rb),MPI_DOUBLE_COMPLEX,partner,8,comm,MPI_STATUS_IGNORE,e)
  if(ns>0)then
   p%messages=p%messages+1;p%bytes=p%bytes+16_int64*size(sb,kind=int64)
  endif
 endif
 do j=1,nr
  idx=p%wave_recv(partner+1)%idx(j)
  tile(idx,:nb)=rb(j,:nb);filled(idx)=.true.
 enddo
 deallocate(sb,rb)
enddo
bad=merge(0,1,all(filled));status=bad;call comm_get_max(status,p%world)
#else
tile=0d0;status=1
#endif
end subroutine
! Transfer one bounded wavefunction batch into the distributed FFT input layout.
! Identical coordinate/copy selection as density delivery; no convolution here.
subroutine stage_waves(p,values,tile,status,packed)
  implicit none
type(group_plan),intent(inout)::p
complex(8),intent(in)::values(:,:)
complex(8),intent(out)::tile(:,:)
integer,intent(out)::status
logical,intent(in)::packed
#ifdef USE_MPI
integer::batch
integer::nb,k,g,z,y,x,hc(3),ii,gr,j,partner,round,e,q(3),rounds,bad
complex(8),allocatable::send(:,:),recv(:,:),halo_data(:,:,:,:),halo_chunk(:,:,:,:)
logical,allocatable::filled(:)
batch=size(tile,2)
nb=size(values,2);status=0;tile=0
allocate(send(product(p%m),batch),recv(product(p%m),batch));send=0
send(:,:nb)=values
rounds=1
do while(rounds<p%peers);rounds=rounds*2;enddo
if(packed)then
 call prepare_wave_maps(p)
endif
if(p%use_halo)then
 allocate(halo_data(0:p%halo_shape(1)-1,0:p%halo_shape(2)-1,0:p%halo_shape(3)-1,batch), &
 filled(product(p%fm)))
 if(.not.packed)allocate(halo_chunk(0:p%halo_shape(1)-1,0:p%halo_shape(2)-1,0:p%halo_shape(3)-1,batch))
endif
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
  if(packed)then
   call packed_delivery(p,values,tile,status,halo_data)
   return
  endif
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
 if(packed)then
  call packed_delivery(p,values,tile,status)
  return
 endif
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
#else
tile=0d0;status=1
#endif
end subroutine
subroutine wave_apply(p,sources,targets,action,status,scratch_bytes,aggregate_return,return_bytes,packed_transfer,batch_size)
  implicit none
type(group_plan),intent(inout)::p
complex(8),intent(in)::sources(:,:),targets(:,:)
complex(8),intent(out)::action(:,:)
integer,intent(out)::status
integer(int64),intent(out),optional::scratch_bytes,return_bytes
logical,intent(in),optional::aggregate_return,packed_transfer
integer,intent(in),optional::batch_size
#ifdef USE_MPI
logical::aggregate,packed
integer::flag,fmin,fmax
integer(int64)::returned
integer::batch
integer::ns,nt,bad,e,hi,lo,first,nb,i,k,g,x,y,z,q(3),j,round,partner
complex(8),allocatable::sw(:,:),tw(:,:),density(:,:),spec(:,:),back(:,:),chunk(:,:),accum(:,:)
real(8)::t,beg
status=1;action=0;bad=0;returned=0;aggregate=.true.
batch=8
if(present(batch_size))batch=batch_size
if(batch/=8.and.batch/=16)bad=1
call MPI_Allreduce(batch,fmin,1,MPI_INTEGER,MPI_MIN,p%world,e)
fmax=batch;call comm_get_max(fmax,p%world)
if(fmin/=fmax)bad=1
if(present(aggregate_return))aggregate=aggregate_return
packed=.true.
if(present(packed_transfer))packed=packed_transfer
flag=merge(1,0,packed)
call MPI_Allreduce(flag,fmin,1,MPI_INTEGER,MPI_MIN,p%world,e)
fmax=flag;call comm_get_max(fmax,p%world)
if(fmin/=fmax)bad=1
if(present(return_bytes))return_bytes=0
flag=merge(1,0,aggregate)
call MPI_Allreduce(flag,fmin,1,MPI_INTEGER,MPI_MIN,p%world,e)
fmax=flag;call comm_get_max(fmax,p%world)
if(fmin/=fmax)bad=1
if(present(scratch_bytes))scratch_bytes=0
if(.not.p%ready)bad=1
ns=size(sources,2);nt=size(targets,2)
if(ns<0.or.ns>batch.or.nt<0)bad=1
if(size(sources,1)/=product(p%m).or.size(targets,1)/=product(p%m))bad=1
if(any(shape(action)/=shape(targets)))bad=1
if(.not.sr_all_finite(sources))bad=1
if(.not.sr_all_finite(targets))bad=1
status=bad;call comm_get_max(status,p%world)
if(status/=0)return
hi=ns;call comm_get_max(hi,p%world)
call MPI_Allreduce(ns,lo,1,MPI_INTEGER,MPI_MIN,p%world,e)
if(hi/=lo)bad=1
hi=nt;call comm_get_max(hi,p%world)
call MPI_Allreduce(nt,lo,1,MPI_INTEGER,MPI_MIN,p%world,e)
if(hi/=lo)bad=1
status=bad;call comm_get_max(status,p%world)
if(status/=0)return
if(ns==0.or.nt==0)then
 status=0;return
endif
allocate(sw(product(p%fm),batch),tw(product(p%fm),batch),density(product(p%fm),batch), &
 spec(product(p%fm),batch),back(product(p%fm),batch),chunk(product(p%fm),batch))
if(aggregate)allocate(accum(product(p%fm),batch))
! Peak explicit complex scratch excludes FFT library/cache and native inputs.
if(packed)call prepare_wave_maps(p)
if(present(scratch_bytes))then
 scratch_bytes=16_int64*batch*(6_int64*product(int(p%fm,int64))+2_int64*product(int(p%m,int64)))
 if(packed)scratch_bytes=scratch_bytes+p%wave_map_bytes
 if(aggregate)scratch_bytes=scratch_bytes+16_int64*batch*product(int(p%fm,int64))
 if(packed)then
  scratch_bytes=scratch_bytes+16_int64*batch*(maxval([(size(p%wave_send(i)%idx),i=1,size(p%wave_send))])+ &
   maxval([(size(p%wave_recv(i)%idx),i=1,size(p%wave_recv))]))
  if(p%use_halo)scratch_bytes=scratch_bytes+16_int64*batch*product(int(p%halo_shape,int64))
 else
  if(p%use_halo)scratch_bytes=scratch_bytes+32_int64*batch*product(int(p%halo_shape,int64))
 endif
endif
beg=sr_walltime();t=beg
call stage_waves(p,sources,sw,status,packed)
if(status/=0)return
p%phase(1)=p%phase(1)+sr_walltime()-t
! Keep these source tiles while visiting all target blocks. Pair products are local.
do first=1,nt,batch
 nb=min(batch,nt-first+1);t=sr_walltime()
 if(aggregate)accum=0
 call stage_waves(p,targets(:,first:first+nb-1),tw,status,packed)
 if(status/=0)return
 p%phase(1)=p%phase(1)+sr_walltime()-t
 do i=1,ns
  density=0
  do k=1,nb;density(:,k)=conjg(sw(:,i))*tw(:,k);enddo
  t=sr_walltime()
  call mesh_transform(p%length,p%gd,p%gc,p%axis,density,spec,-1,bad)
  status=bad;call comm_get_max(status,p%world)
  if(status/=0)return
  p%phase(2)=p%phase(2)+sr_walltime()-t;t=sr_walltime()
  do k=1,nb;spec(:,k)=spec(:,k)*p%filter(:,1);enddo
  p%phase(3)=p%phase(3)+sr_walltime()-t;t=sr_walltime()
  call mesh_transform(p%length,p%gd,p%gc,p%axis,spec,back,1,bad)
  status=bad;call comm_get_max(status,p%world)
  if(status/=0)return
  p%phase(2)=p%phase(2)+sr_walltime()-t;t=sr_walltime()
  if(aggregate)then
   do k=1,nb;accum(:,k)=accum(:,k)-sw(:,i)*back(:,k);enddo
  else
   call return_block(back,i,.false.)
  endif
  p%phase(4)=p%phase(4)+sr_walltime()-t
 enddo
 if(aggregate)then
  t=sr_walltime()
  call return_block(accum,0,.true.)
  p%phase(4)=p%phase(4)+sr_walltime()-t
 endif
enddo
if(present(return_bytes))return_bytes=returned
p%phase(5)=p%phase(5)+sr_walltime()-beg;status=0
contains
subroutine return_block(values,source_index,weighted)
  implicit none
complex(8),intent(in)::values(:,:)
integer,intent(in)::source_index
logical,intent(in)::weighted
  if(all(p%fm==p%m).and.all(p%flo==p%lo).and.all(p%origin==0))then
   do k=1,nb
    if(weighted)then
     action(:,first+k-1)=action(:,first+k-1)+values(:,k)
    else
     action(:,first+k-1)=action(:,first+k-1)-sources(:,source_index)*values(:,k)
    endif
   enddo
  else
   do round=0,p%gpeers-1
    partner=ieor(p%grank,round)
    if(partner==p%grank)then
     chunk=values
    else
     call MPI_Sendrecv(values,size(values),MPI_DOUBLE_COMPLEX,partner,1,chunk,size(chunk), &
      MPI_DOUBLE_COMPLEX,partner,1,p%group,MPI_STATUS_IGNORE,e)
     returned=returned+16_int64*size(values,kind=int64)
     p%messages=p%messages+1;p%bytes=p%bytes+16_int64*size(values,kind=int64)
    endif
    g=0
    do z=0,p%m(3)-1;do y=0,p%m(2)-1;do x=0,p%m(1)-1
     g=g+1;q=modulo(p%lo+[x,y,z]-p%origin,p%n)-p%gcoords(:,partner+1)*p%fm
     if(any(q<0).or.any(q>=p%fm))cycle
     j=1+q(1)+p%fm(1)*(q(2)+p%fm(2)*q(3))
     if(weighted)then
    action(g,first:first+nb-1)=action(g,first:first+nb-1)+chunk(j,:nb)
   else
    action(g,first:first+nb-1)=action(g,first:first+nb-1)-sources(g,source_index)*chunk(j,:nb)
   endif
    enddo;enddo;enddo
   enddo
  endif
end subroutine
#else
action=0d0;status=1
if(present(scratch_bytes))scratch_bytes=0
if(present(return_bytes))return_bytes=0
#endif
end subroutine
! Same-source ACE only: one orientation of each unordered global orbital pair.
! The endpoint selected by parity balances the two orbital owners. Caller must
! reduce reverse_action to the source owner; no replicated global mesh is held.
subroutine wave_apply_pairs(p,sources,targets,source_ids,target_ids,action,reverse_action,status,pairs,scratch_bytes,batch)
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
type(group_plan),intent(inout)::p
complex(8),intent(in)::sources(:,:),targets(:,:)
integer,intent(in)::source_ids(:),target_ids(:),batch
complex(8),intent(out)::action(:,:),reverse_action(:,:)
integer,intent(out)::status
integer(int64),intent(out)::pairs,scratch_bytes
#ifdef USE_MPI
integer::ns,nt,bad,e,hi,lo,first,nb,i,k,j,nsel,selected(16),chosen,g,x,y,z,q(3),idx
integer,allocatable::imin(:),imax(:)
complex(8),allocatable::sw(:,:),tw(:,:),density(:,:),spec(:,:),back(:,:),accum(:,:),raccum(:,:)
real(8)::t,beg,scale
status=1;bad=0;action=0;reverse_action=0;pairs=0;scratch_bytes=0
ns=size(sources,2);nt=size(targets,2)
if(batch/=8.and.batch/=16)bad=1
if(ns>batch.or.size(source_ids)/=ns.or.size(target_ids)/=nt)bad=1
if(.not.p%ready.or.p%gpeers/=1)bad=1
if(size(sources,1)/=product(p%m).or.size(targets,1)/=product(p%m))bad=1
if(any(shape(action)/=shape(targets)).or.any(shape(reverse_action)/=shape(sources)))bad=1
if(.not.sr_all_finite(sources))bad=1
if(.not.sr_all_finite(targets))bad=1
if(any(source_ids<1).or.any(target_ids<1))bad=1
if(ns>1)then
 if(any(source_ids(2:)<=source_ids(:ns-1)))bad=1
endif
if(nt>1)then
 if(any(target_ids(2:)<=target_ids(:nt-1)))bad=1
endif
status=bad;call comm_get_max(status,p%world)
if(status/=0)return
hi=batch;call comm_get_max(hi,p%world)
call MPI_Allreduce(batch,lo,1,MPI_INTEGER,MPI_MIN,p%world,e)
if(hi/=lo)bad=1
hi=ns;call comm_get_max(hi,p%world)
call MPI_Allreduce(ns,lo,1,MPI_INTEGER,MPI_MIN,p%world,e)
if(hi/=lo)bad=1
hi=nt;call comm_get_max(hi,p%world)
call MPI_Allreduce(nt,lo,1,MPI_INTEGER,MPI_MIN,p%world,e)
if(hi/=lo)bad=1
status=bad;call comm_get_max(status,p%world)
if(status/=0)return
allocate(imin(ns+nt),imax(ns+nt))
call MPI_Allreduce([source_ids,target_ids],imin,ns+nt,MPI_INTEGER,MPI_MIN,p%world,e)
call MPI_Allreduce([source_ids,target_ids],imax,ns+nt,MPI_INTEGER,MPI_MAX,p%world,e)
if(any(imin/=imax))bad=1
! A real even spatial kernel has a real Fourier multiplier. Check rather than
! assume this symmetry for arbitrary caller-provided complex kernels.
scale=max(1d0,maxval(abs(p%filter)))
do g=1,size(p%filter,1)
 if(.not.ieee_is_finite(real(p%filter(g,1))).or..not.ieee_is_finite(aimag(p%filter(g,1))))bad=1
enddo
if(maxval(abs(aimag(p%filter)))>128d0*epsilon(1d0)*scale)bad=1
if(any(p%fm/=p%length).or.any(p%flo/=0))bad=1
if(bad==0)then
 do z=0,p%fm(3)-1;do y=0,p%fm(2)-1;do x=0,p%fm(1)-1
  idx=1+x+p%fm(1)*(y+p%fm(2)*z)
  q=modulo(-[x,y,z],p%length)
  j=1+q(1)+p%fm(1)*(q(2)+p%fm(2)*q(3))
  if(abs(p%filter(idx,1)-conjg(p%filter(j,1)))>128d0*epsilon(1d0)*scale)bad=1
 enddo;enddo;enddo
endif
status=bad;call comm_get_max(status,p%world)
if(status/=0)return
if(ns==0.or.nt==0)then
 status=0;return
endif
allocate(sw(product(p%fm),batch),tw(product(p%fm),batch),density(product(p%fm),batch), &
 spec(product(p%fm),batch),back(product(p%fm),batch),accum(product(p%fm),batch),raccum(product(p%fm),batch))
raccum=0
call prepare_wave_maps(p)
scratch_bytes=16_int64*batch*(7_int64*product(int(p%fm,int64))+2_int64*product(int(p%m,int64)))+p%wave_map_bytes
scratch_bytes=scratch_bytes+16_int64*batch*(maxval([(size(p%wave_send(i)%idx),i=1,size(p%wave_send))])+ &
 maxval([(size(p%wave_recv(i)%idx),i=1,size(p%wave_recv))]))
if(p%use_halo)scratch_bytes=scratch_bytes+16_int64*batch*product(int(p%halo_shape,int64))
beg=sr_walltime();t=beg
call stage_waves(p,sources,sw,status,.true.)
if(status/=0)return
p%phase(1)=p%phase(1)+sr_walltime()-t
do first=1,nt,batch
 nb=min(batch,nt-first+1);accum=0;t=sr_walltime()
 call stage_waves(p,targets(:,first:first+nb-1),tw,status,.true.)
 if(status/=0)return
 p%phase(1)=p%phase(1)+sr_walltime()-t
 do i=1,ns
  nsel=0
  do k=1,nb
   j=first+k-1
   chosen=min(source_ids(i),target_ids(j))
   if(modulo(source_ids(i)+target_ids(j),2)==1)chosen=max(source_ids(i),target_ids(j))
   if(target_ids(j)/=chosen)cycle
   nsel=nsel+1;selected(nsel)=k
  enddo
  if(nsel==0)cycle
  pairs=pairs+int(nsel,int64)
  do k=1,nsel;density(:,k)=conjg(sw(:,i))*tw(:,selected(k));enddo
  t=sr_walltime()
  call mesh_transform(p%length,p%gd,p%gc,p%axis,density(:,:nsel),spec(:,:nsel),-1,bad)
  status=bad;call comm_get_max(status,p%world)
  if(status/=0)return
  p%phase(2)=p%phase(2)+sr_walltime()-t;t=sr_walltime()
  do k=1,nsel;spec(:,k)=spec(:,k)*p%filter(:,1);enddo
  p%phase(3)=p%phase(3)+sr_walltime()-t;t=sr_walltime()
  call mesh_transform(p%length,p%gd,p%gc,p%axis,spec(:,:nsel),back(:,:nsel),1,bad)
  status=bad;call comm_get_max(status,p%world)
  if(status/=0)return
  p%phase(2)=p%phase(2)+sr_walltime()-t;t=sr_walltime()
  do k=1,nsel
   j=selected(k)
   accum(:,j)=accum(:,j)-sw(:,i)*back(:,k)
   if(source_ids(i)/=target_ids(first+j-1))raccum(:,i)=raccum(:,i)-tw(:,j)*conjg(back(:,k))
  enddo
  p%phase(4)=p%phase(4)+sr_walltime()-t
 enddo
 t=sr_walltime()
 call crop(accum,action(:,first:first+nb-1))
 p%phase(4)=p%phase(4)+sr_walltime()-t
enddo
t=sr_walltime();call crop(raccum,reverse_action)
p%phase(4)=p%phase(4)+sr_walltime()-t
p%phase(5)=p%phase(5)+sr_walltime()-beg
bad=0
if(.not.sr_all_finite(action))bad=1
if(.not.sr_all_finite(reverse_action))bad=1
status=bad;call comm_get_max(status,p%world)
contains
subroutine crop(values,out)
  implicit none
complex(8),intent(in)::values(:,:)
complex(8),intent(out)::out(:,:)
integer::g,x,y,z,q(3),idx
out=0;g=0
do z=0,p%m(3)-1;do y=0,p%m(2)-1;do x=0,p%m(1)-1
 g=g+1;q=modulo(p%lo+[x,y,z]-p%origin,p%n)-p%flo
 if(any(q<0).or.any(q>=p%fm))cycle
 idx=1+q(1)+p%fm(1)*(q(2)+p%fm(2)*q(3))
 out(g,:)=values(idx,:size(out,2))
enddo;enddo;enddo
end subroutine
#else
action=0d0;reverse_action=0d0;status=1;pairs=0;scratch_bytes=0
#endif
end subroutine

end module
