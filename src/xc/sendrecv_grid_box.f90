! Periodic box-halo extension built on SALMON's unmodified persistent face
! exchange. Independent plans and communicator keep legacy stencils isolated.
module sendrecv_grid_box
use structures,only:s_rgrid,s_sendrecv_grid
use sendrecv_grid,only:init_sendrecv_grid,update_overlap,dealloc_cache
use communication,only:comm_get_groupinfo,comm_get_max,comm_create_group,comm_free_group,comm_proc_null
implicit none
private
public::s_box_halo,init_box_halo,update_box_halo,destroy_box_halo
type s_box_halo
 type(s_sendrecv_grid),allocatable::step(:)
 integer::comm,nb
 integer::is_array(3),ie_array(3)
 logical::ready=.false.
end type
contains
subroutine init_box_halo(p,rg,width,nb,comm,neig,status)
type(s_box_halo),intent(out)::p
type(s_rgrid),intent(in)::rg
integer,intent(in)::width(3),nb,comm,neig(2,3)
integer,intent(out)::status
integer::bad,r,n,m(3),a,d,w,k,count,nn(2,3),v(10),hi,neg,i
bad=0;m=rg%ie-rg%is+1
call comm_get_groupinfo(comm,r,n)
if(any(width<0).or.any(m<=0).or.nb<=0)bad=1
if(any(rg%is_array>rg%is-width).or.any(rg%ie_array<rg%ie+width))bad=1
do a=1,3
 if(width(a)>0.and.any(neig(:,a)==comm_proc_null))bad=1
 if(any((neig(:,a)<0.or.neig(:,a)>=n).and.neig(:,a)/=comm_proc_null))bad=1
enddo
call comm_get_max(bad,comm);status=bad
if(status/=0)return
! All ranks must use the same decomposition size, batch and widths so every
! persistent exchange has an identical count. Reject collectively before split.
v=[m,width,nb,rg%ie_array-rg%is_array+1]
do i=1,size(v)
 hi=v(i);neg=-v(i);call comm_get_max(hi,comm);call comm_get_max(neg,comm)
 if(hi/=-neg)bad=1
enddo
call comm_get_max(bad,comm);status=bad
if(status/=0)return
p%comm=comm_create_group(comm,0,r);p%nb=nb
p%is_array=rg%is_array;p%ie_array=rg%ie_array
count=sum((width+m-1)/m);allocate(p%step(count));k=0
do a=1,3
 d=1
 do while(d<=width(a))
  w=min(m(a),width(a)-d+1);k=k+1
  nn=comm_proc_null;nn(:,a)=neig(:,a)
  call init_sendrecv_grid(p%step(k),rg,nb,p%comm,nn)
  ! Disable unused directions' cache extents; only the current axis communicates.
  p%step(k)%is_block=0;p%step(k)%ie_block=-1
  do i=1,3
   if(i==a)cycle
   p%step(k)%is_block(i,:,:,a)=rg%is(i)
   p%step(k)%ie_block(i,:,:,a)=rg%ie(i)
   if(i<a)then
    p%step(k)%is_block(i,:,:,a)=rg%is(i)-width(i)
    p%step(k)%ie_block(i,:,:,a)=rg%ie(i)+width(i)
   endif
  enddo
  ! Later slabs can source already received opposite halos. w<=local width
  ! guarantees their dependencies were received in an earlier completed slab.
  p%step(k)%is_block(a,1,1,a)=rg%ie(a)-d-w+2
  p%step(k)%ie_block(a,1,1,a)=rg%ie(a)-d+1
  p%step(k)%is_block(a,2,1,a)=rg%ie(a)+d
  p%step(k)%ie_block(a,2,1,a)=rg%ie(a)+d+w-1
  p%step(k)%is_block(a,1,2,a)=rg%is(a)+d-1
  p%step(k)%ie_block(a,1,2,a)=rg%is(a)+d+w-2
  p%step(k)%is_block(a,2,2,a)=rg%is(a)-d-w+1
  p%step(k)%ie_block(a,2,2,a)=rg%is(a)-d
  d=d+w
 enddo
enddo
p%ready=.true.;status=0
end subroutine
subroutine update_box_halo(p,rg,data)
type(s_box_halo),intent(inout)::p
type(s_rgrid),intent(in)::rg
complex(8),intent(inout)::data(rg%is_array(1):rg%ie_array(1),rg%is_array(2):rg%ie_array(2), &
 rg%is_array(3):rg%ie_array(3),p%nb)
integer::k
! Caller retains the grid bounds and batch size used to prepare the plan.
if(.not.p%ready)error stop 'box halo plan not ready'
if(any(rg%is_array/=p%is_array).or.any(rg%ie_array/=p%ie_array))error stop 'box halo grid changed'
do k=1,size(p%step)
 call update_overlap(p%step(k),rg,data)
enddo
end subroutine
subroutine destroy_box_halo(p)
type(s_box_halo),intent(inout)::p
integer::k
if(.not.p%ready)return
do k=1,size(p%step);call dealloc_cache(p%step(k));enddo
deallocate(p%step);call comm_free_group(p%comm);p%ready=.false.
end subroutine
end module
