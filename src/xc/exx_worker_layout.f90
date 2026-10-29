! Integer work model for the existing blocked-cyclic exchange worker loops.
module exx_worker_layout
 use iso_fortran_env,only:int64
 implicit none
 private
 public::choose_worker_layout,choose_balanced_strips
 integer(int64),parameter::density_cap=67108864_int64
contains
 subroutine choose_worker_layout(nfirst,ng,nk,nc,np,requested_groups,requested_divisor, &
  row_groups,column_groups,b,chunk,peak_pairs,active_ranks,ierr)
 implicit none
 integer,intent(in)::nfirst,ng,nk,nc,np,requested_groups,requested_divisor
 integer,intent(out)::row_groups,column_groups,b,chunk,ierr
 integer(int64),intent(out)::peak_pairs,active_ranks
 integer::factor
 integer(int64)::pair_bytes
 logical::found
 ierr=1;row_groups=0;column_groups=0;b=0;chunk=0;peak_pairs=0;active_ranks=0
 if(min(nfirst,ng,nk,np)<1.or.(nc/=1.and.nc/=2))return
 if(requested_groups<0.or.requested_divisor<0)return
 pair_bytes=int(nk,int64)*16_int64*nc*nc
 if(pair_bytes>density_cap)return
 found=.false.
 if(requested_groups>0)then
  if(modulo(np,requested_groups)/=0)return
  call consider(requested_groups)
 else
  factor=1
  do while(factor<=np/factor)
   if(modulo(np,factor)==0)then
    call consider(factor)
    if(factor/=np/factor)call consider(np/factor)
   endif
   factor=factor+1
  enddo
 endif
 if(found)ierr=0
 contains
 subroutine consider(columns)
 implicit none
 integer,intent(in)::columns
 integer::rows,rb,cb,divisor
 integer(int64)::work,active
 rows=np/columns
 rb=min(16,1+(nfirst-1)/rows)
 rb=min(rb,int(density_cap/pair_bytes))
 divisor=requested_divisor
 if(divisor==0)divisor=columns
 cb=min(1+(ng-1)/divisor,int(density_cap/(int(rb,int64)*pair_bytes)))
 work=axis_max(nfirst,rows,rb)*axis_max(ng,columns,cb)
 active=axis_active(nfirst,rows,rb)*axis_active(ng,columns,cb)
 if(found)then
  if(work>peak_pairs)return
  if(work==peak_pairs)then
   if(active<active_ranks)return
   if(active==active_ranks.and.columns>=column_groups)return
  endif
 endif
 found=.true.;row_groups=rows;column_groups=columns;b=rb;chunk=cb
 peak_pairs=work;active_ranks=active
 end subroutine consider
 end subroutine choose_worker_layout

 subroutine choose_balanced_strips(nfirst,ng,nk,nc,np,row_index,col_start,col_count, &
  group_id,leaders,groups,chunk,peak_pairs,active_ranks,ierr)
 implicit none
 integer,intent(in)::nfirst,ng,nk,nc,np
 integer,intent(out)::row_index(:),col_start(:),col_count(:),group_id(:),leaders(:),groups,chunk,ierr
 integer(int64),intent(out)::peak_pairs,active_ranks
 integer::q,extra,r,fanout,j,k,offset,idbase,small,large
 integer(int64)::pair_bytes
 ierr=1;groups=0;chunk=0;peak_pairs=0;active_ranks=0
 row_index=0;col_start=0;col_count=0;group_id=0;leaders=0
 if(min(nfirst,ng,nk,np)<1.or.(nc/=1.and.nc/=2))return
 if(min(size(row_index),size(col_start),size(col_count),size(group_id),size(leaders))<np)return
 if(np<nfirst)return
 q=np/nfirst;extra=modulo(np,nfirst)
 if(q+merge(1,0,extra>0)>ng)return
 pair_bytes=int(nk,int64)*16_int64*nc*nc
 offset=0
 do r=1,nfirst
  fanout=q+merge(1,0,r<=extra)
  idbase=0
  if(r>extra.and.extra>0)idbase=q+1
  do j=0,fanout-1
   k=offset+j+1
   row_index(k)=r
   ! Quotient/remainder partition without default-integer products of grid/rank counts.
   col_count(k)=ng/fanout+merge(1,0,j<modulo(ng,fanout))
   col_start(k)=1+j*(ng/fanout)+min(j,modulo(ng,fanout))
   group_id(k)=idbase+j
   if(r==1.or.r==extra+1)leaders(idbase+j+1)=k-1
  enddo
  offset=offset+fanout
 enddo
 groups=q
 if(extra>0)groups=2*q+1
 if(extra==nfirst)groups=q+1
 small=minval(col_count(:np));large=maxval(col_count(:np))
 if(small<1.or.large-small>1)return
 if(int(large,int64)>density_cap/pair_bytes)return
 chunk=large;peak_pairs=int(large,int64);active_ranks=int(np,int64);ierr=0
 end subroutine choose_balanced_strips

 pure function axis_max(points,groups,block) result(count)
 implicit none
 integer,intent(in)::points,groups,block
 integer(int64)::count,period
 period=int(groups,int64)*block
 count=(int(points,int64)/period)*block+min(int(block,int64),modulo(int(points,int64),period))
 end function axis_max

 pure function axis_active(points,groups,block) result(count)
 implicit none
 integer,intent(in)::points,groups,block
 integer(int64)::count
 count=min(int(groups,int64),1_int64+(int(points,int64)-1)/block)
 end function axis_active
end module exx_worker_layout
