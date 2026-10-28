! Integer work model for the existing blocked-cyclic exchange worker loops.
module exx_worker_layout
 use iso_fortran_env,only:int64
 implicit none
 private
 public::choose_worker_layout
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
