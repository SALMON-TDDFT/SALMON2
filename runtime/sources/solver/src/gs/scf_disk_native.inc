! SCF rank-local disk v1. Host authenticates binary/input/topology and all files.
integer(c_int) function sp_gs_scf_disk_protocol() bind(C)
 sp_gs_scf_disk_protocol=1
end function
integer(c_int) function sp_gs_scf_bootstrap(path,token) bind(C)
 character(c_char),intent(in)::path(*)
 integer(c_int64_t),intent(out)::token
 sp_gs_scf_bootstrap=create_gs(path,token,.false.,0)
end function
integer(c_int) function sp_gs_scf_disk(token,path,mode) bind(C) result(status)
 integer(c_int64_t),value::token
 character(c_char),intent(in)::path(*)
 integer(c_int),value::mode
 type(s_scf_payload),allocatable::value
 character(4096)::filename
 integer::index,i,group,rank,nproc,u,ios,close_status,header(5),expected(5),counters(2),reference
 integer(8)::limit,position
 logical::reading
 call gs_lookup(token,index,status)
 if(status/=0)return
 call comm_get_globalinfo(group,rank,nproc)
 reference=mode;call comm_bcast(reference,group)
 if(mode/=reference.or.(mode/=0.and.mode/=2))status=11
 if(.not.scf_context_ok(gs_handles(index)%context))status=15
 if(allocated(gs_handles(index)%context%scf_pending))status=11
 filename=''
 do i=1,len(filename)
  if(path(i)==c_null_char)exit
  filename(i:i)=path(i)
 enddo
 if(i>len(filename).or.len_trim(filename)==0)status=11
 call comm_get_max(status,group)
 if(status/=0)return
 reading=mode==2;ios=0;u=0
 allocate(value)
 associate(ctx=>gs_handles(index)%context)
 if(.not.reading)call scf_copy(ctx,value,.false.)
 expected=[1,rank,nproc,storage_size(0),storage_size(0d0)]
 header=expected;counters=[ctx%initializations,ctx%executed_iterations]
 if(reading)then
  open(newunit=u,file=trim(filename),status='old',action='read',access='stream',form='unformatted',iostat=ios)
 else
  open(newunit=u,file=trim(filename),status='new',action='write',access='stream',form='unformatted',iostat=ios)
 endif
 if(ios==0)then
  if(reading)then
   inquire(unit=u,size=limit)
   read(u,iostat=ios)header,counters
   if(any(header/=expected).or.any(counters<0))ios=1
  else
   limit=huge(0_8)
   write(u,iostat=ios)header,counters
  endif
  if(ios==0)call scf_io_s_scf_payload(u,reading,value,ios,limit)
  if(reading.and.ios==0)then
   inquire(unit=u,pos=position)
   if(position/=limit+1_8)ios=1
  endif
  close(u,iostat=close_status)
  if(close_status/=0)ios=close_status
 endif
 if(ios/=0)status=4
 call comm_get_max(status,group)
 if(status==0.and.reading)then
  ! Rank-local ACE metric communicators are runtime resources, never disk data.
#ifdef USE_HSE
  value%hse%ace%metric_comm=ctx%info%icomm_r
  value%hse%initial_ace%metric_comm=ctx%info%icomm_r
  value%hse%midpoint_ace%metric_comm=ctx%info%icomm_r
  value%hse%scf_damped_ace%metric_comm=ctx%info%icomm_r
#endif
  call move_alloc(value,ctx%scf_pending)
 endif
 end associate
end function

! Opt-in ordinary-XC full-real-grid k-only migration. Numerical data never
! pass through Python. Existing exact-layout disk/candidate entry is unchanged.
integer(c_int) function sp_gs_scf_k_repartition_protocol() bind(C)
 sp_gs_scf_k_repartition_protocol=1
end function
integer(c_int) function sp_gs_scf_k_repartition(token,path,old_size,kbounds) bind(C) result(status)
 integer(c_int64_t),value::token
 character(c_char),intent(in)::path(*)
 integer(c_int),value::old_size
 integer(c_int),intent(in)::kbounds(2,*)
 type(s_scf_payload),allocatable::value,donor
 character(4096)::directory,filename,expected_directory
 integer::index,group,rank,nproc,i,j,u,ios,cs,header(5),counters(2),ks,ke,lo,hi,iteration,old_reference,pair(2)
 integer(8)::limit,position
 logical::first
 logical,allocatable::covered(:)
 call gs_lookup(token,index,status)
 if(status/=0)return
 call comm_get_globalinfo(group,rank,nproc)
 associate(ctx=>gs_handles(index)%context)
 if(.not.scf_context_ok(ctx))status=15
 if(allocated(ctx%scf_pending))status=11
 if(ctx%info%isize_r/=1.or.ctx%info%isize_o/=1.or.ctx%info%isize_k/=nproc)status=15
 select case(trim(ctx%system%xc_input%xcname))
 case('pz','libxc_pz','libxc_pbe')
 case default
  status=15
 end select
 if(old_size<1.or.old_size>ctx%system%nk)status=11
 if(nproc>ctx%system%nk.or.ctx%info%numk<1)status=15
 directory=''
 do i=1,len(directory)
  if(path(i)==c_null_char)exit
  directory(i:i)=path(i)
 enddo
 if(i>len(directory).or.len_trim(directory)==0)status=11
 old_reference=old_size;call comm_bcast(old_reference,group)
 expected_directory=directory;call comm_bcast(expected_directory,group)
 if(old_size/=old_reference.or.directory/=expected_directory)status=11
 call comm_get_max(status,group)
 if(status/=0)return
 do j=1,old_size
  pair=kbounds(:,j);call comm_bcast(pair,group)
  if(any(pair/=kbounds(:,j)))status=11
 enddo
 call comm_get_max(status,group)
 if(status/=0)return
 ks=ctx%info%ik_s;ke=ctx%info%ik_e
 allocate(covered(ks:ke),stat=ios)
 if(ios/=0)status=4
 call comm_get_max(status,group)
 if(status/=0)return
 covered=.false.
 first=.true.;iteration=-1
 do j=1,old_size
  lo=kbounds(1,j);hi=kbounds(2,j)
  if(lo<1.or.hi<lo.or.hi>ctx%system%nk)then
   status=11;exit
  endif
  if(hi<ks.or.lo>ke)cycle
  if(any(covered(max(lo,ks):min(hi,ke))))then
   status=11;exit
  endif
  allocate(donor,stat=ios)
  if(ios/=0)then
   status=4;exit
  endif
  write(filename,'(a,a,i8.8,a)')trim(directory),'/rank-',j-1,'.bin'
  open(newunit=u,file=trim(filename),status='old',action='read',access='stream',form='unformatted',iostat=ios)
  if(ios==0)then
   inquire(unit=u,size=limit)
   read(u,iostat=ios)header,counters
   if(ios==0)then
    if(any(header/=[1,j-1,old_size,storage_size(0),storage_size(0d0)]).or.any(counters<0))ios=1
   endif
   if(ios==0)call scf_io_s_scf_payload(u,.true.,donor,ios,limit)
   if(ios==0)then
    inquire(unit=u,pos=position)
    if(position/=limit+1_8)ios=1
   endif
   close(u,iostat=cs)
   if(cs/=0)ios=cs
  endif
  if(ios/=0)then
   status=4;exit
  endif
  if(first)then
   allocate(value,stat=ios)
   if(ios/=0)then
    status=4;exit
   endif
   value=donor;iteration=donor%iterations
  endif
  if(donor%iterations/=iteration.or.donor%system%nk/=ctx%system%nk)then
   status=11;exit
  endif
  if(.not.allocated(donor%spsi%zwf))then
   status=15;exit
  endif
  if(lbound(donor%spsi%zwf,6)/=lo.or.ubound(donor%spsi%zwf,6)/=hi)then
   status=11;exit
  endif
  call scf_merge_k_orbital(value%spsi,donor%spsi,first,ks,ke,status)
  call scf_merge_k_orbital(value%shpsi,donor%shpsi,first,ks,ke,status)
  call scf_merge_k_orbital(value%sttpsi,donor%sttpsi,first,ks,ke,status)
  call scf_merge_k_orbital(value%cg%pre_gk,donor%cg%pre_gk,first,ks,ke,status)
  call scf_merge_k_orbital(value%cg%pk,donor%cg%pk,first,ks,ke,status)
  call scf_merge_k_orbital(value%cg%hpk,donor%cg%hpk,first,ks,ke,status)
  if(status/=0)exit
  covered(max(lo,ks):min(hi,ke))=.true.;first=.false.
  deallocate(donor)
 enddo
 if(status==0.and.(.not.all(covered).or..not.allocated(value)))status=11
 call comm_get_max(status,group)
 if(status/=0)return
 call value%continuation_inputs%parallel_from(ctx%inputs)
 call move_alloc(value,ctx%scf_pending)
 end associate
end function
subroutine scf_merge_k_orbital(target,source,first,ks,ke,status)
 type(s_orbital),intent(inout)::target
 type(s_orbital),intent(in)::source
 logical,intent(in)::first
 integer,intent(in)::ks,ke
 integer,intent(inout)::status
 integer::lb(7),ub(7),a,b,istat
 if(status/=0)return
 if(first)target%update_zwf_overlap=source%update_zwf_overlap
 if(first)then
  if(allocated(target%rwf))deallocate(target%rwf)
  if(allocated(source%rwf))then
   lb=lbound(source%rwf);ub=ubound(source%rwf);lb(6)=ks;ub(6)=ke
   allocate(target%rwf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4), &
    lb(5):ub(5),lb(6):ub(6),lb(7):ub(7)),stat=istat)
   if(istat/=0)then
    status=4;return
   endif
   target%rwf=0d0
  endif
 endif
 if(allocated(target%rwf).neqv.allocated(source%rwf))then
  status=11;return
 endif
 if(allocated(source%rwf))then
  lb=lbound(source%rwf);ub=ubound(source%rwf);lb(6)=ks;ub(6)=ke
  if(any(lb/=lbound(target%rwf)).or.any(ub/=ubound(target%rwf)))then
   status=11;return
  endif
  a=max(ks,lbound(source%rwf,6));b=min(ke,ubound(source%rwf,6))
  if(a<=b)target%rwf(:,:,:,:,:,a:b,:)=source%rwf(:,:,:,:,:,a:b,:)
 endif
 if(first)then
  if(allocated(target%zwf))deallocate(target%zwf)
  if(allocated(source%zwf))then
   lb=lbound(source%zwf);ub=ubound(source%zwf);lb(6)=ks;ub(6)=ke
   allocate(target%zwf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4), &
    lb(5):ub(5),lb(6):ub(6),lb(7):ub(7)),stat=istat)
   if(istat/=0)then
    status=4;return
   endif
   target%zwf=0d0
  endif
 endif
 if(allocated(target%zwf).neqv.allocated(source%zwf))then
  status=11;return
 endif
 if(allocated(source%zwf))then
  lb=lbound(source%zwf);ub=ubound(source%zwf);lb(6)=ks;ub(6)=ke
  if(any(lb/=lbound(target%zwf)).or.any(ub/=ubound(target%zwf)))then
   status=11;return
  endif
  a=max(ks,lbound(source%zwf,6));b=min(ke,ubound(source%zwf,6))
  if(a<=b)target%zwf(:,:,:,:,:,a:b,:)=source%zwf(:,:,:,:,:,a:b,:)
 endif
end subroutine
