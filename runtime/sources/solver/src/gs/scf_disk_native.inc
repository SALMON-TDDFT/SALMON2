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
#ifdef USE_HSE
 sp_gs_scf_k_repartition_protocol=2
#else
 sp_gs_scf_k_repartition_protocol=1
#endif
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
#ifdef USE_HSE
 case('hse06')
  if(.not.ctx%hse_bound)status=15
#endif
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
#ifdef USE_HSE
  if(ctx%hse_bound)call scf_merge_k_hse(value%hse,donor%hse,first,lo,hi,ks,ke,status)
#endif
  if(status/=0)exit
  covered(max(lo,ks):min(hi,ke))=.true.;first=.false.
  deallocate(donor)
 enddo
 if(status==0.and.(.not.all(covered).or..not.allocated(value)))status=11
 call comm_get_max(status,group)
 if(status/=0)return
 call value%continuation_inputs%parallel_from(ctx%inputs)
#ifdef USE_HSE
 if(ctx%hse_bound)then
  value%hse%ace%metric_comm=ctx%info%icomm_r
  value%hse%scf_damped_ace%metric_comm=ctx%info%icomm_r
  value%hse%initial_ace%metric_comm=ctx%info%icomm_r
  value%hse%midpoint_ace%metric_comm=ctx%info%icomm_r
  ! Geometric phase tables are rank-specific runtime caches, not exchange history.
  if(allocated(value%hse%native_symmetry%cached_phase))deallocate(value%hse%native_symmetry%cached_phase)
  if(allocated(value%hse%native_symmetry%phase_offset))deallocate(value%hse%native_symmetry%phase_offset)
  value%hse%native_symmetry%cached_first=1;value%hse%native_symmetry%cached_count=0
 endif
#endif
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

#ifdef USE_HSE
subroutine scf_merge_k_complex3(target,source,first,lo,hi,ks,ke,status)
 complex(8),allocatable,intent(inout)::target(:,:,:)
 complex(8),allocatable,intent(in)::source(:,:,:)
 logical,intent(in)::first
 integer,intent(in)::lo,hi,ks,ke
 integer,intent(inout)::status
 integer::lb(3),ub(3),istat,a,b,j
 if(status/=0)return
 if(first)then
  if(allocated(target))deallocate(target)
  if(allocated(source))then
   lb=lbound(source);ub=ubound(source)
   allocate(target(lb(1):ub(1),lb(2):ub(2),1:ke-ks+1),stat=istat)
   if(istat/=0)then
    status=4;return
   endif
   target=0d0
  endif
 endif
 if(allocated(target).neqv.allocated(source))then
  status=11;return
 endif
 if(.not.allocated(source))return
 if(size(source,3)/=hi-lo+1)then
  status=11;return
 endif
 do j=1,2
  if(lbound(target,j)/=lbound(source,j).or.ubound(target,j)/=ubound(source,j))then
   status=11;return
  endif
 enddo
 a=max(ks,lo);b=min(ke,hi)
 if(a<=b)target(:,:,a-ks+1:b-ks+1)=source(:,:,a-lo+lbound(source,3):b-lo+lbound(source,3))
end subroutine
subroutine scf_merge_k_ace(target,source,first,lo,hi,ks,ke,status)
 type(s_exx_ace),intent(inout)::target
 type(s_exx_ace),intent(in)::source
 logical,intent(in)::first
 integer,intent(in)::lo,hi,ks,ke
 integer,intent(inout)::status
 if(status/=0)return
 if(source%packed.or.source%metric_distributed)then
  status=15;return
 endif
 ! condition is diagnostic-only: retain a conservative bound over contributors.
 if(.not.first)target%condition=max(target%condition,source%condition)
 call scf_merge_k_complex3(target%factors,source%factors,first,lo,hi,ks,ke,status)
end subroutine
subroutine scf_merge_k_history(target,source,first,lo,hi,ks,ke,status)
 type(s_factor_history),allocatable,intent(inout)::target(:)
 type(s_factor_history),allocatable,intent(in)::source(:)
 logical,intent(in)::first
 integer,intent(in)::lo,hi,ks,ke
 integer,intent(inout)::status
 integer::j,istat
 if(status/=0)return
 if(first)then
  if(allocated(target))deallocate(target)
  if(allocated(source))then
   allocate(target(ke-ks+1),stat=istat)
   if(istat/=0)then
    status=4;return
   endif
  endif
 endif
 if(allocated(target).neqv.allocated(source))then
  status=11;return
 endif
 if(.not.allocated(source))return
 if(size(source)/=hi-lo+1)then
  status=11;return
 endif
 do j=max(ks,lo),min(ke,hi)
  target(j-ks+1)=source(j-lo+lbound(source,1))
 enddo
end subroutine
subroutine scf_merge_k_hse(target,source,first,lo,hi,ks,ke,status)
 type(s_exx_rt_state),intent(inout)::target
 type(s_exx_rt_state),intent(in)::source
 logical,intent(in)::first
 integer,intent(in)::lo,hi,ks,ke
 integer,intent(inout)::status
 call scf_merge_k_ace(target%ace,source%ace,first,lo,hi,ks,ke,status)
 call scf_merge_k_ace(target%scf_damped_ace,source%scf_damped_ace,first,lo,hi,ks,ke,status)
 call scf_merge_k_ace(target%initial_ace,source%initial_ace,first,lo,hi,ks,ke,status)
 call scf_merge_k_ace(target%midpoint_ace,source%midpoint_ace,first,lo,hi,ks,ke,status)
 call scf_merge_k_complex3(target%cached_source,source%cached_source,first,lo,hi,ks,ke,status)
 call scf_merge_k_complex3(target%cached_action,source%cached_action,first,lo,hi,ks,ke,status)
 call scf_merge_k_complex3(target%full_source,source%full_source,first,lo,hi,ks,ke,status)
 call scf_merge_k_complex3(target%initial_source,source%initial_source,first,lo,hi,ks,ke,status)
 call scf_merge_k_complex3(target%target_work,source%target_work,first,lo,hi,ks,ke,status)
 call scf_merge_k_complex3(target%action_work,source%action_work,first,lo,hi,ks,ke,status)
 call scf_merge_k_history(target%factor_history,source%factor_history,first,lo,hi,ks,ke,status)
 call scf_merge_k_history(target%ptcn_history,source%ptcn_history,first,lo,hi,ks,ke,status)
end subroutine
#endif
