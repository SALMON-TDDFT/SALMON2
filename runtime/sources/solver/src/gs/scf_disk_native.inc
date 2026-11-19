! SCF rank-local disk v1. Host authenticates binary/input/topology and all files.
integer(c_int) function sp_gs_hse_orbital_protocol() bind(C)
#ifdef USE_HSE
 sp_gs_hse_orbital_protocol=2
#else
 sp_gs_hse_orbital_protocol=0
#endif
end function
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
  value%hse%spatial%gauge_tiles%comm=ctx%info%icomm_ro
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
integer(c_int) function sp_gs_scf_o_repartition_protocol() bind(C)
#ifdef USE_HSE
 sp_gs_scf_o_repartition_protocol=2
#else
 sp_gs_scf_o_repartition_protocol=1
#endif
end function
integer(c_int) function sp_gs_scf_k_repartition(token,path,old_size,kbounds) bind(C) result(status)
 integer(c_int64_t),value::token
 character(c_char),intent(in)::path(*)
 integer(c_int),value::old_size
 integer(c_int),intent(in)::kbounds(2,*)
 status=scf_repartition(token,path,old_size,kbounds,6)
end function
integer(c_int) function sp_gs_scf_o_repartition(token,path,old_size,obounds) bind(C) result(status)
 integer(c_int64_t),value::token
 character(c_char),intent(in)::path(*)
 integer(c_int),value::old_size
 integer(c_int),intent(in)::obounds(2,*)
 status=scf_repartition(token,path,old_size,obounds,5)
end function
integer(c_int) function scf_repartition(token,path,old_size,kbounds,axis) result(status)
 integer(c_int64_t),value::token
 character(c_char),intent(in)::path(*)
 integer,intent(in)::axis
 integer(c_int),value::old_size
 integer(c_int),intent(in)::kbounds(2,*)
 type(s_scf_payload),allocatable::value,donor
 character(4096)::directory,filename,expected_directory
 integer::index,group,rank,nproc,i,j,u,ios,cs,header(5),counters(2),ks,ke,lo,hi,iteration,old_reference,pair(2),last,dim
 integer(8)::limit,position
 logical::first
 logical,allocatable::covered(:)
 call gs_lookup(token,index,status)
 if(status/=0)return
 call comm_get_globalinfo(group,rank,nproc)
 associate(ctx=>gs_handles(index)%context)
 if(.not.scf_context_ok(ctx))status=15
 if(allocated(ctx%scf_pending))status=11
 last=ctx%system%nk;dim=6
 if(axis==5)then
  last=ctx%system%no;dim=5
  if(ctx%info%isize_r/=1.or.ctx%info%isize_k/=1.or.ctx%info%isize_o/=nproc)status=15
 else
  if(ctx%info%isize_r/=1.or.ctx%info%isize_o/=1.or.ctx%info%isize_k/=nproc)status=15
 endif
 select case(trim(ctx%system%xc_input%xcname))
 case('pz','libxc_pz','libxc_pbe')
#ifdef USE_HSE
 case('hse06')
  if(.not.ctx%hse_bound)status=15
  if(axis==5.and..not.ctx%hse%orbital_scf)status=15
#endif
 case default
  status=15
 end select
 if(old_size<1.or.old_size>last)status=11
 if(nproc>last)status=15
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
 if(axis==5)then
  ks=ctx%info%io_s;ke=ctx%info%io_e
 endif
 allocate(covered(ks:ke),stat=ios)
 if(ios/=0)status=4
 call comm_get_max(status,group)
 if(status/=0)return
 covered=.false.
 first=.true.;iteration=-1
 do j=1,old_size
  lo=kbounds(1,j);hi=kbounds(2,j)
  if(lo<1.or.hi<lo.or.hi>last)then
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
  if(allocated(donor%spsi%zwf))then
   if(lbound(donor%spsi%zwf,dim)/=lo.or.ubound(donor%spsi%zwf,dim)/=hi)status=11
  else if(axis==5.and.allocated(donor%spsi%rwf))then
   if(lbound(donor%spsi%rwf,dim)/=lo.or.ubound(donor%spsi%rwf,dim)/=hi)status=11
  else
   status=15
  endif
  if(status/=0)exit
  call scf_merge_k_orbital(value%spsi,donor%spsi,first,ks,ke,status,axis)
  call scf_merge_k_orbital(value%shpsi,donor%shpsi,first,ks,ke,status,axis)
  call scf_merge_k_orbital(value%sttpsi,donor%sttpsi,first,ks,ke,status,axis)
  call scf_merge_k_orbital(value%cg%pre_gk,donor%cg%pre_gk,first,ks,ke,status,axis)
  call scf_merge_k_orbital(value%cg%pk,donor%cg%pk,first,ks,ke,status,axis)
  call scf_merge_k_orbital(value%cg%hpk,donor%cg%hpk,first,ks,ke,status,axis)
#ifdef USE_HSE
  if(ctx%hse_bound)then
   if(axis==5)then
    call scf_merge_o_hse(value%hse,donor%hse,first,lo,hi,ks,ke,status)
   else
    call scf_merge_k_hse(value%hse,donor%hse,first,lo,hi,ks,ke,status)
   endif
  endif
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
  value%hse%spatial%gauge_tiles%comm=ctx%info%icomm_ro
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
subroutine scf_merge_k_orbital(target,source,first,ks,ke,status,axis)
 type(s_orbital),intent(inout)::target
 type(s_orbital),intent(in)::source
 logical,intent(in)::first
 integer,intent(in)::ks,ke
 integer,intent(inout)::status
 integer,intent(in)::axis
 integer::lb(7),ub(7),a,b,istat
 if(status/=0)return
 if(first)target%update_zwf_overlap=source%update_zwf_overlap
 if(first)then
  if(allocated(target%rwf))deallocate(target%rwf)
  if(allocated(source%rwf))then
   lb=lbound(source%rwf);ub=ubound(source%rwf);lb(axis)=ks;ub(axis)=ke
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
  lb=lbound(source%rwf);ub=ubound(source%rwf);lb(axis)=ks;ub(axis)=ke
  if(any(lb/=lbound(target%rwf)).or.any(ub/=ubound(target%rwf)))then
   status=11;return
  endif
  a=max(ks,lbound(source%rwf,axis));b=min(ke,ubound(source%rwf,axis))
  if(a<=b)then
   if(axis==5)then
    target%rwf(:,:,:,:,a:b,:,:)=source%rwf(:,:,:,:,a:b,:,:)
   else
    target%rwf(:,:,:,:,:,a:b,:)=source%rwf(:,:,:,:,:,a:b,:)
   endif
  endif
 endif
 if(first)then
  if(allocated(target%zwf))deallocate(target%zwf)
  if(allocated(source%zwf))then
   lb=lbound(source%zwf);ub=ubound(source%zwf);lb(axis)=ks;ub(axis)=ke
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
  lb=lbound(source%zwf);ub=ubound(source%zwf);lb(axis)=ks;ub(axis)=ke
  if(any(lb/=lbound(target%zwf)).or.any(ub/=ubound(target%zwf)))then
   status=11;return
  endif
  a=max(ks,lbound(source%zwf,axis));b=min(ke,ubound(source%zwf,axis))
  if(a<=b)then
   if(axis==5)then
    target%zwf(:,:,:,:,a:b,:,:)=source%zwf(:,:,:,:,a:b,:,:)
   else
    target%zwf(:,:,:,:,:,a:b,:)=source%zwf(:,:,:,:,:,a:b,:)
   endif
  endif
 endif
end subroutine

#ifdef USE_HSE
subroutine scf_merge_k_complex3(target,source,first,lo,hi,ks,ke,status,axis)
 complex(8),allocatable,intent(inout)::target(:,:,:)
 complex(8),allocatable,intent(in)::source(:,:,:)
 logical,intent(in)::first
 integer,intent(in)::lo,hi,ks,ke
 integer,intent(inout)::status
 integer,optional,intent(in)::axis
 integer::lb(3),ub(3),istat,a,b,j,dim
 dim=3
 if(present(axis))dim=axis
 if(status/=0)return
 if(first)then
  if(allocated(target))deallocate(target)
  if(allocated(source))then
   lb=lbound(source);ub=ubound(source)
   lb(dim)=1;ub(dim)=ke-ks+1
   allocate(target(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)),stat=istat)
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
 if(size(source,dim)/=hi-lo+1)then
  status=11;return
 endif
 do j=1,3
  if(j==dim)cycle
  if(lbound(target,j)/=lbound(source,j).or.ubound(target,j)/=ubound(source,j))then
   status=11;return
  endif
 enddo
 a=max(ks,lo);b=min(ke,hi)
 if(a<=b)then
  if(dim==2)then
   target(:,a-ks+1:b-ks+1,:)=source(:,a-lo+lbound(source,2):b-lo+lbound(source,2),:)
  else
   target(:,:,a-ks+1:b-ks+1)=source(:,:,a-lo+lbound(source,3):b-lo+lbound(source,3))
  endif
 endif
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
subroutine scf_merge_o_hse(target,source,first,lo,hi,ks,ke,status)
 type(s_exx_rt_state),intent(inout)::target
 type(s_exx_rt_state),intent(in)::source
 logical,intent(in)::first
 integer,intent(in)::lo,hi,ks,ke
 integer,intent(inout)::status
 complex(8),allocatable::buffer(:,:,:),merged(:,:,:)
 integer::istat
 if(status/=0)return
 if(.not.source%orbital_scf.or.allocated(source%spatial%gauge_tiles%matrix).or. &
    allocated(source%spatial%sparse_source%offset).or.allocated(source%factor_history).or. &
    allocated(source%ptcn_history))then
  status=15;return
 endif
 call merge_ace(target%ace,source%ace)
 call merge_ace(target%scf_damped_ace,source%scf_damped_ace)
 call merge_ace(target%initial_ace,source%initial_ace)
 call merge_ace(target%midpoint_ace,source%midpoint_ace)
 call scf_merge_k_complex3(target%cached_source,source%cached_source,first,lo,hi,ks,ke,status,2)
 call scf_merge_k_complex3(target%cached_action,source%cached_action,first,lo,hi,ks,ke,status,2)
 call scf_merge_k_complex3(target%full_source,source%full_source,first,lo,hi,ks,ke,status,2)
 call scf_merge_k_complex3(target%initial_source,source%initial_source,first,lo,hi,ks,ke,status,2)
 call scf_merge_k_complex3(target%target_work,source%target_work,first,lo,hi,ks,ke,status,2)
 call scf_merge_k_complex3(target%action_work,source%action_work,first,lo,hi,ks,ke,status,2)
 call scf_merge_k_complex3(target%spatial%previous,source%spatial%previous,first,lo,hi,ks,ke,status,2)
 if(status/=0)return
 ! The source matrix has the same orbital-column layout as previous(:,:,1).
 if(allocated(source%spatial%source))then
  allocate(buffer(size(source%spatial%source,1),size(source%spatial%source,2),1),stat=istat)
  if(istat/=0)then
   status=4;return
  endif
  buffer(:,:,1)=source%spatial%source
 endif
 if(.not.first.and.allocated(target%spatial%source))then
  allocate(merged(size(target%spatial%source,1),size(target%spatial%source,2),1),stat=istat)
  if(istat/=0)then
   status=4;return
  endif
  merged(:,:,1)=target%spatial%source
 endif
 call scf_merge_k_complex3(merged,buffer,first,lo,hi,ks,ke,status,2)
 if(status/=0)return
 if(allocated(target%spatial%source))deallocate(target%spatial%source)
 if(allocated(merged))then
  allocate(target%spatial%source(size(merged,1),size(merged,2)),stat=istat)
  if(istat/=0)then
   status=4;return
  endif
  target%spatial%source=merged(:,:,1)
 endif
contains
 subroutine merge_ace(left,right)
  type(s_exx_ace),intent(inout)::left
  type(s_exx_ace),intent(in)::right
  if(right%packed.or.right%metric_distributed)then
   status=15;return
  endif
  if(.not.first)left%condition=max(left%condition,right%condition)
  call scf_merge_k_complex3(left%factors,right%factors,first,lo,hi,ks,ke,status,2)
 end subroutine
end subroutine
#endif

integer(c_int) function sp_gs_scf_r_repartition_protocol() bind(C)
 sp_gs_scf_r_repartition_protocol=1
end function
integer(c_int) function sp_gs_scf_r_repartition(token,path,old_size,boxes) bind(C) result(status)
 integer(c_int64_t),value::token
 character(c_char),intent(in)::path(*)
 integer(c_int),value::old_size
 integer(c_int),intent(in)::boxes(6,*)
 type(s_scf_payload),allocatable::value,donor
 character(4096)::directory,filename,reference_directory
 integer::index,group,rank,nproc,j,i,u,ios,cs,header(5),counters(2),reference,pair(6)
 integer::lo(3),hi(3),a(3),b(3),iteration
 integer(8)::limit,position
 logical::first
 logical,allocatable::covered(:,:,:)
 call gs_lookup(token,index,status)
 if(status/=0)return
 call comm_get_globalinfo(group,rank,nproc)
 associate(ctx=>gs_handles(index)%context)
 if(.not.scf_context_ok(ctx).or.ctx%info%isize_k/=1.or.ctx%info%isize_o/=1)status=15
 if(allocated(ctx%scf_pending))status=11
 select case(trim(ctx%system%xc_input%xcname))
 case('pz','libxc_pz','libxc_pbe')
 case default
  status=15
 end select
 if(old_size<1)status=11
 directory=''
 do i=1,len(directory)
  if(path(i)==c_null_char)exit
  directory(i:i)=path(i)
 enddo
 if(i>len(directory).or.len_trim(directory)==0)status=11
 reference=old_size;call comm_bcast(reference,group)
 reference_directory=directory;call comm_bcast(reference_directory,group)
 if(reference/=old_size.or.reference_directory/=directory)status=11
 call comm_get_max(status,group)
 if(status/=0)return
 do j=1,old_size
  pair=boxes(:,j);call comm_bcast(pair,group)
  if(any(pair/=boxes(:,j)))status=11
 enddo
 call comm_get_max(status,group)
 if(status/=0)return
 allocate(covered(ctx%mg%is(1):ctx%mg%ie(1),ctx%mg%is(2):ctx%mg%ie(2), &
  ctx%mg%is(3):ctx%mg%ie(3)),stat=ios)
 if(ios/=0)status=4
 call comm_get_max(status,group)
 if(status/=0)return
 covered=.false.;first=.true.;iteration=-1
 do j=1,old_size
  lo=boxes(1:6:2,j);hi=boxes(2:6:2,j)
  if(any(lo<1).or.any(hi<lo).or.any(hi>ctx%lg%num))then
   status=11;exit
  endif
  a=max(lo,ctx%mg%is);b=min(hi,ctx%mg%ie)
  if(all(a<=b))then
   if(any(covered(a(1):b(1),a(2):b(2),a(3):b(3))))then
    status=11;exit
   endif
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
   ! Static layout/transform resources and dead FFT/XC scratch belong to the
   ! new host topology. Accepted density, potentials and mixing are mapped below.
   value%fg=ctx%fg;value%poisson=ctx%poisson
   value%system%xc_work=ctx%system%xc_work
   value%system%xc_work%rv_reported=donor%system%xc_work%rv_reported
  endif
  if(donor%iterations/=iteration.or.donor%system%ngrid/=ctx%system%ngrid.or. &
     donor%system%nk/=ctx%system%nk.or.donor%system%no/=ctx%system%no)then
   status=11;exit
  endif
  call scf_merge_r_payload(value,donor,first,lo,hi,ctx%mg,ctx%lg%num,status)
  if(status/=0)exit
  if(all(a<=b))covered(a(1):b(1),a(2):b(2),a(3):b(3))=.true.
  first=.false.;deallocate(donor)
 enddo
 if(status==0.and.(.not.all(covered).or..not.allocated(value)))status=11
 call comm_get_max(status,group)
 if(status/=0)return
 call value%continuation_inputs%parallel_from(ctx%inputs)
 call move_alloc(value,ctx%scf_pending)
 end associate
end function

subroutine scf_merge_r_scalar(target,source,first,lo,hi,mg,status)
 type(s_scalar),intent(inout)::target
 type(s_scalar),intent(in)::source
 logical,intent(in)::first
 integer,intent(in)::lo(3),hi(3)
 type(s_rgrid),intent(in)::mg
 integer,intent(inout)::status
 integer::a(3),b(3),istat
 if(status/=0)return
 if(first)then
  if(allocated(target%f))deallocate(target%f)
  if(allocated(source%f))then
   allocate(target%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)),stat=istat)
   if(istat/=0)then
    status=4;return
   endif
   target%f=0d0
  endif
 endif
 if(allocated(target%f).neqv.allocated(source%f))then
  status=11;return
 endif
 if(.not.allocated(source%f))return
 if(any(lbound(source%f)/=lo).or.any(ubound(source%f)/=hi))then
  status=11;return
 endif
 a=max(lo,mg%is);b=min(hi,mg%ie)
 if(all(a<=b))target%f(a(1):b(1),a(2):b(2),a(3):b(3))=source%f(a(1):b(1),a(2):b(2),a(3):b(3))
end subroutine

subroutine scf_merge_r_orbital(target,source,first,lo,hi,mg,n,status)
 type(s_orbital),intent(inout)::target
 type(s_orbital),intent(in)::source
 logical,intent(in)::first
 integer,intent(in)::lo(3),hi(3),n(3)
 type(s_rgrid),intent(in)::mg
 integer,intent(inout)::status
 integer::lb(7),ub(7),istat,x,y,z,g(3)
 if(status/=0)return
 if(first)target%update_zwf_overlap=source%update_zwf_overlap
 if(first)then
  if(allocated(target%rwf))deallocate(target%rwf)
  if(allocated(source%rwf))then
   lb=lbound(source%rwf);ub=ubound(source%rwf)
   lb(:3)=mg%is_array;ub(:3)=mg%ie_array
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
  lb=lbound(source%rwf);ub=ubound(source%rwf)
  if(any(lb(:3)>lo).or.any(ub(:3)<hi))then
   status=11;return
  endif
  lb(:3)=mg%is_array;ub(:3)=mg%ie_array
  if(any(lb/=lbound(target%rwf)).or.any(ub/=ubound(target%rwf)))then
   status=11;return
  endif
  do z=mg%is_array(3),mg%ie_array(3)
  do y=mg%is_array(2),mg%ie_array(2)
  do x=mg%is_array(1),mg%ie_array(1)
   g=1+modulo([x,y,z]-1,n)
   if(all(g>=lo).and.all(g<=hi))target%rwf(x,y,z,:,:,:,:)=source%rwf(g(1),g(2),g(3),:,:,:,:)
  enddo
  enddo
  enddo
 endif
 if(first)then
  if(allocated(target%zwf))deallocate(target%zwf)
  if(allocated(source%zwf))then
   lb=lbound(source%zwf);ub=ubound(source%zwf)
   lb(:3)=mg%is_array;ub(:3)=mg%ie_array
   allocate(target%zwf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4), &
    lb(5):ub(5),lb(6):ub(6),lb(7):ub(7)),stat=istat)
   if(istat/=0)then
    status=4;return
   endif
   target%zwf=(0d0,0d0)
  endif
 endif
 if(allocated(target%zwf).neqv.allocated(source%zwf))then
  status=11;return
 endif
 if(allocated(source%zwf))then
  lb=lbound(source%zwf);ub=ubound(source%zwf)
  if(any(lb(:3)>lo).or.any(ub(:3)<hi))then
   status=11;return
  endif
  lb(:3)=mg%is_array;ub(:3)=mg%ie_array
  if(any(lb/=lbound(target%zwf)).or.any(ub/=ubound(target%zwf)))then
   status=11;return
  endif
  do z=mg%is_array(3),mg%ie_array(3)
  do y=mg%is_array(2),mg%ie_array(2)
  do x=mg%is_array(1),mg%ie_array(1)
   g=1+modulo([x,y,z]-1,n)
   if(all(g>=lo).and.all(g<=hi))target%zwf(x,y,z,:,:,:,:)=source%zwf(g(1),g(2),g(3),:,:,:,:)
  enddo
  enddo
  enddo
 endif
end subroutine
subroutine scf_merge_r_payload(target,source,first,lo,hi,mg,n,status)
 type(s_scf_payload),intent(inout)::target
 type(s_scf_payload),intent(in)::source
 logical,intent(in)::first
 integer,intent(in)::lo(3),hi(3),n(3)
 type(s_rgrid),intent(in)::mg
 integer,intent(inout)::status
 integer::i,j,istat,a(3),b(3)
 if(status/=0)return
 call scf_merge_r_orbital(target%spsi,source%spsi,first,lo,hi,mg,n,status)
 call scf_merge_r_orbital(target%shpsi,source%shpsi,first,lo,hi,mg,n,status)
 call scf_merge_r_orbital(target%sttpsi,source%sttpsi,first,lo,hi,mg,n,status)
 call scf_merge_r_orbital(target%cg%pre_gk,source%cg%pre_gk,first,lo,hi,mg,n,status)
 call scf_merge_r_orbital(target%cg%pk,source%cg%pk,first,lo,hi,mg,n,status)
 call scf_merge_r_orbital(target%cg%hpk,source%cg%hpk,first,lo,hi,mg,n,status)
 call scf_merge_r_scalar(target%rho,source%rho,first,lo,hi,mg,status)
 call scf_merge_r_scalar(target%rho_jm,source%rho_jm,first,lo,hi,mg,status)
 call scf_merge_r_scalar(target%Vh,source%Vh,first,lo,hi,mg,status)
 call scf_merge_r_scalar(target%Vpsl,source%Vpsl,first,lo,hi,mg,status)
 call scf_merge_r_scalar(target%mixing%scf_progress%rho_old,source%mixing%scf_progress%rho_old,first,lo,hi,mg,status)
 call scf_merge_r_scalar(target%mixing%scf_progress%Vlocal_old,source%mixing%scf_progress%Vlocal_old,first,lo,hi,mg,status)
 if(allocated(target%V_local).neqv.allocated(source%V_local))status=11
 if(status/=0)return
 if(allocated(source%V_local))then
  if(any(shape(target%V_local)/=shape(source%V_local)))then
   status=11;return
  endif
  do i=lbound(source%V_local,1),ubound(source%V_local,1)
   call scf_merge_r_scalar(target%V_local(i),source%V_local(i),first,lo,hi,mg,status)
  enddo
 endif
 if(allocated(target%rho_s).neqv.allocated(source%rho_s))status=11
 if(status/=0)return
 if(allocated(source%rho_s))then
  if(any(shape(target%rho_s)/=shape(source%rho_s)))then
   status=11;return
  endif
  do i=lbound(source%rho_s,1),ubound(source%rho_s,1)
   call scf_merge_r_scalar(target%rho_s(i),source%rho_s(i),first,lo,hi,mg,status)
  enddo
 endif
 if(allocated(target%Vxc).neqv.allocated(source%Vxc))status=11
 if(status/=0)return
 if(allocated(source%Vxc))then
  if(any(shape(target%Vxc)/=shape(source%Vxc)))then
   status=11;return
  endif
  do i=lbound(source%Vxc,1),ubound(source%Vxc,1)
   call scf_merge_r_scalar(target%Vxc(i),source%Vxc(i),first,lo,hi,mg,status)
  enddo
 endif
 if(allocated(target%mixing%rho_in).neqv.allocated(source%mixing%rho_in))status=11
 if(status/=0)return
 if(allocated(source%mixing%rho_in))then
  if(any(shape(target%mixing%rho_in)/=shape(source%mixing%rho_in)))then
   status=11;return
  endif
  do i=lbound(source%mixing%rho_in,1),ubound(source%mixing%rho_in,1)
   call scf_merge_r_scalar(target%mixing%rho_in(i),source%mixing%rho_in(i),first,lo,hi,mg,status)
  enddo
 endif
 if(allocated(target%mixing%rho_out).neqv.allocated(source%mixing%rho_out))status=11
 if(status/=0)return
 if(allocated(source%mixing%rho_out))then
  if(any(shape(target%mixing%rho_out)/=shape(source%mixing%rho_out)))then
   status=11;return
  endif
  do i=lbound(source%mixing%rho_out,1),ubound(source%mixing%rho_out,1)
   call scf_merge_r_scalar(target%mixing%rho_out(i),source%mixing%rho_out(i),first,lo,hi,mg,status)
  enddo
 endif
 if(allocated(target%mixing%rho_s_in).neqv.allocated(source%mixing%rho_s_in))status=11
 if(status/=0)return
 if(allocated(source%mixing%rho_s_in))then
  if(any(shape(target%mixing%rho_s_in)/=shape(source%mixing%rho_s_in)))then
   status=11;return
  endif
  do j=lbound(source%mixing%rho_s_in,2),ubound(source%mixing%rho_s_in,2)
  do i=lbound(source%mixing%rho_s_in,1),ubound(source%mixing%rho_s_in,1)
   call scf_merge_r_scalar(target%mixing%rho_s_in(i,j),source%mixing%rho_s_in(i,j),first,lo,hi,mg,status)
  enddo
  enddo
 endif
 if(allocated(target%mixing%rho_s_out).neqv.allocated(source%mixing%rho_s_out))status=11
 if(status/=0)return
 if(allocated(source%mixing%rho_s_out))then
  if(any(shape(target%mixing%rho_s_out)/=shape(source%mixing%rho_s_out)))then
   status=11;return
  endif
  do j=lbound(source%mixing%rho_s_out,2),ubound(source%mixing%rho_s_out,2)
  do i=lbound(source%mixing%rho_s_out,1),ubound(source%mixing%rho_s_out,1)
   call scf_merge_r_scalar(target%mixing%rho_s_out(i,j),source%mixing%rho_s_out(i,j),first,lo,hi,mg,status)
  enddo
  enddo
 endif
 if(allocated(target%mixing%Vh_in).neqv.allocated(source%mixing%Vh_in))status=11
 if(status/=0)return
 if(allocated(source%mixing%Vh_in))then
  if(any(shape(target%mixing%Vh_in)/=shape(source%mixing%Vh_in)))then
   status=11;return
  endif
  do i=lbound(source%mixing%Vh_in,1),ubound(source%mixing%Vh_in,1)
   call scf_merge_r_scalar(target%mixing%Vh_in(i),source%mixing%Vh_in(i),first,lo,hi,mg,status)
  enddo
 endif
 if(allocated(target%mixing%Vh_out).neqv.allocated(source%mixing%Vh_out))status=11
 if(status/=0)return
 if(allocated(source%mixing%Vh_out))then
  if(any(shape(target%mixing%Vh_out)/=shape(source%mixing%Vh_out)))then
   status=11;return
  endif
  do i=lbound(source%mixing%Vh_out,1),ubound(source%mixing%Vh_out,1)
   call scf_merge_r_scalar(target%mixing%Vh_out(i),source%mixing%Vh_out(i),first,lo,hi,mg,status)
  enddo
 endif
 if(allocated(target%mixing%Vxc_in).neqv.allocated(source%mixing%Vxc_in))status=11
 if(status/=0)return
 if(allocated(source%mixing%Vxc_in))then
  if(any(shape(target%mixing%Vxc_in)/=shape(source%mixing%Vxc_in)))then
   status=11;return
  endif
  do j=lbound(source%mixing%Vxc_in,2),ubound(source%mixing%Vxc_in,2)
  do i=lbound(source%mixing%Vxc_in,1),ubound(source%mixing%Vxc_in,1)
   call scf_merge_r_scalar(target%mixing%Vxc_in(i,j),source%mixing%Vxc_in(i,j),first,lo,hi,mg,status)
  enddo
  enddo
 endif
 if(allocated(target%mixing%Vxc_out).neqv.allocated(source%mixing%Vxc_out))status=11
 if(status/=0)return
 if(allocated(source%mixing%Vxc_out))then
  if(any(shape(target%mixing%Vxc_out)/=shape(source%mixing%Vxc_out)))then
   status=11;return
  endif
  do j=lbound(source%mixing%Vxc_out,2),ubound(source%mixing%Vxc_out,2)
  do i=lbound(source%mixing%Vxc_out,1),ubound(source%mixing%Vxc_out,1)
   call scf_merge_r_scalar(target%mixing%Vxc_out(i,j),source%mixing%Vxc_out(i,j),first,lo,hi,mg,status)
  enddo
  enddo
 endif
 if(status/=0)return
 if(allocated(source%poisson%zrhoG_ele))then
  if(.not.allocated(target%poisson%zrhoG_ele))then
   status=11;return
  endif
  if(any(lbound(source%poisson%zrhoG_ele)/=lo).or.any(ubound(source%poisson%zrhoG_ele)/=hi))then
   status=11;return
  endif
  a=max(lo,mg%is);b=min(hi,mg%ie)
  if(all(a<=b))target%poisson%zrhoG_ele(a(1):b(1),a(2):b(2),a(3):b(3))= &
   source%poisson%zrhoG_ele(a(1):b(1),a(2):b(2),a(3):b(3))
 endif
end subroutine
