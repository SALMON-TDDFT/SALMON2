! Direct molecular repartition. Data movement only; no Hamiltonian or elapsed step.
subroutine isolated_repartition(status)
 use salmonpy_grid_field,only:grid_validate
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 integer,intent(out)::status
 type(isolated_savepoint),allocatable::donor
 integer::u,ios,cs,j,old(6),saved_layout(7),header(8),i,k,lo(3),hi(3),a(3),b(3),ot(4),dst(4),c(4)
 integer::reference_step,grid_status,frame,extra
 real(8)::saved_dt,cell(3,3)
 character(32)::magic
 character(80)::filename
 logical::first,scalar_first,scalar_owner
 logical,allocatable::covered(:,:,:,:),waves(:,:,:,:,:)
 status=0
 open(newunit=u,file='isolated-repartition.dat',status='old',action='read',iostat=ios)
 if(ios==0)then
  read(u,*,iostat=ios)old
  close(u,iostat=cs)
  if(cs/=0)ios=cs
 endif
 if(ios/=0)status=4
 if(status==0)then
  if(any(old<1).or.old(2)/=1.or.product(real(old(2:),8))/=real(old(1),8))status=11
 endif
 if(system%xc_payload%use_tau_operator.or.allocated(rt%tpsi0%zwf).or.allocated(rt%tpsi0%rwf))status=15
 call comm_get_max(status,nproc_group_global)
 if(status/=0)return
 dst=[info%io_s,info%io_e,info%ik_s,info%ik_e]
 allocate(covered(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1), &
  waves(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),dst(1):dst(2),dst(3):dst(4)),stat=ios)
 if(ios/=0)status=4
 call comm_get_max(status,nproc_group_global)
 if(status/=0)return
 covered=.false.;waves=.false.;first=.true.;scalar_first=.true.;reference_step=-1
 worker_checkpoint_limit=max(64_8*int(system%ngrid,8),64_8*(int(nt,8)+1_8))
 if(allocated(rt%near_grid_values)) &
  worker_checkpoint_limit=max(worker_checkpoint_limit,size(rt%near_grid_values,kind=8))
 do j=0,old(1)-1
  allocate(donor,stat=ios)
  if(ios/=0)then
   status=4;exit
  endif
  donor%rt%max_steps=rt%max_steps;donor%rt%time_step=rt%time_step
  donor%rt%recover_errors=rt%recover_errors
  donor%rt%factor_endpoint_before_output=rt%factor_endpoint_before_output
  write(filename,'(a,i8.8,a)')'isolated-state-',j,'.bin'
  open(newunit=u,file=trim(filename),status='old',access='stream',form='unformatted',action='read',iostat=ios)
  if(ios==0)then
   read(u,iostat=ios)magic,header,saved_dt,cell,saved_layout
   if(ios==0)then
    if(magic/='SALMONPY_ISOLATED_MPI_V1'.or.any(header/=[1,nt,system%nspin,system%no,system%nk,lg%num]))ios=1
    if(saved_dt/=dt.or.any(cell/=system%primitive_a).or.any(saved_layout/=[old(1),j,old(2:)]))ios=1
   endif
   if(ios/=0)write(*,*)'Molecular repartition: incompatible checkpoint header',j,ios
   if(ios==0)then
    call isolated_payload_io(u,.true.,donor,ios)
    if(ios/=0)write(*,*)'Molecular repartition: payload read failed',j,ios
   endif
   if(ios==0)then
    read(u,iostat=ios)magic
    if(magic/='SALMONPY_ISOLATED_COMPLETE_V7')ios=1
    if(ios==0)then
     read(u,iostat=ios)extra
     if(is_iostat_end(ios))then
      ios=0
     else
      ios=1
     endif
    endif
    if(ios/=0)write(*,*)'Molecular repartition: incomplete checkpoint',j,ios
   endif
   close(u,iostat=cs)
   if(cs/=0)ios=cs
  endif
  if(ios/=0)then
   status=4;exit
  endif
  if(.not.ieee_is_finite(donor%cval))status=11
  do k=1,4
  do i=1,3
   if(.not.ieee_is_finite(donor%fields(i,k)))status=11
  enddo
  enddo
  if(allocated(donor%rt%near_grid_values))then
   do frame=0,nt+1
    call grid_validate(donor%rt%near_grid_counts,donor%rt%near_grid_origins, &
     donor%rt%near_grid_spacing,donor%rt%near_grid_values(:,:,:,:,:,frame),grid_status)
    if(grid_status/=0)status=11
   enddo
  endif
  if(status/=0)exit
  if(first)then
   reference_step=donor%step
   worker_saved%rt=donor%rt
   ! Retain new topology's static multipole/FFT layouts and dead Poisson scratch.
   worker_saved%poisson=poisson;worker_saved%poisson%iterVh=donor%poisson%iterVh
   worker_saved%energy=donor%energy;worker_saved%fields=donor%fields
   worker_saved%Vbox=donor%Vbox
   worker_saved%cval=donor%cval;worker_saved%step=donor%step
   worker_saved%xc=donor%xc
   if(allocated(worker_saved%rt%near_projectors))deallocate(worker_saved%rt%near_projectors)
   if(allocated(worker_saved%rt%near_Ac%v))deallocate(worker_saved%rt%near_Ac%v)
   if(allocated(worker_saved%rt%near_div%f))deallocate(worker_saved%rt%near_div%f)
  endif
  if(donor%step/=reference_step.or.reference_step<0.or.reference_step>nt)then
   status=11;exit
  endif
  if(.not.allocated(donor%rho%f).or..not.allocated(donor%spsi_in%zwf))then
   status=15;exit
  endif
  lo=lbound(donor%rho%f);hi=ubound(donor%rho%f)
  if(any(lo<lg%is).or.any(hi>lg%ie).or.any(lo>hi))then
   status=11;exit
  endif
  ot=[lbound(donor%spsi_in%zwf,5),ubound(donor%spsi_in%zwf,5), &
      lbound(donor%spsi_in%zwf,6),ubound(donor%spsi_in%zwf,6)]
  a=max(lo,mg%is);b=min(hi,mg%ie)
  c=[max(ot(1),dst(1)),min(ot(2),dst(2)),max(ot(3),dst(3)),min(ot(4),dst(4))]
  if(all(a<=b).and.c(1)<=c(2).and.c(3)<=c(4))then
   if(any(waves(a(1):b(1),a(2):b(2),a(3):b(3),c(1):c(2),c(3):c(4))))then
    status=11;exit
   endif
   waves(a(1):b(1),a(2):b(2),a(3):b(3),c(1):c(2),c(3):c(4))=.true.
  endif
  call molecular_merge_orbital(worker_saved%spsi_in,donor%spsi_in,first,lo,hi,mg,lg%num,status,dst)
  call molecular_merge_orbital(worker_saved%spsi_out,donor%spsi_out,first,lo,hi,mg,lg%num,status,dst)
  call molecular_merge_orbital(worker_saved%tpsi,donor%tpsi,first,lo,hi,mg,lg%num,status,dst)
  scalar_owner=ot(1)==1.and.ot(3)==1
  if(scalar_owner)then
   if(all(a<=b))then
    if(any(covered(a(1):b(1),a(2):b(2),a(3):b(3),1)))then
     status=11;exit
    endif
    covered(a(1):b(1),a(2):b(2),a(3):b(3),1)=.true.
   endif
   call molecular_merge_scalar(worker_saved%rho,donor%rho,scalar_first,lo,hi,mg,status)
   call molecular_merge_scalar(worker_saved%rho_jm,donor%rho_jm,scalar_first,lo,hi,mg,status)
   call molecular_merge_scalar(worker_saved%Vh,donor%Vh,scalar_first,lo,hi,mg,status)
   call molecular_merge_scalar(worker_saved%Vh_stock1,donor%Vh_stock1,scalar_first,lo,hi,mg,status)
   call molecular_merge_scalar(worker_saved%Vh_stock2,donor%Vh_stock2,scalar_first,lo,hi,mg,status)
   call molecular_merge_scalar(worker_saved%Vpsl,donor%Vpsl,scalar_first,lo,hi,mg,status)
   call molecular_merge_scalar(worker_saved%rt%vonf,donor%rt%vonf,scalar_first,lo,hi,mg,status)
   if(allocated(donor%rho_s))then
    do i=lbound(donor%rho_s,1),ubound(donor%rho_s,1)
     call molecular_merge_scalar(worker_saved%rho_s(i),donor%rho_s(i),scalar_first,lo,hi,mg,status)
    enddo
   endif
   if(allocated(donor%V_local))then
    do i=lbound(donor%V_local,1),ubound(donor%V_local,1)
     call molecular_merge_scalar(worker_saved%V_local(i),donor%V_local(i),scalar_first,lo,hi,mg,status)
    enddo
   endif
   if(allocated(donor%Vxc))then
    do i=lbound(donor%Vxc,1),ubound(donor%Vxc,1)
     call molecular_merge_scalar(worker_saved%Vxc(i),donor%Vxc(i),scalar_first,lo,hi,mg,status)
    enddo
   endif
   if(allocated(donor%rt%vloc_t))then
    do i=lbound(donor%rt%vloc_t,1),ubound(donor%rt%vloc_t,1)
     call molecular_merge_scalar(worker_saved%rt%vloc_t(i),donor%rt%vloc_t(i),scalar_first,lo,hi,mg,status)
    enddo
   endif
   if(allocated(donor%rt%vloc_new))then
    do i=lbound(donor%rt%vloc_new,1),ubound(donor%rt%vloc_new,1)
     call molecular_merge_scalar(worker_saved%rt%vloc_new(i),donor%rt%vloc_new(i),scalar_first,lo,hi,mg,status)
    enddo
   endif
   if(allocated(donor%rt%rho0_s))then
    do i=lbound(donor%rt%rho0_s,1),ubound(donor%rt%rho0_s,1)
     call molecular_merge_scalar(worker_saved%rt%rho0_s(i),donor%rt%rho0_s(i),scalar_first,lo,hi,mg,status)
    enddo
   endif
   if(allocated(donor%rt%vloc0))then
    do i=lbound(donor%rt%vloc0,1),ubound(donor%rt%vloc0,1)
     call molecular_merge_scalar(worker_saved%rt%vloc0(i),donor%rt%vloc0(i),scalar_first,lo,hi,mg,status)
    enddo
   endif
   if(allocated(donor%rt%vloc_old))then
    do k=lbound(donor%rt%vloc_old,2),ubound(donor%rt%vloc_old,2)
    do i=lbound(donor%rt%vloc_old,1),ubound(donor%rt%vloc_old,1)
     call molecular_merge_scalar(worker_saved%rt%vloc_old(i,k),donor%rt%vloc_old(i,k),scalar_first,lo,hi,mg,status)
    enddo
    enddo
   endif
   if(allocated(donor%rt%j_e%v))status=15
   scalar_first=.false.
  endif
  if(status/=0)exit
  first=.false.;deallocate(donor)
 enddo
 if(status==0)then
  if(.not.all(covered).or..not.all(waves).or.scalar_first)status=11
 endif
 call comm_get_max(status,nproc_group_global)
 if(status/=0)return
 ! Projector/vector-potential values are derived from the saved field without propagation.
 ! Their layout-specific values are rebuilt on publish, not copied from another rank.
end subroutine

subroutine molecular_merge_scalar(target,source,first,lo,hi,mg,status)
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

subroutine molecular_merge_orbital(target,source,first,lo,hi,mg,n,status,tile)
 type(s_orbital),intent(inout)::target
 type(s_orbital),intent(in)::source
 logical,intent(in)::first
 integer,intent(in)::lo(3),hi(3),n(3)
 type(s_rgrid),intent(in)::mg
 integer,intent(inout)::status
 integer,optional,intent(in)::tile(4)
 integer::lb(7),ub(7),istat,x,y,z,g(3),c(4)
 if(status/=0)return
 if(first)target%update_zwf_overlap=source%update_zwf_overlap
 if(first)then
  if(allocated(target%rwf))deallocate(target%rwf)
  if(allocated(source%rwf))then
   lb=lbound(source%rwf);ub=ubound(source%rwf)
   lb(:3)=mg%is_array;ub(:3)=mg%ie_array
   if(present(tile))then
    lb(5:6)=tile(1:4:2);ub(5:6)=tile(2:4:2)
   endif
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
   if(present(tile))then
    lb(5:6)=tile(1:4:2);ub(5:6)=tile(2:4:2)
   endif
  if(any(lb/=lbound(target%rwf)).or.any(ub/=ubound(target%rwf)))then
   status=11;return
  endif
  do z=mg%is_array(3),mg%ie_array(3)
  do y=mg%is_array(2),mg%ie_array(2)
  do x=mg%is_array(1),mg%ie_array(1)
   g=[x,y,z]
   if(all(g>=lo).and.all(g<=hi))then
    c=[max(lbound(target%rwf,5),lbound(source%rwf,5)), &
       min(ubound(target%rwf,5),ubound(source%rwf,5)), &
       max(lbound(target%rwf,6),lbound(source%rwf,6)), &
       min(ubound(target%rwf,6),ubound(source%rwf,6))]
    if(c(1)<=c(2).and.c(3)<=c(4))target%rwf(x,y,z,:,c(1):c(2),c(3):c(4),:)= &
     source%rwf(g(1),g(2),g(3),:,c(1):c(2),c(3):c(4),:)
   endif
  enddo
  enddo
  enddo
 endif
 if(first)then
  if(allocated(target%zwf))deallocate(target%zwf)
  if(allocated(source%zwf))then
   lb=lbound(source%zwf);ub=ubound(source%zwf)
   lb(:3)=mg%is_array;ub(:3)=mg%ie_array
   if(present(tile))then
    lb(5:6)=tile(1:4:2);ub(5:6)=tile(2:4:2)
   endif
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
   if(present(tile))then
    lb(5:6)=tile(1:4:2);ub(5:6)=tile(2:4:2)
   endif
  if(any(lb/=lbound(target%zwf)).or.any(ub/=ubound(target%zwf)))then
   status=11;return
  endif
  do z=mg%is_array(3),mg%ie_array(3)
  do y=mg%is_array(2),mg%ie_array(2)
  do x=mg%is_array(1),mg%ie_array(1)
   g=[x,y,z]
   if(all(g>=lo).and.all(g<=hi))then
    c=[max(lbound(target%zwf,5),lbound(source%zwf,5)), &
       min(ubound(target%zwf,5),ubound(source%zwf,5)), &
       max(lbound(target%zwf,6),lbound(source%zwf,6)), &
       min(ubound(target%zwf,6),ubound(source%zwf,6))]
    if(c(1)<=c(2).and.c(3)<=c(4))target%zwf(x,y,z,:,c(1):c(2),c(3):c(4),:)= &
     source%zwf(g(1),g(2),g(3),:,c(1):c(2),c(3):c(4),:)
   endif
  enddo
  enddo
  enddo
 endif
end subroutine
