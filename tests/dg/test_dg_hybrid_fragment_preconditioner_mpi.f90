program test_dg_hybrid_fragment_preconditioner_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_get_halting_mode,ieee_invalid
  use dg_hybrid_fragment_preconditioner
  use dg_hybrid_fragment_subspace,only:s_dg_hybrid_fragment_subspace_state,advance_dg_hybrid_fragment_subspace
  use dc_fragment_occupation,only:determine_dc_fragment_occupations
  implicit none
  integer,parameter::n=12,nw=8,m=3
  integer::comm,rank,nproc,ierr,nlocal,i,j,k,a,variant
  integer(int64),allocatable::ids(:),saved_ids(:)
  complex(real64),allocatable::qr(:,:),hr(:,:),sr(:,:),rr(:,:),output(:,:)
  complex(real64)::u(n,n),q(n,n),h(n,n),s(n,n),r(n,m),physical(n,m),expected(n,m)
  real(real64)::hd(n),sd(n),shifts(m),angle,scale,roundoff,floor,denominator
  type(s_dg_hybrid_preconditioner_key)::key,wrong
  type(s_dg_hybrid_fragment_preconditioner)::cache
  integer(int64)::fingerprint,application_fingerprint,minimum_fp,maximum_fp,old_fp
  logical::ok,traps_before,traps_after
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call ieee_get_halting_mode(ieee_invalid,traps_before)
  nlocal=count([(mod(i-1,min(nproc,4))==rank,i=1,n)])
  allocate(ids(nlocal),qr(nlocal,n),hr(nlocal,n),sr(nlocal,n),rr(nlocal,m));k=0
  do i=n,1,-1
    if(mod(i-1,min(nproc,4))/=rank)cycle
    k=k+1;ids(k)=i
  enddo
  saved_ids=ids
  key%fragment_id=7;key%basis_generation=2;key%operator_epoch=3
  key%basis_fingerprint=101_int64;key%metric_fingerprint=103_int64
  key%operator_fingerprint=107_int64;key%reference_fingerprint=109_int64
  hd=[(0.4d0*i,i=1,n)];sd=[(1d0+0.02d0*i,i=1,n)]
  shifts=[-0.1d0,0.6d0,1.1d0]
  do j=1,m
    do i=1,n;r(i,j)=cmplx(sin(real(i*j,real64)),cos(0.2d0*i*j),real64);enddo
  enddo
  do variant=0,3
    call rotation(variant,u);q=conjg(transpose(u));h=0d0;s=0d0
    do i=1,n;h(i,i)=hd(i);s(i,i)=sd(i);enddo
    h=matmul(conjg(transpose(u)),matmul(h,u));s=matmul(conjg(transpose(u)),matmul(s,u))
    call distribute(q,h,s,matmul(conjg(transpose(u)),r))
    call prepare_dg_hybrid_fragment_preconditioner(comm,n,ids,qr,hr,sr,key,1d-10,cache,fingerprint,ok,message)
    call require(ok,'prepare fixed-frame preconditioner: '//trim(message))
    call MPI_Allreduce(fingerprint,minimum_fp,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(fingerprint,maximum_fp,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    call require(fingerprint/=0_int64.and.minimum_fp==maximum_fp,'rank-dependent preparation receipt')
    call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,rr,output,&
      application_fingerprint,ok,message)
    call require(ok,'apply fixed-frame preconditioner: '//trim(message))
    call gather_physical(output,u,physical)
    do j=1,m;expected(:,j)=r(:,j)/(hd-shifts(j)*sd);enddo
    call require(maxval(abs(physical-expected))<1d-10,'physical action changed under WF gauge')
    call require(maxval(abs(physical-r))>1d-2,'production action became identity')
    old_fp=fingerprint
    wrong=key;wrong%operator_epoch=4
    call apply_dg_hybrid_fragment_preconditioner(comm,ids,wrong,cache,shifts,rr,output,&
      application_fingerprint,ok,message)
    call require(.not.ok.and..not.allocated(output),'stale operator epoch accepted')
  enddo
  if(nproc>1)call test_mixed_layouts()
  do a=1,4
    wrong=key
    select case(a)
    case(1);wrong%basis_fingerprint=9_int64
    case(2);wrong%metric_fingerprint=9_int64
    case(3);wrong%operator_fingerprint=9_int64
    case(4);wrong%reference_fingerprint=9_int64
    end select
    call apply_dg_hybrid_fragment_preconditioner(comm,ids,wrong,cache,shifts,rr,output,&
      application_fingerprint,ok,message)
    call require(.not.ok.and..not.allocated(output),'stale provenance accepted')
  enddo
  if(nproc>1)then
    wrong=key;if(rank==0)wrong%operator_epoch=4
    call prepare_dg_hybrid_fragment_preconditioner(comm,n,ids,qr,hr,sr,wrong,1d-10,cache,fingerprint,ok,message)
    call require(.not.ok,'rank-disagreeing key accepted')
    call prepare_dg_hybrid_fragment_preconditioner(comm,n,ids,qr,hr,sr,key,merge(1d-9,1d-10,rank==0),&
      cache,fingerprint,ok,message)
    call require(.not.ok,'rank-disagreeing tolerance accepted')
  endif
  if(nlocal>0)qr(1,:)=2d0*qr(1,:)
  call prepare_dg_hybrid_fragment_preconditioner(comm,n,ids,qr,hr,sr,key,1d-10,cache,fingerprint,ok,message)
  call require(.not.ok,'nonunitary reference frame accepted')
  call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,rr,output,&
    application_fingerprint,ok,message)
  call require(ok,'failed rebuild destroyed previous valid cache')
  call gather_physical(output,u,physical)
  call require(maxval(abs(physical-expected))<1d-10,'failed rebuild changed previous physical action')
  call distribute(q,h,s,matmul(conjg(transpose(u)),r))
  sr=-sr
  call prepare_dg_hybrid_fragment_preconditioner(comm,n,ids,qr,hr,sr,key,1d-10,cache,fingerprint,ok,message)
  call require(.not.ok,'negative reference metric diagonal accepted');sr=-sr
  if(rank==0.and.nlocal>0)hr(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,real64)
  call prepare_dg_hybrid_fragment_preconditioner(comm,n,ids,qr,hr,sr,key,1d-10,cache,fingerprint,ok,message)
  call require(.not.ok,'nonfinite operator accepted')
  call distribute(q,h,s,matmul(conjg(transpose(u)),r))
  if(nlocal>0)hr(1,int(ids(1)))=hr(1,int(ids(1)))+cmplx(0d0,0.2d0,real64)
  call prepare_dg_hybrid_fragment_preconditioner(comm,n,ids,qr,hr,sr,key,1d-10,cache,fingerprint,ok,message)
  call require(.not.ok,'non-Hermitian operator accepted')
  call distribute(q,h,s,matmul(conjg(transpose(u)),r))
  shifts(1)=ieee_value(0d0,ieee_quiet_nan)
  call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,rr,output,&
    application_fingerprint,ok,message)
  call require(.not.ok.and..not.allocated(output),'nonfinite shifts accepted');shifts(1)=-0.1d0
  shifts(1)=huge(1d0)
  call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,rr,output,&
    application_fingerprint,ok,message)
  call require(.not.ok.and..not.allocated(output),'finite denominator overflow accepted');shifts(1)=-0.1d0
  if(nproc>1)then
    shifts(1)=merge(-0.2d0,-0.1d0,rank==0)
    call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,rr,output,&
      application_fingerprint,ok,message)
    call require(.not.ok.and..not.allocated(output),'rank-disagreeing shifts accepted');shifts(1)=-0.1d0
  endif
  if(rank==0.and.nlocal>0)rr(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,real64)
  call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,rr,output,&
    application_fingerprint,ok,message)
  call require(.not.ok.and..not.allocated(output),'nonfinite residual accepted')
  call distribute(q,h,s,matmul(conjg(transpose(u)),r))
  if(nlocal>1)ids=ids(size(ids):1:-1)
  call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,rr,output,&
    application_fingerprint,ok,message)
  call require(.not.ok,'changed row layout accepted');ids=saved_ids
  if(nlocal>0)ids(1)=1_int64
  call prepare_dg_hybrid_fragment_preconditioner(comm,n,ids,qr,hr,sr,key,1d-10,cache,fingerprint,ok,message)
  call require(.not.ok,'duplicate/missing rows accepted');ids=saved_ids

  ! Denominator regularization: exact zero, small resolvable negative, positive.
  call rotation(0,u);q=u;h=0d0;s=0d0
  hd=1d0;sd=1d0;hd(2)=1d0-1d-8;hd(3)=1d0+1d-8
  do i=1,n;h(i,i)=hd(i);s(i,i)=sd(i);enddo
  key%operator_epoch=4;key%operator_fingerprint=211_int64
  call distribute(q,h,s,r)
  call prepare_dg_hybrid_fragment_preconditioner(comm,n,ids,qr,hr,sr,key,1d-6,cache,fingerprint,ok,message)
  call require(ok,'regularized frame preparation failed: '//trim(message));shifts=1d0
  call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,rr,output,&
    application_fingerprint,ok,message)
  call require(ok,'regularized action failed: '//trim(message))
  call gather_physical(output,u,physical)
  scale=max(1d0,maxval(abs(hd))+maxval(abs(sd)));roundoff=64d0*epsilon(1d0)*n*scale
  floor=max(1d-6*scale,roundoff)
  do i=1,n
    denominator=hd(i)-sd(i)
    if(abs(denominator)<=roundoff)then
      denominator=floor
    else
      denominator=sign(max(abs(denominator),floor),denominator)
    endif
    expected(i,:)=r(i,:)/denominator
  enddo
  call require(maxval(abs(physical-expected))<1d-8,'denominator floor lost signs or zero convention')
  ! A larger state inventory is an application extent, not a new frame.
  call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,[1d0,1d0,1d0,1d0],&
    spread(rr(:,1),2,4),output,application_fingerprint,ok,message)
  call require(ok.and.size(output,2)==4,'state-count growth required rebuilding the fixed frame')
  rr=cmplx(huge(1d0)/4d0,0d0,real64)
  call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,rr,output,&
    application_fingerprint,ok,message)
  call require(.not.ok.and..not.allocated(output),'finite action overflow accepted')
  call test_bounded_covariance()
  call test_frame_cancellation()
  call test_rectangular_entry()
  call ieee_get_halting_mode(ieee_invalid,traps_after)
  call require(traps_before.eqv.traps_after,'floating-point trap state changed')
  if(rank==0)write(*,'(a,i0,a)')'PASS fragment preconditioner on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine test_rectangular_entry()
    type(s_dg_hybrid_fragment_preconditioner)::frame_cache,square_cache
    type(s_dg_hybrid_preconditioner_key)::stale
    complex(real64)::fq(2,4),fh(2,2),fs(2,2),res(2,1),oracle(2,1)
    complex(real64)::v(2,2),identity(2,2),rotq(2,4),roth(2,2),rots(2,2),rotr(2,1)
    complex(real64),allocatable::square_answer(:,:)
    complex(real64),allocatable::answer(:,:)
    integer(int64)::fp
    integer::a
    real(real64)::dh,ds
    fq=0d0;fq(1,1:3)=sqrt(2d0/3d0)*[1d0,-0.5d0,-0.5d0]
    fq(2,1:3)=[0d0,sqrt(0.5d0),-sqrt(0.5d0)]
    fh=0d0;fs=0d0;fh(1,1)=2d0;fh(2,2)=4d0;fs(1,1)=1d0;fs(2,2)=2d0
    res(:,1)=[1d0,2d0];oracle=0d0
    do a=1,3
      dh=real(dot_product(fq(:,a),matmul(fh,fq(:,a))),real64)
      ds=real(dot_product(fq(:,a),matmul(fs,fq(:,a))),real64)
      oracle(:,1)=oracle(:,1)+fq(:,a)*dot_product(fq(:,a),res(:,1))/(dh-0.3d0*ds)
    enddo
    call prepare_dg_hybrid_frame_preconditioner(MPI_COMM_SELF,2,[1_int64,2_int64],fq,fh,fs,&
      key,701_int64,1d-10,frame_cache,fp,ok,message)
    call require(ok,'rectangular frame preparation: '//trim(message))
    call apply_dg_hybrid_fragment_preconditioner(MPI_COMM_SELF,[1_int64,2_int64],key,frame_cache,&
      [0.3d0],res,answer,fp,ok,message,selection_fingerprint=701_int64)
    call require(ok,'rectangular frame application: '//trim(message))
    call require(maxval(abs(answer-oracle))<1d-10,'rectangular action differs from physical oracle')
    fq(1,4)=1d-8
    call prepare_dg_hybrid_frame_preconditioner(MPI_COMM_SELF,2,[1_int64,2_int64],fq,fh,fs,&
      key,701_int64,1d-10,frame_cache,fp,ok,message)
    call require(.not.ok.and.index(message,'metric norm')>0,'unresolved nonzero reference column accepted')
    fq(1,4)=0d0
    call apply_dg_hybrid_fragment_preconditioner(MPI_COMM_SELF,[1_int64,2_int64],key,frame_cache,&
      [0.3d0],res,answer,fp,ok,message,selection_fingerprint=701_int64)
    call require(ok,'failed frame rebuild destroyed previous cache')
    call require(maxval(abs(answer-oracle))<1d-10,'failed rebuild changed previous frame action')
    call apply_dg_hybrid_fragment_preconditioner(MPI_COMM_SELF,[1_int64,2_int64],key,frame_cache,&
      [0.3d0],res,answer,fp,ok,message)
    call require(.not.ok.and..not.allocated(answer),'missing selection receipt accepted')
    stale=key;stale%operator_epoch=stale%operator_epoch+1
    call apply_dg_hybrid_fragment_preconditioner(MPI_COMM_SELF,[1_int64,2_int64],stale,frame_cache,&
      [0.3d0],res,answer,fp,ok,message,selection_fingerprint=701_int64)
    call require(.not.ok.and..not.allocated(answer),'stale rectangular operator epoch accepted')
    v(1,:)=[cmplx(1d0,0d0,real64),cmplx(0d0,1d0,real64)]/sqrt(2d0)
    v(2,:)=[cmplx(0d0,1d0,real64),cmplx(1d0,0d0,real64)]/sqrt(2d0)
    rotq=matmul(conjg(transpose(v)),fq)
    roth=matmul(conjg(transpose(v)),matmul(fh,v));rots=matmul(conjg(transpose(v)),matmul(fs,v))
    rotr=matmul(conjg(transpose(v)),res)
    call prepare_dg_hybrid_frame_preconditioner(MPI_COMM_SELF,2,[1_int64,2_int64],rotq,roth,rots,&
      key,701_int64,1d-10,frame_cache,fp,ok,message)
    call require(ok,'rotated rectangular frame rejected: '//trim(message))
    call apply_dg_hybrid_fragment_preconditioner(MPI_COMM_SELF,[1_int64,2_int64],key,frame_cache,&
      [0.3d0],rotr,answer,fp,ok,message,selection_fingerprint=701_int64)
    call require(ok,'rotated rectangular action rejected')
    call require(maxval(abs(matmul(v,answer)-oracle))<1d-10,'rectangular physical covariance failed')
    identity=0d0;identity(1,1)=1d0;identity(2,2)=1d0
    call prepare_dg_hybrid_fragment_preconditioner(MPI_COMM_SELF,2,[1_int64,2_int64],identity,fh,fs,&
      key,1d-10,square_cache,fp,ok,message)
    call require(ok,'square reference preparation failed')
    call apply_dg_hybrid_fragment_preconditioner(MPI_COMM_SELF,[1_int64,2_int64],key,square_cache,&
      [0.3d0],res,square_answer,fp,ok,message)
    call require(ok,'square reference application failed')
    call prepare_dg_hybrid_frame_preconditioner(MPI_COMM_SELF,2,[1_int64,2_int64],identity,fh,fs,&
      key,701_int64,1d-10,frame_cache,fp,ok,message)
    call require(ok,'square-limit frame preparation failed')
    call apply_dg_hybrid_fragment_preconditioner(MPI_COMM_SELF,[1_int64,2_int64],key,frame_cache,&
      [0.3d0],res,answer,fp,ok,message,selection_fingerprint=701_int64)
    call require(ok,'square-limit frame application failed')
    call require(maxval(abs(answer-square_answer))<1d-12,'rectangular square limit changed old action')
    ! Matrix validation must not introduce positivity tests in an arbitrary
    ! identity frame: only the actual reference Q defines these norms.
    v(1,:)=[1d0,1d0]/sqrt(2d0);v(2,:)=[1d0,-1d0]/sqrt(2d0)
    fh=identity;fs=identity;fs(1,1)=1d-9
    call prepare_dg_hybrid_fragment_preconditioner(MPI_COMM_SELF,2,[1_int64,2_int64],v,fh,fs,&
      key,1d-6,square_cache,fp,ok,message)
    call require(ok,'square anisotropic reference failed')
    call prepare_dg_hybrid_frame_preconditioner(MPI_COMM_SELF,2,[1_int64,2_int64],v,fh,fs,&
      key,701_int64,1d-6,frame_cache,fp,ok,message)
    call require(ok,'identity validation imposed spurious reference norm rejection: '//trim(message))
    call apply_dg_hybrid_fragment_preconditioner(MPI_COMM_SELF,[1_int64,2_int64],key,frame_cache,&
      [0.3d0],res,answer,fp,ok,message,selection_fingerprint=702_int64)
    call require(.not.ok.and..not.allocated(answer),'changed frame selection accepted')
    fh(1,1)=1d0;fh(2,2)=-1d0;fs=identity;res(:,1)=[1d0,0d0]
    call prepare_dg_hybrid_frame_preconditioner(MPI_COMM_SELF,2,[1_int64,2_int64],fq,fh,fs,&
      key,701_int64,1d-10,frame_cache,fp,ok,message)
    call require(ok,'signed frame preparation failed')
    call apply_dg_hybrid_fragment_preconditioner(MPI_COMM_SELF,[1_int64,2_int64],key,frame_cache,&
      [0d0],res,answer,fp,ok,message,selection_fingerprint=701_int64)
    call require(.not.ok.and..not.allocated(answer).and.index(message,'cancellation')>0,&
      'rectangular entry published an annihilated residual')
    fq(1,:)=0d0
    call prepare_dg_hybrid_frame_preconditioner(MPI_COMM_SELF,2,[1_int64,2_int64],fq,fh,fs,&
      key,701_int64,1d-10,frame_cache,fp,ok,message)
    call require(.not.ok,'rank-deficient rectangular frame accepted')
    if(nproc>1)then
      if(rank==0)then
        call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,rr,answer,fp,ok,message,&
          selection_fingerprint=701_int64)
      else
        call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,rr,answer,fp,ok,message)
      endif
      call require(.not.ok.and..not.allocated(answer),'rank-disagreeing selection presence accepted')
    endif
  end subroutine
  subroutine test_frame_cancellation()
    complex(real64)::frame(2,3),rotated(2,3),v(2,2),amplitudes(3,1)
    complex(real64),allocatable::local_frame(:,:)
    real(real64)::diagonal(3),tol
    integer::p,rows
    frame(1,:)=sqrt(2d0/3d0)*[1d0,-0.5d0,-0.5d0]
    frame(2,:)=[0d0,sqrt(0.5d0),-sqrt(0.5d0)]
    diagonal=abs(frame(1,:))**2-abs(frame(2,:))**2
    rows=0;if(rank==0)rows=2
    allocate(local_frame(rows,3));local_frame=frame(:rows,:)
    amplitudes(:,1)=conjg(frame(1,:))/diagonal
    call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,1d-10,ok,message)
    call require(.not.ok.and.index(message,'cancellation')>0,'signed rectangular cancellation was accepted')
    v(1,:)=[cmplx(1d0,0d0,real64),cmplx(0d0,1d0,real64)]/sqrt(2d0)
    v(2,:)=[cmplx(0d0,1d0,real64),cmplx(1d0,0d0,real64)]/sqrt(2d0)
    rotated=matmul(conjg(transpose(v)),frame);local_frame=rotated(:rows,:)
    call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,1d-10,ok,message)
    call require(.not.ok,'rotated cancellation was accepted')
    amplitudes(:,1)=(conjg(frame(1,:))+1d-12*conjg(frame(2,:)))/diagonal
    call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,1d-10,ok,message)
    call require(.not.ok,'near-cancellation below tolerance was accepted')
    amplitudes(:,1)=(conjg(frame(1,:))+1d-6*conjg(frame(2,:)))/diagonal
    call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,1d-10,ok,message)
    call require(ok,'resolved action above cancellation tolerance was rejected')
    do p=1,2
      if(p==1)local_frame=frame(:rows,:)
      if(p==2)local_frame=rotated(:rows,:)
      amplitudes(:,1)=conjg(frame(2,:))/diagonal
      call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,1d-10,ok,message)
      call require(ok,'resolved signed action rejected: '//trim(message))
    enddo
    amplitudes=0d0
    call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,1d-10,ok,message)
    call require(ok,'zero residual was treated as cancellation')
    amplitudes(:,1)=conjg(frame(2,:))/diagonal*1d200
    call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,1d-10,ok,message)
    call require(ok,'finite scaled action overflowed the cancellation diagnostic')
    amplitudes(:,1)=conjg(frame(2,:))/diagonal*1d-200
    call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,1d-10,ok,message)
    call require(ok,'tiny finite action underflowed the cancellation diagnostic')
    if(nproc>1)then
      tol=1d-10;if(rank==0)tol=2d-10
      call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,tol,ok,message)
      call require(.not.ok,'different cancellation controls accepted')
      if(rank==0)amplitudes(1,1)=1d0
      call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,1d-10,ok,message)
      call require(.not.ok,'different replicated amplitudes accepted')
    endif
    amplitudes=0d0
    if(rank==0)amplitudes(1,1)=ieee_value(0d0,ieee_quiet_nan)
    call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,1d-10,ok,message)
    call require(.not.ok,'nonfinite amplitudes accepted')
    ! Exercise nonzero partial sums on multiple ranks, not only idle ranks.
    if(nproc>1)then
      deallocate(local_frame);rows=merge(1,0,rank<2)
      allocate(local_frame(rows,3))
      if(rows==1)local_frame(1,:)=rotated(rank+1,:)
      amplitudes(:,1)=conjg(frame(1,:))/diagonal
      call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,1d-10,ok,message)
      call require(.not.ok,'distributed signed cancellation accepted')
      amplitudes(:,1)=conjg(frame(2,:))/diagonal
      call check_dg_hybrid_frame_cancellation(comm,2,local_frame,amplitudes,1d-10,ok,message)
      call require(ok,'distributed resolved sum rejected')
    endif
  end subroutine
  subroutine test_mixed_layouts()
    type(s_dg_hybrid_fragment_preconditioner)::other_cache
    integer(int64),allocatable::other_ids(:)
    integer(int64)::other_fp
    integer::row,owner,position
    allocate(other_ids(nlocal));position=0
    do row=n,1,-1
      owner=mod(row-1,min(nproc,4))
      if(owner==0)then
        owner=1
      else if(owner==1)then
        owner=0
      endif
      if(owner/=rank)cycle
      position=position+1;other_ids(position)=row
    enddo
    ids=other_ids;call distribute(q,h,s,matmul(conjg(transpose(u)),r))
    call prepare_dg_hybrid_fragment_preconditioner(comm,n,ids,qr,hr,sr,key,1d-10,&
      other_cache,other_fp,ok,message)
    call require(ok,'alternate ownership cache could not be prepared')
    if(rank==0)then
      ids=saved_ids;call distribute(q,h,s,matmul(conjg(transpose(u)),r))
      call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,rr,output,&
        application_fingerprint,ok,message)
    else
      call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,other_cache,shifts,rr,output,&
        application_fingerprint,ok,message)
    endif
    call require(.not.ok.and..not.allocated(output),'mixed ownership caches caused duplicate/missing residual rows')
    ids=saved_ids;call distribute(q,h,s,matmul(conjg(transpose(u)),r))
  end subroutine
  subroutine test_bounded_covariance()
    type(s_dg_hybrid_fragment_subspace_state)::state
    complex(real64)::href(n,n),sref(n,n),projector(n,n),reference_projector(n,n)
    real(real64)::spectrum(m),reference_spectrum(m),residual,reference_residual
    real(real64)::table(m,2),weights(m,2),density(n),reference_density(n),mu,electrons
    real(real64),allocatable::occupations(:,:)
    integer::variant,a,b,steps
    integer(int64)::workspace,receipt
    logical::advanced,converged,representative(2),tail(2)
    character(64)::reason
    href=0d0;sref=0d0
    do a=1,n
      href(a,a)=0.15d0*a;sref(a,a)=1d0+0.02d0*a
      if(a==n)cycle
      href(a,a+1)=cmplx(-0.03d0,0.01d0,real64);href(a+1,a)=conjg(href(a,a+1))
    enddo
    allocate(state%vectors(nlocal,m),state%directions(nlocal,m))
    state%fragment_id=key%fragment_id;state%basis_generation=key%basis_generation
    state%basis_fingerprint=key%basis_fingerprint;state%metric_fingerprint=key%metric_fingerprint
    state%state_count=m
    do variant=0,3
      call rotation(variant,u);q=conjg(transpose(u))
      h=matmul(conjg(transpose(u)),matmul(href,u));s=matmul(conjg(transpose(u)),matmul(sref,u))
      call distribute(q,h,s,matmul(conjg(transpose(u)),r))
      call prepare_dg_hybrid_fragment_preconditioner(comm,n,ids,qr,hr,sr,key,1d-10,cache,receipt,ok,message)
      call require(ok,'bounded covariance frame preparation: '//trim(message))
      state%vectors=rr;state%directions=0d0
      call advance_dg_hybrid_fragment_subspace(comm,n,ids,key%fragment_id,key%basis_generation,&
        key%basis_fingerprint,key%metric_fingerprint,apply_h,apply_s,maximum_steps=3,&
        intermediate_tolerance=1d-12,orthogonality_tolerance=1d-10,allowed_residual_growth=2d0,&
        state=state,eigenvalues=spectrum,iterations=steps,relative_residual=residual,&
        eigensolver_converged=converged,advanced=advanced,stop_reason=reason,&
        workspace_peak_bytes=workspace,fingerprint=receipt,ok=ok,message=message,&
        apply_shifted_preconditioner=apply_shifted)
      call require(ok.and.advanced.and.steps<=3,'bounded fixed-frame update failed: '//trim(message))
      call gather_physical(state%vectors,u,physical)
      projector=matmul(physical,conjg(transpose(physical)))
      if(variant==0)then
        reference_projector=projector;reference_spectrum=spectrum;reference_residual=residual
      endif
      call require(maxval(abs(projector-reference_projector))<1d-8.and.&
        maxval(abs(spectrum-reference_spectrum))<1d-9.and.abs(residual-reference_residual)<1d-9,&
        'three-step fixed-frame update depends on WF gauge')
      table(:,1)=spectrum;table(:,2)=[-0.2d0,0.5d0,1.5d0];weights(:,2)=0.4d0
      density=0d0
      do a=1,m
        weights(a,1)=0.6d0*sum([(real(sref(b,b),real64)*abs(physical(b,a))**2,b=1,n)])
      enddo
      representative=[rank==0,rank==mod(1,nproc)]
      call determine_dc_fragment_occupations(comm,table,weights,representative,0d0,2d0,2d0,1d-8,&
        mu,occupations,electrons,ok,message,tail)
      call require(ok.and..not.any(tail),'common-mu fixed-frame fixture failed: '//trim(message))
      do a=1,m
        do b=1,n
          density(b)=density(b)+0.6d0*occupations(a,1)*real(sref(b,b),real64)*abs(physical(b,a))**2
        enddo
      enddo
      density(:m)=density(:m)+0.4d0*occupations(:,2)
      if(variant==0)reference_density=density
      call require(abs(sum(density)-2d0)<1d-8.and.maxval(abs(density-reference_density))<1d-8,&
        'common-mu density depends on fixed-frame WF gauge')
    enddo
  end subroutine
  subroutine apply_h(input,output,valid)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    call matrix_action(h,input,output);valid=.true.
  end subroutine
  subroutine apply_s(input,output,valid)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    call matrix_action(s,input,output);valid=.true.
  end subroutine
  subroutine matrix_action(matrix,input,output)
    complex(real64),intent(in)::matrix(:,:),input(:,:)
    complex(real64),intent(out)::output(:,:)
    complex(real64)::all_rows(n,size(input,2))
    integer::a
    all_rows=0d0
    do a=1,nlocal;all_rows(int(ids(a)),:)=input(a,:);enddo
    call MPI_Allreduce(MPI_IN_PLACE,all_rows,size(all_rows),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do a=1,nlocal;output(a,:)=matmul(matrix(int(ids(a)),:),all_rows);enddo
  end subroutine
  subroutine apply_shifted(input,shifts,output,valid)
    complex(real64),intent(in)::input(:,:)
    real(real64),intent(in)::shifts(:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    complex(real64),allocatable::preconditioned(:,:)
    integer(int64)::receipt
    character(256)::diagnostic
    call apply_dg_hybrid_fragment_preconditioner(comm,ids,key,cache,shifts,input,preconditioned,receipt,valid,diagnostic)
    output=0d0
    if(valid)output=preconditioned
  end subroutine
  subroutine rotation(variant,u)
    integer,intent(in)::variant
    complex(real64),intent(out)::u(n,n)
    integer::a,b
    real(real64)::angle
    u=0d0
    do a=1,n;u(a,a)=1d0;enddo
    select case(variant)
    case(1)
      do a=1,nw;u(a,a)=cmplx(cos(0.3d0*a),sin(0.3d0*a),real64);enddo
    case(2)
      u(:nw,:nw)=0d0
      do a=1,nw;u(a,nw-a+1)=1d0;enddo
    case(3)
      do a=1,nw
        do b=1,nw
          angle=2d0*acos(-1d0)*real((a-1)*(b-1),real64)/nw
          u(a,b)=cmplx(cos(angle),sin(angle),real64)/sqrt(real(nw,real64))
        enddo
      enddo
    end select
  end subroutine
  subroutine distribute(q,h,s,r)
    complex(real64),intent(in)::q(n,n),h(n,n),s(n,n),r(n,m)
    integer::a
    do a=1,nlocal
      qr(a,:)=q(int(ids(a)),:);hr(a,:)=h(int(ids(a)),:);sr(a,:)=s(int(ids(a)),:);rr(a,:)=r(int(ids(a)),:)
    enddo
  end subroutine
  subroutine gather_physical(local,u,physical)
    complex(real64),intent(in)::local(:,:),u(n,n)
    complex(real64),intent(out)::physical(n,m)
    complex(real64)::all_rows(n,m)
    integer::a
    all_rows=0d0
    do a=1,nlocal;all_rows(int(ids(a)),:)=local(a,:);enddo
    call MPI_Allreduce(MPI_IN_PLACE,all_rows,n*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    physical=matmul(u,all_rows)
  end subroutine
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::bad
    call MPI_Allreduce(merge(0,1,condition),bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(bad/=0)then
      if(rank==0)write(*,'(a)')trim(label)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine
end program
