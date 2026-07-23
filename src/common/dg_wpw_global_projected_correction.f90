module dg_wpw_global_projected_correction
  use mpi,only:MPI_Allreduce,MPI_INTEGER,MPI_MAX,MPI_MIN,MPI_SUCCESS,MPI_DOUBLE_PRECISION
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private

  integer,parameter,public::DG_WPW_GLOBAL_STATE_ZERO_RHS=1
  integer,parameter,public::DG_WPW_GLOBAL_STATE_CONVERGED=2
  integer,parameter,public::DG_WPW_GLOBAL_STATE_NONCONVERGED=3
  integer,parameter,public::DG_WPW_GLOBAL_STATE_CALLBACK_FAILURE=4
  integer,parameter,public::DG_WPW_GLOBAL_STATE_STALE_OPERATOR=5
  integer,parameter,public::DG_WPW_GLOBAL_STATE_INVALID_RESULT=6
  integer,parameter,public::DG_WPW_GLOBAL_BREAKDOWN_NONE=0
  integer,parameter,public::DG_WPW_GLOBAL_BREAKDOWN_HAPPY=1
  integer,parameter,public::DG_WPW_GLOBAL_BREAKDOWN_FAILED=2

  type,public::s_dg_wpw_global_correction_controls
    integer::restart=8
    integer::max_iterations=32
    integer::state_batch=8
    real(8)::relative_tolerance=1d-2
  end type

  type,public::s_dg_wpw_global_correction_diagnostics
    real(8),allocatable::initial_residual(:),final_residual(:)
    real(8),allocatable::relative_residual(:),s_orthogonality(:)
    real(8),allocatable::equation_defect(:),projected_fraction(:)
    real(8),allocatable::correction_norm(:),amplification(:)
    integer,allocatable::iterations(:),restart_count(:)
    integer,allocatable::breakdown_status(:),state_status(:)
    logical,allocatable::zero_rhs(:),converged(:)
  end type

  abstract interface
    subroutine apply_batch_interface(context,xw,xp,yw,yp,info)
      class(*),intent(inout)::context
      complex(8),intent(in)::xw(:,:),xp(:,:)
      complex(8),intent(out)::yw(:,:),yp(:,:)
      integer,intent(out)::info
    end subroutine
    subroutine global_gram_interface(left,right,nrow,nleft,nright,gram,info)
      integer,intent(in)::nrow,nleft,nright
      complex(8),intent(in)::left(nrow,nleft),right(nrow,nright)
      complex(8),intent(out)::gram(nleft,nright)
      integer,intent(out)::info
    end subroutine
    subroutine validate_operator_state_interface(context,expected_epoch,expected_fingerprint,info)
      class(*),intent(inout)::context
      integer,intent(in)::expected_epoch
      integer(8),intent(in)::expected_fingerprint
      integer,intent(out)::info
    end subroutine
  end interface
  interface
    subroutine zlartg(f,g,c,s,r)
      complex(8),intent(in)::f,g
      real(8),intent(out)::c
      complex(8),intent(out)::s,r
    end subroutine
  end interface

  public::solve_dg_wpw_global_projected_correction
  public::release_dg_wpw_global_correction_diagnostics
contains

  subroutine solve_dg_wpw_global_projected_correction(context,comm,apply_h,apply_s,global_gram,&
    validate_operator_state,expected_epoch,expected_fingerprint,qw,qp,eigenvalues,rw,rp,controls,&
    zw,zp,diagnostics,info)
    class(*),intent(inout)::context
    integer,intent(in)::comm,expected_epoch
    integer(8),intent(in)::expected_fingerprint
    procedure(apply_batch_interface)::apply_h,apply_s
    procedure(global_gram_interface)::global_gram
    procedure(validate_operator_state_interface)::validate_operator_state
    complex(8),intent(in)::qw(:,:),qp(:,:),rw(:,:),rp(:,:)
    real(8),intent(in)::eigenvalues(:)
    type(s_dg_wpw_global_correction_controls),intent(in)::controls
    complex(8),intent(out)::zw(:,:),zp(:,:)
    type(s_dg_wpw_global_correction_diagnostics),intent(inout)::diagnostics
    integer,intent(out)::info
    type(s_dg_wpw_global_correction_diagnostics)::candidate
    complex(8),allocatable::sqw(:,:),sqp(:,:)
    integer::nstate,nq,first,last,astat,local_bad,global_bad,ierr,state_info

    call release_dg_wpw_global_correction_diagnostics(diagnostics)
    info=1;zw=(0d0,0d0);zp=(0d0,0d0);nstate=size(eigenvalues);nq=size(qw,2)
    local_bad=0
    if(controls%restart<1.or.controls%restart>16.or.controls%max_iterations<1.or.&
      controls%max_iterations>64.or.controls%state_batch<1.or.controls%state_batch>16.or.&
      .not.ieee_is_finite(controls%relative_tolerance).or.controls%relative_tolerance<=0d0.or.&
      controls%relative_tolerance>=1d0.or.nstate<1.or.nq<1)local_bad=1
    if(size(qp,2)/=nq.or.size(rw,2)/=nstate.or.size(rp,2)/=nstate.or.&
      size(zw,2)/=nstate.or.size(zp,2)/=nstate.or.size(qw,1)/=size(rw,1).or.&
      size(qp,1)/=size(rp,1).or.any(shape(zw)/=shape(rw)).or.any(shape(zp)/=shape(rp)))local_bad=1
    if(.not.finite2(qw).or..not.finite2(qp).or..not.finite2(rw).or..not.finite2(rp).or.&
      .not.all(ieee_is_finite(eigenvalues)))local_bad=1
    if(nstate>huge(1)/17.or.size(rw,1)>huge(1)/max(1,nstate).or.&
      size(rp,1)>huge(1)/max(1,nstate))local_bad=1
    if(size(rw,1)>huge(1)-size(rp,1).or.size(qw,1)>huge(1)-size(qp,1).or.&
      nq>huge(1)/16.or.size(rw,1)>huge(1)/(17*16).or.&
      size(rp,1)>huge(1)/(17*16))local_bad=1
    call MPI_Allreduce(nstate,astat,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.astat/=nstate)local_bad=1
    call MPI_Allreduce(nstate,astat,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.astat/=nstate)local_bad=1
    call MPI_Allreduce(nq,astat,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.astat/=nq)local_bad=1
    call MPI_Allreduce(nq,astat,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.astat/=nq)local_bad=1
    call collective_bad(comm,local_bad,global_bad)
    if(global_bad/=0)return
    call allocate_diagnostics(candidate,nstate,astat)
    call collective_bad(comm,merge(0,1,astat==0),global_bad)
    if(global_bad/=0)return
    allocate(sqw(size(qw,1),nq),sqp(size(qp,1),nq),stat=astat)
    call collective_bad(comm,merge(0,1,astat==0),global_bad)
    if(global_bad/=0)then;call move_diagnostics(candidate,diagnostics);return;endif
    call checked_apply(context,comm,apply_s,validate_operator_state,expected_epoch,&
      expected_fingerprint,qw,qp,sqw,sqp,state_info)
    if(state_info/=0)then
      candidate%state_status=merge(DG_WPW_GLOBAL_STATE_STALE_OPERATOR,&
        DG_WPW_GLOBAL_STATE_CALLBACK_FAILURE,state_info==2)
      call synchronize_diagnostics(comm,candidate,astat)
      call move_diagnostics(candidate,diagnostics);return
    endif

    do first=1,nstate,controls%state_batch
      last=min(nstate,first+controls%state_batch-1)
      ! Deterministic contiguous state batching; each state has an independent
      ! shifted operator but no callback order depends on rank-local data.
      call solve_state_batch(context,comm,apply_h,apply_s,global_gram,validate_operator_state,&
        expected_epoch,expected_fingerprint,qw,qp,sqw,sqp,eigenvalues(first:last),rw(:,first:last),&
        rp(:,first:last),controls,zw(:,first:last),zp(:,first:last),candidate,first,state_info)
      if(state_info/=0)then
        zw=0;zp=0
        call synchronize_diagnostics(comm,candidate,astat)
        call move_diagnostics(candidate,diagnostics);return
      endif
    enddo
    call move_diagnostics(candidate,diagnostics);info=0
  end subroutine

  subroutine solve_state_batch(context,comm,apply_h,apply_s,global_gram,validate_operator_state,&
    expected_epoch,expected_fingerprint,qw,qp,sqw,sqp,shifts,rw,rp,controls,zw,zp,diag,first_state,info)
    class(*),intent(inout)::context
    integer,intent(in)::comm,expected_epoch,first_state
    integer(8),intent(in)::expected_fingerprint
    procedure(apply_batch_interface)::apply_h,apply_s
    procedure(global_gram_interface)::global_gram
    procedure(validate_operator_state_interface)::validate_operator_state
    complex(8),intent(in)::qw(:,:),qp(:,:),sqw(:,:),sqp(:,:),rw(:,:),rp(:,:)
    real(8),intent(in)::shifts(:)
    type(s_dg_wpw_global_correction_controls),intent(in)::controls
    complex(8),intent(out)::zw(:,:),zp(:,:)
    type(s_dg_wpw_global_correction_diagnostics),intent(inout)::diag
    integer,intent(out)::info
    complex(8),allocatable::bw(:,:),bp(:,:),resw(:,:),resp(:,:),aw(:,:),ap(:,:),&
      workw(:,:),workp(:,:),vw(:,:,:),vp(:,:,:),hess(:,:,:),g(:,:),y(:,:),&
      gram(:,:),qgram(:,:),sn(:,:)
    complex(8)::temp
    real(8),allocatable::cs(:,:),raw_norm(:),initial_norm(:),residual_norm(:),&
      arnoldi_norm(:),target(:),absolute_floor(:),z_norm(:)
    integer,allocatable::used(:),remaining(:),stop_mask(:),mask_min(:),mask_max(:),decision(:)
    logical,allocatable::active(:),cycle_active(:),breakdown_event(:)
    integer::nb,astat,cycle,j,k,s,idx,max_steps,apply_info
    real(8)::breakdown_floor,hscale

    info=1;zw=0;zp=0;nb=size(shifts)
    allocate(bw(size(rw,1),nb),bp(size(rp,1),nb),resw(size(rw,1),nb),resp(size(rp,1),nb),&
      aw(size(rw,1),nb),ap(size(rp,1),nb),workw(size(rw,1),nb),workp(size(rp,1),nb),&
      vw(size(rw,1),controls%restart+1,nb),vp(size(rp,1),controls%restart+1,nb),&
      hess(controls%restart+1,controls%restart,nb),g(controls%restart+1,nb),&
      y(controls%restart,nb),gram(nb,nb),qgram(size(qw,2),nb),&
      cs(controls%restart,nb),sn(controls%restart,nb),raw_norm(nb),initial_norm(nb),&
      residual_norm(nb),arnoldi_norm(nb),target(nb),absolute_floor(nb),z_norm(nb),&
      used(nb),remaining(nb),stop_mask(nb),mask_min(nb),mask_max(nb),active(nb),&
      cycle_active(nb),breakdown_event(nb),decision(nb),stat=astat)
    call collective_bad(comm,merge(0,1,astat==0),apply_info)
    if(apply_info/=0)then
      diag%state_status(first_state:first_state+nb-1)=DG_WPW_GLOBAL_STATE_INVALID_RESULT;return
    endif
    call vector_norms(comm,global_gram,rw,rp,raw_norm,apply_info)
    if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
    call apply_left_projector(comm,global_gram,qw,qp,sqw,sqp,-rw,-rp,bw,bp,apply_info)
    if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
    call vector_norms(comm,global_gram,bw,bp,initial_norm,apply_info)
    if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
    do s=1,nb
      idx=first_state+s-1
      diag%initial_residual(idx)=initial_norm(s)
      diag%final_residual(idx)=initial_norm(s);diag%relative_residual(idx)=1d0
      diag%equation_defect(idx)=initial_norm(s)
      diag%projected_fraction(idx)=initial_norm(s)/max(100d0*epsilon(1d0),raw_norm(s))
      absolute_floor(s)=100d0*epsilon(1d0)*max(1d0,initial_norm(s))
      target(s)=max(absolute_floor(s),controls%relative_tolerance*initial_norm(s))
      stop_mask(s)=merge(1,0,initial_norm(s)<=absolute_floor(s))
    enddo
    call agree_masks(comm,stop_mask,mask_min,mask_max,apply_info)
    if(apply_info/=0)then
      diag%state_status(first_state:first_state+nb-1)=DG_WPW_GLOBAL_STATE_INVALID_RESULT;return
    endif
    active=.true.;remaining=controls%max_iterations;resw=bw;resp=bp
    do s=1,nb
      if(stop_mask(s)==1)then
        idx=first_state+s-1;active(s)=.false.;diag%zero_rhs(idx)=.true.;diag%converged(idx)=.true.
        diag%state_status(idx)=DG_WPW_GLOBAL_STATE_ZERO_RHS;diag%final_residual(idx)=initial_norm(s)
        diag%relative_residual(idx)=0;diag%equation_defect(idx)=initial_norm(s)
        diag%s_orthogonality(idx)=0
      endif
    enddo
    cycle=0
    do while(any(active))
      cycle=cycle+1;hess=0;g=0;cs=0;sn=0;vw=0;vp=0;y=0;used=0
      cycle_active=active;breakdown_event=.false.
      do s=1,nb
        if(.not.cycle_active(s))then;resw(:,s)=0;resp(:,s)=0;endif
      enddo
      call vector_norms(comm,global_gram,resw,resp,residual_norm,apply_info)
      if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
      do s=1,nb
        if(cycle_active(s))then
          vw(:,1,s)=resw(:,s)/residual_norm(s);vp(:,1,s)=resp(:,s)/residual_norm(s)
          g(1,s)=residual_norm(s)
        endif
      enddo
      max_steps=min(controls%restart,maxval(remaining,mask=cycle_active))
      do j=1,max_steps
        workw=0;workp=0
        do s=1,nb
          if(cycle_active(s).and.remaining(s)>0)then
            workw(:,s)=vw(:,j,s);workp(:,s)=vp(:,j,s)
          endif
        enddo
        call apply_projected_operator(context,comm,apply_h,apply_s,global_gram,&
          validate_operator_state,expected_epoch,expected_fingerprint,qw,qp,sqw,sqp,shifts,&
          workw,workp,aw,ap,apply_info)
        if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
        do k=1,j
          call pair_gram(comm,global_gram,vw(:,k,:),vp(:,k,:),aw,ap,gram,apply_info)
          if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
          do s=1,nb
            if(cycle_active(s).and.remaining(s)>0)then
              hess(k,j,s)=gram(s,s);aw(:,s)=aw(:,s)-vw(:,k,s)*gram(s,s)
              ap(:,s)=ap(:,s)-vp(:,k,s)*gram(s,s)
            endif
          enddo
        enddo
        do k=1,j
          call pair_gram(comm,global_gram,vw(:,k,:),vp(:,k,:),aw,ap,gram,apply_info)
          if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
          do s=1,nb
            if(cycle_active(s).and.remaining(s)>0)then
              hess(k,j,s)=hess(k,j,s)+gram(s,s);aw(:,s)=aw(:,s)-vw(:,k,s)*gram(s,s)
              ap(:,s)=ap(:,s)-vp(:,k,s)*gram(s,s)
            endif
          enddo
        enddo
        call vector_norms(comm,global_gram,aw,ap,arnoldi_norm,apply_info)
        if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
        stop_mask=0
        do s=1,nb
          if(.not.cycle_active(s).or.remaining(s)<=0)cycle
          hess(j+1,j,s)=cmplx(arnoldi_norm(s),0d0,8)
          hscale=max(tiny(1d0),maxval(abs(hess(1:j+1,j,s))))
          breakdown_floor=100d0*epsilon(1d0)*hscale
          breakdown_event(s)=arnoldi_norm(s)<=breakdown_floor
          do k=1,j-1
            temp=cs(k,s)*hess(k,j,s)+sn(k,s)*hess(k+1,j,s)
            hess(k+1,j,s)=-conjg(sn(k,s))*hess(k,j,s)+cs(k,s)*hess(k+1,j,s)
            hess(k,j,s)=temp
          enddo
          call zlartg(hess(j,j,s),hess(j+1,j,s),cs(j,s),sn(j,s),temp)
          hess(j,j,s)=temp;hess(j+1,j,s)=0
          temp=cs(j,s)*g(j,s)+sn(j,s)*g(j+1,s)
          g(j+1,s)=-conjg(sn(j,s))*g(j,s)+cs(j,s)*g(j+1,s);g(j,s)=temp
          idx=first_state+s-1;diag%iterations(idx)=diag%iterations(idx)+1
          remaining(s)=remaining(s)-1;used(s)=j
          stop_mask(s)=merge(1,0,abs(g(j+1,s))<=target(s).or.breakdown_event(s).or.remaining(s)==0)
          if(stop_mask(s)==0)then
            vw(:,j+1,s)=aw(:,s)/arnoldi_norm(s);vp(:,j+1,s)=ap(:,s)/arnoldi_norm(s)
          endif
        enddo
        call agree_masks(comm,stop_mask,mask_min,mask_max,apply_info)
        if(apply_info/=0)then
          diag%state_status(first_state:first_state+nb-1)=DG_WPW_GLOBAL_STATE_INVALID_RESULT;return
        endif
        do s=1,nb;if(stop_mask(s)==1)cycle_active(s)=.false.;enddo
        if(.not.any(cycle_active))exit
      enddo
      stop_mask=0
      do s=1,nb
        if(.not.active(s))cycle
        idx=first_state+s-1
        call triangular_solve(hess(:,:,s),g(:,s),used(s),y(:,s),apply_info)
        if(apply_info/=0)then
          stop_mask(s)=1
        else
          do k=1,used(s)
            zw(:,s)=zw(:,s)+vw(:,k,s)*y(k,s);zp(:,s)=zp(:,s)+vp(:,k,s)*y(k,s)
          enddo
        endif
      enddo
      call agree_masks(comm,stop_mask,mask_min,mask_max,apply_info)
      if(apply_info/=0)then
        diag%state_status(first_state:first_state+nb-1)=DG_WPW_GLOBAL_STATE_INVALID_RESULT;return
      endif
      if(any(stop_mask/=0))then
        do s=1,nb
          if(stop_mask(s)/=0)then
            idx=first_state+s-1;diag%breakdown_status(idx)=DG_WPW_GLOBAL_BREAKDOWN_FAILED
            diag%state_status(idx)=DG_WPW_GLOBAL_STATE_INVALID_RESULT
          endif
        enddo
        return
      endif
      workw=0;workp=0
      do s=1,nb
        if(active(s))then;workw(:,s)=zw(:,s);workp(:,s)=zp(:,s);endif
      enddo
      call explicit_final_residual(context,comm,apply_h,apply_s,global_gram,validate_operator_state,&
        expected_epoch,expected_fingerprint,qw,qp,sqw,sqp,shifts,bw,bp,workw,workp,resw,resp,&
        residual_norm,active,apply_info)
      if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
      call vector_norms(comm,global_gram,workw,workp,z_norm,apply_info)
      if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
      decision=0
      do s=1,nb
        if(.not.active(s))cycle
        idx=first_state+s-1;diag%final_residual(idx)=residual_norm(s)
        diag%relative_residual(idx)=residual_norm(s)/max(initial_norm(s),absolute_floor(s))
        diag%equation_defect(idx)=residual_norm(s)
        diag%correction_norm(idx)=z_norm(s)
        diag%amplification(idx)=z_norm(s)/max(raw_norm(s),absolute_floor(s))
        if(breakdown_event(s))then
          if(residual_norm(s)<=target(s))then
            decision(s)=1
          else
            decision(s)=2
          endif
        elseif(residual_norm(s)<=target(s))then
          decision(s)=1
        elseif(remaining(s)<=0)then
          decision(s)=3
        else
          diag%restart_count(idx)=diag%restart_count(idx)+1
        endif
      enddo
      call agree_masks(comm,decision,mask_min,mask_max,apply_info)
      if(apply_info/=0)then
        diag%state_status(first_state:first_state+nb-1)=DG_WPW_GLOBAL_STATE_INVALID_RESULT;return
      endif
      do s=1,nb
        idx=first_state+s-1
        select case(decision(s))
        case(1)
          if(breakdown_event(s))diag%breakdown_status(idx)=DG_WPW_GLOBAL_BREAKDOWN_HAPPY
          active(s)=.false.;diag%converged(idx)=.true.;diag%state_status(idx)=DG_WPW_GLOBAL_STATE_CONVERGED
        case(2)
          diag%breakdown_status(idx)=DG_WPW_GLOBAL_BREAKDOWN_FAILED
          diag%state_status(idx)=DG_WPW_GLOBAL_STATE_NONCONVERGED
        case(3)
          diag%state_status(idx)=DG_WPW_GLOBAL_STATE_NONCONVERGED
        end select
      enddo
      if(any(decision==2).or.any(decision==3))return
      stop_mask=merge(1,0,active)
      call agree_masks(comm,stop_mask,mask_min,mask_max,apply_info)
      if(apply_info/=0)then
        diag%state_status(first_state:first_state+nb-1)=DG_WPW_GLOBAL_STATE_INVALID_RESULT;return
      endif
    enddo
    call apply_right_projector(context,comm,apply_s,global_gram,validate_operator_state,&
      expected_epoch,expected_fingerprint,qw,qp,zw,zp,workw,workp,apply_info)
    if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
    zw=workw;zp=workp
    call checked_apply(context,comm,apply_s,validate_operator_state,expected_epoch,&
      expected_fingerprint,zw,zp,workw,workp,apply_info)
    if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
    call pair_gram(comm,global_gram,qw,qp,workw,workp,qgram,apply_info)
    if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
    call vector_norms(comm,global_gram,zw,zp,z_norm,apply_info)
    if(apply_info/=0)then;call mark_batch_failure(diag,first_state,nb,apply_info);return;endif
    stop_mask=0
    do s=1,nb
      idx=first_state+s-1;diag%s_orthogonality(idx)=maxval(abs(qgram(:,s)))
      diag%correction_norm(idx)=z_norm(s)
      diag%amplification(idx)=z_norm(s)/max(raw_norm(s),absolute_floor(s))
      if(.not.finite2(zw(:,s:s)).or..not.finite2(zp(:,s:s)).or.diag%s_orthogonality(idx)>1d-11)then
        stop_mask(s)=1
      endif
    enddo
    call agree_masks(comm,stop_mask,mask_min,mask_max,apply_info)
    if(apply_info/=0)then
      diag%state_status(first_state:first_state+nb-1)=DG_WPW_GLOBAL_STATE_INVALID_RESULT;return
    endif
    if(any(stop_mask/=0))then
      do s=1,nb
        if(stop_mask(s)/=0)diag%state_status(first_state+s-1)=DG_WPW_GLOBAL_STATE_INVALID_RESULT
      enddo
      return
    endif
    info=0
  end subroutine

  subroutine apply_projected_operator(context,comm,apply_h,apply_s,global_gram,validator,&
    epoch,fingerprint,qw,qp,sqw,sqp,shift,xw,xp,yw,yp,info)
    class(*),intent(inout)::context
    integer,intent(in)::comm,epoch
    integer(8),intent(in)::fingerprint
    procedure(apply_batch_interface)::apply_h,apply_s
    procedure(global_gram_interface)::global_gram
    procedure(validate_operator_state_interface)::validator
    complex(8),intent(in)::qw(:,:),qp(:,:),sqw(:,:),sqp(:,:),xw(:,:),xp(:,:)
    real(8),intent(in)::shift(:);complex(8),intent(out)::yw(:,:),yp(:,:);integer,intent(out)::info
    complex(8),allocatable::rightw(:,:),rightp(:,:),hw(:,:),hp(:,:),sw(:,:),sp(:,:)
    integer::astat
    allocate(rightw(size(xw,1),size(xw,2)),rightp(size(xp,1),size(xp,2)),&
      hw(size(xw,1),size(xw,2)),hp(size(xp,1),size(xp,2)),&
      sw(size(xw,1),size(xw,2)),sp(size(xp,1),size(xp,2)),stat=astat)
    call collective_bad(comm,merge(0,1,astat==0),info);if(info/=0)return
    call apply_right_projector(context,comm,apply_s,global_gram,validator,epoch,fingerprint,&
      qw,qp,xw,xp,rightw,rightp,info);if(info/=0)return
    call checked_apply(context,comm,apply_h,validator,epoch,fingerprint,rightw,rightp,hw,hp,info)
    if(info/=0)return
    call checked_apply(context,comm,apply_s,validator,epoch,fingerprint,rightw,rightp,sw,sp,info)
    if(info/=0)return
    do astat=1,size(shift)
      hw(:,astat)=hw(:,astat)-shift(astat)*sw(:,astat)
      hp(:,astat)=hp(:,astat)-shift(astat)*sp(:,astat)
    enddo
    call apply_left_projector(comm,global_gram,qw,qp,sqw,sqp,hw,hp,yw,yp,info)
  end subroutine

  subroutine apply_right_projector(context,comm,apply_s,global_gram,validator,epoch,fingerprint,&
    qw,qp,xw,xp,rightw,rightp,info)
    class(*),intent(inout)::context
    integer,intent(in)::comm,epoch
    integer(8),intent(in)::fingerprint
    procedure(apply_batch_interface)::apply_s
    procedure(global_gram_interface)::global_gram
    procedure(validate_operator_state_interface)::validator
    complex(8),intent(in)::qw(:,:),qp(:,:),xw(:,:),xp(:,:)
    complex(8),intent(out)::rightw(:,:),rightp(:,:);integer,intent(out)::info
    complex(8),allocatable::sxw(:,:),sxp(:,:),global_gram_q_sx(:,:)
    integer::astat,global_bad
    allocate(sxw(size(xw,1),size(xw,2)),sxp(size(xp,1),size(xp,2)),&
      global_gram_q_sx(size(qw,2),size(xw,2)),stat=astat)
    call collective_bad(comm,merge(0,1,astat==0),global_bad)
    if(global_bad/=0)then;info=1;return;endif
    call checked_apply(context,comm,apply_s,validator,epoch,fingerprint,xw,xp,sxw,sxp,info)
    if(info/=0)return
    call pair_gram(comm,global_gram,qw,qp,sxw,sxp,global_gram_q_sx,info);if(info/=0)return
    ! right=x-matmul(q,global_gram_q_sx)
    rightw=xw-matmul(qw,global_gram_q_sx)
    rightp=xp-matmul(qp,global_gram_q_sx)
  end subroutine

  subroutine apply_left_projector(comm,global_gram,qw,qp,sqw,sqp,yw,yp,leftw,leftp,info)
    integer,intent(in)::comm
    procedure(global_gram_interface)::global_gram
    complex(8),intent(in)::qw(:,:),qp(:,:),sqw(:,:),sqp(:,:),yw(:,:),yp(:,:)
    complex(8),intent(out)::leftw(:,:),leftp(:,:);integer,intent(out)::info
    complex(8),allocatable::global_gram_q_y(:,:)
    integer::astat,global_bad
    allocate(global_gram_q_y(size(qw,2),size(yw,2)),stat=astat)
    call collective_bad(comm,merge(0,1,astat==0),global_bad)
    if(global_bad/=0)then;info=1;return;endif
    call pair_gram(comm,global_gram,qw,qp,yw,yp,global_gram_q_y,info);if(info/=0)return
    ! left=y-matmul(sq,global_gram_q_y)
    leftw=yw-matmul(sqw,global_gram_q_y)
    leftp=yp-matmul(sqp,global_gram_q_y)
  end subroutine

  subroutine explicit_final_residual(context,comm,apply_h,apply_s,global_gram,validator,epoch,&
    fingerprint,qw,qp,sqw,sqp,shift,bw,bp,zw,zp,resw,resp,norm,active,info)
    class(*),intent(inout)::context
    integer,intent(in)::comm,epoch
    integer(8),intent(in)::fingerprint
    procedure(apply_batch_interface)::apply_h,apply_s
    procedure(global_gram_interface)::global_gram
    procedure(validate_operator_state_interface)::validator
    complex(8),intent(in)::qw(:,:),qp(:,:),sqw(:,:),sqp(:,:),bw(:,:),bp(:,:),zw(:,:),zp(:,:)
    real(8),intent(in)::shift(:)
    logical,intent(in)::active(:)
    complex(8),intent(out)::resw(:,:),resp(:,:);real(8),intent(out)::norm(:);integer,intent(out)::info
    complex(8),allocatable::aw(:,:),ap(:,:)
    integer::astat,global_bad
    allocate(aw(size(zw,1),size(zw,2)),ap(size(zp,1),size(zp,2)),stat=astat)
    call collective_bad(comm,merge(0,1,astat==0),global_bad)
    if(global_bad/=0)then;info=1;return;endif
    call apply_projected_operator(context,comm,apply_h,apply_s,global_gram,validator,epoch,&
      fingerprint,qw,qp,sqw,sqp,shift,zw,zp,aw,ap,info);if(info/=0)return
    resw=bw-aw;resp=bp-ap
    do astat=1,size(active)
      if(.not.active(astat))then;resw(:,astat)=0;resp(:,astat)=0;endif
    enddo
    call vector_norms(comm,global_gram,resw,resp,norm,info)
  end subroutine

  subroutine checked_apply(context,comm,apply,validator,epoch,fingerprint,xw,xp,yw,yp,info)
    class(*),intent(inout)::context
    integer,intent(in)::comm,epoch
    integer(8),intent(in)::fingerprint
    procedure(apply_batch_interface)::apply
    procedure(validate_operator_state_interface)::validator
    complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(out)::yw(:,:),yp(:,:)
    integer,intent(out)::info
    integer::local_info,global_info
    call apply(context,xw,xp,yw,yp,local_info)
    call MPI_Allreduce(merge(0,1,local_info==0),global_info,1,MPI_INTEGER,MPI_MAX,comm,info)
    if(info/=MPI_SUCCESS.or.global_info/=0)then;info=1;return;endif
    call validator(context,epoch,fingerprint,local_info)
    call MPI_Allreduce(merge(0,1,local_info==0),global_info,1,MPI_INTEGER,MPI_MAX,comm,info)
    if(info/=MPI_SUCCESS.or.global_info/=0)then;info=2;return;endif
    info=0
  end subroutine

  subroutine pair_gram(comm,global_gram,lw,lp,rw,rp,gram,info)
    integer,intent(in)::comm
    procedure(global_gram_interface)::global_gram
    complex(8),intent(in)::lw(:,:),lp(:,:),rw(:,:),rp(:,:)
    complex(8),intent(out)::gram(:,:);integer,intent(out)::info
    complex(8),allocatable::left(:,:),right(:,:)
    integer::local_bad,global_bad,astat
    local_bad=merge(0,1,size(lw,1)<=huge(1)-size(lp,1).and.size(rw,1)<=huge(1)-size(rp,1))
    call collective_bad(comm,local_bad,global_bad);if(global_bad/=0)then;info=1;return;endif
    allocate(left(size(lw,1)+size(lp,1),size(lw,2)),right(size(rw,1)+size(rp,1),size(rw,2)),stat=astat)
    call collective_bad(comm,merge(0,1,astat==0),global_bad)
    if(global_bad/=0)then;info=1;return;endif
    left(1:size(lw,1),:)=lw;left(size(lw,1)+1:,:)=lp
    right(1:size(rw,1),:)=rw;right(size(rw,1)+1:,:)=rp
    call global_gram(left,right,size(left,1),size(left,2),size(right,2),gram,info)
    if(info/=0)then
      local_bad=1
    else
      local_bad=merge(0,1,finite2(gram))
    endif
    call collective_bad(comm,local_bad,global_bad)
    info=merge(0,1,global_bad==0)
  end subroutine

  subroutine vector_norms(comm,global_gram,w,p,norms,info)
    integer,intent(in)::comm
    procedure(global_gram_interface)::global_gram
    complex(8),intent(in)::w(:,:),p(:,:)
    real(8),intent(out)::norms(:)
    integer,intent(out)::info
    complex(8),allocatable::gram(:,:)
    integer::n,s,astat,local_bad,global_bad
    n=size(w,2);allocate(gram(n,n),stat=astat)
    call collective_bad(comm,merge(0,1,astat==0),global_bad)
    if(global_bad/=0)then;info=1;return;endif
    call pair_gram(comm,global_gram,w,p,w,p,gram,info);if(info/=0)return
    local_bad=0
    do s=1,n
      if(.not.ieee_is_finite(real(gram(s,s),8)).or.&
        abs(aimag(gram(s,s)))>1d-10*max(1d0,abs(real(gram(s,s),8))).or.&
        real(gram(s,s),8)<-1d-13)local_bad=1
    enddo
    call collective_bad(comm,local_bad,global_bad)
    if(global_bad/=0)then;info=1;return;endif
    do s=1,n;norms(s)=sqrt(max(0d0,real(gram(s,s),8)));enddo
    info=0
  end subroutine

  subroutine agree_masks(comm,mask,mask_min,mask_max,info)
    integer,intent(in)::comm,mask(:)
    integer,intent(out)::mask_min(:),mask_max(:),info
    integer::ierr
    call MPI_Allreduce(mask,mask_min,size(mask),MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;info=1;return;endif
    call MPI_Allreduce(mask,mask_max,size(mask),MPI_INTEGER,MPI_MAX,comm,ierr)
    info=merge(0,1,ierr==MPI_SUCCESS.and.all(mask_min==mask_max))
  end subroutine

  subroutine mark_batch_failure(diag,first_state,nb,failure_kind)
    type(s_dg_wpw_global_correction_diagnostics),intent(inout)::diag
    integer,intent(in)::first_state,nb,failure_kind
    integer::status
    status=merge(DG_WPW_GLOBAL_STATE_STALE_OPERATOR,DG_WPW_GLOBAL_STATE_CALLBACK_FAILURE,&
      failure_kind==2)
    diag%state_status(first_state:first_state+nb-1)=status
  end subroutine

  subroutine synchronize_diagnostics(comm,diag,info)
    integer,intent(in)::comm
    type(s_dg_wpw_global_correction_diagnostics),intent(inout)::diag
    integer,intent(out)::info
    real(8),allocatable::rwork(:)
    integer,allocatable::iwork(:),logical_as_int(:)
    integer::n,i,ierr,astat,global_bad
    info=1;if(.not.allocated(diag%state_status))return
    n=size(diag%state_status);allocate(rwork(n),iwork(n),logical_as_int(n),stat=astat)
    call collective_bad(comm,merge(0,1,astat==0),global_bad);if(global_bad/=0)return
    call sync_real(diag%initial_residual);call sync_real(diag%final_residual)
    call sync_real(diag%relative_residual);call sync_real(diag%s_orthogonality)
    call sync_real(diag%equation_defect);call sync_real(diag%projected_fraction)
    call sync_real(diag%correction_norm);call sync_real(diag%amplification)
    call sync_int(diag%iterations);call sync_int(diag%restart_count)
    call sync_int(diag%breakdown_status);call sync_int(diag%state_status)
    if(info==2)return
    do i=1,n;logical_as_int(i)=merge(1,0,diag%zero_rhs(i));enddo
    call MPI_Allreduce(logical_as_int,iwork,n,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    diag%zero_rhs=iwork/=0
    do i=1,n;logical_as_int(i)=merge(1,0,diag%converged(i));enddo
    call MPI_Allreduce(logical_as_int,iwork,n,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    diag%converged=iwork/=0
    info=0
  contains
    subroutine sync_real(values)
      real(8),intent(inout)::values(:)
      if(info==2)return
      call MPI_Allreduce(values,rwork,n,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;info=2;return;endif
      values=rwork
    end subroutine
    subroutine sync_int(values)
      integer,intent(inout)::values(:)
      if(info==2)return
      call MPI_Allreduce(values,iwork,n,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;info=2;return;endif
      values=iwork
    end subroutine
  end subroutine

  subroutine triangular_solve(hess,g,n,y,info)
    complex(8),intent(in)::hess(:,:),g(:)
    integer,intent(in)::n
    complex(8),intent(out)::y(:)
    integer,intent(out)::info
    integer::i
    y=0;info=1
    do i=n,1,-1
      if(abs(hess(i,i))<=100d0*epsilon(1d0)*&
        max(tiny(1d0),maxval(abs(hess(1:i,i)))))return
      y(i)=(g(i)-sum(hess(i,i+1:n)*y(i+1:n)))/hess(i,i)
    enddo
    info=0
  end subroutine

  subroutine allocate_diagnostics(diag,n,info)
    type(s_dg_wpw_global_correction_diagnostics),intent(out)::diag
    integer,intent(in)::n;integer,intent(out)::info
    integer::astat
    allocate(diag%initial_residual(n),diag%final_residual(n),diag%relative_residual(n),&
      diag%s_orthogonality(n),diag%equation_defect(n),diag%projected_fraction(n),&
      diag%correction_norm(n),diag%amplification(n),diag%iterations(n),diag%restart_count(n),&
      diag%breakdown_status(n),diag%state_status(n),diag%zero_rhs(n),diag%converged(n),stat=astat)
    info=astat;if(astat/=0)return
    diag%initial_residual=0;diag%final_residual=huge(1d0);diag%relative_residual=huge(1d0)
    diag%s_orthogonality=huge(1d0);diag%equation_defect=huge(1d0);diag%projected_fraction=0
    diag%correction_norm=0;diag%amplification=0;diag%iterations=0;diag%restart_count=0
    diag%breakdown_status=DG_WPW_GLOBAL_BREAKDOWN_NONE
    diag%state_status=DG_WPW_GLOBAL_STATE_INVALID_RESULT
    diag%zero_rhs=.false.;diag%converged=.false.
  end subroutine

  subroutine move_diagnostics(source,destination)
    type(s_dg_wpw_global_correction_diagnostics),intent(inout)::source,destination
    call release_dg_wpw_global_correction_diagnostics(destination)
    call move_alloc(source%initial_residual,destination%initial_residual)
    call move_alloc(source%final_residual,destination%final_residual)
    call move_alloc(source%relative_residual,destination%relative_residual)
    call move_alloc(source%s_orthogonality,destination%s_orthogonality)
    call move_alloc(source%equation_defect,destination%equation_defect)
    call move_alloc(source%projected_fraction,destination%projected_fraction)
    call move_alloc(source%correction_norm,destination%correction_norm)
    call move_alloc(source%amplification,destination%amplification)
    call move_alloc(source%iterations,destination%iterations)
    call move_alloc(source%restart_count,destination%restart_count)
    call move_alloc(source%breakdown_status,destination%breakdown_status)
    call move_alloc(source%state_status,destination%state_status)
    call move_alloc(source%zero_rhs,destination%zero_rhs)
    call move_alloc(source%converged,destination%converged)
  end subroutine

  subroutine release_dg_wpw_global_correction_diagnostics(diag)
    type(s_dg_wpw_global_correction_diagnostics),intent(inout)::diag
    if(allocated(diag%initial_residual))deallocate(diag%initial_residual)
    if(allocated(diag%final_residual))deallocate(diag%final_residual)
    if(allocated(diag%relative_residual))deallocate(diag%relative_residual)
    if(allocated(diag%s_orthogonality))deallocate(diag%s_orthogonality)
    if(allocated(diag%equation_defect))deallocate(diag%equation_defect)
    if(allocated(diag%projected_fraction))deallocate(diag%projected_fraction)
    if(allocated(diag%correction_norm))deallocate(diag%correction_norm)
    if(allocated(diag%amplification))deallocate(diag%amplification)
    if(allocated(diag%iterations))deallocate(diag%iterations)
    if(allocated(diag%restart_count))deallocate(diag%restart_count)
    if(allocated(diag%breakdown_status))deallocate(diag%breakdown_status)
    if(allocated(diag%state_status))deallocate(diag%state_status)
    if(allocated(diag%zero_rhs))deallocate(diag%zero_rhs)
    if(allocated(diag%converged))deallocate(diag%converged)
  end subroutine

  subroutine collective_bad(comm,local_bad,global_bad)
    integer,intent(in)::comm,local_bad;integer,intent(out)::global_bad
    integer::ierr
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)global_bad=1
  end subroutine

  logical function finite2(a)
    complex(8),intent(in)::a(:,:)
    finite2=all(ieee_is_finite(real(a,8))).and.all(ieee_is_finite(aimag(a)))
  end function
end module
