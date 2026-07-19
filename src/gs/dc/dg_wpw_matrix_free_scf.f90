module dg_wpw_matrix_free_scf
  use mpi,only:MPI_Allreduce,MPI_Comm_rank,MPI_INTEGER,MPI_MAX,MPI_MIN,MPI_DOUBLE_PRECISION,MPI_SUCCESS
  use dg_wpw_matrix_free_operator,only:apply_h_batch,apply_s_batch,global_gram_batch
  use dg_generalized_algebra,only:dg_generalized_eigh,dg_reduced_generalized_eigh,dg_metric_factor
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private

  type,public::s_dg_wpw_matrix_free_scf_result
    logical::converged=.false.
    integer::iterations=0,retained_dimension=0
    real(8)::density_residual=huge(1d0),potential_residual=huge(1d0)
    real(8)::energy_residual=huge(1d0),projector_residual=huge(1d0)
    real(8)::generalized_residual=huge(1d0),metric_orthonormality=huge(1d0)
    real(8)::charge_error=huge(1d0),gap=0d0,fixed_point_residual=huge(1d0)
    real(8)::total_energy=huge(1d0)
  end type

  abstract interface
    subroutine dg_wpw_scf_map(context,cw,cp,occ,n_occ,density_in,mix,density_out,&
        potential_residual,total_energy,charge_error,info)
      class(*),intent(inout)::context
      integer,intent(in)::n_occ
      complex(8),intent(in)::cw(:,:),cp(:,:)
      real(8),intent(in)::occ(n_occ),density_in(:),mix
      real(8),intent(out)::density_out(:),potential_residual,total_energy,charge_error
      integer,intent(out)::info
    end subroutine
    subroutine dg_wpw_preconditioner(context,eigenvalues,rw,rp,zw,zp,info)
      class(*),intent(inout)::context
      real(8),intent(in)::eigenvalues(:)
      complex(8),intent(in)::rw(:,:),rp(:,:)
      complex(8),intent(out)::zw(:,:),zp(:,:)
      integer,intent(out)::info
    end subroutine
  end interface

  public::run_dg_wpw_matrix_free_scf,dg_wpw_scf_map
  public::run_dg_wpw_matrix_free_algebra_step,dg_wpw_preconditioner
  public::initialize_dg_wpw_deterministic_seed
  public::initialize_dg_wpw_projected_occupied
  public::complete_dg_wpw_projected_subspace
  public::solve_dg_wpw_metric_projection
  public::initialize_dg_wpw_metric_projected_occupied
contains
  subroutine initialize_dg_wpw_metric_projected_occupied(context,comm,apply_s,global_gram,nw,np,nocc,&
      tolerance,max_iterations,diagonal_w,diagonal_p,bw,bp,qw,qp,relative_residual,rhs_residuals,&
      rhs_residual_history,iterations,diagonal_spread,effective_rank,orthogonality,info)
    class(*),intent(inout)::context
    integer,intent(in)::comm,nw,np,nocc,max_iterations
    procedure(apply_s_batch)::apply_s
    procedure(global_gram_batch)::global_gram
    real(8),intent(in)::tolerance,diagonal_w(nw),diagonal_p(np)
    complex(8),intent(in)::bw(nw,nocc),bp(np,nocc)
    complex(8),intent(out)::qw(nw,nocc),qp(np,nocc)
    real(8),intent(out)::relative_residual,rhs_residuals(nocc),rhs_residual_history(max_iterations,nocc),&
      diagonal_spread,orthogonality
    integer,intent(out)::iterations,effective_rank,info

    qw=(0d0,0d0);qp=(0d0,0d0);effective_rank=0;orthogonality=huge(1d0)
    call solve_dg_wpw_metric_projection(context,comm,apply_s,global_gram,nw,np,nocc,&
      tolerance,max_iterations,diagonal_w,diagonal_p,bw,bp,qw,qp,relative_residual,rhs_residuals,&
      rhs_residual_history,iterations,diagonal_spread,info)
    if(info/=0)return
    call initialize_dg_wpw_projected_occupied(context,comm,apply_s,global_gram,nw,np,nocc,&
      tolerance,qw,qp,effective_rank,orthogonality,info)
    if(info/=0.or.effective_rank/=nocc)then
      qw=(0d0,0d0);qp=(0d0,0d0);info=1
    endif
  end subroutine initialize_dg_wpw_metric_projected_occupied

  subroutine solve_dg_wpw_metric_projection(context,comm,apply_s,global_gram,nw,np,nrhs,&
      tolerance,max_iterations,diagonal_w,diagonal_p,bw,bp,cw,cp,relative_residual,rhs_residuals,&
      rhs_residual_history,iterations,diagonal_spread,info)
    class(*),intent(inout)::context
    integer,intent(in)::comm,nw,np,nrhs,max_iterations
    procedure(apply_s_batch)::apply_s
    procedure(global_gram_batch)::global_gram
    real(8),intent(in)::tolerance,diagonal_w(nw),diagonal_p(np)
    complex(8),intent(in)::bw(nw,nrhs),bp(np,nrhs)
    complex(8),intent(out)::cw(nw,nrhs),cp(np,nrhs)
    real(8),intent(out)::relative_residual,rhs_residuals(nrhs),rhs_residual_history(max_iterations,nrhs),&
      diagonal_spread
    integer,intent(out)::iterations,info
    complex(8),allocatable::b(:,:),x(:,:),r(:,:),z(:,:),p(:,:),ap(:,:),gram(:,:)
    real(8),allocatable::rho_old(:),rho_new(:),rr(:),denom(:),rhs_norm2(:),diagonal(:),best(:)
    logical,allocatable::active(:)
    integer,allocatable::stagnant(:)
    real(8)::bnorm,local_min,local_max,global_min,global_max
    integer::i,iter,astat,local_bad,global_bad,ierr,apply_info,gram_info

    cw=(0d0,0d0);cp=(0d0,0d0);relative_residual=huge(1d0);rhs_residuals=huge(1d0)
    rhs_residual_history=huge(1d0);diagonal_spread=huge(1d0);iterations=0;info=1
    local_bad=merge(0,1,nw>=0.and.np>=0.and.nw+np>0.and.nrhs>0.and.&
      max_iterations>0.and.tolerance>0d0.and.finite_complex(bw).and.finite_complex(bp).and.&
      all(ieee_is_finite(diagonal_w)).and.all(ieee_is_finite(diagonal_p)).and.&
      all(diagonal_w>0d0).and.all(diagonal_p>0d0))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(b(nw+np,nrhs),x(nw+np,nrhs),r(nw+np,nrhs),z(nw+np,nrhs),p(nw+np,nrhs),&
      ap(nw+np,nrhs),gram(nrhs,nrhs),rho_old(nrhs),rho_new(nrhs),rr(nrhs),denom(nrhs),&
      rhs_norm2(nrhs),diagonal(nw+np),best(nrhs),active(nrhs),stagnant(nrhs),stat=astat)
    local_bad=merge(0,1,astat==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    diagonal(1:nw)=diagonal_w;diagonal(nw+1:nw+np)=diagonal_p
    local_min=minval(diagonal);local_max=maxval(diagonal)
    call MPI_Allreduce(local_min,global_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_max,global_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    diagonal_spread=global_max/global_min
    b(1:nw,:)=bw;b(nw+1:nw+np,:)=bp;x=(0d0,0d0);r=b
    do i=1,nrhs;z(:,i)=r(:,i)/diagonal;enddo
    p=z
    call global_gram(b,b,nw+np,nrhs,nrhs,gram,gram_info)
    if(gram_info/=0)return
    do i=1,nrhs;rhs_norm2(i)=real(gram(i,i),8);enddo
    call global_gram(r,z,nw+np,nrhs,nrhs,gram,gram_info)
    if(gram_info/=0)return
    do i=1,nrhs;rho_old(i)=real(gram(i,i),8);enddo
    local_bad=merge(0,1,all(ieee_is_finite(rhs_norm2)).and.all(rhs_norm2>tolerance*tolerance).and.&
      all(ieee_is_finite(rho_old)).and.all(rho_old>0d0))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    active=.true.;bnorm=sqrt(sum(rhs_norm2));relative_residual=1d0;rhs_residuals=1d0
    best=1d0;stagnant=0
    do iter=1,max_iterations
      call apply_s(context,p(1:nw,:),p(nw+1:nw+np,:),ap(1:nw,:),ap(nw+1:nw+np,:),apply_info)
      if(apply_info/=0)return
      call global_gram(p,ap,nw+np,nrhs,nrhs,gram,gram_info)
      if(gram_info/=0)return
      do i=1,nrhs;denom(i)=real(gram(i,i),8);enddo
      local_bad=0
      do i=1,nrhs
        if(active(i).and.(.not.ieee_is_finite(denom(i)).or.denom(i)<=0d0))local_bad=1
      enddo
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
      do i=1,nrhs
        if(active(i))then
          x(:,i)=x(:,i)+(rho_old(i)/denom(i))*p(:,i)
          r(:,i)=r(:,i)-(rho_old(i)/denom(i))*ap(:,i)
        endif
      enddo
      call global_gram(r,r,nw+np,nrhs,nrhs,gram,gram_info)
      if(gram_info/=0)return
      do i=1,nrhs
        rr(i)=max(0d0,real(gram(i,i),8));rhs_residuals(i)=sqrt(rr(i)/rhs_norm2(i))
      enddo
      rhs_residual_history(iter,:)=rhs_residuals
      relative_residual=sqrt(sum(rr))/bnorm;iterations=iter
      if(.not.ieee_is_finite(relative_residual))return
      if(all(rhs_residuals<=tolerance))then
        cw=x(1:nw,:);cp=x(nw+1:nw+np,:);info=0;return
      endif
      do i=1,nrhs
        if(rhs_residuals(i)<=tolerance)then
          active(i)=.false.;p(:,i)=(0d0,0d0)
        else if(active(i))then
          if(rhs_residuals(i)<best(i)*(1d0-1d-12))then;best(i)=rhs_residuals(i);stagnant(i)=0
          else;stagnant(i)=stagnant(i)+1;endif
        endif
      enddo
      local_bad=merge(1,0,any(stagnant>=12))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
      do i=1,nrhs;z(:,i)=r(:,i)/diagonal;enddo
      call global_gram(r,z,nw+np,nrhs,nrhs,gram,gram_info)
      if(gram_info/=0)return
      do i=1,nrhs;rho_new(i)=real(gram(i,i),8);enddo
      local_bad=0
      do i=1,nrhs
        if(active(i).and.(.not.ieee_is_finite(rho_new(i)).or.rho_new(i)<=0d0.or.rho_old(i)<=0d0))local_bad=1
      enddo
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
      do i=1,nrhs
        if(active(i))p(:,i)=z(:,i)+(rho_new(i)/rho_old(i))*p(:,i)
      enddo
      rho_old=rho_new
    enddo
  end subroutine solve_dg_wpw_metric_projection

  subroutine complete_dg_wpw_projected_subspace(context,comm,apply_s,global_gram,nw,np,nocc,nretain,&
      metric_cutoff,qw,qp,effective_rank,orthogonality,info)
    class(*),intent(inout)::context
    integer,intent(in)::comm,nw,np,nocc,nretain
    procedure(apply_s_batch)::apply_s
    procedure(global_gram_batch)::global_gram
    real(8),intent(in)::metric_cutoff
    complex(8),intent(inout)::qw(nw,nretain),qp(np,nretain)
    integer,intent(out)::effective_rank
    real(8),intent(out)::orthogonality
    integer,intent(out)::info
    complex(8),allocatable::qocc(:,:),qextra(:,:),sextra(:,:),overlap(:,:),qall(:,:),sall(:,:),full_gram(:,:)
    integer::nextra,i

    info=1;effective_rank=0;orthogonality=huge(1d0);nextra=nretain-nocc
    if(nocc<1.or.nextra<1)return
    allocate(qocc(nw+np,nocc),qextra(nw+np,nextra),sextra(nw+np,nextra),overlap(nocc,nextra))
    qocc(1:nw,:)=qw(:,1:nocc);qocc(nw+1:,:)=qp(:,1:nocc)
    qextra(1:nw,:)=qw(:,nocc+1:nretain);qextra(nw+1:,:)=qp(:,nocc+1:nretain)
    call apply_s(context,qextra(1:nw,:),qextra(nw+1:,:),sextra(1:nw,:),sextra(nw+1:,:),info)
    if(info/=0)return
    call global_gram(qocc,sextra,nw+np,nocc,nextra,overlap,info);if(info/=0)return
    qextra=qextra-matmul(qocc,overlap)
    qw(:,nocc+1:nretain)=qextra(1:nw,:);qp(:,nocc+1:nretain)=qextra(nw+1:,:)
    call initialize_dg_wpw_projected_occupied(context,comm,apply_s,global_gram,nw,np,nextra,&
      metric_cutoff,qw(:,nocc+1:nretain),qp(:,nocc+1:nretain),effective_rank,orthogonality,info)
    if(info/=0.or.effective_rank/=nextra)return
    allocate(qall(nw+np,nretain),sall(nw+np,nretain),full_gram(nretain,nretain))
    qall(1:nw,:)=qw;qall(nw+1:,:)=qp
    call apply_s(context,qall(1:nw,:),qall(nw+1:,:),sall(1:nw,:),sall(nw+1:,:),info)
    if(info/=0)return
    call global_gram(qall,sall,nw+np,nretain,nretain,full_gram,info);if(info/=0)return
    do i=1,nretain;full_gram(i,i)=full_gram(i,i)-(1d0,0d0);enddo
    orthogonality=maxval(abs(full_gram))
    if(.not.ieee_is_finite(orthogonality))then;info=1;effective_rank=0;return;endif
    effective_rank=nretain;info=0
  end subroutine complete_dg_wpw_projected_subspace

  subroutine initialize_dg_wpw_projected_occupied(context,comm,apply_s,global_gram,nw,np,nocc,&
      metric_cutoff,qw,qp,effective_rank,orthogonality,info)
    class(*),intent(inout)::context
    integer,intent(in)::comm,nw,np,nocc
    procedure(apply_s_batch)::apply_s
    procedure(global_gram_batch)::global_gram
    real(8),intent(in)::metric_cutoff
    complex(8),intent(inout)::qw(nw,nocc),qp(np,nocc)
    integer,intent(out)::effective_rank
    real(8),intent(out)::orthogonality
    integer,intent(out)::info
    complex(8),allocatable::q(:,:),sq(:,:),metric(:,:),inverse_sqrt(:,:),check(:,:)
    real(8)::emin,emax,condition
    integer::discarded,local_bad,global_bad,ierr,i

    info=1;effective_rank=0;orthogonality=huge(1d0)
    local_bad=merge(0,1,nw>=0.and.np>=0.and.nocc>=1.and.metric_cutoff>0d0.and.&
      finite_complex(qw).and.finite_complex(qp))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(q(nw+np,nocc),sq(nw+np,nocc),metric(nocc,nocc),inverse_sqrt(nocc,nocc),check(nocc,nocc))
    q(1:nw,:)=qw;q(nw+1:nw+np,:)=qp
    call apply_s(context,q(1:nw,:),q(nw+1:nw+np,:),sq(1:nw,:),sq(nw+1:nw+np,:),info)
    if(info/=0)return
    call global_gram(q,sq,nw+np,nocc,nocc,metric,info);if(info/=0)return
    call dg_metric_factor(metric,nocc,metric_cutoff,inverse_sqrt,emin,emax,condition,discarded,info)
    effective_rank=nocc-discarded
    if(info/=0)return
    q=matmul(q,inverse_sqrt)
    call apply_s(context,q(1:nw,:),q(nw+1:nw+np,:),sq(1:nw,:),sq(nw+1:nw+np,:),info)
    if(info/=0)return
    call global_gram(q,sq,nw+np,nocc,nocc,check,info);if(info/=0)return
    do i=1,nocc;check(i,i)=check(i,i)-1d0;enddo
    orthogonality=maxval(abs(check))
    if(.not.ieee_is_finite(orthogonality))then;info=2;return;endif
    qw=q(1:nw,:);qp=q(nw+1:nw+np,:);info=0
  end subroutine initialize_dg_wpw_projected_occupied

  subroutine initialize_dg_wpw_deterministic_seed(w_ids,p_ids,qw,qp,info)
    integer,intent(in)::w_ids(:),p_ids(:)
    complex(8),intent(out)::qw(:,:),qp(:,:)
    integer,intent(out)::info
    integer::i,j,nretain
    info=1;qw=(0d0,0d0);qp=(0d0,0d0)
    nretain=size(qw,2)
    if(size(qw,1)/=size(w_ids).or.size(qp,1)/=size(p_ids).or.size(qp,2)/=nretain.or.&
       nretain<1.or.any(w_ids<=0).or.any(p_ids<=0))return
    do j=1,nretain
      do i=1,size(w_ids)
        qw(i,j)=deterministic_seed_value(1,w_ids(i),j)
      enddo
      do i=1,size(p_ids)
        qp(i,j)=deterministic_seed_value(2,p_ids(i),j)
      enddo
    enddo
    info=0
  end subroutine

  complex(8) function deterministic_seed_value(namespace,stable_id,state_index)result(value)
    integer,intent(in)::namespace,stable_id,state_index
    integer(8),parameter::prime=2147483647_8
    integer(8)::basis,state,real_key,imag_key,real_residue,imag_residue
    basis=int(stable_id,8);state=int(state_index,8)
    real_key=modulo(1009_8*basis+9176_8*state+32452843_8*int(namespace,8),prime)
    imag_key=modulo(2003_8*basis+6113_8*state+49979687_8*int(namespace,8),prime)
    real_residue=mix_seed_key(real_key,104729_8)
    imag_residue=mix_seed_key(imag_key,224737_8)
    value=cmplx((dble(real_residue)+0.5d0)/dble(prime)-0.5d0,&
      (dble(imag_residue)+0.5d0)/dble(prime)-0.5d0,8)
  contains
    integer(8) function mix_seed_key(key,salt)result(mixed)
      integer(8),intent(in)::key,salt
      mixed=modulo(key*key+48271_8*key+salt,prime)
      mixed=modulo(mixed*mixed+69621_8*mixed+3_8*salt,prime)
    end function
  end function

  subroutine run_dg_wpw_matrix_free_algebra_step(context,comm,apply_h,apply_s,global_gram,nw,np,nocc,nretain,&
      iteration,metric_cutoff,residual_tolerance,qw,qp,q_old_occ,eigenvalues,gap,residual,orth,&
      projector_residual,info,precondition)
    class(*),intent(inout)::context
    integer,intent(in)::comm
    procedure(apply_h_batch)::apply_h;procedure(apply_s_batch)::apply_s
    procedure(global_gram_batch)::global_gram
    procedure(dg_wpw_preconditioner),optional::precondition
    integer,intent(in)::nw,np,nocc,nretain,iteration
    real(8),intent(in)::metric_cutoff,residual_tolerance
    complex(8),intent(inout)::qw(nw,nretain),qp(np,nretain)
    complex(8),intent(in)::q_old_occ(nw+np,nocc)
    real(8),intent(out)::eigenvalues(nretain),gap,residual,orth,projector_residual
    integer,intent(out)::info
    complex(8),allocatable::q(:,:)
    integer::local_bad,global_bad,ierr
    info=1;gap=0d0;residual=huge(1d0);orth=huge(1d0);projector_residual=huge(1d0)
    local_bad=merge(0,1,nw>=0.and.np>=0.and.nw+np>=1.and.nocc>=1.and.nretain>nocc.and.iteration>=1.and.&
       metric_cutoff>0d0.and.residual_tolerance>0d0.and.finite_complex(qw).and.&
       finite_complex(qp).and.finite_complex(q_old_occ))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(q(nw+np,nretain));q(1:nw,:)=qw;q(nw+1:nw+np,:)=qp
    call solve_window(context,comm,apply_h,apply_s,global_gram,nw,np,nretain,metric_cutoff,&
      residual_tolerance,q,eigenvalues,residual,orth,info,precondition)
    if(info/=0)return
    gap=eigenvalues(nocc+1)-eigenvalues(nocc)
    if(.not.ieee_is_finite(gap))then;info=2;return;endif
    projector_residual=projector_change(context,apply_s,q_old_occ,q(:,1:nocc),nw,np,iteration,&
      global_gram,info)
    if(info/=0)return
    qw=q(1:nw,:);qp=q(nw+1:nw+np,:);info=0
  end subroutine

  subroutine run_dg_wpw_matrix_free_scf(context,comm,apply_h,apply_s,global_gram,potential_map,&
      n_w_local,n_p_local,n_occ,extra_states,occ,density,mix,max_iter,metric_cutoff,&
      gap_threshold,residual_tolerance,qw,qp,eigenvalues,result,info,precondition)
    class(*),intent(inout)::context
    integer,intent(in)::comm
    procedure(apply_h_batch)::apply_h
    procedure(apply_s_batch)::apply_s
    procedure(global_gram_batch)::global_gram
    procedure(dg_wpw_preconditioner),optional::precondition
    procedure(dg_wpw_scf_map)::potential_map
    integer,intent(in)::n_w_local,n_p_local,n_occ,extra_states,max_iter
    real(8),intent(in)::occ(n_occ),mix,metric_cutoff,gap_threshold,residual_tolerance
    real(8),intent(inout)::density(:)
    complex(8),intent(inout)::qw(:,:),qp(:,:)
    real(8),intent(out)::eigenvalues(:)
    type(s_dg_wpw_matrix_free_scf_result),intent(out)::result
    integer,intent(out)::info
    integer::nretain,nlocal,iter,map_info,local_bad,global_bad,ierr
    real(8)::energy,energy_old,potential_residual,charge_error
    real(8),allocatable::density_new(:),density_check(:)
    complex(8),allocatable::q(:,:),q_old_occ(:,:)

    result=s_dg_wpw_matrix_free_scf_result();info=0
    nretain=n_occ + extra_states
    result%retained_dimension=nretain
    nlocal=n_w_local+n_p_local
    local_bad=0
    if(n_w_local<0.or.n_p_local<0.or.nlocal<1.or.n_occ<1.or.extra_states<1.or.&
       max_iter<1.or.mix<=0d0.or.mix>1d0.or.metric_cutoff<=0d0.or.gap_threshold<=0d0.or.&
       residual_tolerance<=0d0.or.size(density)<1.or.any(occ<=0d0).or.&
       any(shape(qw)/=[n_w_local,nretain]).or.any(shape(qp)/=[n_p_local,nretain]).or.&
       size(eigenvalues)/=nretain)then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=1;return;endif
    allocate(q(nlocal,nretain),q_old_occ(nlocal,n_occ),density_new(size(density)),density_check(size(density)))
    q(1:n_w_local,:)=qw;q(n_w_local+1:nlocal,:)=qp
    if(.not.finite_complex(q).or..not.all(ieee_is_finite(density)))then;info=2;return;endif
    energy_old=huge(1d0);q_old_occ=(0d0,0d0)

    do iter=1,max_iter
      call solve_window(context,comm,apply_h,apply_s,global_gram,n_w_local,n_p_local,nretain,&
        metric_cutoff,residual_tolerance,q,eigenvalues,result%generalized_residual,&
        result%metric_orthonormality,info,precondition)
      if(info/=0)return
      result%gap=eigenvalues(n_occ+1)-eigenvalues(n_occ)
      if(.not.ieee_is_finite(result%gap).or.result%gap<gap_threshold)then;info=3;return;endif
      result%projector_residual=projector_change(context,apply_s,q_old_occ,q(:,1:n_occ),&
        n_w_local,n_p_local,iter,global_gram,info)
      if(info/=0)return
      call potential_map(context,q(1:n_w_local,1:n_occ),q(n_w_local+1:nlocal,1:n_occ),&
        occ,n_occ,density,mix,density_new,potential_residual,energy,charge_error,map_info)
      if(map_info/=0.or..not.all(ieee_is_finite(density_new)).or.&
         .not.all(ieee_is_finite([potential_residual,energy,charge_error])))then
        info=20+max(1,map_info);return
      endif
      result%density_residual=relative_real(density_new,density)
      result%potential_residual=potential_residual
      if(iter>1)result%energy_residual=abs(energy-energy_old)
      result%charge_error=charge_error;result%total_energy=energy;result%iterations=iter
      if(iter>1.and.max(result%density_residual,result%potential_residual)<residual_tolerance.and.&
         result%energy_residual<residual_tolerance.and.result%projector_residual<residual_tolerance.and.&
         result%generalized_residual<residual_tolerance.and.result%metric_orthonormality<residual_tolerance.and.&
         abs(result%charge_error)<residual_tolerance)then
        result%converged=.true.;exit
      endif
      density=density_new;energy_old=energy;q_old_occ=q(:,1:n_occ)
    enddo
    if(.not.result%converged)then;info=30;return;endif

    call potential_map(context,q(1:n_w_local,1:n_occ),q(n_w_local+1:nlocal,1:n_occ),&
      occ,n_occ,density_new,1d0,density_check,potential_residual,energy,charge_error,map_info)
    if(map_info/=0.or..not.all(ieee_is_finite(density_check)).or.&
       .not.all(ieee_is_finite([potential_residual,energy,charge_error])))then;info=31;return;endif
    result%fixed_point_residual=relative_real(density_check,density_new)
    if(result%fixed_point_residual>=residual_tolerance.or.abs(potential_residual)>=residual_tolerance.or.&
       abs(energy-result%total_energy)>=residual_tolerance.or.abs(charge_error)>=residual_tolerance)then
      result%converged=.false.;info=32;return
    endif
    density=density_check;qw=q(1:n_w_local,:);qp=q(n_w_local+1:nlocal,:)
  end subroutine

  subroutine solve_window(context,comm,apply_h,apply_s,gram,nw,np,nretain,cutoff,tol,q,eval,residual,orth,info,&
      precondition)
    class(*),intent(inout)::context
    procedure(apply_h_batch)::apply_h;procedure(apply_s_batch)::apply_s;procedure(global_gram_batch)::gram
    procedure(dg_wpw_preconditioner),optional::precondition
    integer,intent(in)::comm,nw,np,nretain;real(8),intent(in)::cutoff,tol
    complex(8),intent(inout)::q(nw+np,nretain);real(8),intent(out)::eval(nretain),residual,orth
    integer,intent(out)::info
    complex(8)::hq(nw+np,nretain),sq(nw+np,nretain),a(nretain,nretain),b(nretain,nretain)
    complex(8)::u(nretain,nretain),r(nw+np,nretain),g(nretain,nretain)
    complex(8)::search(nw+np,nretain),preconditioned(nw+np,nretain),z(nw+np,3*nretain)
    complex(8)::hz(nw+np,3*nretain),sz(nw+np,3*nretain)
    complex(8)::reduced_h(3*nretain,3*nretain),reduced_s(3*nretain,3*nretain)
    complex(8)::reduced_c(3*nretain,nretain)
    real(8)::condition,reduced_eval(nretain),reduced_residual,reduced_orth
    integer::discarded,effective_rank,i,inner
    real(8)::stage_clock
    info=0;residual=huge(1d0);orth=huge(1d0);search=(0d0,0d0);call cpu_time(stage_clock)
    do inner=1,160
      call trace_solve_window('initial_apply_begin',comm,stage_clock,inner,0,residual,orth)
      call apply_s(context,q(1:nw,:),q(nw+1:nw+np,:),sq(1:nw,:),sq(nw+1:nw+np,:),info);if(info/=0)return
      call apply_h(context,q(1:nw,:),q(nw+1:nw+np,:),hq(1:nw,:),hq(nw+1:nw+np,:),info);if(info/=0)return
      call trace_solve_window('initial_apply_end',comm,stage_clock,inner,0,residual,orth)
      call gram(q,hq,nw+np,nretain,nretain,a,info);if(info/=0)return
      call gram(q,sq,nw+np,nretain,nretain,b,info);if(info/=0)return
      call dg_generalized_eigh(a,b,nretain,cutoff,eval,u,residual,orth,condition,discarded,info);if(info/=0)return
      call trace_solve_window('initial_eigh_end',comm,stage_clock,inner,nretain-discarded,residual,orth)
      q=matmul(q,u);hq=matmul(hq,u);sq=matmul(sq,u)
      do i=1,nretain;r(:,i)=hq(:,i)-eval(i)*sq(:,i);enddo
      call gram(r,r,nw+np,nretain,nretain,g,info);if(info/=0)return
      residual=sqrt(max(0d0,maxval(real([(g(i,i),i=1,nretain)],8))))
      call gram(q,sq,nw+np,nretain,nretain,g,info);if(info/=0)return
      do i=1,nretain;g(i,i)=g(i,i)-1d0;enddo
      orth=maxval(abs(g))
      if(residual<tol.and.orth<tol)return
      if(present(precondition))then
        call precondition(context,eval,r(1:nw,:),r(nw+1:nw+np,:),preconditioned(1:nw,:),&
          preconditioned(nw+1:nw+np,:),info)
        if(info/=0)return
      else
        preconditioned=r
      endif
      z(:,1:nretain)=q;z(:,nretain+1:2*nretain)=preconditioned;z(:,2*nretain+1:)=search
      call apply_s(context,z(1:nw,:),z(nw+1:nw+np,:),sz(1:nw,:),sz(nw+1:nw+np,:),info)
      if(info/=0)return
      call apply_h(context,z(1:nw,:),z(nw+1:nw+np,:),hz(1:nw,:),hz(nw+1:nw+np,:),info)
      if(info/=0)return
      call trace_solve_window('expanded_apply_end',comm,stage_clock,inner,0,residual,orth)
      call gram(z,hz,nw+np,3*nretain,3*nretain,reduced_h,info);if(info/=0)return
      call gram(z,sz,nw+np,3*nretain,3*nretain,reduced_s,info);if(info/=0)return
      call dg_reduced_generalized_eigh(reduced_h,reduced_s,3*nretain,nretain,cutoff,&
        reduced_eval,reduced_c,effective_rank,reduced_residual,reduced_orth,info)
      if(info/=0)return
      call trace_solve_window('reduced_eigh_end',comm,stage_clock,inner,effective_rank,&
        reduced_residual,reduced_orth)
      q=matmul(z,reduced_c);eval=reduced_eval
      search=matmul(preconditioned,reduced_c(nretain+1:2*nretain,:))+&
        matmul(search,reduced_c(2*nretain+1:3*nretain,:))
      call trace_solve_window('inner_end',comm,stage_clock,inner,effective_rank,residual,orth)
    enddo
    info=40
  end subroutine

  subroutine trace_solve_window(stage,comm,stage_clock,inner,effective_rank,residual,orth)
    character(*),intent(in)::stage
    integer,intent(in)::comm,inner,effective_rank
    real(8),intent(inout)::stage_clock
    real(8),intent(in)::residual,orth
    integer::rank,ierr
    real(8)::now
    call cpu_time(now);call MPI_Comm_rank(comm,rank,ierr)
    if(ierr==MPI_SUCCESS.and.rank==0)then
      write(*,'(1x,a,a,a,i0,a,i0,3(a,es12.4))')'[DG-WPW-SOLVE-STAGE] stage=',trim(stage),&
        ' inner=',inner,' rank=',effective_rank,' residual=',residual,' orth=',orth,&
        ' cpu_seconds=',now-stage_clock
      flush(6)
    endif
    stage_clock=now
  end subroutine trace_solve_window

  real(8) function projector_change(context,apply_s,old,new,nw,np,iter,gram,info) result(value)
    class(*),intent(inout)::context
    procedure(apply_s_batch)::apply_s
    complex(8),intent(in)::old(:,:),new(:,:);integer,intent(in)::nw,np,iter
    procedure(global_gram_batch)::gram;integer,intent(out)::info
    complex(8)::overlap(size(old,2),size(new,2))
    complex(8)::snew(size(new,1),size(new,2))
    real(8)::deficit
    if(iter==1)then;value=huge(1d0);info=0;return;endif
    call apply_s(context,new(1:nw,:),new(nw+1:nw+np,:),snew(1:nw,:),snew(nw+1:nw+np,:),info)
    if(info/=0)then;value=huge(1d0);return;endif
    call gram(old,snew,size(old,1),size(old,2),size(new,2),overlap,info)
    if(info/=0)then;value=huge(1d0);return;endif
    deficit=1d0-sum(abs(overlap)**2)/dble(size(new,2))
    if(abs(deficit)<=100d0*epsilon(1d0))deficit=0d0
    value=sqrt(max(0d0,deficit))
  end function
  real(8) function relative_real(a,b) result(value)
    real(8),intent(in)::a(:),b(:);value=sqrt(sum((a-b)**2))/max(1d-30,sqrt(sum(a**2)))
  end function
  logical function finite_complex(a) result(ok)
    complex(8),intent(in)::a(:,:);ok=all(ieee_is_finite(real(a,8))).and.all(ieee_is_finite(aimag(a)))
  end function
end module
