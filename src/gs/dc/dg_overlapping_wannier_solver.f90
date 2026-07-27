#include "config.h"
module dg_overlapping_wannier_solver
  use,intrinsic::iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::solve_dg_overlapping_wannier_coefficients
contains
  subroutine solve_dg_overlapping_wannier_coefficients(comm,row_ids,hrows,srows,nstate,max_iterations,&
      tolerance,metric_tolerance,coefficients,eigenvalues,maximum_residual,orthogonality_defect,&
      metric_condition,ok,message)
    integer,intent(in)::comm,nstate,max_iterations
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::hrows(:,:),srows(:,:)
    real(8),intent(in)::tolerance,metric_tolerance
    complex(8),intent(out)::coefficients(:,:)
    real(8),intent(out)::eigenvalues(:),maximum_residual,orthogonality_defect,metric_condition
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(8),allocatable::qlocal(:,:),vector(:),hvector(:),svector(:),residual(:),direction(:),&
      hdirection(:),sdirection(:),gram(:,:)
    complex(8)::value,reduced_h(2,2),reduced_s(2,2),ritz_vector(2)
    real(8)::ritz_value,pivot,residual_norm
    integer,allocatable::ownership(:)
    integer::n,nlocal,i,j,band,iteration,ierr,local_bad,global_bad
    integer::scalar_local(3),scalar_min(3),scalar_max(3)
    real(8)::real_local(2),real_min(2),real_max(2)

    ok=.false.;message='';coefficients=(0d0,0d0);eigenvalues=0d0
    maximum_residual=huge(1d0);orthogonality_defect=huge(1d0);metric_condition=huge(1d0)
    n=size(hrows,2);nlocal=size(row_ids)
    scalar_local=[n,nstate,max_iterations]
    call MPI_Allreduce(scalar_local,scalar_min,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(scalar_local,scalar_max,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    real_local=[tolerance,metric_tolerance]
    call MPI_Allreduce(real_local,real_min,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(real_local,real_max,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    local_bad=0
    if(any(scalar_min/=scalar_max).or.any(real_min/=real_max))local_bad=1
    if(n<1.or.nstate<1.or.nstate>n.or.max_iterations<1)local_bad=1
    if(tolerance<=0d0.or.metric_tolerance<=0d0.or.metric_tolerance>=1d0)local_bad=1
    if(.not.all(ieee_is_finite(real_local)))local_bad=1
    if(size(hrows,1)/=nlocal.or.any(shape(srows)/=shape(hrows)))local_bad=1
    if(any(shape(coefficients)/=[n,nstate]).or.size(eigenvalues)/=nstate)local_bad=1
    if(any(row_ids<1_8).or.any(row_ids>int(n,8)))local_bad=1
    if(.not.finite_matrix(hrows).or..not.finite_matrix(srows))local_bad=1
    if(int(n,int64)>0_int64.and.int(n,int64)>int(huge(1),int64)/int(n,int64))local_bad=1
    if(int(nstate,int64)>0_int64.and.int(nstate,int64)>int(huge(1),int64)/int(nstate,int64))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid distributed overlapping-Wannier coefficient solve contract';return;endif

    allocate(ownership(n));ownership=0
    do i=1,nlocal;ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(any(ownership/=1))then;message='coefficient rows do not form a unique global partition';return;endif
    call distributed_hermiticity(comm,row_ids,hrows,tolerance,local_bad)
    call distributed_hermiticity(comm,row_ids,srows,metric_tolerance,global_bad)
    if(local_bad/=0.or.global_bad/=0)then;message='non-Hermitian distributed H or S';return;endif

    allocate(qlocal(nlocal,n),vector(n),hvector(nlocal),svector(nlocal),residual(nlocal),&
      direction(nlocal),hdirection(nlocal),sdirection(nlocal),gram(nstate,nstate))
    qlocal=(0d0,0d0)
    do j=1,n
      do i=1,nlocal
        if(row_ids(i)==int(j,8))qlocal(i,j)=1d0
      enddo
      call s_orthogonalize(comm,row_ids,srows,qlocal(:,j),qlocal(:,1:j-1),pivot,ok)
      if(.not.ok)then;message='overlap metric is indefinite or numerically singular';return;endif
    enddo
    call estimate_metric_condition(comm,row_ids,srows,qlocal,metric_condition,ok)
    if(.not.ok.or.metric_condition*metric_tolerance>=1d0)then
      ok=.false.;message='overlap metric condition gate failed';return
    endif

    do band=1,nstate
      call s_orthogonalize(comm,row_ids,srows,qlocal(:,band),qlocal(:,1:band-1),pivot,ok)
      if(.not.ok)then;message='starting coefficient vector lost S rank';return;endif
      do iteration=1,max_iterations
        call gather_vector(comm,row_ids,qlocal(:,band),vector)
        hvector=matmul(hrows,vector);svector=matmul(srows,vector)
        call inner(comm,qlocal(:,band),hvector,value)
        eigenvalues(band)=real(value,8)
        residual=hvector-eigenvalues(band)*svector
        call relative_residual(comm,residual,hvector,svector,eigenvalues(band),residual_norm)
        if(.not.ieee_is_finite(residual_norm))then;ok=.false.;message='nonfinite coefficient residual';return;endif
        if(residual_norm<=tolerance)exit
        direction=-residual
        call s_project_out(comm,row_ids,srows,direction,qlocal(:,1:band-1),ok)
        if(.not.ok)return
        call s_orthogonalize(comm,row_ids,srows,direction,qlocal(:,1:band-1),pivot,ok)
        if(.not.ok)then;message='coefficient search direction lost metric rank';return;endif
        call gather_vector(comm,row_ids,direction,vector)
        hdirection=matmul(hrows,vector);sdirection=matmul(srows,vector)
        if(.not.finite_matrix(reshape(hdirection,[size(hdirection),1])).or.&
            .not.finite_matrix(reshape(sdirection,[size(sdirection),1])))then
          ok=.false.;message='nonfinite coefficient search action';return
        endif
        call inner(comm,qlocal(:,band),hvector,reduced_h(1,1))
        call inner(comm,qlocal(:,band),hdirection,reduced_h(1,2))
        call inner(comm,direction,hvector,reduced_h(2,1))
        call inner(comm,direction,hdirection,reduced_h(2,2))
        call inner(comm,qlocal(:,band),svector,reduced_s(1,1))
        call inner(comm,qlocal(:,band),sdirection,reduced_s(1,2))
        call inner(comm,direction,svector,reduced_s(2,1))
        call inner(comm,direction,sdirection,reduced_s(2,2))
        reduced_h=0.5d0*(reduced_h+conjg(transpose(reduced_h)))
        reduced_s=0.5d0*(reduced_s+conjg(transpose(reduced_s)))
        if(.not.finite_matrix(reduced_h).or..not.finite_matrix(reduced_s))then
          ok=.false.;message='nonfinite reduced coefficient problem';return
        endif
        if(real(reduced_s(2,2),8)<=metric_tolerance.or.real(reduced_s(1,1),8)<=metric_tolerance.or.&
            real(reduced_s(2,2),8)-abs(reduced_s(1,2))**2/real(reduced_s(1,1),8)<=metric_tolerance)then
          ok=.false.;message='coefficient search space lost positive metric rank';return
        endif
        call lowest_generalized_2x2(reduced_h,reduced_s,ritz_value,ritz_vector,ok)
        if(.not.ok)then;message='coefficient Rayleigh-Ritz step failed';return;endif
        qlocal(:,band)=qlocal(:,band)*ritz_vector(1)+direction*ritz_vector(2)
        call s_orthogonalize(comm,row_ids,srows,qlocal(:,band),qlocal(:,1:band-1),pivot,ok)
        if(.not.ok)return
      enddo
      if(iteration>max_iterations)then;ok=.false.;message='coefficient iteration did not converge';return;endif
    enddo
    do band=1,nstate
      call gather_vector(comm,row_ids,qlocal(:,band),coefficients(:,band))
    enddo
    call block_rayleigh_ritz(comm,row_ids,hrows,srows,qlocal(:,1:nstate),coefficients,eigenvalues,ok)
    if(.not.ok)then;message='final block Rayleigh-Ritz failed';return;endif
    call coefficient_diagnostics(comm,row_ids,hrows,srows,coefficients,eigenvalues,&
      maximum_residual,orthogonality_defect)
    if(maximum_residual>10d0*tolerance.or.orthogonality_defect>10d0*tolerance)then
      ok=.false.;message='coefficient residual or S-orthonormality gate failed';return
    endif
    ok=.true.;message=''
#else
    ok=.false.;message='overlapping-Wannier coefficient solve requires MPI'
#endif
  end subroutine

#ifdef USE_MPI
  subroutine gather_vector(comm,row_ids,local,global)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::local(:)
    complex(8),intent(out)::global(:)
    complex(8),allocatable::staged(:)
    integer::i,ierr
    allocate(staged(size(global)));staged=(0d0,0d0)
    do i=1,size(local);staged(int(row_ids(i)))=local(i);enddo
    call MPI_Allreduce(staged,global,size(global),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  end subroutine

  subroutine inner(comm,left,right,value)
    integer,intent(in)::comm
    complex(8),intent(in)::left(:),right(:)
    complex(8),intent(out)::value
    complex(8)::local
    integer::ierr
    local=sum(conjg(left)*right)
    call MPI_Allreduce(local,value,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  end subroutine

  subroutine s_orthogonalize(comm,row_ids,srows,vector,previous,pivot,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::srows(:,:),previous(:,:)
    complex(8),intent(inout)::vector(:)
    real(8),intent(out)::pivot
    logical,intent(out)::ok
    complex(8),allocatable::global(:),svector(:)
    complex(8)::value
    integer::j
    allocate(global(size(srows,2)),svector(size(vector)));ok=.true.
    do j=1,size(previous,2)
      call gather_vector(comm,row_ids,vector,global);svector=matmul(srows,global)
      call inner(comm,previous(:,j),svector,value);vector=vector-previous(:,j)*value
    enddo
    call gather_vector(comm,row_ids,vector,global);svector=matmul(srows,global)
    call inner(comm,vector,svector,value);pivot=real(value,8)
    ok=ieee_is_finite(pivot).and.pivot>0d0.and.abs(aimag(value))<=1d-10*max(1d0,abs(value))
    if(ok)vector=vector/sqrt(pivot)
  end subroutine

  subroutine s_project_out(comm,row_ids,srows,vector,previous,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::srows(:,:),previous(:,:)
    complex(8),intent(inout)::vector(:)
    logical,intent(out)::ok
    complex(8),allocatable::global(:),svector(:)
    complex(8)::value
    integer::j
    allocate(global(size(srows,2)),svector(size(vector)));ok=.true.
    do j=1,size(previous,2)
      call gather_vector(comm,row_ids,vector,global);svector=matmul(srows,global)
      call inner(comm,previous(:,j),svector,value);vector=vector-previous(:,j)*value
    enddo
    ok=finite_matrix(reshape(vector,[size(vector),1]))
  end subroutine

  subroutine distributed_hermiticity(comm,row_ids,rows,tolerance,bad)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::rows(:,:)
    real(8),intent(in)::tolerance
    integer,intent(out)::bad
    complex(8),allocatable::full(:,:)
    real(8)::scale
    integer::i,ierr,n
    n=size(rows,2);allocate(full(n,n));full=(0d0,0d0)
    do i=1,size(row_ids);full(int(row_ids(i)),:)=rows(i,:);enddo
    call MPI_Allreduce(MPI_IN_PLACE,full,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    scale=max(1d0,maxval(abs(full)))
    bad=merge(1,0,maxval(abs(full-conjg(transpose(full))))>tolerance*scale)
  end subroutine

  subroutine relative_residual(comm,residual,hvector,svector,eigenvalue,value)
    integer,intent(in)::comm
    complex(8),intent(in)::residual(:),hvector(:),svector(:)
    real(8),intent(in)::eigenvalue
    real(8),intent(out)::value
    complex(8)::z
    real(8)::nr,nh,ns
    call inner(comm,residual,residual,z);nr=sqrt(max(0d0,real(z,8)))
    call inner(comm,hvector,hvector,z);nh=sqrt(max(0d0,real(z,8)))
    call inner(comm,svector,svector,z);ns=sqrt(max(0d0,real(z,8)))
    value=nr/max(tiny(1d0),nh+abs(eigenvalue)*ns)
  end subroutine

  subroutine estimate_metric_condition(comm,row_ids,srows,q,condition,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::srows(:,:),q(:,:)
    real(8),intent(out)::condition
    logical,intent(out)::ok
    complex(8),allocatable::local_gram(:,:),gram(:,:)
    real(8)::squares_s,squares_inverse,local_squares,local_scale,scale_s,scale_inverse,&
      norm_s,norm_inverse
    real(8)::factor_s,factor_inverse
    integer::ierr,n
    n=size(q,2);allocate(local_gram(n,n),gram(n,n))
    local_scale=maxval(abs(srows))
    call MPI_Allreduce(local_scale,scale_s,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(scale_s<=0d0.or..not.ieee_is_finite(scale_s))then;ok=.false.;return;endif
    local_squares=sum((abs(srows)/scale_s)**2)
    call MPI_Allreduce(local_squares,squares_s,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    local_gram=matmul(conjg(transpose(q)),q)
    call MPI_Allreduce(local_gram,gram,n*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    scale_inverse=maxval(abs(gram))
    if(scale_inverse<=0d0.or..not.ieee_is_finite(scale_inverse))then;ok=.false.;return;endif
    squares_inverse=sum((abs(gram)/scale_inverse)**2)
    factor_s=sqrt(squares_s);factor_inverse=sqrt(squares_inverse)
    if(factor_s<=0d0)then;ok=.false.;return;endif
    if(factor_inverse<=0d0)then;ok=.false.;return;endif
    if(factor_s>1d0)then
      if(scale_s>huge(1d0)/factor_s)then
        condition=huge(1d0);ok=.false.;return
      endif
    endif
    if(factor_inverse>1d0)then
      if(scale_inverse>huge(1d0)/factor_inverse)then
        condition=huge(1d0);ok=.false.;return
      endif
    endif
    norm_s=scale_s*factor_s;norm_inverse=scale_inverse*factor_inverse
    if(norm_inverse>1d0)then
      if(norm_s>huge(1d0)/norm_inverse)then
        condition=huge(1d0);ok=.false.;return
      else
        condition=norm_s*norm_inverse
      endif
    else
      condition=norm_s*norm_inverse
    endif
    ok=ieee_is_finite(condition).and.condition>=1d0
  end subroutine

  subroutine block_rayleigh_ritz(comm,row_ids,hrows,srows,qlocal,coefficients,eigenvalues,ok)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::hrows(:,:),srows(:,:)
    complex(8),intent(inout)::qlocal(:,:)
    complex(8),intent(out)::coefficients(:,:)
    real(8),intent(out)::eigenvalues(:)
    logical,intent(out)::ok
    complex(8),allocatable::hq(:,:),projected(:,:),local_projected(:,:),rotation(:,:),&
      rotated_projected(:,:),rotated_q(:,:)
    complex(8)::subh(2,2),subs(2,2),v(2),tmpcol(size(qlocal,1))
    real(8)::ev,offdiag
    integer::p,q,sweep,i,ierr,n
    n=size(eigenvalues)
    hq=matmul(hrows,coefficients)
    local_projected=matmul(conjg(transpose(qlocal)),hq)
    allocate(projected(size(eigenvalues),size(eigenvalues)))
    call MPI_Allreduce(local_projected,projected,size(eigenvalues)**2,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    projected=0.5d0*(projected+conjg(transpose(projected)))
    if(.not.finite_matrix(projected))then;ok=.false.;return;endif
    allocate(rotation(n,n),rotated_projected(n,n),rotated_q(size(qlocal,1),n))
    subs=(0d0,0d0);subs(1,1)=1d0;subs(2,2)=1d0
    do sweep=1,100
      offdiag=0d0
      do q=2,n;do p=1,q-1
        offdiag=max(offdiag,abs(projected(p,q)))
        if(abs(projected(p,q))<=1d-13*max(1d0,maxval(abs(projected))))cycle
        subh=reshape([projected(p,p),projected(q,p),projected(p,q),projected(q,q)],[2,2])
        call lowest_generalized_2x2(subh,subs,ev,v,ok);if(.not.ok)return
        rotation=(0d0,0d0);do i=1,n;rotation(i,i)=1d0;enddo
        rotation(p,p)=v(1);rotation(q,p)=v(2)
        rotation(p,q)=-conjg(v(2));rotation(q,q)=conjg(v(1))
        rotated_projected=matmul(conjg(transpose(rotation)),matmul(projected,rotation))
        rotated_q=matmul(qlocal,rotation)
        projected=0.5d0*(rotated_projected+conjg(transpose(rotated_projected)))
        qlocal=rotated_q
      enddo;enddo
      if(offdiag<=1d-12*max(1d0,maxval(abs(projected))))exit
    enddo
    if(sweep>100.or..not.finite_matrix(projected))then;ok=.false.;return;endif
    eigenvalues=[(real(projected(p,p),8),p=1,size(eigenvalues))]
    do p=1,size(eigenvalues)-1
      q=minloc(eigenvalues(p:),dim=1)+p-1
      if(q/=p)then
        ev=eigenvalues(p);eigenvalues(p)=eigenvalues(q);eigenvalues(q)=ev
        tmpcol=qlocal(:,p);qlocal(:,p)=qlocal(:,q);qlocal(:,q)=tmpcol
      endif
    enddo
    do p=1,size(eigenvalues);call gather_vector(comm,row_ids,qlocal(:,p),coefficients(:,p));enddo
    ok=all(ieee_is_finite(eigenvalues)).and.finite_matrix(coefficients)
  end subroutine

  subroutine coefficient_diagnostics(comm,row_ids,hrows,srows,c,e,residual,orthogonality)
    integer,intent(in)::comm
    integer(8),intent(in)::row_ids(:)
    complex(8),intent(in)::hrows(:,:),srows(:,:),c(:,:)
    real(8),intent(in)::e(:)
    real(8),intent(out)::residual,orthogonality
    complex(8),allocatable::rlocal(:,:),gram(:,:),local_gram(:,:),hv(:),sv(:)
    real(8)::band_residual
    integer::j,ierr
    allocate(rlocal(size(row_ids),size(e)),gram(size(e),size(e)),local_gram(size(e),size(e)))
    allocate(hv(size(row_ids)),sv(size(row_ids)));residual=0d0
    do j=1,size(e)
      hv=matmul(hrows,c(:,j));sv=matmul(srows,c(:,j));rlocal(:,j)=hv-e(j)*sv
      call relative_residual(comm,rlocal(:,j),hv,sv,e(j),band_residual)
      residual=max(residual,band_residual)
    enddo
    local_gram=matmul(conjg(transpose(c(int(row_ids),:))),matmul(srows,c))
    call MPI_Allreduce(local_gram,gram,size(e)**2,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do j=1,size(e);gram(j,j)=gram(j,j)-1d0;enddo
    orthogonality=maxval(abs(gram))
    if(.not.all(ieee_is_finite([residual,orthogonality])))then
      residual=huge(1d0);orthogonality=huge(1d0)
    endif
  end subroutine

  subroutine lowest_generalized_2x2(h,s,eigenvalue,eigenvector,ok)
    complex(8),intent(in)::h(2,2),s(2,2)
    real(8),intent(out)::eigenvalue
    complex(8),intent(out)::eigenvector(2)
    logical,intent(out)::ok
    complex(8)::l(2,2),linv(2,2),a(2,2),u(2),phase,hs(2,2),temporary(2,2)
    real(8)::l11,l22,s22_root,disc_scaled,norm,scale,delta,t,c,absb_scaled,&
      difference_scaled,eigenvalue_scaled,xdisc,ydisc,metric_ratio,phase_threshold
    real(8)::hamiltonian_scale,total_scale
    complex(8)::quotient
    integer::i,j,k
    ok=.false.;eigenvalue=0d0;eigenvector=(0d0,0d0)
    if(.not.finite_matrix(h).or..not.finite_matrix(s))return
    if(real(s(1,1),8)<=0d0.or.real(s(2,2),8)<=0d0)return
    l11=sqrt(real(s(1,1),8))
    l=(0d0,0d0);l(1,1)=l11
    call safe_complex_divide(s(2,1),l11,l(2,1),ok);if(.not.ok)return
    s22_root=sqrt(real(s(2,2),8))
    if(safe_complex_abs(l(2,1))>=s22_root)return
    metric_ratio=safe_complex_abs(l(2,1))/s22_root
    l22=s22_root*sqrt((1d0-metric_ratio)*(1d0+metric_ratio));l(2,2)=l22
    if(.not.all(ieee_is_finite([l11,l22])).or.min(l11,l22)<=0d0)return
    if(l11<1d0/huge(1d0).or.l22<1d0/huge(1d0))return
    linv=(0d0,0d0);linv(1,1)=1d0/l11;linv(2,2)=1d0/l22
    call safe_complex_divide(l(2,1),l11,quotient,ok);if(.not.ok)return
    call safe_complex_divide(-quotient,l22,linv(2,1),ok);if(.not.ok)return
    hamiltonian_scale=max(tiny(1d0),max(maxval(abs(real(h,8))),maxval(abs(aimag(h)))))
    hs=h/hamiltonian_scale
    temporary=(0d0,0d0);a=(0d0,0d0)
    do j=1,2;do i=1,2;do k=1,2
      temporary(i,j)=temporary(i,j)+linv(i,k)*hs(k,j)
    enddo;enddo;enddo
    do j=1,2;do i=1,2;do k=1,2
      a(i,j)=a(i,j)+temporary(i,k)*conjg(linv(j,k))
    enddo;enddo;enddo
    a=0.5d0*(a+conjg(transpose(a)))
    if(.not.finite_matrix(a))return
    scale=max(tiny(1d0),max(maxval(abs(real(a,8))),maxval(abs(aimag(a)))))
    delta=real(a(1,1),8)/scale-real(a(2,2),8)/scale
    xdisc=abs(delta)
    absb_scaled=safe_real_hypot(real(a(1,2),8)/scale,aimag(a(1,2))/scale)
    ydisc=2d0*absb_scaled
    disc_scaled=sqrt(xdisc*xdisc+ydisc*ydisc)
    eigenvalue_scaled=0.5d0*(real(a(1,1),8)/scale+real(a(2,2),8)/scale-disc_scaled)
    if(scale>1d0)then
      if(hamiltonian_scale>huge(1d0)/scale)return
    endif
    total_scale=scale*hamiltonian_scale
    if(abs(eigenvalue_scaled)>1d0)then
      if(total_scale>huge(1d0)/abs(eigenvalue_scaled))return
    endif
    eigenvalue=total_scale*eigenvalue_scaled
    phase_threshold=epsilon(1d0)*max(1d0,1d0/scale)
    if(absb_scaled>phase_threshold)then
      phase=(a(1,2)/scale)/absb_scaled
      if(real(a(1,1),8)<=real(a(2,2),8))then
        difference_scaled=real(a(2,2),8)/scale-real(a(1,1),8)/scale
        t=2d0*absb_scaled/(difference_scaled+disc_scaled);c=1d0/sqrt(1d0+t*t)
        u=[cmplx(c,0d0,8),-t*c*conjg(phase)]
      else
        difference_scaled=real(a(1,1),8)/scale-real(a(2,2),8)/scale
        t=2d0*absb_scaled/(difference_scaled+disc_scaled);c=1d0/sqrt(1d0+t*t)
        u=[-t*c*phase,cmplx(c,0d0,8)]
      endif
    else
      u=(0d0,0d0);u(merge(1,2,real(a(1,1),8)<=real(a(2,2),8)))=1d0
    endif
    scale=max(safe_complex_abs(u(1)),safe_complex_abs(u(2)))
    if(scale<=0d0.or..not.ieee_is_finite(scale))return
    norm=scale*safe_real_hypot(safe_complex_abs(u(1))/scale,safe_complex_abs(u(2))/scale)
    if(norm<=0d0.or..not.ieee_is_finite(norm))return
    u=u/norm;eigenvector=matmul(conjg(transpose(linv)),u)
    ok=all(ieee_is_finite(real(eigenvector,8))).and.all(ieee_is_finite(aimag(eigenvector)))
  end subroutine
#endif

  real(8) function safe_complex_abs(value)
    complex(8),intent(in)::value
    safe_complex_abs=safe_real_hypot(real(value,8),aimag(value))
  end function

  real(8) function safe_real_hypot(x,y)
    real(8),intent(in)::x,y
    real(8)::larger,ratio,factor
    larger=max(abs(x),abs(y))
    if(larger==0d0)then
      safe_real_hypot=0d0
      return
    endif
    ratio=min(abs(x),abs(y))/larger
    factor=sqrt(1d0+ratio*ratio)
    if(larger>huge(1d0)/factor)then
      safe_real_hypot=huge(1d0)
    else
      safe_real_hypot=larger*factor
    endif
  end function

  subroutine safe_complex_divide(value,divisor,result,ok)
    complex(8),intent(in)::value
    real(8),intent(in)::divisor
    complex(8),intent(out)::result
    logical,intent(out)::ok
    real(8)::magnitude
    result=(0d0,0d0);ok=.false.
    if(divisor<=0d0.or..not.ieee_is_finite(divisor))return
    magnitude=safe_complex_abs(value)
    if(divisor<1d0)then
      if(magnitude>huge(1d0)*divisor)return
    endif
    result=value/divisor
    ok=ieee_is_finite(real(result,8))
    if(ok)ok=ieee_is_finite(aimag(result))
  end subroutine

  logical function finite_matrix(matrix)
    complex(8),intent(in)::matrix(:,:)
    integer::i,j
    real(8)::real_part,imaginary_part
    finite_matrix=.true.
    do j=1,size(matrix,2)
      do i=1,size(matrix,1)
        real_part=real(matrix(i,j),8)
        imaginary_part=aimag(matrix(i,j))
        if(real_part/=real_part)then
          finite_matrix=.false.
          return
        endif
        if(imaginary_part/=imaginary_part)then
          finite_matrix=.false.
          return
        endif
        if(abs(real_part)>huge(1d0))then
          finite_matrix=.false.
          return
        endif
        if(abs(imaginary_part)>huge(1d0))then
          finite_matrix=.false.
          return
        endif
      enddo
    enddo
  end function
end module
