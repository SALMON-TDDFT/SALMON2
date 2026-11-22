program test_frozen_response
  use hsex_frozen_response, only: frozen_source,frozen_sources,frozen_solve,frozen_trace,s_frozen_solution,s_frozen_stats
  use hsex_frozen_response, only: frozen_density
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  integer,parameter :: n=6,nv=2,ns=2
  complex(8) :: hl(n,n),hr(n,n),pl(n,nv),pr(n,nv),vertex(n),probe(n),z(ns),expected(ns),actual(ns)
  complex(8) :: left_vectors(n,n),right_vectors(n,n),weighted(n,n)
  complex(8),allocatable :: u(:,:),v(:,:),c(:,:)
  complex(8),allocatable :: extra(:,:,:)
  complex(8) :: batch_values(3*ns),second(ns)
  complex(8) :: batch_vertices(n,3),batch_expected(3*ns)
  complex(8) :: density_expected(n,ns),density_actual(n,ns)
  real(8) :: occ(nv),angle
  type(s_frozen_solution) :: solution
  type(s_frozen_stats) :: stats
  integer :: i,j,status,coherent,precondition_fault
  hl=0d0;hr=0d0;occ=2d0
  do i=1,n
    vertex(i)=cmplx(cos(.37d0*i),sin(.37d0*i),8)
    probe(i)=cmplx(cos(.21d0*i),-.2d0*sin(.21d0*i),8)
    do j=1,n
      left_vectors(i,j)=exp(cmplx(0d0,2d0*acos(-1d0)*(i-1)*(j-1)/n,8))/sqrt(dble(n))
      right_vectors(i,j)=left_vectors(i,j)*exp(cmplx(0d0,.17d0*i,8))
    enddo
  enddo
  do j=1,n
    weighted(:,j)=left_vectors(:,j)*(.17d0*j)
  enddo
  hl=matmul(weighted,conjg(transpose(left_vectors)))
  do j=1,n
    weighted(:,j)=right_vectors(:,j)*(.13d0*j+.02d0)
  enddo
  hr=matmul(weighted,conjg(transpose(right_vectors)))
  z=[cmplx(0d0,.08d0,8),cmplx(.21d0,.04d0,8)]
  do coherent=0,1
    pl=left_vectors(:,:nv);pr=right_vectors(:,:nv)
    if(coherent==1)then
      angle=.31d0
      pl(:,1)=cos(angle)*left_vectors(:,1)+cmplx(0d0,sin(angle),8)*left_vectors(:,4)
      pr(:,2)=cos(.19d0)*right_vectors(:,2)+sin(.19d0)*right_vectors(:,5)
    endif
    call frozen_source(pr,pl,occ,occ,vertex,u,v,c,status)
    if(status/=0)error stop 'valid source rejected'
    call direct_reference(u,c,v,expected,density_expected)
    if(maxval(abs(expected))<1d-3)error stop 'model accidentally has zero response'
    call frozen_solve(left_action,right_action,u,c,v,z,n,1d-10,solution,stats,status)
    if(status/=0)error stop 'full-space response not converged'
    call check_explicit_residual()
    if(maxval(solution%residual)>1d-10)error stop 'false convergence'
    call frozen_trace(solution,probe,actual,status)
    if(status/=0.or.maxval(abs(actual-expected))>2d-10)error stop 'Liouvillian reference mismatch'
    call frozen_density(solution,density_actual,status)
    if(status/=0.or.maxval(abs(density_actual-density_expected))>2d-10)error stop 'full induced density mismatch'
    print '(a,i2,2es16.7,3i6)','FROZEN_MODEL coherent/error/residual/bases/actions ',coherent, &
      maxval(abs(actual-expected)),maxval(solution%residual),stats%left_size,stats%right_size,stats%h_actions
    if(allocated(extra))deallocate(extra)
    allocate(extra(size(c,1),size(c,2),2));extra=0d0
    extra(:,:,1)=.7d0*c;extra(1,2,1)=cmplx(.13d0,-.07d0,8)
    call direct_reference(u,extra(:,:,1),v,second)
    call frozen_solve(left_action,right_action,u,c,v,z,n,1d-10,solution,stats,status,extra_sources=extra)
    if(status/=0.or.stats%h_actions/=2*n)error stop 'batch basis not shared'
    call frozen_trace(solution,probe,batch_values,status)
    if(maxval(abs(batch_values(:ns)-expected))>2d-10.or. &
       maxval(abs(batch_values(ns+1:2*ns)-second))>2d-10.or. &
       maxval(abs(batch_values(2*ns+1:)))>1d-13)error stop 'shared sources mismatch'
    batch_vertices(:,1)=vertex;batch_vertices(:,2)=probe;batch_vertices(:,3)=1d0
    call frozen_sources(pr,pl,occ,occ,batch_vertices,u,v,extra,status)
    if(status/=0)error stop 'joint source factory rejected'
    do j=1,3
      call direct_reference(u,extra(:,:,j),v,batch_expected((j-1)*ns+1:j*ns))
    enddo
    call frozen_solve(left_action,right_action,u,extra(:,:,1),v,z,n,1d-10,solution,stats,status, &
      extra_sources=extra(:,:,2:))
    if(status/=0.or.stats%h_actions/=2*n)error stop 'joint sources failed'
    call frozen_trace(solution,probe,batch_values,status)
    if(maxval(abs(batch_values-batch_expected))>2d-10)error stop 'joint factor response mismatch'
    call frozen_source(pr,pl,occ,occ,vertex,u,v,c,status)
    call frozen_solve(left_action,right_action,u,c,v,z,n,1d-10,solution,stats,status,residual_bound=.true.)
    if(status/=0)error stop 'conservative spectral residual bound did not converge'
    call frozen_trace(solution,probe,actual,status)
    if(status/=0.or.maxval(abs(actual-expected))>2d-10)error stop 'spectral-bound response mismatch'
    call check_explicit_residual()
    call frozen_solve(left_action,right_action,u,c,v,z,n,1d-10,solution,stats,status, &
      residual_bound=.true.,new_vectors_only=.true.)
    if(status/=0)error stop 'one-pass expansion response rejected'
    call check_explicit_residual()
    call frozen_trace(solution,probe,actual,status)
    if(status/=0.or.maxval(abs(actual-expected))>2d-10)error stop 'one-pass expansion mismatch'
    call frozen_solve(left_action,right_action,u,c,v,z,n,1d-10,solution,stats,status, &
      residual_bound=.true.,precondition_left=left_action,precondition_right=right_action)
    if(status/=0)error stop 'preconditioned response rejected'
    call check_explicit_residual()
    call frozen_trace(solution,probe,actual,status)
    if(status/=0.or.maxval(abs(actual-expected))>2d-10)error stop 'preconditioned response mismatch'
    do precondition_fault=0,2
      call frozen_solve(left_action,right_action,u,c,v,z,n,1d-10,solution,stats,status, &
        precondition_left=bad_precondition,precondition_right=right_action)
      if(status/=3)error stop 'failed/NaN/Inf precondition silently accepted'
    enddo
    call frozen_solve(left_action,right_action,u,c,v,z,n,1d-10,solution,stats,status, &
      precondition_left=left_action)
    if(status/=1)error stop 'unpaired precondition silently accepted'

    call frozen_solve(left_action,right_action,u,c,v,z,2,1d-12,solution,stats,status,residual_bound=.true.)
    if(status/=2)error stop 'spectral-bound capacity exhaustion silently accepted'
    call check_explicit_residual(.true.)
    call frozen_solve(left_action,right_action,u,c,v,z,2,1d-12,solution,stats,status)
    if(status/=2)error stop 'capacity exhaustion silently accepted'
    call check_explicit_residual()
  enddo
  ! A scalar potential commutes with the same initial density, including coherent density.
  call frozen_source(pl,pl,occ,occ,cmplx(1d0,0d0,8)*[(1d0,i=1,n)],u,v,c,status)
  call frozen_solve(left_action,left_action,u,c,v,z,n,1d-10,solution,stats,status)
  if(status/=0)error stop 'zero source rejected'
  call frozen_trace(solution,probe,actual,status)
  if(maxval(abs(actual))>1d-12)error stop 'uniform-potential response'
  vertex(1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call frozen_source(pr,pl,occ,occ,vertex,u,v,c,status)
  if(status==0)error stop 'NaN source accepted'
  vertex(1)=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
  call frozen_source(pr,pl,occ,occ,vertex,u,v,c,status)
  if(status==0)error stop 'Infinity source accepted'
  vertex=probe
  call frozen_source(pr,pl,occ,occ,vertex,u,v,c,status)
  call frozen_solve(bad_action,right_action,u,c,v,z,n,1d-10,solution,stats,status)
  if(status/=3)error stop 'callback failure not propagated'
  call frozen_solve(nan_action,right_action,u,c,v,z,n,1d-10,solution,stats,status)
  if(status/=3)error stop 'NaN Hamiltonian action accepted'
  z(1)=cmplx(0d0,-.1d0,8)
  call frozen_solve(left_action,right_action,u,c,v,z,n,1d-10,solution,stats,status)
  if(status/=1)error stop 'advanced shift accepted'
  call cancellation_regression()
  print *,'FROZEN_RESPONSE_DENSE_COHERENT_BATCH_SHIFT_CAP_IEEE_PASS'
contains
  subroutine bad_precondition(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    code=1
    if(precondition_fault==2)return
    y=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
    if(precondition_fault==1)y=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
    code=0
  end subroutine
  subroutine check_explicit_residual(loose)
    implicit none
    logical,optional,intent(in) :: loose
    complex(8) :: x(n,n),r(n,n),source(n,n)
    real(8) :: norm,measured
    integer :: s
    logical :: tight
    tight=.true.
    if(present(loose))tight=.not.loose
    source=matmul(u,matmul(c,conjg(transpose(v))))
    norm=sqrt(sum(abs(source)**2))
    do s=1,ns
      x=matmul(solution%left,matmul(solution%coefficient(:,:,s),conjg(transpose(solution%right))))
      r=z(s)*x-matmul(hl,x)+matmul(x,hr)-source
      measured=sqrt(sum(abs(r)**2))/norm
      if(solution%residual(s)<measured-5d-13)error stop 'residual underestimated'
      if(tight.and.solution%residual(s)>measured+1d-8)error stop 'residual contraction mismatch'
    enddo
  end subroutine
  subroutine cancellation_regression()
    implicit none
    complex(8) :: su(n,2),sv(n,2),sc(2,2),small_shift(1)
    su=0d0;sv=0d0;sc=0d0
    su(1,:)=1d0;su(2,2)=1d-14;sv(1,:)=1d0;sc(1,1)=-1d0;sc(2,2)=1d0
    small_shift=cmplx(0d0,.1d0,8)
    call frozen_solve(left_action,right_action,su,sc,sv,small_shift,n,1d-10,solution,stats,status)
    if(status==0)error stop 'unresolved near-cancelled source falsely accepted'
    su(2,2)=0d0
    call frozen_solve(bad_action,bad_action,su,sc,sv,small_shift,n,1d-10,solution,stats,status)
    if(status/=0.or.stats%h_actions/=0)error stop 'exactly zero original source rejected'
    call frozen_trace(solution,probe,actual(:1),status)
    if(status/=0.or.any(actual(:1)/=(0d0,0d0)))error stop 'exact-zero source has a response'
    su=0d0;sv=0d0;sc=0d0
    su(:,1)=pl(:,1);sv(:,1)=pl(:,1);sc(1,1)=1d0
    small_shift=cmplx(0d0,1d-15,8)
    call frozen_solve(left_action,left_action,su,sc,sv,small_shift,n,1d-10,solution,stats,status)
    if(status==0)error stop 'resonant roundoff falsely accepted'
  end subroutine
  subroutine bad_action(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    y=0d0;code=1
  end subroutine
  subroutine nan_action(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    code=1
    if(precondition_fault==2)return
    y=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
    if(precondition_fault==1)y=cmplx(ieee_value(0d0,ieee_positive_inf),0d0,8)
    code=0
  end subroutine
  subroutine left_action(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    y=matmul(hl,x);code=0
  end subroutine
  subroutine right_action(x,y,code)
    implicit none
    complex(8),intent(in) :: x(:)
    complex(8),intent(out) :: y(:)
    integer,intent(out) :: code
    y=matmul(hr,x);code=0
  end subroutine
  subroutine direct_reference(su,sc,sv,values,densities)
    implicit none
    complex(8),intent(in) :: su(:,:),sc(:,:),sv(:,:)
    complex(8),intent(out) :: values(:)
    complex(8),optional,intent(out) :: densities(:,:)
    complex(8) :: a(n*n,n*n),b(n*n,1),source(n,n)
    integer :: p,q,r,s,shift,info,piv(n*n)
    source=matmul(su,matmul(sc,conjg(transpose(sv))))
    do shift=1,ns
      a=0d0;b(:,1)=reshape(source,[n*n])
      do q=1,n
        do p=1,n
          r=p+n*(q-1);a(r,r)=z(shift)
          do s=1,n
            a(r,s+n*(q-1))=a(r,s+n*(q-1))-hl(p,s)
            a(r,p+n*(s-1))=a(r,p+n*(s-1))+hr(s,q)
          enddo
        enddo
      enddo
      call zgesv(n*n,1,a,n*n,piv,b,n*n,info)
      if(info/=0)error stop 'independent reference solve failed'
      values(shift)=(0d0,0d0)
      do p=1,n
        values(shift)=values(shift)+conjg(probe(p))*b(p+n*(p-1),1)
        if(present(densities))densities(p,shift)=b(p+n*(p-1),1)
      enddo
    enddo
  end subroutine
end program
