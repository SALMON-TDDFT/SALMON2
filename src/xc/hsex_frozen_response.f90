! Low-rank finite-q resolvent. H actions are supplied by the frozen real-space adapter.
module hsex_frozen_response
  implicit none
  private
  public :: frozen_source,frozen_sources,frozen_solve,frozen_trace,s_frozen_solution,s_frozen_stats
  public :: frozen_background_rate
  public :: frozen_density
  type :: s_frozen_solution
    integer :: nshifts=0,nsources=0
    complex(8),allocatable :: left(:,:),right(:,:),coefficient(:,:,:)
    real(8),allocatable :: residual(:)
    ! Window-channel bound: ||response error||_F <= residual*error_scale.
    ! For window channels error_scale=integral(w)*||source||, not a resolvent residual.
    integer,allocatable :: window_kind(:)
    real(8),allocatable :: error_scale(:)
  end type
  type :: s_frozen_stats
    integer :: left_size=0,right_size=0,h_actions=0,iterations=0,precondition_actions=0
  end type
  abstract interface
    subroutine h_action(x,y,status)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  function product(a,b,ac,bc) result(c)
    implicit none
    complex(8),intent(in) :: a(:,:),b(:,:)
    logical,optional,intent(in) :: ac,bc
    complex(8),allocatable :: c(:,:)
    character :: ta,tb
    integer :: m,n,k
    ta='N';tb='N'
    if(present(ac))then
      if(ac)ta='C'
    endif
    if(present(bc))then
      if(bc)tb='C'
    endif
    m=size(a,1);k=size(a,2);n=size(b,2)
    if(ta=='C')then
      m=size(a,2);k=size(a,1)
    endif
    if(tb=='C')n=size(b,1)
    allocate(c(m,n));c=0d0
    if(m==0.or.n==0.or.k==0)return
    call zgemm(ta,tb,m,n,k,cmplx(1d0,0d0,8),a,size(a,1),b,size(b,1),cmplx(0d0,0d0,8),c,m)
  end function

  subroutine frozen_background_rate(psi,hpsi,rate,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    complex(8),intent(in) :: psi(:,:),hpsi(:,:)
    real(8),intent(out) :: rate
    integer,intent(out) :: status
    complex(8),allocatable :: gram(:,:),outside(:,:)
    integer :: a,r
    status=1;rate=0d0;r=size(psi,2)
    if(r<1.or.any(shape(psi)/=shape(hpsi)))return
    if(.not.screen_finite(psi).or..not.screen_finite(hpsi))return
    gram=product(psi,psi,.true.)
    do a=1,r
      gram(a,a)=gram(a,a)-1d0
    enddo
    if(maxval(abs(gram))>1d-8)return
    outside=hpsi-product(psi,product(psi,hpsi,.true.))
    rate=sqrt(2d0*sum(abs(outside)**2)/r)
    if(.not.screen_finite([rate]))return
    status=0
  end subroutine

  subroutine frozen_sources(source,target,fs,ft,vertices,u,v,c,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    complex(8),intent(in) :: source(:,:),target(:,:),vertices(:,:)
    real(8),intent(in) :: fs(:),ft(:)
    complex(8),allocatable,intent(out) :: u(:,:),v(:,:),c(:,:,:)
    integer,intent(out) :: status
    integer :: n,rs,rt,tests,j,a
    n=size(vertices,1);rs=size(source,2);rt=size(target,2);tests=size(vertices,2);status=1
    if(n<1.or.tests<1.or.size(source,1)/=n.or.size(target,1)/=n)return
    if(size(fs)/=rs.or.size(ft)/=rt)return
    if(.not.screen_finite(source).or..not.screen_finite(target).or..not.screen_finite(vertices))return
    if(.not.screen_finite(fs).or..not.screen_finite(ft))return
    if(any(fs<0d0).or.any(ft<0d0).or.any(fs>2d0).or.any(ft>2d0))return
    allocate(u(n,rt+tests*rs),v(n,rs+tests*rt),c(rt+tests*rs,rs+tests*rt,tests));c=0d0
    u(:,:rt)=target;v(:,:rs)=source
    do j=1,tests
      do a=1,rs
        u(:,rt+(j-1)*rs+a)=vertices(:,j)*source(:,a)
        c(rt+(j-1)*rs+a,a,j)=fs(a)
      enddo
      do a=1,rt
        v(:,rs+(j-1)*rt+a)=conjg(vertices(:,j))*target(:,a)
        c(a,rs+(j-1)*rt+a,j)=-ft(a)
      enddo
      if(rs==rt)then
        if(all(source==target).and.all(fs==ft).and.all(vertices(:,j)==vertices(1,j)))c(:,:,j)=0d0
      endif
    enddo
    status=0
  end subroutine

  subroutine frozen_source(source,target,fs,ft,vertex,u,v,c,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    complex(8),intent(in) :: source(:,:),target(:,:),vertex(:)
    real(8),intent(in) :: fs(:),ft(:)
    complex(8),allocatable,intent(out) :: u(:,:),v(:,:),c(:,:)
    integer,intent(out) :: status
    integer :: n,a,rs,rt,r
    status=1;n=size(vertex);rs=size(source,2);rt=size(target,2);r=rs+rt
    if(n<1.or.size(source,1)/=n.or.size(target,1)/=n.or.size(fs)/=rs.or.size(ft)/=rt)return
    if(.not.screen_finite(source).or..not.screen_finite(target).or..not.screen_finite(vertex))return
    if(.not.screen_finite(fs).or..not.screen_finite(ft))return
    if(any(fs<0d0).or.any(ft<0d0).or.any(fs>2d0).or.any(ft>2d0))return
    if(rs==rt)then
      if(all(source==target).and.all(fs==ft).and.all(vertex==vertex(1)))then
        allocate(u(n,0),v(n,0),c(0,0));status=0;return
      endif
    endif
    allocate(u(n,r),v(n,r),c(r,r));c=0d0
    do a=1,rs
      u(:,a)=vertex*source(:,a);v(:,a)=source(:,a);c(a,a)=fs(a)
    enddo
    do a=1,rt
      u(:,rs+a)=target(:,a);v(:,rs+a)=conjg(vertex)*target(:,a);c(rs+a,rs+a)=-ft(a)
    enddo
    status=0
  end subroutine

  subroutine append_vector(q,count,x,accepted)
    implicit none
    complex(8),intent(inout) :: q(:,:)
    integer,intent(inout) :: count
    complex(8),intent(in) :: x(:)
    logical,intent(out) :: accepted
    complex(8) :: y(size(x))
    real(8) :: norm,scale
    integer :: pass,j
    accepted=.false.
    if(count==size(q,2))return
    y=x;scale=sqrt(sum(abs(x)**2))
    do pass=1,2
      do j=1,count
        y=y-q(:,j)*dot_product(q(:,j),y)
      enddo
    enddo
    norm=sqrt(sum(abs(y)**2))
    if(norm<=1d-13*max(scale,tiny(scale)))return
    count=count+1;q(:,count)=y/norm;accepted=.true.
  end subroutine

  subroutine factor(a,r,status)
    implicit none
    complex(8),intent(in) :: a(:,:)
    complex(8),allocatable,intent(out) :: r(:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: matrix(:,:),tau(:),work(:)
    integer :: rows,cols,count,j,i
    rows=size(a,1);cols=size(a,2);count=min(rows,cols)
    allocate(r(count,cols),stat=status)
    if(status/=0)return
    r=0d0
    if(count==0)return
    allocate(matrix(rows,cols),tau(count),work(max(1,64*cols)),stat=status)
    if(status/=0)return
    matrix=a
    call zgeqrf(rows,cols,matrix,rows,tau,work,size(work),status)
    if(status/=0)return
    do j=1,cols
      do i=1,min(j,count)
        r(i,j)=matrix(i,j)
      enddo
    enddo
    ! Only R is needed: ||A C B^H||_F = ||R_A C R_B^H||_F.
    ! Do not form unused unitary Q or truncate nearly dependent columns.
  end subroutine

  subroutine diagonalize(a,values,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    complex(8),intent(inout) :: a(:,:)
    real(8),intent(out) :: values(:)
    integer,intent(out) :: status
    complex(8) :: work(max(1,2*size(a,1)-1))
    real(8) :: rw(max(1,3*size(a,1)-2))
    integer :: n
    n=size(a,1);status=1
    if(.not.screen_finite(a))return
    if(maxval(abs(a-conjg(transpose(a))))>1d-10*max(1d0,maxval(abs(a))))return
    a=.5d0*(a+conjg(transpose(a)))
    call zheev('V','L',n,a,n,values,work,size(work),rw,status)
    if(status==0.and..not.screen_finite(values))status=1
  end subroutine

  logical function source_cancels_exactly(u,c,v)
    implicit none
    complex(8),intent(in) :: u(:,:),c(:,:),v(:,:)
    complex(8),allocatable :: pending(:,:)
    integer :: i,j,k,l,code
    logical :: matched
    source_cancels_exactly=.false.
    allocate(pending(size(c,1),size(c,2)),stat=code)
    if(code/=0)return
    pending=c
    do j=1,size(v,2)
      do i=1,size(u,2)
        if(all(u(:,i)==(0d0,0d0)).or.all(v(:,j)==(0d0,0d0)))pending(i,j)=0d0
      enddo
    enddo
    do j=1,size(v,2)
      do i=1,size(u,2)
        if(pending(i,j)==(0d0,0d0))cycle
        matched=.false.
        search: do l=1,size(v,2)
          if(.not.all(v(:,l)==v(:,j)))cycle
          do k=1,size(u,2)
            if(pending(k,l)/=-pending(i,j))cycle
            if(.not.all(u(:,k)==u(:,i)))cycle
            pending(k,l)=0d0;pending(i,j)=0d0;matched=.true.
            exit search
          enddo
        enddo search
        if(.not.matched)return
      enddo
    enddo
    source_cancels_exactly=.true.
  end function

  subroutine frozen_solve(apply_left,apply_right,u,c,v,z,capacity,tolerance,solution,statistics,status,extra_sources, &
    residual_bound,precondition_left,precondition_right,new_vectors_only,window_kind,window_time)
    use hsex_screen_math, only: screen_finite
    implicit none
    procedure(h_action) :: apply_left,apply_right
    complex(8),intent(in) :: u(:,:),c(:,:),v(:,:),z(:)
    integer,intent(in) :: capacity
    real(8),intent(in) :: tolerance
    type(s_frozen_solution),intent(out) :: solution
    type(s_frozen_stats),intent(out) :: statistics
    integer,intent(out) :: status
    complex(8),optional,intent(in) :: extra_sources(:,:,:)
    procedure(h_action),optional :: precondition_left,precondition_right
    complex(8),allocatable :: expanded(:)
    integer :: expanded_l,expanded_r
    logical,optional,intent(in) :: residual_bound,new_vectors_only
    integer,optional,intent(in) :: window_kind(:)
    real(8),optional,intent(in) :: window_time(:)
    real(8) :: first_moment
    integer :: channel
    logical :: windowed
    complex(8),allocatable :: ql(:,:),qr(:,:),hql(:,:),hqr(:,:),al(:,:),ar(:,:),rhs(:,:),rotated(:,:)
    complex(8),allocatable :: ru(:,:),rv(:,:),rc(:,:),source_ru(:,:),source_rv(:,:)
    complex(8),allocatable :: residual_ru(:,:),residual_rv(:,:),source_c(:,:),small(:,:)
    real(8),allocatable :: el(:),er(:),source_norm(:)
    real(8) :: residual,roundoff,gamma,source_bound
    integer :: nl,nr,cl,cr,done_l,done_r,a,b,s,code,old_l,old_r,cases,case_index,index
    logical :: accepted,bounded,once
    statistics=s_frozen_stats();status=1
    expanded_l=0;expanded_r=0
    if(present(precondition_left).neqv.present(precondition_right))return
    once=.false.
    if(present(new_vectors_only))once=new_vectors_only
    bounded=.false.
    if(present(residual_bound))bounded=residual_bound
    if(size(u,1)<1.or.size(v,1)<1.or.size(z)<1.or.capacity<1)return
    if(size(c,1)/=size(u,2).or.size(c,2)/=size(v,2))return
    if(.not.screen_finite(u).or..not.screen_finite(v).or..not.screen_finite(c))return
    if(.not.screen_finite(z).or..not.screen_finite([tolerance]))return
    if(tolerance<=0d0.or.any(aimag(z)<=0d0))return
    windowed=.false.
    if(present(window_kind).neqv.present(window_time))return
    allocate(solution%window_kind(size(z)),stat=code)
    if(code/=0)return
    solution%window_kind=0
    if(present(window_kind))then
      if(size(window_kind)/=size(z).or.size(window_time)/=size(z))return
      if(any(window_kind<0).or.any(window_kind>2))return
      if(.not.screen_finite(window_time))return
      if(any(window_time<0d0))return
      solution%window_kind=window_kind
      windowed=any(window_kind/=0)
      if(windowed.and..not.bounded)return
      do channel=1,size(z)
        if(window_kind(channel)==0)cycle
        call validate_window(window_kind(channel),window_time(channel),first_moment,code)
        if(code/=0)return
      enddo
    endif
    cases=1
    if(present(extra_sources))then
      if(size(extra_sources,1)/=size(c,1).or.size(extra_sources,2)/=size(c,2))return
      if(.not.screen_finite(extra_sources))return
      cases=cases+size(extra_sources,3)
    endif
    solution%nshifts=size(z);solution%nsources=cases
    allocate(source_norm(cases),stat=code)
    if(code/=0)return
    call factor(u,source_ru,code)
    if(code/=0)return
    call factor(v,source_rv,code)
    if(code/=0)return
    do case_index=1,cases
      source_c=c
      if(case_index>1)source_c=extra_sources(:,:,case_index-1)
      small=product(source_ru,product(source_c,source_rv,bc=.true.))
      source_norm(case_index)=sqrt(sum(abs(small)**2))
      source_bound=sqrt(sum(abs(source_ru)**2)*sum(abs(source_c)**2)*sum(abs(source_rv)**2))
      if(source_bound>0d0.and.source_norm(case_index)<=64d0*epsilon(source_bound)*source_bound)then
        ! A small compressed source remains unresolved unless its ORIGINAL
        ! factors cancel in exact opposite pairs. Never infer zero from QR
        ! rank truncation or a small numerical norm.
        if(.not.source_cancels_exactly(u,source_c,v))return
        source_norm(case_index)=0d0
      endif
    enddo
    if(.not.screen_finite(source_norm))return
    allocate(solution%error_scale(size(z)*cases),stat=code)
    if(code/=0)return
    do case_index=1,cases
      do channel=1,size(z)
        index=(case_index-1)*size(z)+channel
        solution%error_scale(index)=source_norm(case_index)
        if(solution%window_kind(channel)/=0)then
          solution%error_scale(index)=source_norm(case_index)*window_time(channel)
        endif
      enddo
    enddo
    if(.not.screen_finite(solution%error_scale))return
    if(maxval(source_norm)==0d0)then
      allocate(solution%left(size(u,1),0),solution%right(size(v,1),0),solution%coefficient(0,0,size(z)*cases))
      allocate(solution%residual(size(z)*cases));solution%residual=0d0;status=0;return
    endif
    cl=min(capacity,size(u,1));cr=min(capacity,size(v,1))
    allocate(ql(size(u,1),cl),qr(size(v,1),cr),hql(size(u,1),cl),hqr(size(v,1),cr),stat=code)
    if(code/=0)return
    ql=0d0;qr=0d0;nl=0;nr=0;done_l=0;done_r=0
    do a=1,size(u,2)
      call append_vector(ql,nl,u(:,a),accepted)
    enddo
    do a=1,size(v,2)
      call append_vector(qr,nr,v(:,a),accepted)
    enddo
    if(nl==0.or.nr==0)return
    do
      statistics%iterations=statistics%iterations+1
      do a=done_l+1,nl
        call apply_left(ql(:,a),hql(:,a),code)
        statistics%h_actions=statistics%h_actions+1
        if(code/=0)then
          status=3;return
        endif
        if(.not.screen_finite(hql(:,a)))then
          status=3;return
        endif
      enddo
      do a=done_r+1,nr
        call apply_right(qr(:,a),hqr(:,a),code)
        statistics%h_actions=statistics%h_actions+1
        if(code/=0)then
          status=3;return
        endif
        if(.not.screen_finite(hqr(:,a)))then
          status=3;return
        endif
      enddo
      done_l=nl;done_r=nr;statistics%left_size=nl;statistics%right_size=nr
      al=product(ql(:,:nl),hql(:,:nl),.true.)
      ar=product(qr(:,:nr),hqr(:,:nr),.true.)
      allocate(el(nl),er(nr),stat=code)
      if(code/=0)return
      call diagonalize(al,el,code)
      if(code/=0)return
      call diagonalize(ar,er,code)
      if(code/=0)return
      solution%left=ql(:,:nl);solution%right=qr(:,:nr)
      if(allocated(solution%coefficient))deallocate(solution%coefficient,solution%residual)
      allocate(solution%coefficient(nl,nr,size(z)*cases),solution%residual(size(z)*cases),stat=code)
      if(code/=0)return
      if(bounded)then
        call spectral_residual_bound(ql(:,:nl),hql(:,:nl),qr(:,:nr),hqr(:,:nr),al,ar,el,er, &
          u,c,v,z,source_norm,solution,code,extra_sources,window_time)
        if(code/=0)return
      else
      allocate(ru(size(u,1),2*nl+size(u,2)),rv(size(v,1),2*nr+size(v,2)),stat=code)
      if(code/=0)return
      ru(:,:nl)=ql(:,:nl);ru(:,nl+1:2*nl)=hql(:,:nl);ru(:,2*nl+1:)=u
      rv(:,:nr)=qr(:,:nr);rv(:,nr+1:2*nr)=hqr(:,:nr);rv(:,2*nr+1:)=v
      allocate(rc(size(ru,2),size(rv,2)),stat=code)
      if(code/=0)return
      ! Residual factor spaces are independent of both source and spectral shift.
      call factor(ru,residual_ru,code)
      if(code/=0)return
      call factor(rv,residual_rv,code)
      if(code/=0)return
      do case_index=1,cases
        source_c=c
        if(case_index>1)source_c=extra_sources(:,:,case_index-1)
        rhs=product(product(ql(:,:nl),u,.true.),product(source_c,product(v,qr(:,:nr),.true.)))
        rotated=product(al,product(rhs,ar),.true.)
        do s=1,size(z)
          index=(case_index-1)*size(z)+s
          if(source_norm(case_index)==0d0)then
            solution%coefficient(:,:,index)=0d0;solution%residual(index)=0d0
            cycle
          endif
          rhs=rotated
          do b=1,nr
            do a=1,nl
              rhs(a,b)=rhs(a,b)/(z(s)-el(a)+er(b))
            enddo
          enddo
          solution%coefficient(:,:,index)=product(al,product(rhs,ar,bc=.true.))
          rc=0d0
          rc(:nl,:nr)=z(s)*solution%coefficient(:,:,index)
          rc(nl+1:2*nl,:nr)=-solution%coefficient(:,:,index)
          rc(:nl,nr+1:2*nr)=solution%coefficient(:,:,index)
          rc(2*nl+1:,2*nr+1:)=-source_c
          small=product(residual_ru,product(rc,residual_rv,bc=.true.))
          residual=sqrt(sum(abs(small)**2))/source_norm(case_index)
          ! Retain the previous, conservative bound despite the smaller factors.
          gamma=16d0*epsilon(residual)*max(size(u,1),size(v,1),3*nl+size(u,2),3*nr+size(v,2))
          roundoff=gamma*(sqrt(sum(abs(solution%coefficient(:,:,index))**2))* &
            (abs(z(s))+sqrt(sum(abs(hql(:,:nl))**2))+sqrt(sum(abs(hqr(:,:nr))**2)))+ &
            sqrt(sum(abs(u)**2)*sum(abs(source_c)**2)*sum(abs(v)**2)))/source_norm(case_index)
          solution%residual(index)=residual+roundoff
        enddo
      enddo
      deallocate(ru,rv,rc,residual_ru,residual_rv)
      endif
      deallocate(el,er)
      if(.not.screen_finite(solution%residual).or..not.screen_finite(solution%coefficient))then
        status=1;return
      endif
      if(maxval(solution%residual)<=tolerance)then
        status=0;return
      endif
      old_l=nl;old_r=nr
      if(present(precondition_left))then
        allocate(expanded(size(ql,1)))
        do a=expanded_l+1,done_l
          call precondition_left(hql(:,a),expanded,code)
          statistics%precondition_actions=statistics%precondition_actions+1
          if(code/=0)then
            status=3;return
          endif
          if(.not.screen_finite(expanded))then
            status=3;return
          endif
          call append_vector(ql,nl,expanded,accepted)
        enddo
        deallocate(expanded);allocate(expanded(size(qr,1)))
        do a=expanded_r+1,done_r
          call precondition_right(hqr(:,a),expanded,code)
          statistics%precondition_actions=statistics%precondition_actions+1
          if(code/=0)then
            status=3;return
          endif
          if(.not.screen_finite(expanded))then
            status=3;return
          endif
          call append_vector(qr,nr,expanded,accepted)
        enddo
        deallocate(expanded)
      else
        ! A dependent old candidate cannot become independent as the space
        ! grows. Pair workers consume each H-basis vector once; the legacy
        ! order remains available for reproducibility comparisons.
        do a=merge(expanded_l+1,1,once),done_l
          call append_vector(ql,nl,hql(:,a),accepted)
        enddo
        do a=merge(expanded_r+1,1,once),done_r
          call append_vector(qr,nr,hqr(:,a),accepted)
        enddo
      endif
      expanded_l=done_l;expanded_r=done_r
      if(old_l==nl.and.old_r==nr)then
        status=2;return
      endif
    enddo
  end subroutine

  subroutine validate_window(kind,area,first,status)
    use hsex_frozen_window, only: window_moments
    implicit none
    integer,intent(in) :: kind
    real(8),intent(in) :: area
    real(8),intent(out) :: first
    integer,intent(out) :: status
    call window_moments(kind,area,first,status)
  end subroutine

  subroutine spectral_residual_bound(ql,hql,qr,hqr,al,ar,el,er,u,c,v,z,source_norm,solution,status,extra,window_time)
    use hsex_frozen_window, only: window_filter_matrix,window_moments
    use hsex_screen_math, only: screen_finite
    implicit none
    complex(8),intent(in) :: ql(:,:),hql(:,:),qr(:,:),hqr(:,:),al(:,:),ar(:,:),u(:,:),c(:,:),v(:,:),z(:)
    real(8),intent(in) :: el(:),er(:),source_norm(:)
    type(s_frozen_solution),intent(inout) :: solution
    integer,intent(out) :: status
    complex(8),optional,intent(in) :: extra(:,:,:)
    real(8),optional,intent(in) :: window_time(:)
    complex(8),allocatable :: filters(:,:,:)
    real(8) :: area,first,rhsnorm,numerator
    complex(8),allocatable :: defect(:,:),ls(:,:),rs(:,:),sc(:,:),rhs(:,:),d(:,:)
    real(8) :: lnorm,rnorm,ldef,rdef,uout,vout,vnorm,lnorm_source,snorm,cnorm,source_error,gamma,roundoff
    integer :: i,j,a,b,s,index,nl,nr
    status=0;nl=size(ql,2);nr=size(qr,2)
    solution%left=product(ql,al);solution%right=product(qr,ar)
    defect=product(solution%left,solution%left,.true.)
    do a=1,nl
      defect(a,a)=defect(a,a)-1d0
    enddo
    lnorm=sqrt(1d0+sqrt(sum(abs(defect)**2)))
    defect=product(solution%right,solution%right,.true.)
    do a=1,nr
      defect(a,a)=defect(a,a)-1d0
    enddo
    rnorm=sqrt(1d0+sqrt(sum(abs(defect)**2)))
    defect=product(hql,al)
    do a=1,nl
      defect(:,a)=defect(:,a)-solution%left(:,a)*el(a)
    enddo
    ldef=sqrt(sum(abs(defect)**2))
    defect=product(hqr,ar)
    do a=1,nr
      defect(:,a)=defect(:,a)-solution%right(:,a)*er(a)
    enddo
    rdef=sqrt(sum(abs(defect)**2))
    ls=product(solution%left,u,.true.);rs=product(v,solution%right,.true.)
    defect=u-product(solution%left,ls);uout=sqrt(sum(abs(defect)**2))
    defect=v-product(solution%right,rs,bc=.true.);vout=sqrt(sum(abs(defect)**2))
    vnorm=sqrt(sum(abs(v)**2));lnorm_source=sqrt(sum(abs(ls)**2))
    gamma=16d0*epsilon(1d0)*max(size(u,1),size(v,1),3*nl+size(u,2),3*nr+size(v,2))
    if(any(solution%window_kind/=0))then
      allocate(filters(nl,nr,size(z)),stat=status)
      if(status/=0)return
      do s=1,size(z)
        if(solution%window_kind(s)==0)cycle
        call window_filter_matrix(solution%window_kind(s),window_time(s),el,er,filters(:,:,s),status)
        if(status/=0)return
      enddo
    endif
    do i=1,size(source_norm)
      sc=c
      if(i>1)sc=extra(:,:,i-1)
      snorm=sqrt(sum(abs(sc)**2))
      rhs=product(ls,product(sc,rs))
      source_error=snorm*(uout*vnorm+lnorm*lnorm_source*vout)
      do s=1,size(z)
        index=(i-1)*size(z)+s
        if(source_norm(i)==0d0)then
          solution%coefficient(:,:,index)=0d0;solution%residual(index)=0d0
          cycle
        endif
        if(solution%window_kind(s)/=0)then
          area=window_time(s)
          call window_moments(solution%window_kind(s),area,first,status)
          if(status/=0)return
          solution%coefficient(:,:,index)=rhs*filters(:,:,s)
          rhsnorm=sqrt(sum(abs(rhs)**2))
          ! Duhamel bound for the two unitary propagators, integrated over
          ! the causal window. The first time moment controls H-basis defects.
          numerator=area*source_error+first*rhsnorm*(ldef*rnorm+rdef*lnorm)
          roundoff=gamma*(area*(lnorm*rnorm*rhsnorm+ &
            sqrt(sum(abs(u)**2)*sum(abs(sc)**2)*sum(abs(v)**2)))+ &
            first*rhsnorm*(sqrt(sum(abs(hql)**2))+sqrt(sum(abs(hqr)**2))))
          ! Reserve an additional per-entry filter-evaluation budget. This
          ! covers quadrature/asymptotic floating-point error in addition to
          ! the basis/BLAS roundoff allowance above.
          roundoff=roundoff+1024d0*epsilon(1d0)*area*sqrt(dble(nl)*nr)*rhsnorm*lnorm*rnorm
          if(.not.screen_finite([numerator,roundoff]).or.solution%error_scale(index)<=0d0)then
            status=1;return
          endif
          solution%residual(index)=(numerator+roundoff)/solution%error_scale(index)
          cycle
        endif
        do b=1,nr
          do a=1,nl
            solution%coefficient(a,b,index)=rhs(a,b)/(z(s)-el(a)+er(b))
          enddo
        enddo
        d=solution%coefficient(:,:,index)
        cnorm=sqrt(sum(abs(d)**2))
        do b=1,nr
          do a=1,nl
            d(a,b)=(z(s)-el(a)+er(b))*d(a,b)-rhs(a,b)
          enddo
        enddo
        ! R = L D R^H - (HL-L eL) C R^H + L C (HR-R eR)^H - source_out.
        ! Triangle/submultiplicative inequalities give an upper bound, never
        ! cancellation of large Gram terms. No source or H-defect is truncated.
        roundoff=gamma*(cnorm*(abs(z(s))+sqrt(sum(abs(hql)**2))+sqrt(sum(abs(hqr)**2)))+ &
          sqrt(sum(abs(u)**2)*sum(abs(sc)**2)*sum(abs(v)**2)))/source_norm(i)
        solution%residual(index)=(lnorm*rnorm*sqrt(sum(abs(d)**2))+ &
          cnorm*(ldef*rnorm+rdef*lnorm)+source_error)/source_norm(i)+roundoff
      enddo
    enddo
  end subroutine

  subroutine frozen_trace(solution,vertex,values,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    type(s_frozen_solution),intent(in) :: solution
    complex(8),intent(in) :: vertex(:)
    complex(8),intent(out) :: values(:)
    integer,intent(out) :: status
    complex(8),allocatable :: weighted(:,:),matrix(:,:)
    integer :: a,s
    status=1;values=0d0
    if(.not.allocated(solution%coefficient))return
    if(size(solution%left,1)/=size(vertex).or.size(solution%right,1)/=size(vertex))return
    if(size(values)/=size(solution%coefficient,3).or..not.screen_finite(vertex))return
    weighted=solution%left
    do a=1,size(weighted,2)
      weighted(:,a)=conjg(vertex)*weighted(:,a)
    enddo
    matrix=product(solution%right,weighted,.true.)
    do s=1,size(values)
      values(s)=sum(transpose(matrix)*solution%coefficient(:,:,s))
    enddo
    if(.not.screen_finite(values))return
    status=0
  end subroutine

  subroutine frozen_density(solution,density,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    type(s_frozen_solution),intent(in) :: solution
    complex(8),intent(out) :: density(:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: weighted(:,:)
    integer :: s,b
    status=1;density=0d0
    if(.not.allocated(solution%coefficient))return
    if(size(density,1)/=size(solution%left,1).or.size(solution%right,1)/=size(density,1))return
    if(size(density,2)/=size(solution%coefficient,3))return
    do s=1,size(density,2)
      weighted=product(solution%left,solution%coefficient(:,:,s))
      do b=1,size(weighted,2)
        density(:,s)=density(:,s)+weighted(:,b)*conjg(solution%right(:,b))
      enddo
    enddo
    if(.not.screen_finite(density))return
    status=0
  end subroutine
end module
