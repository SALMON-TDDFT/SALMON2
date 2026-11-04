! One spin, equal weights, primitive-cell normalization. i=(band,k).
module kli_diagonal
  implicit none
  private
  public :: solve_kli
contains
  subroutine solve_kli(p,pq,nk,dv,tol,d,v,rowerr,compat,resid,status)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in):: p(:,:),pq(:,:),dv,tol ! pq = Re(conjg(psi_i)*q_i)
    integer,intent(in):: nk
    real(8),intent(out):: d(:),v(:),rowerr,compat,resid
    integer,intent(out):: status
    real(8),allocatable:: a(:,:),b(:),rho(:),sl(:),weight(:),basis(:,:),eig(:),work(:)
    real(8),allocatable:: weighted_p(:,:),scale(:)
    real(8):: lambda,gap_bound
    integer:: ng,n,g,i,j,info
    status=1
    ng=size(p,1); n=size(p,2)
    if(ng<1.or.n<1.or.nk<1)return
    if(mod(n,nk)/=0.or.any(shape(pq)/=shape(p)).or.size(d)/=n.or.size(v)/=ng)return
    if(.not.ieee_is_finite(dv).or..not.ieee_is_finite(tol))return
    if(dv<=0d0.or.tol<=0d0)return
    if(n>2*ng)then
      call solve_in_grid_space
      return
    endif
    allocate(a(n,n),b(n),rho(ng),sl(ng),weight(ng),basis(n,n),eig(n),work(max(1,3*n)))
    do i=1,n
      do g=1,ng
        if(.not.ieee_is_finite(p(g,i)).or..not.ieee_is_finite(pq(g,i)))return
        if(p(g,i)<0d0)return
      enddo
    enddo
    rho=sum(p,dim=2)/nk
    sl=0d0;weight=0d0
    do g=1,ng
      if(rho(g)>tiny(1d0))then
        weight(g)=1d0/(nk*rho(g))
        sl(g)=sum(pq(g,:))*weight(g)
      else if(any(pq(g,:)/=0d0))then
        return
      endif
    enddo
    ! A = I - X^T X; form each symmetric overlap once with blocked BLAS.
    allocate(weighted_p(ng,n),scale(ng))
    scale=sqrt(dv*weight)
    do i=1,n
      weighted_p(:,i)=scale*p(:,i)
      b(i)=dv*sum(p(:,i)*sl-pq(:,i))
    enddo
    call dsyrk('U','T',n,ng,-1d0,weighted_p,ng,0d0,a,n)
    do j=1,n
      a(j+1:n,j)=a(j,j+1:n)
    enddo
    gap_bound=-n*maxval(a)
    do i=1,n
      a(i,i)=a(i,i)+1d0
    enddo
    rowerr=maxval(abs(sum(a,dim=2)));compat=abs(sum(b)/nk)
    if(rowerr>tol.or.compat>tol)return
    ! On the constant-free subspace the weighted graph Laplacian has
    ! Rayleigh quotient >= n * min_ij (X^T X)_ij. Row-sum defects
    ! perturb that bound by at most rowerr. If the bound is inconclusive,
    ! use the spectral test to distinguish weakly connected and singular cases.
    if(gap_bound-rowerr>max(tol,100*epsilon(dv)))then
      basis=a+1d0/n
      d=b
      call dposv('U',n,1,basis,n,d,n,info)
      if(info/=0)return
    else
      basis=a
      call dsyev('V','U',n,basis,n,eig,work,size(work),info)
      if(info/=0)return
      ! A single constant is arbitrary. Reject disconnected occupied components.
      if(abs(eig(1))>tol)return
      d=0d0
      do i=2,n
        lambda=eig(i)
        if(lambda<=max(tol,100*epsilon(lambda)))return
        d=d+basis(:,i)*dot_product(basis(:,i),b)/lambda
      enddo
    endif
    d=d-sum(d)/n
    resid=maxval(abs(matmul(a,d)-b))
    if(.not.ieee_is_finite(resid).or.resid>tol)return
    v=sl+weight*matmul(p,d)
    do g=1,ng
      if(.not.ieee_is_finite(v(g)))return
    enddo
    status=0
  contains
    subroutine solve_in_grid_space
      implicit none
      ! M = X^T X. The nonzero eigenvalues are also those of X X^T.
      ! D = b + X^T (I-X X^T)^+ X b on the constant-free subspace.
      real(8),allocatable::xmat(:,:),gram(:,:),vec(:),projected(:),ones(:)
      allocate(xmat(ng,n),gram(ng,ng),vec(ng),projected(ng),ones(ng))
      allocate(b(n),rho(ng),sl(ng),weight(ng),eig(ng),work(3*ng))
      do j=1,n
        do g=1,ng
          if(.not.ieee_is_finite(p(g,j)).or..not.ieee_is_finite(pq(g,j)).or.p(g,j)<0)return
        enddo
      enddo
      rho=sum(p,dim=2)/nk;sl=0;weight=0
      do g=1,ng
        if(rho(g)>tiny(dv))then
          weight(g)=1d0/(nk*rho(g))
          sl(g)=sum(pq(g,:))*weight(g)
        else if(any(pq(g,:)/=0d0))then
          return
        endif
      enddo
      do j=1,n
        xmat(:,j)=sqrt(dv*weight)*p(:,j)
        b(j)=dv*sum(p(:,j)*sl-pq(:,j))
      enddo
      ones=sum(xmat,dim=2)
      rowerr=maxval(abs(matmul(ones,xmat)-1d0));compat=abs(sum(b)/nk)
      if(rowerr>tol.or.compat>tol)return
      call dsyrk('U','N',ng,n,-1d0,xmat,ng,0d0,gram,ng)
      do g=1,ng
        gram(g,g)=gram(g,g)+1
      enddo
      call dsyev('V','U',ng,gram,ng,eig,work,size(work),info)
      if(info/=0.or.abs(eig(1))>tol)return
      vec=matmul(xmat,b);projected=0
      do i=2,ng
        if(eig(i)<=max(tol,100*epsilon(dv)))return
        projected=projected+gram(:,i)*(dot_product(gram(:,i),vec)/eig(i))
      enddo
      d=b+matmul(projected,xmat);d=d-sum(d)/n
      resid=maxval(abs(d-matmul(matmul(xmat,d),xmat)-b))
      if(.not.ieee_is_finite(resid).or.resid>tol)return
      v=sl+weight*matmul(p,d)
      do g=1,ng
        if(.not.ieee_is_finite(v(g)))return
      enddo
      status=0
    end subroutine
  end subroutine
end module
