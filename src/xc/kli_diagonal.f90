! One spin, equal weights, primitive-cell normalization. i=(band,k).
module kli_diagonal
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: solve_kli
contains
  subroutine solve_kli(p,pq,nk,dv,tol,d,v,rowerr,compat,resid,status)
    real(8),intent(in):: p(:,:),pq(:,:),dv,tol ! pq = Re(conjg(psi_i)*q_i)
    integer,intent(in):: nk
    real(8),intent(out):: d(:),v(:),rowerr,compat,resid
    integer,intent(out):: status
    real(8),allocatable:: a(:,:),b(:),rho(:),sl(:),weight(:),basis(:,:),eig(:),work(:)
    real(8):: lambda
    integer:: ng,n,g,i,j,info
    status=1
    ng=size(p,1); n=size(p,2)
    if(ng<1.or.n<1.or.nk<1)return
    if(mod(n,nk)/=0.or.any(shape(pq)/=shape(p)).or.size(d)/=n.or.size(v)/=ng)return
    if(.not.ieee_is_finite(dv).or..not.ieee_is_finite(tol))return
    if(dv<=0d0.or.tol<=0d0)return
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
    do i=1,n
      b(i)=dv*sum(p(:,i)*sl-pq(:,i))
      do j=1,n
        a(i,j)=-dv*sum(p(:,i)*(weight*p(:,j)))
      enddo
      a(i,i)=a(i,i)+1d0
    enddo
    rowerr=maxval(abs(sum(a,dim=2)));compat=abs(sum(b)/nk)
    if(rowerr>tol.or.compat>tol)return
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
    d=d-sum(d)/n
    resid=maxval(abs(matmul(a,d)-b))
    if(.not.ieee_is_finite(resid).or.resid>tol)return
    v=sl+weight*matmul(p,d)
    do g=1,ng
      if(.not.ieee_is_finite(v(g)))return
    enddo
    status=0
  end subroutine
end module
