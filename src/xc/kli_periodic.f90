! Primitive cell: integral |psi_nk|^2 = 1, per spin. Cell-periodic SALMON orbitals.
module kli_periodic
  use exx_wannier, only: s_exx_wannier,wannier_init,wannier_localize,wannier_forward,wannier_backward
  use kli_self_hartree
  use kli_diagonal
  use builtin_pz_sp, only: exc_cor_pz_sp
  implicit none
  private
  public :: kli_frame,evaluate_kli
  type kli_frame
    type(s_exx_wannier)::w
    type(isolated_box)::box,check_box
    real(8)::audit(7)=0d0
    integer,allocatable::image_shift(:,:)
  end type
contains
  subroutine evaluate_kli(f,psi,n,mesh,cell,k,dv,ltol,ktol,htol,maxiter,v,energy,locerr,edge,iter,nfft,status)
    type(kli_frame),intent(inout)::f
    complex(8),intent(in)::psi(:,:,:)
    integer,intent(in)::n(3),mesh(3),maxiter
    real(8),intent(in)::cell(3,3),k(:,:),dv,ltol,ktol,htol
    real(8),intent(out)::v(:),energy,locerr,edge
    integer,intent(out)::iter,nfft,status
    integer:: ng,nb,nk,ns(3),sz,ik,i,j,ic,ios,ls
    real(8)::step(3,3),alpha,eold,err,compat,resid,eh,bnd
    real(8),allocatable::vcheck(:,:,:),rho(:,:,:)
    real(8),allocatable:: kred(:,:),occ(:,:),eps(:),exc(:),vxc(:,:),p(:,:),pq(:,:),d(:)
    complex(8),allocatable::w(:,:),uw(:,:),rot(:,:,:),uwk(:,:,:),gen(:,:,:),old(:,:,:),h(:,:),t(:,:),eye(:,:),q(:,:)
    integer,allocatable::ipiv(:)
    complex(8),allocatable::prior(:,:,:),gu(:,:,:)
    ng=size(psi,1);nb=size(psi,2);nk=size(psi,3);ns=n*mesh;sz=product(ns)
    status=1;nfft=0;iter=0;locerr=huge(1d0);edge=0;energy=0;f%audit=0
    if(product(n)/=ng.or.product(mesh)/=nk.or.size(v)/=ng.or.nb<1)return
    allocate(kred(3,nk),occ(sz,2),eps(sz),exc(sz),vxc(sz,2),p(ng,nb*nk),pq(ng,nb*nk),d(nb*nk))
    allocate(w(sz,nb),uw(sz,nb),rot(ng,nb,nk),uwk(ng,nb,nk),gen(nb,nb,nk),old(nb,nb,nk))
    allocate(h(nb,nb),t(nb,nb),eye(nb,nb),q(ng,nb),ipiv(nb))
    eye=0
    do i=1,nb
      eye(i,i)=1
    enddo
    do ik=1,nk
      if(maxval(abs(dv*matmul(conjg(transpose(psi(:,:,ik))),psi(:,:,ik))-eye))>ktol)return
    enddo
    if(f%w%ng==0)then
      ! Use fractional grid coordinates in the existing discrete Bloch transform.
      ! The physical volume and Coulomb metric come from the actual lattice.
      kred=matmul(transpose(cell),k)
      call move_alloc(f%w%previous,prior);call move_alloc(f%w%gauge,gu)
      call wannier_init(f%w,n,mesh,1d0/n,kred,0d0,ios)
      if(ios/=0)return
      call move_alloc(prior,f%w%previous);call move_alloc(gu,f%w%gauge)
      f%w%dv=dv
      if(.not.allocated(f%image_shift))then
        allocate(f%image_shift(3,nb));f%image_shift=huge(1)
      endif
      do i=1,3
        step(:,i)=cell(:,i)/n(i)
      enddo
      call prepare_box(f%box,ns,step,ios)
      if(ios/=0)return
      call prepare_box(f%check_box,ns,step,ios,1.25d0)
      if(ios/=0)return
    endif
    call wannier_localize(f%w,psi,0,ltol,ios)
    if(ios/=0)return
    call orbitals(.true.)
    if(ios/=0)return
    do
      do ik=1,nk
        h=dv*matmul(conjg(transpose(rot(:,:,ik))),uwk(:,:,ik))
        gen(:,:,ik)=h-conjg(transpose(h))
      enddo
      ! Fourier coefficients are the independent pair/translation stationarity defects.
      locerr=0
      do ic=1,nk
        h=0
        do ik=1,nk
          h=h+gen(:,:,ik)*exp(cmplx(0d0,sum((f%w%k(:,ik)-f%w%k(:,1))*f%w%cell(:,ic)/n),8))/nk
        enddo
        locerr=max(locerr,maxval(abs(h)))
      enddo
      if(locerr<=ltol)exit
      if(iter>=maxiter)return
      old=f%w%gauge;eold=energy;alpha=4d0
      do ls=1,24
        do ik=1,nk
          ! Unitary Cayley step along negative energy gradient.
          h=eye+.5d0*alpha*gen(:,:,ik);t=eye-.5d0*alpha*gen(:,:,ik)
          call zgesv(nb,nb,h,nb,ipiv,t,nb,ios)
          if(ios/=0)return
          f%w%gauge(:,:,ik)=matmul(old(:,:,ik),t)
        enddo
        call orbitals(.false.)
        if(ios/=0)return
        if(energy<eold)exit
        alpha=alpha/2
      enddo
      if(ls>24)return
      iter=iter+1
    enddo
    do ik=1,nk
      q=matmul(uwk(:,:,ik),conjg(transpose(f%w%gauge(:,:,ik))))
      do j=1,nb
        i=j+(ik-1)*nb
        pq(:,i)=real(conjg(psi(:,j,ik))*q(:,j),8)
        p(:,i)=abs(psi(:,j,ik))**2
      enddo
    enddo
    ! All orbital pairs and all BvK translations, independently of the KLI system.
    do ic=1,nk
      h=0
      do ik=1,nk
        h=h+dv*matmul(conjg(transpose(rot(:,:,ik))),rot(:,:,ik))* &
          exp(cmplx(0d0,sum((f%w%k(:,ik)-f%w%k(:,1))*f%w%cell(:,ic)/n),8))/nk
      enddo
      if(ic==1)h=h-eye
      f%audit(1)=max(f%audit(1),maxval(abs(h)))
    enddo
    if(f%audit(1)>ktol)return
    ! Compare the final physical potential with a larger auxiliary FFT volume;
    ! no new orbital samples enter this independent box-size test.
    allocate(vcheck(ns(1),ns(2),ns(3)),rho(ns(1),ns(2),ns(3)))
    do j=1,nb
      occ(:,1)=abs(w(:,j))**2;occ(:,2)=0d0;rho=reshape(occ(:,1),ns)
      call self_hartree(f%check_box,rho,vcheck,eh,bnd,ios,f%image_shift(:,j))
      if(ios/=0)return
      call exc_cor_pz_sp(sz,occ,eps,exc,vxc)
      ! RMS potential difference weighted with the unit orbital density.
      f%audit(7)=max(f%audit(7),sqrt(dv*sum(abs(uw(:,j)+(reshape(vcheck,[sz])+vxc(:,1))*w(:,j))**2)))
      nfft=nfft+2
    enddo
    if(f%audit(7)>htol)return
    ! Sum of translated orbital densities must match the primitive spin density.
    eps(:ng)=0d0
    do i=1,sz
      j=f%w%primitive_point(i)
      eps(j)=eps(j)+sum(abs(w(i,:))**2)
    enddo
    f%audit(2)=maxval(abs(eps(:ng)-sum(p,dim=2)/nk))
    eps(:ng)=0d0
    do i=1,sz
      j=f%w%primitive_point(i)
      eps(j)=eps(j)+real(dot_product(w(i,:),uw(i,:)),8)
    enddo
    f%audit(3)=maxval(abs(eps(:ng)-sum(pq,dim=2)/nk))
    if(maxval(f%audit(2:3))>ktol)return
    call solve_kli(p,pq,nk,dv,ktol,d,v,err,compat,resid,ios)
    if(ios/=0)return
    f%audit(4:6)=[err,compat,resid]
    f%w%previous=rot
    status=0
  contains
    subroutine orbitals(track)
      logical,intent(in)::track
      integer::a,l
      real(8)::eh,bnd
      real(8),allocatable::rho(:,:,:),vh(:,:,:)
      allocate(rho(ns(1),ns(2),ns(3)),vh(ns(1),ns(2),ns(3)))
      do l=1,nk
        rot(:,:,l)=matmul(psi(:,:,l),f%w%gauge(:,:,l))
      enddo
      call wannier_forward(f%w,rot,w)
      edge=0;energy=0
      do a=1,nb
        occ(:,1)=abs(w(:,a))**2;occ(:,2)=0d0
        rho=reshape(occ(:,1),ns)
        call self_hartree(f%box,rho,vh,eh,bnd,ios,f%image_shift(:,a),track)
        if(ios/=0)return
        nfft=nfft+2;edge=max(edge,bnd)
        call exc_cor_pz_sp(sz,occ,eps,exc,vxc)
        energy=energy-eh-dv*sum(exc)
        uw(:,a)=-(reshape(vh,[sz])+vxc(:,1))*w(:,a)
      enddo
      call wannier_backward(f%w,uw,uwk)
    end subroutine
  end subroutine
end module
