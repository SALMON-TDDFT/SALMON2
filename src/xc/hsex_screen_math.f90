! Pure numerical pieces of the microscopic static screening calculation.
module hsex_screen_math
  implicit none
  private
  public :: screen_sites,screen_geometry,screen_forms,screen_assemble,screen_channels
  public :: screen_metal_head,screen_transform,screen_population,screen_finite
  interface screen_finite
    module procedure finite_real,finite_complex,finite_real2,finite_complex2,finite_complex3
  end interface
contains
  logical function finite_real(a)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in)::a(:)
    integer::i
    finite_real=.false.
    do i=1,size(a)
      if(.not.ieee_is_finite(a(i)))return
    enddo
    finite_real=.true.
  end function

  logical function finite_complex(a)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    complex(8),intent(in)::a(:)
    integer::i
    finite_complex=.false.
    do i=1,size(a)
      if(.not.ieee_is_finite(real(a(i),8)).or..not.ieee_is_finite(aimag(a(i))))return
    enddo
    finite_complex=.true.
  end function

  logical function finite_real2(a)
    implicit none
    real(8),intent(in)::a(:,:)
    integer::j
    finite_real2=.false.
    do j=1,size(a,2)
      if(.not.finite_real(a(:,j)))return
    enddo
    finite_real2=.true.
  end function

  logical function finite_complex2(a)
    implicit none
    complex(8),intent(in)::a(:,:)
    integer::j
    finite_complex2=.false.
    do j=1,size(a,2)
      if(.not.finite_complex(a(:,j)))return
    enddo
    finite_complex2=.true.
  end function

  logical function finite_complex3(a)
    implicit none
    complex(8),intent(in)::a(:,:,:)
    integer::j,k
    finite_complex3=.false.
    do k=1,size(a,3)
      do j=1,size(a,2)
        if(.not.finite_complex(a(:,j,k)))return
      enddo
    enddo
    finite_complex3=.true.
  end function

  subroutine screen_sites(mesh,sites)
    implicit none
    integer,intent(in)::mesh(3)
    integer,intent(out)::sites(3,product(mesh))
    integer::i,a
    do i=0,product(mesh)-1
      sites(:,i+1)=[mod(i,mesh(1)),mod(i/mesh(1),mesh(2)),i/(mesh(1)*mesh(2))]
      do a=1,3
        if(2*sites(a,i+1)>=mesh(a))sites(a,i+1)=sites(a,i+1)-mesh(a)
      enddo
    enddo
  end subroutine

  subroutine screen_geometry(cell,rec,vol,status)
    implicit none
    real(8),intent(in)::cell(3,3)
    real(8),intent(out)::rec(3,3),vol
    integer,intent(out)::status
    real(8)::det
    status=1
    if(.not.screen_finite(cell))return
    rec(:,1)=[cell(2,2)*cell(3,3)-cell(3,2)*cell(2,3), &
      cell(3,2)*cell(1,3)-cell(1,2)*cell(3,3),cell(1,2)*cell(2,3)-cell(2,2)*cell(1,3)]
    rec(:,2)=[cell(2,3)*cell(3,1)-cell(3,3)*cell(2,1), &
      cell(3,3)*cell(1,1)-cell(1,3)*cell(3,1),cell(1,3)*cell(2,1)-cell(2,3)*cell(1,1)]
    rec(:,3)=[cell(2,1)*cell(3,2)-cell(3,1)*cell(2,2), &
      cell(3,1)*cell(1,2)-cell(1,1)*cell(3,2),cell(1,1)*cell(2,2)-cell(2,1)*cell(1,2)]
    det=sum(cell(:,1)*rec(:,1));vol=abs(det)
    if(vol<1d-20)return
    rec=rec*(2d0*acos(-1d0)/det)
    status=0
  end subroutine

  subroutine screen_forms(n,boxes,q,forms)
    implicit none
    integer,intent(in)::n(3),boxes
    real(8),intent(in)::q(:,:)
    complex(8),intent(out)::forms(boxes**3,size(q,2))
    integer::i,g,a,b,j,m
    complex(8)::factor
    real(8)::angle
    forms=1d0
    do g=1,size(q,2)
      do i=0,boxes**3-1
        do a=1,3
          b=mod(i/boxes**(a-1),boxes);m=n(a)/boxes
          factor=0d0
          do j=b*m,(b+1)*m-1
            angle=-2d0*acos(-1d0)*j*q(a,g)/n(a)
            factor=factor+cmplx(cos(angle),sin(angle),8)
          enddo
          forms(i+1,g)=forms(i+1,g)*factor/m
        enddo
      enddo
    enddo
  end subroutine

  subroutine screen_assemble(n,boxes,q,rec,vol,nk,chi,pump,mat,head,qhead,status)
    implicit none
    integer,intent(in)::n(3),boxes,nk
    real(8),intent(in)::q(:,:),rec(3,3),vol
    complex(8),intent(in)::chi(:,:)
    logical,intent(in)::pump
    complex(8),intent(out)::mat(boxes**3,boxes**3)
    real(8),intent(out)::head,qhead
    integer,intent(out)::status
    complex(8),allocatable::eps(:,:),inv(:,:),forms(:,:),tmp(:,:),work(:)
    real(8),allocatable::b(:),q2(:),levels(:),rwork(:)
    integer,allocatable::piv(:)
    integer::g,i,j,l,na,info
    g=size(q,2);na=boxes**3;status=1;head=0d0;qhead=0d0
    if(g<1.or.vol<=0d0.or.nk<1.or.boxes<1)return
    if(any(mod(n,boxes)/=0))return
    if(.not.screen_finite(chi))return
    if(.not.screen_finite(q))return
    if(.not.screen_finite(rec).or..not.screen_finite([vol]))return
    allocate(eps(g,g),inv(g,g),forms(na,g),tmp(na,g),b(g),q2(g),piv(g))
    do i=1,g
      q2(i)=sum(matmul(rec,q(:,i))**2)
      b(i)=0d0
      if(q2(i)>=1d-14)b(i)=sqrt(4d0*acos(-1d0)/q2(i))
    enddo
    do j=1,g
      do i=1,g
        eps(i,j)=2d0/(nk*vol)*b(i)*chi(i,j)*b(j)
      enddo
      eps(j,j)=eps(j,j)+1d0
    enddo
    if(.not.screen_finite(eps))return
    inv=eps
    if(pump)then
      allocate(levels(g),work(max(1,2*g)),rwork(max(1,3*g-2)))
      call zheev('N','L',g,inv,g,levels,work,size(work),rwork,info)
      if(info/=0)return
      if(minval(levels)<=1d-10)return
    endif
    inv=0d0
    do i=1,g
      inv(i,i)=1d0
    enddo
    call zgesv(g,g,eps,g,piv,inv,g,info)
    if(info/=0)return
    if(all(q2>=1d-14))then
      l=minloc(q2,dim=1);head=real(inv(l,l),8);qhead=q2(l)
    endif
    do j=1,g
      inv(:,j)=inv(:,j)*b*b(j)
    enddo
    call screen_forms(n,boxes,q,forms)
    call zgemm('N','N',na,g,g,cmplx(1d0,0d0,8),conjg(forms),na,inv,g,cmplx(0d0,0d0,8),tmp,na)
    call zgemm('N','T',na,na,g,cmplx(1d0/vol,0d0,8),tmp,na,forms,na,cmplx(0d0,0d0,8),mat,na)
    if(.not.screen_finite(mat))return
    status=0
  end subroutine

  subroutine screen_population(levels,pop,nv,insulating,status)
    implicit none
    real(8),intent(in)::levels(:,:),pop(:,:)
    integer,intent(in)::nv
    logical,intent(out)::insulating
    integer,intent(out)::status
    integer::nb
    status=1;insulating=.false.;nb=size(levels,1)
    if(any(shape(levels)/=shape(pop)).or.nv<1.or.nv>nb)return
    if(.not.screen_finite(levels))return
    if(.not.screen_finite(pop))return
    if(any(levels(2:,:)-levels(:nb-1,:)<-1d-10))return
    if(minval(pop)<-1d-10.or.maxval(pop)>2d0+1d-7)return
    if(any(abs(sum(pop,dim=1)-2*nv)>1d-4))return
    if(nb>nv)then
      insulating=all(abs(pop(:nv,:)-2d0)<=1d-7).and.all(abs(pop(nv+1:,:))<=1d-7)
      insulating=insulating.and.minval(levels(nv+1:,:))-maxval(levels(:nv,:))>0d0
    endif
    status=0
  end subroutine

  subroutine screen_channels(a,b,error,aggregate)
    implicit none
    complex(8),intent(in)::a(:,:,:),b(:,:,:)
    real(8),intent(out)::error
    logical,intent(in)::aggregate
    real(8)::norma(4),normb(4),diff(4),floor,local
    complex(8),allocatable::ar(:),al(:),br(:),bl(:)
    complex(8)::ah,bh
    integer::n,j,lo,hi,block,c
    n=size(a,1);error=huge(1d0)
    if(any(shape(a)/=shape(b)).or.size(a,2)/=n)return
    if(.not.screen_finite(a).or..not.screen_finite(b))return
    allocate(ar(n),al(n),br(n),bl(n));error=0d0
    do block=1,merge(1,size(a,3),aggregate)
      lo=block;hi=block
      if(aggregate)then
        lo=1;hi=size(a,3)
      endif
      norma=0d0;normb=0d0;diff=0d0
      do j=lo,hi
        ar=sum(a(:,:,j),dim=2)/sqrt(dble(n));al=sum(a(:,:,j),dim=1)/sqrt(dble(n))
        br=sum(b(:,:,j),dim=2)/sqrt(dble(n));bl=sum(b(:,:,j),dim=1)/sqrt(dble(n))
        ah=sum(ar)/sqrt(dble(n));bh=sum(br)/sqrt(dble(n))
        ar=ar-ah/sqrt(dble(n));al=al-ah/sqrt(dble(n))
        br=br-bh/sqrt(dble(n));bl=bl-bh/sqrt(dble(n))
        norma=norma+[sum(abs(a(:,:,j))**2),abs(ah)**2,sum(abs(ar)**2),sum(abs(al)**2)]
        normb=normb+[sum(abs(b(:,:,j))**2),abs(bh)**2,sum(abs(br)**2),sum(abs(bl)**2)]
        diff=diff+[sum(abs(a(:,:,j)-b(:,:,j))**2),abs(ah-bh)**2,sum(abs(ar-br)**2),sum(abs(al-bl)**2)]
      enddo
      norma=sqrt(norma);normb=sqrt(normb);diff=sqrt(diff)
      floor=max(256d0*epsilon(1d0)*max(norma(1),normb(1)),1d-14)
      do c=1,4
        local=0d0
        if(diff(c)>floor)local=diff(c)/max(normb(c),floor)
        error=max(error,local)
      enddo
    enddo
  end subroutine

  subroutine screen_metal_head(values,qmesh,rec,mesh,tolerance,point,average,status,sample_indices)
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    complex(8),intent(in)::values(:,:,:)
    integer,intent(in)::qmesh(3),mesh(3)
    real(8),intent(in)::rec(3,3),tolerance
    complex(8),intent(out)::point(:,:),average(:,:)
    integer,intent(out)::status
    integer,optional,intent(in)::sample_indices(:)
    integer,allocatable::sites(:,:),valid(:),idx(:),lookup(:)
    real(8),allocatable::q(:,:),r2(:),sorted(:),design(:,:),s(:),rwork(:),eval(:)
    complex(8),allocatable::a(:,:),b(:,:),work(:),pred(:,:,:),target(:,:,:),av(:,:,:),fit(:,:,:),check(:,:)
    real(8)::radius,scale,cov(3,3),moment(10),error,errmax,norm,detrec
    integer::nq,na,nvalid,i,j,l,m,pass,rank,info,counts(2),lw
    status=1;nq=product(qmesh);na=size(values,1)
    if(any(qmesh<6).or.any(mesh<qmesh).or.tolerance<=0d0.or.tolerance>=1d0)return
    if(size(values,2)/=na.or.na<1)return
    if(.not.screen_finite(rec).or..not.screen_finite([tolerance]))return
    detrec=rec(1,1)*(rec(2,2)*rec(3,3)-rec(3,2)*rec(2,3)) &
      -rec(1,2)*(rec(2,1)*rec(3,3)-rec(3,1)*rec(2,3)) &
      +rec(1,3)*(rec(2,1)*rec(3,2)-rec(3,1)*rec(2,2))
    if(.not.screen_finite([detrec]))return
    if(abs(detrec)<1d-20)return
    allocate(lookup(nq));lookup=0
    if(present(sample_indices))then
      if(size(sample_indices)/=size(values,3))return
      do i=1,size(sample_indices)
        if(sample_indices(i)<1.or.sample_indices(i)>nq)return
        if(lookup(sample_indices(i))/=0)return
        lookup(sample_indices(i))=i
      enddo
    else
      if(size(values,3)/=nq)return
      do i=1,nq
        lookup(i)=i
      enddo
    endif
    if(.not.screen_finite(values))return
    allocate(sites(3,nq),q(3,nq),r2(nq),valid(nq),sorted(nq),design(nq,10),av(na,na,2),fit(na,na,2))
    call screen_sites(qmesh,sites);nvalid=0
    do i=1,nq
      q(:,i)=matmul(rec,dble(sites(:,i))/qmesh);r2(i)=sum(q(:,i)**2)
      if(r2(i)>0d0.and.all(2*abs(sites(:,i))<qmesh))then
        nvalid=nvalid+1;valid(nvalid)=i;sorted(nvalid)=r2(i)
      endif
    enddo
    if(nvalid<52)return
    do i=2,nvalid
      radius=sorted(i);j=i-1
      do while(j>=1)
        if(sorted(j)<=radius)exit
        sorted(j+1)=sorted(j);j=j-1
      enddo
      sorted(j+1)=radius
    enddo
    radius=sorted(52)*(1d0+1d-10);scale=0d0
    do i=1,nvalid
      if(r2(valid(i))<=radius)scale=max(scale,sqrt(r2(valid(i))))
    enddo
    q=q/scale
    do i=1,nq
      design(i,:)=[1d0,q(:,i),q(:,i)**2,q(1,i)*q(2,i),q(1,i)*q(3,i),q(2,i)*q(3,i)]
    enddo
    cov=0d0
    do i=1,3
      do j=1,3
        cov(i,j)=sum(rec(i,:)*rec(j,:)/dble(mesh)**2)/(12d0*scale**2)
      enddo
    enddo
    moment=[1d0,0d0,0d0,0d0,cov(1,1),cov(2,2),cov(3,3),cov(1,2),cov(1,3),cov(2,3)]
    counts=[26,52];errmax=0d0
    do pass=1,2
      radius=sorted(counts(pass))*(1d0+1d-10)
      m=count(r2(valid(:nvalid))<=radius)
      allocate(idx(m),a(m,10),b(max(m,10),na*na),s(10),rwork(5*10))
      l=0
      do i=1,nvalid
        if(r2(valid(i))>radius)cycle
        l=l+1;idx(l)=valid(i)
        if(lookup(idx(l))==0)return
      enddo
      a=cmplx(design(idx,:),0d0,8);b=0d0
      do i=1,m
        b(i,:)=reshape(values(:,:,lookup(idx(i))),[na*na])
      enddo
      allocate(work(1))
      status=2
      call zgelss(m,10,na*na,a,m,b,max(m,10),s,-1d0,rank,work,-1,rwork,info)
      if(info/=0)return
      lw=max(1,int(real(work(1),8)));deallocate(work);allocate(work(lw))
      call zgelss(m,10,na*na,a,m,b,max(m,10),s,-1d0,rank,work,lw,rwork,info)
      if(info/=0)return
      status=3
      if(rank/=10)return
      if(s(1)/s(10)>1d6)return
      fit(:,:,pass)=reshape(b(1,:),[na,na])
      av(:,:,pass)=reshape(matmul(moment,b(:10,:)),[na,na])
      allocate(pred(na,na,m),target(na,na,m))
      do i=1,m
        pred(:,:,i)=reshape(matmul(design(idx(i),:),b(:10,:)),[na,na]);target(:,:,i)=values(:,:,lookup(idx(i)))
      enddo
      call screen_channels(pred,target,error,.true.);errmax=max(errmax,error)
      write(error_unit,'(a,i0,a,i0,a,es16.8)')'HSEX_HEAD_FIT pass=',pass,' samples=',m,' relative_error=',error
      deallocate(idx,a,b,s,rwork,work,pred,target)
    enddo
    call screen_channels(av(:,:,2:2),av(:,:,1:1),error,.true.);errmax=max(errmax,error)
    status=4
    write(error_unit,'(a,3es16.8)')'HSEX_HEAD_ERRORS nested/max/tolerance ',error,errmax,tolerance
    if(errmax>tolerance)return
    allocate(check(na,na),eval(na),work(max(1,2*na)),rwork(max(1,3*na-2)))
    do pass=1,2
      check=fit(:,:,1)
      if(pass==2)check=av(:,:,1)
      norm=max(sqrt(sum(abs(check)**2)),1d-14)
      status=5
      if(sqrt(sum(aimag(check)**2))>1d-8*norm)return
      if(sqrt(sum(abs(check-transpose(conjg(check)))**2))>1d-8*norm)return
      status=6
      call zheev('N','L',na,check,na,eval,work,size(work),rwork,info)
      if(info/=0)return
      status=7
      if(minval(eval)<-1d-8*norm)then
        write(error_unit,'(a,i0,a,2es16.8)')'HSEX_HEAD_NEGATIVE pass=',pass,' eigenvalue/scale ',minval(eval),norm
        return
      endif
    enddo
    point=real(fit(:,:,1),8);average=real(av(:,:,1),8);status=0
  end subroutine

  subroutine fft_line(x)
    implicit none
    complex(8),intent(inout)::x(:)
    complex(8)::tmp,factor,w
    complex(8),allocatable::y(:)
    integer::n,i,j,m,l,k
    n=size(x)
    if(iand(n,n-1)/=0)then
      allocate(y(n));y=0d0
      do j=0,n-1
        do i=0,n-1
          y(j+1)=y(j+1)+x(i+1)*exp(cmplx(0d0,-2d0*acos(-1d0)*i*j/n,8))
        enddo
      enddo
      x=y;return
    endif
    j=0
    do i=0,n-1
      if(i<j)then
        tmp=x(i+1);x(i+1)=x(j+1);x(j+1)=tmp
      endif
      m=n/2
      do while(m>=1)
        if(j<m)exit
        j=j-m;m=m/2
      enddo
      j=j+m
    enddo
    l=2
    do while(l<=n)
      factor=exp(cmplx(0d0,-2d0*acos(-1d0)/l,8))
      do k=1,n,l
        w=1d0
        do j=0,l/2-1
          tmp=w*x(k+j+l/2);x(k+j+l/2)=x(k+j)-tmp;x(k+j)=x(k+j)+tmp;w=w*factor
        enddo
      enddo
      l=l*2
    enddo
  end subroutine

  subroutine screen_transform(values,mesh)
    implicit none
    integer,intent(in)::mesh(3)
    complex(8),intent(inout)::values(:,:)
    complex(8),allocatable::grid(:,:,:)
    integer::p,i,j
    !$omp parallel do default(shared) private(grid,p,i,j) schedule(static)
    do p=1,size(values,1)
      allocate(grid(mesh(1),mesh(2),mesh(3)));grid=reshape(values(p,:),mesh)
      do j=1,mesh(3)
        do i=1,mesh(2)
          call fft_line(grid(:,i,j))
        enddo
      enddo
      do j=1,mesh(3)
        do i=1,mesh(1)
          call fft_line(grid(i,:,j))
        enddo
      enddo
      do j=1,mesh(2)
        do i=1,mesh(1)
          call fft_line(grid(i,j,:))
        enddo
      enddo
      values(p,:)=reshape(grid,[product(mesh)])/product(mesh);deallocate(grid)
    enddo
    !$omp end parallel do
  end subroutine
end module
