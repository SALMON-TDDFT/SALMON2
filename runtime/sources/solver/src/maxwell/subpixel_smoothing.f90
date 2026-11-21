! Reference-frequency reciprocal subpixel material coefficients.
! Time evolution is unchanged: freeze interface correctors, not epsilon(omega).
module subpixel_smoothing
  use iso_c_binding
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::subpixel_mix,plane_fraction,shape_cut,subpixel_spd_bound
  public::sp_subpixel_mix,sp_subpixel_plane_fraction,sp_subpixel_shape_cut
  public::s_subpixel_profile,read_subpixel_profile,allocate_subpixel_profile,subpixel_bulk,subpixel_at
  type s_subpixel_profile
    integer::nm=0,ns=0,nq=0,maxiter=200
    real(c_double)::omega=0d0,rtol=1d-10
    real(c_double),allocatable::background(:,:,:),strength(:,:,:),resonance(:),damping(:)
    real(c_double),allocatable::shape_size(:,:),center(:,:)
    integer,allocatable::owner(:),kind(:),material(:),ld_phase(:)
  end type
contains
  logical function real_matrix_ok(a,strict)result(ok)
    real(c_double),intent(in)::a(3,3)
    logical,intent(in)::strict
    real(c_double)::b(3,3),scale,tol,det
    integer::i,j
    ok=.false.;scale=0d0;tol=128d0*epsilon(1d0)
    do j=1,3
      do i=1,3
        if(.not.ieee_is_finite(a(i,j)))return
        scale=max(scale,abs(a(i,j)))
      enddo
    enddo
    if(scale==0d0)then
      ok=.not.strict;return
    endif
    b=a/scale
    do i=1,3
      if(b(i,i)<0d0)return
      if(strict.and.b(i,i)<=0d0)return
      do j=1,3
        if(abs(b(i,j)-b(j,i))>tol)return
      enddo
    enddo
    do i=1,2
      do j=i+1,3
        if(b(i,i)*b(j,j)-b(i,j)*b(j,i)<-tol)return
      enddo
    enddo
    det=b(1,1)*(b(2,2)*b(3,3)-b(2,3)*b(3,2)) &
       -b(1,2)*(b(2,1)*b(3,3)-b(2,3)*b(3,1)) &
       +b(1,3)*(b(2,1)*b(3,2)-b(2,2)*b(3,1))
    if(det < -tol)return
    if(strict.and.(det<=0d0.or.b(1,1)*b(2,2)-b(1,2)*b(2,1)<=0d0))return
    ok=.true.
  end function

  real(c_double) function subpixel_spd_bound(a)result(bound)
    real(c_double),intent(in)::a(3,3)
    real(c_double)::b(3,3),adj(3,3),det,scale
    bound=0d0
    if(.not.real_matrix_ok(a,.true.))return
    scale=maxval(abs(a));b=a/scale
    adj(1,:)=[b(2,2)*b(3,3)-b(2,3)*b(3,2),b(1,3)*b(3,2)-b(1,2)*b(3,3), &
      b(1,2)*b(2,3)-b(1,3)*b(2,2)]
    adj(2,:)=[b(2,3)*b(3,1)-b(2,1)*b(3,3),b(1,1)*b(3,3)-b(1,3)*b(3,1), &
      b(1,3)*b(2,1)-b(1,1)*b(2,3)]
    adj(3,:)=[b(2,1)*b(3,2)-b(2,2)*b(3,1),b(1,2)*b(3,1)-b(1,1)*b(3,2), &
      b(1,1)*b(2,2)-b(1,2)*b(2,1)]
    det=sum(b(1,:)*adj(:,1))
    if(det<=0d0)return
    bound=scale*det/maxval(sum(abs(adj),dim=2))
    if(.not.ieee_is_finite(bound).or.bound<=0d0)bound=0d0
  end function

  integer function subpixel_mix(nq,omega,fraction,normal,background,strength,resonance,damping,owner, &
      effective_background,effective_strength,target,matched)result(status)
    integer,intent(in)::nq,owner(nq)
    real(c_double),intent(in)::omega,fraction,normal(3),background(3,3,2),strength(3,3,nq)
    real(c_double),intent(in)::resonance(nq),damping(nq)
    real(c_double),intent(out)::effective_background(3,3),effective_strength(3,3,nq)
    complex(c_double_complex),intent(out)::target(3,3),matched(3,3)
    real(c_double)::f(2),n(3),p(3,3),identity(3,3),scale,length,a(3,3)
    complex(c_double_complex)::eps(3,3,2),t(3,3,2),den,an(2),bn(3,2),dn(3),inverse_mean
    complex(c_double_complex)::work(3,3),field(3,3)
    integer::i,j,m,q
    status=1;effective_background=0d0;effective_strength=0d0;target=0d0;matched=0d0
    if(nq<0.or..not.ieee_is_finite(omega).or.omega<=0d0)return
    if(.not.ieee_is_finite(fraction).or.fraction<0d0.or.fraction>1d0)return
    do i=1,3
      if(.not.ieee_is_finite(normal(i)))return
    enddo
    scale=maxval(abs(normal))
    if(scale==0d0)return
    length=sqrt(sum((normal/scale)**2));n=(normal/scale)/length
    f=[fraction,1d0-fraction];identity=0d0
    do i=1,3
      identity(i,i)=1d0
    enddo
    do j=1,3
      do i=1,3
        p(i,j)=identity(i,j)-n(i)*n(j)
      enddo
    enddo
    do m=1,2
      if(.not.real_matrix_ok(background(:,:,m),.true.))return
      eps(:,:,m)=cmplx(background(:,:,m),0d0,c_double)
    enddo
    do q=1,nq
      if(owner(q)<1.or.owner(q)>2)return
      if(.not.ieee_is_finite(resonance(q)).or.resonance(q)<0d0)return
      if(.not.ieee_is_finite(damping(q)).or.damping(q)<0d0)return
      if(.not.real_matrix_ok(strength(:,:,q),.false.))return
      if(maxval(abs(strength(:,:,q)))==0d0)cycle
      if(resonance(q)==omega.and.damping(q)==0d0)then
        status=2;return
      endif
      den=cmplx(resonance(q)**2-omega**2,-damping(q)*omega,c_double)
      if(abs(den)==0d0)then
        status=2;return
      endif
      eps(:,:,owner(q))=eps(:,:,owner(q))+strength(:,:,q)/den
    enddo
    status=2
    do m=1,2
      do j=1,3
        do i=1,3
          if(.not.ieee_is_finite(real(eps(i,j,m),c_double)))return
          if(.not.ieee_is_finite(aimag(eps(i,j,m))))return
        enddo
      enddo
    enddo
    if(fraction==0d0.or.fraction==1d0)then
      m=merge(1,2,fraction==1d0)
      effective_background=background(:,:,m);target=eps(:,:,m);matched=target
      do q=1,nq
        if(owner(q)==m)effective_strength(:,:,q)=strength(:,:,q)
      enddo
      status=0;return
    endif
    inverse_mean=0d0;dn=0d0
    do m=1,2
      an(m)=sum(n*matmul(eps(:,:,m),n))
      if(abs(an(m))<=tiny(1d0))return
      bn(:,m)=matmul(n,matmul(eps(:,:,m),p))
      inverse_mean=inverse_mean+f(m)/an(m)
      dn=dn+f(m)*bn(:,m)/an(m)
    enddo
    if(abs(inverse_mean)<=tiny(1d0))return
    dn=(cmplx(n,0d0,c_double)+dn)/inverse_mean
    do m=1,2
      do j=1,3
        do i=1,3
          t(i,j,m)=p(i,j)+n(i)*(dn(j)-bn(j,m))/an(m)
        enddo
      enddo
      target=target+f(m)*matmul(eps(:,:,m),t(:,:,m))
      field=matmul(cmplx(background(:,:,m),0d0,c_double),t(:,:,m))
      work=matmul(transpose(conjg(t(:,:,m))),field)
      a=real(work,c_double)
      effective_background=effective_background+f(m)*.5d0*(a+transpose(a))
    enddo
    do q=1,nq
      m=owner(q)
      field=matmul(cmplx(strength(:,:,q),0d0,c_double),t(:,:,m))
      work=matmul(transpose(conjg(t(:,:,m))),field)
      a=real(work,c_double)
      effective_strength(:,:,q)=f(m)*.5d0*(a+transpose(a))
    enddo
    if(.not.real_matrix_ok(effective_background,.true.))return
    matched=cmplx(effective_background,0d0,c_double)
    do q=1,nq
      if(.not.real_matrix_ok(effective_strength(:,:,q),.false.))return
      if(maxval(abs(effective_strength(:,:,q)))==0d0)cycle
      den=cmplx(resonance(q)**2-omega**2,-damping(q)*omega,c_double)
      matched=matched+effective_strength(:,:,q)/den
    enddo
    do j=1,3
      do i=1,3
        if(.not.ieee_is_finite(real(target(i,j),c_double)))return
        if(.not.ieee_is_finite(aimag(target(i,j))))return
        if(.not.ieee_is_finite(real(matched(i,j),c_double)))return
        if(.not.ieee_is_finite(aimag(matched(i,j))))return
      enddo
    enddo
    status=0
  end function subpixel_mix

  integer(c_int) function sp_subpixel_mix(nq,omega,fraction,normal,background,strength,resonance,damping,owner, &
      effective_background,effective_strength,target,matched)bind(C)result(status)
    integer(c_int),value::nq
    real(c_double),value::omega,fraction
    real(c_double),intent(in)::normal(3),background(3,3,2),strength(3,3,nq),resonance(nq),damping(nq)
    integer(c_int),intent(in)::owner(nq)
    real(c_double),intent(out)::effective_background(3,3),effective_strength(3,3,nq)
    complex(c_double_complex),intent(out)::target(3,3),matched(3,3)
    status=subpixel_mix(nq,omega,fraction,normal,background,strength,resonance,damping,owner, &
      effective_background,effective_strength,target,matched)
  end function

  pure function cross(a,b)result(c)
    real(c_double),intent(in)::a(3),b(3)
    real(c_double)::c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function

  integer function plane_fraction(normal,distance,spacing,fraction)result(status)
    real(c_double),intent(in)::normal(3),distance,spacing(3)
    real(c_double),intent(out)::fraction
    integer,parameter::faces(4,6)=reshape([1,4,3,2,5,6,7,8,1,5,8,4,2,3,7,6,1,2,6,5,4,8,7,3],[4,6])
    real(c_double)::vertex(3,8),polygon(3,8),clipped(3,8),cap(3,24),angles(24),center(3)
    real(c_double)::a(3),b(3),point(3),normal_unit(3),u(3),v(3),tmp(3),angle
    real(c_double)::h(3),d,scale,length,da,db,t,volume,tol
    integer::i,j,k,face,np,nclip,ncap,unique,axis
    logical::inside_a,inside_b,duplicate
    status=1;fraction=0d0
    if(.not.ieee_is_finite(distance))return
    do i=1,3
      if(.not.ieee_is_finite(normal(i)))return
      if(.not.ieee_is_finite(spacing(i)).or.spacing(i)<=0d0)return
    enddo
    scale=maxval(abs(normal));if(scale==0d0)return
    length=sqrt(sum((normal/scale)**2));normal_unit=(normal/scale)/length;d=(distance/scale)/length
    h=spacing/2d0
    if(d>=sum(abs(normal_unit)*h))then
      fraction=1d0;status=0;return
    endif
    if(d<=-sum(abs(normal_unit)*h))then
      fraction=0d0;status=0;return
    endif
    vertex(:,1)=[-h(1),-h(2),-h(3)];vertex(:,2)=[h(1),-h(2),-h(3)]
    vertex(:,3)=[h(1),h(2),-h(3)];vertex(:,4)=[-h(1),h(2),-h(3)]
    vertex(:,5)=[-h(1),-h(2),h(3)];vertex(:,6)=[h(1),-h(2),h(3)]
    vertex(:,7)=[h(1),h(2),h(3)];vertex(:,8)=[-h(1),h(2),h(3)]
    volume=0d0;ncap=0;tol=64d0*epsilon(1d0)*maxval(spacing)
    do face=1,6
      polygon(:,1:4)=vertex(:,faces(:,face));nclip=0
      do i=1,4
        j=mod(i,4)+1;a=polygon(:,i);b=polygon(:,j)
        da=sum(normal_unit*a)-d;db=sum(normal_unit*b)-d
        inside_a=da<=0d0;inside_b=db<=0d0
        if(inside_a)then
          nclip=nclip+1;clipped(:,nclip)=a
        endif
        if(inside_a.neqv.inside_b)then
          t=da/(da-db);point=a+t*(b-a)
          nclip=nclip+1;clipped(:,nclip)=point
          duplicate=.false.
          do k=1,ncap
            if(maxval(abs(cap(:,k)-point))<=tol)duplicate=.true.
          enddo
          if(.not.duplicate)then
            ncap=ncap+1;cap(:,ncap)=point
          endif
        endif
      enddo
      do i=2,nclip-1
        volume=volume+sum(clipped(:,1)*cross(clipped(:,i),clipped(:,i+1)))/6d0
      enddo
    enddo
    if(ncap>=3)then
      center=sum(cap(:,1:ncap),dim=2)/real(ncap,c_double)
      axis=minloc(abs(normal_unit),dim=1);a=0d0;a(axis)=1d0
      u=cross(a,normal_unit);u=u/sqrt(sum(u*u));v=cross(normal_unit,u)
      do i=1,ncap
        tmp=cap(:,i)-center;angles(i)=atan2(sum(tmp*v),sum(tmp*u))
      enddo
      do i=2,ncap
        point=cap(:,i);angle=angles(i);j=i-1
        do while(j>=1)
          if(angles(j)<=angle)exit
          cap(:,j+1)=cap(:,j);angles(j+1)=angles(j);j=j-1
        enddo
        cap(:,j+1)=point;angles(j+1)=angle
      enddo
      do i=2,ncap-1
        volume=volume+sum(cap(:,1)*cross(cap(:,i),cap(:,i+1)))/6d0
      enddo
    endif
    fraction=volume/product(spacing)
    if(.not.ieee_is_finite(fraction))return
    if(fraction< -1d-12.or.fraction>1d0+1d-12)return
    fraction=max(0d0,min(1d0,fraction));status=0
  end function plane_fraction

  integer(c_int) function sp_subpixel_plane_fraction(normal,distance,spacing,fraction)bind(C)result(status)
    real(c_double),intent(in)::normal(3),spacing(3)
    real(c_double),value::distance
    real(c_double),intent(out)::fraction
    status=plane_fraction(normal,distance,spacing,fraction)
  end function

  subroutine classify_shape(kind,size,center,point,spacing,classification,normal,distance,corner)
    integer,intent(in)::kind
    real(c_double),intent(in)::size(3),center(3),point(3),spacing(3)
    integer,intent(out)::classification
    real(c_double),intent(out)::normal(3),distance
    logical,intent(out)::corner
    real(c_double)::lo(3),hi(3),delta(3),radius(3),minsum,maxsum,g(3),length
    integer::a,hits
    delta=point-center;radius=size/2d0;lo=delta-spacing/2d0;hi=delta+spacing/2d0
    normal=0d0;distance=0d0;corner=.false.;classification=0
    if(kind==1)then
      if(any(lo>=radius).or.any(hi<=-radius))then
        classification=-1;return
      endif
      if(all(lo>=-radius).and.all(hi<=radius))then
        classification=1;return
      endif
      hits=0
      do a=1,3
        if(lo(a)<-radius(a).and.hi(a)>-radius(a))then
          hits=hits+1;normal=0d0;normal(a)=-1d0;distance=delta(a)+radius(a)
        endif
        if(lo(a)<radius(a).and.hi(a)>radius(a))then
          hits=hits+1;normal=0d0;normal(a)=1d0;distance=radius(a)-delta(a)
        endif
      enddo
      corner=hits/=1
    else
      minsum=sum((max(0d0,abs(delta)-spacing/2d0)/radius)**2)
      maxsum=sum(((abs(delta)+spacing/2d0)/radius)**2)
      if(minsum>=1d0)then
        classification=-1;return
      endif
      if(maxsum<=1d0)then
        classification=1;return
      endif
      g=2d0*delta/radius**2;length=sqrt(sum(g*g))
      if(length==0d0)then
        corner=.true.;return
      endif
      normal=g/length;distance=(1d0-sum((delta/radius)**2))/length
    endif
  end subroutine

  subroutine shape_cut(ns,kind,size,center,material,point,spacing,pair,fraction,normal,flag)
    integer,intent(in)::ns,kind(ns),material(ns)
    real(c_double),intent(in)::size(3,ns),center(3,ns),point(3),spacing(3)
    integer,intent(out)::pair(2),flag
    real(c_double),intent(out)::fraction,normal(3)
    integer::s,classification,partial,point_material,status
    real(c_double)::n(3),distance,dummy,delta(3),cut_distance
    logical::corner
    point_material=0;pair=0;fraction=1d0;normal=[1d0,0d0,0d0];flag=0;partial=0;cut_distance=0d0
    do s=1,ns
      delta=2d0*(point-center(:,s))/size(:,s)
      if(kind(s)==1)then
        if(all(abs(delta)<=1d0))point_material=material(s)
      else
        if(sum(delta*delta)<=1d0)point_material=material(s)
      endif
    enddo
    do s=ns,1,-1
      call classify_shape(kind(s),size(:,s),center(:,s),point,spacing,classification,n,distance,corner)
      if(classification==-1)cycle
      if(classification==1)then
        if(partial==0)then
          pair=material(s);return
        endif
        pair(2)=material(s);exit
      endif
      if(partial/=0.and..not.corner)then
        if(kind(s)==1.and.kind(partial)==1)then
          if(maxval(abs(n+normal))<64d0*epsilon(1d0).and. &
             abs(distance+cut_distance)<64d0*epsilon(1d0)*maxval(spacing))then
            pair(2)=material(s);exit
          endif
          if(maxval(abs(n-normal))<64d0*epsilon(1d0).and. &
             abs(distance-cut_distance)<64d0*epsilon(1d0)*maxval(spacing))cycle
        endif
      endif
      if(partial/=0.or.corner)then
        pair=point_material;flag=2;return
      endif
      partial=s;pair(1)=material(s);normal=n;cut_distance=distance
    enddo
    if(partial==0)then
      pair=point_material;return
    endif
    if(pair(1)==pair(2))return
    status=plane_fraction(normal,cut_distance,spacing,fraction)
    if(status/=0)then
      pair=point_material;flag=3;return
    endif
    flag=1
  end subroutine shape_cut

  integer(c_int) function sp_subpixel_shape_cut(ns,kind,size,center,material,point,spacing,pair,fraction,normal,flag) &
      bind(C)result(status)
    integer(c_int),value::ns
    integer(c_int),intent(in)::kind(ns),material(ns)
    real(c_double),intent(in)::size(3,ns),center(3,ns),point(3),spacing(3)
    integer(c_int),intent(out)::pair(2),flag
    real(c_double),intent(out)::fraction,normal(3)
    integer::s,a
    status=1
    if(ns<0)return
    do a=1,3
      if(.not.ieee_is_finite(point(a)))return
      if(.not.ieee_is_finite(spacing(a)).or.spacing(a)<=0d0)return
    enddo
    do s=1,ns
      if(kind(s)<1.or.kind(s)>2.or.material(s)<1)return
      do a=1,3
        if(.not.ieee_is_finite(size(a,s)).or.size(a,s)<=0d0)return
        if(.not.ieee_is_finite(center(a,s)))return
      enddo
    enddo
    call shape_cut(ns,kind,size,center,material,point,spacing,pair,fraction,normal,flag)
    status=0
  end function
  subroutine read_subpixel_profile(filename,profile,status)
    character(*),intent(in)::filename
    type(s_subpixel_profile),intent(out)::profile
    integer,intent(out)::status
    integer::unit,ios,version,nactual,m,id,np,q,k,total,a,b
    character(64)::magic
    status=1
    open(newunit=unit,file=filename,status='old',action='read',iostat=ios)
    if(ios/=0)return
    read(unit,*,iostat=ios)magic,version
    if(ios==0)then
      if(magic/='SALMONPY_SUBPIXEL_REFERENCE'.or.version/=2)ios=1
    endif
    if(ios==0)read(unit,*,iostat=ios)profile%omega,profile%nm,profile%ns,nactual
    if(ios==0)then
      if(profile%nm<0.or.profile%nm>200.or.profile%ns<0.or.profile%ns>1000)ios=1
      if(nactual<0.or.nactual>20000)ios=1
      if(.not.ieee_is_finite(profile%omega).or.profile%omega<=0d0)ios=1
    endif
    if(ios==0)read(unit,*,iostat=ios)profile%rtol,profile%maxiter
    if(ios==0)then
      if(.not.ieee_is_finite(profile%rtol).or.profile%rtol<=0d0.or.profile%rtol>=1d0)ios=1
      if(profile%maxiter<1)ios=1
    endif
    if(ios/=0)then
      close(unit);return
    endif
    profile%nq=max(1,nactual)
    call allocate_subpixel_profile(profile)
    profile%owner=0;profile%strength=0d0;profile%resonance=0d0;profile%damping=0d0
    total=0
    do m=0,profile%nm
      read(unit,*,iostat=ios)id,np,profile%ld_phase(m)
      if(ios/=0)exit
      if(id/=m.or.np<0.or.np>100.or.total+np>nactual)then
        ios=1;exit
      endif
      if(profile%ld_phase(m)<0.or.profile%ld_phase(m)>1)then
        ios=1;exit
      endif
      if((m==0.or.np>0).and.profile%ld_phase(m)/=merge(0,1,m==0))then
        ios=1;exit
      endif
      read(unit,*,iostat=ios)profile%background(:,:,m)
      if(ios/=0)exit
      if(.not.real_matrix_ok(profile%background(:,:,m),.true.))then
        ios=1;exit
      endif
      do k=1,np
        q=total+k
        read(unit,*,iostat=ios)profile%resonance(q),profile%damping(q)
        if(ios/=0)exit
        read(unit,*,iostat=ios)profile%strength(:,:,q)
        if(ios/=0)exit
        if(.not.ieee_is_finite(profile%resonance(q)).or.profile%resonance(q)<0d0)ios=1
        if(.not.ieee_is_finite(profile%damping(q)).or.profile%damping(q)<0d0)ios=1
        if(.not.real_matrix_ok(profile%strength(:,:,q),.false.))ios=1
        if(ios/=0)exit
        profile%owner(q)=m
      enddo
      if(ios/=0)exit
      total=total+np
    enddo
    if(total/=nactual)ios=1
    if(ios==0)then
      do k=1,profile%ns
        read(unit,*,iostat=ios)profile%kind(k),profile%material(k)
        if(ios/=0)exit
        read(unit,*,iostat=ios)profile%shape_size(:,k),profile%center(:,k)
        if(ios/=0)exit
        if(profile%kind(k)<1.or.profile%kind(k)>2)ios=1
        if(profile%material(k)<1.or.profile%material(k)>profile%nm)ios=1
        do a=1,3
          if(.not.ieee_is_finite(profile%shape_size(a,k)).or.profile%shape_size(a,k)<=0d0)ios=1
          if(.not.ieee_is_finite(profile%center(a,k)))ios=1
        enddo
        if(ios/=0)exit
      enddo
    endif
    close(unit)
    if(ios==0)status=0
  end subroutine

  subroutine allocate_subpixel_profile(profile)
    type(s_subpixel_profile),intent(inout)::profile
    allocate(profile%background(3,3,0:profile%nm),profile%strength(3,3,profile%nq), &
      profile%resonance(profile%nq),profile%damping(profile%nq),profile%owner(profile%nq), &
      profile%shape_size(3,profile%ns),profile%center(3,profile%ns), &
      profile%kind(profile%ns),profile%material(profile%ns),profile%ld_phase(0:profile%nm))
  end subroutine

  subroutine subpixel_bulk(profile,material,background,strength)
    type(s_subpixel_profile),intent(in)::profile
    integer,intent(in)::material
    real(c_double),intent(out)::background(3,3),strength(3,3,profile%nq)
    integer::q
    background=profile%background(:,:,material);strength=0d0
    do q=1,profile%nq
      if(profile%owner(q)==material)strength(:,:,q)=profile%strength(:,:,q)
    enddo
  end subroutine

  subroutine subpixel_at(profile,point,spacing,background,strength,flag,error)
    type(s_subpixel_profile),intent(in)::profile
    real(c_double),intent(in)::point(3),spacing(3)
    real(c_double),intent(out)::background(3,3),strength(3,3,profile%nq),error
    integer,intent(out)::flag
    integer::pair(2),owners(profile%nq),q,status,material
    real(c_double)::fraction,normal(3),input_background(3,3,2),input_strength(3,3,profile%nq)
    complex(c_double_complex)::target(3,3),matched(3,3)
    call shape_cut(profile%ns,profile%kind,profile%shape_size,profile%center,profile%material, &
      point,spacing,pair,fraction,normal,flag)
    error=0d0
    if(flag/=1.or.fraction==0d0.or.fraction==1d0)then
      material=pair(1)
      if(flag==1.and.fraction==0d0)material=pair(2)
      call subpixel_bulk(profile,material,background,strength)
      return
    endif
    if(profile%ld_phase(pair(1))==0.and.profile%ld_phase(pair(2))==0)then
      flag=5;material=pair(merge(1,2,fraction>=.5d0))
      call subpixel_bulk(profile,material,background,strength)
      return
    endif
    input_background(:,:,1)=profile%background(:,:,pair(1))
    input_background(:,:,2)=profile%background(:,:,pair(2))
    input_strength=0d0;owners=1
    do q=1,profile%nq
      if(profile%owner(q)==pair(1))then
        input_strength(:,:,q)=profile%strength(:,:,q)
      else if(profile%owner(q)==pair(2))then
        input_strength(:,:,q)=profile%strength(:,:,q);owners(q)=2
      endif
    enddo
    status=subpixel_mix(profile%nq,profile%omega,fraction,normal,input_background,input_strength, &
      profile%resonance,profile%damping,owners,background,strength,target,matched)
    if(status/=0)then
      flag=4;material=pair(merge(1,2,fraction>=.5d0))
      call subpixel_bulk(profile,material,background,strength)
    else
      error=maxval(abs(target-matched))
    endif
  end subroutine

end module subpixel_smoothing
