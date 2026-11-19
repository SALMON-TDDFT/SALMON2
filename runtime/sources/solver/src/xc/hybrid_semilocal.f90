! Libxc C ABI avoids compiler-specific libxcf03 module dependencies.
module hybrid_semilocal
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: hse_semilocal_evaluate,pbeh_semilocal_evaluate,hybrid_density_evaluate
  interface
    function xc_func_alloc() bind(C) result(p)
      import c_ptr
      implicit none
      type(c_ptr) :: p
    end function
    function xc_func_init(p,id,spin) bind(C) result(ierr)
      import c_ptr,c_int
      implicit none
      type(c_ptr),value :: p
      integer(c_int),value :: id,spin
      integer(c_int) :: ierr
    end function
    subroutine xc_func_end(p) bind(C)
      import c_ptr
      implicit none
      type(c_ptr),value :: p
    end subroutine
    subroutine xc_func_free(p) bind(C)
      import c_ptr
      implicit none
      type(c_ptr),value :: p
    end subroutine
    subroutine xc_func_set_ext_params(p,values) bind(C)
      import c_ptr,c_double
      implicit none
      type(c_ptr),value :: p
      real(c_double),intent(in) :: values(*)
    end subroutine
    subroutine xc_hyb_cam_coef(p,omega,alpha,beta) bind(C)
      import c_ptr,c_double
      implicit none
      type(c_ptr),value :: p
      real(c_double) :: omega,alpha,beta
    end subroutine
    subroutine xc_gga_exc_vxc(p,n,rho,sigma,eps,vrho,vsigma) bind(C)
      import c_ptr,c_double,c_size_t
      implicit none
      type(c_ptr),value :: p
      integer(c_size_t),value :: n
      real(c_double) :: rho(*),sigma(*),eps(*),vrho(*),vsigma(*)
    end subroutine
  end interface
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_1d
  interface salmon_all_finite
    module procedure finite_real_1d
  end interface
contains
  subroutine pbeh_semilocal_evaluate(rho,sigma,eps,vrho,vsigma,ierr,exchange_fraction)
    implicit none
    real(c_double),intent(in) :: rho(:),sigma(:)
    real(c_double),intent(out) :: eps(:),vrho(:),vsigma(:)
    integer,intent(out) :: ierr
    real(c_double),intent(in),optional :: exchange_fraction
    type(c_ptr) :: func
    real(c_double) :: e(size(rho)),v(size(rho)),s(size(rho)),weight
    integer(c_int) :: ids(2)=[101_c_int,130_c_int]
    integer :: n,j
    ierr=1;n=size(rho)
    if(size(sigma)/=n.or.size(eps)/=n.or.size(vrho)/=n.or.size(vsigma)/=n)return
    if(any(rho<0d0).or.any(sigma<0d0))return
    if(.not.salmon_all_finite(rho).or..not.salmon_all_finite(sigma))return
    eps=0d0;vrho=0d0;vsigma=0d0
    do j=1,2
      ierr=1
      func=xc_func_alloc();if(.not.c_associated(func))return
      ierr=xc_func_init(func,ids(j),1_c_int)
      if(ierr/=0)then
        call xc_func_free(func);return
      endif
      call xc_gga_exc_vxc(func,int(n,c_size_t),rho,sigma,e,v,s)
      call xc_func_end(func);call xc_func_free(func)
      weight=1d0
      if(j==1)then
        weight=.6d0
        if(present(exchange_fraction))weight=1d0-exchange_fraction
      endif
      eps=eps+weight*e;vrho=vrho+weight*v;vsigma=vsigma+weight*s
    enddo
    where(rho==0d0)
      eps=0d0;vrho=0d0;vsigma=0d0
    endwhere
    ierr=0
    if(.not.salmon_all_finite(eps).or..not.salmon_all_finite(vrho).or. &
       .not.salmon_all_finite(vsigma))ierr=1
  end subroutine
  subroutine hse_semilocal_evaluate(rho,sigma,eps,vrho,vsigma,ierr,screening)
    implicit none
    real(c_double),intent(in) :: rho(:),sigma(:)
    real(c_double),intent(out) :: eps(:),vrho(:),vsigma(:)
    integer,intent(out) :: ierr
    real(c_double),intent(in),optional :: screening
    type(c_ptr) :: func
    integer :: n
    real(c_double) :: omega,alpha,beta
    real(c_double) :: requested
    ierr=1;n=size(rho)
    requested=.11d0
    if(present(screening))requested=screening
    if(.not.ieee_is_finite(requested).or.requested<=0d0)return
    if(size(sigma)/=n.or.size(eps)/=n.or.size(vrho)/=n.or.size(vsigma)/=n)return
    if(any(rho<0).or.any(sigma<0))return
    if(.not.salmon_all_finite(rho).or..not.salmon_all_finite(sigma))return
    func=xc_func_alloc();if(.not.c_associated(func))return
    ierr=xc_func_init(func,428_c_int,1_c_int)
    if(ierr/=0)then
      call xc_func_free(func);return
    endif
    ! Set beta, omega_HF and omega_PBE together: old Libxc named setters
    ! restore the other parameters to their defaults on each call.
    call xc_func_set_ext_params(func,[.25d0,requested,requested])
    call xc_hyb_cam_coef(func,omega,alpha,beta)
    if(abs(omega-requested)>1d-12.or.abs(alpha)>1d-12.or.abs(alpha+beta-.25d0)>1d-12)then
      ierr=1
    else
      call xc_gga_exc_vxc(func,int(n,c_size_t),rho,sigma,eps,vrho,vsigma)
      call stable_hse_small_gradient(1,n,rho,sigma,eps,vrho,vsigma,requested,ierr)
      where(rho==0d0)
        eps=0d0;vrho=0d0;vsigma=0d0
      endwhere
      if(.not.salmon_all_finite(eps).or..not.salmon_all_finite(vrho).or. &
         .not.salmon_all_finite(vsigma))ierr=1
    endif
    call xc_func_end(func);call xc_func_free(func)
  end subroutine

  ! Fixed-core NLCC: differentiate E_xc[n_val+n_core] with respect to n_val.
  ! The core is never a Fock source. gradient_derivative is +dE/d(grad n).
  subroutine hybrid_density_evaluate(rho,gradient,core,core_gradient,energy,potential, &
                                     gradient_derivative,ierr,screening,exchange_fraction)
    implicit none
    real(c_double),intent(in) :: rho(:,:),gradient(:,:,:),core(:),core_gradient(:,:)
    real(c_double),intent(out) :: energy(:),potential(:,:),gradient_derivative(:,:,:)
    integer,intent(out) :: ierr
    real(c_double),intent(in),optional :: screening,exchange_fraction
    real(c_double) :: r(size(rho,1),size(core)),g(size(rho,1),3,size(core))
    real(c_double) :: sigma(2*size(rho,1)-1,size(core)),ep(size(core)),v(size(rho,1),size(core))
    real(c_double) :: vs(2*size(rho,1)-1,size(core)),weight,omega
    type(c_ptr) :: func
    integer :: j,n,ns,part,id,axis
    ierr=1;n=size(core);ns=size(rho,1)
    if(ns/=1.and.ns/=2)return
    if(size(rho,2)/=n.or.any(shape(gradient)/=[ns,3,n]).or.any(shape(core_gradient)/=[3,n]))return
    if(size(energy)/=n.or.any(shape(potential)/=[ns,n]).or.any(shape(gradient_derivative)/=[ns,3,n]))return
    if(.not.salmon_all_finite(core).or.any(core<0d0))return
    do j=1,ns
      if(.not.salmon_all_finite(rho(j,:)).or.any(rho(j,:)<-1d-12))return
      r(j,:)=max(rho(j,:),0d0)+core/real(ns,c_double)
      g(j,:,:)=gradient(j,:,:)+core_gradient/real(ns,c_double)
    enddo
    sigma(1,:)=sum(g(1,:,:)**2,dim=1)
    if(ns==2)then
      sigma(2,:)=sum(g(1,:,:)*g(2,:,:),dim=1)
      sigma(3,:)=sum(g(2,:,:)**2,dim=1)
    endif
    do j=1,size(sigma,1)
      if(.not.salmon_all_finite(sigma(j,:)))return
    enddo
    omega=.11d0
    if(present(screening))omega=screening
    if(.not.ieee_is_finite(omega).or.omega<=0d0)return
    energy=0d0;potential=0d0;gradient_derivative=0d0
    do part=1,merge(2,1,present(exchange_fraction))
      id=428;weight=1d0
      if(present(exchange_fraction))then
        id=130
        if(part==1)then
          id=101;weight=1d0-exchange_fraction
        endif
      endif
      func=xc_func_alloc();if(.not.c_associated(func))return
      ierr=xc_func_init(func,int(id,c_int),int(ns,c_int))
      if(ierr/=0)then
        call xc_func_free(func);return
      endif
      if(id==428)call xc_func_set_ext_params(func,[.25d0,omega,omega])
      call xc_gga_exc_vxc(func,int(n,c_size_t),r,sigma,ep,v,vs)
      if(id==428)call stable_hse_small_gradient(ns,n,r,sigma,ep,v,vs,omega,ierr)
      call xc_func_end(func);call xc_func_free(func)
      if(ierr/=0)return
      energy=energy+weight*sum(r,dim=1)*ep
      potential=potential+weight*v
      do j=1,3
        gradient_derivative(1,j,:)=gradient_derivative(1,j,:)+weight*2d0*vs(1,:)*g(1,j,:)
        if(ns==2)then
          gradient_derivative(1,j,:)=gradient_derivative(1,j,:)+weight*vs(2,:)*g(2,j,:)
          gradient_derivative(2,j,:)=gradient_derivative(2,j,:)+ &
            weight*(2d0*vs(3,:)*g(2,j,:)+vs(2,:)*g(1,j,:))
        endif
      enddo
    enddo
    ierr=1
    if(.not.salmon_all_finite(energy))return
    do j=1,ns
      if(.not.salmon_all_finite(potential(j,:)))return
      do axis=1,3
        if(.not.salmon_all_finite(gradient_derivative(j,axis,:)))return
      enddo
    enddo
    ierr=0
  end subroutine hybrid_density_evaluate


  ! Stable evaluation of the SAME WPBEh exchange formula used by Libxc HSE06.
  ! Formula provenance: Libxc 5.2.3 maple/gga_exc/gga_x_wpbeh.mpl,
  ! Copyright (C) 2017 M.A.L. Marques, MPL-2.0. See runtime/NOTICE.md.
  ! For omega=0, combine exp(x)*E1(x)+log(y/(D+y)) BEFORE differentiating:
  ! x=c*y, y=u*H(sqrt(u)), u=s*s. The cancelling 1/y terms disappear.
  ! The omitted E1-series terms are O(x*x*abs(log(x))); at u<=1e-8,
  ! x<3e-18, far below double precision. This is not a sigma floor.
  subroutine wpbeh_zero_small(r,sigma,e,v,vs)
    real(c_double),intent(in) :: r,sigma
    real(c_double),intent(out) :: e,v,vs
    real(c_double),parameter :: a=1.0161144d0,b=-.37170836d0,c=-.077215461d0
    real(c_double),parameter :: d=.57786348d0,ee=-.051955731d0,gamma=.5772156649015328606d0
    real(c_double) :: pi,k,u,h,dh,den,num,y,dy,z,x,l,t1,dt1,t2,dt2,f,df,gg,dg,enh,de
    pi=acos(-1d0);k=(3d0*pi*pi*r)**(1d0/3d0);u=sigma/(4d0*k*k*r*r)
    den=1d0+.187440d0*u*u+.00120824d0*u**2.5d0+.0347188d0*u**3
    num=.00979681d0*u+.0410834d0*u*u
    h=num/den
    dh=(.00979681d0+2d0*.0410834d0*u-h*(2d0*.187440d0*u+ &
        2.5d0*.00120824d0*u**1.5d0+3d0*.0347188d0*u*u))/den
    y=u*h;dy=h+u*dh;z=d+y;x=9d0*y/(4d0*a)
    t1=.5d0*a*(-gamma-log(9d0/(4d0*a))-log(z));dt1=0d0
    if(x>tiny(1d0))then
      l=-gamma-log(x)
      t1=t1+.5d0*a*x*(l+1d0)
      dt1=.5d0*a*(9d0/(4d0*a)*(l+x*(l+1d0))-1d0/z)*dy
    endif
    f=6.4753871d0*h+.47965830d0;df=6.4753871d0*dh
    gg=-.02628417880d0-.07117647788d0*u+.08534541323d0*u*u
    dg=-.07117647788d0+2d0*.08534541323d0*u
    num=z*z*b+z*c+2d0*ee+z*u*c*f+2d0*u*gg
    t2=num/(2d0*z**3)
    dt2=(2d0*z*dy*b+dy*c+c*((dy*u+z)*f+z*u*df)+ &
         2d0*gg+2d0*u*dg-3d0*num*dy/z)/(2d0*z**3)
    enh=-8d0/9d0*(t1+t2);de=-8d0/9d0*(dt1+dt2)
    z=-.75d0*(3d0/pi)**(1d0/3d0)
    e=z*r**(1d0/3d0)*enh
    v=z*r**(1d0/3d0)*(4d0/3d0*enh-8d0/3d0*u*de)
    vs=z*r**(1d0/3d0)*de/(4d0*k*k*r)
  end subroutine

  pure real(c_double) function wpbeh_screened_quadratic(w) result(f2)
    real(c_double),intent(in) :: w
    real(c_double),parameter :: d=.57786348d0,c=-.077215461d0,f=.47965830d0,g=-.02628417880d0
    real(c_double) :: q
    q=d+w*w
    f2=-8d0/9d0*((d*c*f+2d0*g)/(2d0*d**3) &
       -w*(6d0*q*c*f+15d0*g)/(8d0*d*q**2.5d0) &
       -w**3*(q*c*f+5d0*g)/(2d0*d*d*q**2.5d0)-w**5*g/(d**3*q**2.5d0))
  end function

  subroutine stable_hse_small_gradient(ns,n,r,sigma,e,v,vs,omega,ierr)
    integer,intent(in) :: ns,n
    integer,intent(out) :: ierr
    real(c_double),intent(in) :: r(ns,n),sigma(2*ns-1,n),omega
    real(c_double),intent(inout) :: e(n),v(ns,n),vs(2*ns-1,n)
    type(c_ptr) :: corr,screen,unscreen
    real(c_double) :: rr(1),ss(1),ex(1),vx(1),sx(1),es(1),vv(1),sv(1)
    real(c_double) :: ec(1),vc(ns,1),sc(2*ns-1,1),e0,v0,s0,k,u,density,scale,total,pi
    integer :: i,j,idx,status
    logical :: selected,ready_c,ready_s,ready_u
    ierr=0;ready_c=.false.;ready_s=.false.;ready_u=.false.
    corr=c_null_ptr;screen=c_null_ptr;unscreen=c_null_ptr
    pi=acos(-1d0);scale=real(ns,c_double)
    do i=1,n
      selected=.false.
      do j=1,ns
        density=scale*r(j,i)
        if(density<=2d-15)cycle ! Retain Libxc's low-density screening.
        k=(3d0*pi*pi*density)**(1d0/3d0)
        u=scale**2*max(0d0,sigma(2*j-1,i))/(4d0*k*k*density*density)
        if(u<=1d-8)selected=.true.
      enddo
      if(.not.selected)cycle
      if(.not.c_associated(corr))then
        corr=xc_func_alloc();screen=xc_func_alloc();unscreen=xc_func_alloc()
        ierr=1
        if(.not.c_associated(corr).or..not.c_associated(screen).or..not.c_associated(unscreen))goto 900
        status=xc_func_init(corr,130_c_int,int(ns,c_int));if(status/=0)goto 900
        ready_c=.true.
        status=xc_func_init(screen,524_c_int,1_c_int);if(status/=0)goto 900
        ready_s=.true.
        status=xc_func_init(unscreen,524_c_int,1_c_int);if(status/=0)goto 900
        ready_u=.true.;ierr=0
        call xc_func_set_ext_params(screen,[omega])
        call xc_func_set_ext_params(unscreen,[0d0])
      endif
      call xc_gga_exc_vxc(corr,1_c_size_t,r(:,i),sigma(:,i),ec,vc,sc)
      total=sum(r(:,i));e(i)=ec(1)*total;v(:,i)=vc(:,1);vs(:,i)=sc(:,1)
      do j=1,ns
        idx=2*j-1;rr(1)=scale*r(j,i);ss(1)=scale**2*max(0d0,sigma(idx,i))
        if(rr(1)<=2d-15)cycle
        k=(3d0*pi*pi*rr(1))**(1d0/3d0);u=ss(1)/(4d0*k*k*rr(1)**2)
        if(u<=1d-8)then
          call wpbeh_zero_small(rr(1),ss(1),e0,v0,s0)
          ex(1)=e0;vx(1)=v0;sx(1)=s0
        else
          call xc_gga_exc_vxc(unscreen,1_c_size_t,rr,ss,ex,vx,sx)
        endif
        call xc_gga_exc_vxc(screen,1_c_size_t,rr,ss,es,vv,sv)
        ! At the library's exact-zero s clamp, its derivative is zero rather
        ! than the finite quadratic limit. Higher terms are below roundoff.
        if(u<=epsilon(1d0))sv(1)=-.75d0*(3d0/pi)**(1d0/3d0)* &
          wpbeh_screened_quadratic(omega/k)/(4d0*k*k*rr(1)**(2d0/3d0))
        e(i)=e(i)+r(j,i)*(ex(1)-.25d0*es(1))
        v(j,i)=v(j,i)+vx(1)-.25d0*vv(1)
        vs(idx,i)=vs(idx,i)+scale*(sx(1)-.25d0*sv(1))
      enddo
      if(total>0d0)e(i)=e(i)/total
    enddo
900 continue
    if(c_associated(corr))then
      if(ready_c)call xc_func_end(corr)
      call xc_func_free(corr)
    endif
    if(c_associated(screen))then
      if(ready_s)call xc_func_end(screen)
      call xc_func_free(screen)
    endif
    if(c_associated(unscreen))then
      if(ready_u)call xc_func_end(unscreen)
      call xc_func_free(unscreen)
    endif
  end subroutine

  pure logical function finite_real_1d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:)
    real(8) :: value
    integer :: i
    finite=.false.
    do i=1,size(values,1)
      value=values(i)
      if(.not.ieee_is_finite(value))return
    enddo
    finite=.true.
  end function

end module
