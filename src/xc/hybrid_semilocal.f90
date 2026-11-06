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
      call xc_func_end(func);call xc_func_free(func)
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
