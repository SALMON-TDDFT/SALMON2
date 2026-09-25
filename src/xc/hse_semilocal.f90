! Libxc C ABI avoids compiler-specific libxcf03 module dependencies.
module hse_semilocal
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: hse_semilocal_evaluate
  interface
    function xc_func_alloc() bind(C) result(p)
      import c_ptr
      type(c_ptr) :: p
    end function
    function xc_func_init(p,id,spin) bind(C) result(ierr)
      import c_ptr,c_int
      type(c_ptr),value :: p
      integer(c_int),value :: id,spin
      integer(c_int) :: ierr
    end function
    subroutine xc_func_end(p) bind(C)
      import c_ptr
      type(c_ptr),value :: p
    end subroutine
    subroutine xc_func_free(p) bind(C)
      import c_ptr
      type(c_ptr),value :: p
    end subroutine
    subroutine xc_func_set_ext_params(p,values) bind(C)
      import c_ptr,c_double
      type(c_ptr),value :: p
      real(c_double),intent(in) :: values(*)
    end subroutine
    subroutine xc_hyb_cam_coef(p,omega,alpha,beta) bind(C)
      import c_ptr,c_double
      type(c_ptr),value :: p
      real(c_double) :: omega,alpha,beta
    end subroutine
    subroutine xc_gga_exc_vxc(p,n,rho,sigma,eps,vrho,vsigma) bind(C)
      import c_ptr,c_double,c_size_t
      type(c_ptr),value :: p
      integer(c_size_t),value :: n
      real(c_double) :: rho(*),sigma(*),eps(*),vrho(*),vsigma(*)
    end subroutine
  end interface
contains
  subroutine hse_semilocal_evaluate(rho,sigma,eps,vrho,vsigma,ierr,screening)
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
    if(.not.all(ieee_is_finite(rho)).or..not.all(ieee_is_finite(sigma)))return
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
      if(.not.all(ieee_is_finite(eps)).or..not.all(ieee_is_finite(vrho)).or. &
         .not.all(ieee_is_finite(vsigma)))ierr=1
    endif
    call xc_func_end(func);call xc_func_free(func)
  end subroutine
end module
