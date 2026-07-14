module dg_wpw_fixed_operator
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  type,public :: s_dg_wpw_fixed_operator
    private
    integer :: n_wannier=0,n_pw=0,n_basis=0
    complex(8),allocatable :: s(:,:),h_fixed(:,:)
  end type
  public :: initialize_dg_wpw_fixed_operator,compose_dg_wpw_hamiltonian
  public :: dg_wpw_fixed_dimensions,copy_dg_wpw_metric,copy_dg_wpw_fixed_hamiltonian
contains
  subroutine initialize_dg_wpw_fixed_operator(op,sww,swp,spp,hww,hwp,hpp, &
      hfaceww,hfacewp,nw,np,info)
    type(s_dg_wpw_fixed_operator),intent(inout)::op
    integer,intent(in)::nw,np
    complex(8),intent(in)::sww(nw,nw),swp(nw,np),spp(np,np)
    complex(8),intent(in)::hww(nw,nw),hwp(nw,np),hpp(np,np)
    complex(8),intent(in)::hfaceww(nw,nw),hfacewp(nw,np)
    integer,intent(out)::info
    integer::n
    info=0;n=nw+np
    if(nw<1.or.np<0)then;info=1;return;endif
    if(.not.finite_matrix(sww).or..not.finite_matrix(swp).or..not.finite_matrix(spp).or. &
       .not.finite_matrix(hww).or..not.finite_matrix(hwp).or..not.finite_matrix(hpp).or. &
       .not.finite_matrix(hfaceww).or..not.finite_matrix(hfacewp))then
      info=2;return
    endif
    if(allocated(op%s))deallocate(op%s,op%h_fixed)
    allocate(op%s(n,n),op%h_fixed(n,n));op%s=0;op%h_fixed=0
    op%s(1:nw,1:nw)=sww
    op%h_fixed(1:nw,1:nw)=hww+hfaceww
    if(np>0)then
      op%s(1:nw,nw+1:n)=swp;op%s(nw+1:n,1:nw)=conjg(transpose(swp));op%s(nw+1:n,nw+1:n)=spp
      op%h_fixed(1:nw,nw+1:n)=hwp+hfacewp
      op%h_fixed(nw+1:n,1:nw)=conjg(transpose(hwp+hfacewp))
      ! P columns are global H1 functions: [[P]]=0, so PP face terms vanish.
      op%h_fixed(nw+1:n,nw+1:n)=hpp
    endif
    if(.not.hermitian(op%s).or..not.hermitian(op%h_fixed))then
      deallocate(op%s,op%h_fixed);info=3;return
    endif
    op%n_wannier=nw;op%n_pw=np;op%n_basis=n
  end subroutine

  subroutine dg_wpw_fixed_dimensions(op,nw,np,n)
    type(s_dg_wpw_fixed_operator),intent(in)::op
    integer,intent(out)::nw,np,n
    nw=op%n_wannier;np=op%n_pw;n=op%n_basis
  end subroutine
  subroutine copy_dg_wpw_metric(op,s,info)
    type(s_dg_wpw_fixed_operator),intent(in)::op;complex(8),intent(out)::s(:,:);integer,intent(out)::info
    info=0;s=0;if(.not.allocated(op%s).or.any(shape(s)/=[op%n_basis,op%n_basis]))then;info=1;return;endif;s=op%s
  end subroutine
  subroutine copy_dg_wpw_fixed_hamiltonian(op,h,info)
    type(s_dg_wpw_fixed_operator),intent(in)::op;complex(8),intent(out)::h(:,:);integer,intent(out)::info
    info=0;h=0;if(.not.allocated(op%h_fixed).or.any(shape(h)/=[op%n_basis,op%n_basis]))then;info=1;return;endif;h=op%h_fixed
  end subroutine

  subroutine compose_dg_wpw_hamiltonian(op,v_local,h,info)
    type(s_dg_wpw_fixed_operator),intent(in)::op
    complex(8),intent(in)::v_local(:,:)
    complex(8),intent(out)::h(:,:)
    integer,intent(out)::info
    info=0;h=0
    if(.not.allocated(op%h_fixed).or.size(v_local,1)/=op%n_basis.or.size(v_local,2)/=op%n_basis.or. &
       size(h,1)/=op%n_basis.or.size(h,2)/=op%n_basis)then;info=1;return;endif
    if(.not.finite_matrix(v_local).or..not.hermitian(v_local))then;info=2;return;endif
    h=op%h_fixed+v_local
  end subroutine

  logical function finite_matrix(a)
    complex(8),intent(in)::a(:,:)
    finite_matrix=all(ieee_is_finite(real(a,8))).and.all(ieee_is_finite(aimag(a)))
  end function
  logical function hermitian(a)
    complex(8),intent(in)::a(:,:)
    hermitian=size(a,1)==size(a,2)
    if(hermitian)hermitian=maxval(abs(a-conjg(transpose(a))))<=1d-11*max(1d0,maxval(abs(a)))
  end function
end module
