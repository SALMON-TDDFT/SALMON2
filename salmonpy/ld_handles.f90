! Copyright 2026 SALMON developers. Apache-2.0.
! Private persistent LD adapter. Shared kernel arithmetic remains unchanged.
module ld_handles
  use iso_c_binding, only: c_int,c_int64_t,c_double
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use lorentz_drude, only: s_ld_state,ld_coefficients,ld_store_old,ld_update
  use lorentz_drude_tensor,only:sp_ld_tensor_validate
  implicit none
  private
  public :: sp_ld_create_tensor,sp_ld_create,sp_ld_advance,sp_ld_export,sp_ld_set,sp_ld_destroy
  type entry
    integer(c_int64_t) :: id=0
    logical :: tensor=.false.
    real(c_double),allocatable :: drive(:,:,:)
    integer :: n=0
    real(c_double) :: dt=0,field(3)=0
    type(s_ld_state) :: ld
    real(c_double),allocatable :: saved(:,:,:)
    type(entry),pointer :: next=>null()
  end type
  type(entry),pointer,save :: first=>null()
  integer(c_int64_t),save :: next_id=1
contains
  function lookup(id) result(s)
    integer(c_int64_t),intent(in) :: id
    type(entry),pointer :: s
    s=>first
    do while(associated(s))
      if(s%id==id)return
      s=>s%next
    end do
  end function

  integer(c_int) function sp_ld_create(n,dt,wp,f,w,g,e,id) bind(C) result(status)
    integer(c_int),value :: n
    real(c_double),value :: dt,wp
    real(c_double),intent(in) :: f(n),w(n),g(n),e(3)
    integer(c_int64_t),intent(out) :: id
    type(entry),pointer :: s
    real(c_double),allocatable :: ff(:,:),ww(:,:),gg(:,:),c(:,:)
    real(c_double) :: c4(1),c5(1)
    integer :: q,k,err
    status=1;id=0
    if(n<1 .or. n>huge(n)/3)return
    if(.not.ieee_is_finite(dt) .or. .not.ieee_is_finite(wp))return
    if(dt<=0 .or. wp<0)return
    do k=1,3
      if(.not.ieee_is_finite(e(k)))return
    end do
    do q=1,n
      if(.not.ieee_is_finite(f(q)) .or. .not.ieee_is_finite(w(q)) .or. .not.ieee_is_finite(g(q)))return
      if(f(q)<0 .or. w(q)<0 .or. g(q)<0)return
    end do
    status=3
    if(next_id==huge(next_id))return
    allocate(s,stat=err)
    if(err/=0)return
    s%n=n;s%dt=dt;s%field=e
    s%ld%enabled=.true.;s%ld%num_ld=1;s%ld%max_pole_num_ld=n
    allocate(s%saved(3,n,2),ff(0:0,n),ww(0:0,n),gg(0:0,n),c(n,4),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    allocate(s%ld%rjx_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%rjx_ld=0
    allocate(s%ld%rjx_old_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%rjx_old_ld=0
    allocate(s%ld%rjy_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%rjy_ld=0
    allocate(s%ld%rjy_old_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%rjy_old_ld=0
    allocate(s%ld%rjz_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%rjz_ld=0
    allocate(s%ld%rjz_old_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%rjz_old_ld=0
    allocate(s%ld%px_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%px_ld=0
    allocate(s%ld%py_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%py_ld=0
    allocate(s%ld%pz_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%pz_ld=0
    allocate(s%ld%c1_jx_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%c1_jx_ld=0
    allocate(s%ld%c2_jx_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%c2_jx_ld=0
    allocate(s%ld%c3_jx_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%c3_jx_ld=0
    allocate(s%ld%c1_jy_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%c1_jy_ld=0
    allocate(s%ld%c2_jy_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%c2_jy_ld=0
    allocate(s%ld%c3_jy_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%c3_jy_ld=0
    allocate(s%ld%c1_jz_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%c1_jz_ld=0
    allocate(s%ld%c2_jz_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%c2_jz_ld=0
    allocate(s%ld%c3_jz_ld(1,1,1,n),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%c3_jz_ld=0
    allocate(s%ld%rjx_fdtd_ld(1,1,1),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%rjx_fdtd_ld=0
    allocate(s%ld%rjy_fdtd_ld(1,1,1),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%rjy_fdtd_ld=0
    allocate(s%ld%rjz_fdtd_ld(1,1,1),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%rjz_fdtd_ld=0
    allocate(s%ld%ex_old_ld(1,1,1),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%ex_old_ld=0
    allocate(s%ld%ey_old_ld(1,1,1),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%ey_old_ld=0
    allocate(s%ld%ez_old_ld(1,1,1),stat=err)
    if(err/=0)then
      deallocate(s);return
    end if
    s%ld%ez_old_ld=0
    ff(0,:)=f;ww(0,:)=w;gg(0,:)=g;c4=0;c5=0;c=0
    call ld_coefficients([0],[n],gg,ww,ff,[wp],dt,acos(-1.d0),c(:,1:1),c(:,2:2),c(:,3:3),c(:,4:4),c4,c5)
    status=4
    do q=1,n
      do k=1,4
        if(.not.ieee_is_finite(c(q,k)))then
          deallocate(s);return
        end if
      end do
      s%ld%c1_jx_ld(1,1,1,q)=c(q,2)
      s%ld%c2_jx_ld(1,1,1,q)=c(q,3)
      s%ld%c3_jx_ld(1,1,1,q)=c(q,4)
      s%ld%c1_jy_ld(1,1,1,q)=c(q,2)
      s%ld%c2_jy_ld(1,1,1,q)=c(q,3)
      s%ld%c3_jy_ld(1,1,1,q)=c(q,4)
      s%ld%c1_jz_ld(1,1,1,q)=c(q,2)
      s%ld%c2_jz_ld(1,1,1,q)=c(q,3)
      s%ld%c3_jz_ld(1,1,1,q)=c(q,4)
    end do
    s%id=next_id;next_id=next_id+1;s%next=>first;first=>s;id=s%id;status=0
  end function

  integer(c_int) function sp_ld_create_tensor(n,dt,wp,f,w,g,eps,e,id) bind(C) result(status)
    integer(c_int),value::n
    real(c_double),value::dt,wp
    real(c_double),intent(in)::f(3,3,n),w(n),g(n),eps(3,3),e(3)
    integer(c_int64_t),intent(out)::id
    real(c_double)::den,v
    real(c_double),allocatable::zero_f(:)
    integer::i,j,q,err,cleanup
    type(entry),pointer::s
    status=sp_ld_tensor_validate(n,wp,f,w,g,eps);id=0
    if(status/=0)return
    status=3
    allocate(zero_f(n),stat=err)
    if(err/=0)return
    zero_f=0d0
    status=sp_ld_create(n,dt,1d0,zero_f,w,g,e,id)
    if(status/=0)return
    s=>lookup(id)
    allocate(s%drive(3,3,n),stat=err)
    if(err/=0)then
      cleanup=sp_ld_destroy(id);id=0;status=3;return
    endif
    status=4
    do q=1,n
      den=1d0+g(q)*dt/2d0+(w(q)*dt)**2/4d0
      do j=1,3
        do i=1,3
          v=f(i,j,q)*wp**2*dt/(8d0*acos(-1d0)*den)
          if(.not.ieee_is_finite(v))then
            cleanup=sp_ld_destroy(id);id=0;return
          endif
          s%drive(i,j,q)=v
        enddo
      enddo
    enddo
    s%tensor=.true.;status=0
  end function

  subroutine tensor_update(s,old,e)
    type(entry),intent(inout)::s
    real(c_double),intent(in)::old(3),e(3)
    real(c_double)::p(3,s%n),j(3,s%n),v(3),c1,c3
    integer::q,i,k
    call transfer_poles(s%ld,s%n,p,j,.false.)
    do q=1,s%n
      c1=s%ld%c1_jx_ld(1,1,1,q);c3=s%ld%c3_jx_ld(1,1,1,q)
      do i=1,3
        v(i)=c1*j(i,q)-c3*p(i,q)
        do k=1,3
          v(i)=v(i)+s%drive(i,k,q)*(old(k)+e(k))
        enddo
      enddo
      do i=1,3
        p(i,q)=p(i,q)+.5d0*s%dt*(v(i)+j(i,q))
        j(i,q)=v(i)
      enddo
    enddo
    call transfer_poles(s%ld,s%n,p,j,.true.)
  end subroutine

  subroutine transfer_poles(ld,n,p,j,put)
    type(s_ld_state),intent(inout) :: ld
    integer,intent(in) :: n
    real(c_double),intent(inout) :: p(3,n),j(3,n)
    logical,intent(in) :: put
    integer :: q
    do q=1,n
      if(put)then
        ld%px_ld(1,1,1,q)=p(1,q);ld%rjx_ld(1,1,1,q)=j(1,q)
        ld%py_ld(1,1,1,q)=p(2,q);ld%rjy_ld(1,1,1,q)=j(2,q)
        ld%pz_ld(1,1,1,q)=p(3,q);ld%rjz_ld(1,1,1,q)=j(3,q)
      else
        p(1,q)=ld%px_ld(1,1,1,q);j(1,q)=ld%rjx_ld(1,1,1,q)
        p(2,q)=ld%py_ld(1,1,1,q);j(2,q)=ld%rjy_ld(1,1,1,q)
        p(3,q)=ld%pz_ld(1,1,1,q);j(3,q)=ld%rjz_ld(1,1,1,q)
      end if
    end do
  end subroutine

  subroutine totals(s,p,j)
    type(entry),intent(in) :: s
    real(c_double),intent(out) :: p(3),j(3)
    p(1)=sum(s%ld%px_ld(1,1,1,:));j(1)=sum(s%ld%rjx_ld(1,1,1,:))
    p(2)=sum(s%ld%py_ld(1,1,1,:));j(2)=sum(s%ld%rjy_ld(1,1,1,:))
    p(3)=sum(s%ld%pz_ld(1,1,1,:));j(3)=sum(s%ld%rjz_ld(1,1,1,:))
  end subroutine

  logical function finite_state(s)
    type(entry),intent(in) :: s
    integer :: q
    finite_state=.false.
    do q=1,s%n
      if(.not.ieee_is_finite(s%ld%px_ld(1,1,1,q)))return
      if(.not.ieee_is_finite(s%ld%rjx_ld(1,1,1,q)))return
      if(.not.ieee_is_finite(s%ld%py_ld(1,1,1,q)))return
      if(.not.ieee_is_finite(s%ld%rjy_ld(1,1,1,q)))return
      if(.not.ieee_is_finite(s%ld%pz_ld(1,1,1,q)))return
      if(.not.ieee_is_finite(s%ld%rjz_ld(1,1,1,q)))return
    end do
    finite_state=.true.
  end function

  integer(c_int) function sp_ld_advance(id,nt,e,p,j) bind(C) result(status)
    integer(c_int64_t),value :: id
    integer(c_int),value :: nt
    real(c_double),intent(in) :: e(3,nt)
    real(c_double),intent(out) :: p(3,nt+1),j(3,nt+1)
    type(entry),pointer :: s
    real(c_double) :: old(3),x(1,1,1),y(1,1,1),z(1,1,1),zero(1,1,1)
    integer :: k,b
    status=2;s=>lookup(id)
    if(.not.associated(s))return
    status=1
    if(nt<0 .or. nt==huge(nt))return
    do k=1,nt
      do b=1,3
        if(.not.ieee_is_finite(e(b,k)))return
      end do
    end do
    call transfer_poles(s%ld,s%n,s%saved(:,:,1),s%saved(:,:,2),.false.)
    call totals(s,p(:,1),j(:,1));old=s%field;zero=0
    status=4
    do k=1,nt
      x=old(1);y=old(2);z=old(3)
      call ld_store_old(s%ld,[1,1,1],[1,1,1],[1,1,1],[1,1,1],s%dt,x,zero,y,zero,z,zero)
      x=e(1,k);y=e(2,k);z=e(3,k)
      if(s%tensor)then
        call tensor_update(s,old,e(:,k))
      else
        call ld_update(s%ld,[1,1,1],[1,1,1],[1,1,1],[1,1,1],s%dt,x,zero,y,zero,z,zero)
      endif
      call totals(s,p(:,k+1),j(:,k+1))
      old=e(:,k)
    end do
    ! Match the stateless contract: validate final poles and every returned sample.
    ! Intermediate pole histories are neither returned nor committed.
    if(.not.finite_state(s))then
      call transfer_poles(s%ld,s%n,s%saved(:,:,1),s%saved(:,:,2),.true.)
      return
    end if
    do k=1,nt+1
      do b=1,3
        if(.not.ieee_is_finite(p(b,k)) .or. .not.ieee_is_finite(j(b,k)))then
          call transfer_poles(s%ld,s%n,s%saved(:,:,1),s%saved(:,:,2),.true.)
          return
        end if
      end do
    end do
    s%field=old;status=0
  end function

  integer(c_int) function sp_ld_export(id,p,j) bind(C) result(status)
    integer(c_int64_t),value :: id
    real(c_double),intent(out) :: p(*),j(*)
    type(entry),pointer :: s
    status=2;s=>lookup(id)
    if(.not.associated(s))return
    call transfer_poles(s%ld,s%n,p(1:3*s%n),j(1:3*s%n),.false.);status=0
  end function

  integer(c_int) function sp_ld_set(id,e,p,j,pt,jt) bind(C) result(status)
    integer(c_int64_t),value :: id
    real(c_double),intent(in) :: e(3),p(*),j(*)
    real(c_double),intent(out) :: pt(3),jt(3)
    type(entry),pointer :: s
    integer :: b,q
    status=2;s=>lookup(id)
    if(.not.associated(s))return
    status=1
    do b=1,3
      if(.not.ieee_is_finite(e(b)))return
    end do
    do b=1,3*s%n
      if(.not.ieee_is_finite(p(b)) .or. .not.ieee_is_finite(j(b)))return
    end do
    ! Copy incoming data into checked persistent scratch, without changing accepted state.
    do q=1,s%n
      do b=1,3
        s%saved(b,q,1)=p(3*(q-1)+b)
        s%saved(b,q,2)=j(3*(q-1)+b)
      end do
    end do
    pt=sum(s%saved(:,:,1),dim=2);jt=sum(s%saved(:,:,2),dim=2)
    status=4
    do b=1,3
      if(.not.ieee_is_finite(pt(b)) .or. .not.ieee_is_finite(jt(b)))return
    end do
    call transfer_poles(s%ld,s%n,s%saved(:,:,1),s%saved(:,:,2),.true.)
    s%field=e;status=0
  end function

  integer(c_int) function sp_ld_destroy(id) bind(C) result(status)
    integer(c_int64_t),value :: id
    type(entry),pointer :: s,previous
    status=2;s=>first;nullify(previous)
    do while(associated(s))
      if(s%id==id)then
        if(associated(previous))then
          previous%next=>s%next
        else
          first=>s%next
        end if
        deallocate(s);status=0;return
      end if
      previous=>s;s=>s%next
    end do
  end function
end module
