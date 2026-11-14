! Persistent SBE ownership adapter. Shared kernels remain unchanged.
module sbe_handles
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use sbe_types
  use sbe_velocity_kernel
  implicit none
  private
  type model_node
    integer(c_int64_t) :: id=0,live=0
    integer :: order=0
    logical :: vnl=.false.
    type(s_sbe_gs_info) :: gs
    type(model_node),pointer :: next=>null()
  end type
  type state_node
    integer(c_int64_t) :: id=0
    real(c_double) :: dt=0d0,A(3)=0d0,J(3)=0d0
    type(s_sbe_bloch_solver) :: solver,trial
    type(model_node),pointer :: model=>null()
    type(state_node),pointer :: next=>null()
  end type
  type(model_node),pointer,save :: models=>null()
  type(state_node),pointer,save :: states=>null()
  integer(c_int64_t),save :: next_id=1
contains
  function find_model(id) result(node)
    integer(c_int64_t),intent(in) :: id
    type(model_node),pointer :: node
    node=>models
    do while(associated(node))
      if(node%id==id)return
      node=>node%next
    end do
  end function
  function find_state(id) result(node)
    integer(c_int64_t),intent(in) :: id
    type(state_node),pointer :: node
    node=>states
    do while(associated(node))
      if(node%id==id)return
      node=>node%next
    end do
  end function
  logical function vector_finite(v) result(ok)
    real(c_double),intent(in) :: v(3)
    integer :: k
    ok=.true.
    do k=1,3
      if(.not.ieee_is_finite(v(k)))ok=.false.
    end do
  end function
  logical function density_finite(rho) result(ok)
    complex(c_double_complex),intent(in) :: rho(:,:,:)
    integer :: i,j,k
    ok=.true.
    do k=1,size(rho,3)
      do j=1,size(rho,2)
        do i=1,size(rho,1)
          if(.not.ieee_is_finite(real(rho(i,j,k),c_double)))ok=.false.
          if(.not.ieee_is_finite(aimag(rho(i,j,k))))ok=.false.
        end do
      end do
    end do
  end function

  integer(c_int32_t) function sp_sbe_model_create(nb,nk,ne,order,vnl,volume,eigen,weight,p,rvnl,id) &
      bind(C) result(status)
    integer(c_int32_t),value :: nb,nk,ne,order,vnl
    real(c_double),value :: volume
    real(c_double),intent(in) :: eigen(nb,nk),weight(nk)
    complex(c_double_complex),intent(in) :: p(nb,nb,3,nk),rvnl(nb,nb,3,nk)
    integer(c_int64_t),intent(out) :: id
    type(model_node),pointer :: m
    integer(c_int64_t) :: cells
    integer :: i,j,k,d,ierr
    logical :: ok
    id=0;status=1
    if(nb<1.or.nk<1.or.ne<0.or.mod(ne,2)/=0.or.order<0.or.order>3.or.vnl<0.or.vnl>1)return
    if(int(ne,c_int64_t)>2_c_int64_t*int(nb,c_int64_t))return
    if(.not.ieee_is_finite(volume))return
    if(volume<=0d0)return
    status=3;cells=int(nb,c_int64_t)*int(nb,c_int64_t)
    if(cells>huge(cells)/int(nk,c_int64_t))return
    cells=cells*int(nk,c_int64_t)
    ! Checked address arithmetic, not a machine-wide memory budget.
    if(cells>huge(cells)/128_c_int64_t.or.next_id==huge(next_id))return
    ok=.true.
    do k=1,nk
      if(.not.ieee_is_finite(weight(k)))ok=.false.
      if(weight(k)<0d0)ok=.false.
      do i=1,nb
        if(.not.ieee_is_finite(eigen(i,k)))ok=.false.
      end do
      do d=1,3
        do j=1,nb
          do i=1,nb
            if(.not.ieee_is_finite(real(p(i,j,d,k),c_double)))ok=.false.
            if(.not.ieee_is_finite(aimag(p(i,j,d,k))))ok=.false.
            if(.not.ieee_is_finite(real(rvnl(i,j,d,k),c_double)))ok=.false.
            if(.not.ieee_is_finite(aimag(rvnl(i,j,d,k))))ok=.false.
          end do
        end do
      end do
    end do
    status=4
    if(.not.ok)return
    if(.not.ieee_is_finite(sum(weight)))return
    if(sum(weight)<=0d0)return
    status=3;allocate(m,stat=ierr)
    if(ierr/=0)return
    allocate(m%gs%eigen(nb,nk),m%gs%kweight(nk),m%gs%p_tm_matrix(nb,nb,3,nk), &
             m%gs%rvnl_tm_matrix(nb,nb,3,nk),m%gs%delta_omega(nb,nb,nk),m%gs%occup(nb,nk),stat=ierr)
    if(ierr/=0)then
      deallocate(m);return
    end if
    m%gs%nb=nb;m%gs%nk=nk;m%gs%ne=ne;m%gs%volume=volume;m%order=order;m%vnl=vnl/=0
    m%gs%eigen=eigen;m%gs%kweight=weight;m%gs%p_tm_matrix=p;m%gs%rvnl_tm_matrix=rvnl
    m%gs%occup=0d0;m%gs%occup(1:ne/2,:)=2d0
    ok=.true.
    do k=1,nk
      do j=1,nb
        do i=1,nb
          m%gs%delta_omega(i,j,k)=eigen(i,k)-eigen(j,k)
          if(.not.ieee_is_finite(m%gs%delta_omega(i,j,k)))ok=.false.
        end do
      end do
    end do
    status=4
    if(.not.ok)then
      deallocate(m);return
    end if
    m%id=next_id;next_id=next_id+1;m%next=>models;models=>m;id=m%id;status=0
  end function

  integer(c_int32_t) function sp_sbe_model_destroy(id) bind(C) result(status)
    integer(c_int64_t),value :: id
    type(model_node),pointer :: m,previous
    status=2;m=>models;nullify(previous)
    do while(associated(m))
      if(m%id==id)exit
      previous=>m;m=>m%next
    end do
    if(.not.associated(m))return
    status=5
    if(m%live/=0)return
    if(associated(previous))then
      previous%next=>m%next
    else
      models=>m%next
    end if
    deallocate(m);status=0
  end function

  subroutine evaluate(sbe,m,A,J)
    type(s_sbe_bloch_solver),intent(in) :: sbe
    type(model_node),intent(in) :: m
    real(c_double),intent(in) :: A(3)
    real(c_double),intent(out) :: J(3)
    complex(8) :: tmp(3)
    real(8) :: trace
    integer :: k,i
    call velocity_current_local(sbe,m%gs,A,m%order,tmp)
    J=real(tmp)/sum(m%gs%kweight)/m%gs%volume
    if(m%order<1)then
      trace=0d0
      do k=1,m%gs%nk
        do i=1,m%gs%nb
          trace=trace+real(sbe%rho(i,i,k))*m%gs%kweight(k)
        end do
      end do
      J=J+A*(trace/sum(m%gs%kweight))/m%gs%volume
    end if
  end subroutine

  integer(c_int32_t) function sp_sbe_state_create(model,dt,A,id,J) bind(C) result(status)
    integer(c_int64_t),value :: model
    real(c_double),value :: dt
    real(c_double),intent(in) :: A(3)
    integer(c_int64_t),intent(out) :: id
    real(c_double),intent(out) :: J(3)
    type(model_node),pointer :: m
    type(state_node),pointer :: s
    integer :: nb,nk,i,ierr
    id=0;status=2;m=>find_model(model)
    if(.not.associated(m))return
    status=1
    if(.not.ieee_is_finite(dt))return
    if(dt<=0d0.or..not.vector_finite(A))return
    status=3
    if(next_id==huge(next_id).or.m%live==huge(m%live))return
    allocate(s,stat=ierr)
    if(ierr/=0)return
    nb=m%gs%nb;nk=m%gs%nk
    allocate(s%solver%rho(nb,nb,nk),s%trial%rho(nb,nb,nk),stat=ierr)
    if(ierr/=0)then
      deallocate(s);return
    end if
    s%model=>m;s%dt=dt;s%A=A
    s%solver%nb=nb;s%solver%nk=nk;s%solver%ik_min=1;s%solver%ik_max=nk
    s%solver%flag_vnl_correction=m%vnl
    s%trial%nb=nb;s%trial%nk=nk;s%trial%ik_min=1;s%trial%ik_max=nk
    s%trial%flag_vnl_correction=m%vnl
    s%solver%rho=0d0
    do i=1,m%gs%ne/2
      s%solver%rho(i,i,:)=2d0
    end do
    call evaluate(s%solver,m,A,J)
    status=4
    if(.not.vector_finite(J))then
      deallocate(s);return
    end if
    s%J=J;s%id=next_id;next_id=next_id+1;s%next=>states;states=>s;m%live=m%live+1
    id=s%id;status=0
  end function

  integer(c_int32_t) function sp_sbe_step(id,A,J) bind(C) result(status)
    integer(c_int64_t),value :: id
    real(c_double),intent(in) :: A(3)
    real(c_double),intent(out) :: J(3)
    type(state_node),pointer :: s
    complex(c_double_complex),allocatable :: old(:,:,:)
    real(c_double) :: midpoint(3)
    s=>find_state(id);status=2
    if(.not.associated(s))return
    status=4
    if(.not.vector_finite(A))return
    midpoint=.5d0*(s%A+A)
    if(.not.vector_finite(midpoint))return
    s%trial%rho=s%solver%rho
    call velocity_evolve(s%trial,s%model%gs,midpoint,s%dt)
    call evaluate(s%trial,s%model,A,J)
    if(.not.density_finite(s%trial%rho).or..not.vector_finite(J))return
    call move_alloc(s%solver%rho,old)
    call move_alloc(s%trial%rho,s%solver%rho)
    call move_alloc(old,s%trial%rho)
    s%A=A;s%J=J;status=0
  end function

  integer(c_int32_t) function sp_sbe_density(id,rho) bind(C) result(status)
    integer(c_int64_t),value :: id
    complex(c_double_complex),intent(out) :: rho(*)
    type(state_node),pointer :: s
    integer :: i,j,k
    integer(c_int64_t) :: n
    s=>find_state(id);status=2
    if(.not.associated(s))return
    n=0
    do k=1,s%model%gs%nk
      do j=1,s%model%gs%nb
        do i=1,s%model%gs%nb
          n=n+1;rho(n)=s%solver%rho(i,j,k)
        end do
      end do
    end do
    status=0
  end function

  integer(c_int32_t) function sp_sbe_set(id,rho,A,J) bind(C) result(status)
    integer(c_int64_t),value :: id
    complex(c_double_complex),intent(in) :: rho(*)
    real(c_double),intent(in) :: A(3),J(3)
    type(state_node),pointer :: s
    integer :: i,b,k
    integer(c_int64_t) :: n
    logical :: ok
    s=>find_state(id);status=2
    if(.not.associated(s))return
    status=4;ok=vector_finite(A).and.vector_finite(J);n=0
    do k=1,s%model%gs%nk
      do b=1,s%model%gs%nb
        do i=1,s%model%gs%nb
          n=n+1
          if(.not.ieee_is_finite(real(rho(n),c_double)))ok=.false.
          if(.not.ieee_is_finite(aimag(rho(n))))ok=.false.
        end do
      end do
    end do
    if(.not.ok)return
    n=0
    do k=1,s%model%gs%nk
      do b=1,s%model%gs%nb
        do i=1,s%model%gs%nb
          n=n+1;s%solver%rho(i,b,k)=rho(n)
        end do
      end do
    end do
    s%A=A;s%J=J;status=0
  end function

  integer(c_int32_t) function sp_sbe_state_destroy(id) bind(C) result(status)
    integer(c_int64_t),value :: id
    type(state_node),pointer :: s,previous
    status=2;s=>states;nullify(previous)
    do while(associated(s))
      if(s%id==id)exit
      previous=>s;s=>s%next
    end do
    if(.not.associated(s))return
    if(associated(previous))then
      previous%next=>s%next
    else
      states=>s%next
    end if
    s%model%live=s%model%live-1;deallocate(s);status=0
  end function
end module
