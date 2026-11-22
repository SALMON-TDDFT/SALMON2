! Electric-impulse response of an evolving occupied density under a frozen H.
! This diagnostic is not a stationary Gamma-cell screened interaction.
module hsex_coherent_optical
  implicit none
  private
  public :: coherent_optical,stationary_current_reference
  abstract interface
    subroutine operator_action(x,y,status)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine stationary_current_reference(apply,velocity,acceleration,psi,occupation,volume,tol,dc,alpha,status)
    use hsex_frozen_response, only: frozen_background_rate
    use hsex_longwave, only: longwave_band
    use hsex_screen_math, only: screen_finite
    implicit none
    procedure(operator_action) :: apply,velocity,acceleration
    complex(8),intent(in) :: psi(:,:)
    real(8),intent(in) :: occupation,volume,tol
    real(8),intent(out) :: dc,alpha
    integer,intent(out) :: status
    complex(8),allocatable :: hp(:,:),matrix(:,:),orbitals(:,:),vel(:,:),body(:,:),work(:),mass(:)
    real(8),allocatable :: energy(:),rwork(:)
    complex(8) :: part(4,4)
    real(8) :: rate
    integer :: n,nv,a,code
    dc=0d0;alpha=0d0;status=1;n=size(psi,1);nv=size(psi,2)
    if(nv<1.or.n<nv.or..not.screen_finite(psi))return
    if(.not.screen_finite([occupation,volume,tol]).or.volume<=0d0.or.tol<=0d0)return
    if(occupation<0d0.or.occupation>2d0)return
    allocate(hp(n,nv),matrix(nv,nv),orbitals(n,nv),vel(n,3),body(n,1),mass(n), &
      work(max(1,4*nv)),energy(nv),rwork(max(1,3*nv-2)),stat=code)
    if(code/=0)return
    do a=1,nv
      call apply(psi(:,a),hp(:,a),code)
      if(code/=0)return
      if(.not.screen_finite(hp(:,a)))return
    enddo
    call frozen_background_rate(psi,hp,rate,code)
    if(code/=0)return
    status=2
    if(rate>tol)return
    matrix=matmul(conjg(transpose(psi)),hp)
    if(maxval(abs(matrix-conjg(transpose(matrix))))>tol)return
    matrix=.5d0*(matrix+conjg(transpose(matrix)))
    call zheev('V','U',nv,matrix,nv,energy,work,size(work),rwork,code)
    if(code/=0)return
    orbitals=matmul(psi,matrix);status=3
    do a=1,nv
      vel=0d0
      call velocity(orbitals(:,a),vel(:,3),code)
      if(code/=0)return
      if(.not.screen_finite(vel))return
      body(:,1)=vel(:,3)
      call longwave_band(apply,orbitals,energy(a),vel,body,part,code)
      if(code/=0)return
      call acceleration(orbitals(:,a),mass,code)
      if(code/=0)return
      if(.not.screen_finite(mass))return
      dc=dc+occupation*(real(dot_product(orbitals(:,a),mass),8)-2d0*real(part(4,4),8))/volume
      alpha=alpha+2d0*occupation*real(part(3,3),8)/volume
    enddo
    if(.not.screen_finite([dc,alpha]))return
    status=0
  end subroutine

  subroutine coherent_optical(apply,velocity,acceleration,psi,occ,vertices,z,capacity,tol,volume, &
      conductivity,density,residual,actions,status,seed)
    use hsex_screen_math, only: screen_finite
    use hsex_frozen_response, only: frozen_solve,frozen_trace,s_frozen_solution,s_frozen_stats
    implicit none
    procedure(operator_action) :: apply,velocity,acceleration
    complex(8),intent(in) :: psi(:,:),vertices(:,:),z(:)
    real(8),intent(in) :: occ(:),tol,volume
    integer,intent(in) :: capacity
    complex(8),intent(out) :: conductivity(:),density(:,:)
    real(8),intent(out) :: residual(:,:)
    integer,intent(out) :: actions,status
    complex(8),optional,intent(in) :: seed(:,:)
    type(s_frozen_solution) :: background,response
    type(s_frozen_stats) :: statistics
    complex(8),allocatable :: c(:,:),left_factor(:,:),right_factor(:,:),vl(:,:),vr(:,:),contact(:),current(:)
    complex(8),allocatable :: flat(:),occupation(:,:),initial(:,:)
    integer :: n,nv,nz,nl,nr,j,s,code,ni
    status=1;actions=0;conductivity=0d0;density=0d0;residual=huge(1d0)
    n=size(psi,1);nv=size(psi,2);nz=size(z)
    if(n<1.or.nv<1.or.nz<1.or.size(occ)/=nv.or.size(vertices,1)/=n)return
    if(size(conductivity)/=nz.or.any(shape(density)/=[size(vertices,2),nz]))return
    if(any(shape(residual)/=[2,nz]).or.capacity<1)return
    if(.not.screen_finite(psi).or..not.screen_finite(occ).or..not.screen_finite(vertices))return
    if(.not.screen_finite(z).or..not.screen_finite([tol,volume]))return
    if(volume<=0d0.or.tol<=0d0.or.any(aimag(z)<=0d0).or.any(occ<0d0).or.any(occ>2d0))return
    ni=nv
    if(present(seed))then
      if(size(seed,1)/=n.or..not.screen_finite(seed))return
      ni=ni+size(seed,2)
    endif
    allocate(occupation(ni,ni),initial(n,ni),contact(nz),current(1),flat(1),stat=code)
    if(code/=0)return
    occupation=0d0
    initial(:,:nv)=psi
    if(present(seed))initial(:,nv+1:)=seed
    do j=1,nv
      occupation(j,j)=occ(j)
    enddo
    ! R rho, not rho/z: the latter silently removes coherent background motion.
    ! Extra seed columns have ZERO coefficients: the source is exactly rho,
    ! with no change of occupations or accepted residual tolerance.
    call frozen_solve(apply,apply,initial,occupation,initial,z,capacity,tol,background,statistics,code, &
      residual_bound=.true.)
    actions=statistics%h_actions;status=2
    if(code/=0)return
    residual(1,:)=background%residual
    call traces(background,acceleration,contact,code)
    status=3
    if(code/=0)return
    nl=size(background%left,2);nr=size(background%right,2)
    allocate(vl(n,nl),vr(n,nr),left_factor(n,2*nl),right_factor(n,2*nr),c(2*nl,2*nr),stat=code)
    if(code/=0)return
    status=4
    do j=1,nl
      call velocity(background%left(:,j),vl(:,j),code)
      if(code/=0)return
      if(.not.screen_finite(vl(:,j)))return
    enddo
    do j=1,nr
      call velocity(background%right(:,j),vr(:,j),code)
      if(code/=0)return
      if(.not.screen_finite(vr(:,j)))return
    enddo
    left_factor(:,:nl)=vl;left_factor(:,nl+1:)=background%left
    right_factor(:,:nr)=background%right;right_factor(:,nr+1:)=vr
    do s=1,nz
      current=0d0;residual(2,s)=0d0
      ! An identically zero velocity has an exactly zero commutator, while
      ! its contact term can remain finite (free-electron Drude response).
      if(any(vl/=(0d0,0d0)).or.any(vr/=(0d0,0d0)))then
        c=0d0;c(:nl,:nr)=background%coefficient(:,:,s)
        c(nl+1:,nr+1:)=-background%coefficient(:,:,s)
        call frozen_solve(apply,apply,left_factor,c,right_factor,z(s:s),capacity,tol,response,statistics,code, &
          residual_bound=.true.)
        actions=actions+statistics%h_actions;status=5
        if(code/=0)return
        residual(2,s)=response%residual(1)
        call traces(response,velocity,current,code)
        status=6
        if(code/=0)return
        do j=1,size(vertices,2)
          call frozen_trace(response,vertices(:,j),flat,code)
          if(code/=0)return
          density(j,s)=cmplx(0d0,-1d0,8)*flat(1)
        enddo
      endif
      conductivity(s)=cmplx(0d0,1d0,8)*(current(1)+contact(s))/volume
    enddo
    if(.not.screen_finite(conductivity).or..not.screen_finite(density))return
    status=0
  end subroutine

  subroutine traces(solution,apply,values,status)
    use hsex_frozen_response, only: s_frozen_solution
    use hsex_screen_math, only: screen_finite
    implicit none
    type(s_frozen_solution),intent(in) :: solution
    procedure(operator_action) :: apply
    complex(8),intent(out) :: values(:)
    integer,intent(out) :: status
    complex(8),allocatable :: matrix(:,:),vector(:)
    integer :: j,s,code
    status=1;values=0d0
    if(size(values)/=size(solution%coefficient,3))return
    allocate(matrix(size(solution%right,2),size(solution%left,2)),vector(size(solution%left,1)),stat=code)
    if(code/=0)return
    do j=1,size(solution%left,2)
      call apply(solution%left(:,j),vector,code)
      if(code/=0)return
      if(.not.screen_finite(vector))return
      matrix(:,j)=matmul(conjg(transpose(solution%right)),vector)
    enddo
    do s=1,size(values)
      values(s)=sum(transpose(matrix)*solution%coefficient(:,:,s))
    enddo
    if(.not.screen_finite(values))return
    status=0
  end subroutine
end module
