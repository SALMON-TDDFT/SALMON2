! Unitary space-group star reconstruction on an exact signed-permutation fractional grid.
! Sources average little-group projectors, not individual orbital gauges.
module exx_symmetry
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 public :: s_exx_symmetry_map,symmetry_init,symmetry_expand,symmetry_validate_atoms,symmetry_transform
 type s_exx_symmetry_map
   integer :: nfull=0,nrep=0,nsym=0,max_little=0
   integer :: n(3)=0
   real(8) :: h(3)=0
   real(8) :: cell(3,3)=0d0,inverse_cell(3,3)=0d0
   real(8),allocatable :: cart_rotation(:,:,:)
   real(8),allocatable :: full_k(:,:),rep_k(:,:),rotation(:,:,:),translation(:,:)
   integer,allocatable :: first(:),owner(:),operations(:,:),multiplicity(:),grid(:,:)
 end type
  ! Rank-specific loops keep IEEE inquiries scalar on Fujitsu compilers.
  private :: salmon_all_finite,finite_real_1d,finite_real_2d,finite_real_3d
  interface salmon_all_finite
    module procedure finite_real_1d,finite_real_2d,finite_real_3d
  end interface
 interface symmetry_init
 module procedure symmetry_init_scalar,symmetry_init_vector
 end interface
 public :: symmetry_spin_rotation
contains
 subroutine symmetry_init_scalar(map,n,h,krep,weights,SymMatA,SymMatB,mesh,ierr,cell)
 implicit none
 type(s_exx_symmetry_map),intent(out)::map
 integer,intent(in)::n(3),mesh
 real(8),intent(in)::h(3),krep(:,:),weights(:),SymMatA(:,:,:),SymMatB(:,:,:)
 integer,intent(out)::ierr
 real(8),intent(in),optional::cell(3,3)
 call symmetry_init_vector(map,n,h,krep,weights,SymMatA,SymMatB,[mesh,mesh,mesh],ierr,cell)
 end subroutine
 subroutine symmetry_init_vector(map,n,h,krep,weights,SymMatA,SymMatB,mesh,ierr,cell)
 implicit none
 type(s_exx_symmetry_map),intent(out)::map
 integer,intent(in):: n(3),mesh(3)
 real(8),intent(in)::h(3),krep(:,:),weights(:),SymMatA(:,:,:),SymMatB(:,:,:)
 integer,intent(out)::ierr
 real(8),intent(in),optional::cell(3,3)
 real(8):: L(3),pi2,ident(3,3),rot(3,3),v(3),f(3),offset(3),xyz(3),cart(3,3),det,cell_scale
 integer:: nr,ns,nf,ng,r,s,d,g,j,i,id,slot,base,c(3),c0(3),nstar,matches
 integer,allocatable::lookup(:),order(:),used(:)
 ierr=1
 if(any(n<1).or.any(mesh<1).or.any(h<=0).or..not.salmon_all_finite(h))return
 if(size(krep,1)/=3.or.size(SymMatA,1)/=3.or.size(SymMatA,2)/=4)return
 if(any(shape(SymMatA)/=shape(SymMatB)))return
 nr=size(krep,2);ns=size(SymMatA,3);nf=product(mesh);ng=product(n)
 if(nr<1.or.ns<1.or.size(weights)/=nr)return
 if(.not.salmon_all_finite(krep).or..not.salmon_all_finite(weights))return
 if(.not.salmon_all_finite(SymMatA).or..not.salmon_all_finite(SymMatB))return
 L=n*h;pi2=2*acos(-1d0);ident=0
 do d=1,3
 ident(d,d)=1
 enddo

 map%cell=0d0
 do d=1,3
 map%cell(d,d)=L(d)
 enddo
 if(present(cell))map%cell=cell
 if(.not.salmon_all_finite(map%cell))return
 cell_scale=maxval(abs(map%cell))
 if(cell_scale==0d0)return
 det=determinant3(map%cell/cell_scale)
 if(.not.ieee_is_finite(det).or.abs(det)<=1d-14)return
 call inverse3(map%cell/cell_scale,map%inverse_cell,det)
 map%inverse_cell=map%inverse_cell/cell_scale
 if(.not.salmon_all_finite(map%inverse_cell))return
 allocate(map%cart_rotation(3,3,ns))
 allocate(order(ns));id=0
 do s=1,ns
 rot=SymMatA(:,1:3,s)
 if(maxval(abs(rot-anint(rot)))>1d-10)return
 if(maxval(abs(matmul(transpose(rot),rot)-ident))>1d-10)return
 ! Fractional signed permutations have inverse transpose equal to themselves.
 if(maxval(abs(rot-SymMatB(:,1:3,s)))>1d-10)return
 cart=matmul(map%cell,matmul(rot,map%inverse_cell))
 if(.not.salmon_all_finite(cart))return
 if(maxval(abs(matmul(transpose(cart),cart)-ident))>1d-10)return
 map%cart_rotation(:,:,s)=cart
 if(maxval(abs(rot-ident))<1d-10.and.maxval(abs(SymMatA(:,4,s)))<1d-10)id=s
 enddo
 ! A supplied operation list must be a finite group modulo lattice translations.
 ! Duplicate operations would bias the projector average and are rejected.
 do r=1,ns
 do s=1,ns
 rot=matmul(SymMatA(:,1:3,r),SymMatA(:,1:3,s))
 v=matmul(SymMatA(:,1:3,r),SymMatA(:,4,s))+SymMatA(:,4,r)
 matches=0
 do j=1,ns
 f=v-SymMatA(:,4,j);f=f-anint(f)
 if(maxval(abs(rot-SymMatA(:,1:3,j)))<1d-10.and.maxval(abs(f))<1d-10)matches=matches+1
 enddo
 if(matches/=1)return
 enddo
 enddo
 if(id==0)return
 order(1)=id;j=1
 do s=1,ns
 if(s==id)cycle
 j=j+1;order(j)=s
 enddo
 map%n=n;map%h=h;map%nrep=nr;map%nsym=ns;map%nfull=nf
 allocate(map%full_k(3,nf),map%rep_k(3,nr),map%rotation(3,3,ns),map%translation(3,ns))
 allocate(map%first(nr+1),map%owner(nf),map%operations(ns,nf),map%multiplicity(nf),map%grid(ng,ns))
 map%rep_k=krep;map%rotation=SymMatB(:,1:3,:);map%translation=SymMatA(:,4,:)
 map%multiplicity=0;map%operations=0
 allocate(lookup(nf),used(ng));lookup=0
 ! SALMON periodic grid is indexed from one, coordinate (index-1)*h.
 do s=1,ns
 used=0
 do g=1,ng
 c0=[mod(g-1,n(1)),mod((g-1)/n(1),n(2)),(g-1)/(n(1)*n(2))]
 xyz=matmul(transpose(SymMatA(:,1:3,s)),real(c0,8)/n-SymMatA(:,4,s))*n
 if(maxval(abs(xyz-anint(xyz)))>1d-8)return
 c=modulo(nint(xyz),n)
 j=1+c(1)+n(1)*(c(2)+n(2)*c(3))
 if(used(j)/=0)return
 used(j)=1;map%grid(g,s)=j
 enddo
 enddo
 offset=matmul(transpose(map%cell),krep(:,1))/pi2*mesh;offset=offset-floor(offset)
 i=0
 do r=1,nr
 base=i;map%first(r)=i+1
 do slot=1,ns
 s=order(slot)
 v=matmul(SymMatB(:,1:3,s),matmul(transpose(map%cell),krep(:,r))/pi2)
 f=v*mesh-offset
 if(maxval(abs(f-anint(f)))>1d-8)return
 c=modulo(nint(f),mesh)
 j=1+c(1)+mesh(1)*(c(2)+mesh(2)*c(3))
 g=lookup(j)
 if(g==0)then
 i=i+1
 if(i>nf)return
 g=i;lookup(j)=g;map%owner(g)=r
 if(slot==1)then
 map%full_k(:,g)=krep(:,r)
 else
 map%full_k(:,g)=pi2*matmul(transpose(map%inverse_cell),v-floor(v+0.5d0))
 endif
 else
 if(map%owner(g)/=r)return
 endif
 d=map%multiplicity(g)+1
 map%multiplicity(g)=d;map%operations(d,g)=s
 enddo
 nstar=i-base
 if(abs(weights(r)-real(nstar,8)/nf)>1d-10)return
 enddo
 if(i/=nf.or.any(lookup==0))return
 map%first(nr+1)=nf+1
 map%max_little=maxval(map%multiplicity)
 ierr=0
 end subroutine

 ! Verify species-preserving ionic permutations before using a user symmetry file.
 subroutine symmetry_validate_atoms(map,rion,kion,ierr)
 implicit none
 type(s_exx_symmetry_map),intent(in)::map
 real(8),intent(in)::rion(:,:)
 integer,intent(in)::kion(:)
 integer,intent(out)::ierr
 integer::na,s,i,j,match,count_matches
 integer,allocatable::used(:)
 real(8)::v(3),delta(3)
 ierr=1;na=size(kion)
 if(size(rion,1)/=3.or.size(rion,2)/=na.or.na<1)return
 if(.not.salmon_all_finite(rion))return
 if(.not.allocated(map%rotation).or..not.allocated(map%translation))return
 if(.not.salmon_all_finite(map%cell).or..not.salmon_all_finite(map%inverse_cell))return
 allocate(used(na))
 do s=1,map%nsym
 used=0
 do i=1,na
 v=matmul(map%rotation(:,:,s),matmul(map%inverse_cell,rion(:,i)))+map%translation(:,s)
 count_matches=0;match=0
 do j=1,na
 if(kion(i)/=kion(j))cycle
 delta=v-matmul(map%inverse_cell,rion(:,j));delta=delta-anint(delta)
 if(sqrt(sum(matmul(map%cell,delta)**2))<1d-7)then
 count_matches=count_matches+1;match=j
 endif
 enddo
 if(count_matches/=1)return
 if(used(match)/=0)return
 used(match)=1
 enddo
 enddo
 ierr=0
 end subroutine

 ! Transform one representative orbital block; no little-group normalization.
 subroutine symmetry_transform(map,full_index,operation_slot,orbitals,transformed,ierr)
 implicit none
 type(s_exx_symmetry_map),intent(in)::map
 integer,intent(in)::full_index,operation_slot
 complex(8),intent(in)::orbitals(:,:)
 complex(8),intent(out)::transformed(:,:)
 integer,intent(out)::ierr
 integer::ng,r,s,g,c(3),nc
 complex(8)::spin(2,2)
 real(8)::rk(3),Gvec(3),pos(3),angle
 complex(8)::phase
 ierr=1;ng=product(map%n)
 if(.not.allocated(map%multiplicity))return
 if(full_index<1.or.full_index>map%nfull)return
 if(operation_slot<1.or.operation_slot>map%multiplicity(full_index))return
 nc=size(orbitals,1)/ng
 if((nc/=1.and.nc/=2).or.size(orbitals,1)/=nc*ng.or.any(shape(orbitals)/=shape(transformed)))return
 if(.not.salmon_all_finite(real(orbitals)).or..not.salmon_all_finite(aimag(orbitals)))return
 r=map%owner(full_index);s=map%operations(operation_slot,full_index)
 if(nc==2)then
 call symmetry_spin_rotation(map%cart_rotation(:,:,s),spin,ierr)
 if(ierr/=0)return
 endif
 rk=matmul(map%cart_rotation(:,:,s),map%rep_k(:,r))
 Gvec=rk-map%full_k(:,full_index)
 do g=1,ng
 c=[mod(g-1,map%n(1)),mod((g-1)/map%n(1),map%n(2)),(g-1)/(map%n(1)*map%n(2))]
 pos=matmul(map%cell,real(c,8)/map%n)
 angle=dot_product(Gvec,pos)-dot_product(rk,matmul(map%cell,map%translation(:,s)))
 phase=cmplx(cos(angle),sin(angle),8)
 if(nc==1)then
 transformed(g,:)=phase*orbitals(map%grid(g,s),:)
 else
 transformed(g,:)=phase*(spin(1,1)*orbitals(map%grid(g,s),:)+spin(1,2)*orbitals(ng+map%grid(g,s),:))
 transformed(ng+g,:)=phase*(spin(2,1)*orbitals(map%grid(g,s),:)+spin(2,2)*orbitals(ng+map%grid(g,s),:))
 endif
 enddo
 ierr=0
 end subroutine

 subroutine symmetry_expand(map,rep_start,source,target,expanded_source,expanded_target,ierr)
 type(s_exx_symmetry_map),intent(in)::map
 integer,intent(in)::rep_start
 complex(8),intent(in)::source(:,:,:),target(:,:,:)
 complex(8),allocatable,intent(out)::expanded_source(:,:,:),expanded_target(:,:,:)
 integer,intent(out)::ierr
 integer::nr,ng,no,nv,lo,hi,j,l,r,m
 real(8)::scale
 ierr=1;ng=size(source,1);no=size(source,2);nv=size(target,2);nr=size(source,3)
 if(.not.allocated(map%first))return
 if(rep_start<1.or.rep_start+nr>map%nrep+1)return
 if(size(target,1)/=ng.or.size(target,3)/=nr)return
 if(no<1.or.nv<1)return
 lo=map%first(rep_start);hi=map%first(rep_start+nr)-1
 allocate(expanded_source(ng,no*map%max_little,hi-lo+1),expanded_target(ng,nv,hi-lo+1))
 expanded_source=0;expanded_target=0
 do j=lo,hi
 r=map%owner(j)-rep_start+1;l=j-lo+1;scale=1/sqrt(real(map%multiplicity(j),8))
 do m=1,map%multiplicity(j)
 call symmetry_transform(map,j,m,source(:,:,r),expanded_source(:,(m-1)*no+1:m*no,l),ierr)
 if(ierr/=0)return
 expanded_source(:,(m-1)*no+1:m*no,l)=expanded_source(:,(m-1)*no+1:m*no,l)*scale
 enddo
 call symmetry_transform(map,j,1,target(:,:,r),expanded_target(:,:,l),ierr)
 if(ierr/=0)return
 enddo
 ierr=0
 end subroutine


 ! Axial SO(3) lift for proper/improper Cartesian isometries. The signed
 ! diagonal branch retains the established gauge exactly; other operations use
 ! a stable largest-diagonal quaternion. Overall SU(2) signs cancel in projectors.
 subroutine symmetry_spin_rotation(rotation,u,ierr)
 implicit none
 real(8),intent(in)::rotation(3,3)
 complex(8),intent(out)::u(2,2)
 integer,intent(out)::ierr
 real(8)::r(3),off(3,3),det,a(3,3),ident(3,3),q(4),t,norm
 integer::i,j,k
 ierr=1;u=0
 if(.not.salmon_all_finite(rotation))return
 ident=0d0
 do i=1,3
 ident(i,i)=1d0
 enddo
 if(maxval(abs(matmul(transpose(rotation),rotation)-ident))>1d-10)return
 det=determinant3(rotation)
 if(abs(abs(det)-1d0)>1d-10)return
 off=rotation
 do i=1,3
 off(i,i)=0d0;r(i)=rotation(i,i)
 enddo
 if(maxval(abs(off))<1d-12.and.all(abs(abs(r)-1d0)<1d-12))then
 r=sign(1d0,det)*r
 if(all(r>0))then
 u(1,1)=1;u(2,2)=1
 else if(r(1)>0)then
 u(1,2)=(0d0,-1d0);u(2,1)=(0d0,-1d0)
 else if(r(2)>0)then
 u(1,2)=-1;u(2,1)=1
 else
 u(1,1)=(0d0,-1d0);u(2,2)=(0d0,1d0)
 endif
 else
 a=sign(1d0,det)*rotation
 t=a(1,1)+a(2,2)+a(3,3)
 if(t>0d0)then
 t=2d0*sqrt(1d0+t)
 q(1)=t/4d0
 q(2)=(a(3,2)-a(2,3))/t
 q(3)=(a(1,3)-a(3,1))/t
 q(4)=(a(2,1)-a(1,2))/t
 else
 i=maxloc([a(1,1),a(2,2),a(3,3)],dim=1);j=mod(i,3)+1;k=mod(j,3)+1
 t=2d0*sqrt(max(0d0,1d0+a(i,i)-a(j,j)-a(k,k)))
 if(t<=1d-14)return
 q(1)=(a(k,j)-a(j,k))/t;q(i+1)=t/4d0
 q(j+1)=(a(j,i)+a(i,j))/t;q(k+1)=(a(k,i)+a(i,k))/t
 endif
 norm=sqrt(sum(q*q))
 if(norm<=1d-14)return
 q=q/norm
 u(1,1)=cmplx(q(1),-q(4),8);u(2,2)=cmplx(q(1),q(4),8)
 u(1,2)=cmplx(-q(3),-q(2),8);u(2,1)=cmplx(q(3),-q(2),8)
 endif
 ierr=0
 end subroutine

 pure real(8) function determinant3(a) result(det)
 implicit none
 real(8),intent(in)::a(3,3)
 det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
    -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
    +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
 end function
 subroutine inverse3(a,b,det)
 implicit none
 real(8),intent(in)::a(3,3)
 real(8),intent(out)::b(3,3),det
 integer::i,j,k,l
 det=determinant3(a);b=0d0
 if(det==0d0)return
 do i=1,3
 j=mod(i,3)+1;k=mod(j,3)+1
 do l=1,3
 b(i,l)=(a(mod(l,3)+1,j)*a(mod(l+1,3)+1,k)-a(mod(l+1,3)+1,j)*a(mod(l,3)+1,k))/det
 enddo
 enddo
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

  pure logical function finite_real_2d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:)
    real(8) :: value
    integer :: i,j
    finite=.false.
    do j=1,size(values,2)
      do i=1,size(values,1)
        value=values(i,j)
        if(.not.ieee_is_finite(value))return
      enddo
    enddo
    finite=.true.
  end function

  pure logical function finite_real_3d(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(8),intent(in) :: values(:,:,:)
    real(8) :: value
    integer :: i,j,k
    finite=.false.
    do k=1,size(values,3)
      do j=1,size(values,2)
        do i=1,size(values,1)
          value=values(i,j,k)
          if(.not.ieee_is_finite(value))return
        enddo
      enddo
    enddo
    finite=.true.
  end function
end module
