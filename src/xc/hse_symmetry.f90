! Unitary space-group star reconstruction for a complete cubic k mesh.
! Sources average little-group projectors, not individual orbital gauges.
module hse_symmetry
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 public :: hse_symmetry_map,symmetry_init,symmetry_expand,symmetry_validate_atoms,symmetry_transform
 type hse_symmetry_map
   integer :: nfull=0,nrep=0,nsym=0,max_little=0
   integer :: n(3)=0
   real(8) :: h(3)=0
   real(8),allocatable :: full_k(:,:),rep_k(:,:),rotation(:,:,:),translation(:,:)
   integer,allocatable :: first(:),owner(:),operations(:,:),multiplicity(:),grid(:,:)
 end type
contains
 subroutine symmetry_init(map,n,h,krep,weights,SymMatA,SymMatB,mesh,ierr)
 type(hse_symmetry_map),intent(out)::map
 integer,intent(in):: n(3),mesh
 real(8),intent(in)::h(3),krep(:,:),weights(:),SymMatA(:,:,:),SymMatB(:,:,:)
 integer,intent(out)::ierr
 real(8):: L(3),pi2,ident(3,3),rot(3,3),v(3),f(3),offset(3),xyz(3)
 integer:: nr,ns,nf,ng,r,s,d,g,j,i,id,slot,base,c(3),c0(3),nstar,matches
 integer,allocatable::lookup(:),order(:),used(:)
 ierr=1
 if(any(n<1).or.mesh<1.or.any(h<=0).or..not.all(ieee_is_finite(h)))return
 if(size(krep,1)/=3.or.size(SymMatA,1)/=3.or.size(SymMatA,2)/=4)return
 if(any(shape(SymMatA)/=shape(SymMatB)))return
 nr=size(krep,2);ns=size(SymMatA,3);nf=mesh**3;ng=product(n)
 if(nr<1.or.ns<1.or.size(weights)/=nr)return
 if(.not.all(ieee_is_finite(krep)).or..not.all(ieee_is_finite(weights)))return
 if(.not.all(ieee_is_finite(SymMatA)).or..not.all(ieee_is_finite(SymMatB)))return
 L=n*h;pi2=2*acos(-1d0);ident=0
 do d=1,3
 ident(d,d)=1
 enddo
 if(maxval(abs(L-L(1)))>1d-10)return
 allocate(order(ns));id=0
 do s=1,ns
 rot=SymMatA(:,1:3,s)
 if(maxval(abs(rot-anint(rot)))>1d-10)return
 if(maxval(abs(matmul(transpose(rot),rot)-ident))>1d-10)return
 if(maxval(abs(rot-SymMatB(:,1:3,s)))>1d-10)return
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
 offset=krep(:,1)*L/pi2*mesh;offset=offset-floor(offset)
 i=0
 do r=1,nr
 base=i;map%first(r)=i+1
 do slot=1,ns
 s=order(slot)
 v=matmul(SymMatB(:,1:3,s),krep(:,r)*L/pi2)
 f=v*mesh-offset
 if(maxval(abs(f-anint(f)))>1d-8)return
 c=modulo(nint(f),mesh)
 j=1+c(1)+mesh*(c(2)+mesh*c(3))
 g=lookup(j)
 if(g==0)then
 i=i+1
 if(i>nf)return
 g=i;lookup(j)=g;map%owner(g)=r
 if(slot==1)then
 map%full_k(:,g)=krep(:,r)
 else
 map%full_k(:,g)=(v-floor(v+0.5d0))*pi2/L
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
 type(hse_symmetry_map),intent(in)::map
 real(8),intent(in)::rion(:,:)
 integer,intent(in)::kion(:)
 integer,intent(out)::ierr
 integer::na,s,i,j,match,count_matches
 integer,allocatable::used(:)
 real(8)::L(3),v(3),delta(3)
 ierr=1;na=size(kion)
 if(size(rion,1)/=3.or.size(rion,2)/=na.or.na<1)return
 if(.not.all(ieee_is_finite(rion)))return
 if(.not.allocated(map%rotation).or..not.allocated(map%translation))return
 L=map%n*map%h
 if(any(L<=0))return
 allocate(used(na))
 do s=1,map%nsym
 used=0
 do i=1,na
 v=matmul(map%rotation(:,:,s),rion(:,i)/L)+map%translation(:,s)
 count_matches=0;match=0
 do j=1,na
 if(kion(i)/=kion(j))cycle
 delta=v-rion(:,j)/L;delta=delta-anint(delta)
 if(maxval(abs(delta*L))<1d-7)then
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
 type(hse_symmetry_map),intent(in)::map
 integer,intent(in)::full_index,operation_slot
 complex(8),intent(in)::orbitals(:,:)
 complex(8),intent(out)::transformed(:,:)
 integer,intent(out)::ierr
 integer::ng,r,s,g,c(3)
 real(8)::rk(3),Gvec(3),pos(3),angle
 complex(8)::phase
 ierr=1;ng=product(map%n)
 if(.not.allocated(map%multiplicity))return
 if(full_index<1.or.full_index>map%nfull)return
 if(operation_slot<1.or.operation_slot>map%multiplicity(full_index))return
 if(size(orbitals,1)/=ng.or.any(shape(orbitals)/=shape(transformed)))return
 if(.not.all(ieee_is_finite(real(orbitals))).or..not.all(ieee_is_finite(aimag(orbitals))))return
 r=map%owner(full_index);s=map%operations(operation_slot,full_index)
 rk=matmul(map%rotation(:,:,s),map%rep_k(:,r))
 Gvec=rk-map%full_k(:,full_index)
 do g=1,ng
 c=[mod(g-1,map%n(1)),mod((g-1)/map%n(1),map%n(2)),(g-1)/(map%n(1)*map%n(2))]
 pos=c*map%h
 angle=dot_product(Gvec,pos)-dot_product(rk,map%translation(:,s)*map%n*map%h)
 phase=cmplx(cos(angle),sin(angle),8)
 transformed(g,:)=phase*orbitals(map%grid(g,s),:)
 enddo
 ierr=0
 end subroutine

 subroutine symmetry_expand(map,rep_start,source,target,expanded_source,expanded_target,ierr)
 type(hse_symmetry_map),intent(in)::map
 integer,intent(in)::rep_start
 complex(8),intent(in)::source(:,:,:),target(:,:,:)
 complex(8),allocatable,intent(out)::expanded_source(:,:,:),expanded_target(:,:,:)
 integer,intent(out)::ierr
 integer::nr,ng,no,nv,lo,hi,j,l,r,s,g,m,c(3),q
 real(8)::rk(3),Gvec(3),pos(3),angle,scale
 complex(8)::phase
 ierr=1;ng=product(map%n);no=size(source,2);nv=size(target,2);nr=size(source,3)
 if(.not.allocated(map%first))return
 if(rep_start<1.or.rep_start+nr>map%nrep+1)return
 if(size(source,1)/=ng.or.size(target,1)/=ng.or.size(target,3)/=nr)return
 if(no<1.or.nv<1)return
 if(.not.all(ieee_is_finite(real(source))).or..not.all(ieee_is_finite(aimag(source))))return
 if(.not.all(ieee_is_finite(real(target))).or..not.all(ieee_is_finite(aimag(target))))return
 lo=map%first(rep_start);hi=map%first(rep_start+nr)-1
 allocate(expanded_source(ng,no*map%max_little,hi-lo+1),expanded_target(ng,nv,hi-lo+1))
 expanded_source=0;expanded_target=0
 do j=lo,hi
 r=map%owner(j);l=j-lo+1;scale=1/sqrt(real(map%multiplicity(j),8))
 do m=1,map%multiplicity(j)
 s=map%operations(m,j)
 rk=matmul(map%rotation(:,:,s),map%rep_k(:,r))
 Gvec=rk-map%full_k(:,j)
 do g=1,ng
 c=[mod(g-1,map%n(1)),mod((g-1)/map%n(1),map%n(2)),(g-1)/(map%n(1)*map%n(2))]
 pos=c*map%h
 angle=dot_product(Gvec,pos)-dot_product(rk,map%translation(:,s)*map%n*map%h)
 phase=cmplx(cos(angle),sin(angle),8);q=map%grid(g,s)
 expanded_source(g,(m-1)*no+1:m*no,l)=phase*source(q,:,r-rep_start+1)*scale
 if(m==1)expanded_target(g,:,l)=phase*target(q,:,r-rep_start+1)
 enddo
 enddo
 enddo
 ierr=0
 end subroutine
end module
