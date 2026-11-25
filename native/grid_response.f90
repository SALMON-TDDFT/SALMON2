! Transpose of the electric line-integral coupling, on E-component grids.
! Signed charge density is supplied by the caller, including its quadrature volumes.
! This stateless operator is not a live TDDFT/FDTD coupling driver.
module salmonpy_grid_response
 use iso_c_binding
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::sp_grid_electric_response,sp_grid_difference,sp_grid_magnetization_current,sp_grid_current_sum
contains
 subroutine stencil(counts,origin,spacing,point,index,weights,status)
  integer,intent(in)::counts(3)
  real(8),intent(in)::origin(3),spacing(3),point(3)
  integer,intent(out)::index(3),status
  real(8),intent(out)::weights(0:1,3)
  real(8)::q,f,tolerance
  integer::k
  status=1
  do k=1,3
   q=(point(k)-origin(k))/spacing(k)
   if(.not.ieee_is_finite(q))return
   tolerance=32d0*epsilon(1d0)*max(1d0,abs(q))
   if(q< -tolerance.or.q>counts(k)-1+tolerance)return
   q=min(max(q,0d0),real(counts(k)-1,8))
   index(k)=min(floor(q)+1,counts(k)-1);f=q-(index(k)-1)
   weights(1,k)=f;weights(0,k)=1d0-f
  enddo
  status=0
 end subroutine
 integer(c_int) function sp_grid_electric_response(mode,n,counts,origins,spacing,position,center, &
                              points,volumes,old,new,dt,result) bind(C) result(status)
  integer(c_int),value::mode,n
  integer(c_int),intent(in)::counts(3)
  real(c_double),intent(in)::origins(3,3),spacing(3),position(3),center(3),points(3,n),volumes(n),old(n),new(n)
  real(c_double),value::dt
  real(c_double),intent(inout)::result(3,counts(1),counts(2),counts(3))
  real(8),parameter::nodes(3)=[-.774596669241483377d0,0d0,.774596669241483377d0]
  real(8),parameter::gauss(3)=[5d0/9d0,8d0/9d0,5d0/9d0]
  real(8),allocatable::buffer(:,:,:,:)
  real(8)::r(3),next_face(3,3),increment(3),q,tol,left,right,t,w,basis(0:1,3),charge,cell_volume,value
  integer::a,d,k,j,plane,index(3),code,ios,ix,iy,iz
  integer(8)::elements,iterations
  status=1
  if(n<1.or.any(counts<2).or.mode<0.or.mode>1)return
  if(.not.ieee_is_finite(dt).or.dt<=0d0)return
  elements=24_8;cell_volume=1d0
  do k=1,3
   if(elements>huge(elements)/counts(k))return
   elements=elements*counts(k)
   if(.not.ieee_is_finite(spacing(k)).or.spacing(k)<=0d0)return
   if(.not.ieee_is_finite(position(k)).or..not.ieee_is_finite(center(k)))return
   do d=1,3
    if(.not.ieee_is_finite(origins(k,d)))return
   enddo
   cell_volume=cell_volume*spacing(k)
  enddo
  if(.not.ieee_is_finite(cell_volume).or.cell_volume<=0d0)return
  do a=1,n
   if(.not.ieee_is_finite(volumes(a)).or.volumes(a)<0d0)return
   if(.not.ieee_is_finite(old(a)).or..not.ieee_is_finite(new(a)))return
   do k=1,3
    if(.not.ieee_is_finite(points(k,a)))return
   enddo
  enddo
  status=2
  allocate(buffer(3,counts(1),counts(2),counts(3)),stat=ios)
  if(ios/=0)return
  buffer=0d0
  do a=1,n
   r=points(:,a)-center
   charge=new(a)*volumes(a)
   if(mode==1)charge=((new(a)-old(a))/dt)*volumes(a)
   status=4
   if(.not.ieee_is_finite(charge))return
   do k=1,3
    if(.not.ieee_is_finite(r(k)))return
   enddo
   do d=1,3
    call stencil(counts,origins(:,d),spacing,position,index,basis,code)
    if(code/=0)then
     status=code;return
    endif
    call stencil(counts,origins(:,d),spacing,position+r,index,basis,code)
    if(code/=0)then
     status=code;return
    endif
   enddo
   ! Split at all E-component faces. Three-point Gauss integrates each
   ! trilinear basis along the ray exactly, including paths across cells.
   next_face=huge(1d0);increment=0d0
   do k=1,3
    if(r(k)==0d0)cycle
    increment(k)=spacing(k)/abs(r(k))
    do d=1,3
     q=(position(k)-origins(k,d))/spacing(k)
     tol=32d0*epsilon(1d0)*max(1d0,abs(q))
     if(r(k)>0d0)then
      plane=floor(q+tol)+1
     else
      plane=ceiling(q-tol)-1
     endif
     next_face(k,d)=(origins(k,d)+plane*spacing(k)-position(k))/r(k)
    enddo
   enddo
   left=0d0;iterations=0
   do while(left<1d0)
    right=min(1d0,minval(next_face));status=4
    if(right<left)return
    do j=1,3
     t=.5d0*((right-left)*nodes(j)+right+left);w=.5d0*(right-left)*gauss(j)
     do d=1,3
      call stencil(counts,origins(:,d),spacing,position+t*r,index,basis,code)
      if(code/=0)then
       status=code;return
      endif
      do iz=0,1
      do iy=0,1
      do ix=0,1
       value=charge*r(d)*w*basis(ix,1)*basis(iy,2)*basis(iz,3)/cell_volume
       buffer(d,index(1)+ix,index(2)+iy,index(3)+iz)=buffer(d,index(1)+ix,index(2)+iy,index(3)+iz)+value
      enddo
      enddo
      enddo
     enddo
    enddo
    do d=1,3
    do k=1,3
     if(next_face(k,d)<=right+32d0*epsilon(1d0))next_face(k,d)=next_face(k,d)+increment(k)
    enddo
    enddo
    left=right;iterations=iterations+1
    if(iterations>3_8*sum(int(counts,8))+1_8)return
   enddo
  enddo
  status=4
  do iz=1,counts(3)
  do iy=1,counts(2)
  do ix=1,counts(1)
  do d=1,3
   if(.not.ieee_is_finite(buffer(d,ix,iy,iz)))return
  enddo
  enddo
  enddo
  enddo
  result=buffer;status=0
 end function
 integer(c_int) function sp_grid_difference(n,total,own,result) bind(C) result(status)
  integer(c_int),value::n
  real(c_double),intent(in)::total(n),own(n)
  real(c_double),intent(inout)::result(n)
  real(8),allocatable::buffer(:)
  integer::i,ios
  status=1
  if(n<1)return
  do i=1,n
   if(.not.ieee_is_finite(total(i)).or..not.ieee_is_finite(own(i)))return
  enddo
  status=2;allocate(buffer(n),stat=ios);if(ios/=0)return
  status=4
  do i=1,n
   buffer(i)=total(i)-own(i)
   if(.not.ieee_is_finite(buffer(i)))return
  enddo
  result=buffer;status=0
 end function
 integer(c_int) function sp_grid_current_sum(n,electric,magnetic,total) bind(C) result(status)
  integer(c_int),value::n
  real(c_double),intent(in)::electric(n),magnetic(n)
  real(c_double),intent(inout)::total(n)
  real(8),allocatable::buffer(:)
  integer::i,ios
  status=1;if(n<1)return
  do i=1,n
   if(.not.ieee_is_finite(electric(i)).or..not.ieee_is_finite(magnetic(i)))return
  enddo
  status=2;allocate(buffer(n),stat=ios);if(ios/=0)return
  status=4
  do i=1,n
   buffer(i)=electric(i)+magnetic(i)
   if(.not.ieee_is_finite(buffer(i)))return
  enddo
  total=buffer;status=0
 end function
 integer(c_int) function sp_grid_magnetization_current(counts,spacing,M,current) bind(C) result(status)
  ! Vacuum Gaussian/atomic units: j_M=c curl M. H components at Yee faces,
  ! E components at Yee edges. Zero extension includes the patch surface current.
  integer(c_int),intent(in)::counts(3)
  real(c_double),intent(in)::spacing(3),M(3,counts(1),counts(2),counts(3))
  real(c_double),intent(inout)::current(3,counts(1),counts(2),counts(3))
  real(8),allocatable::buffer(:,:,:,:)
  real(8)::back(3,3),v(3)
  integer::i,j,k,d,a,ios,idx(3),lower(3)
  status=1
  if(any(counts<2))return
  do a=1,3
   if(.not.ieee_is_finite(spacing(a)).or.spacing(a)<=0d0)return
  enddo
  do k=1,counts(3)
  do j=1,counts(2)
  do i=1,counts(1)
  do d=1,3
   if(.not.ieee_is_finite(M(d,i,j,k)))return
  enddo
  enddo
  enddo
  enddo
  ! The positive-face support must lie inside the patch; truncating a curl loses surface current.
  if(any(M(:,counts(1),:,:)/=0d0).or.any(M(:,:,counts(2),:)/=0d0).or. &
     any(M(:,:,:,counts(3))/=0d0))return
  status=2;allocate(buffer(3,counts(1),counts(2),counts(3)),stat=ios);if(ios/=0)return
  do k=1,counts(3)
  do j=1,counts(2)
  do i=1,counts(1)
   idx=[i,j,k];v=M(:,i,j,k);back=0d0
   do a=1,3
    lower=idx;lower(a)=lower(a)-1
    if(lower(a)>=1)back(:,a)=M(:,lower(1),lower(2),lower(3))
    back(:,a)=(v-back(:,a))/spacing(a)
   enddo
   buffer(:,i,j,k)=137.035999084d0*[back(3,2)-back(2,3),back(1,3)-back(3,1),back(2,1)-back(1,2)]
  enddo
  enddo
  enddo
  status=4
  do k=1,counts(3)
  do j=1,counts(2)
  do i=1,counts(1)
  do d=1,3
   if(.not.ieee_is_finite(buffer(d,i,j,k)))return
  enddo
  enddo
  enddo
  enddo
  current=buffer;status=0
 end function
end module
