! Uniform Cartesian E/H field provider; interpolation and molecular coupling stay in Fortran.
module salmonpy_grid_field
 use iso_c_binding
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 use salmonpy_near_field_coupling,only:near_coupling_map
 implicit none
 private
 public::grid_validate,grid_values,grid_potentials,grid_line_phase
 public::grid_magnetic_add,grid_phase_magnetic_add
 public::sp_grid_protocol,sp_grid_fields,sp_grid_potentials,sp_grid_pack,sp_grid_average
contains
 integer(c_int) function sp_grid_protocol() bind(C)
  sp_grid_protocol=1
 end function
 subroutine grid_metadata(counts,origins,spacing,status)
  integer,intent(in)::counts(3)
  real(8),intent(in)::origins(3,3,2),spacing(3)
  integer,intent(out)::status
  integer::i,j,k
  integer(8)::elements
  status=1
  if(any(counts<2))return
  elements=48_8
  do i=1,3
   if(elements>huge(elements)/counts(i))return
   elements=elements*counts(i)
  enddo
  do i=1,3
   if(.not.ieee_is_finite(spacing(i)))return
   if(spacing(i)<=0d0)return
   do j=1,3
   do k=1,2
    if(.not.ieee_is_finite(origins(i,j,k)))return
   enddo
   enddo
  enddo
  status=0
 end subroutine
 subroutine grid_validate(counts,origins,spacing,data,status)
  integer,intent(in)::counts(3)
  real(8),intent(in)::origins(3,3,2),spacing(3),data(3,counts(1),counts(2),counts(3),2)
  integer,intent(out)::status
  integer::i,j,k,d,which
  call grid_metadata(counts,origins,spacing,status)
  if(status/=0)return
  status=4
  do which=1,2
  do k=1,counts(3)
  do j=1,counts(2)
  do i=1,counts(1)
  do d=1,3
   if(.not.ieee_is_finite(data(d,i,j,k,which)))return
  enddo
  enddo
  enddo
  enddo
  enddo
  status=0
 end subroutine
 subroutine grid_values(counts,origins,spacing,data,point,E,H,status)
  integer,intent(in)::counts(3)
  real(8),intent(in)::origins(3,3,2),spacing(3),data(3,counts(1),counts(2),counts(3),2),point(3)
  real(8),intent(out)::E(3),H(3)
  integer,intent(out)::status
  real(8)::coordinate(3),fraction(3),w(0:1,3),v(3,2),item,tolerance
  integer::index(3),which,d,k,ix,iy,iz
  status=1
  do k=1,3
   if(.not.ieee_is_finite(point(k)))return
  enddo
  do which=1,2
  do d=1,3
   coordinate=(point-origins(:,d,which))/spacing
   do k=1,3
    if(.not.ieee_is_finite(coordinate(k)))return
    tolerance=32d0*epsilon(1d0)*max(1d0,abs(coordinate(k)))
    if(coordinate(k)<-tolerance.or.coordinate(k)>counts(k)-1+tolerance)return
    coordinate(k)=min(max(coordinate(k),0d0),real(counts(k)-1,8))
    index(k)=min(floor(coordinate(k))+1,counts(k)-1)
    fraction(k)=coordinate(k)-(index(k)-1)
    w(1,k)=fraction(k);w(0,k)=1d0-fraction(k)
   enddo
   item=0d0
   do iz=0,1
   do iy=0,1
   do ix=0,1
    if(.not.ieee_is_finite(data(d,index(1)+ix,index(2)+iy,index(3)+iz,which)))then
     status=4;return
    endif
    item=item+w(ix,1)*w(iy,2)*w(iz,3)*data(d,index(1)+ix,index(2)+iy,index(3)+iz,which)
   enddo
   enddo
   enddo
   if(.not.ieee_is_finite(item))then
    status=4;return
   endif
   v(d,which)=item
  enddo
  enddo
  E=v(:,1);H=v(:,2);status=0
 end subroutine
 subroutine grid_potentials(counts,origins,spacing,data,field_center,r,potential,Ac,status)
  integer,intent(in)::counts(3)
  real(8),intent(in)::origins(3,3,2),spacing(3),data(3,counts(1),counts(2),counts(3),2)
  real(8),intent(in)::field_center(3),r(3)
  real(8),intent(out)::potential,Ac(3)
  integer,intent(out)::status
  real(8),parameter::nodes(3)=[-.774596669241483377d0,0d0,.774596669241483377d0]
  real(8),parameter::weights(3)=[5d0/9d0,8d0/9d0,5d0/9d0]
  real(8)::next_face(3,3,2),increment(3),q,tol,left,right,t,w,E(3),H(3),ei(3),hi(3)
  integer::which,d,k,j,plane
  integer(8)::iterations
  call grid_values(counts,origins,spacing,data,field_center,E,H,status)
  if(status/=0)return
  call grid_values(counts,origins,spacing,data,field_center+r,E,H,status)
  if(status/=0)return
  next_face=huge(1d0);increment=0d0
  do k=1,3
   if(r(k)==0d0)cycle
   increment(k)=spacing(k)/abs(r(k))
   do which=1,2
   do d=1,3
    q=(field_center(k)-origins(k,d,which))/spacing(k)
    tol=32d0*epsilon(1d0)*max(1d0,abs(q))
    if(r(k)>0d0)then
     plane=floor(q+tol)+1
    else
     plane=ceiling(q-tol)-1
    endif
    next_face(k,d,which)=(origins(k,d,which)+plane*spacing(k)-field_center(k))/r(k)
   enddo
   enddo
  enddo
  left=0d0;ei=0d0;hi=0d0;iterations=0
  do while(left<1d0)
   right=min(1d0,minval(next_face))
   if(right<left)then
    status=4;return
   endif
   ! Within a cell: E(s) is cubic and s*H(s) quartic. Three Gauss points are exact.
   do j=1,3
    t=.5d0*((right-left)*nodes(j)+right+left);w=.5d0*(right-left)*weights(j)
    call grid_values(counts,origins,spacing,data,field_center+t*r,E,H,status)
    if(status/=0)return
    ei=ei+w*E;hi=hi+w*t*H
   enddo
   do which=1,2
   do d=1,3
   do k=1,3
    if(next_face(k,d,which)<=right+32d0*epsilon(1d0)) &
     next_face(k,d,which)=next_face(k,d,which)+increment(k)
   enddo
   enddo
   enddo
   left=right;iterations=iterations+1
   if(iterations>6_8*sum(int(counts,8))+1_8)then
    status=4;return
   endif
  enddo
  call near_coupling_map(r,ei,hi,potential,Ac)
  status=4
  if(.not.ieee_is_finite(potential))return
  do k=1,3
   if(.not.ieee_is_finite(Ac(k)))return
  enddo
  status=0
 end subroutine
 subroutine grid_line_phase(counts,origins,spacing,data,field_center,center,start,finish,phase,status)
  integer,intent(in)::counts(3)
  real(8),intent(in)::origins(3,3,2),spacing(3),data(3,counts(1),counts(2),counts(3),2)
  real(8),intent(in)::field_center(3),center(3),start(3),finish(3)
  real(8),intent(out)::phase
  integer,intent(out)::status
  real(8),parameter::nodes(8)=[.095012509837637440d0,.281603550779258913d0, &
    .458016777657227386d0,.617876244402643748d0,.755404408355003034d0, &
    .865631202387831744d0,.944575023073232576d0,.989400934991649933d0]
  real(8),parameter::weights(8)=[.189450610455068496d0,.182603415044923589d0, &
    .169156519395002538d0,.149595988816576732d0,.124628971255533872d0, &
    .095158511682492785d0,.062253523938647893d0,.027152459411754095d0]
  real(8)::r(3),potential,Ac(3),t,length_in_cells
  integer::segments,segment,k,side
  status=4;length_in_cells=2d0*maxval(abs(finish-start)/spacing)
  if(.not.ieee_is_finite(length_in_cells))return
  if(length_in_cells>real(huge(segments)-1,8))return
  segments=max(1,ceiling(length_in_cells));phase=0d0
  do segment=1,segments
  do k=1,8
  do side=-1,1,2
   t=(segment-.5d0+.5d0*side*nodes(k))/segments
   r=start+t*(finish-start)-center
   call grid_potentials(counts,origins,spacing,data,field_center,r,potential,Ac,status)
   if(status/=0)return
   phase=phase+.5d0*weights(k)*dot_product(Ac,finish-start)/segments
  enddo
  enddo
  enddo
  if(.not.ieee_is_finite(phase))status=4
 end subroutine
 subroutine grid_magnetic_add(counts,origins,spacing,position,r,q,M,status)
  ! Add -d(q.A/c)/dH / cell_volume, transpose of grid_potentials.
  integer,intent(in)::counts(3)
  real(8),intent(in)::origins(3,3),spacing(3),position(3),r(3),q(3)
  real(8),intent(inout)::M(3,counts(1),counts(2),counts(3))
  integer,intent(out)::status
  real(8),parameter::nodes(3)=[-.774596669241483377d0,0d0,.774596669241483377d0]
  real(8),parameter::weights(3)=[5d0/9d0,8d0/9d0,5d0/9d0]
  real(8)::next_face(3,3),increment(3),x,f,tol,left,right,t,w,point(3),basis(0:1,3),v(3),dv
  integer::d,k,j,plane,index(3),ix,iy,iz
  integer(8)::iterations
  status=4;dv=product(spacing)
  if(.not.ieee_is_finite(dv).or.dv<=0d0.or.any(counts<2))return
  do k=1,3
   if(.not.ieee_is_finite(r(k)).or..not.ieee_is_finite(q(k)))return
  enddo
  v=-[r(2)*q(3)-r(3)*q(2),r(3)*q(1)-r(1)*q(3),r(1)*q(2)-r(2)*q(1)]/(137.035999084d0*dv)
  next_face=huge(1d0);increment=0d0
  do k=1,3
   if(r(k)==0d0)cycle
   increment(k)=spacing(k)/abs(r(k))
   do d=1,3
    x=(position(k)-origins(k,d))/spacing(k)
    tol=32d0*epsilon(1d0)*max(1d0,abs(x))
    if(r(k)>0d0)then
     plane=floor(x+tol)+1
    else
     plane=ceiling(x-tol)-1
    endif
    next_face(k,d)=(origins(k,d)+plane*spacing(k)-position(k))/r(k)
   enddo
  enddo
  left=0d0;iterations=0
  do while(left<1d0)
   right=min(1d0,minval(next_face));if(right<left)return
   do j=1,3
    t=.5d0*((right-left)*nodes(j)+right+left);w=.5d0*(right-left)*weights(j)*t
    point=position+t*r
    do d=1,3
     do k=1,3
      x=(point(k)-origins(k,d))/spacing(k)
      if(.not.ieee_is_finite(x))return
      tol=32d0*epsilon(1d0)*max(1d0,abs(x))
      if(x< -tol.or.x>counts(k)-1+tol)return
      x=min(max(x,0d0),real(counts(k)-1,8));index(k)=min(floor(x)+1,counts(k)-1)
      f=x-(index(k)-1);basis(1,k)=f;basis(0,k)=1d0-f
     enddo
     do iz=0,1
     do iy=0,1
     do ix=0,1
      M(d,index(1)+ix,index(2)+iy,index(3)+iz)= &
       M(d,index(1)+ix,index(2)+iy,index(3)+iz)+v(d)*w*basis(ix,1)*basis(iy,2)*basis(iz,3)
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
  status=0
 end subroutine
 subroutine grid_phase_magnetic_add(counts,origins,spacing,field_center,center,start,finish,gradient,M,status)
  integer,intent(in)::counts(3)
  real(8),intent(in)::origins(3,3),spacing(3)
  real(8),intent(in)::field_center(3),center(3),start(3),finish(3)
  real(8),intent(in)::gradient
  real(8),intent(inout)::M(3,counts(1),counts(2),counts(3))
  integer,intent(out)::status
  real(8),parameter::nodes(8)=[.095012509837637440d0,.281603550779258913d0, &
    .458016777657227386d0,.617876244402643748d0,.755404408355003034d0, &
    .865631202387831744d0,.944575023073232576d0,.989400934991649933d0]
  real(8),parameter::weights(8)=[.189450610455068496d0,.182603415044923589d0, &
    .169156519395002538d0,.149595988816576732d0,.124628971255533872d0, &
    .095158511682492785d0,.062253523938647893d0,.027152459411754095d0]
  real(8)::r(3),q(3),t,length_in_cells
  integer::segments,segment,k,side
  status=4;length_in_cells=2d0*maxval(abs(finish-start)/spacing)
  if(.not.ieee_is_finite(length_in_cells))return
  if(length_in_cells>real(huge(segments)-1,8))return
  segments=max(1,ceiling(length_in_cells))
  do segment=1,segments
  do k=1,8
  do side=-1,1,2
   t=(segment-.5d0+.5d0*side*nodes(k))/segments
   r=start+t*(finish-start)-center
   q=gradient*.5d0*weights(k)*(finish-start)/segments
   call grid_magnetic_add(counts,origins,spacing,field_center,r,q,M,status)
   if(status/=0)return
  enddo
  enddo
  enddo

 end subroutine
 integer(c_int) function sp_grid_fields(n,counts,origins,spacing,field_center,center,points,data,E,H) &
  bind(C) result(status)
  integer(c_int),value::n
  integer(c_int),intent(in)::counts(3)
  real(c_double),intent(in)::origins(3,3,2),spacing(3),field_center(3),center(3),points(3,n)
  real(c_double),intent(in)::data(3,counts(1),counts(2),counts(3),2)
  real(c_double),intent(inout)::E(3,n),H(3,n)
  real(8),allocatable::buffer(:,:,:)
  integer::i,code,allocation_status
  status=1
  if(n<1)return
  call grid_validate(counts,origins,spacing,data,code)
  if(code/=0)then
   status=code;return
  endif
  allocate(buffer(3,n,2),stat=allocation_status);status=2
  if(allocation_status/=0)return
  do i=1,n
   call grid_values(counts,origins,spacing,data,field_center+points(:,i)-center, &
                    buffer(:,i,1),buffer(:,i,2),code)
   if(code/=0)then
    status=code;return
   endif
  enddo
  E=buffer(:,:,1);H=buffer(:,:,2);status=0
 end function
 integer(c_int) function sp_grid_potentials(n,counts,origins,spacing,field_center,center,points,data,V,Ac) &
  bind(C) result(status)
  integer(c_int),value::n
  integer(c_int),intent(in)::counts(3)
  real(c_double),intent(in)::origins(3,3,2),spacing(3),field_center(3),center(3),points(3,n)
  real(c_double),intent(in)::data(3,counts(1),counts(2),counts(3),2)
  real(c_double),intent(inout)::V(n),Ac(3,n)
  real(8),allocatable::buffer(:,:)
  integer::i,code,allocation_status
  status=1
  if(n<1)return
  call grid_validate(counts,origins,spacing,data,code)
  if(code/=0)then
   status=code;return
  endif
  allocate(buffer(4,n),stat=allocation_status);status=2
  if(allocation_status/=0)return
  do i=1,n
   call grid_potentials(counts,origins,spacing,data,field_center,points(:,i)-center, &
                        buffer(1,i),buffer(2:4,i),code)
   if(code/=0)then
    status=code;return
   endif
  enddo
  V=buffer(1,:);Ac=buffer(2:4,:);status=0
 end function
 integer(c_int) function sp_grid_pack(n,counts,origin,spacing,points,fields,data) bind(C) result(status)
  integer(c_int),value::n
  integer(c_int),intent(in)::counts(3)
  real(c_double),intent(in)::origin(3),spacing(3),points(3,n),fields(3,n,2)
  real(c_double),intent(inout)::data(3,counts(1),counts(2),counts(3),2)
  real(8),allocatable::buffer(:,:,:,:,:)
  logical,allocatable::assigned(:,:,:)
  real(8)::q(3),origins(3,3,2),tol
  integer::i,k,d,which,index(3),allocation_status
  status=1
  do which=1,2
  do d=1,3
   origins(:,d,which)=origin
  enddo
  enddo
  call grid_metadata(counts,origins,spacing,k)
  if(k/=0)return
  if(int(n,8)/=product(int(counts,8)))return
  allocate(buffer(3,counts(1),counts(2),counts(3),2), &
           assigned(counts(1),counts(2),counts(3)),stat=allocation_status);status=2
  if(allocation_status/=0)return
  assigned=.false.;status=1
  do i=1,n
   do k=1,3
    if(.not.ieee_is_finite(points(k,i)))return
   enddo
   q=(points(:,i)-origin)/spacing
   do k=1,3
    if(.not.ieee_is_finite(q(k)))return
    tol=64d0*epsilon(1d0)*max(1d0,abs(q(k)))
    if(q(k)<-tol.or.q(k)>counts(k)-1+tol)return
    index(k)=nint(q(k))+1
    if(abs(q(k)-(index(k)-1))>tol)return
   enddo
   if(assigned(index(1),index(2),index(3)))return
   assigned(index(1),index(2),index(3))=.true.
   do which=1,2
   do d=1,3
    if(.not.ieee_is_finite(fields(d,i,which)))return
    buffer(d,index(1),index(2),index(3),which)=fields(d,i,which)
   enddo
   enddo
  enddo
  data=buffer;status=0
 end function
 integer(c_int) function sp_grid_average(n,left,right,centered) bind(C) result(status)
  integer(c_int),value::n
  real(c_double),intent(in)::left(n),right(n)
  real(c_double),intent(inout)::centered(n)
  real(8),allocatable::buffer(:)
  integer::i,allocation_status
  status=1
  if(n<1)return
  allocate(buffer(n),stat=allocation_status);status=2
  if(allocation_status/=0)return
  status=4
  do i=1,n
   if(.not.ieee_is_finite(left(i)).or..not.ieee_is_finite(right(i)))return
   buffer(i)=.5d0*left(i)+.5d0*right(i)
   if(.not.ieee_is_finite(buffer(i)))return
  enddo
  centered=buffer;status=0
 end function
end module
