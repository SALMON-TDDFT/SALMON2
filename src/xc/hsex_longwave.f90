! Insulating long-wavelength response. Optical vertices are R_c (dH/dk) u_v,
! where R_c=Q(H-e_v)^(-1)Q. No conduction eigenstates enter these solves.
module hsex_longwave
  implicit none
  private
  public :: longwave_band, longwave_cell
  abstract interface
    subroutine action(x,y,status)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine longwave_band(apply,occupied,energy,velocity,body,response,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    procedure(action) :: apply
    complex(8),intent(in) :: occupied(:,:),velocity(:,:),body(:,:)
    real(8),intent(in) :: energy
    complex(8),intent(out) :: response(:,:)
    integer,intent(out) :: status
    complex(8) :: vertices(size(occupied,1),3+size(body,2)),solution(size(occupied,1))
    integer :: axis,j
    status=1
    response=0d0
    if(size(velocity,2)/=3.or.size(velocity,1)/=size(occupied,1))return
    if(.not.screen_finite(velocity).or..not.screen_finite(body))return
    do axis=1,3
      call solve(velocity(:,axis),vertices(:,axis))
      if(status/=0)return
    enddo
    vertices(:,4:)=body
    do j=1,size(vertices,2)
      call solve(vertices(:,j),solution)
      if(status/=0)return
      response(:,j)=matmul(conjg(transpose(vertices)),solution)
    enddo
    if(.not.screen_finite(response))status=1
  contains
    subroutine solve(rhs,x)
      implicit none
      complex(8),intent(in) :: rhs(:)
      complex(8),intent(out) :: x(:)
      complex(8) :: b(size(rhs)),r(size(rhs)),p(size(rhs)),ap(size(rhs))
      real(8) :: rr,next_rr,curvature,tol
      integer :: iteration,ierr
      b=rhs-matmul(occupied,matmul(conjg(transpose(occupied)),rhs))
      x=0d0;r=b;p=r
      rr=real(dot_product(r,r),8)
      tol=max(1d-12,1d-10*sqrt(rr))
      status=1
      do iteration=1,2000
        if(sqrt(rr)<=tol)then
          call apply(x,ap,ierr)
          if(ierr/=0.or..not.screen_finite(ap))return
          ap=ap-energy*x
          ap=ap-matmul(occupied,matmul(conjg(transpose(occupied)),ap))
          r=b-ap
          rr=real(dot_product(r,r),8)
          if(sqrt(rr)<=tol)then
            status=0
            return
          endif
          p=r
        endif
        call apply(p,ap,ierr)
        if(ierr/=0.or..not.screen_finite(ap))return
        ap=ap-energy*p
        ap=ap-matmul(occupied,matmul(conjg(transpose(occupied)),ap))
        curvature=real(dot_product(p,ap),8)
        if(curvature<=0d0)return
        x=x+(rr/curvature)*p
        r=r-(rr/curvature)*ap
        next_rr=real(dot_product(r,r),8)
        if(.not.screen_finite([next_rr]))return
        p=r+(next_rr/rr)*p
        rr=next_rr
      enddo
    end subroutine
  end subroutine longwave_band

  subroutine longwave_cell(n,boxes,g,rec,vol,nk,chi,macro,point,average,status)
    use hsex_screen_math, only: screen_forms,screen_finite,screen_inverse
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    integer,intent(in) :: n(3),boxes,g(:,:),nk
    real(8),intent(in) :: rec(3,3),vol
    complex(8),intent(in) :: chi(:,:)
    real(8),intent(out) :: macro(3,3)
    complex(8),intent(out) :: point(:,:),average(:,:)
    integer,intent(out) :: status
    complex(8),allocatable :: dielectric(:,:),inverse(:,:),coupling(:,:),screened(:,:),forms(:,:),base(:,:),projected(:,:)
    complex(8),allocatable :: previous(:,:),previous_regular(:,:),regular(:,:),vector(:)
    real(8),allocatable :: bare(:),nodes(:),weights(:),radial(:),radial_weight(:)
    real(8) :: pi,factor,qc,mu,phi,direction(3),denom,weight,radius,error
    complex(8) :: tensor(3,3)
    real(8) :: coords(3,size(g,2)+1)
    real(8) :: rec_inverse(3,3),rec_work(3,3)
    integer :: rec_piv(3)
    integer :: nb,na,i,j,info,order,azimuth,polar,rad,iteration,nrad
    status=1;macro=0d0;point=0d0;average=0d0
    nb=size(g,2);na=boxes**3
    if(any(shape(chi)/=[nb+3,nb+3]).or.any(shape(point)/=[na,na]).or.any(shape(average)/=[na,na]))return
    if(any(n<1).or.boxes<1.or.nk<1.or.size(g,1)/=3)return
    if(any(mod(n,boxes)/=0))return
    if(.not.screen_finite(chi).or..not.screen_finite(rec).or..not.screen_finite([vol]).or.vol<=0d0)return
    status=2
    pi=acos(-1d0);factor=2d0/(nk*vol)
    rec_work=rec;rec_inverse=0d0
    do i=1,3
      rec_inverse(i,i)=1d0
    enddo
    call dgesv(3,3,rec_work,3,rec_piv,rec_inverse,3,info)
    if(info/=0)return
    allocate(bare(nb),dielectric(nb,nb),inverse(nb,nb),coupling(nb,3),screened(nb,3))
    do i=1,nb
      radius=sum(matmul(rec,dble(g(:,i)))**2)
      if(radius<=1d-14)return
      bare(i)=sqrt(4d0*pi/radius)
    enddo
    do j=1,nb
      do i=1,nb
        dielectric(i,j)=factor*bare(i)*chi(i+3,j+3)*bare(j)
      enddo
      dielectric(j,j)=dielectric(j,j)+1d0
    enddo
    inverse=0d0
    do i=1,nb
      inverse(i,i)=1d0
    enddo
    status=3
    if(nb>0)then
      call screen_inverse(dielectric,inverse,info)
      if(info/=0)return
    endif
    do i=1,nb
      coupling(i,:)=factor*sqrt(4d0*pi)*bare(i)*chi(i+3,:3)
    enddo
    screened=0d0
    tensor=4d0*pi*factor*chi(:3,:3)
    if(nb>0)then
      screened=matmul(inverse,coupling)
      tensor=tensor-matmul(conjg(transpose(coupling)),screened)
    endif
    do i=1,3
      tensor(i,i)=tensor(i,i)+1d0
    enddo
    status=4
    if(maxval(abs(aimag(tensor)))>1d-8.or.maxval(abs(tensor-transpose(tensor)))>1d-8)return
    macro=real(tensor,8)
    ! Positive definiteness is required in every direction, including ones
    ! absent from the spherical quadrature.
    status=5
    block
      real(8) :: eigen(3),matrix(3,3),workspace(16)
      matrix=macro
      call dsyev('N','L',3,matrix,3,eigen,workspace,size(workspace),info)
      if(info/=0.or.minval(eigen)<=0d0)return
    end block
    allocate(forms(na,nb+1),base(na,na),projected(na,3),previous(na,na),previous_regular(na,na),regular(na,na),vector(na))
    coords(:,1)=0d0;coords(:,2:)=dble(g)
    call screen_forms(n,boxes,coords,forms)
    do i=1,nb
      forms(:,i+1)=forms(:,i+1)*bare(i)
    enddo
    base=0d0;projected=0d0
    if(nb>0)then
      base=matmul(matmul(conjg(forms(:,2:)),inverse),transpose(forms(:,2:)))/vol
      projected=matmul(conjg(forms(:,2:)),screened)
    endif
    qc=(6d0*pi*pi/(vol*nk))**(1d0/3d0)
    previous=0d0;previous_regular=0d0
    status=6
    do iteration=1,6
      order=4*2**(iteration-1)
      nrad=min(order,32)
      allocate(nodes(order),weights(order),radial(nrad),radial_weight(nrad))
      call legendre(nodes,weights)
      call legendre(radial,radial_weight)
      radial=(1d0+radial)/2d0;radial_weight=radial_weight/2d0
      average=base;regular=base
      do polar=1,order
        mu=nodes(polar)
        do azimuth=1,2*order
          phi=pi*(azimuth-.5d0)/order
          direction=[sqrt(1d0-mu*mu)*cos(phi),sqrt(1d0-mu*mu)*sin(phi),mu]
          denom=dot_product(direction,matmul(macro,direction))
          weight=weights(polar)/(4d0*order*vol*denom)
          vector=matmul(projected,direction)
          do j=1,na
            regular(:,j)=regular(:,j)+weight*vector*conjg(vector(j))
          enddo
          do rad=1,nrad
            radius=qc*radial(rad)
            ! Only the analytic head/wings are continued inside this mini BZ.
            ! The finite-G body uses its q=0 limit. Box form factors retain
            ! their phases, so the integrable wing contribution is included.
            coords(:,1)=radius*matmul(rec_inverse,direction)
            do i=1,nb
              coords(:,i+1)=g(:,i)+coords(:,1)
            enddo
            call screen_forms(n,boxes,coords,forms)
            do i=1,nb
              forms(:,i+1)=forms(:,i+1)*bare(i)
            enddo
            vector=sqrt(4d0*pi)/radius*conjg(forms(:,1))
            if(nb>0)vector=vector-matmul(matmul(conjg(forms(:,2:)),screened),direction)
            do j=1,na
              average(:,j)=average(:,j)+weight*3d0*radial(rad)**2*radial_weight(rad)*vector*conjg(vector(j))
            enddo
          enddo
        enddo
      enddo
      deallocate(nodes,weights,radial,radial_weight)
      if(iteration>1)then
        error=max(maxval(abs(average-previous)),maxval(abs(regular-previous_regular)))/ &
          max(1d-14,maxval(abs(average)))
        write(error_unit,'(a,i0,a,es13.5)')'HSEX_W_LONGWAVE_QUADRATURE order=',order,' error=',error
        if(error<1d-7)exit
      endif
      previous=average
      previous_regular=regular
    enddo
    if(iteration>6)return
    if(.not.screen_finite(average).or..not.screen_finite(regular))return
    point=regular
    status=0
  end subroutine longwave_cell

  subroutine legendre(nodes,weights)
    implicit none
    real(8),intent(out) :: nodes(:),weights(:)
    real(8) :: z,p1,p2,p3,deriv,change
    integer :: i,j,l,n
    n=size(nodes)
    do i=1,n
      z=cos(acos(-1d0)*(i-.25d0)/(n+.5d0))
      do j=1,40
        p1=1d0;p2=0d0
        do l=1,n
          p3=p2;p2=p1
          p1=((2*l-1)*z*p2-(l-1)*p3)/l
        enddo
        deriv=n*(z*p1-p2)/(z*z-1d0)
        change=p1/deriv
        z=z-change
        if(abs(change)<1d-15)exit
      enddo
      nodes(i)=z;weights(i)=2d0/((1d0-z*z)*deriv*deriv)
    enddo
  end subroutine legendre
end module hsex_longwave
