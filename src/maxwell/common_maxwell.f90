!
!  Copyright 2019-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!-----------------------------------------------------------------------------------------
module common_maxwell
  use structures, only: s_vector, s_scalar
  implicit none
  
  type ls_fdtd_work
    !weyl
    type(s_vector) :: vec_Ac_tmp            ! Tempolary variable used in weyl_calc
    type(s_scalar) :: edensity_emfield      ! Field energy: (E^2 + H^2)
    type(s_scalar) :: edensity_absorb       ! Loss energy: integrate(E*J)
  end type ls_fdtd_work
  
  public :: set_coo_em
  public :: find_point_em
  public :: make_shape_em
  public :: input_i3d_em
  public :: output_i3d_em
  
contains
  
  !=========================================================================================
  != set coordinate system in fdtd =========================================================
  subroutine set_coo_em(Nd,ioe,ista,iend,hgs,coo,peri)
    implicit none
    integer,intent(in)      :: Nd
    integer,intent(in)      :: ioe(3),ista(3),iend(3)
    real(8),intent(in)      :: hgs(3)
    real(8),intent(out)     :: coo(minval(ista(:))-Nd:maxval(iend(:))+Nd,3)
    character(1),intent(in) :: peri
    integer :: ii,ij
    
    do ii=1,3
      select case(peri)
      case('n')
        select case(ioe(ii))
        case(1)
          do ij=ista(ii)-Nd,iend(ii)+Nd
            coo(ij,ii)=dble(ij)*hgs(ii)
          end do
        case(2)
          do ij=ista(ii)-Nd,iend(ii)+Nd
            coo(ij,ii)=(dble(ij)-0.5d0)*hgs(ii)
          end do
        end select
      case('y')
        do ij=ista(ii)-Nd,iend(ii)+Nd
          coo(ij,ii)=dble(ij-1)*hgs(ii)
        end do
      end select
    end do
    
    return
  end subroutine set_coo_em
  
  !=========================================================================================
  != find point, line, and plane in fdtd ===================================================
  subroutine find_point_em(rloc,id,ipo,ili,ipl,ista,iend,icoo_sta,icoo_end,coo)
    use parallelization, only: nproc_id_global,nproc_size_global,nproc_group_global
    use communication,   only: comm_summation
    implicit none
    real(8),intent(in)    :: rloc(3)
    integer,intent(inout) :: id(3)
    integer,intent(out)   :: ipo
    integer,intent(out)   :: ili(3),ipl(3)
    integer,intent(in)    :: ista(3),iend(3)
    integer,intent(in)    :: icoo_sta,icoo_end
    real(8),intent(in)    :: coo(icoo_sta:icoo_end,3)
    integer               :: ii,ix,iy,iz,ipe_tmp,i1,i1s,i2,i2s
    integer               :: id_tmp(3),id_tmp2(3)
    real(8)               :: err(0:nproc_size_global-1),err2(0:nproc_size_global-1)
    real(8)               :: err_tmp
    
    !set initial value
    err(:)              =0.0d0
    err(nproc_id_global)=1.0d9
    err2(:)             =0.0d0
    id_tmp(:)           =0
    id_tmp2(:)          =0
    
    !find observation point in each PE
    do iz=ista(3),iend(3)
    do iy=ista(2),iend(2)
    do ix=ista(1),iend(1)
      err_tmp=sqrt( (coo(ix,1)-rloc(1))**2.0d0 &
                   +(coo(iy,2)-rloc(2))**2.0d0 &
                   +(coo(iz,3)-rloc(3))**2.0d0 )
      if(err(nproc_id_global)>err_tmp) then
        err(nproc_id_global)=err_tmp
        id_tmp(1)=ix; id_tmp(2)=iy; id_tmp(3)=iz;
      end if
    end do
    end do
    end do
    call comm_summation(err,err2,nproc_size_global,nproc_group_global)
    
    !determine and share id + determine pe including the point
    ipe_tmp=-1; err_tmp=1.0d9;
    do ii=0,nproc_size_global-1
      if(err_tmp>err2(ii)) then
        err_tmp=err2(ii)
        ipe_tmp=ii
      end if
    end do
    if(nproc_id_global==ipe_tmp) then
      ipo=1;
    else
      ipo=0; id_tmp(:)=0;
    end if
    call comm_summation(id_tmp,id_tmp2,3,nproc_group_global)
    id(:)=id_tmp2(:)
    
    !determine pe including the line
    do ii=1,3
      if(ii==1) then     !x-line(searching at yz-plane)
        i1s=3; i2s=2;
      elseif(ii==2) then !x-line(searching at xz-plane)
        i1s=3; i2s=1;
      elseif(ii==3) then !z-line(searching at xy-plane)
        i1s=2; i2s=1;
      end if
      do i2=ista(i2s),iend(i2s)
      do i1=ista(i1s),iend(i1s)
        if( (i1==id(i1s)).and.(i2==id(i2s)) ) ili(ii)=1
      end do
      end do
    end do
    
    !determine pe including the plane
    do ii=1,3
      if(ii==1) then     !xy-plane(searching at z-line)
        i1s=3;
      elseif(ii==2) then !yz-plane(searching at x-line)
        i1s=1;
      elseif(ii==3) then !xz-plane(searching at y-line)
        i1s=2;
      end if
      do i1=ista(i1s),iend(i1s)
        if(i1==id(i1s)) ipl(ii)=1
      end do
    end do
    
    return
  end subroutine find_point_em
  
  !=========================================================================================
  != make shape ============================================================================
  subroutine make_shape_em(mg,lg,rlsize,imat,Nd,coo)
    use salmon_global,   only: n_s,typ_s,id_s,inf_s,ori_s,rot_s,&
                               yn_copy_x,yn_copy_y,yn_copy_z,rot_type
    use structures,      only: s_rgrid
    use math_constants,  only: pi
    implicit none
    type(s_rgrid), intent(in)  :: mg,lg
    real(8),       intent(in)  :: rlsize(3)
    integer,       intent(out) :: imat(mg%is_array(1):mg%ie_array(1),&
                                       mg%is_array(2):mg%ie_array(2),&
                                       mg%is_array(3):mg%ie_array(3))
    integer,       intent(in)  :: Nd
    real(8),       intent(in)  :: coo(minval(lg%is(:))-Nd:maxval(lg%ie(:))+Nd,3)
    integer             :: icopy_num(3)
    real(8)             :: rot_s_d(1000,3)
    real(8),allocatable :: rmove_x(:,:),rmove_y(:,:),rmove_z(:,:)
    integer :: ii,ix,iy,iz,ip_x,ip_y,ip_z,il,l_max
    real(8) :: x,y,z,x_o,y_o,z_o,x_tmp,y_tmp,z_tmp,adj_err,cal_tmp
    
    !convert from degree to radian or keep radian
    if(trim(rot_type)=='degree') then
      rot_s_d(:,:) = (rot_s(:,:)/360.0d0) * (2.0d0*pi)
    else
      rot_s_d(:,:) = rot_s(:,:)
    end if
    
    !make icopy_num and l_max
    if    (yn_copy_x=='y') then
      icopy_num(1) = 3
    elseif(yn_copy_x=='n') then
      icopy_num(1) = 1     
    end if
    if    (yn_copy_y=='y') then
      icopy_num(2) = 3
    elseif(yn_copy_y=='n') then
      icopy_num(2) = 1     
    end if
    if    (yn_copy_z=='y') then
      icopy_num(3) = 3
    elseif(yn_copy_z=='n') then
      icopy_num(3) = 1     
    end if
    l_max = icopy_num(1)*icopy_num(2)*icopy_num(3)
    
    !make move matrix
    allocate(rmove_x(n_s,l_max),rmove_y(n_s,l_max),rmove_z(n_s,l_max))
    do ii=1,n_s
      il = 1
      do iz=1,icopy_num(3)
      do iy=1,icopy_num(2)
      do ix=1,icopy_num(1)
        ip_x = ix-1; ip_y = iy-1; ip_z = iz-1;
        if(yn_copy_x=='y') ip_x = ip_x-1;
        if(yn_copy_y=='y') ip_y = ip_y-1;
        if(yn_copy_z=='y') ip_z = ip_z-1;
        x_tmp = rlsize(1) * dble(ip_x)
        y_tmp = rlsize(2) * dble(ip_y)
        z_tmp = rlsize(3) * dble(ip_z)
        call rotate_x(x_tmp,y_tmp,z_tmp,rot_s_d(ii,:),rmove_x(ii,il))
        call rotate_y(x_tmp,y_tmp,z_tmp,rot_s_d(ii,:),rmove_y(ii,il))
        call rotate_z(x_tmp,y_tmp,z_tmp,rot_s_d(ii,:),rmove_z(ii,il))
        il = il+1
      end do
      end do
      end do
    end do
    
    !set adjust parameter
    adj_err = 1.0d-6
    
    !make shape
    cal_tmp = 0.0d0
    do ii=1,n_s
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        !move origin
        x_tmp = coo(ix,1) - ori_s(ii,1)
        y_tmp = coo(iy,2) - ori_s(ii,2)
        z_tmp = coo(iz,3) - ori_s(ii,3)
        
        !rotation
        call rotate_x(x_tmp,y_tmp,z_tmp,rot_s_d(ii,:),x_o)
        call rotate_y(x_tmp,y_tmp,z_tmp,rot_s_d(ii,:),y_o)
        call rotate_z(x_tmp,y_tmp,z_tmp,rot_s_d(ii,:),z_o)
        
        !copy loop
        do il=1,l_max
          !determine point
          x = x_o + rmove_x(ii,il)
          y = y_o + rmove_y(ii,il)
          z = z_o + rmove_z(ii,il)
          
          !determine shape
          if    (trim(typ_s(ii))=='ellipsoid')            then
            cal_tmp = (x/(inf_s(ii,1)/2.0d0))**2.0d0 + (y/(inf_s(ii,2)/2.0d0))**2.0d0 + (z/(inf_s(ii,3)/2.0d0))**2.0d0
            if(cal_tmp<=1.0d0) imat(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='half-ellipsoid')       then
            cal_tmp = (x/(inf_s(ii,1)/2.0d0))**2.0d0 + (y/(inf_s(ii,2)/2.0d0))**2.0d0 + (z/(inf_s(ii,3)      ))**2.0d0
            if((cal_tmp<=1.0d0).and.(z>=-adj_err)) imat(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='elliptic-cylinder')    then
            cal_tmp = (x/(inf_s(ii,1)/2.0d0))**2.0d0 + (y/(inf_s(ii,2)/2.0d0))**2.0d0
            if((cal_tmp<=1.0d0).and.(z>=-inf_s(ii,3)/2.0d0).and.(z<=inf_s(ii,3)/2.0d0)) imat(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='triangular-cylinder')  then
            if( (x>= -inf_s(ii,1)/2.0d0).and.(x<=inf_s(ii,1)/2.0d0).and.                  &
                (y>= -inf_s(ii,2)/3.0d0).and.                                             &
                (y<=( inf_s(ii,2)/(inf_s(ii,1)/2.0d0)*x + inf_s(ii,2)*2.0d0/3.0d0 )).and. &
                (y<=(-inf_s(ii,2)/(inf_s(ii,1)/2.0d0)*x + inf_s(ii,2)*2.0d0/3.0d0 )).and. &
                (z>= -inf_s(ii,3)/2.0d0).and.(z<=inf_s(ii,3)/2.0d0) )                     &
              imat(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='rectangular-cylinder') then
            if( (x>=-inf_s(ii,1)/2.0d0).and.(x<=inf_s(ii,1)/2.0d0).and. &
                (y>=-inf_s(ii,2)/2.0d0).and.(y<=inf_s(ii,2)/2.0d0).and. &
                (z>=-inf_s(ii,3)/2.0d0).and.(z<=inf_s(ii,3)/2.0d0) )    &
              imat(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='elliptic-cone')        then
            if(inf_s(ii,3)-z/=0.0d0) then 
              cal_tmp= (x/( inf_s(ii,1)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3) ))**2.0d0 &
                     + (y/( inf_s(ii,2)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3) ))**2.0d0
            else
              cal_tmp=10
            end if
            if( (cal_tmp<=1.0d0).and.(z>=-adj_err).and.(z<=inf_s(ii,3)) ) imat(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='triangular-cone')      then
            if( (x>= -inf_s(ii,1)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3))        .and. &
                (x<=  inf_s(ii,1)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3))        .and. &
                (y>= -inf_s(ii,2)/3.0d0*(inf_s(ii,3)-z)/inf_s(ii,3))        .and. &
                (y<=( inf_s(ii,2)/(inf_s(ii,1)/2.0d0)*x                           &
                     +inf_s(ii,2)*2.0d0/3.0d0*(inf_s(ii,3)-z)/inf_s(ii,3) )).and. &
                (y<=(-inf_s(ii,2)/(inf_s(ii,1)/2.0d0)*x                           &
                     +inf_s(ii,2)*2.0d0/3.0d0*(inf_s(ii,3)-z)/inf_s(ii,3) )).and. &
                (z>=-adj_err).and.(z<=inf_s(ii,3)) )                               &
              imat(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='rectangular-cone')     then
            if( (x>= -inf_s(ii,1)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3)).and. &
                (x<=  inf_s(ii,1)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3)).and. &
                (y>= -inf_s(ii,2)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3)).and. &
                (y<=  inf_s(ii,2)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3)).and. &
                (z>=-adj_err).and.(z<=inf_s(ii,3)) )                      &
              imat(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='elliptic-ring')        then
            cal_tmp=(x/(inf_s(ii,1)/2.0d0))**2.0d0 + (y/(inf_s(ii,2)/2.0d0))**2.0d0
            if((cal_tmp<=1.0d0).and.(z>=-inf_s(ii,3)/2.0d0).and.(z<=inf_s(ii,3)/2.0d0)) then
              cal_tmp=(x/(inf_s(ii,4)/2.0d0))**2.0d0 + (y/(inf_s(ii,5)/2.0d0))**2.0d0
              if(cal_tmp>=1.0d0) imat(ix,iy,iz) = id_s(ii)
            end if
          end if
        end do
      end do
      end do
      end do
    end do
    
    return
  contains
    
    !+ CONTAINED IN make_shape_em ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ rotation for x ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine rotate_x(x_in,y_in,z_in,rot_in,x_out)
      implicit none
      real(8),intent(in)  :: x_in,y_in,z_in
      real(8),intent(in)  :: rot_in(3)
      real(8),intent(out) :: x_out
      
      x_out = cos(rot_in(3))*( cos(rot_in(2))*x_in - sin(rot_in(2))*(cos(rot_in(1))*z_in - sin(rot_in(1))*y_in) ) &
             +sin(rot_in(3))*( sin(rot_in(1))*z_in + cos(rot_in(1))*y_in )
      
      return
    end subroutine rotate_x
    
    !+ CONTAINED IN make_shape_em ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ rotation for y ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine rotate_y(x_in,y_in,z_in,rot_in,y_out)
      implicit none
      real(8),intent(in)  :: x_in,y_in,z_in
      real(8),intent(in)  :: rot_in(3)
      real(8),intent(out) :: y_out
      
      y_out = -sin(rot_in(3))*( cos(rot_in(2))*x_in - sin(rot_in(2))*( cos(rot_in(1))*z_in - sin(rot_in(1))*y_in) ) &
              +cos(rot_in(3))*( sin(rot_in(1))*z_in + cos(rot_in(1))*y_in )
      
      return
    end subroutine rotate_y
    
    !+ CONTAINED IN make_shape_em ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ rotation for z ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine rotate_z(x_in,y_in,z_in,rot_in,z_out)
      implicit none
      real(8),intent(in)  :: x_in,y_in,z_in
      real(8),intent(in)  :: rot_in(3)
      real(8),intent(out) :: z_out
      
      z_out = sin(rot_in(2))*x_in + cos(rot_in(2))*( cos(rot_in(1))*z_in - sin(rot_in(1))*y_in )
      
      return
    end subroutine rotate_z
    
  end subroutine make_shape_em
  
  !=========================================================================================
  != input i3d file ========================================================================
  subroutine input_i3d_em(file_name,ifn,ng_is,ng_ie,lg_is,lg_ie,Nd,imat,format)
    use parallelization, only: nproc_id_global
    use communication,   only: comm_is_root
    implicit none
    character(256),intent(in)  :: file_name
    integer,       intent(in)  :: ifn,Nd
    integer,       intent(in)  :: ng_is(3),ng_ie(3),lg_is(3),lg_ie(3)
    integer,       intent(out) :: imat(ng_is(1)-Nd:ng_ie(1)+Nd,&
                                       ng_is(2)-Nd:ng_ie(2)+Nd,&
                                       ng_is(3)-Nd:ng_ie(3)+Nd)
    character(2),  intent(in)  :: format
    real(8),allocatable        :: rtmp1d(:)
    integer                    :: inum(3),inum_check(3)
    integer(8)                 :: ii
    integer                    :: ij,ix,iy,iz,iflag_x,iflag_y,iflag_z
    
    !open file
    open(ifn,file=trim(file_name), status='old')
    
    if(trim(format)=='cu') then
      !check grid information
      inum(:)=lg_ie(:)-lg_is(:)+1
      read(ifn,*); read (ifn,*); read (ifn,*); !skip
      allocate(rtmp1d(4))
      read (ifn,*) rtmp1d; inum_check(1)=int(rtmp1d(1)+1d-3);
      read (ifn,*) rtmp1d; inum_check(2)=int(rtmp1d(1)+1d-3);
      read (ifn,*) rtmp1d; inum_check(3)=int(rtmp1d(1)+1d-3);
      deallocate(rtmp1d)
      do ii=1,3
        if(inum(ii)/=inum_check(ii)) then
          if(comm_is_root(nproc_id_global)) write(*,*) "al_em, dl_em or num_rgrid_em do not mutch shape file."
          stop
        end if
      end do
      read (ifn,*); !skip
      
      !input data(general case)
      allocate(rtmp1d(6))
      ix=lg_is(1); iy=lg_is(2); iz=lg_is(3);
      do ii=1,int(inum(1)*inum(2)*inum(3)/6)
        read (ifn,*) rtmp1d
        do ij=1,6
          !check flag and write imat
          iflag_x=0; iflag_y=0; iflag_z=0;
          if(ix>=ng_is(1) .and. ix<=ng_ie(1)) iflag_x=1
          if(iy>=ng_is(2) .and. iy<=ng_ie(2)) iflag_y=1
          if(iz>=ng_is(3) .and. iz<=ng_ie(3)) iflag_z=1
          if(iflag_x==1 .and. iflag_y==1 .and. iflag_z==1) then
            imat(ix,iy,iz)=int(rtmp1d(ij)+1d-3)
          end if        
          
          !update iz, iy, ix 
          iz=iz+1                                          !iz
          if(iz>lg_ie(3))                      iz=lg_is(3) !iz
          if(iz==lg_is(3))                     iy=iy+1     !iy
          if(iy>lg_ie(2))                      iy=lg_is(2) !iy
          if(iz==lg_is(3) .and. iy==lg_is(2))  ix=ix+1     !ix
        end do
      end do
      deallocate(rtmp1d)
      
      !input data(special case)
      if(mod(inum(1)*inum(2)*inum(3),6)>0) then
        allocate(rtmp1d(mod(inum(1)*inum(2)*inum(3),6)))
        read (ifn,*) rtmp1d
        do ij=1,mod(inum(1)*inum(2)*inum(3),6)
          !check flag and write imat
          iflag_x=0; iflag_y=0; iflag_z=0;
          if(ix>=ng_is(1) .and. ix<=ng_ie(1)) iflag_x=1
          if(iy>=ng_is(2) .and. iy<=ng_ie(2)) iflag_y=1
          if(iz>=ng_is(3) .and. iz<=ng_ie(3)) iflag_z=1
          if(iflag_x==1 .and. iflag_y==1 .and. iflag_z==1) then
            imat(ix,iy,iz)=int(rtmp1d(ij)+1d-3)
          end if        
          
          !update iz, iy, ix 
          iz=iz+1                                         !iz
          if(iz> lg_ie(3))                    iz=lg_is(3) !iz
          if(iz==lg_is(3))                    iy=iy+1     !iy
          if(iy> lg_ie(2))                    iy=lg_is(2) !iy
          if(iz==lg_is(3) .and. iy==lg_is(2)) ix=ix+1     !ix
        end do      
        deallocate(rtmp1d)
      end if
    elseif(trim(format)=='mp') then
      !This format is still unimplemented
    end if
    
    !close file
    close(ifn)
    
    return
  end subroutine input_i3d_em
  
  !=========================================================================================
  != output i3d file =======================================================================
  subroutine output_i3d_em(mg,lg,ifn,imessage,imat,format,file_name,Nd,coo)
    use inputoutput,     only: ulength_from_au
    use parallelization, only: nproc_id_global,nproc_group_global
    use communication,   only: comm_is_root,comm_summation
    use structures,      only: s_rgrid
    implicit none
    type(s_rgrid), intent(in) :: lg,mg
    integer,       intent(in) :: ifn,imessage
    integer,       intent(in) :: imat(mg%is_array(1):mg%ie_array(1),&
                                      mg%is_array(2):mg%ie_array(2),&
                                      mg%is_array(3):mg%ie_array(3))
    character(2),  intent(in) :: format
    character(256),intent(in) :: file_name
    integer,       intent(in) :: Nd
    real(8),       intent(in) :: coo(minval(lg%is(:))-Nd:maxval(lg%ie(:))+Nd,3)
    integer(8) :: ii
    integer    :: ij,ik,ix,iy,iz
    integer    :: i1d_tmp1(6),i1d_tmp2(6)
    real(8)    :: dx,dy,dz,conv
    
    !output
    if(trim(format)=='cu') then
      !set spacing
      dx = coo(2,1) - coo(1,1)
      dy = coo(2,2) - coo(1,2)
      dz = coo(2,3) - coo(1,3)
      
      !open cube file and write basic information
      if(comm_is_root(nproc_id_global)) then
        open(ifn,file=trim(file_name)//".cube")
        if(imessage==1) then
          write(ifn,*) "An input shape described in a cube file."
          conv = ulength_from_au
        else
          write(ifn,*) "None"
          conv = 1.0d0
        end if
        write(ifn,'(1X,A)') "All values here are in a.u. and A hydrogen atom is used to set the origin of the model."
        write(ifn,'(i5,3f12.6)') 1,coo(lg%is(1),1)*conv,coo(lg%is(2),2)*conv,coo(lg%is(3),3)*conv
        write(ifn,'(i5,3f12.6)') lg%num(1),dx*conv,0.0d0,0.0d0
        write(ifn,'(i5,3f12.6)') lg%num(2),0.0d0,dy*conv,0.0d0
        write(ifn,'(i5,3f12.6)') lg%num(3),0.0d0,0.0d0,dz*conv
        write(ifn,'(i5,4f12.6)') 1,1.0d0,0.0d0,0.0d0,0.0d0
      end if
      
      !write 3d data
      ix=lg%is(1); iy=lg%is(2); iz=lg%is(3);
      ij=1; i1d_tmp1(:)=0; i1d_tmp2(:)=0;
      do ii=1,lg%num(1)*lg%num(2)*lg%num(3)
        !collect data
        if( mg%is(1)<=ix .and. ix<=mg%ie(1) .and. &
            mg%is(2)<=iy .and. iy<=mg%ie(2) .and. &
            mg%is(3)<=iz .and. iz<=mg%ie(3) ) then
          i1d_tmp1(ij)=imat(ix,iy,iz);
        end if
        
        !(write data and reset i1d & ij) or (update ij)
        if(mod(ij,6)==0) then
          call comm_summation(i1d_tmp1(:),i1d_tmp2(:),6,nproc_group_global)
          if(comm_is_root(nproc_id_global)) write(ifn,'(6(1X,I2))', advance="yes") i1d_tmp2(:)
          i1d_tmp1(:)=0; i1d_tmp2(:)=0;
          ij=1;
        else
          ij=ij+1;
        end if
        
        !update iz
        iz=iz+1;
        if(iz>lg%ie(3)) iz=lg%is(3);
        
        !update iy
        if(iz==lg%is(3)) iy=iy+1;
        if(iy>lg%ie(2)) iy=lg%is(2);
        
        !update ix
        if(iz==lg%is(3) .and. iy==lg%is(2)) ix=ix+1;
        
        !final output for special case
        if(ii==lg%num(1)*lg%num(2)*lg%num(3) .and. ij>1) then
          call comm_summation(i1d_tmp1(:),i1d_tmp2(:),6,nproc_group_global)
          do ik=1,(ij-1)
            if(comm_is_root(nproc_id_global)) write(ifn,'(1X,I2)', advance="no") i1d_tmp2(ik)
          end do
        end if
      end do
    elseif(trim(format)=='mp') then
      !This format is still unimplemented
    end if
    
    !close file
    if(comm_is_root(nproc_id_global)) close(ifn)
    
    return
  end subroutine output_i3d_em
  
end module common_maxwell
