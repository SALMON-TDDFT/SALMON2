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
  
  !===========================================================================================
  != input shape data ========================================================================
  subroutine input_shape_em(shape_file,ifn,ng_is,ng_ie,lg_is,lg_ie,Nd,imat,format)
    use parallelization, only: nproc_id_global
    use communication,   only: comm_is_root
    implicit none
    character(256),intent(in)  :: shape_file
    integer,       intent(in)  :: ifn,Nd
    integer,       intent(in)  :: ng_is(3),ng_ie(3),lg_is(3),lg_ie(3)
    integer,       intent(out) :: imat(ng_is(1)-Nd:ng_ie(1)+Nd,&
                                       ng_is(2)-Nd:ng_ie(2)+Nd,&
                                       ng_is(3)-Nd:ng_ie(3)+Nd)
    character(2),  intent(in)  :: format
    real(8),allocatable        :: rtmp1d(:)
    integer                    :: inum(3),inum_check(3)
    integer                    :: ii,ij,ix,iy,iz,iflag_x,iflag_y,iflag_z
    
    !open file
    open(ifn,file=trim(shape_file), status='old')
    
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
          if(comm_is_root(nproc_id_global)) write(*,*) "al_em or dl_em does not mutch shape file."
          stop
        end if
      end do
      read (ifn,*); !skip
      
      !input shape(general case)
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
      
      !input shape(special case)
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
          iz=iz+1                                          !iz
          if(iz>lg_ie(3))                      iz=lg_is(3) !iz
          if(iz==lg_is(3))                     iy=iy+1     !iy
          if(iy>lg_ie(2))                      iy=lg_is(2) !iy
          if(iz==lg_is(3) .and. iy==lg_is(2)) ix=ix+1      !ix
        end do      
        deallocate(rtmp1d)
      end if
    elseif(trim(format)=='mp') then
    end if
    
    !close file
    close(ifn)
    
    return
  end subroutine input_shape_em
  
end module common_maxwell
