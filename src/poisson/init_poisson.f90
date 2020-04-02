!
!  Copyright 2019 SALMON developers
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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------
module init_poisson_sub
  implicit none

contains

!=====================================================================
subroutine make_corr_pole(lg,mg,poisson)
  use salmon_global, only: natom,Rion,layout_multipole,num_multipole_xyz,al
  use inputoutput, only: au_length_aa
  use structures, only: s_rgrid,s_poisson
  implicit none
  type(s_rgrid), intent(in) :: lg
  type(s_rgrid), intent(in) :: mg
  type(s_poisson),intent(inout) :: poisson
  integer :: a,i
  integer :: ix,iy,iz
  integer :: ibox
  integer :: j1,j2,j3
  integer,allocatable :: ista_Mxin_pole(:,:)
  integer,allocatable :: iend_Mxin_pole(:,:)
  integer,allocatable :: inum_Mxin_pole(:,:)
  integer,allocatable :: iflag_pole(:)
  integer :: amin,amax
  real(8) :: rmin,r
  real(8),allocatable :: Rion2(:,:)
  integer,allocatable :: nearatomnum(:,:,:)
  integer,allocatable :: inv_icorr_polenum(:)
  integer :: maxval_ig_num
  real(8) :: dip_spacing

  if(layout_multipole==2)then
  
    amax=natom
    allocate(Rion2(3,natom))
    Rion2(:,:)=Rion(:,:)
  
    allocate(poisson%ig_num(1:amax))
    allocate(nearatomnum(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    poisson%ig_num=0
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rmin=1.d6
      do a=1,amax
        r=sqrt( (lg%coordinate(ix,1)-Rion2(1,a))**2      &
              + (lg%coordinate(iy,2)-Rion2(2,a))**2      &
              + (lg%coordinate(iz,3)-Rion2(3,a))**2 )
        if ( r < rmin ) then
          rmin=r ; amin=a
        end if
      end do
      poisson%ig_num(amin)=poisson%ig_num(amin)+1
      nearatomnum(ix,iy,iz)=amin
    end do
    end do
    end do
  
    allocate(poisson%ipole_tbl(1:amax))
    allocate(inv_icorr_polenum(1:amax))
    poisson%ipole_tbl=0
    inv_icorr_polenum=0
    ibox=0
    do a=1,amax
      if(poisson%ig_num(a)>=1)then
        ibox=ibox+1
        poisson%ipole_tbl(ibox)=a
        inv_icorr_polenum(a)=ibox
      end if
    end do
    poisson%npole_partial=ibox
  
    maxval_ig_num=maxval(poisson%ig_num(:)) 
    allocate(poisson%ig(3,maxval(poisson%ig_num(:)),poisson%npole_partial))
  
    poisson%ig_num(:)=0
  
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      ibox=inv_icorr_polenum(nearatomnum(ix,iy,iz))
      poisson%ig_num(ibox)=poisson%ig_num(ibox)+1
      poisson%ig(1,poisson%ig_num(ibox),ibox)=ix
      poisson%ig(2,poisson%ig_num(ibox),ibox)=iy
      poisson%ig(3,poisson%ig_num(ibox),ibox)=iz
    end do
    end do
    end do
  
    deallocate(Rion2)
    deallocate(nearatomnum)
    deallocate(inv_icorr_polenum)
  
  else if(layout_multipole==3)then
  
    if(num_multipole_xyz(1)==0.and.num_multipole_xyz(2)==0.and.num_multipole_xyz(3)==0)then
      dip_spacing = 8.d0/au_length_aa  ! approximate spacing of multipoles 
      poisson%n_multipole_xyz(1:3)=int((al(1:3)+dip_spacing)/dip_spacing-1.d-8)
    else
      poisson%n_multipole_xyz(1:3)=num_multipole_xyz(1:3)
    end if
    poisson%npole_total=poisson%n_multipole_xyz(1)*poisson%n_multipole_xyz(2)*poisson%n_multipole_xyz(3)

    allocate(ista_Mxin_pole(3,0:poisson%npole_total-1))
    allocate(iend_Mxin_pole(3,0:poisson%npole_total-1))
    allocate(inum_Mxin_pole(3,0:poisson%npole_total-1))
    allocate(iflag_pole(1:poisson%npole_total))
  
    do j3=0,poisson%n_multipole_xyz(3)-1
    do j2=0,poisson%n_multipole_xyz(2)-1
    do j1=0,poisson%n_multipole_xyz(1)-1
      ibox = j1 + poisson%n_multipole_xyz(1)*j2 + poisson%n_multipole_xyz(1)*poisson%n_multipole_xyz(2)*j3 
      ista_Mxin_pole(1,ibox)=j1*lg%num(1)/poisson%n_multipole_xyz(1)+lg%is(1)
      iend_Mxin_pole(1,ibox)=(j1+1)*lg%num(1)/poisson%n_multipole_xyz(1)+lg%is(1)-1
      ista_Mxin_pole(2,ibox)=j2*lg%num(2)/poisson%n_multipole_xyz(2)+lg%is(2)
      iend_Mxin_pole(2,ibox)=(j2+1)*lg%num(2)/poisson%n_multipole_xyz(2)+lg%is(2)-1
      ista_Mxin_pole(3,ibox)=j3*lg%num(3)/poisson%n_multipole_xyz(3)+lg%is(3)
      iend_Mxin_pole(3,ibox)=(j3+1)*lg%num(3)/poisson%n_multipole_xyz(3)+lg%is(3)-1
    end do
    end do
    end do
  
    iflag_pole=0
  
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      do i=1,poisson%npole_total
        if(ista_Mxin_pole(3,i-1)<=iz.and.iend_Mxin_pole(3,i-1)>=iz.and.   &
           ista_Mxin_pole(2,i-1)<=iy.and.iend_Mxin_pole(2,i-1)>=iy.and.   &
           ista_Mxin_pole(1,i-1)<=ix.and.iend_Mxin_pole(1,i-1)>=ix)then
          iflag_pole(i)=1
        end if
      end do
    end do
    end do
    end do
  
    poisson%npole_partial=0
    do i=1,poisson%npole_total
      if(iflag_pole(i)==1)then
        poisson%npole_partial=poisson%npole_partial+1
      end if
    end do
  
    allocate(poisson%ipole_tbl(1:poisson%npole_partial))
    allocate(poisson%ig_num(1:poisson%npole_partial))
  
    ibox=1
    do i=1,poisson%npole_total
      if(iflag_pole(i)==1)then
        poisson%ipole_tbl(ibox)=i
        ibox=ibox+1
      end if
    end do
  
    poisson%ig_num=0
  
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      do i=1,poisson%npole_partial
        if(ista_Mxin_pole(3,poisson%ipole_tbl(i)-1)<=iz.and.iend_Mxin_pole(3,poisson%ipole_tbl(i)-1)>=iz.and.   &
           ista_Mxin_pole(2,poisson%ipole_tbl(i)-1)<=iy.and.iend_Mxin_pole(2,poisson%ipole_tbl(i)-1)>=iy.and.   &
           ista_Mxin_pole(1,poisson%ipole_tbl(i)-1)<=ix.and.iend_Mxin_pole(1,poisson%ipole_tbl(i)-1)>=ix)then
          poisson%ig_num(i)=poisson%ig_num(i)+1
        end if
      end do
    end do
    end do
    end do
  
    maxval_ig_num=maxval(poisson%ig_num(:)) 
    allocate(poisson%ig(3,maxval(poisson%ig_num(:)),poisson%npole_partial))
   
    poisson%ig_num=0
  
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      do i=1,poisson%npole_partial
        if(ista_Mxin_pole(3,poisson%ipole_tbl(i)-1)<=iz.and.iend_Mxin_pole(3,poisson%ipole_tbl(i)-1)>=iz.and.   &
           ista_Mxin_pole(2,poisson%ipole_tbl(i)-1)<=iy.and.iend_Mxin_pole(2,poisson%ipole_tbl(i)-1)>=iy.and.   &
           ista_Mxin_pole(1,poisson%ipole_tbl(i)-1)<=ix.and.iend_Mxin_pole(1,poisson%ipole_tbl(i)-1)>=ix)then
          poisson%ig_num(i)=poisson%ig_num(i)+1
          poisson%ig(1,poisson%ig_num(i),i)=ix
          poisson%ig(2,poisson%ig_num(i),i)=iy
          poisson%ig(3,poisson%ig_num(i),i)=iz
        end if
      end do
    end do
    end do
    end do
  
    deallocate(ista_Mxin_pole,iend_Mxin_pole,inum_Mxin_pole)
    deallocate(iflag_pole)
  
  end if

  return
 
end subroutine make_corr_pole

!=====================================================================

subroutine set_ig_bound(lg,mg,poisson)
  use structures, only: s_rgrid,s_poisson
  implicit none
  type(s_rgrid), intent(in)     :: lg, mg
  type(s_poisson),intent(inout) :: poisson
  integer :: ix,iy,iz
  integer :: ibox
  integer :: icount
  integer,parameter :: ndh=4
  
  ibox=mg%num(1)*mg%num(2)*mg%num(3)/minval(mg%num(1:3))*2*ndh
  allocate( poisson%ig_bound(3,ibox,3) )
  
  icount=0
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=lg%is(1)-ndh,lg%is(1)-1
    icount=icount+1
    poisson%ig_bound(1,icount,1)=ix
    poisson%ig_bound(2,icount,1)=iy
    poisson%ig_bound(3,icount,1)=iz
  end do
  end do
  end do
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=lg%ie(1)+1,lg%ie(1)+ndh
    icount=icount+1
    poisson%ig_bound(1,icount,1)=ix
    poisson%ig_bound(2,icount,1)=iy
    poisson%ig_bound(3,icount,1)=iz
  end do
  end do
  end do
  icount=0
  do iz=mg%is(3),mg%ie(3)
  do iy=lg%is(2)-ndh,lg%is(2)-1
  do ix=mg%is(1),mg%ie(1)
    icount=icount+1
    poisson%ig_bound(1,icount,2)=ix
    poisson%ig_bound(2,icount,2)=iy
    poisson%ig_bound(3,icount,2)=iz
  end do
  end do
  end do
  do iz=mg%is(3),mg%ie(3)
  do iy=lg%ie(2)+1,lg%ie(2)+ndh
  do ix=mg%is(1),mg%ie(1)
    icount=icount+1
    poisson%ig_bound(1,icount,2)=ix
    poisson%ig_bound(2,icount,2)=iy
    poisson%ig_bound(3,icount,2)=iz
  end do
  end do
  end do
  icount=0
  do iz=lg%is(3)-ndh,lg%is(3)-1
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    icount=icount+1
    poisson%ig_bound(1,icount,3)=ix
    poisson%ig_bound(2,icount,3)=iy
    poisson%ig_bound(3,icount,3)=iz
  end do
  end do
  end do
  do iz=lg%ie(3)+1,lg%ie(3)+ndh
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    icount=icount+1
    poisson%ig_bound(1,icount,3)=ix
    poisson%ig_bound(2,icount,3)=iy
    poisson%ig_bound(3,icount,3)=iz
  end do
  end do
  end do
  
  return

end subroutine set_ig_bound

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------

end module init_poisson_sub
