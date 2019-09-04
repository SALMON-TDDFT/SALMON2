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
subroutine calcVpsl_periodic_FFTE(lg,ng,info_field,poisson)
  use structures,      only: s_rgrid,s_field_parallel,s_poisson
  use salmon_parallel, only: nproc_group_global, nproc_size_global, nproc_id_global
  use salmon_communication, only: comm_bcast, comm_summation, comm_is_root
  use scf_data
  use new_world_sub
  use allocate_psl_sub
  use allocate_mat_sub
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: ng
  type(s_field_parallel),intent(in) :: info_field
  type(s_poisson),intent(inout) :: poisson
  
  integer :: ii,ix,iy,iz,ak
  integer :: iix,iiy,iiz
  integer :: iatom
  
  integer :: n
  real(8) :: bLx,bLy,bLz
  real(8) :: aLxyz
  integer :: NG_s,NG_e
  integer :: NG_l_s_para,NG_l_e_para
  integer :: numtmp
  complex(8),parameter :: zI=(0.d0,1.d0)
  real(8) :: G2sq,G2
  real(8) :: Gd
  real(8) :: s
  real(8) :: r
  integer :: imax
  integer :: iy_sta,iy_end,iz_sta,iz_end
  integer :: i,iix2,iiy2,iiz2
  
  integer :: npuy,npuz

  npuy=info_field%isize_ffte(2)
  npuz=info_field%isize_ffte(3)
  
  if(.not.allocated(poisson%a_ffte))then
    allocate(poisson%a_ffte(lg%num(1),lg%num(2)/npuy,lg%num(3)/npuz))
    allocate(poisson%b_ffte(lg%num(1),lg%num(2)/npuy,lg%num(3)/npuz))
  end if

!calculate reciprocal lattice vector
  bLx=2.d0*Pi/(Hgs(1)*dble(lg%num(1)))
  bLy=2.d0*Pi/(Hgs(2)*dble(lg%num(2)))
  bLz=2.d0*Pi/(Hgs(3)*dble(lg%num(3)))

  iz_sta=1
  iz_end=lg%num(3)/npuz
  iy_sta=1
  iy_end=lg%num(2)/npuy
 
  NG_s=1
  NG_e=lg%num(1)*lg%num(2)*lg%num(3)
  
  numtmp=(NG_e-NG_s+1)/nproc_size_global
  
  NG_l_s_para = nproc_id_global*numtmp+1
  NG_l_e_para = (nproc_id_global+1)*numtmp
  if(nproc_id_global==nproc_size_global-1) NG_l_e_para=NG_e
  
  nGzero=-1
  
  do ak=1,MKI
    do ii=1,Mr(ak)
      vloctbl(ii,ak)=pp%vpp_f(ii,Lref(ak),ak)
    enddo
  end do

  do iz=1,lg%num(3)/npuz
  do iy=1,lg%num(2)/npuy
  do ix=1,lg%num(1)
    n=(iz-1)*lg%num(2)/npuy*lg%num(1)+(iy-1)*lg%num(1)+ix
    iix2=ix-1+lg%is(1)
    iiy2=iy-1+info_field%id_ffte(2)*lg%num(2)/npuy+lg%is(2)
    iiz2=iz-1+info_field%id_ffte(3)*lg%num(3)/npuz+lg%is(3)
    if(ix==1.and.iy==1.and.iz==1.and.info_field%id_ffte(3)==0.and.info_field%id_ffte(2)==0) nGzero=n
    iix=ix-1-lg%num(1)*(1+sign(1,(iix2-1-(lg%num(1)+1)/2)))/2
    iiy=iy-1+info_field%id_ffte(2)*lg%num(2)/npuy-lg%num(2)*(1+sign(1,(iiy2-1-(lg%num(2)+1)/2)))/2
    iiz=iz-1+info_field%id_ffte(3)*lg%num(3)/npuz-lg%num(3)*(1+sign(1,(iiz2-1-(lg%num(3)+1)/2)))/2
    Gx(n)=dble(iix)*bLx
    Gy(n)=dble(iiy)*bLy
    Gz(n)=dble(iiz)*bLz
  enddo
  enddo
  enddo

  dVloc_G(:,:)=0.d0
  do ak=1,MKI
    imax=min(Mr(ak),Nr-1)
    do iz=1,lg%num(3)/npuz
    do iy=1,lg%num(2)/npuy
    do ix=1,lg%num(1)
      n=(iz-1)*lg%num(2)/npuy*lg%num(1)+(iy-1)*lg%num(1)+ix
      G2sq=sqrt(Gx(n)**2+Gy(n)**2+Gz(n)**2)
      s=0.d0
      if (n == nGzero) then
        do i=2,imax
          r=pp%rad(i+1,ak) !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
          s=s+4*Pi*r**2*(vloctbl(i,ak)+Zps(ak)/r)*(pp%rad(i+2,ak)-pp%rad(i+1,ak))
        enddo
      else
        do i=2,imax
          r=pp%rad(i+1,ak) !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
          s=s+4*Pi*r**2*sin(G2sq*r)/(G2sq*r)*(vloctbl(i,ak)+Zps(ak)/r)*(pp%rad(i+2,ak)-pp%rad(i+1,ak))
        enddo
      endif
      dVloc_G(n,ak)=s
    enddo
    enddo
    enddo
  enddo
 
  aLxyz=Hvol*dble(lg%num(1)*lg%num(2)*lg%num(3))
  rhoion_G=0.d0
  do iatom=1,MI
    do iz=1,lg%num(3)/npuz
    do iy=1,lg%num(2)/npuy
    do ix=1,lg%num(1)
      n=(iz-1)*lg%num(2)/npuy*lg%num(1)+(iy-1)*lg%num(1)+ix
      rhoion_G(n)=rhoion_G(n)+Zps(Kion(iatom))/aLxyz*exp(-zI*(Gx(n)*Rion(1,iatom)+Gy(n)*Rion(2,iatom)+Gz(n)*Rion(3,iatom)))
    enddo
    enddo
    enddo
  enddo

  Vion_G=0.d0
  do iatom=1,MI
    ak=Kion(iatom)
    do iz=1,lg%num(3)/npuz
    do iy=1,lg%num(2)/npuy
    do ix=1,lg%num(1)
      n=(iz-1)*lg%num(2)/npuy*lg%num(1)+(iy-1)*lg%num(1)+ix
      G2=Gx(n)**2+Gy(n)**2+Gz(n)**2
      Gd=Gx(n)*Rion(1,iatom)+Gy(n)*Rion(2,iatom)+Gz(n)*Rion(3,iatom)
      Vion_G(n)=Vion_G(n)+dVloc_G(n,ak)*exp(-zI*Gd)/aLxyz
      if(n == nGzero) cycle
      Vion_G(n)=Vion_G(n)-4*Pi/G2*Zps(ak)*exp(-zI*Gd)/aLxyz
    enddo
    enddo
    enddo
  enddo

  CALL PZFFT3DV_MOD(poisson%a_ffte,poisson%b_ffte,lg%num(1),lg%num(2),lg%num(3),npuy,npuz,0 &
                   ,info_field%icomm_ffte(2),info_field%icomm_ffte(3))

  do iz=1,lg%num(3)/npuz
  do iy=1,lg%num(2)/npuy
  do ix=1,lg%num(1)
    n=(iz-1)*lg%num(2)/npuy*lg%num(1)+(iy-1)*lg%num(1)+ix
    poisson%b_ffte(ix,iy,iz)=Vion_G(n)*dble(lg%num(1)*lg%num(2)*lg%num(3))
  enddo
  enddo
  enddo

  CALL PZFFT3DV_MOD(poisson%b_ffte,poisson%a_ffte,lg%num(1),lg%num(2),lg%num(3),npuy,npuz,1 &
                   ,info_field%icomm_ffte(2),info_field%icomm_ffte(3))

!$OMP parallel do
  do iz = lg%is(3),lg%ie(3)
  do iy = lg%is(2),lg%ie(2)
  do ix = lg%is(1),lg%ie(1)
    matbox_l(ix,iy,iz)=0.d0
  end do
  end do
  end do
  if(info_field%isize_ffte(1)==1)then
!$OMP parallel do private(iiz,iiy)
    do iz=iz_sta,iz_end
      iiz=iz+info_field%id_ffte(3)*lg%num(3)/npuz
      do iy=iy_sta,iy_end
        iiy=iy+info_field%id_ffte(2)*lg%num(2)/npuy
        matbox_l(1:lg%ie(1),iiy,iiz)=poisson%a_ffte(1:lg%ie(1),iy,iz)
      end do
    end do
  else
!$OMP parallel do private(iiz,iiy,ix)
    do iz=iz_sta,iz_end
      iiz=iz+info_field%id_ffte(3)*lg%num(3)/npuz
      do iy=iy_sta,iy_end
        iiy=iy+info_field%id_ffte(2)*lg%num(2)/npuy
        do iix=ng%is(1),ng%ie(1)
          ix=iix-lg%is(1)+1
          matbox_l(iix,iiy,iiz)=poisson%a_ffte(ix,iy,iz)
        end do
      end do
    end do
  end if
  call comm_summation(matbox_l,matbox_l2,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)
!$OMP parallel do
  do iz = mg_sta(3),mg_end(3)
  do iy = mg_sta(2),mg_end(2)
  do ix = mg_sta(1),mg_end(1)
    Vpsl(ix,iy,iz)=matbox_l2(ix,iy,iz)
  end do
  end do
  end do

  return

end subroutine calcVpsl_periodic_FFTE
