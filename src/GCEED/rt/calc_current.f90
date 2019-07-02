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
subroutine calc_current(mg,tpsi)
use salmon_parallel, only: nproc_group_global, nproc_group_korbital
use salmon_communication, only: comm_summation
use timer
use calc_allob_sub
use scf_data
use sendrecv_groupob_tmp_sub
use allocate_psl_sub
use structures, only: s_rgrid
implicit none
type(s_rgrid),intent(in) :: mg
complex(8) :: tpsi(mg%is_overlap(1):mg%ie_overlap(1) &
&                 ,mg%is_overlap(2):mg%ie_overlap(2) &
&                 ,mg%is_overlap(3):mg%ie_overlap(3), 1:iobnum, k_sta:k_end)
integer :: ix,iy,iz,iob,iik
complex(8),parameter :: zi=(0.d0,1.d0)
complex(8) :: ekr(maxMps,MI,k_sta:k_end)
real(8) :: r(3)
integer :: iatom,ik,jj,lm
complex(8) :: uVpsi0,uVpsix,uVpsiy,uVpsiz
real(8) :: jxt,jyt,jzt
real(8) :: curr1(3),curr2(3)
integer :: p_allob

call timer_begin(LOG_CALC_CURRENT)

iwk_size=2
call make_iwksta_iwkend

curr1(1:3)=0.d0


call timer_begin(LOG_CUR_SENDRECV)
call sendrecv_groupob_tmp(mg,tpsi)
call timer_end(LOG_CUR_SENDRECV)


call timer_begin(LOG_CUR_LOCAL)
do iik=k_sta,k_end
do iob=1,iobnum
  call calc_allob(iob,p_allob,iparaway_ob,itotmst,mst,iobnum)
  jxt=0.d0
!$OMP parallel do private(iz,iy,ix) reduction(+ : jxt)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    jxt=jxt+real(tpsi(ix,iy,iz,iob,iik))**2+aimag(tpsi(ix,iy,iz,iob,iik))**2
  end do
  end do
  end do
  jyt=jxt ; jzt=jxt
  curr1(1)=curr1(1)+rocc(p_allob,iik)*wtk(iik)*k_rd(1,iik)*jxt
  curr1(2)=curr1(2)+rocc(p_allob,iik)*wtk(iik)*k_rd(2,iik)*jyt
  curr1(3)=curr1(3)+rocc(p_allob,iik)*wtk(iik)*k_rd(3,iik)*jzt

  jxt=0.d0; jyt=0.d0; jzt=0.d0
!$OMP parallel do private(iz,iy,ix) reduction(+ : jxt,jyt,jzt)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    jxt=jxt+aimag(conjg(tpsi(ix,iy,iz,iob,iik))*   &
                            ( bN1*( tpsi(ix+1,iy,iz,iob,iik) - tpsi(ix-1,iy,iz,iob,iik) )    &
                             +bN2*( tpsi(ix+2,iy,iz,iob,iik) - tpsi(ix-2,iy,iz,iob,iik) )    &
                             +bN3*( tpsi(ix+3,iy,iz,iob,iik) - tpsi(ix-3,iy,iz,iob,iik) )    &
                             +bN4*( tpsi(ix+4,iy,iz,iob,iik) - tpsi(ix-4,iy,iz,iob,iik) ))/Hgs(1))
    jyt=jyt+aimag(conjg(tpsi(ix,iy,iz,iob,iik))*   &
                            ( bN1*( tpsi(ix,iy+1,iz,iob,iik) - tpsi(ix,iy-1,iz,iob,iik) )    &
                             +bN2*( tpsi(ix,iy+2,iz,iob,iik) - tpsi(ix,iy-2,iz,iob,iik) )    &
                             +bN3*( tpsi(ix,iy+3,iz,iob,iik) - tpsi(ix,iy-3,iz,iob,iik) )    &
                             +bN4*( tpsi(ix,iy+4,iz,iob,iik) - tpsi(ix,iy-4,iz,iob,iik) ))/Hgs(2))
    jzt=jzt+aimag(conjg(tpsi(ix,iy,iz,iob,iik))*   &
                            ( bN1*( tpsi(ix,iy,iz+1,iob,iik) - tpsi(ix,iy,iz-1,iob,iik) )    &
                             +bN2*( tpsi(ix,iy,iz+2,iob,iik) - tpsi(ix,iy,iz-2,iob,iik) )    &
                             +bN3*( tpsi(ix,iy,iz+3,iob,iik) - tpsi(ix,iy,iz-3,iob,iik) )    &
                             +bN4*( tpsi(ix,iy,iz+4,iob,iik) - tpsi(ix,iy,iz-4,iob,iik) ))/Hgs(3))
  end do
  end do
  end do
  curr1(1)=curr1(1)+rocc(p_allob,iik)*wtk(iik)*jxt
  curr1(2)=curr1(2)+rocc(p_allob,iik)*wtk(iik)*jyt
  curr1(3)=curr1(3)+rocc(p_allob,iik)*wtk(iik)*jzt
end do
end do
curr1(1:3)=curr1(1:3)*Hvol/(dble(lg_num(1)*lg_num(2)*lg_num(3))*Hvol)
call timer_end(LOG_CUR_LOCAL)


call timer_begin(LOG_CUR_NONLOCAL1)
jxt=0.d0;jyt=0.d0;jzt=0.d0
do iik=k_sta,k_end
!$OMP parallel do private(iatom,jj,ik,r)
  do iatom=1,MI
    ik=Kion(iatom)
    do jj=1,Mps(iatom)
      r(1)=(dble(Jxyz(1,jj,iatom)-1)-dble(Jxxyyzz(1,jj,iatom)*lg_num(1)))*Hgs(1)
      r(2)=(dble(Jxyz(2,jj,iatom)-1)-dble(Jxxyyzz(2,jj,iatom)*lg_num(2)))*Hgs(2)
      r(3)=(dble(Jxyz(3,jj,iatom)-1)-dble(Jxxyyzz(3,jj,iatom)*lg_num(3)))*Hgs(3)
      ekr(jj,iatom,iik)=exp(zi*(k_rd(1,iik)*r(1)   &
                          +k_rd(2,iik)*r(2)   &
                          +k_rd(3,iik)*r(3)))
    end do
  end do
  do iob=1,iobnum
    call calc_allob(iob,p_allob,iparaway_ob,itotmst,mst,iobnum)
!$OMP parallel do private(iatom,lm,ik,uVpsi0,uVpsix,uVpsiy,uVpsiz,r)
    do iatom=1,MI
      ik=Kion(iatom)
      do lm=1,(Mlps(ik)+1)**2
        if ( abs(uVu(lm,iatom))<1.d-5 ) then
          uVpsibox1_j(1:4,lm,iatom,iob,iik)=0.d0
        else
          uVpsi0=0.d0; uVpsix=0.d0; uVpsiy=0.d0; uVpsiz=0.d0
          do jj=1,Mps(iatom)
            r(1)=(dble(Jxyz(1,jj,iatom)-1)-dble(Jxxyyzz(1,jj,iatom)*lg_num(1)))*Hgs(1)
            r(2)=(dble(Jxyz(2,jj,iatom)-1)-dble(Jxxyyzz(2,jj,iatom)*lg_num(2)))*Hgs(2)
            r(3)=(dble(Jxyz(3,jj,iatom)-1)-dble(Jxxyyzz(3,jj,iatom)*lg_num(3)))*Hgs(3)
            uVpsi0=uVpsi0+uV(jj,lm,iatom)*ekr(jj,iatom,iik)      &
                     *tpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)
            uVpsix=uVpsix+uV(jj,lm,iatom)*ekr(jj,iatom,iik)*r(1) &
                     *tpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)
            uVpsiy=uVpsiy+uV(jj,lm,iatom)*ekr(jj,iatom,iik)*r(2) &
                     *tpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)
            uVpsiz=uVpsiz+uV(jj,lm,iatom)*ekr(jj,iatom,iik)*r(3) &
                     *tpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)
          end do  
          uVpsibox1_j(1,lm,iatom,iob,iik)=uVpsi0*Hvol/uVu(lm,iatom)
          uVpsibox1_j(2,lm,iatom,iob,iik)=uVpsix*Hvol
          uVpsibox1_j(3,lm,iatom,iob,iik)=uVpsiy*Hvol
          uVpsibox1_j(4,lm,iatom,iob,iik)=uVpsiz*Hvol
        end if
      end do
    end do
  end do
end do
call timer_end(LOG_CUR_NONLOCAL1)


call timer_begin(LOG_ALLREDUCE_CURRENT)
call timer_begin(LOG_CUR_NONLOCAL1_ALLREDUCE)
call comm_summation(uVpsibox1_j,uVpsibox2_j,4*maxlm*MI*iobnum*k_num,nproc_group_korbital)
call timer_end(LOG_CUR_NONLOCAL1_ALLREDUCE)
call timer_end(LOG_ALLREDUCE_CURRENT)


call timer_begin(LOG_CUR_NONLOCAL2)
do iik=k_sta,k_end
  do iob=1,iobnum
    call calc_allob(iob,p_allob,iparaway_ob,itotmst,mst,iobnum)
    do iatom=1,MI
      ik=Kion(iatom)
      do lm=1,(Mlps(ik)+1)**2
        jxt=jxt+rocc(p_allob,iik)*wtk(iik)/(lg_num(1)*lg_num(2)*lg_num(3)*Hvol)   &
                          *2.d0*aimag(conjg(uVpsibox2_j(2,lm,iatom,iob,iik))  &
                          *uVpsibox2_j(1,lm,iatom,iob,iik))/dble(nproc_Mxin_mul) 
        jyt=jyt+rocc(p_allob,iik)*wtk(iik)/(lg_num(1)*lg_num(2)*lg_num(3)*Hvol)   &
                          *2.d0*aimag(conjg(uVpsibox2_j(3,lm,iatom,iob,iik))  &
                          *uVpsibox2_j(1,lm,iatom,iob,iik))/dble(nproc_Mxin_mul)
        jzt=jzt+rocc(p_allob,iik)*wtk(iik)/(lg_num(1)*lg_num(2)*lg_num(3)*Hvol)   &
                          *2.d0*aimag(conjg(uVpsibox2_j(4,lm,iatom,iob,iik))  &
                          *uVpsibox2_j(1,lm,iatom,iob,iik))/dble(nproc_Mxin_mul)
      end do
    end do
  end do
end do

curr1(1)=curr1(1)+jxt
curr1(2)=curr1(2)+jyt
curr1(3)=curr1(3)+jzt
call timer_end(LOG_CUR_NONLOCAL2)


call timer_begin(LOG_CUR_NONLOCAL2_ALLREDUCE)
call comm_summation(curr1,curr2,3,nproc_group_global)
call timer_end(LOG_CUR_NONLOCAL2_ALLREDUCE)


curr(1:3,itt)=curr2(1:3)

if(iflag_indA==1)then
  A_ind(:,itt+1)=2.d0*A_ind(:,itt)-A_ind(:,itt-1)-4.d0*Pi*curr(:,itt)*dt**2
else if(iflag_indA==0)then
  A_ind(:,itt+1)=0.d0
end if

A_tot(:,itt+1)=A_ext(:,itt+1)+A_ind(:,itt+1)

E_ext(:,itt)=-(A_ext(:,itt+1)-A_ext(:,itt-1))/(2.d0*dt)
E_ind(:,itt)=-(A_ind(:,itt+1)-A_ind(:,itt-1))/(2.d0*dt)
E_tot(:,itt)=-(A_tot(:,itt+1)-A_tot(:,itt-1))/(2.d0*dt)

call timer_end(LOG_CALC_CURRENT)

end subroutine calc_current

subroutine calc_current_ion(system,j_ion)
  use structures, only: s_system
  use salmon_global, only: MI,Kion
  use scf_data, only: pp,lg_num,Hvol
  implicit none
  type(s_system) :: system
  integer :: ia
  real(8) :: j_ion(3)

  !AY memo
  !current of ion: defined by positive charge-->minus sign
  !This is NOT matter current NOR electric current.... strange definition....
  !This is defined so as to the total electric current = -(curr + curr_ion)
  !Should change this ion current but if you change, 
  !please change all part in ARTED, multiscale ..... 
  j_ion(:)=0d0
  do ia=1,MI
     j_ion(:) = j_ion(:) - pp%Zps(Kion(ia)) * system%Velocity(:,ia)
  enddo
  j_ion(:) = j_ion(:)/(dble(lg_num(1)*lg_num(2)*lg_num(3))*Hvol)

end subroutine calc_current_ion
