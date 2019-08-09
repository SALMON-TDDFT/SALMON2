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
module hpsi_sub
  implicit none
  integer,private,parameter :: Nd = 4

contains

!===================================================================================================================================

SUBROUTINE hpsi(tpsi,htpsi,info,mg,V_local,system,stencil,srg,ppg,ttpsi)
  use structures
  use stencil_sub
  use pseudo_pt_sub
  use sendrecv_grid, only: s_sendrecv_grid, update_overlap_real8, update_overlap_complex8
  use timer
  implicit none
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_rgrid)  ,intent(in) :: mg
  type(s_scalar) ,intent(in) :: V_local(system%Nspin)
  type(s_stencil),intent(in) :: stencil
  type(s_sendrecv_grid),intent(inout) :: srg
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital)            :: tpsi,htpsi
  type(s_orbital),optional   :: ttpsi
  !
  integer :: nspin,ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb,ix,iy,iz
  real(8) :: k_nabt(Nd,3),k_lap0,kAc(3)
  logical :: if_kAc,if_singlescale

  call timer_begin(LOG_UHPSI_ALL)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  nspin = system%nspin
  norb = Nspin* info%numo * info%numk * info%numm
  
  if_kAc = allocated(stencil%vec_kAc)
  if_singlescale = allocated(system%Ac_micro%v)

  if(allocated(tpsi%rwf)) then

  ! overlap region communication
    call timer_begin(LOG_UHPSI_UPDATE_OVERLAP)
    if(info%if_divide_rspace) then
      call update_overlap_real8(srg, mg, tpsi%rwf)
    end if
    call timer_end(LOG_UHPSI_UPDATE_OVERLAP)

  ! stencil
    call timer_begin(LOG_UHPSI_STENCIL)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      call stencil_R(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz &
                    ,tpsi%rwf(:,:,:,ispin,io,ik,im),htpsi%rwf(:,:,:,ispin,io,ik,im) &
                    ,V_local(ispin)%f,stencil%coef_lap0,stencil%coef_lap)
    end do
    end do
    end do
    end do
    call timer_end(LOG_UHPSI_STENCIL)

  ! pseudopotential
    call pseudo_R(tpsi,htpsi,info,Nspin,ppg)

  else

  ! overlap region communication
    call timer_begin(LOG_UHPSI_UPDATE_OVERLAP)
    if(info%if_divide_rspace) then
      call update_overlap_complex8(srg, mg, tpsi%zwf)
    end if
    call timer_end(LOG_UHPSI_UPDATE_OVERLAP)

  ! stencil
    call timer_begin(LOG_UHPSI_STENCIL)
    if(stencil%if_orthogonal .and. .not.if_singlescale) then
    ! orthogonal lattice (general)
      do im=im_s,im_e
      do ik=ik_s,ik_e
        if(if_kAc) then
          kAc(1:3) = stencil%vec_kAc(1:3,ik)
          k_lap0 = stencil%coef_lap0 + 0.5d0* sum(kAc(1:3)**2)
          k_nabt(:,1) = kAc(1) * stencil%coef_nab(:,1)
          k_nabt(:,2) = kAc(2) * stencil%coef_nab(:,2)
          k_nabt(:,3) = kAc(3) * stencil%coef_nab(:,3)
        else
          k_lap0 = stencil%coef_lap0
          k_nabt = 0d0
        end if
        do io=io_s,io_e
        do ispin=1,Nspin
          call stencil_C(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz &
                        ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                        ,V_local(ispin)%f,k_lap0,stencil%coef_lap,k_nabt)
        end do
        end do
      end do
      end do
    else if(stencil%if_orthogonal .and. if_singlescale) then
    ! orthogonal lattice, sigle-scale Maxwell-TDDFT
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin
        call stencil_microAc(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz &
                      ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                      ,V_local(ispin)%f,system%Ac_micro%v,system%div_Ac%f,stencil%coef_lap0 &
                      ,stencil%coef_lap,stencil%coef_nab,system%vec_k(1:3,ik))
      end do
      end do
      end do
      end do
    else if(.not.stencil%if_orthogonal) then
    ! non-orthogonal lattice
      if(.not.allocated(htpsi%ztmp)) allocate(htpsi%ztmp(mg%is_array(1):mg%ie_array(1) &
                                                        ,mg%is_array(2):mg%ie_array(2) &
                                                        ,mg%is_array(3):mg%ie_array(3),2) )
      do im=im_s,im_e
      do ik=ik_s,ik_e
        kAc = 0d0
        k_lap0 = 0d0
        if(if_kAc) then
          kAc(1:3) = stencil%vec_kAc(1:3,ik) ! Cartesian vector k+A/c
          k_lap0 = stencil%coef_lap0 + 0.5d0* sum(kAc(1:3)**2)
          kAc(1:3) = matmul(stencil%rmatrix_B,kAc) ! B* (k+A/c)
        end if
        do io=io_s,io_e
        do ispin=1,Nspin
          call stencil_nonorthogonal(mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz,htpsi%ztmp &
                                    ,tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                                    ,V_local(ispin)%f,k_lap0,stencil%coef_lap,stencil%coef_nab,kAc,stencil%coef_F)
        end do
        end do
      end do
      end do
    end if
    call timer_end(LOG_UHPSI_STENCIL)

  ! subtraction
    call timer_begin(LOG_UHPSI_SUBTRACTION)
    if(present(ttpsi)) then
!$omp parallel do collapse(6) default(none) &
!$omp          private(im,ik,io,ispin,iz,iy,ix) &
!$omp          shared(im_s,im_e,ik_s,ik_e,io_s,io_e,nspin,mg) &
!$omp          shared(ttpsi,htpsi,V_local,tpsi)
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          ttpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) &
                                             - V_local(ispin)%f(ix,iy,iz) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
        end do
        end do
        end do
      end do
      end do
      end do
      end do
!$omp end parallel do
    end if
    call timer_end(LOG_UHPSI_SUBTRACTION)

  ! pseudopotential
    call pseudo_C(tpsi,htpsi,info,nspin,ppg)

  end if

  call timer_end(LOG_UHPSI_ALL)

  return
end subroutine hpsi

!===================================================================================================================================

subroutine update_kvector_nonlocalpt(ppg,kAc,ik_s,ik_e)
  use math_constants,only : zi
  use structures
  implicit none
  type(s_pp_grid)    :: ppg
  integer,intent(in) :: ik_s,ik_e
  real(8),intent(in) :: kAc(3,ik_s:ik_e)
  !
  integer :: ilma,iatom,j,ik
  real(8) :: x,y,z,k(3)
  complex(8) :: ekr
  if(.not.allocated(ppg%zekr_uV)) allocate(ppg%zekr_uV(ppg%nps,ppg%nlma,ik_s:ik_e))
  do ik=ik_s,ik_e
    k = kAc(:,ik)
    do ilma=1,ppg%nlma
      iatom = ppg%ia_tbl(ilma)
      do j=1,ppg%mps(iatom)
        x = ppg%rxyz(1,j,iatom)
        y = ppg%rxyz(2,j,iatom)
        z = ppg%rxyz(3,j,iatom)
        ekr = exp(zi*(k(1)*x+k(2)*y+k(3)*z))
        ppg%zekr_uV(j,ilma,ik) = conjg(ekr) * ppg%uv(j,ilma)
      end do
    end do
  end do
  return
end subroutine update_kvector_nonlocalpt

subroutine update_kvector_nonlocalpt_microAc(ik_s,ik_e,system,ppg)
  use math_constants,only : zi
  use structures
  implicit none
  integer           ,intent(in) :: ik_s,ik_e
  type(s_dft_system),intent(in) :: system
  type(s_pp_grid)               :: ppg
  !
  integer :: ilma,iatom,j,ik,ix,iy,iz,nj
  real(8) :: r(3),r_i(3),k(3),Ac(3),integral
  complex(8) :: ekr
  if(.not.allocated(ppg%zekr_uV)) allocate(ppg%zekr_uV(ppg%nps,ppg%nlma,ik_s:ik_e))
  do ik=ik_s,ik_e
    k = system%vec_k(:,ik)
    do ilma=1,ppg%nlma
      iatom = ppg%ia_tbl(ilma)
      nj = 0
      Ac = 0d0
      do j=1,ppg%mps(iatom)
        nj = nj + 1
        ix = ppg%jxyz(1,j,iatom)
        iy = ppg%jxyz(2,j,iatom)
        iz = ppg%jxyz(3,j,iatom)
        Ac = Ac + system%Ac_micro%v(:,ix,iy,iz)
      end do
      Ac = Ac/dble(nj) ! Ac averaged in the cutoff radius of the atom "iatom"
      do j=1,ppg%mps(iatom)
        r(1) = ppg%rxyz(1,j,iatom)
        r(2) = ppg%rxyz(2,j,iatom)
        r(3) = ppg%rxyz(3,j,iatom)
        r_i = system%rion(:,iatom)
        integral = Ac(1)*(r(1)-r_i(1)) + Ac(2)*(r(2)-r_i(2)) + Ac(3)*(r(3)-r_i(3))
        ekr = exp(zi*( k(1)*r(1)+k(2)*r(2)+k(3)*r(3) + integral ))
        ppg%zekr_uV(j,ilma,ik) = conjg(ekr) * ppg%uv(j,ilma)
      end do
    end do
  end do
  return
end subroutine update_kvector_nonlocalpt_microAc

end module hpsi_sub
