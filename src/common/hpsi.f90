!
!  Copyright 2017 SALMON developers
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
  integer,private,parameter :: Nd = 4 !????????

contains

!===================================================================================================================================

SUBROUTINE hpsi(tpsi,htpsi,info,mg,V_local,Nspin,stencil,srg,ppg,ttpsi)
  use structures
  use update_overlap_sub
  use stencil_sub
  use sendrecv_grid, only: s_sendrecv_grid, update_overlap_real8, update_overlap_complex8
  implicit none
  integer        ,intent(in) :: Nspin
  type(s_wf_info),intent(in) :: info
  type(s_rgrid)  ,intent(in) :: mg
  type(s_scalar) ,intent(in) :: V_local(Nspin)
  type(s_stencil),intent(in) :: stencil
  type(s_sendrecv_grid),intent(inout) :: srg
  type(s_pp_grid),intent(in) :: ppg
  type(s_wavefunction)       :: tpsi,htpsi
  type(s_wavefunction),optional :: ttpsi
  !
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb,ix,iy,iz
  real(8) :: k_nabt(Nd,3),k_lap0,kAc(3) !?????
  logical :: if_kAc

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  norb = Nspin* info%numo * info%numk * info%numm
  
  if_kAc = allocated(stencil%kAc)

  if(allocated(tpsi%rwf)) then

  ! overlap region communication
    if(info%if_divide_rspace) then
      call update_overlap_real8(srg, mg, tpsi%rwf)
!      call update_overlap_R(tpsi%rwf,mg%is_array,mg%ie_array,norb,Nd & !?????????
!                           ,mg%is,mg%ie,info%irank_r,info%icomm_r)
    end if
  ! stencil
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      call stencil_R(tpsi%rwf(:,:,:,ispin,io,ik,im),htpsi%rwf(:,:,:,ispin,io,ik,im),mg%is_array,mg%ie_array &
                    ,V_local(ispin)%f,mg%is,mg%ie &
                    ,mg%idx,mg%idy,mg%idz,stencil%lap0,stencil%lapt)
    end do
    end do
    end do
    end do
  ! pseudopotential
    call pseudo_R(tpsi,htpsi,info,Nspin,ppg)

  else

  ! overlap region communication
    if(info%if_divide_rspace) then
      !call update_overlap_C(tpsi%zwf,mg%is_array,mg%ie_array,norb,Nd & !????????
      !                     ,mg%is,mg%ie,info%irank_r,info%icomm_r)
      call update_overlap_complex8(srg, mg, tpsi%zwf)
    end if
  ! stencil
    if(stencil%if_orthogonal) then
    ! orthogonal lattice
      do im=im_s,im_e
      do ik=ik_s,ik_e
        if(if_kAc) then
          kAc(1:3) = stencil%kAc(ik,1:3)
          k_lap0 = stencil%lap0 + 0.5d0* sum(kAc(1:3)**2)
          k_nabt(:,1) = kAc(1) * stencil%nabt(:,1)
          k_nabt(:,2) = kAc(2) * stencil%nabt(:,2)
          k_nabt(:,3) = kAc(3) * stencil%nabt(:,3)
        else
          k_lap0 = stencil%lap0
          k_nabt = 0d0
        end if
        do io=io_s,io_e
        do ispin=1,Nspin
          call stencil_C(tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im),mg%is_array,mg%ie_array &
                        ,V_local(ispin)%f,mg%is,mg%ie &
                        ,mg%idx,mg%idy,mg%idz,k_lap0,stencil%lapt,k_nabt)
        end do
        end do
      end do
      end do
    else
      if(mg%ndir==3) then
      ! non-orthogonal lattice (general)
        if(.not.allocated(htpsi%wrk)) allocate(htpsi%wrk(mg%is_array(1):mg%ie_array(1) &
                                                        ,mg%is_array(2):mg%ie_array(2) &
                                                        ,mg%is_array(3):mg%ie_array(3),2) )
        do im=im_s,im_e
        do ik=ik_s,ik_e
          kAc = 0d0
          k_lap0 = 0d0
          if(if_kAc) then
            kAc(1:3) = stencil%kAc(ik,1:3) ! Cartesian vector
            k_lap0 = stencil%lap0 + 0.5d0* sum(kAc(1:3)**2)
            kAc(1:3) = matmul(stencil%matrix_B,kAc) ! B* (k+A/c)
          end if
          do io=io_s,io_e
          do ispin=1,Nspin
            call stencil_nonorthogonal(tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                                          ,mg%is_array,mg%ie_array,V_local(ispin)%f,mg%is,mg%ie &
                                          ,mg%idx,mg%idy,mg%idz,k_lap0,stencil%lapt,stencil%nabt &
                                          ,kAc,stencil%coef_F,htpsi%wrk)
          end do
          end do
        end do
        end do
      else if(mg%ndir > 3) then
      ! non-orthogonal lattice (high symmetry)
        do im=im_s,im_e
        do ik=ik_s,ik_e
          if(if_kAc) then
            kAc(1:3) = stencil%kAc(ik,1:3) ! Cartesian vector
            k_lap0 = stencil%lap0 + 0.5d0* sum(kAc(1:3)**2)
            k_nabt(:,1) = kAc(1) * stencil%nabt(:,1)
            k_nabt(:,2) = kAc(2) * stencil%nabt(:,2)
            k_nabt(:,3) = kAc(3) * stencil%nabt(:,3)
          else
            k_lap0 = stencil%lap0
            k_nabt = 0d0
          end if
          do io=io_s,io_e
          do ispin=1,Nspin
            call stencil_nonorthogonal_highsymmetry(tpsi%zwf(:,:,:,ispin,io,ik,im),htpsi%zwf(:,:,:,ispin,io,ik,im) &
                                          ,mg%is_array,mg%ie_array,V_local(ispin)%f,mg%is,mg%ie &
                                          ,mg%idx,mg%idy,mg%idz,k_lap0,stencil%coef_lap,k_nabt,mg%ndir,stencil%sign)
          end do
          end do
        end do
        end do
      end if
    end if
  ! subtraction
    if(present(ttpsi)) then
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
    end if
  ! pseudopotential
    call pseudo_C(tpsi,htpsi,info,nspin,ppg)

  end if

  return
end subroutine hpsi

!===================================================================================================================================

subroutine pseudo_R(tpsi,htpsi,info,nspin,ppg)
  use structures
  use salmon_communication, only: comm_summation
  implicit none
  integer,intent(in) :: nspin
  type(s_wf_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_wavefunction),intent(in) :: tpsi
  type(s_wavefunction) :: htpsi
  !
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,iorb,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma
  real(8) :: uVpsi,wrk
  real(8),allocatable :: uVpsibox(:,:),uVpsibox2(:,:)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  norb = Nspin* info%numo * info%numk * info%numm

  Nlma = ppg%Nlma

  if(info%if_divide_rspace) then

    allocate(uVpsibox(Nlma,Norb),uVpsibox2(Nlma,Norb))

    iorb = 0
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      iorb = iorb + 1
      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%rwf(ix,iy,iz,ispin,io,ik,im)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
        uVpsibox(ilma,iorb) = uVpsi
      end do
    end do
    end do
    end do
    end do
    call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,info%icomm_r)
    iorb = 0
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      iorb = iorb + 1
      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsibox2(ilma,iorb) * ppg%uV(j,ilma)
          htpsi%rwf(ix,iy,iz,ispin,io,ik,im) = htpsi%rwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do
    end do
    end do
    end do
    end do

    deallocate(uVpsibox,uVpsibox2)

  else

    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%rwf(ix,iy,iz,ispin,io,ik,im)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%uV(j,ilma)
          htpsi%rwf(ix,iy,iz,ispin,io,ik,im) = htpsi%rwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do
    end do
    end do
    end do
    end do

  end if

  return
end subroutine pseudo_R

!-----------------------------------------------------------------------------------------------------------------------------------

subroutine pseudo_C(tpsi,htpsi,info,nspin,ppg)
  use structures
  use salmon_communication, only: comm_summation
  implicit none
  integer,intent(in) :: nspin
  type(s_wf_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_wavefunction),intent(in) :: tpsi
  type(s_wavefunction) :: htpsi
  !
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,iorb,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma
  logical :: if_zproj
  complex(8) :: uVpsi,wrk
  complex(8),allocatable :: uVpsibox(:,:),uVpsibox2(:,:)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  norb = Nspin* info%numo * info%numk * info%numm

  Nlma = ppg%Nlma
  if_zproj = allocated(ppg%zproj)

!      ! gather (load) pseudo potential point
!      do i=1,NPI
!        spseudo(i) = tpsi(idx_proj(i))
!      end do

  if(info%if_divide_rspace) then

    allocate(uVpsibox(Nlma,Norb),uVpsibox2(Nlma,Norb))
    if(if_zproj) then

      iorb = 0
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin
        iorb = iorb + 1
        do ilma=1,Nlma
          ia = ppg%ia_tbl(ilma)
          uVpsi=0.d0
          do j=1,ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            uVpsi = uVpsi + conjg(ppg%zproj(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
          end do
          uVpsi = uVpsi * ppg%rinv_uvu(ilma)
          uVpsibox(ilma,iorb) = uVpsi
        end do
      end do
      end do
      end do
      end do
      call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,info%icomm_r)
      iorb = 0
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin
        iorb = iorb + 1
        do ilma=1,Nlma
          ia = ppg%ia_tbl(ilma)
          do j=1,ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            wrk = uVpsibox2(ilma,iorb) * ppg%zproj(j,ilma,ik)
            htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
          end do
        end do
      end do
      end do
      end do
      end do

    else

      iorb = 0
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin
        iorb = iorb + 1
        do ilma=1,Nlma
          ia = ppg%ia_tbl(ilma)
          uVpsi=0.d0
          do j=1,ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
          end do
          uVpsi = uVpsi * ppg%rinv_uvu(ilma)
          uVpsibox(ilma,iorb) = uVpsi
        end do
      end do
      end do
      end do
      end do
      call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,info%icomm_r)
      iorb = 0
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin
        iorb = iorb + 1
        do ilma=1,Nlma
          ia = ppg%ia_tbl(ilma)
          do j=1,ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            wrk = uVpsibox2(ilma,iorb) * ppg%uV(j,ilma)
            htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
          end do
        end do
      end do
      end do
      end do
      end do

    end if
    deallocate(uVpsibox,uVpsibox2)

  else

    if(if_zproj) then

      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin
        do ilma=1,Nlma
          ia = ppg%ia_tbl(ilma)
          uVpsi=0.d0
          do j=1,ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            uVpsi = uVpsi + conjg(ppg%zproj(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
          end do
          uVpsi = uVpsi * ppg%rinv_uvu(ilma)
          do j=1,ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            wrk = uVpsi * ppg%zproj(j,ilma,ik)
            htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
          end do
        end do
      end do
      end do
      end do
      end do

    else

      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin
        do ilma=1,Nlma
          ia = ppg%ia_tbl(ilma)
          uVpsi=0.d0
          do j=1,ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
          end do
          uVpsi = uVpsi * ppg%rinv_uvu(ilma)
          do j=1,ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            wrk = uVpsi * ppg%uV(j,ilma)
            htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
          end do
        end do
      end do
      end do
      end do
      end do

    end if

  end if

  return
end subroutine pseudo_C

!===================================================================================================================================

subroutine update_kvector_nonlocalpt(ppg,kAc,ik_s,ik_e)
  use structures
  implicit none
  type(s_pp_grid)    :: ppg
  integer,intent(in) :: ik_s,ik_e
  real(8),intent(in) :: kAc(ik_s:ik_e,3)
  !
  integer :: ilma,iatom,j,ik
  real(8) :: x,y,z,k(3)
  complex(8) :: ekr
  complex(8),parameter :: zi = (0d0,1d0)
  if(.not.allocated(ppg%zproj)) allocate(ppg%zproj(ppg%nps,ppg%nlma,ik_s:ik_e))
  do ik=ik_s,ik_e
    k = kAc(ik,:)
    do ilma=1,ppg%nlma
      iatom = ppg%ia_tbl(ilma)
      do j=1,ppg%mps(iatom)
        x = ppg%rxyz(1,j,iatom)
        y = ppg%rxyz(2,j,iatom)
        z = ppg%rxyz(3,j,iatom)
        ekr = exp(zi*(k(1)*x+k(2)*y+k(3)*z))
        ppg%zproj(j,ilma,ik) = conjg(ekr) * ppg%uv(j,ilma)
      end do
    end do
  end do
  return
end subroutine update_kvector_nonlocalpt

end module hpsi_sub
