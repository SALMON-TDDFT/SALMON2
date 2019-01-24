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

contains

!===================================================================================================================================

SUBROUTINE hpsi(tpsi,htpsi,rg_wf &
                 ,V_local,Nspin &
                 ,stencil &
                 ,ppg,ttpsi &
                 ,nproc_Mxin_mul,irank_overlap,icomm_overlap,icomm_pseudo)
  use structures
  use update_overlap_sub, only: update_overlap_R
  use stencil_sub, only: stencil_R
  implicit none
  integer,intent(in)  :: Nspin &
                        ,nproc_Mxin_mul,irank_overlap(6),icomm_overlap,icomm_pseudo
  type(s_rgrid)  ,intent(in) :: rg_wf
  type(s_scalar) ,intent(in) :: V_local(Nspin)
  type(s_stencil),intent(in) :: stencil
  type(s_pp_grid),intent(in) :: ppg
  type(s_wavefunction),intent(in)  :: tpsi
  type(s_wavefunction),intent(out) :: htpsi
  type(s_wavefunction),optional,intent(out) :: ttpsi
  !
  integer :: ispin,io,ik,i1,i1_s,i1_e,ik_s,ik_e,io_s,io_e,norb
  real(8) :: k_nabt(4,3),k_lap0
  logical :: if_kAc

  i1_s = tpsi%i1_s
  i1_e = tpsi%i1_e
  ik_s = tpsi%ik_s
  ik_e = tpsi%ik_e
  io_s = tpsi%io_s
  io_e = tpsi%io_e
  norb = Nspin* tpsi%numo * tpsi%numk * tpsi%num1

  if(allocated(tpsi%rwf)) then ! real or complex

  ! overlap region communication
    if(nproc_Mxin_mul.ne.1) then
      call update_overlap_R(tpsi%rwf,rg_wf%is_array,rg_wf%ie_array,norb,4 & !?????????
                           ,rg_wf%is,rg_wf%ie,irank_overlap,icomm_overlap)
    end if
  ! stencil
    do i1=i1_s,i1_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      call stencil_R(tpsi%rwf(:,:,:,ispin,io,ik,i1),htpsi%rwf(:,:,:,ispin,io,ik,i1),rg_wf%is_array,rg_wf%ie_array &
                    ,V_local(ispin)%f,rg_wf%is,rg_wf%ie &
                    ,rg_wf%idx,rg_wf%idy,rg_wf%idz,stencil%lap0,stencil%lapt)
    end do
    end do
    end do
    end do
  ! pseudopotential
    call pseudo_R(tpsi,htpsi,Nspin,ppg &
                 ,nproc_Mxin_mul,icomm_pseudo)

  else

  if_kAc = allocated(stencil%kAc) !?????

! overlap region communication
    if(nproc_Mxin_mul.ne.1) then
      call update_overlap_C(tpsi%zwf,rg_wf%is_array,rg_wf%ie_array,norb,4 & !????????
                           ,rg_wf%is,rg_wf%ie,irank_overlap,icomm_overlap)
    end if
  ! stencil
    do i1=i1_s,i1_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      if(if_kAc) then
        k_lap0 = stencil%lap0 + 0.5d0* sum(stencil%kAc(ik,:)**2)
        k_nabt(:,1) = stencil%kAc(ik,1) * stencil%nabt(:,1)
        k_nabt(:,2) = stencil%kAc(ik,2) * stencil%nabt(:,2)
        k_nabt(:,3) = stencil%kAc(ik,3) * stencil%nabt(:,3)
      else
        k_lap0 = stencil%lap0
        k_nabt = 0d0
      end if
      ! spin collinear
      call stencil_C(tpsi%zwf(:,:,:,ispin,io,ik,i1),htpsi%zwf(:,:,:,ispin,io,ik,i1),rg_wf%is_array,rg_wf%ie_array &
                    ,V_local(ispin)%f,rg_wf%is,rg_wf%ie &
                    ,rg_wf%idx,rg_wf%idy,rg_wf%idz,k_lap0,stencil%lapt,k_nabt)
    end do
    end do
    end do
    end do
  ! subtraction
    if(present(ttpsi)) then
      do i1=i1_s,i1_e
      do ik=ik_s,ik_e
      do io=io_s,io_e
      do ispin=1,Nspin
        ttpsi%zwf(:,:,:,ispin,io,ik,i1) = htpsi%zwf(:,:,:,ispin,io,ik,i1) &
          - V_local(ispin)%f(:,:,:) * tpsi%zwf(:,:,:,ispin,io,ik,i1)
      end do
      end do
      end do
      end do
    end if
  ! pseudopotential
    call pseudo_C(tpsi,htpsi,nspin,ppg &
                  ,nproc_Mxin_mul,icomm_pseudo)

  end if

  return
end subroutine hpsi

!===================================================================================================================================

subroutine pseudo_R(tpsi,htpsi,nspin,ppg &
                   ,nproc_Mxin_mul,icomm_pseudo)
  use structures
  use salmon_communication, only: comm_summation
  implicit none
  integer,intent(in) :: nspin,nproc_Mxin_mul,icomm_pseudo
  type(s_pp_grid),intent(in) :: ppg
  type(s_wavefunction)   ,intent(in) :: tpsi
  type(s_wavefunction),intent(inout) :: htpsi
  !
  integer :: ispin,io,ik,i1,i1_s,i1_e,ik_s,ik_e,io_s,io_e,iorb,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma
  real(8) :: uVpsi,wrk
  real(8),allocatable :: uVpsibox(:,:),uVpsibox2(:,:)

  i1_s = tpsi%i1_s
  i1_e = tpsi%i1_e
  ik_s = tpsi%ik_s
  ik_e = tpsi%ik_e
  io_s = tpsi%io_s
  io_e = tpsi%io_e
  norb = Nspin* tpsi%numo * tpsi%numk * tpsi%num1

  Nlma = ppg%Nlma

  if(nproc_Mxin_mul.ne.1) then

    allocate(uVpsibox(Nlma,Norb),uVpsibox2(Nlma,Norb))

    iorb = 0
    do i1=i1_s,i1_e
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
          uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%rwf(ix,iy,iz,ispin,io,ik,i1)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
        uVpsibox(ilma,iorb) = uVpsi
      end do
    end do
    end do
    end do
    end do
    call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,icomm_pseudo)
    iorb = 0
    do i1=i1_s,i1_e
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
          htpsi%rwf(ix,iy,iz,ispin,io,ik,i1) = htpsi%rwf(ix,iy,iz,ispin,io,ik,i1) + wrk
        end do
      end do
    end do
    end do
    end do
    end do

    deallocate(uVpsibox,uVpsibox2)

  else

    do i1=i1_s,i1_e
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
          uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%rwf(ix,iy,iz,ispin,io,ik,i1)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%uV(j,ilma)
          htpsi%rwf(ix,iy,iz,ispin,io,ik,i1) = htpsi%rwf(ix,iy,iz,ispin,io,ik,i1) + wrk
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

subroutine pseudo_C(tpsi,htpsi,nspin,ppg &
                   ,nproc_Mxin_mul,icomm_pseudo)
  use structures
  use salmon_communication, only: comm_summation
  implicit none
  integer,intent(in) :: nspin,nproc_Mxin_mul,icomm_pseudo
  type(s_pp_grid),intent(in) :: ppg
  type(s_wavefunction)   ,intent(in) :: tpsi
  type(s_wavefunction),intent(inout) :: htpsi
  !
  integer :: ispin,io,ik,i1,i1_s,i1_e,ik_s,ik_e,io_s,io_e,iorb,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma
  logical :: if_zproj
  complex(8) :: uVpsi,wrk
  complex(8),allocatable :: uVpsibox(:,:),uVpsibox2(:,:)

  i1_s = tpsi%i1_s
  i1_e = tpsi%i1_e
  ik_s = tpsi%ik_s
  ik_e = tpsi%ik_e
  io_s = tpsi%io_s
  io_e = tpsi%io_e
  norb = Nspin* tpsi%numo * tpsi%numk * tpsi%num1

  Nlma = ppg%Nlma
  if_zproj = allocated(ppg%zproj)

!      ! gather (load) pseudo potential point
!      do i=1,NPI
!        spseudo(i) = tpsi(idx_proj(i))
!      end do

  if(nproc_Mxin_mul.ne.1) then

    allocate(uVpsibox(Nlma,Norb),uVpsibox2(Nlma,Norb))
    if(if_zproj) then

      iorb = 0
      do i1=i1_s,i1_e
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
            uVpsi = uVpsi + conjg(ppg%zproj(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,i1)
          end do
          uVpsi = uVpsi * ppg%rinv_uvu(ilma)
          uVpsibox(ilma,iorb) = uVpsi
        end do
      end do
      end do
      end do
      end do
      call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,icomm_pseudo)
      iorb = 0
      do i1=i1_s,i1_e
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
            htpsi%zwf(ix,iy,iz,ispin,io,ik,i1) = htpsi%zwf(ix,iy,iz,ispin,io,ik,i1) + wrk
          end do
        end do
      end do
      end do
      end do
      end do

    else

      iorb = 0
      do i1=i1_s,i1_e
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
            uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%zwf(ix,iy,iz,ispin,io,ik,i1)
          end do
          uVpsi = uVpsi * ppg%rinv_uvu(ilma)
          uVpsibox(ilma,iorb) = uVpsi
        end do
      end do
      end do
      end do
      end do
      call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,icomm_pseudo)
      iorb = 0
      do i1=i1_s,i1_e
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
            htpsi%zwf(ix,iy,iz,ispin,io,ik,i1) = htpsi%zwf(ix,iy,iz,ispin,io,ik,i1) + wrk
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

      iorb = 0
      do i1=i1_s,i1_e
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
            uVpsi = uVpsi + conjg(ppg%zproj(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,i1)
          end do
          uVpsi = uVpsi * ppg%rinv_uvu(ilma)
          do j=1,ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            wrk = uVpsi * ppg%zproj(j,ilma,ik)
            htpsi%zwf(ix,iy,iz,ispin,io,ik,i1) = htpsi%zwf(ix,iy,iz,ispin,io,ik,i1) + wrk
          end do
        end do
      end do
      end do
      end do
      end do

    else

      iorb = 0
      do i1=i1_s,i1_e
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
            uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%zwf(ix,iy,iz,ispin,io,ik,i1)
          end do
          uVpsi = uVpsi * ppg%rinv_uvu(ilma)
          do j=1,ppg%mps(ia)
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            wrk = uVpsi * ppg%uV(j,ilma)
            htpsi%zwf(ix,iy,iz,ispin,io,ik,i1) = htpsi%zwf(ix,iy,iz,ispin,io,ik,i1) + wrk
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

end module hpsi_sub
