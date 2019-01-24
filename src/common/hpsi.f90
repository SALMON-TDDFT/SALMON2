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
                 ,system,stencil &
                 ,pp,ttpsi &
                 ,nproc_Mxin_mul,irank_overlap,icomm_overlap,icomm_pseudo)
  use structures
  use salmon_structure
  use update_overlap_sub, only: update_overlap_R
  use stencil_sub, only: stencil_R
  implicit none
  integer,intent(in)  :: Nspin &
                        ,nproc_Mxin_mul,irank_overlap(6),icomm_overlap,icomm_pseudo
  type(s_rgrid)  ,intent(in) :: rg_wf
  type(struct_field)  ,intent(in) :: V_local(Nspin)
  type(struct_stencil),intent(in) :: stencil
  type(struct_pseudopotential),intent(in)  :: ppg !???????? --> type pp_grid
  type(struct_wavefunction)   ,intent(in)  :: tpsi
  type(struct_wavefunction)   ,intent(out) :: htpsi
  type(struct_wavefunction),optional,intent(out) :: ttpsi
  !
  integer :: is,ik,iorb
  real(8) :: k_nabt(4,3),k2
  logical :: if_kAc

  i1_s = tpsi%i1_s
  i1_e = tpsi%i1_e
  n1 = tpsi%n1

  i2_s = tpsi%i2_s
  i2_e = tpsi%i2_e
  n2 = tpsi%n2

  i3_s = tpsi%i3_s
  i3_e = tpsi%i3_e
  n3 = tpsi%n3

  if(allocated(tpsi%rwf)) then ! real or complex

  ! overlap region communication
    if(nproc_Mxin_mul.ne.1) then
      call update_overlap_R(tpsi%rwf,rg_wf%is_array,rg_wf%ie_array,Nspin*n1*n2*n3,4 & !?????????
                           ,rg_wf%is,rg_wf%ie,irank_overlap,icomm_overlap)
    end if
  ! stencil
    do i1=i1_s,i1_e!?????
    do i2=i2_s,i2_e
    do i3=i3_s,i3_e
    do ispin=1,Nspin
      call stencil_R(tpsi%rwf(:,:,:,ispin,i3,i2,i1),htpsi%rwf(:,:,:,ispin,i3,i2,i1),rg_wf%is_array,rg_wf%ie_array &
                    ,V_local(ispin)%f(:,:,:),rg_wf%is,rg_wf%ie &
                    ,rg_wf%idx,rg_wf%idy,rg_wf%idz,stencil%lap0,stencil%lapt)
    end do
    end do
    end do
    end do
  ! pseudopotential
    call pseudo_R(tpsi,htpsi,rg_wf &
                 ,ppg &
                 ,nproc_Mxin_mul,icomm_pseudo)

  else

  if_kAc = allocated(stencil%kAc) !?????

! overlap region communication
    if(nproc_Mxin_mul.ne.1) then
      call update_overlap_C(tpsi%zwf,rg_wf%is_array,rg_wf%ie_array,Nspin*n1*n2*n3,4 & !????????
                           ,rg_wf%is,rg_wf%ie,irank_overlap,icomm_overlap)
    end if
  ! stencil
    do i1=i1_s,i1_e !??????
    do i2=i2_s,i2_e
    do i3=i3_s,i3_e
    do ispin=1,Nspin
      if(if_kAc) then
        ik = table&kpoint(i3,i2,i1) !?????
        k2 = 0.5d0* sum(stencil%kAc(ik,:)**2)
        k_nabt(:,1) = stencil%kAc(ik,1) * stencil%nabt(:,1)
        k_nabt(:,2) = stencil%kAc(ik,2) * stencil%nabt(:,2)
        k_nabt(:,3) = stencil%kAc(ik,3) * stencil%nabt(:,3)
      else
        k2 = 0d0
        k_nabt = 0d0
      end if
      ! spin collinear
      call stencil_C(tpsi%zwf(:,:,:,ispin,i3,i2,i1),htpsi%zwf(:,:,:,ispin,i3,i2,i1),rg_wf%is_array,rg_wf%ie_array &
                    ,V_local(ispin)%f(:,:,:),rg_wf%is,rg_wf%ie &
                    ,rg_wf%idx,rg_wf%idy,rg_wf%idz,stencil%lap0+k2,stencil%lapt,k_nabt)
    end do
    end do
    end do
    end do
  ! subtraction
    if(present(ttpsi)) then
      do iorb=1,Norb
        is = table%ispin
        do iz=rg%ib(3),rg%ie(3)
          do iy=rg%ib(2),rg%ie(2)
            do ix=rg%ib(1),rg%ie(1)
            ! spin collinear
              ttpsi%zwf(ix,iy,iz,1,iorb) = htpsi%zwf(ix,iy,iz,1,iorb) - V_local(ix,iy,iz,is) * tpsi%zwf(ix,iy,iz,1,iorb)
            end do
          end do
        end do
      end do
    end if
  ! pseudopotential
    call pseudo_C(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                 ,NI,Nps,Nlma,ia_table,Mps,Jxyz,uV,uVu &
                 ,Nk,nproc_Mxin_mul,icomm_pseudo &
                 ,ik_table,exp_ikr)

  end if

  return
end subroutine hpsi

!===================================================================================================================================

subroutine pseudo_R(tpsi,htpsi,rg,Norb &
                   ,ppg &
                   ,nproc_Mxin_mul,icomm_pseudo)
  use salmon_pp
  use salmon_communication, only: comm_summation
  implicit none
  integer,intent(in)  :: Norb,nproc_Mxin_mul,icomm_pseudo
  type(struct_rgrid),intent(in) :: rg
  type(pp_grid),intent(in) :: ppg !??????
  type(struct_wavefunction) ,intent(in) :: tpsi
  type(struct_wavefunction),intent(out) :: htpsi
  !
  integer :: ilma,ia,j,ix,iy,iz,iorb,Nlma
  real(8) :: uVpsi,wrk
  real(8),allocatable :: uVpsibox(:,:),uVpsibox2(:,:)

  integer,pointer :: ia_table(:),Mps(:),Jxyz(:,:,:)
  real(8),pointer :: tpsi_(:,:,:,:,:),htpsi_(:,:,:,:,:) !???????
  real(8),pointer :: uV(:,:),uVu(:)

  Nlma = ppg%Nlma
  ia_table => ppg%ia_tbl
  Mps => ppg%Mps
  Jxyz => ppg%Jxyz
  uV => ppg%uV
  uVu => ppg%rinv_uvu

  tpsi_ => tpsi%rwf !?????
  htpsi_ => htpsi%rwf

  if(nproc_Mxin_mul.ne.1) then

    allocate(uVpsibox(Nlma,Norb),uVpsibox2(Nlma,Norb))

    do iorb=1,Norb !??????
      do ilma=1,Nlma
        ia = ia_table(ilma)
        uVpsi=0.d0
        do j=1,Mps(ia)
          ix = Jxyz(1,j,ia)
          iy = Jxyz(2,j,ia)
          iz = Jxyz(3,j,ia)
          uVpsi = uVpsi + uV(j,ilma) * tpsi_(ix,iy,iz,iorb)
        end do
        uVpsi = uVpsi * uVu(ilma)
        uVpsibox(ilma,iorb) = uVpsi
      end do
    end do
    call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,icomm_pseudo)
    do iorb=1,Norb
      do ilma=1,Nlma
        ia = ia_table(ilma)
        do j=1,Mps(ia)
          ix = Jxyz(1,j,ia)
          iy = Jxyz(2,j,ia)
          iz = Jxyz(3,j,ia)
          wrk = uVpsibox2(ilma,iorb) * uV(j,ilma)
          htpsi_(ix,iy,iz,iorb) = htpsi_(ix,iy,iz,iorb) + wrk
        end do
      end do
    end do

    deallocate(uVpsibox,uVpsibox2)

  else

    do iorb=1,Norb
      do ilma=1,Nlma
        ia = ia_table(ilma)
        uVpsi=0.d0
        do j=1,Mps(ia)
          ix = Jxyz(1,j,ia)
          iy = Jxyz(2,j,ia)
          iz = Jxyz(3,j,ia)
          uVpsi = uVpsi + uV(j,ilma) * tpsi_(ix,iy,iz,iorb)
        end do
        uVpsi = uVpsi * uVu(ilma) 
        do j=1,Mps(ia)
          ix = Jxyz(1,j,ia)
          iy = Jxyz(2,j,ia)
          iz = Jxyz(3,j,ia)
          wrk = uVpsi * uV(j,ilma)
          htpsi_(ix,iy,iz,iorb) = htpsi_(ix,iy,iz,iorb) + wrk
        end do
      end do
    end do

  end if

  return
end subroutine pseudo_R

!-----------------------------------------------------------------------------------------------------------------------------------

subroutine pseudo_C(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                   ,NI,Nps,Nlma,ia_table,Mps,Jxyz,uV,uVu &
                   ,Nk,nproc_Mxin_mul,icomm_pseudo &
                   ,ik_table,exp_ikr)
  use salmon_communication, only: comm_summation
  implicit none
  integer   ,intent(in) :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb &
                          ,NI,Nps,Nlma,ia_table(Nlma),Mps(NI),Jxyz(3,Nps,NI) &
                          ,Nk,nproc_Mxin_mul,icomm_pseudo
  complex(8),intent(in) :: tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  real(8)   ,intent(in) :: uV(Nps,Nlma),uVu(Nlma)
  !
  integer   ,optional,intent(in) :: ik_table(Norb)
  complex(8),optional,intent(in) :: exp_ikr(Nps,NI,Nk)
  !
  complex(8),intent(out) :: htpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  !
  integer    :: ilma,ia,j,ix,iy,iz,ik,iorb
  complex(8) :: uVpsi,wrk
  complex(8),allocatable :: uVpsibox(:,:),uVpsibox2(:,:)


  if(nproc_Mxin_mul.ne.1) then

    allocate(uVpsibox(Nlma,Norb),uVpsibox2(Nlma,Norb))
    if(present(exp_ikr)) then

      do iorb=1,Norb
        ik = ik_table(iorb)
        do ilma=1,Nlma
          ia = ia_table(ilma)
          uVpsi=0.d0
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            uVpsi = uVpsi + uV(j,ilma) * exp_ikr(j,ia,ik) * tpsi(ix,iy,iz,iorb)
          end do
          uVpsi = uVpsi * uVu(ilma)
          uVpsibox(ilma,iorb) = uVpsi
        end do
      end do
      call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,icomm_pseudo)
      do iorb=1,Norb
        ik = ik_table(iorb)
        do ilma=1,Nlma
          ia = ia_table(ilma)
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            wrk = conjg(exp_ikr(j,ia,ik)) * uVpsibox2(ilma,iorb) * uV(j,ilma)
            htpsi(ix,iy,iz,iorb) = htpsi(ix,iy,iz,iorb) + wrk
          end do
        end do
      end do

    else

      do iorb=1,Norb
        do ilma=1,Nlma
          ia = ia_table(ilma)
          uVpsi=0.d0
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            uVpsi = uVpsi + uV(j,ilma) * tpsi(ix,iy,iz,iorb)
          end do
          uVpsi = uVpsi * uVu(ilma)
          uVpsibox(ilma,iorb) = uVpsi
        end do
      end do
      call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,icomm_pseudo)
      do iorb=1,Norb
        do ilma=1,Nlma
          ia = ia_table(ilma)
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            wrk = uVpsibox2(ilma,iorb) * uV(j,ilma)
            htpsi(ix,iy,iz,iorb) = htpsi(ix,iy,iz,iorb) + wrk
          end do
        end do
      end do

    end if
    deallocate(uVpsibox,uVpsibox2)

  else

    if(present(exp_ikr)) then

      do iorb=1,Norb
        ik = ik_table(iorb)
        do ilma=1,Nlma
          ia = ia_table(ilma)
          uVpsi=0.d0
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            uVpsi = uVpsi + uV(j,ilma) * exp_ikr(j,ia,ik) * tpsi(ix,iy,iz,iorb)
          end do
          uVpsi = uVpsi * uVu(ilma)
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            wrk = conjg(exp_ikr(j,ia,ik)) * uVpsi * uV(j,ilma)
            htpsi(ix,iy,iz,iorb) = htpsi(ix,iy,iz,iorb) + wrk
          end do
        end do
      end do

    else

      do iorb=1,Norb
        do ilma=1,Nlma
          ia = ia_table(ilma)
          uVpsi=0.d0
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            uVpsi = uVpsi + uV(j,ilma) * tpsi(ix,iy,iz,iorb)
          end do
          uVpsi = uVpsi * uVu(ilma) 
          do j=1,Mps(ia)
            ix = Jxyz(1,j,ia)
            iy = Jxyz(2,j,ia)
            iz = Jxyz(3,j,ia)
            wrk = uVpsi * uV(j,ilma)
            htpsi(ix,iy,iz,iorb) = htpsi(ix,iy,iz,iorb) + wrk
          end do
        end do
      end do

    end if

  end if

  return
end subroutine pseudo_C

!===================================================================================================================================

end module hpsi_sub
