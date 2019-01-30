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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
#define LOG_BEG(id) call timer_thread_begin(id)
#define LOG_END(id) call timer_thread_end(id)

#ifdef ARTED_USE_NVTX
#define NVTX_BEG(name,id)  call nvtxStartRange(name,id)
#define NVTX_END()         call nvtxEndRange()
#else
#define NVTX_BEG(name,id)
#define NVTX_END()
#endif

module hpsi
  use timer
  implicit none

contains
  subroutine hpsi_omp_KB_GS(ik,tpsi,ttpsi,htpsi)
    use Global_Variables, only: NL,NLz,NLy,NLx
    use opt_variables, only: zhtpsi,zttpsi,PNLx,PNLy,PNLz
    use omp_lib
    implicit none
    integer,intent(in)     :: ik
    complex(8),intent(in)  :: tpsi(NL)
    complex(8),intent(out) :: ttpsi(NL),htpsi(NL)
    integer :: tid

    LOG_BEG(LOG_HPSI)

    tid = omp_get_thread_num()
    call init(tpsi,zhtpsi(:,1,tid))
!    call hpsi_omp_KB_base(ik,zhtpsi(:,1,tid),zhtpsi(:,2,tid),zttpsi(:,tid))
    call hpsi_test3(ik,zhtpsi(:,1,tid),zhtpsi(:,2,tid),zttpsi(:,tid))
    call copyout(zhtpsi(:,2,tid),zttpsi(:,tid),htpsi,ttpsi)

   LOG_END(LOG_HPSI)

  contains
      subroutine init(zu,tpsi)
      implicit none
      complex(8),intent(in)  :: zu(0:NLz-1,0:NLy-1,0:NLx-1)
      complex(8),intent(out) :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
      integer :: ix,iy,iz

!dir$ vector aligned
      do ix=0,NLx-1
      do iy=0,NLy-1
      do iz=0,NLz-1
        tpsi(iz,iy,ix)=zu(iz,iy,ix)
      end do
      end do
      end do
    end subroutine

    subroutine copyout(zhtpsi,zttpsi,htpsi,ttpsi)
      implicit none
      complex(8), intent(in)  :: zhtpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
      complex(8), intent(in)  :: zttpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
      complex(8), intent(out) :: htpsi(0:NLz-1,0:NLy-1,0:NLx-1)
      complex(8), intent(out) :: ttpsi(0:NLz-1,0:NLy-1,0:NLx-1)
      integer :: ix,iy,iz

!dir$ vector aligned
      do ix=0,NLx-1
      do iy=0,NLy-1
      do iz=0,NLz-1
        htpsi(iz,iy,ix) = zhtpsi(iz,iy,ix)
        ttpsi(iz,iy,ix) = zttpsi(iz,iy,ix)
      end do
      end do
      end do
    end subroutine
  end subroutine

  subroutine hpsi_omp_KB_RT(ik,tpsi,htpsi)
    use opt_variables, only: PNLx,PNLy,PNLz
    implicit none
    integer,intent(in)     :: ik
    complex(8),intent(in)  ::  tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8),intent(out) :: htpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
!    call hpsi_omp_KB_base(ik,tpsi,htpsi)
    call hpsi_test3(ik,tpsi,htpsi)
  end subroutine

!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine hpsi_test3(ik,tpsi0,htpsi0,ttpsi0)
    use structures
    use hpsi_sub, only: hpsi
    use timer
    use Global_Variables, only: NLx,NLy,NLz,NL,kAc,lapx,lapy,lapz,nabx,naby,nabz,Vloc,ppg,zproj
    use opt_variables, only: PNLx,PNLy,PNLz
!    use Global_Variables, only: NLx,NLy,NLz,kAc,lapx,lapy,lapz,nabx,naby,nabz,Vloc,Mps,uV,iuV,Hxyz,ekr_omp,Nlma,a_tbl,Nps,NI,Jxyz
!    use opt_variables, only: lapt,PNLx,PNLy,PNLz,PNL
#ifdef ARTED_USE_NVTX
    use nvtx
#endif
    implicit none
    integer,intent(in)              :: ik
    complex(8),intent(in)           ::  tpsi0(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8),intent(out)          :: htpsi0(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8),intent(out),optional :: ttpsi0(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    !
    integer :: i,ix,iy,iz,irank_overlap(6),icomm_overlap,icomm_pseudo
    type(s_rgrid) :: rg
    type(s_stencil) :: stencil
    type(s_wavefunction) :: tpsi, htpsi, ttpsi
    type(s_scalar) :: V_local(1)

    stencil%lap0 = -(lapx(0)+lapy(0)+lapz(0))*0.5d0

    stencil%lapt(1:4,1) = lapz(1:4) ! x <--> z
    stencil%lapt(1:4,2) = lapy(1:4)
    stencil%lapt(1:4,3) = lapx(1:4) ! x <--> z

    stencil%nabt(1:4,1) = nabz(1:4) ! x <--> z
    stencil%nabt(1:4,2) = naby(1:4)
    stencil%nabt(1:4,3) = nabx(1:4) ! x <--> z

    allocate(stencil%kAc(1,3))
    stencil%kAc(1,1) = kAc(ik,3) ! x <--> z
    stencil%kAc(1,2) = kAc(ik,2)
    stencil%kAc(1,3) = kAc(ik,1) ! x <--> z

    rg%is = 0
    rg%ie(1) = NLz-1 ! x <--> z
    rg%ie(2) = NLy-1
    rg%ie(3) = NLx-1 ! x <--> z

    rg%is_overlap = rg%is - 4
    rg%ie_overlap = rg%ie + 4

    rg%is_array = 0
    rg%ie_array(1) = PNLz-1 ! x <--> z
    rg%ie_array(2) = PNLy-1
    rg%ie_array(3) = PNLx-1 ! x <--> z

    allocate(rg%idx(rg%is_overlap(1):rg%ie_overlap(1)) &
            ,rg%idy(rg%is_overlap(2):rg%ie_overlap(2)) &
            ,rg%idz(rg%is_overlap(3):rg%ie_overlap(3)))
    do i=rg%is_overlap(1),rg%ie_overlap(1)
      rg%idx(i) = mod(NLz+i,NLz) ! x <--> z
    end do
    do i=rg%is_overlap(2),rg%ie_overlap(2)
      rg%idy(i) = mod(NLy+i,NLy)
    end do
    do i=rg%is_overlap(3),rg%ie_overlap(3)
      rg%idz(i) = mod(NLx+i,NLx) ! x <--> z
    end do

    allocate(tpsi%zwf(rg%is_array(1):rg%ie_array(1),rg%is_array(2):rg%ie_array(2),rg%is_array(3):rg%ie_array(3),1,1,1,1) &
           ,htpsi%zwf(rg%is_array(1):rg%ie_array(1),rg%is_array(2):rg%ie_array(2),rg%is_array(3):rg%ie_array(3),1,1,1,1))
    if(present(ttpsi0)) then
      allocate(ttpsi%zwf(rg%is_array(1):rg%ie_array(1),rg%is_array(2):rg%ie_array(2),rg%is_array(3):rg%ie_array(3),1,1,1,1))
    end if
    tpsi%io_s = 1
    tpsi%io_e = 1
    tpsi%numo = 1
    tpsi%ik_s = 1
    tpsi%ik_e = 1
    tpsi%numk = 1
    tpsi%i1_s = 1
    tpsi%i1_e = 1
    tpsi%num1 = 1
    tpsi%zwf(:,:,:,1,1,1,1) = tpsi0

    allocate(V_local(1)%f(0:NLz-1,0:NLy-1,0:NLx-1))
    do i=1,NL
      iz = mod(i-1,NLz)
      iy = mod((i-1-iz)/NLz,NLy)
      ix = (i-1-iz-iy*NLy)/(NLy*NLz)
      V_local(1)%f(iz,iy,ix) = Vloc(i)
    end do

    if(.not.allocated(ppg%zproj)) allocate(ppg%zproj(ppg%nps,ppg%nlma,1))
    ppg%zproj(:,:,1) = zproj(:,:,ik)

    call hpsi(tpsi,htpsi,rg,V_local,1,stencil,ppg &
                 ,1,irank_overlap,icomm_overlap,icomm_pseudo,ttpsi)

    htpsi0 = htpsi%zwf(:,:,:,1,1,1,1)
    if(present(ttpsi0)) ttpsi0 = ttpsi%zwf(:,:,:,1,1,1,1)

    call deallocate_rgrid(rg)
    call deallocate_stencil(stencil)
    call deallocate_wavefunction(tpsi)
    call deallocate_wavefunction(htpsi)
    call deallocate_wavefunction(ttpsi)
    call deallocate_scalar(V_local(1))

    return
  end subroutine hpsi_test3
!-----------------------------------------------------------------------------------------------------------------------------------

end module
