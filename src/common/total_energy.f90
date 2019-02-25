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
MODULE Total_Energy

CONTAINS

  SUBROUTINE Total_Energy_eigensum(system,info,mg,V_local,Nspin,stencil,ppg)
    use structures
    use salmon_parallel, only: nproc_group_global
    use salmon_communication, only: comm_summation
    implicit none
    integer        ,intent(in) :: Nspin
    type(s_wf_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_scalar) ,intent(in) :: V_local(Nspin)
    type(s_stencil),intent(in) :: stencil
    type(s_pp_grid),intent(in) :: ppg
    type(s_system)             :: system
    !
    type(s_wavefunction)       :: htpsi

    integer :: io,ia,ib,ik,im,ispin
    integer :: ix,iy,iz,is(3),ie(3)
    real(8) :: rab
    real(8) :: sum1,sum2
    complex(8) :: cbox
    integer :: iob_allob

!??????
!    Etot=0.d0
!    do iik=1,num_kpoints_rd
!      Etot = Etot + sum( rocc(:itotMST,iik)*esp(:itotMST,iik) )*wtk(iik)
!    end do
!
!    do ia=1,MI
!    do ib=1,ia-1
!      rab=sqrt((Rion(1,ia)-Rion(1,ib))**2      &
!               +(Rion(2,ia)-Rion(2,ib))**2      &
!               +(Rion(3,ia)-Rion(3,ib))**2)
!      Etot=Etot+Zps(Kion(ia))*Zps(Kion(ib))/rab
!    end do
!    end do
!
!    if(ilsda == 0)then
!      sum1=0.d0
!      do iz=ng_sta(3),ng_end(3)
!      do iy=ng_sta(2),ng_end(2)
!      do ix=ng_sta(1),ng_end(1)
!        sum1=sum1+rho(ix,iy,iz)*(-0.5d0*Vh(ix,iy,iz)-Vxc(ix,iy,iz))
!      end do
!      end do
!      end do
!      call comm_summation(sum1,sum2,nproc_group_h)
!      Etot=Etot+sum2*Hvol+Exc
!    else if(ilsda == 1)then
!      sum1=0.d0
!      do iz=ng_sta(3),ng_end(3)
!      do iy=ng_sta(2),ng_end(2)
!      do ix=ng_sta(1),ng_end(1)
!        sum1=sum1+rho(ix,iy,iz)*(-0.5d0*Vh(ix,iy,iz))   &
!                 -Vxc_s(ix,iy,iz,1)*rho_s(ix,iy,iz,1)   &
!                 -Vxc_s(ix,iy,iz,2)*rho_s(ix,iy,iz,2)
!      end do
!      end do
!      end do
!      call comm_summation(sum1,sum2,nproc_group_korbital)
!      Etot=Etot+sum2*Hvol+Exc
!    end if

  end SUBROUTINE Total_Energy_eigensum

  SUBROUTINE Total_Energy_periodic()
  END SUBROUTINE Total_Energy_periodic

  Subroutine calc_eigen_energy(system,psi,hpsi_tmp,info,mg,V_local,Nspin,stencil,ppg)
    use structures
    use salmon_communication, only: comm_summation
    use hpsi_sub
    implicit none
    integer        ,intent(in) :: Nspin
    type(s_wf_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_scalar) ,intent(in) :: V_local(Nspin)
    type(s_stencil),intent(in) :: stencil
    type(s_pp_grid),intent(in) :: ppg
    type(s_wavefunction)       :: psi,hpsi_tmp
    type(s_system)             :: system
    !
    integer :: ik,io,ispin,im,nk,no,is(3),ie(3)
    real(8),allocatable :: wrk1(:,:),wrk2(:,:)

    if(info%im_s/=1 .or. info%im_e/=1) then
      write(*,*) "error: calc_eigen_energy"
      stop
    end if
    im = 1
    is = mg%is
    ie = mg%ie
    no = system%no
    nk = system%nk
    if(.not. allocated(system%esp)) allocate(system%esp(no,nk,Nspin))
    allocate(wrk1(no,nk),wrk2(no,nk))
    wrk1 = 0d0

    call hpsi(psi,hpsi_tmp,info,mg,V_local,Nspin,stencil,ppg)

    if(allocated(psi%rwf)) then
      do ispin=1,Nspin
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          wrk1(io,ik) = sum( psi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
                      * hpsi_tmp%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) * system%Hvol
        end do
        end do
        call comm_summation(wrk1,wrk2,no*nk,info%icomm_rko)
        system%esp(:,:,ispin) = wrk2
      end do
    else
      do ispin=1,Nspin
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          wrk1(io,ik) = sum( conjg( psi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) &
                             * hpsi_tmp%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) * system%Hvol
        end do
        end do
        call comm_summation(wrk1,wrk2,no*nk,info%icomm_rko)
        system%esp(:,:,ispin) = wrk2
      end do
    end if

    deallocate(wrk1,wrk2)
    return
  End Subroutine calc_eigen_energy

  Subroutine calc_ion_energy
    implicit none

!    if(iperiodic==0) then

!      do ia=1,MI
!      do ib=1,ia-1
!        rab=sqrt((Rion(1,ia)-Rion(1,ib))**2      &
!                 +(Rion(2,ia)-Rion(2,ib))**2      &
!                 +(Rion(3,ia)-Rion(3,ib))**2)
!        Etot=Etot+Zps(Kion(ia))*Zps(Kion(ib))/rab
!      end do
!      end do

!    else if(iperiodic==3) then
    ! Ewald sum

! ARTED
!      thr_id=0
!      Eion_tmp1=0.d0
!      Eion_l=0.d0
!      do ia=1,NI
!      do ix=-NEwald,NEwald
!      do iy=-NEwald,NEwald
!      do iz=-NEwald,NEwald
!      do ib=1,NI
!        if (ix**2+iy**2+iz**2 == 0 .and. ia == ib) then
!          cycle
!        end if
!        rab(1)=Rion(1,ia)-ix*aLx-Rion(1,ib)
!        rab(2)=Rion(2,ia)-iy*aLy-Rion(2,ib)
!        rab(3)=Rion(3,ia)-iz*aLz-Rion(3,ib)
!        rab2=sum(rab(:)**2)
!        Eion_tmp1=Eion_tmp1 + 0.5d0*Zps(Kion(ia))*Zps(Kion(ib))*erfc_salmon(sqrt(aEwald*rab2))/sqrt(rab2)
!      end do
!      end do
!      end do
!      end do
!      end do
!      do n=NG_s,NG_e
!        if(n == nGzero) cycle
!        G2=Gx(n)**2+Gy(n)**2+Gz(n)**2
!        Eion_l=Eion_l+aLxyz*(4*Pi/G2)*(abs(rhoion_G(n))**2*exp(-G2/(4*aEwald))*0.5d0)
!      end do
!      Eion_tmp1=Eion_tmp1-Pi*sum(Zps(Kion(:)))**2/(2*aEwald*aLxyz)-sqrt(aEwald/Pi)*sum(Zps(Kion(:))**2)

!    end if
    return
  end Subroutine calc_ion_energy

END MODULE Total_Energy
