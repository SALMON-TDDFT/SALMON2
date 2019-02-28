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
MODULE salmon_Total_Energy

CONTAINS

  SUBROUTINE calc_Total_Energy(Etot,Exc,system,info,mg,ng,stencil,pp,ppg,rho,Vh,Vxc)
    use structures
    use salmon_communication, only: comm_summation
    implicit none
    real(8)        ,intent(in) :: Exc
    type(s_system) ,intent(in) :: system
    type(s_wf_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg,ng
    type(s_stencil),intent(in) :: stencil
    type(s_pp_info),intent(in) :: pp
    type(s_pp_grid),intent(in) :: ppg
    type(s_scalar) ,intent(in) :: rho(system%Nspin),Vh,Vxc(system%Nspin)
    real(8)                    :: Etot
    !
    integer :: io,ik,ispin,Nspin
    integer :: ix,iy,iz,is
    real(8) :: sum1,sum2,Eion
    Nspin = system%Nspin

    call calc_ion_energy(Eion,system,pp)

    Etot = 0d0
    do ispin=1,Nspin
    do ik=1,system%nk
    do io=1,system%no
      Etot = Etot + system%rocc(io,ik,ispin) * system%esp(io,ik,ispin)
    end do
    end do
    end do

    sum1 = 0d0
    do is=1,Nspin
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        sum1 = sum1 - 0.5d0* Vh%f(ix,iy,iz) * rho(is)%f(ix,iy,iz)    &
                    - ( Vxc(is)%f(ix,iy,iz) * rho(is)%f(ix,iy,iz) )
      end do
      end do
      end do
    end do
    
    call comm_summation(sum1,sum2,info%icomm_r)

    Etot = Etot + sum2*system%Hvol + Exc + Eion

    return
  end SUBROUTINE calc_Total_Energy

  Subroutine calc_eigen_energy(system,psi,hpsi_tmp,info,mg,V_local,stencil,ppg)
    use structures
    use salmon_communication, only: comm_summation
    use hpsi_sub
    implicit none
    type(s_system)             :: system
    type(s_wavefunction)       :: psi,hpsi_tmp
    type(s_wf_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_scalar) ,intent(in) :: V_local(system%Nspin)
    type(s_stencil),intent(in) :: stencil
    type(s_pp_grid),intent(in) :: ppg
    !
    integer :: ik,io,jo,ispin,im,nk,no,is(3),ie(3),Nspin
    real(8),allocatable :: wrk1(:,:),wrk2(:,:)

    Nspin = system%Nspin
    if(info%im_s/=1 .or. info%im_e/=1) then
      write(*,*) "error: calc_eigen_energy"
      stop
    end if
    im = 1
    is = mg%is
    ie = mg%ie
    no = system%no
    nk = system%nk
    allocate(wrk1(no,nk),wrk2(no,nk))
    wrk1 = 0d0

    call hpsi(psi,hpsi_tmp,info,mg,V_local,Nspin,stencil,ppg)

    if(allocated(psi%rwf)) then
      do ispin=1,Nspin
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          jo = info%io_tbl(io)
          wrk1(jo,ik) = sum( psi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) &
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
          jo = info%io_tbl(io)
          wrk1(jo,ik) = sum( conjg( psi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,ik,im) ) &
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

  Subroutine calc_ion_energy(Eion,system,pp)
    use structures
    use salmon_global, only: kion
    implicit none
    type(s_system) ,intent(in) :: system
    type(s_pp_info),intent(in) :: pp
    real(8)                    :: Eion
    !
    integer :: ia1,ia2
    real(8) :: r

!    if(iperiodic==0) then
    do ia1=1,system%nion
    do ia2=1,ia1-1
      r = sqrt((system%Rion(1,ia1)-system%Rion(1,ia2))**2      &
              +(system%Rion(2,ia1)-system%Rion(2,ia2))**2      &
              +(system%Rion(3,ia1)-system%Rion(3,ia2))**2)
      Eion = Eion + pp%Zps(Kion(ia1)) * pp%Zps(Kion(ia2)) /r !???????? zps
    end do
    end do
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

END MODULE salmon_Total_Energy
