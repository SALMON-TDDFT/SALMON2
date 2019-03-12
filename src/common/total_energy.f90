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
MODULE salmon_Total_Energy

CONTAINS

  SUBROUTINE calc_Total_Energy(Etot,Exc,system,info,ng,pp,rho,Vh,Vxc,fg)
    use structures
    use salmon_communication, only: comm_summation
    implicit none
    real(8)        ,intent(in) :: Exc
    type(s_system) ,intent(in) :: system
    type(s_wf_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: ng
    type(s_pp_info),intent(in) :: pp
    type(s_scalar) ,intent(in) :: rho(system%Nspin),Vh,Vxc(system%Nspin)
    type(s_fourier_grid),intent(in) :: fg
    real(8)                    :: Etot
    !
    integer :: io,ik,ispin,Nspin
    integer :: ix,iy,iz,is
    real(8) :: sum1,sum2,Eion
    Nspin = system%Nspin

!    if (Rion_update) then
    call calc_ion_energy(Eion,system,pp,fg) ! ion-ion energy
!    end if

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

  Subroutine calc_eigen_energy(system,psi,hpsi_tmp,info,mg,V_local,stencil,srg,ppg)
    use structures
    use salmon_communication, only: comm_summation
    use hpsi_sub
    use sendrecv_grid, only: s_sendrecv_grid
    implicit none
    type(s_system)             :: system
    type(s_wavefunction)       :: psi,hpsi_tmp
    type(s_wf_info),intent(in) :: info
    type(s_rgrid)  ,intent(in) :: mg
    type(s_scalar) ,intent(in) :: V_local(system%Nspin)
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
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

    call hpsi(psi,hpsi_tmp,info,mg,V_local,Nspin,stencil,srg,ppg)

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

  Subroutine calc_ion_energy(Eion,system,pp,fg)
    use structures
    use salmon_math
    use salmon_global, only: kion,NEwald,aEwald
    use salmon_communication, only: comm_summation
    implicit none
    type(s_system) ,intent(in) :: system
    type(s_pp_info),intent(in) :: pp
    type(s_fourier_grid),intent(in) :: fg
    real(8)                    :: Eion
    !
    integer :: ia,ib,ix,iy,iz,n
    real(8) :: r,rab(3),Eion_wrk,Eion_wrk2,G2
    real(8),parameter :: Pi=3.141592653589793d0 !??????????? salmon_math ? global parameter ?

    Eion = 0d0
    if(system%iperiodic==0) then
      do ia=1,system%nion
      do ib=1,ia-1
        r = sqrt((system%Rion(1,ia)-system%Rion(1,ib))**2      &
                +(system%Rion(2,ia)-system%Rion(2,ib))**2      &
                +(system%Rion(3,ia)-system%Rion(3,ib))**2)
        Eion = Eion + pp%Zps(Kion(ia)) * pp%Zps(Kion(ib)) /r
      end do
      end do
    else if(system%iperiodic==3) then
    ! Ewald sum
      do ia=1,system%nion
        do ix=-NEwald,NEwald
        do iy=-NEwald,NEwald
        do iz=-NEwald,NEwald
          do ib=1,system%nion
            if (ix**2+iy**2+iz**2 == 0 .and. ia == ib) then
              cycle
            end if
            rab(1) = system%Rion(1,ia)-ix*system%al(1,1) - system%Rion(1,ib)
            rab(2) = system%Rion(2,ia)-iy*system%al(2,2) - system%Rion(2,ib)
            rab(3) = system%Rion(3,ia)-iz*system%al(3,3) - system%Rion(3,ib)
            r=sum(rab(:)**2)
            Eion = Eion + 0.5d0*pp%Zps(Kion(ia))*pp%Zps(Kion(ib))*erfc_salmon(sqrt(aEwald*r))/sqrt(r)
          end do
        end do
        end do
        end do
      end do

      Eion_wrk = 0d0
      do n=fg%NG_s,fg%NG_e
        if(n == fg%nGzero) cycle
        G2 = fg%Gx(n)**2+fg%Gy(n)**2+fg%Gz(n)**2
        Eion_wrk = Eion_wrk + system%det_al*(4*Pi/G2)*(abs(fg%rhoion_G(n))**2*exp(-G2/(4*aEwald))*0.5d0)
      end do
      call comm_summation(Eion_wrk,Eion_wrk2,fg%icomm_fourier)

      Eion = Eion + Eion_wrk2 - Pi*sum(pp%Zps(Kion(:)))**2/(2*aEwald*system%det_al)-sqrt(aEwald/Pi)*sum(pp%Zps(Kion(:))**2)
    end if
    return
  end Subroutine calc_ion_energy

END MODULE salmon_Total_Energy
