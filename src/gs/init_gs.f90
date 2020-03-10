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

module init_gs
  implicit none

contains

!===================================================================================================================================

SUBROUTINE init_wf(lg,mg,system,info,spsi,pinfo)
  use structures
  use inputoutput, only: au_length_aa
  use salmon_global, only: yn_periodic,natom,rion
  use gram_schmidt_orth
  implicit none

  type(s_rgrid)           ,intent(in) :: lg,mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital)                     :: spsi
  type(s_process_info)    ,intent(in) :: pinfo
  !
  integer :: ik,io,is,a,ix,iy,iz,ip
  real(8) :: xx,yy,zz,x1,y1,z1,rr,Xmax,Ymax,Zmax,q(3)
  real(8) :: Xzero,Yzero,Zzero

  call init_wf_rand

  ! get offset (0-th element)
  if (yn_periodic == 'y') then
    Xzero = lg%coordinate(lg%num(1)/2+mod(lg%num(1),2),1)
    Yzero = lg%coordinate(lg%num(2)/2+mod(lg%num(2),2),2)
    Zzero = lg%coordinate(lg%num(3)/2+mod(lg%num(3),2),3)
  else
    Xzero = 0d0
    Yzero = 0d0
    Zzero = 0d0
  end if

  if(system%if_real_orbital) then

      Xmax=0.d0 ; Ymax=0.d0 ; Zmax=0.d0
      do a=1,natom
        if ( abs(Rion(1,a)) > Xmax ) Xmax=abs(Rion(1,a))
        if ( abs(Rion(2,a)) > Ymax ) Ymax=abs(Rion(2,a))
        if ( abs(Rion(3,a)) > Zmax ) Zmax=abs(Rion(3,a))
      end do

      Xmax=Xmax-Xzero+1.d0/au_length_aa
      Ymax=Ymax-Yzero+1.d0/au_length_aa
      Zmax=Zmax-Zzero+1.d0/au_length_aa

      do is=1,system%nspin
      do io=1,system%no
        call random_number(q)
        x1=Xmax*(2.d0*q(1)-1.d0)
        y1=Ymax*(2.d0*q(2)-1.d0)
        z1=Zmax*(2.d0*q(3)-1.d0)
        if(info%io_s <= io .and. io <= info%io_e) then
!$OMP parallel do collapse(2) private(iz,iy,ix,xx,yy,zz,rr)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            xx=lg%coordinate(ix,1)-Xzero
            yy=lg%coordinate(iy,2)-Yzero
            zz=lg%coordinate(iz,3)-Zzero
            rr=sqrt((xx-x1)**2+(yy-y1)**2+(zz-z1)**2)
            spsi%rwf(ix,iy,iz,is,io,1,1) = exp(-0.5d0*(rr*au_length_aa)**2)*(au_length_aa)**(3/2)
          end do
          end do
          end do
!$omp end parallel do
        end if
      end do
      end do

  else

      Xmax = sqrt(sum(system%primitive_a(1:3,1)**2))-Xzero
      Ymax = sqrt(sum(system%primitive_a(1:3,2)**2))-Yzero
      Zmax = sqrt(sum(system%primitive_a(1:3,3)**2))-Zzero

      do is=1,system%nspin
      do ik=1,system%nk
      do io=1,system%no
        call random_number(q)
        x1=Xmax*q(1)
        y1=Ymax*q(2)
        z1=Zmax*q(3)
        if(info%ik_s <= ik .and. ik <= info%ik_e .and.   &
           info%io_s <= io .and. io <= info%io_e) then
!$OMP parallel do collapse(2) private(iz,iy,ix,xx,yy,zz,rr)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            xx=lg%coordinate(ix,1)-Xzero
            yy=lg%coordinate(iy,2)-Yzero
            zz=lg%coordinate(iz,3)-Zzero
            rr=sqrt((xx-x1)**2+(yy-y1)**2+(zz-z1)**2)
            spsi%zwf(ix,iy,iz,is,io,ik,1) = exp(-0.5d0*rr**2)
          end do
          end do
          end do
!$omp end parallel do
        end if
      end do
      end do
      end do

  end if

  call gram_schmidt(system, mg, info, spsi, pinfo)

  return

CONTAINS

  ! cf. RSDFT
  subroutine init_wf_rand
    implicit none
    integer :: s,k,n,i,llen
    integer,allocatable :: iseed(:)

    call random_seed(size = n)
    allocate(iseed(n))
    llen = product(lg%num)
    iseed(:) = (info%ik_s * system%no + info%io_s) * llen &
             + (mg%is(3) - lg%is(3)) * lg%num(2) * lg%num(1) &
             + (mg%is(2) - lg%is(2)) * lg%num(1) &
             + (mg%is(1) - lg%is(1))
    call random_seed(put = iseed)
    deallocate(iseed)
  end subroutine

  subroutine gen_rwf_periodic
    implicit none
    real(8) :: u

    do ip=lbound(spsi%rwf,7),ubound(spsi%rwf,7)
    do ik=lbound(spsi%rwf,6),ubound(spsi%rwf,6)
    do io=lbound(spsi%rwf,5),ubound(spsi%rwf,5)
    do is=lbound(spsi%rwf,4),ubound(spsi%rwf,4)
    do iz=lbound(spsi%rwf,3),ubound(spsi%rwf,3)
    do iy=lbound(spsi%rwf,2),ubound(spsi%rwf,2)
    do ix=lbound(spsi%rwf,1),ubound(spsi%rwf,1)
      call random_number(u)
      spsi%rwf(ix,iy,iz,is,io,ik,ip) = u
    end do
    end do
    end do
    end do
    end do
    end do
    end do
  end subroutine

  subroutine gen_zwf_periodic
    implicit none
    real(8) :: u(2)

    do ip=lbound(spsi%zwf,7),ubound(spsi%zwf,7)
    do ik=lbound(spsi%zwf,6),ubound(spsi%zwf,6)
    do io=lbound(spsi%zwf,5),ubound(spsi%zwf,5)
    do is=lbound(spsi%zwf,4),ubound(spsi%zwf,4)
    do iz=lbound(spsi%zwf,3),ubound(spsi%zwf,3)
    do iy=lbound(spsi%zwf,2),ubound(spsi%zwf,2)
    do ix=lbound(spsi%zwf,1),ubound(spsi%zwf,1)
      call random_number(u)
      spsi%zwf(ix,iy,iz,is,io,ik,ip) = cmplx(u(1),u(2))
    end do
    end do
    end do
    end do
    end do
    end do
    end do
  end subroutine

END SUBROUTINE init_wf

end module init_gs
