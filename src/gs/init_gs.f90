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
  integer :: ik,io,is,iseed,a,ix,iy,iz
  real(8) :: xx,yy,zz,x1,y1,z1,rr,rnd,Xmax,Ymax,Zmax

  if(system%if_real_orbital) then

    call init_wf_rand

!  else
!
!    Xmax=0.d0 ; Ymax=0.d0 ; Zmax=0.d0
!    do a=1,natom
!      if ( abs(Rion(1,a)) > Xmax ) Xmax=abs(Rion(1,a))
!      if ( abs(Rion(2,a)) > Ymax ) Ymax=abs(Rion(2,a))
!      if ( abs(Rion(3,a)) > Zmax ) Zmax=abs(Rion(3,a))
!    end do
!
!    Xmax=Xmax+1.d0/au_length_aa ; Ymax=Ymax+1.d0/au_length_aa ; Zmax=Zmax+1.d0/au_length_aa
!
!    iseed=123
!    do is=1,system%nspin
!    do io=1,system%no
!      call quickrnd_ns ; x1=Xmax*(2.d0*rnd-1.d0)
!      call quickrnd_ns ; y1=Ymax*(2.d0*rnd-1.d0)
!      call quickrnd_ns ; z1=Zmax*(2.d0*rnd-1.d0)
!      if(info%io_s <= io .and. io <= info%io_e) then
!!$OMP parallel do collapse(2) private(iz,iy,ix,xx,yy,zz,rr)
!        do iz=mg%is(3),mg%ie(3)
!        do iy=mg%is(2),mg%ie(2)
!        do ix=mg%is(1),mg%ie(1)
!          xx=lg%coordinate(ix,1) ; yy=lg%coordinate(iy,2) ; zz=lg%coordinate(iz,3)
!          rr=sqrt((xx-x1)**2+(yy-y1)**2+(zz-z1)**2)
!          spsi%rwf(ix,iy,iz,is,io,1,1) = exp(-0.5d0*(rr*au_length_aa)**2)*(au_length_aa)**(3/2)
!        end do
!        end do
!        end do
!!$omp end parallel do
!      end if
!    end do
!    end do
!
!  end if

  else

    Xmax = sqrt(sum(system%primitive_a(1:3,1)**2))
    Ymax = sqrt(sum(system%primitive_a(1:3,2)**2))
    Zmax = sqrt(sum(system%primitive_a(1:3,3)**2))

    iseed=123
    do is=1,system%nspin
    do ik=1,system%nk
    do io=1,system%no
      call quickrnd_ns ; x1=Xmax*rnd
      call quickrnd_ns ; y1=Ymax*rnd
      call quickrnd_ns ; z1=Zmax*rnd
      if(info%ik_s <= ik .and. ik <= info%ik_e .and.   &
         info%io_s <= io .and. io <= info%io_e) then
!$OMP parallel do collapse(2) private(iz,iy,ix,xx,yy,zz,rr)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          xx=lg%coordinate(ix,1) ; yy=lg%coordinate(iy,2) ; zz=lg%coordinate(iz,3)
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

  subroutine quickrnd_ns
  implicit none
  integer,parameter :: im=6075,ia=106,ic=1283
  iseed=mod(iseed*ia+ic,im) ; rnd=real(iseed,8)/real(im,8)
  end subroutine quickrnd_ns

  subroutine init_wf_rand
    implicit none
    integer :: s,k,n,i
    integer,allocatable :: ir(:)
    real(8) :: u

    if (.not. allocated(spsi%rwf)) stop 'spsi%rwf not allocated.'

    call random_seed( size=n )
    allocate( ir(n) )
    ir(:) = info%io_s + mg%is(1)
    call random_seed( put=ir )
    deallocate( ir )

    do io=lbound(spsi%rwf,5),ubound(spsi%rwf,5)
    do is=lbound(spsi%rwf,4),ubound(spsi%rwf,4)
    do iz=lbound(spsi%rwf,3),ubound(spsi%rwf,3)
    do iy=lbound(spsi%rwf,2),ubound(spsi%rwf,2)
    do ix=lbound(spsi%rwf,1),ubound(spsi%rwf,1)
      call random_number(u)
      spsi%rwf(ix,iy,iz,is,io,1,1) = u
    end do
    end do
    end do
    end do
    end do
  end subroutine

END SUBROUTINE init_wf

end module init_gs
