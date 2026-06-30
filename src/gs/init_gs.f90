!
!  Copyright 2019-2020 SALMON developers
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

SUBROUTINE init_wf(lg,mg,system,info,spsi)
  use structures
  use inputoutput, only: au_length_aa, method_init_wf
  use salmon_global, only: yn_periodic,natom,Rion,yn_jm
  use gram_schmidt_orth
  implicit none

  type(s_rgrid)           ,intent(in) :: lg,mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_parallel_info)   ,intent(in) :: info
  type(s_orbital)                     :: spsi
  !
  integer :: ik,io,is,a,ix,iy,iz,ip, ig,ngauss
  real(8) :: xx,yy,zz,x1,y1,z1,rr,Xmax,Ymax,Zmax,q(3)
  real(8) :: Xzero,Yzero,Zzero

  call init_wf_rand

  select case(method_init_wf)
  case ('gauss'  ) ; ngauss=1
  case ('gauss2' ) ; ngauss=2
  case ('gauss3' ) ; ngauss=3
  case ('gauss4' ) ; ngauss=4
  case ('gauss5' ) ; ngauss=5
  case ('gauss10') ; ngauss=10
  end select

  ! get offset (0-th element) : Xzero means center
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

    select case(method_init_wf)

    case ('random')
      call gen_random_rwf

    case ('gauss','gauss2','gauss3','gauss4','gauss5','gauss10')
      if(yn_jm=='n') then
         Xmax=0.d0 ; Ymax=0.d0 ; Zmax=0.d0
         do a=1,natom
           if ( abs(system%Rion(1,a)) > Xmax ) Xmax=abs(system%Rion(1,a))
           if ( abs(system%Rion(2,a)) > Ymax ) Ymax=abs(system%Rion(2,a))
           if ( abs(system%Rion(3,a)) > Zmax ) Zmax=abs(system%Rion(3,a))
         end do

         Xmax=Xmax-Xzero+1.d0/au_length_aa
         Ymax=Ymax-Yzero+1.d0/au_length_aa
         Zmax=Zmax-Zzero+1.d0/au_length_aa
      else
         Xmax=( lg%coordinate(lg%ie(1),1)-Xzero )*0.8d0
         Ymax=( lg%coordinate(lg%ie(2),2)-Yzero )*0.8d0
         Zmax=( lg%coordinate(lg%ie(3),3)-Zzero )*0.8d0
      end if

      do is=1,system%nspin
      do io=1,system%no
      do ig=1,ngauss
        call random_number(q)
        x1=Xmax*(2.d0*q(1)-1.d0)
        y1=Ymax*(2.d0*q(2)-1.d0)
        z1=Zmax*(2.d0*q(3)-1.d0)
        if(info%io_s <= io .and. io <= info%io_e) then
#ifdef USE_OPENACC
!$acc parallel loop collapse(2) private(iz,iy,ix,xx,yy,zz,rr)
#else
!$OMP parallel do collapse(2) private(iz,iy,ix,xx,yy,zz,rr)
#endif
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            xx=lg%coordinate(ix,1)-Xzero
            yy=lg%coordinate(iy,2)-Yzero
            zz=lg%coordinate(iz,3)-Zzero
            rr=sqrt((xx-x1)**2+(yy-y1)**2+(zz-z1)**2)
            if(ig==1) then
              spsi%rwf(ix,iy,iz,is,io,1,1) = exp(-0.5d0*(rr*au_length_aa)**2)*(au_length_aa)**(3/2)
            else
              spsi%rwf(ix,iy,iz,is,io,1,1) = spsi%rwf(ix,iy,iz,is,io,1,1) &
                   + exp(-0.5d0*(rr*au_length_aa)**2)*(au_length_aa)**(3/2)
            endif
          end do
          end do
          end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$omp end parallel do
#endif
        end if
      end do !ig
      end do
      end do

    end select

  else

    select case(method_init_wf)

    case ('random')
      call gen_random_zwf

    case ('gauss','gauss2','gauss3','gauss4','gauss5','gauss10')
      if(yn_jm=='n') then
        Xmax = sqrt(sum(system%primitive_a(1:3,1)**2))-Xzero
        Ymax = sqrt(sum(system%primitive_a(1:3,2)**2))-Yzero
        Zmax = sqrt(sum(system%primitive_a(1:3,3)**2))-Zzero
      else
         Xmax=( lg%coordinate(lg%ie(1),1)-Xzero )*0.8d0
         Ymax=( lg%coordinate(lg%ie(2),2)-Yzero )*0.8d0
         Zmax=( lg%coordinate(lg%ie(3),3)-Zzero )*0.8d0
      endif

      do is=1,system%nspin
      do ik=1,system%nk
      do io=1,system%no
      do ig=1,ngauss
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
            if(ig==1) then
               spsi%zwf(ix,iy,iz,is,io,ik,1) = exp(-0.5d0*rr**2)
            else
               spsi%zwf(ix,iy,iz,is,io,ik,1) = spsi%zwf(ix,iy,iz,is,io,ik,1) + exp(-0.5d0*rr**2)
            endif
          end do
          end do
          end do
!$omp end parallel do
        end if
      end do !ig
      end do
      end do
      end do

    end select

  end if

  call gram_schmidt(system, mg, info, spsi)

  return

CONTAINS

  ! cf. RSDFT
  subroutine init_wf_rand
    use salmon_global, only: iseed_number_change
    implicit none
    integer :: n, j
    integer(8) :: llen, seed64, seed64_mixed, seed_max, seed_value
    integer,allocatable :: iseed(:)
    logical :: duplicated

    call random_seed(size = n)
    allocate(iseed(n))
    llen = int(lg%num(1),8) * int(lg%num(2),8) * int(lg%num(3),8)
    seed_max = int(huge(iseed(1)),8)
    seed64 = (int(info%ik_s,8) * int(system%no,8) + int(info%io_s,8) - 1_8) * llen &
           + int(mg%is(3) - lg%is(3) + 1,8) * int(lg%num(2),8) * int(lg%num(1),8) &
           + int(mg%is(2) - lg%is(2) + 1,8) * int(lg%num(1),8) &
           + int(mg%is(1) - lg%is(1) + 1,8) + int(iseed_number_change,8)
    do j = 1, n
      seed64_mixed = splitmix64_next(seed64)
      seed_value = modulo(seed64_mixed, seed_max) + 1_8
      do
        duplicated = .false.
        if (j > 1) duplicated = any(iseed(1:j-1) == int(seed_value))
        if (.not. duplicated) exit
        seed_value = modulo(seed_value, seed_max) + 1_8
      end do
      iseed(j) = int(seed_value)
    end do
    call random_seed(put = iseed)
    deallocate(iseed)
  end subroutine

  integer(8) function splitmix64_next(state) result(mixed)
    implicit none
    integer(8), intent(inout) :: state
    integer(8), parameter :: splitmix_gamma(4) = [31765_8, 32586_8, 31161_8, 40503_8]
    integer(8), parameter :: splitmix_mul1(4) = [58809_8,  7396_8, 18285_8, 48984_8]
    integer(8), parameter :: splitmix_mul2(4) = [ 4587_8,  4913_8, 18875_8, 38096_8]

    ! SplitMix64 constants above are low-to-high 16-bit limbs of
    ! z'9E3779B97F4A7C15', z'BF58476D1CE4E5B9', z'94D049BB133111EB'.
    state = add_mod_2_64(state, splitmix_gamma)
    mixed = ieor(state, ishft(state, -30))
    mixed = mul_mod_2_64(mixed, splitmix_mul1)
    mixed = ieor(mixed, ishft(mixed, -27))
    mixed = mul_mod_2_64(mixed, splitmix_mul2)
    mixed = ieor(mixed, ishft(mixed, -31))
  end function

  integer(8) function add_mod_2_64(a, b_limb) result(res)
    implicit none
    integer(8), intent(in) :: a, b_limb(4)
    integer(8) :: carry, limb_sum, r_limb(4)
    integer :: k

    carry = 0_8
    do k = 1, 4
      limb_sum = limb16(a, k) + b_limb(k) + carry
      r_limb(k) = iand(limb_sum, 65535_8)
      carry = ishft(limb_sum, -16)
    end do
    res = pack_limbs(r_limb)
  end function

  integer(8) function mul_mod_2_64(a, b_limb) result(res)
    implicit none
    integer(8), intent(in) :: a, b_limb(4)
    integer(8) :: a_limb(4), r_limb(4), accum, carry
    integer :: i, j

    do i = 1, 4
      a_limb(i) = limb16(a, i)
    end do

    r_limb = 0_8
    do i = 1, 4
      carry = 0_8
      do j = 1, 5 - i
        accum = r_limb(i+j-1) + a_limb(i) * b_limb(j) + carry
        r_limb(i+j-1) = iand(accum, 65535_8)
        carry = ishft(accum, -16)
      end do
    end do
    res = pack_limbs(r_limb)
  end function

  integer(8) function limb16(value, index) result(limb)
    implicit none
    integer(8), intent(in) :: value
    integer, intent(in) :: index

    limb = iand(ishft(value, -16*(index - 1)), 65535_8)
  end function

  integer(8) function pack_limbs(limb) result(value)
    implicit none
    integer(8), intent(in) :: limb(4)

    value = ior(ior(limb(1), ishft(limb(2), 16)), &
              & ior(ishft(limb(3), 32), ishft(limb(4), 48)))
  end function

  subroutine gen_random_rwf
    implicit none
    real(8) :: u(3),v(3)

    do ip=lbound(spsi%rwf,7),ubound(spsi%rwf,7)
    do ik=lbound(spsi%rwf,6),ubound(spsi%rwf,6)
    do io=lbound(spsi%rwf,5),ubound(spsi%rwf,5)
    do is=lbound(spsi%rwf,4),ubound(spsi%rwf,4)
    do iz=lbound(spsi%rwf,3),ubound(spsi%rwf,3)
    do iy=lbound(spsi%rwf,2),ubound(spsi%rwf,2)
    do ix=lbound(spsi%rwf,1),ubound(spsi%rwf,1)
      v = dble([ix, iy, iz])
      call random_number(u)
      spsi%rwf(ix,iy,iz,is,io,ik,ip) = product(sign(u(1:3),v(1:3)))
    end do
    end do
    end do
    end do
    end do
    end do
    end do
  end subroutine

  subroutine gen_random_zwf
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
      spsi%zwf(ix,iy,iz,is,io,ik,ip) = dcmplx(u(1),u(2))
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
