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

#include "config.h"

module opt_variables
  implicit none

  integer :: NUMBER_THREADS_POW2

  real(8) :: lapt(12)

  integer                :: PNLx,PNLy,PNLz,PNL
  complex(8),allocatable :: zhtpsi(:,:,:),zttpsi(:,:)
  complex(8),allocatable :: ghtpsi(:,:,:)
  complex(8),allocatable :: spseudo(:,:),dpseudo(:,:)  ! (NPI, # of threads)

  real(8),allocatable :: zrhotmp(:,:)

  integer,allocatable :: zJxyz(:,:),zKxyz(:,:)

  integer,allocatable :: modx(:),mody(:),modz(:)

  real(8),allocatable :: zcx(:,:),zcy(:,:),zcz(:,:)

  integer,allocatable    :: nprojector(:)       ! # of projector
  integer                :: NPI                 ! size of pseudo-vector (packed vector)
  integer,allocatable    :: idx_proj(:)         ! projector element index
  integer,allocatable    :: idx_lma(:)          ! start index of lma
  integer,allocatable    :: pseudo_start_idx(:) ! start index of pseudo-vector

#if defined(__KNC__) || defined(__AVX512F__)
# define MEM_ALIGNED 64
#else
# define MEM_ALIGNED 32
#endif

!dir$ attributes align:MEM_ALIGNED :: lapt
!dir$ attributes align:MEM_ALIGNED :: zhtpsi,zttpsi
!dir$ attributes align:MEM_ALIGNED :: zrhotmp
!dir$ attributes align:MEM_ALIGNED :: zJxyz,zKxyz
!dir$ attributes align:MEM_ALIGNED :: zcx,zcy,zcz
!dir$ attributes align:MEM_ALIGNED :: modx,mody,modz

contains
  subroutine opt_vars_initialize_p1
    use global_variables
    use misc_routines, only: ceiling_pow2
    implicit none

    select case(functional)
      case('TPSS','VS98')
        call err_finalize('functional: TPSS/VS98 versions not implemented.')
    end select

    NUMBER_THREADS_POW2 = ceiling_pow2(NUMBER_THREADS)
    allocate(zrhotmp(0:NL-1,0:NUMBER_THREADS_POW2-1))
    zrhotmp(:,:) = 0.0d0
  end subroutine

  subroutine opt_vars_initialize_p2
    use global_variables
    use code_optimization, only: switch_stencil_optimization, &
    &                            stencil_is_parallelized_by_omp, &
    &                            set_modulo_tables
    implicit none

    integer :: ix,iy,iz
    integer :: num(3)

    PNLx = NLx
#ifdef USE_OPT_ARRAY_PADDING
    PNLy = NLy + 1
#else
    PNLy = NLy
#endif
    PNLz = NLz
    PNL  = PNLx * PNLy * PNLz

    allocate(zhtpsi(0:PNL-1,4,0:NUMBER_THREADS-1))
    allocate(zttpsi(0:PNL-1,0:NUMBER_THREADS-1))

    allocate(zcx(NBoccmax,NK_s:NK_e))
    allocate(zcy(NBoccmax,NK_s:NK_e))
    allocate(zcz(NBoccmax,NK_s:NK_e))

    lapt( 1: 4)=lapx(1:4)
    lapt( 5: 8)=lapy(1:4)
    lapt( 9:12)=lapz(1:4)

    if(.not. allocated(zJxyz)) allocate(zJxyz(Nps,NI))  !AY see subroutine prep_ps_periodic
    !allocate(zJxyz(Nps,NI))

    zJxyz(1:Nps,1:NI) = Jxyz(1:Nps,1:NI) - 1

    allocate(modx(0:NLx*2+Nd-1))
    allocate(mody(0:NLy*2+Nd-1))
    allocate(modz(0:NLz*2+Nd-1))

    do ix=0,NLx*2+Nd-1
      modx(ix) = mod(ix,NLx)
    end do
    do iy=0,NLy*2+Nd-1
      mody(iy) = mod(iy,NLy)
    end do
    do iz=0,NLz*2+Nd-1
      modz(iz) = mod(iz,NLz)
    end do

#ifdef USE_OPT_ARRAY_PADDING
    call init_for_padding
    call init_projector(zKxyz)
#else
    call init_projector(zJxyz)
#endif

    num(1) = NLz
    num(2) = NLy
    num(3) = NLx
    call switch_stencil_optimization(num)
    ! In ARTED, we must disable openmp parallelization in the stencil.
    stencil_is_parallelized_by_omp = .false.
    call set_modulo_tables(num)
  end subroutine

  subroutine init_for_padding
    use global_variables
    implicit none
    integer :: a,ik,ix,iy,iz,jx,jy,jz,i,j,k
    real(8) :: x,y,z,r

    if(allocated(zKxyz)) deallocate(zKxyz)
    allocate(zKxyz(Nps,NI))

    do a=1,NI
      ik=Kion(a)
      j=0
      do ix=-2,2
      do iy=-2,2
      do iz=-2,2
        do jx=0,NLx-1
        do jy=0,NLy-1
        do jz=0,NLz-1
          i=jx* NLy* NLz + jy* NLz + jz + 1
          k=jx*PNLy*PNLz + jy*PNLz + jz
          x=Lx(i)*Hx-(Rion(1,a)+ix*aLx)
          y=Ly(i)*Hy-(Rion(2,a)+iy*aLy)
          z=Lz(i)*Hz-(Rion(3,a)+iz*aLz)
          r=sqrt(x*x+y*y+z*z)
          if (r<Rps(ik)) then
            j=j+1
            if(j<=Nps) then
              zKxyz(j,a)=k
            end if
          end if
        end do
        end do
        end do
      end do
      end do
      end do
    end do
  end subroutine

  function count_if_integer(vec, val) result(n)
    implicit none
    integer,intent(in) :: vec(:)
    integer,intent(in) :: val
    integer :: n, i
    n = 0
    do i=1,size(vec)
      if (vec(i) == val) then
        n = n + 1
      end if
    end do
  end function

  subroutine init_projector(zJxyz)
    use global_variables
    implicit none
    integer,intent(in) :: zJxyz(Nps,NI)
    integer :: i,ioffset

    NPI = sum(Mps)

    allocate(nprojector(NI),idx_lma(NI))

    do i=1,NI
      nprojector(i) = count_if_integer(a_tbl, i)
    end do

    idx_lma(1) = 0
    do i=2,NI
      idx_lma(i) = idx_lma(i-1) + nprojector(i-1)
    end do

    allocate(idx_proj(NPI))
    allocate(pseudo_start_idx(NI))

    ioffset = 0
    do i=1,NI
      pseudo_start_idx(i) = ioffset
      idx_proj(ioffset+1:ioffset+Mps(i)) = zJxyz(1:Mps(i),i)
      ioffset = ioffset + Mps(i)
    end do

    if(allocated(spseudo)) deallocate(spseudo,dpseudo)
    allocate(spseudo(NPI,0:NUMBER_THREADS-1))
    allocate(dpseudo(NPI,0:NUMBER_THREADS-1))
  end subroutine


  subroutine opt_vars_reinitialize
    implicit none

    call opt_vars_finalize
    call opt_vars_initialize_p1
    call opt_vars_initialize_p2
  end subroutine

  subroutine opt_vars_finalize
    implicit none

#define SAFE_DEALLOCATE(var) if(allocated(var)) deallocate(var)

    SAFE_DEALLOCATE(zhtpsi)
    SAFE_DEALLOCATE(zttpsi)

    SAFE_DEALLOCATE(zrhotmp)
    SAFE_DEALLOCATE(zJxyz)
    SAFE_DEALLOCATE(zKxyz)

    SAFE_DEALLOCATE(modx)
    SAFE_DEALLOCATE(mody)
    SAFE_DEALLOCATE(modz)

    SAFE_DEALLOCATE(zcx)
    SAFE_DEALLOCATE(zcy)
    SAFE_DEALLOCATE(zcz)

    SAFE_DEALLOCATE(nprojector)
    SAFE_DEALLOCATE(idx_proj)
    SAFE_DEALLOCATE(idx_lma)
    SAFE_DEALLOCATE(pseudo_start_idx)
  end subroutine
end module opt_variables
