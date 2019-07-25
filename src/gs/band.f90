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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module band
    implicit none
  
  contains

    !> Calculate inner-product to orbits in neightboring k-grid.
    subroutine calc_kgrid_prod(sys, lg, mg, wfi, wf, nk1, nk2, nk3, ndk, ik3d_tbl, prod_dk)
        use structures, only: s_dft_system, s_rgrid, s_orbital_parallel, s_orbital
        use salmon_parallel, only: nproc_group_global, nproc_id_global
        use math_constants, only: zI
        implicit none
        type(s_dft_system), intent(in) :: sys
        type(s_rgrid), intent(in) :: lg, mg
        type(s_orbital_parallel), intent(in) :: wfi
        type(s_orbital), intent(in) :: wf
        integer, intent(in) :: nk1, nk2, nk3, ndk
        integer, intent(out) :: ik3d_tbl(1:3, nk1*nk2*nk3)
        complex(8), intent(out) :: prod_dk( &
            & nk1*nk2*nk3, &
            & -ndk:ndk, -ndk:ndk, -ndk:ndk, &
            & sys%no, sys%no)
        
        integer, parameter :: nrep = 1

        integer :: im, ik, jk
        integer :: ik1, ik2, ik3
        integer :: jdk1, jdk2, jdk3
        complex(8) :: zwf_all( &
            & lg%is(1):lg%ie(1), &
            & lg%is(2):lg%ie(2), &
            & lg%is(3):lg%ie(3), &
            & sys%nspin, sys%no, sys%nk)
        integer :: ik_tbl( &
            & (1-nrep*nk1):(nrep+1)*nk1, &
            & (1-nrep*nk2):(nrep+1)*nk2, &
            & (1-nrep*nk3):(nrep+1)*nk3)

        ! Create ik_tbl with periodic boundary condition:
        call create_ik_tbl()

        ! Retrieve entire wf:
        call retrieve_entire_zwf()

        ! Calculate production <k,io|k+dk,jo> for all k,io,jo:
        !$omp parallel do collapse(6) default(none) &
        !$omp private(ik1, ik2, ik3, jdk1, jdk2, jdk3, ik, jk) &
        !$omp shared(nk1, nk2, nk3, ndk, ik_tbl, prod_dk)
        do ik3 = 1, nk3
            do ik2 = 1, nk2
                do ik1 = 1, nk1
                    do jdk3 = -ndk, ndk
                        do jdk2 = -ndk, ndk
                            do jdk1 = -ndk, ndk
                                ik = ik_tbl(ik1, ik2, ik3)
                                jk = ik_tbl(ik1+jdk1, ik2+jdk2, ik3+jdk3)
                                call calc_prod(ik, jk, prod_dk(ik, jdk1, jdk2, jdk3, :, :))
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        return

    contains

    subroutine create_ik_tbl()
        implicit none
        integer :: ik_count
        integer :: ik1_o, ik2_o, ik3_o
        integer :: ik1_r, ik2_r, ik3_r
        integer :: ir1, ir2, ir3

        ik_count = 0

        do ik3_o = 1, nk3
            do ik2_o = 1, nk2
                do ik1_o = 1, nk1
                    ik_count = ik_count + 1
                    ! Assign to ik3d_tbl:
                    ik3d_tbl(1, ik_count) = ik1_o
                    ik3d_tbl(2, ik_count) = ik2_o
                    ik3d_tbl(3, ik_count) = ik3_o
                    ! Assign to ik_tbl with replicated coordinates:
                    do ir3 = -nrep, nrep
                        ik3_r = ik3_o + ir3 * nk3
                        do ir2 = -nrep, nrep
                            ik2_r = ik2_o + ir2 * nk2
                            do ir1 = -nrep, nrep
                                ik1_r = ik1_o + ir1 * nk1
                                ik_tbl(ik1_r, ik2_r, ik3_r) = ik_count
                            end do
                        end do
                    end do
                end do
            end do
        end do

        return
    end subroutine create_ik_tbl


    subroutine retrieve_entire_zwf()
        use pack_unpack, only: copy_data
        use salmon_communication, only: comm_summation, comm_is_root
      
        implicit none
        integer :: io
        integer, parameter :: im = 1, ispin = 1
        complex(8) :: zwf_all_tmp( &
            & lg%is(1):lg%ie(1), &
            & lg%is(2):lg%ie(2), &
            & lg%is(3):lg%ie(3), &
            & sys%nspin, sys%no, sys%nk)

        zwf_all_tmp = 0d0
        
        call copy_data( &
            wf%zwf( &
                & mg%is(1):mg%ie(1), &
                & mg%is(2):mg%ie(2), &
                & mg%is(3):mg%ie(3), &
                & 1:sys%nspin, &
                & wfi%io_s:wfi%io_e, &
                & wfi%ik_s:wfi%ik_e, &
                & wfi%im_s), &
            zwf_all_tmp( &
                & mg%is(1):mg%ie(1), &
                & mg%is(2):mg%ie(2), &
                & mg%is(3):mg%ie(3), &
                & 1:sys%nspin, &
                & wfi%io_s:wfi%io_e, &
                & wfi%ik_s:wfi%ik_e))
        
        call comm_summation( &
            & zwf_all_tmp, zwf_all, &
            & sys%ngrid*sys%nspin*sys%no*sys%nk, &
            & wfi%icomm_rko)

        return
    end subroutine retrieve_entire_zwf


    subroutine calc_prod(iik, jjk, prod_ij)
        implicit none
        integer, intent(in) :: iik, jjk
        complex(8), intent(out) :: prod_ij(sys%no, sys%no)

        integer :: ir1, ir2, ir3
        integer :: iio, jjo
        real(8) :: dk_ij(3), r(3), r_red(3)
        complex(8) :: phase_ij( &
            & mg%is(1):mg%ie(1), &
            & mg%is(2):mg%ie(2), &
            & mg%is(3):mg%ie(3), &
            & sys%nspin)

        complex(8) ZDOTC ! From BLAS

        ! Bloch wavenumber difference: dk_ij = ki - kj:
        dk_ij(1:3) = sys%vec_k(1:3, iik) - sys%vec_k(1:3, jjk)
        
        do ir3 = mg%is(3), mg%ie(3)
            do ir2 = mg%is(2), mg%ie(2)
                do ir1 = mg%is(1), mg%ie(1)
                    ! ii-th grid-point in cartesian coordinate: r
                    r(1:3) = sys%primitive_a(1:3, 1) * sys%hgs(1) * ir1 &
                        &  + sys%primitive_a(1:3, 2) * sys%hgs(2) * ir2 &
                        &  + sys%primitive_a(1:3, 3) * sys%hgs(3) * ir3 
                    ! phase factor: exp(-i(ki-kj)r)
                    phase_ij(ir1, ir2, ir3, :) = exp(-zI * dot_product(dk_ij, r))
                end do
            end do
        end do

        do jjo = 1, sys%no
            do iio = 1, sys%no
                ! Compute dot-product of two eigenstate psi(ik,io) and psi(jk,jo):
                ! <psi(ik,io)|psi(jk,jo)> = <u(ik,io)|exp(-i(ki-kj))*u(jk,j)>
                prod_ij(iio, jjo) = sys%Hvol * ZDOTC( &
                    & sys%ngrid * sys%nspin, &
                    & zwf_all(:, :, :, :, iio, iik), 1, &
                    & zwf_all(:, :, :, :, jjo, jjk) * phase_ij(:, :, :, :), 1)
            end do
        end do

        return
    end subroutine calc_prod


    end subroutine calc_kgrid_prod
end module














  
