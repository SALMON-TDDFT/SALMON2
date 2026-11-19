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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module band
    implicit none
  
  contains

    !> Calculate inner-product between bloch states at neightboring k-grid.
    subroutine calc_kgrid_prod(sys, lg, mg, par, wf, nk1, nk2, nk3, ndk, ik3d_tbl, prod_dk)
        use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
        use pack_unpack, only: copy_data
        use communication, only: comm_summation, comm_is_root
        use math_constants, only: zI
        implicit none
        type(s_dft_system), intent(in) :: sys
        type(s_rgrid), intent(in) :: lg, mg
        type(s_parallel_info), intent(in) :: par
        type(s_orbital), intent(in) :: wf
        integer, intent(in) :: nk1, nk2, nk3, ndk
        integer, intent(out) :: ik3d_tbl(1:3, nk1*nk2*nk3)
        complex(8), intent(out) :: prod_dk( &
            & sys%no, sys%no, &
            & -ndk:ndk, -ndk:ndk, -ndk:ndk, &
            & nk1*nk2*nk3)

        integer, parameter :: nrep = 1

        integer :: ik
        integer :: ik1, ik2, ik3
        integer :: jk1, jk2, jk3
        integer :: jdk1, jdk2, jdk3
        complex(8) :: zwf_all( &
            & lg%is(1):lg%ie(1), &
            & lg%is(2):lg%ie(2), &
            & lg%is(3):lg%ie(3), &
            & sys%nspin, sys%no, sys%nk)
        integer :: ik_tbl( &
            & -nrep*nk1:(nrep+1)*nk1, &
            & -nrep*nk2:(nrep+1)*nk2, &
            & -nrep*nk3:(nrep+1)*nk3)
        integer :: ik_rep_tbl(1:3, &
            & -nrep*nk1:(nrep+1)*nk1, &
            & -nrep*nk2:(nrep+1)*nk2, &
            & -nrep*nk3:(nrep+1)*nk3)
        complex(8) :: prod_dk_tmp( &
            & lbound(prod_dk,1):ubound(prod_dk,1), &
            & lbound(prod_dk,2):ubound(prod_dk,2), &
            & lbound(prod_dk,3):ubound(prod_dk,3), &
            & lbound(prod_dk,4):ubound(prod_dk,4), &
            & lbound(prod_dk,5):ubound(prod_dk,5), &
            & lbound(prod_dk,6):ubound(prod_dk,6))
        
        prod_dk_tmp(:,:,:,:,:,:) = 0d0

        ! Create ik_tbl with periodic boundary condition:
        call create_ik_tbl()

        ! Retrieve entire wf:
        call retrieve_entire_zwf()

        ! Calculate production <k,io|k+dk,jo> for all k,io,jo:
!$omp parallel do collapse(4) default(none) &
!$omp private(ik,jdk1,jdk2,jdk3,ik1,ik2,ik3,jk1,jk2,jk3) &
!$omp shared(par,ndk,ik3d_tbl,prod_dk_tmp)
        do ik = par%ik_s, par%ik_e
            do jdk3 = -ndk, ndk
                do jdk2 = -ndk, ndk
                    do jdk1 = -ndk, ndk
                        ! i-th k-point grid: ik
                        ik1 = ik3d_tbl(1, ik)
                        ik2 = ik3d_tbl(2, ik)
                        ik3 = ik3d_tbl(3, ik)
                        ! Neighboring k-point: jk
                        jk1 = ik1 + jdk1
                        jk2 = ik2 + jdk2
                        jk3 = ik3 + jdk3
                        
                        call calc_prod( &
                            & ik1, ik2, ik3, &
                            & jk1, jk2, jk3, & 
                            & prod_dk_tmp(:, :, jdk1, jdk2, jdk3, ik))
                    end do
                end do
            end do
        end do
!$omp end parallel do

        ! Summarize results of prod_dk:
        call comm_summation( &
            & prod_dk_tmp, prod_dk, size(prod_dk), &
            & par%icomm_k)

        return

    contains

    !> Create ik3d_tbl and ik_tbl which represents 3D coordinates to 1D index in k-grid.
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
                                ik_rep_tbl(1, ik1_r, ik2_r, ik3_r) = ir1
                                ik_rep_tbl(2, ik1_r, ik2_r, ik3_r) = ir2
                                ik_rep_tbl(3, ik1_r, ik2_r, ik3_r) = ir3
                            end do
                        end do
                    end do
                end do
            end do
        end do

        return
    end subroutine create_ik_tbl


    !> Retrieve entire wavefunction by gathering the distributed data in MPI.
    subroutine retrieve_entire_zwf()
        implicit none
        integer :: io
        integer, parameter :: im = 1, ispin = 1
        complex(8) :: zwf_all_tmp( &
            & lbound(zwf_all,1):ubound(zwf_all,1), &
            & lbound(zwf_all,2):ubound(zwf_all,2), &
            & lbound(zwf_all,3):ubound(zwf_all,3), &
            & lbound(zwf_all,4):ubound(zwf_all,4), &
            & lbound(zwf_all,5):ubound(zwf_all,5), &
            & lbound(zwf_all,6):ubound(zwf_all,6))

        zwf_all_tmp = 0d0
        
        call copy_data( &
            wf%zwf( &
                & mg%is(1):mg%ie(1), &
                & mg%is(2):mg%ie(2), &
                & mg%is(3):mg%ie(3), &
                & 1:sys%nspin, &
                & par%io_s:par%io_e, &
                & par%ik_s:par%ik_e, &
                & par%im_s), &
            zwf_all_tmp( &
                & mg%is(1):mg%ie(1), &
                & mg%is(2):mg%ie(2), &
                & mg%is(3):mg%ie(3), &
                & 1:sys%nspin, &
                & par%io_s:par%io_e, &
                & par%ik_s:par%ik_e))
        
        call comm_summation( &
            & zwf_all_tmp, zwf_all, size(zwf_all), &
            & par%icomm_rko)

        return
    end subroutine retrieve_entire_zwf


    !> Calculate dot-product of Bloch WFs between two k-grid points.
    subroutine calc_prod(iik1, iik2, iik3, jjk1, jjk2, jjk3, prod_k)
        implicit none
        integer, intent(in) :: iik1, iik2, iik3
        integer, intent(in) :: jjk1, jjk2, jjk3
        complex(8), intent(out) :: prod_k(sys%no, sys%no)

        integer :: ir1, ir2, ir3
        integer :: iik, iio, jjk, jjo
        real(8) :: qi(3), qj(3), r(3)
        complex(8) :: phase( &
            & mg%is(1):mg%ie(1), &
            & mg%is(2):mg%ie(2), &
            & mg%is(3):mg%ie(3), &
            & sys%nspin)

        complex(8) ZDOTC ! from BLAS library

        iik = ik_tbl(iik1, iik2, iik3)
        jjk = ik_tbl(jjk1, jjk2, jjk3)

        ! Calculate phase shift at j-th k-point coordinate:
        qi(1:3) = ik_rep_tbl(1, iik1, iik2, iik3) * sys%primitive_b(1:3, 1) &
            &   + ik_rep_tbl(2, iik1, iik2, iik3) * sys%primitive_b(1:3, 2) &
            &   + ik_rep_tbl(3, iik1, iik2, iik3) * sys%primitive_b(1:3, 3) 
        qj(1:3) = ik_rep_tbl(1, jjk1, jjk2, jjk3) * sys%primitive_b(1:3, 1) &
            &   + ik_rep_tbl(2, jjk1, jjk2, jjk3) * sys%primitive_b(1:3, 2) &
            &   + ik_rep_tbl(3, jjk1, jjk2, jjk3) * sys%primitive_b(1:3, 3) 

        ! Calculate phase factor:
        do ir3 = mg%is(3), mg%ie(3)
            do ir2 = mg%is(2), mg%ie(2)
                do ir1 = mg%is(1), mg%ie(1)
                    ! Grid-point coordinate in cartesian: r
                    r(1:3) =  ir1 * (sys%primitive_a(1:3, 1) / mg%num(1)) &
                        &  +  ir2 * (sys%primitive_a(1:3, 2) / mg%num(2)) &
                        &  +  ir3 * (sys%primitive_a(1:3, 3) / mg%num(3)) 
                    ! Phase factor:
                    phase(ir1, ir2, ir3, :) = exp(zI * dot_product(qi - qj, r))
                end do
            end do
        end do

        do jjo = 1, sys%no
            do iio = 1, sys%no
                ! Compute dot-product of two eigenstate u(ik,io) and u(jk,jo):
                ! <u(ik,io)|u(jk,jo)> = <u(ik-iq,io)|u(jk-jq,jo)*exp(zI(iq-jq)r)> 
                prod_k(iio, jjo) = sys%Hvol * ZDOTC( &
                    & sys%ngrid * sys%nspin, &
                    & zwf_all(:, :, :, :, iio, iik), 1, &
                    & zwf_all(:, :, :, :, jjo, jjk) * phase(:, :, :, :), 1)
            end do
        end do

        return
    end subroutine calc_prod

end subroutine calc_kgrid_prod
end module




  
