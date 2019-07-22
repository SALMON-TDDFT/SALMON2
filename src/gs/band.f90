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
    subroutine calc_kgrid_prod(system, rgrid_lg, rgrid_mg, wf_info, wavefunction, nk1, nk2, nk3, ndk, ik3d_tbl, prod_dk)
        use structures
        implicit none
        type(s_system), intent(in) :: system
        type(s_rgrid), intent(in) :: rgrid_lg, rgrid_mg
        type(s_wf_info), intent(in) :: wf_info
        type(s_wavefunction), intent(in) :: wavefunction
        integer, intent(in) :: nk1, nk2, nk3, ndk
        integer, intent(out) :: ik3d_tbl(1:3, nk1*nk2*nk3)
        complex(8), intent(out) :: prod_dk( &
            & nk1*nk2*nk3, &
            & -ndk:ndk, -ndk:ndk, -ndk:ndk, &
            & system%no, system%no)
        
        integer, parameter :: nrep = 1

        integer :: im, ik, jk
        integer :: ik1, ik2, ik3
        integer :: jdk1, jdk2, jdk3
        complex(8) :: zwf_all( &
            & rgrid_lg%is(1):rgrid_lg%ie(1), &
            & rgrid_lg%is(2):rgrid_lg%ie(2), &
            & rgrid_lg%is(3):rgrid_lg%ie(3), &
            & system%nspin, system%no, system%nk)
        integer :: ik_tbl( &
            & -nrep*nk1:(nrep+1)*nk1, &
            & -nrep*nk2:(nrep+1)*nk2, &
            & -nrep*nk3:(nrep+1)*nk3)

        ! Create ik_tbl with periodic boundary condition:
        call create_ik_tbl()

        ! Retrieve entire wavefunction:
        call retrieve_entire_zwf()

        ! Calculate production <k,io|k+dk,jo> for all k,io,jo:
        do ik3 = 1, nk3
            do ik2 = 1, nk2
                do ik1 = 1, nk1
                    ik = ik_tbl(ik1, ik2, ik3)
                    do jdk3 = -ndk, ndk
                        do jdk2 = -ndk, ndk
                            do jdk1 = -ndk, ndk
                                jk = ik_tbl(ik1+jdk1, ik2+jdk2, ik3+jdk3)
                                call calc_prod(ik, jk, prod_dk(ik, jdk1, jdk2, jdk3, :, :))
                            end do
                        end do
                    end do
                end do
            end do
        end do

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
        use salmon_parallel, only: nproc_group_global, nproc_id_global
      
        implicit none
        integer :: io
        integer, parameter :: im = 1, ispin = 1
        complex(8) :: zwf_all_tmp( &
            & rgrid_lg%is(1):rgrid_lg%ie(1), &
            & rgrid_lg%is(2):rgrid_lg%ie(2), &
            & rgrid_lg%is(3):rgrid_lg%ie(3), &
            & system%nspin, system%no, system%nk)

        zwf_all_tmp = 0d0
        
        call copy_data( &
            wavefunction%zwf( &
                & rgrid_mg%is(1):rgrid_mg%ie(1), &
                & rgrid_mg%is(2):rgrid_mg%ie(2), &
                & rgrid_mg%is(3):rgrid_mg%ie(3), &
                & 1:system%nspin, &
                & wf_info%io_s:wf_info%io_e, &
                & wf_info%ik_s:wf_info%ik_e, &
                & wf_info%im_s), &
            zwf_all_tmp( &
                & rgrid_mg%is(1):rgrid_mg%ie(1), &
                & rgrid_mg%is(2):rgrid_mg%ie(2), &
                & rgrid_mg%is(3):rgrid_mg%ie(3), &
                & 1:system%nspin, &
                & wf_info%io_s:wf_info%io_e, &
                & wf_info%ik_s:wf_info%ik_e))
        
        call comm_summation( &
            & zwf_all_tmp, zwf_all, &
            & system%ngrid*system%nspin*system%no*system%nk, &
            & wf_info%icomm_rko)

        return
    end subroutine retrieve_entire_zwf


    subroutine calc_prod(iik, jjk, prod_ij)
        implicit none
        integer, intent(in) :: iik, jjk
        complex(8), intent(out) :: prod_ij(system%no, system%no)
        integer :: iio, jjo

        complex(8) ZDOTC ! From BLAS

        do jjo = 1, system%no
            do iio = 1, system%no
                ! Compute dot-products: <iik,iio|jjk,jjo>
                prod_ij(iio, jjo) = system%Hvol * ZDOTC( &
                    & system%ngrid * system%nspin, &
                    & zwf_all(:, :, :, :, iio, iik), 1, &
                    & zwf_all(:, :, :, :, jjo, jjk), 1)
            end do
        end do

        return
    end subroutine calc_prod


    end subroutine calc_kgrid_prod
end module














  