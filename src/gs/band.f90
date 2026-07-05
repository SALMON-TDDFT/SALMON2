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

    !> Suggested number of consecutive bra k-points per block for the streaming
    !> overlap export (calc_kgrid_prod_block).
    !> The old full-nk export replicated ngrid*nspin*no*nk wavefunction buffers
    !> on every rank (O(nk) memory; SIGBUS at large nk).  The streaming variant
    !> only gathers the k-slices referenced by one block of bra k-points: the
    !> block itself plus its (2*ndk+1)^3 neighbor shells.  This routine picks a
    !> block size such that that gathered buffer stays within a fixed budget:
    !>   - plane-aligned blocks (multiples of nk1*nk2): p planes reference at
    !>     most (p + 2*ndk) planes,
    !>   - row-aligned blocks (multiples of nk1): r rows reference at most
    !>     (2*ndk+1)**2 * r rows,
    !>   - generic blocks: n points reference at most (2*ndk+1)**3 * n slices.
    !> Peak memory is ~2x the budget while the gather allreduce is in flight
    !> (result buffer + zero-padded contribution buffer).
    integer function kgrid_prod_block_size(sys, nk1, nk2, nk3, ndk) result(nblk)
        use structures, only: s_dft_system
        implicit none
        type(s_dft_system), intent(in) :: sys
        integer, intent(in) :: nk1, nk2, nk3, ndk

        ! Memory budget (bytes) for the gathered neighbor-wavefunction buffer:
        integer(8), parameter :: budget_bytes = 1073741824_8 ! 1 GiB

        integer(8) :: slice_bytes, budget_slices, prod_bytes_per_k, nblk8, ncap8
        integer :: nk, nsh, nplane

        nk = nk1 * nk2 * nk3
        nsh = 2 * ndk + 1
        nplane = nk1 * nk2

        ! Bytes of one full-grid, all-orbital wavefunction k-slice:
        slice_bytes = 16_8 * int(sys%ngrid, 8) * int(sys%nspin, 8) * int(sys%no, 8)
        budget_slices = budget_bytes / max(slice_bytes, 1_8)

        ! Cap (in bra k-points) from the per-block product-table memory (two
        ! copies live during the k-communicator reduction).  The cap is
        ! applied INSIDE each branch, in whole planes/rows, so it never
        ! breaks the block alignment that the union bounds above rely on:
        prod_bytes_per_k = 16_8 * int(sys%no, 8)**2 * int(nsh, 8)**3
        ncap8 = max(1_8, budget_bytes / max(2_8 * prod_bytes_per_k, 1_8))

        if (budget_slices >= int(nplane, 8) * int(nsh, 8) &
            & .and. ncap8 >= int(nplane, 8)) then
            ! Plane-aligned blocks: p full (ik1,ik2)-planes per block.
            nblk8 = min(budget_slices / int(nplane, 8) - int(2 * ndk, 8), &
                &       int(nk3, 8), ncap8 / int(nplane, 8)) * int(nplane, 8)
        else if (budget_slices >= int(nk1, 8) * int(nsh, 8)**2 &
            & .and. ncap8 >= int(nk1, 8)) then
            ! Row-aligned blocks: r full ik1-rows per block.
            nblk8 = min(budget_slices / (int(nk1, 8) * int(nsh, 8)**2), &
                &       ncap8 / int(nk1, 8)) * int(nk1, 8)
        else
            ! Generic blocks of n consecutive bra k-points:
            nblk8 = min(max(1_8, budget_slices / int(nsh, 8)**3), ncap8)
        end if

        ! (nk is a multiple of both nplane and nk1, so this keeps alignment.)
        nblk8 = min(nblk8, int(nk, 8))
        nblk = max(1, int(nblk8))

        return
    end function kgrid_prod_block_size


    !> Calculate inner-product between bloch states at neightboring k-grid
    !> for the contiguous block of bra k-points ik = ikb_s..ikb_e (global
    !> 1..nk1*nk2*nk3 indices).  Streaming replacement for the old full-nk
    !> calc_kgrid_prod: instead of replicating the entire wavefunction
    !> (ngrid*nspin*no*nk) on every rank, only the k-slices referenced by
    !> this block (bra points and their (2*ndk+1)^3 neighbor shells, with
    !> periodic wrap) are gathered, so memory scales with the block size and
    !> not with the total number of k-points.  On return prod_dk_blk is
    !> replicated on all ranks (allreduce over the k-communicator), exactly
    !> as the old version but restricted to columns ikb_s..ikb_e.
    subroutine calc_kgrid_prod_block(sys, lg, mg, par, wf, nk1, nk2, nk3, ndk, &
        & ikb_s, ikb_e, ik3d_tbl, prod_dk_blk)
        use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
        use pack_unpack, only: copy_data
        use communication, only: comm_summation
        use math_constants, only: zI
        implicit none
        type(s_dft_system), intent(in) :: sys
        type(s_rgrid), intent(in) :: lg, mg
        type(s_parallel_info), intent(in) :: par
        type(s_orbital), intent(in) :: wf
        integer, intent(in) :: nk1, nk2, nk3, ndk
        integer, intent(in) :: ikb_s, ikb_e
        integer, intent(out) :: ik3d_tbl(1:3, nk1*nk2*nk3)
        complex(8), intent(out) :: prod_dk_blk( &
            & sys%no, sys%no, &
            & -ndk:ndk, -ndk:ndk, -ndk:ndk, &
            & ikb_s:ikb_e)

        integer, parameter :: nrep = 1

        integer :: ik
        integer :: ik1, ik2, ik3
        integer :: jk1, jk2, jk3
        integer :: jdk1, jdk2, jdk3
        integer :: nslot
        ! Gathered wavefunction slices (only the k-points this block needs):
        complex(8), allocatable :: zwf_blk(:, :, :, :, :, :)
        ! Periodic k-index tables (as in the old full-nk version):
        integer, allocatable :: ik_tbl(:, :, :)
        integer, allocatable :: ik_rep_tbl(:, :, :, :)
        ! Global k-index <-> gathered slot maps (0 = slice not needed):
        integer, allocatable :: slot_tbl(:)
        integer, allocatable :: slot_ik(:)
        complex(8), allocatable :: prod_dk_tmp(:, :, :, :, :, :)

        allocate(ik_tbl( &
            & -nrep*nk1:(nrep+1)*nk1, &
            & -nrep*nk2:(nrep+1)*nk2, &
            & -nrep*nk3:(nrep+1)*nk3))
        allocate(ik_rep_tbl(1:3, &
            & -nrep*nk1:(nrep+1)*nk1, &
            & -nrep*nk2:(nrep+1)*nk2, &
            & -nrep*nk3:(nrep+1)*nk3))
        allocate(slot_tbl(nk1*nk2*nk3))
        allocate(slot_ik(nk1*nk2*nk3))

        ! Create ik_tbl with periodic boundary condition:
        call create_ik_tbl()

        ! Mark the k-slices referenced by this block (bra + neighbor shells):
        call create_slot_tbl()

        ! Retrieve the wavefunction slices needed by this block:
        call retrieve_block_zwf()

        allocate(prod_dk_tmp( &
            & sys%no, sys%no, &
            & -ndk:ndk, -ndk:ndk, -ndk:ndk, &
            & ikb_s:ikb_e))
        prod_dk_tmp(:,:,:,:,:,:) = 0d0

        ! Calculate production <k,io|k+dk,jo> for the locally-owned part of
        ! the block (the k-distribution ranges are contiguous, so the
        ! intersection with the block is a contiguous ik range):
!$omp parallel do collapse(4) default(none) &
!$omp private(ik,jdk1,jdk2,jdk3,ik1,ik2,ik3,jk1,jk2,jk3) &
!$omp shared(par,ndk,ikb_s,ikb_e,ik3d_tbl,prod_dk_tmp)
        do ik = max(ikb_s, par%ik_s), min(ikb_e, par%ik_e)
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

        ! Summarize results of prod_dk_blk:
        call comm_summation( &
            & prod_dk_tmp, prod_dk_blk, size(prod_dk_blk), &
            & par%icomm_k)

        deallocate(prod_dk_tmp)
        deallocate(zwf_blk)
        deallocate(slot_ik, slot_tbl, ik_rep_tbl, ik_tbl)

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


    !> Enumerate the (deduplicated) set of global k-points whose wavefunction
    !> slices are referenced by this block: every bra point ikb_s..ikb_e and
    !> its (2*ndk+1)^3 neighbor shell (periodic wrap through ik_tbl).  The
    !> enumeration order is deterministic, so every rank builds the same maps.
    subroutine create_slot_tbl()
        implicit none
        integer :: iik, jjk
        integer :: iik1, iik2, iik3
        integer :: iid1, iid2, iid3

        slot_tbl(:) = 0
        slot_ik(:) = 0
        nslot = 0

        do iik = ikb_s, ikb_e
            iik1 = ik3d_tbl(1, iik)
            iik2 = ik3d_tbl(2, iik)
            iik3 = ik3d_tbl(3, iik)
            do iid3 = -ndk, ndk
                do iid2 = -ndk, ndk
                    do iid1 = -ndk, ndk
                        jjk = ik_tbl(iik1+iid1, iik2+iid2, iik3+iid3)
                        if (slot_tbl(jjk) == 0) then
                            nslot = nslot + 1
                            slot_tbl(jjk) = nslot
                            slot_ik(nslot) = jjk
                        end if
                    end do
                end do
            end do
        end do

        return
    end subroutine create_slot_tbl


    !> Retrieve the wavefunction slices referenced by this block by gathering
    !> the distributed data in MPI (allreduce of the zero-padded local
    !> contributions, as in the old full-nk version, but restricted to the
    !> nslot needed k-slices).
    subroutine retrieve_block_zwf()
        implicit none
        integer :: islot, iik
        complex(8), allocatable :: zwf_blk_tmp(:, :, :, :, :, :)

        allocate(zwf_blk( &
            & lg%is(1):lg%ie(1), &
            & lg%is(2):lg%ie(2), &
            & lg%is(3):lg%ie(3), &
            & sys%nspin, sys%no, nslot))
        allocate(zwf_blk_tmp( &
            & lg%is(1):lg%ie(1), &
            & lg%is(2):lg%ie(2), &
            & lg%is(3):lg%ie(3), &
            & sys%nspin, sys%no, nslot))

        zwf_blk_tmp = 0d0

        do islot = 1, nslot
            iik = slot_ik(islot)
            if (par%ik_s <= iik .and. iik <= par%ik_e) then
                call copy_data( &
                    wf%zwf( &
                        & mg%is(1):mg%ie(1), &
                        & mg%is(2):mg%ie(2), &
                        & mg%is(3):mg%ie(3), &
                        & 1:sys%nspin, &
                        & par%io_s:par%io_e, &
                        & iik, par%im_s), &
                    zwf_blk_tmp( &
                        & mg%is(1):mg%ie(1), &
                        & mg%is(2):mg%ie(2), &
                        & mg%is(3):mg%ie(3), &
                        & 1:sys%nspin, &
                        & par%io_s:par%io_e, &
                        & islot))
            end if
        end do

        call comm_summation( &
            & zwf_blk_tmp, zwf_blk, size(zwf_blk), &
            & par%icomm_rko)

        deallocate(zwf_blk_tmp)

        return
    end subroutine retrieve_block_zwf


    !> Calculate dot-product of Bloch WFs between two k-grid points.
    subroutine calc_prod(iik1, iik2, iik3, jjk1, jjk2, jjk3, prod_k)
        implicit none
        integer, intent(in) :: iik1, iik2, iik3
        integer, intent(in) :: jjk1, jjk2, jjk3
        complex(8), intent(out) :: prod_k(sys%no, sys%no)

        integer :: ir1, ir2, ir3
        integer :: iik, iio, jjk, jjo
        integer :: iislot, jjslot
        real(8) :: qi(3), qj(3), r(3)
        complex(8) :: phase( &
            & mg%is(1):mg%ie(1), &
            & mg%is(2):mg%ie(2), &
            & mg%is(3):mg%ie(3), &
            & sys%nspin)

        complex(8) ZDOTC ! from BLAS library

        iik = ik_tbl(iik1, iik2, iik3)
        jjk = ik_tbl(jjk1, jjk2, jjk3)
        iislot = slot_tbl(iik)
        jjslot = slot_tbl(jjk)

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
                    & zwf_blk(:, :, :, :, iio, iislot), 1, &
                    & zwf_blk(:, :, :, :, jjo, jjslot) * phase(:, :, :, :), 1)
            end do
        end do

        return
    end subroutine calc_prod

end subroutine calc_kgrid_prod_block
end module





