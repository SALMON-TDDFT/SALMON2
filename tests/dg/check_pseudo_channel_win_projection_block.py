#!/usr/bin/env python3
from pathlib import Path

src = Path("src/gs/dc/lcfo_flux.f90").read_text()

if "subroutine write_pseudo_channel_projection_block" not in src:
    raise SystemExit("missing pseudo-channel Wannier90 projection-block writer")

base_writer_start = src.find("subroutine write_wannier_base_win_atomic")
base_writer_end = src.find("end subroutine write_wannier_base_win_atomic", base_writer_start)
if base_writer_start < 0 or base_writer_end < 0:
    raise SystemExit("missing atomic Wannier90 .win writer")

base_writer = src[base_writer_start:base_writer_end]
pseudo_branch_start = base_writer.find("else if(io_status == 0 .and. is_pseudo_channel_projection")
pseudo_branch_end = base_writer.find("else if(io_status == 0)", pseudo_branch_start + 1)
if pseudo_branch_start < 0 or pseudo_branch_end < 0:
    raise SystemExit("missing pseudo_channels branch in atomic Wannier90 .win writer")

pseudo_branch = base_writer[pseudo_branch_start:pseudo_branch_end]
if "call write_pseudo_channel_projection_block" not in pseudo_branch:
    raise SystemExit("pseudo_channels .win branch must write symmetry-aware projections")

for token in ["write_sawf_projection_block", "wannier_site_symmetry"]:
    if token not in src:
        raise SystemExit(f"missing token: {token}")

writer_start = src.find("subroutine write_pseudo_channel_projection_block")
writer_end = src.find("end subroutine write_pseudo_channel_projection_block", writer_start)
writer = src[writer_start:writer_end]
if "trim(wannier_site_symmetry) == 'off'" not in writer:
    raise SystemExit("pseudo_channels .win writer must preserve an explicit off path")
if "write_sawf_projection_block" not in writer:
    raise SystemExit("pseudo_channels .win writer must use the tested off/on formatter")
if "pseudo_channel_win_projection_lmax" in writer:
    raise SystemExit("off pseudo_channels .win behavior must not re-enable shell expansion")
