#!/usr/bin/env python3
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text()

seed = source.split("subroutine write_wannier_seed_files", 1)[1].split(
    "end subroutine write_wannier_seed_files", 1
)[0]
prepare = seed.find("call prepare_sawf_closed_seed_eigensystem")
eig_write = seed.find("open(iunit,file=eigfile")
amn_write = seed.find("call write_wannier_amn_file")
if min(prepare, eig_write, amn_write) < 0 or not prepare < eig_write < amn_write:
    raise SystemExit("closed SAWF eigensystem is not prepared before eig/amn publication")

if "sawf_explicit_basis_active" not in source:
    raise SystemExit("explicit SAWF basis mode is absent")
hpsi = source.split("subroutine hpsi_basis", 1)[1].split("end subroutine hpsi_basis", 1)[0]
if "sawf_explicit_basis_active" not in hpsi or "sawf_explicit_buffer" not in hpsi:
    raise SystemExit("hpsi does not receive the explicit closed buffer basis")
value = source.split("function local_basis_value", 1)[1].split(
    "end function local_basis_value", 1
)[0]
if "sawf_explicit_basis_active" not in value or "sawf_explicit_buffer" not in value:
    raise SystemExit("Flux traces do not use the explicit closed buffer basis")

prepare_body = source.split("subroutine prepare_sawf_closed_seed_eigensystem", 1)[1].split(
    "end subroutine prepare_sawf_closed_seed_eigensystem", 1
)[0]
for token in [
    "call hpsi_basis",
    "call calc_hamiltonian_matrix",
    "call diag_eigenexa",
    "call write_sawf_closed_seed_file",
    "nbasis_closed>dc%nstate_frag",
]:
    if token.replace(" ", "") not in prepare_body.replace(" ", ""):
        raise SystemExit(f"closed SAWF eigensystem preparation is missing {token}")

print("PASS SAWF closed basis is wired to a second physical Flux eigensystem")
