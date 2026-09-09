#!/usr/bin/env python3
"""Reject dense/global Hybrid GS-to-RT data paths in the distributed-v4 route."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[2]
checkpoint = (root / "src/rt/dg/rt_dg_hybrid_checkpoint.f90").read_text().lower()
checkpoint_v4 = (root / "src/rt/dg/rt_dg_hybrid_checkpoint_v4.f90").read_text().lower()
initialization = (root / "src/rt/dg/rt_dg_hybrid_initialization_v4.f90").read_text().lower()
density = (root / "src/rt/dg/rt_dg_hybrid_density_update.f90").read_text().lower()
exchange = (root / "src/rt/dg/rt_dg_hybrid_sparse_exchange.f90").read_text().lower()
main_gs = (root / "src/gs/main_dft.f90").read_text().lower()
main_rt = (root / "src/rt/main_tddft.f90").read_text().lower()
wf_checkpoint = (root / "src/gs/dc/dg_overlapping_wannier_checkpoint.f90").read_text().lower()
assert not (root / "src/rt/dg/rt_dg_hybrid_initialization.f90").exists(), (
    "RED: dense-v3 initializer implementation remains in the production source tree")
assert "publish_dg_hybrid_divided_v3_legacy_unreachable" not in main_gs, (
    "RED: dead dense-v3 divided publisher remains in production source")

assert re.search(r"schema_version\s*=\s*5", checkpoint_v4), "RED: authenticated distributed Hybrid handoff is not schema v5"
assert "salmon_hybrid_dg_manifest_v5" in checkpoint_v4 and "salmon_hybrid_dg_rank_shard_v5" in checkpoint_v4, (
    "RED: schema-v5 manifest/shard magic is missing")
assert "unsupported distributed checkpoint schema v4" in checkpoint_v4, (
    "RED: legacy schema-v4 input is not rejected with a named migration diagnostic")
for digest in ("system_fingerprint(4)", "pseudopotential_digest(4)", "shard_digest(4)"):
    assert digest in checkpoint_v4, f"RED: schema-v5 does not retain full SHA-256 {digest}"
publisher = main_gs.split("subroutine publish_dg_hybrid_divided_v4", 1)[1].split(
    "end subroutine publish_dg_hybrid_divided_v4", 1)[0]
for forbidden in ("full_coefficients", "collect_dg_hybrid_full_rows", "symmetry_representation(n,n", "real(n,8)*real(n,8)"):
    assert forbidden not in publisher, f"RED: v4 publisher retains dense/global construct {forbidden}"
for required in ("metric_values", "operator_values", "basis_point_offsets", "basis_support_ids",
                 "initial_occupied_amplitudes"):
    assert required in checkpoint_v4, f"RED: distributed-v4 payload lacks {required}"

initializer = initialization.split("subroutine initialize_rt_dg_hybrid_from_checkpoint", 1)[1].split(
    "end subroutine initialize_rt_dg_hybrid_from_checkpoint", 1)[0]
assert "read_rt_dg_hybrid_checkpoint_v4" in initializer,"RED: formal RT initializer does not consume v4"
for forbidden in ("full_u", "full_s", "full_h", "local_matrix(r,r)", "projected_s(r,r)",
                  "global_coefficients(hybrid_state%certified_rank"):
    assert forbidden not in initializer + main_rt, f"RED: formal RT startup/main retains {forbidden}"
assert "apply_rt_dg_sparse_rows_tiled" in exchange, "RED: no unique-row tiled sparse multi-RHS action"
assert "exchange_rt_dg_sparse_matrix" not in exchange, "RED: obsolete edge x Nocc halo remains"
reconstruction = density.split("subroutine reconstruct_rt_dg_hybrid_density", 1)[1].split(
    "end subroutine reconstruct_rt_dg_hybrid_density", 1)[0]
for forbidden in ("do owner=0,nproc-1", "mpi_bcast", "mpi_reduce(", "basis_batch(state%certified_rank"):
    assert forbidden not in reconstruction, f"RED: density reconstruction retains {forbidden}"
for required in ("basis_point_offsets", "basis_support_ids", "reconstruct_rt_dg_point_csr_density"):
    if required == "basis_support_ids": required = "basis_support_slots"
    assert required in reconstruction, f"RED: local point-CSR density path lacks {required}"
projection_callback = main_rt.split("subroutine project_salmon_local_rows",1)[1].split(
    "end subroutine project_salmon_local_rows",1)[0]
assert "project_rt_dg_hybrid_point_csr_edges" in projection_callback, (
    "RED: production local-potential projection still scans dense global basis rows")

writer = checkpoint_v4.split("subroutine write_rt_dg_hybrid_checkpoint_v4", 1)[1].split(
    "end subroutine write_rt_dg_hybrid_checkpoint_v4", 1)[0]
reader = checkpoint_v4.split("subroutine read_rt_dg_hybrid_checkpoint_v4", 1)[1].split(
    "end subroutine read_rt_dg_hybrid_checkpoint_v4", 1)[0]
assert "do owner=0,nproc-1" not in writer + reader, "RED: v4 checkpoint I/O loops over all owners"
assert "manifest" in writer and "shard" in writer, "RED: v4 does not publish rank shards through an atomic manifest"
assert "file_nproc/=nproc" in reader, "RED: v4 reader no longer enforces exact MPI rank reuse"

# The expensive fragment-WF checkpoint remains its established independent shard/manifest format.
for token in ("versioned_shard_name", "manifest", "mlwf_backend"):
    assert token in wf_checkpoint, f"fragment-WF checkpoint contract was accidentally removed: {token}"

print("PASS distributed-native Hybrid v4 static architecture contract")
