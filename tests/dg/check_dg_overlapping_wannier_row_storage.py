#!/usr/bin/env python3
"""Source-level production storage contract for overlapping-Wannier matrices."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]


def source(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8").lower()


operators = source("src/gs/dc/dg_overlapping_wannier_operators.f90")
metric = source("src/gs/dc/dg_overlapping_wannier_metric.f90")
checkpoint = source("src/gs/dc/dg_overlapping_wannier_checkpoint.f90")
density = source("src/gs/dc/dg_overlapping_wannier_density.f90")
scf = source("src/gs/dc/dg_overlapping_wannier_scf.f90")
solver = source("src/gs/dc/dg_overlapping_wannier_solver.f90")
main = source("src/gs/main_dft.f90")

assert "assemble_dg_overlapping_wannier_weak_operator_rows" in operators, (
    "RED: weak operators do not expose a row-owned production interface"
)

row_routine = operators.split(
    "subroutine assemble_dg_overlapping_wannier_weak_operator_rows", 1
)[1].split("end subroutine", 1)[0]
assert re.search(r"\brow_ids\s*\(\s*:\s*\)", row_routine), (
    "row-owned weak assembly must consume explicit global row IDs"
)
assert not re.search(
    r"allocate\s*\([^)]*\(\s*nwann\s*,\s*nwann\s*\)", row_routine
), "row-owned weak assembly must not allocate a full Nw-by-Nw matrix"
assert not re.search(
    r"mpi_allreduce\s*\([^,\n]+,[^,\n]+,\s*matrix_count\b", row_routine
), "row-owned weak assembly must not all-reduce a full Nw-by-Nw payload"
assert "row_batch_size" in row_routine, (
    "row-owned weak assembly must bound temporary destination-row batches"
)
assert "nwann_min" in row_routine and "expected_min" in row_routine, (
    "row-owned weak assembly must collectively validate scalar contracts"
)

builder = main.split("subroutine ow_build_hamiltonian", 1)[1].split(
    "end subroutine", 1
)[0]
for forbidden in (
    "global_h(:,:)",
    "allocate(global_h",
    "kinetic(:,:)",
    "local_matrix(:,:)",
):
    assert forbidden not in builder, (
        f"production Hamiltonian builder still owns replicated storage: {forbidden}"
    )
assert "assemble_dg_overlapping_wannier_weak_operator_rows" in builder, (
    "production Hamiltonian builder does not use row-owned weak assembly"
)
assert "int(size(hrows),8)*16_8" in re.sub(r"\s+", "", builder), (
    "Hamiltonian byte evidence must derive from the actual row-owned allocation"
)

hermiticity = main.split("subroutine ow_distributed_hermiticity", 1)[1].split(
    "end subroutine", 1
)[0]
assert "row_batch_size" in hermiticity, (
    "distributed Hermiticity diagnostics must use bounded row broadcasts"
)
assert "max(1,maxval(counts)),size(rows,2)" not in re.sub(r"\s+", "", hermiticity), (
    "distributed Hermiticity diagnostics must not allocate the largest owner block"
)

assert "assemble_dg_overlapping_wannier_metric_rows" in metric, (
    "RED: metric does not expose a row-owned production interface"
)
metric_rows = metric.split(
    "subroutine assemble_dg_overlapping_wannier_metric_rows", 1
)[1].split("\n  end subroutine", 1)[0]
assert "row_batch_size" in metric_rows
assert "row_ids(:)" in re.sub(r"\s+", "", metric_rows)
assert "if(rank==0)allocate(root_metric" in re.sub(r"\s+", "", metric_rows), (
    "only rank zero may allocate the temporary full metric"
)
assert "mpi_allreduce(local_metric,metric,matrix_count" not in metric_rows

assert re.search(r"\bow_srows\s*\(\s*:\s*,\s*:\s*\)", main), (
    "production state does not retain row-owned overlap storage"
)
assert not re.search(r"\bow_s\s*\(\s*:\s*,\s*:\s*\)", main), (
    "production state still declares a persistent full overlap"
)

assert "overlap_row_ids(:)" in re.sub(r"\s+", "", checkpoint), (
    "route checkpoint does not record global IDs for overlap row shards"
)
assert "checkpoint%overlap_row_ids" in checkpoint
assert not re.search(
    r"mpi_bcast\s*\(\s*checkpoint%overlap\b", checkpoint
), "route checkpoint must not broadcast a full overlap"
assert "allocate(ow_checkpoint%overlap,source=ow_srows)" in re.sub(r"\s+", "", main)
assert "allocate(ow_checkpoint%overlap_row_ids,source=ow_row_ids)" in re.sub(
    r"\s+", "", main
)
evidence = main.split("subroutine write_ow_ground_state_evidence", 1)[1].split(
    "end subroutine", 1
)[0]
assert "int(size(ow_srows),8)*16_8" in re.sub(r"\s+", "", evidence), (
    "overlap byte evidence must derive from the actual row-owned allocation"
)
assert "replicated_full_matrix_count" not in evidence, (
    "full-matrix absence must be enforced by this source contract, not self-reported"
)
assert "overlap(nwann,nwann)" not in re.sub(r"\s+", "", scf), (
    "SCF must not reconstruct a replicated full overlap"
)
for forbidden in ("density_matrix(nwann,nwann)", "reference_overlap", "matmul(overlap,density_matrix)"):
    assert forbidden not in re.sub(r"\s+", "", density), (
        f"density reconstruction retains a full projected matrix: {forbidden}"
    )
assert "row_ids(:)" in re.sub(r"\s+", "", density)
assert "srows(:,:)" in re.sub(r"\s+", "", density)
for forbidden in ("full(n,n)", "gram(n,n)", "local_gram(n,n)"):
    assert forbidden not in re.sub(r"\s+", "", solver), (
        f"coefficient solver retains a full projected matrix: {forbidden}"
    )
condition_estimator = solver.split("subroutine estimate_metric_condition", 1)[1].split(
    "end subroutine", 1
)[0]
hermiticity_checker = solver.split("subroutine distributed_hermiticity", 1)[1].split(
    "end subroutine", 1
)[0]
assert "size(local_batch)" in hermiticity_checker, (
    "Hermiticity reduction must include the full contiguous batch allocation"
)
assert "mpi_gatherv" in condition_estimator, (
    "metric condition estimate must gather row-owned inverse factors to root"
)
assert "mpi_allreduce" not in condition_estimator.split("squares_inverse=0d0", 1)[1], (
    "metric inverse norm must not use one collective per matrix element"
)

print("overlapping-Wannier row-storage source contract: PASS")
