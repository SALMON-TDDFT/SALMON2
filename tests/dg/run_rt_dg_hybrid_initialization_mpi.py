#!/usr/bin/env python3
"""Distributed-v4 initializer architecture and collective dense-v3 rejection."""
from pathlib import Path
import os
import re
import shlex
import shutil
import subprocess
import tempfile

root=Path(__file__).resolve().parents[2]
source=(root/"src/rt/dg/rt_dg_hybrid_initialization_v4.f90").read_text()
lower=source.lower()
body=lower.split("subroutine initialize_rt_dg_hybrid_from_checkpoint",1)[1].split(
  "end subroutine initialize_rt_dg_hybrid_from_checkpoint",1)[0]
compact=re.sub(r"\s+|&","",body)
assert "read_rt_dg_hybrid_checkpoint_v4" in body
assert "dense hybrid v3 checkpoint is unsupported; regenerate distributed-native v4" in body
for forbidden in (
  "collect_complex_rows", "compute_construction_projection", "build_certified_rt_state",
  "allocate(full_", "mpi_allgatherv", "mpi_bcast(c_buffer", "mpi_bcast(b_buffer",
):
  assert forbidden not in body, f"formal initializer retains dense-v3 route: {forbidden}"
for token in (
  "payload%metric_offsets", "payload%operator_offsets", "payload%basis_point_offsets",
  "state%coefficients", "apply_rt_dg_sparse_rows_tiled", "reconstruct_rt_dg_point_csr_density",
):
  assert token in body, f"distributed-v4 initializer omits {token}"

with tempfile.TemporaryDirectory(prefix="hybrid-v3-reject-") as name:
  build=Path(name);(build/"config.h").write_text("")
  exe=build/"reject-v3"
  sources=[
    "src/common/dg_hybrid_sparse_metric.f90", "src/common/dg_hybrid_sparse_operators.f90",
    "src/rt/dg/rt_dg_hybrid_sparse_exchange.f90", "src/rt/dg/rt_dg_hybrid_point_density.f90",
    "src/rt/dg/rt_dg_hybrid_checkpoint_v4.f90", "src/rt/dg/rt_dg_hybrid_checkpoint.f90",
    "src/rt/dg/rt_dg_hybrid_structural_graph.f90", "src/rt/dg/rt_dg_hybrid_sparse_projection.f90",
    "src/rt/dg/rt_dg_hybrid_initialization_v4.f90", "tests/dg/test_rt_dg_hybrid_v3_rejection_mpi.f90",
  ]
  if os.environ.get("SALMON_LAPACK_LIBS"):
    libs=shlex.split(os.environ["SALMON_LAPACK_LIBS"])
  elif shutil.which("brew") and subprocess.run(["brew","--prefix","openblas"],capture_output=True).returncode==0:
    prefix=subprocess.check_output(["brew","--prefix","openblas"],text=True).strip();libs=[f"-L{prefix}/lib","-lopenblas"]
  else: libs=["-llapack","-lblas"]
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",*[str(root/s) for s in sources],
    *libs,"-o",str(exe)],check=True,capture_output=True,text=True)
  (build/"absent-dense-v3.chk").write_bytes(b"SALMON_DG_GS001 stale dense payload")
  env={**os.environ,"OMP_NUM_THREADS":"1","OMPI_MCA_rmaps_base_oversubscribe":"1"}
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],cwd=build,env=env,
      capture_output=True,text=True,timeout=30)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS dense v3 early rejection on {nrank} ranks" in run.stdout

print("PASS distributed-v4 initializer contract and collective dense-v3 rejection on 1/2/4/8 ranks")
