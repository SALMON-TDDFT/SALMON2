#!/usr/bin/env python3
"""Behavioral sparse-graph/projection oracle and communication mutations."""
from pathlib import Path
import os,shlex,shutil,subprocess,tempfile

root=Path(__file__).resolve().parents[2]
graph_source=root/"src/rt/dg/rt_dg_hybrid_structural_graph.f90"
projection_source=root/"src/rt/dg/rt_dg_hybrid_sparse_projection.f90"
assert graph_source.exists(),"RED: missing exact structural Hybrid operator graph"
assert projection_source.exists(),"RED: missing testable production sparse projection kernel"
projection_text=projection_source.read_text()
partner_scan="do edge=1,total" in projection_text and "do i=1,total" in projection_text
assert not partner_scan,"RED: Hermitian partner validation still performs O(nnz^2) global edge scans"
assert "sort_directed_edges" in projection_text and "find_directed_edge" in projection_text, \
  "RED: Hermitian partner validation lacks deterministic O(nnz log nnz) lookup"

if os.environ.get("SALMON_LAPACK_LIBS"):
  libs=shlex.split(os.environ["SALMON_LAPACK_LIBS"])
elif shutil.which("brew") and subprocess.run(["brew","--prefix","openblas"],capture_output=True).returncode==0:
  prefix=subprocess.check_output(["brew","--prefix","openblas"],text=True).strip();libs=[f"-L{prefix}/lib","-lopenblas"]
else: libs=["-llapack","-lblas"]

def compile_and_run(build,projection,label,expect_success):
  exe=build/f"probe-{label}"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    str(graph_source),str(projection),str(root/"tests/dg/test_rt_dg_hybrid_sparse_projection_mpi.f90"),
    *libs,"-o",str(exe)],check=True,capture_output=True,text=True)
  mutation_detected=False
  for nrank in (1,2,4):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,
      env={**os.environ,"OMP_NUM_THREADS":"1","OMPI_MCA_rmaps_base_oversubscribe":"1"},timeout=45)
    if expect_success:
      assert run.returncode==0,(label,nrank,run.stdout,run.stderr)
      assert f"PASS structural Hybrid sparse projection on {nrank} ranks" in run.stdout
    else:
      mutation_detected = mutation_detected or run.returncode!=0 or \
        f"PASS structural Hybrid sparse projection on {nrank} ranks" not in run.stdout
  if not expect_success: assert mutation_detected,label

with tempfile.TemporaryDirectory(prefix="hybrid-sparse-projection-") as name:
  build=Path(name);(build/"config.h").write_text("")
  compile_and_run(build,projection_source,"baseline",True)
  source=projection_text
  mutations={
    "counts":("call MPI_Allgather(local_edge_count,1,MPI_INTEGER,edge_counts", "call MPI_Allgather(0,1,MPI_INTEGER,edge_counts"),
    "order":("local_edge_rows(row_offsets(p):row_offsets(p+1)-1)=row_ids(p)", "local_edge_rows(row_offsets(p):row_offsets(p+1)-1)=row_ids(size(row_ids)-p+1)"),
    "displacements":("edge_displacements(p)=edge_displacements(p-1)+edge_counts(p-1)", "edge_displacements(p)=0"),
    "permutation":("basis_values(global_edge_columns(edge),p)", "basis_values(global_edge_columns(global_edge_count-edge+1),p)"),
  }
  for label,(old,new) in mutations.items():
    assert source.count(old)==1,(label,source.count(old));mutated=build/f"projection-{label}.f90"
    mutated.write_text(source.replace(old,new,1));compile_and_run(build,mutated,label,False)

print("PASS structural/support and production sparse projection on 1, 2, and 4 ranks; mutations rejected")
