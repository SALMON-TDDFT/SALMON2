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
graph_text=graph_source.read_text()
fixture_text=(root/"tests/dg/test_rt_dg_hybrid_sparse_projection_mpi.f90").read_text()
for label,text in (("graph",graph_text),("projection",projection_text)):
  assert "MPI_Allgatherv" not in text, f"RED: {label} replicates the global sparse edge catalog"
assert "operator_raw_count=nactive*(nactive+1)/2" not in graph_text.replace(" ",""), \
  "RED: structural graph accumulates point-by-point duplicate pairs before unique"
assert "peak_workspace_keys" in graph_text and "integer(int64)" in graph_text, \
  "RED: structural graph has no auditable int64 bounded-workspace receipt"
assert "make_checked_displacements" in graph_text and "integer(int64)::running" in graph_text, \
  "RED: structural graph cumulative MPI extents are not checked in int64"
assert "make_displacements" in projection_text and "integer(int64)::running" in projection_text, \
  "RED: sparse projection cumulative MPI extents are not checked in int64"
assert "partial_values(global_edge_count)" not in projection_text.replace(" ",""), \
  "RED: sparse projection allocates a global-edge-sized replicated buffer"
partner_scan="do edge=1,total" in projection_text and "do i=1,total" in projection_text
assert not partner_scan,"RED: Hermitian partner validation still performs O(nnz^2) global edge scans"
assert "MPI_Alltoallv" in graph_text and "MPI_Alltoallv" in projection_text, \
  "RED: sparse support/projection is not routed directly to row owners"
assert "fragment*nproc/8" in fixture_text and "nproc==8" in fixture_text and "nowned==50" in fixture_text, \
  "RED: scaling fixture does not model one 50-basis fragment per rank at eight ranks"
assert "call project_rt_dg_hybrid_sparse_edges(comm,n,rows" in fixture_text, \
  "RED: production sparse projection is not exercised by the 400-basis scaling fixture"
assert "if(global_bad/=0)then;ierr=-1;return;endif" in projection_text, \
  "RED: owner-only missing-edge failures are not propagated collectively"

if os.environ.get("SALMON_LAPACK_LIBS"):
  libs=shlex.split(os.environ["SALMON_LAPACK_LIBS"])
elif shutil.which("brew") and subprocess.run(["brew","--prefix","openblas"],capture_output=True).returncode==0:
  prefix=subprocess.check_output(["brew","--prefix","openblas"],text=True).strip();libs=[f"-L{prefix}/lib","-lopenblas"]
else: libs=["-llapack","-lblas"]

def compile_and_run(build,projection,label,expect_success,graph=graph_source):
  exe=build/f"probe-{label}"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    str(graph),str(projection),str(root/"tests/dg/test_rt_dg_hybrid_sparse_projection_mpi.f90"),
    *libs,"-o",str(exe)],check=True,capture_output=True,text=True)
  mutation_detected=False
  for nrank in (1,2,4,8):
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
    "counts":("count=merge(size(keys),0,rank==sender)", "count=0"),
    "order":("call MPI_Send(keys,count,MPI_INTEGER8,destination-1,1800+sender,comm,ierr)",
      "call MPI_Send(keys(size(keys):1:-1),count,MPI_INTEGER8,destination-1,1800+sender,comm,ierr)"),
    "displacements":("running=running+int(counts(p),int64)", "running=0_int64"),
    "permutation":("call MPI_Send(values,count,MPI_DOUBLE_COMPLEX,destination-1,1900+sender,comm,ierr)",
      "call MPI_Send(-values,count,MPI_DOUBLE_COMPLEX,destination-1,1900+sender,comm,ierr)"),
    "missing-collective":("if(global_bad/=0)then;ierr=-1;return;endif",
      "if(.false.)then;ierr=-1;return;endif"),
    "stale-success":("if(ierr/=MPI_SUCCESS.or..not.ok)then;message='sparse projection rows do not have unique owners';return;endif\n    ok=.false.",
      "if(ierr/=MPI_SUCCESS.or..not.ok)then;message='sparse projection rows do not have unique owners';return;endif\n    ok=.true."),
  }
  for label,(old,new) in mutations.items():
    assert source.count(old)==1,(label,source.count(old));mutated=build/f"projection-{label}.f90"
    mutated.write_text(source.replace(old,new,1));compile_and_run(build,mutated,label,False)
  dedup_old="elseif(set%slot(position)==key)then;return"
  dedup_new="elseif(set%slot(position)==key)then;set%count=set%count+1_int64;return"
  assert graph_text.count(dedup_old)==1
  mutated_graph=build/"graph-duplicate-accumulation.f90"
  mutated_graph.write_text(graph_text.replace(dedup_old,dedup_new,1))
  compile_and_run(build,projection_source,"duplicate-accumulation",False,mutated_graph)

print("PASS structural/support and production sparse projection on 1, 2, 4, and 8 ranks; mutations rejected")
