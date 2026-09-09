#!/usr/bin/env python3
"""Behavioral sparse-graph/projection oracle and communication mutations."""
from pathlib import Path
import os,signal,shlex,shutil,subprocess,tempfile

root=Path(__file__).resolve().parents[2]
graph_source=root/"src/rt/dg/rt_dg_hybrid_structural_graph.f90"
projection_source=root/"src/rt/dg/rt_dg_hybrid_sparse_projection.f90"
assert graph_source.exists(),"RED: missing exact structural Hybrid operator graph"
assert projection_source.exists(),"RED: missing testable production sparse projection kernel"
projection_text=projection_source.read_text()
graph_text=graph_source.read_text()
fixture_text=(root/"tests/dg/test_rt_dg_hybrid_sparse_projection_mpi.f90").read_text()
owner_route=projection_text.split("subroutine route_contributions_to_row_owners",1)[1].split(
  "end subroutine route_contributions_to_row_owners",1)[0]
for label,text in (("graph",graph_text),("projection",projection_text)):
  assert "MPI_Allgatherv" not in text, f"RED: {label} replicates the global sparse edge catalog"
  for forbidden in ("MPI_Bcast", "MPI_Send", "MPI_Recv", "1700+sender", "1800+sender", "1900+sender"):
    assert forbidden not in text, f"RED: {label} retains rank-quadratic sender scheduling: {forbidden}"
assert "do destination=1,nproc" not in graph_text and "do destination=1,nproc" not in projection_text, \
  "RED: sparse graph/projection rescans once per rank owner"
assert graph_text.count("MPI_Alltoallv") <= 3 and projection_text.count("MPI_Alltoallv") <= 3, \
  "RED: sparse collective count grows with rank count instead of remaining constant"
assert owner_route.count("MPI_Alltoallv")==1, \
  "RED: each sparse local-potential update must use one packed owner-directed MPI_Alltoallv"
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
  "RED: sparse support/projection is not packed and routed directly to row owners"
assert "fragment*nproc/8" in fixture_text and "nproc==8" in fixture_text and "nowned==50" in fixture_text, \
  "RED: scaling fixture does not model one 50-basis fragment per rank at eight ranks"
assert "call project_rt_dg_hybrid_sparse_edges(comm,n,rows" in fixture_text, \
  "RED: production sparse projection is not exercised by the 400-basis scaling fixture"
assert "if(global_bad/=0)ierr=-1" in projection_text, \
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
    command=[shutil.which("mpiexec"),"-n",str(nrank),str(exe)]
    process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True,
      env={**os.environ,"OMP_NUM_THREADS":"1","OMPI_MCA_rmaps_base_oversubscribe":"1"})
    try:
      stdout,stderr=process.communicate(timeout=45)
    except subprocess.TimeoutExpired:
      os.killpg(process.pid,signal.SIGTERM);stdout,stderr=process.communicate(timeout=5)
      if expect_success: raise AssertionError((label,nrank,"MPI test timed out",stdout,stderr))
      mutation_detected=True;break
    if expect_success:
      assert process.returncode==0,(label,nrank,stdout,stderr)
      assert f"PASS structural Hybrid sparse projection on {nrank} ranks" in stdout
    else:
      mutation_detected = mutation_detected or process.returncode!=0 or \
        f"PASS structural Hybrid sparse projection on {nrank} ranks" not in stdout
      if mutation_detected: break
  if not expect_success: assert mutation_detected,label

with tempfile.TemporaryDirectory(prefix="hybrid-sparse-projection-") as name:
  build=Path(name);(build/"config.h").write_text("")
  compile_and_run(build,projection_source,"baseline",True)
  source=projection_text
  mutations={
    "counts":("words=3_int64*int(counts(p),int64)", "words=4_int64*int(counts(p),int64)"),
    "routing":("row=int((keys(q)-1_int64)/int(n,int64))+1;destination=owners(row)",
      "row=int((keys(q)-1_int64)/int(n,int64))+1;destination=1",2),
    "order":("send_payload(base+1)=keys(q)", "send_payload(base+1)=directed_key(1,1,n)"),
    "displacements":("word_counts(p)=int(words);word_displacements(p)=int(running);running=running+words",
      "word_counts(p)=int(words);word_displacements(p)=0;running=running+words"),
    "permutation":("send_payload(base+2)=transfer(real(values(q),real64),0_int64)",
      "send_payload(base+2)=transfer(-real(values(q),real64),0_int64)"),
    "missing-collective":("if(global_bad/=0)ierr=-1", "if(.false.)ierr=-1"),
    "stale-success":("if(ierr/=MPI_SUCCESS.or..not.ok)then;message='sparse projection rows do not have unique owners';return;endif\n    ok=.false.",
      "if(ierr/=MPI_SUCCESS.or..not.ok)then;message='sparse projection rows do not have unique owners';return;endif\n    ok=.true."),
  }
  for label,replacement in mutations.items():
    old,new,*expected=replacement;occurrences=expected[0] if expected else 1
    assert source.count(old)==occurrences,(label,source.count(old));mutated=build/f"projection-{label}.f90"
    mutated.write_text(source.replace(old,new,occurrences));compile_and_run(build,mutated,label,False)
  dedup_old="elseif(set%slot(position)==key)then;return"
  dedup_new="elseif(set%slot(position)==key)then;set%count=set%count+1_int64;return"
  assert graph_text.count(dedup_old)==1
  mutated_graph=build/"graph-duplicate-accumulation.f90"
  mutated_graph.write_text(graph_text.replace(dedup_old,dedup_new,1))
  compile_and_run(build,projection_source,"duplicate-accumulation",False,mutated_graph)

print("PASS structural/support and production sparse projection on 1, 2, 4, and 8 ranks; mutations rejected")
