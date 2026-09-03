#!/usr/bin/env python3
from pathlib import Path
import os,re,shlex,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]

source=(root/"src/common/dg_hybrid_wannier_complement.f90").read_text(errors="replace")
prepare_match=re.search(
  r"\bsubroutine\s+prepare_dg_hybrid_generalized_wannier_metric\b(?P<body>.*?)"
  r"\bend\s+subroutine\s+prepare_dg_hybrid_generalized_wannier_metric\b",
  source,re.IGNORECASE|re.DOTALL)
assert prepare_match,"missing generalized metric preparation routine"
prepare_body=re.sub(r"&\s*\n\s*&?","",prepare_match.group("body").lower())
prepare_body=re.sub(r"\s+","",prepare_body)
assert "negative_limit=roundoff_floor" in prepare_body, \
  "negative Gram rejection must use a roundoff-sized limit independent of the positive rank cutoff"
assert "eigenvalues<-negative_limit" in prepare_body, \
  "negative Gram modes beyond roundoff are not rejected"
assert "eigenvalues<-metric_cutoff" not in prepare_body, \
  "positive metric rank cutoff is still being reused as the negative-mode acceptance limit"

def lapack_libraries():
  if os.environ.get("SALMON_LAPACK_LIBS"):
    return shlex.split(os.environ["SALMON_LAPACK_LIBS"])
  if shutil.which("pkg-config"):
    probe=subprocess.run(["pkg-config","--libs","openblas"],capture_output=True,text=True)
    if probe.returncode==0:return shlex.split(probe.stdout)
    probe=subprocess.run(["pkg-config","--libs","lapack"],capture_output=True,text=True)
    if probe.returncode==0:return shlex.split(probe.stdout)
  return ["-llapack","-lblas"]

with tempfile.TemporaryDirectory(prefix="hybrid-complement-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_complement"
  compile_result=subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI",
    "-ffree-line-length-none","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/common/dg_hybrid_wannier_complement.f90"),
    str(root/"tests/dg/test_dg_hybrid_wannier_complement_mpi.f90"),*lapack_libraries(),"-o",str(exe)],
    capture_output=True,text=True)
  if compile_result.returncode!=0:
    diagnostic=compile_result.stdout+compile_result.stderr
    expected="compute_dg_hybrid_generalized_wannier_projection_tile"
    if expected in diagnostic.lower() and "syntax error" not in diagnostic.lower():
      raise RuntimeError("EXPECTED RED: production is missing "+expected+"\n"+diagnostic)
    raise RuntimeError("hybrid complement compile failed:\n"+diagnostic)
  env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1");fingerprints=[]
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid Wannier complement on {nrank} ranks" in run.stdout
    match=re.search(r"HYBRID_COMPLEMENT ranks=\d+ fingerprint=(-?\d+)",run.stdout);assert match,run.stdout
    fingerprints.append(int(match.group(1)))
  assert len(set(fingerprints))==1,fingerprints
print("PASS hybrid Wannier complement on 1, 2, 4, and 8 ranks")
