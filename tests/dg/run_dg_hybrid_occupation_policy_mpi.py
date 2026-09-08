#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile

root=Path(__file__).resolve().parents[2]
occupation_source=(root/"src/gs/occupation.f90").read_text().lower()
assert "use occupation_kernel" in occupation_source
assert "call solve_spectrum_occupations" in occupation_source
assert "subroutine mu2ne" not in occupation_source
with tempfile.TemporaryDirectory(prefix="hybrid-occupation-policy-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_occupation_policy"
  subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-std=f2008","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
    str(root/"src/gs/occupation_kernel.f90"),
    str(root/"src/gs/dc/dg_hybrid_occupation_policy.f90"),
    str(root/"tests/dg/test_dg_hybrid_occupation_policy_mpi.f90"),"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
  fingerprints=[]
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS hybrid occupation policy on {nrank} ranks" in run.stdout
    match=re.search(r"HYBRID_OCCUPATION_POLICY ranks=\d+ fingerprint=(-?\d+)",run.stdout);assert match,run.stdout
    fingerprints.append(int(match.group(1)))
  assert len(set(fingerprints))==1,fingerprints
print("PASS hybrid occupation policy on 1, 2, 4, and 8 ranks")
