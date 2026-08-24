#!/usr/bin/env python3
from pathlib import Path
import os,shutil,subprocess,tempfile

root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="hybrid-production-pw-") as name:
    build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_production_pw"
    subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-std=f2008","-I",str(build),"-J",str(build),
      "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
      str(root/"src/common/dg_hybrid_windowed_pw_types.f90"),
      str(root/"src/common/dg_hybrid_reciprocal_catalog.f90"),
      str(root/"src/common/dg_hybrid_window_distribution.f90"),
      str(root/"src/common/dg_hybrid_windowed_pw_basis.f90"),
      str(root/"src/common/dg_hybrid_production_pw_basis.f90"),
      str(root/"tests/dg/test_dg_hybrid_production_pw_basis_mpi.f90"),"-o",str(exe)],check=True)
    env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
    for nrank in (1,2,4,8):
        run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env)
        assert run.returncode==0,(nrank,run.stdout,run.stderr)
        assert f"PASS hybrid production PW basis on {nrank} ranks" in run.stdout
print("PASS hybrid production PW basis on 1, 2, 4, and 8 ranks")
