#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile

root=Path(__file__).resolve().parents[2]
overlay=Path(os.environ.get("SALMON_GP5_OVERLAY","/tmp/salmon-w4-overlay.OScR7T/build"))
eigen_build=overlay/"eigenexa/src/eigenexa-project-build/src"
assert (eigen_build/"libEigenExa.a").exists(),"clean full-feature EigenExa overlay is required"
with tempfile.TemporaryDirectory(prefix="ow-eigenexa-") as name:
    build=Path(name);(build/"config.h").write_text("#define USE_MPI\n#define USE_EIGENEXA\n")
    common=[shutil.which("mpifort"),"-cpp","-fopenmp","-I",str(build),"-I",str(eigen_build),"-J",str(build)]
    support=build/"support.o";metric=build/"metric.o";construction=build/"construction.o"
    solver=build/"solver.o";test=build/"test.o";exe=build/"eigenexa_solver"
    subprocess.run(common+["-c",str(root/"tests/dg/dg_eigenexa_test_support.f90"),"-o",str(support)],check=True)
    subprocess.run(common+["-c",str(root/"src/gs/dc/dg_overlapping_wannier_metric.f90"),"-o",str(metric)],check=True)
    subprocess.run(common+["-c",str(root/"src/gs/dc/dg_overlapping_wannier_construction.f90"),"-o",str(construction)],check=True)
    subprocess.run(common+["-c",str(root/"src/gs/dc/dg_overlapping_wannier_solver.f90"),"-o",str(solver)],check=True)
    subprocess.run(common+["-c",str(root/"tests/dg/test_dg_overlapping_wannier_eigenexa_mpi.f90"),"-o",str(test)],check=True)
    subprocess.run([shutil.which("mpifort"),"-fopenmp",str(support),str(metric),str(construction),str(solver),str(test),
      str(eigen_build/"libEigenExa.a"),"/opt/homebrew/lib/libscalapack.2.2.3.dylib",
      "-L/opt/homebrew/opt/openblas/lib","-lopenblas","-o",str(exe)],check=True)
    env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1");signatures=[]
    for nproc in (1,2,4,8):
        result=subprocess.run([shutil.which("mpiexec"),"-n",str(nproc),str(exe)],
          capture_output=True,text=True,env=env,timeout=60)
        assert result.returncode==0,(nproc,result.stdout,result.stderr)
        match=re.search(r"EIGENEXA ranks=\d+ signature=(-?\d+)",result.stdout);assert match,result.stdout
        signatures.append(int(match.group(1)))
    assert len(set(signatures))==1,signatures
    averaged_signatures=[]
    for nproc in (1,2,4,8):
      result=subprocess.run([shutil.which("mpiexec"),"-n",str(nproc),str(exe),"average"],
        capture_output=True,text=True,env=env,timeout=60)
      assert result.returncode==0,("average",nproc,result.stdout,result.stderr)
      match=re.search(r"AVERAGE ranks=\d+ signature=(-?\d+)",result.stdout);assert match,result.stdout
      averaged_signatures.append(int(match.group(1)))
    assert len(set(averaged_signatures))==1,averaged_signatures
    for case_name in ("degenerate","nonreal","illmetric","residual"):
      for nproc in (1,2,4,8):
        result=subprocess.run([shutil.which("mpiexec"),"-n",str(nproc),str(exe),case_name],
          capture_output=True,text=True,env=env,timeout=60)
        assert result.returncode==0,(case_name,nproc,result.stdout,result.stderr)
        assert f"REJECT {case_name} ranks={nproc}" in result.stdout,result.stdout
print("PASS generalized EigenExa solver on 1, 2, 4, and 8 ranks")
