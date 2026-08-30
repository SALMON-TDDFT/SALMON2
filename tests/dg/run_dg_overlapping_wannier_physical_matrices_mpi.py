#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
main_dft=(root/"src/gs/main_dft.f90").read_text()
divided_nonlocal=main_dft.split("subroutine assemble_dg_hybrid_divided_nonlocal_rows",1)[1].split(
    "end subroutine assemble_dg_hybrid_divided_nonlocal_rows",1)[0]
assert "matrix_strength" in divided_nonlocal
assert "action_strength" in divided_nonlocal
assert "apply_dg_overlapping_wannier_nonlocal_action" in divided_nonlocal
assert "local_strength(ilma)=system%hvol*ppg%rinv_uvu(ilma)" not in divided_nonlocal
assert "if(position<=0)then;message='divided basis omits nonlocal projector support';return;endif" not in divided_nonlocal
assert "if(q==0)then;message='complete nonlocal projector payload lacks local identity';ok=.false.;return;endif" not in divided_nonlocal
nonlocal_source=(root/"src/gs/dc/dg_overlapping_wannier_nonlocal.f90").read_text().lower()
collector=nonlocal_source.split("subroutine collect_dg_overlapping_wannier_projector_overlaps",1)[1].split(
    "end subroutine collect_dg_overlapping_wannier_projector_overlaps",1)[0]
assert "mpi_allgatherv(partial_overlap" not in collector
assert "do q=1,p-1" not in collector
with tempfile.TemporaryDirectory(prefix="ow-physical-") as name:
    build=Path(name);(build/"config.h").write_text("")
    env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
    fixtures=[
      ("nonlocal","dg_overlapping_wannier_nonlocal.f90","test_dg_overlapping_wannier_nonlocal_mpi.f90"),
      ("observables","dg_overlapping_wannier_observables.f90","test_dg_overlapping_wannier_observables_mpi.f90"),
    ]
    for label,module,test in fixtures:
      exe=build/label
      subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
        "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
        str(root/"src/gs/dc"/module),str(root/"tests/dg"/test),"-o",str(exe)],check=True)
      signatures=[]
      for n in (1,2,4,8):
        p=subprocess.run([shutil.which("mpiexec"),"-n",str(n),str(exe)],capture_output=True,text=True,env=env)
        assert p.returncode==0,(label,n,p.stdout,p.stderr)
        assert f"PASS overlapping-Wannier {label} on {n} ranks" in p.stdout
        marker="NONLOCAL" if label=="nonlocal" else "OBSERVABLES"
        match=re.search(rf"{marker} ranks=\d+ values=([^\n]+)",p.stdout)
        assert match,p.stdout
        signatures.append([float(value) for value in match.group(1).split()])
      reference=signatures[0]
      assert all(max(abs(a-b) for a,b in zip(reference,row))<1e-13 for row in signatures[1:])
print("PASS overlapping-Wannier physical matrices on 1, 2, 4, and 8 ranks")
