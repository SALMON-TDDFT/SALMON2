#!/usr/bin/env python3
"""Check occupied-state persistence and formal v4 publication preconditions."""
from pathlib import Path
import os, re, shutil, subprocess, tempfile

root=Path(__file__).resolve().parents[2]
source_path=root/"src/rt/dg/rt_dg_hybrid_checkpoint.f90"
source=source_path.read_text().lower()
main=(root/"src/gs/main_dft.f90").read_text().lower()
continuation=main.split("subroutine run_dg_hybrid_concrete_continuation",1)[1].split(
    "end subroutine run_dg_hybrid_concrete_continuation",1)[0]
publisher=main.split("subroutine publish_dg_hybrid_divided_v4",1)[1].split(
    "end subroutine publish_dg_hybrid_divided_v4",1)[0]

for forbidden in (
    "s_rt_dg_hybrid_ground_state_payload", "write_rt_dg_hybrid_ground_state_checkpoint",
    "read_rt_dg_hybrid_ground_state_checkpoint", "ground_state_version", "full_coefficients",
):
    assert forbidden not in source, f"obsolete dense checkpoint implementation remains: {forbidden}"
assert "call publish_dg_hybrid_divided_v4" in continuation
assert "collect_dg_hybrid_full_rows" not in continuation
assert "write_rt_dg_hybrid_ground_state_checkpoint" not in continuation
assert continuation.index("if(.not.final_refresh_performed)") < continuation.index(
    "call publish_dg_hybrid_divided_v4")
assert continuation.count("call publish_dg_hybrid_divided_v4")==1
assert publisher.count("call publish_rt_dg_hybrid_checkpoint_v4")==1
assert "call write_rt_dg_hybrid_checkpoint_v4" not in publisher
assert "collect_dg_hybrid_full_rows" not in publisher and "full_metric" not in publisher
assert re.search(r"occupied_version\s*=\s*2",source)
for required in ("collective_rt_dg_hybrid_publication_precondition",
                 "collective_rt_dg_hybrid_publication_mapping_precondition",
                 "write_rt_dg_hybrid_occupied_checkpoint",
                 "read_rt_dg_hybrid_occupied_checkpoint"):
    assert required in source

with tempfile.TemporaryDirectory(prefix="hybrid-v4-publication-") as name:
    build=Path(name);(build/"config.h").write_text("")
    precondition=build/"publication";occupied=build/"occupied";checkpoint=build/"occupied.chk"
    common=[shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
            "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
            str(root/"src/rt/dg/rt_dg_hybrid_checkpoint_v4.f90"),str(source_path)]
    subprocess.run([*common,str(root/"tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90"),"-o",str(precondition)],check=True)
    subprocess.run([*common,str(root/"tests/dg/test_rt_dg_hybrid_occupied_checkpoint_mpi.f90"),"-o",str(occupied)],check=True)
    env={**os.environ,"OMP_NUM_THREADS":"1","OMPI_MCA_rmaps_base_oversubscribe":"1"}
    for nrank in (1,2,4,8):
        good=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(precondition),"valid"],
                            capture_output=True,text=True,env=env,timeout=30)
        assert good.returncode==0,(nrank,good.stdout,good.stderr)
        assert f"PASS v4 publication preconditions on {nrank} ranks" in good.stdout
        for mode,diagnostic in (
            ("bad_local","PASS named collective v4 publication precondition rejection"),
            ("bad_mapping","PASS named collective v4 mapping rejection before publication"),
        ):
            if mode=="bad_mapping" and nrank==1: continue
            bad=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(precondition),mode],
                               capture_output=True,text=True,env=env,timeout=30)
            assert bad.returncode!=0,(nrank,mode,bad.stdout,bad.stderr)
            assert diagnostic in bad.stdout+bad.stderr,(nrank,mode,bad.stdout,bad.stderr)
        written=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(occupied),"write",str(checkpoint)],
                               capture_output=True,text=True,env=env,timeout=30)
        assert written.returncode==0,(nrank,written.stdout,written.stderr)
        read=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(occupied),"read",str(checkpoint)],
                            capture_output=True,text=True,env=env,timeout=30)
        assert read.returncode==0,(nrank,read.stdout,read.stderr)
        assert f"PASS occupied checkpoint on {nrank} ranks" in read.stdout
print("PASS v4 publication and occupied checkpoints on 1, 2, 4, and 8 ranks")
