#!/usr/bin/env python3
from pathlib import Path
import os,shutil,struct,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="hybrid-refinement-receipt-") as name:
    build=Path(name);(build/"config.h").write_text("");exe=build/"refinement_receipt"
    subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-std=f2008","-ffree-line-length-none",
      "-I",str(build),"-J",str(build),"-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
      str(root/"src/common/dg_portable_sha256.f90"),str(root/"src/rt/dg/rt_dg_hybrid_refinement_receipt.f90"),
      str(root/"tests/dg/test_rt_dg_hybrid_refinement_receipt_mpi.f90"),"-o",str(exe)],check=True)
    env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
    for ranks in (1,2,4,8):
        case=build/f"case-{ranks}";case.mkdir();prefix=case/"hybrid.chk.refinement"
        command=[shutil.which("mpiexec"),"-n",str(ranks),str(exe),"roundtrip",str(prefix)]
        run=subprocess.run(command,capture_output=True,text=True,env=env,timeout=30)
        assert run.returncode==0,(ranks,run.stdout,run.stderr)
        assert f"PASS refinement receipt round trip ranks={ranks}" in run.stdout
        manifest=Path(str(prefix)+".manifest");assert manifest.is_file()
        shards=sorted(case.glob("hybrid.chk.refinement.transaction-*.rank-*"),key=lambda p:p.stat().st_mtime_ns)
        assert len(shards)>=ranks;target=shards[-ranks];original=target.read_bytes();target.write_bytes(original[:-1])
        corrupt=subprocess.run([shutil.which("mpiexec"),"-n",str(ranks),str(exe),"read-corrupt",str(prefix)],
          capture_output=True,text=True,env=env,timeout=30)
        assert corrupt.returncode==0,(ranks,corrupt.stdout,corrupt.stderr)
        run=subprocess.run(command,capture_output=True,text=True,env=env,timeout=30);assert run.returncode==0
        shards=sorted(case.glob("hybrid.chk.refinement.transaction-*.rank-*"),key=lambda p:p.stat().st_mtime_ns)
        target=shards[-ranks];payload=bytearray(target.read_bytes());needle=struct.pack("=d",1.25e-6)
        offset=payload.find(needle);assert offset>=0;payload[offset]^=1;target.write_bytes(payload)
        corrupt=subprocess.run([shutil.which("mpiexec"),"-n",str(ranks),str(exe),"read-corrupt",str(prefix)],
          capture_output=True,text=True,env=env,timeout=30)
        assert corrupt.returncode==0,(ranks,corrupt.stdout,corrupt.stderr)
print("PASS authenticated refinement receipt on 1, 2, 4, and 8 ranks")
