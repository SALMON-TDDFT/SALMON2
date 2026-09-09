#!/usr/bin/env python3
from pathlib import Path
import os, shutil, subprocess, tempfile

root=Path(__file__).resolve().parents[2]
source=(root/"src/rt/dg/rt_dg_hybrid_sparse_exchange.f90").read_text().lower()
assert "subroutine exchange_rt_dg_sparse_matrix" in source, "RED: packed coefficient-matrix halo API is absent"
body=source.split("subroutine exchange_rt_dg_sparse_matrix",1)[1].split(
  "end subroutine exchange_rt_dg_sparse_matrix",1)[0]
assert body.count("mpi_alltoallv")==1,"RED: coefficient payload exchange is not one packed Alltoallv"
assert "mpi_bcast" not in body and "mpi_allreduce(mpi_in_place" not in body
main=(root/"src/rt/main_tddft.f90").read_text().lower()
metric_body=main.split("function apply_hybrid_metric_to_coefficients",1)[1].split(
  "end function apply_hybrid_metric_to_coefficients",1)[0]
assert "global_coefficients" not in metric_body,"RED: metric action still replicates R x Nocc coefficients"
assert "exchange_rt_dg_sparse_matrix" in metric_body,"RED: metric action does not use the packed coefficient halo"
invariant_body=main.split("subroutine evaluate_hybrid_rt_physical_invariants",1)[1].split(
  "end subroutine evaluate_hybrid_rt_physical_invariants",1)[0]
assert "global_coefficients" not in invariant_body,"RED: invariant evaluation still replicates R x Nocc coefficients"
with tempfile.TemporaryDirectory(prefix="hybrid-v4-distributed-") as name:
  build=Path(name);(build/"config.h").write_text("");exe=build/"hybrid_v4_distributed"
  flags=[shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
    "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace"]
  test_source=str(root/"tests/dg/test_rt_dg_hybrid_distributed_v4_mpi.f90")
  projection_source=root/"src/rt/dg/rt_dg_hybrid_sparse_projection.f90"
  exchange_source=root/"src/rt/dg/rt_dg_hybrid_sparse_exchange.f90"
  density_source=str(root/"src/rt/dg/rt_dg_hybrid_point_density.f90")
  subprocess.run(flags+[str(exchange_source),density_source,str(projection_source),test_source,"-o",str(exe)],check=True)
  env=os.environ.copy();env["OMP_NUM_THREADS"]="1";env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
  for nrank in (1,2,4,8):
    run=subprocess.run([shutil.which("mpiexec"),"-n",str(nrank),str(exe)],capture_output=True,text=True,env=env,timeout=60)
    assert run.returncode==0,(nrank,run.stdout,run.stderr)
    assert f"PASS distributed-v4 coefficient halo ranks={nrank}" in run.stdout
  mutations=(
    (exchange_source,"send_values((i-1)*nrhs+j)=local_values(plan%send_positions(i),j)",
     "send_values((i-1)*nrhs+j)=conjg(local_values(plan%send_positions(i),j))","packed coefficient permutation"),
    (projection_source,"factor*conjg(support_values(i))*support_values(j)",
     "factor*support_values(i)*support_values(j)","point-CSR bra conjugation"),
  )
  for index,(target,old,new,label) in enumerate(mutations):
    mutated=target.read_text().replace(old,new,1)
    assert mutated!=target.read_text(),f"mutation target missing: {label}"
    mutated_path=build/f"mutation-{index}-{target.name}";mutated_path.write_text(mutated)
    mutation_exe=build/f"mutation-{index}"
    sources=[str(exchange_source),density_source,str(projection_source),test_source]
    sources[0 if target==exchange_source else 2]=str(mutated_path)
    subprocess.run(flags+sources+["-o",str(mutation_exe)],check=True)
    run=subprocess.run([shutil.which("mpiexec"),"-n","2",str(mutation_exe)],capture_output=True,text=True,env=env,timeout=30)
    assert run.returncode!=0 or "PASS distributed-v4 coefficient halo ranks=2" not in run.stdout, \
      f"mutation survived: {label}"
print("PASS distributed-v4 packed coefficient halo on 1, 2, 4, and 8 ranks")
