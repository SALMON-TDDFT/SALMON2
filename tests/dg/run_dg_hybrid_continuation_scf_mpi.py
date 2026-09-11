#!/usr/bin/env python3
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
obsolete = root / "src/gs/dc/dg_hybrid_continuation_scf.f90"
assert not obsolete.exists(), "obsolete callback continuation module remains"

source = (root / "src/gs/main_dft.f90").read_text().lower()
assert "use dg_hybrid_continuation_scf" not in source, (
    "production still imports the obsolete callback continuation module"
)
assert "run_dg_hybrid_continuation_scf_fixture" not in source, (
    "production still calls the obsolete fixture-only continuation entry point"
)
# Production no longer uses the global continuation driver. Retain the
# standalone schedule tests while its shared controller module is supported.
assert "run_dg_hybrid_concrete_continuation" not in source

with tempfile.TemporaryDirectory(prefix="dg-continuation-schedule-") as name:
    build = Path(name)
    (build / "config.h").write_text("")
    fixture = build / "test_schedule.f90"
    fixture.write_text(r'''program test_schedule
  use mpi
  use dg_hybrid_continuation_controller,only:s_dg_hybrid_stage_schedule,initialize_dg_hybrid_stage_schedule,&
    begin_dg_hybrid_stage_solve,schedule_dg_hybrid_candidate_checks,complete_dg_hybrid_stage_solve
  implicit none
  integer,parameter::nscf=4
  integer::i,ierr,solve_count,action_count,symmetry_count,hermiticity_count
  logical::cheap_candidate,stage_converged,run_solve,run_expensive,refresh_scheduled,final_refresh_performed
  type(s_dg_hybrid_stage_schedule)::schedule
  call MPI_Init(ierr)
  solve_count=0;action_count=0;symmetry_count=0;hermiticity_count=0
  call initialize_dg_hybrid_stage_schedule(nscf,schedule)
  do
    call begin_dg_hybrid_stage_solve(schedule,run_solve,i)
    if(.not.run_solve)exit
    solve_count=solve_count+1;stage_converged=.false.
    cheap_candidate=i>=nscf
    call schedule_dg_hybrid_candidate_checks(schedule,cheap_candidate,run_expensive)
    if(run_expensive)then
      action_count=action_count+1;symmetry_count=symmetry_count+1;hermiticity_count=hermiticity_count+1
      stage_converged=.true.
    endif
    call complete_dg_hybrid_stage_solve(schedule,stage_converged,1d0,refresh_scheduled,final_refresh_performed)
  enddo
  if(solve_count/=nscf+1)error stop 'final refresh consumed the ordinary iteration budget'
  if(action_count/=2.or.symmetry_count/=2.or.hermiticity_count/=2)error stop 'refresh call counts'
  if(.not.final_refresh_performed)error stop 'final refresh was not completed'
  call MPI_Finalize(ierr)
end program test_schedule
''')
    exe = build / "test_schedule"
    subprocess.run([
        shutil.which("mpifort"), "-cpp", "-DUSE_MPI", "-I", str(build), "-J", str(build),
        str(root / "src/gs/dc/dg_hybrid_continuation_controller.f90"), str(fixture), "-o", str(exe),
    ], check=True)
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_rmaps_base_oversubscribe", "1")
    for nrank in (1, 2, 4):
        subprocess.run([shutil.which("mpiexec"), "-n", str(nrank), str(exe)], check=True, env=env)
print("PASS obsolete callback continuation layer removed")
