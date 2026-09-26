#!/usr/bin/env python3
"""Sequential MPI2 native integration test. Needs a USE_HSE binary and H_rps.dat.
Usage: python3 test_native.py --binary /abs/salmon --pseudo /abs/H_rps.dat
Each run owns a fresh temporary directory; no production data is modified.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--binary", type=Path, required=True)
parser.add_argument("--pseudo", type=Path, required=True)
parser.add_argument("--mpirun", default="mpirun")
parser.add_argument("--work-root", type=Path, default=Path(tempfile.gettempdir()))
args = parser.parse_args()
root = Path(tempfile.mkdtemp(prefix="lcfo-native-test-", dir=args.work_root))
repo = Path(__file__).resolve().parents[2]
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
env.pop("SALMON_LCFO_RT", None)
command = [args.mpirun, "-np", "2", str(args.binary.resolve())]

def run(name, text, rt=False, reject=False, extra_env=None):
    folder = root / name
    folder.mkdir()
    (folder / "inputfile").write_text(text)
    shutil.copy(args.pseudo, folder / "H_rps.dat")
    if rt:
        (folder / "data_dcdft").symlink_to(root / "gs/data_dcdft", target_is_directory=True)
    local_env = dict(env)
    if rt:
        local_env["SALMON_LCFO_RT"] = "1"
        local_env["SALMON_LCFO_RT_CONTINUITY"] = "1"
    local_env.update(extra_env or {})
    with (folder / "inputfile").open("rb") as inp, (folder / "run.log").open("wb") as log:
        status = subprocess.run(command, cwd=folder, env=local_env, stdin=inp,
                                stdout=log, stderr=subprocess.STDOUT, timeout=120)
    log = (folder / "run.log").read_text()
    if reject:
        assert "LCFO RT: restart/checkpoint output is not supported yet" in log, folder
        assert "Native LCFO RT active" not in log, folder
    else:
        assert status.returncode == 0 and "end SALMON" in log, folder
        if rt:
            assert "Native LCFO RT active" in log and "LCFO HSE ACE build" in log, folder
    return folder

def rows(file):
    data = [[float(x) for x in line.split()] for line in file.read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    assert all(math.isfinite(x) for row in data for x in row), file
    return data

source = (repo / "testsuites/422_H_dcdft_hse/inputfile").read_text()
source = source.replace("nproc_rgrid_tot=4,1,1", "nproc_rgrid_tot=2,1,1")
source = source.replace("nproc_k=2", "nproc_k=1").replace("num_kgrid=1,2,1", "num_kgrid=1,1,1")
run("gs", source)
rt = source.replace("theory='dft'", "theory='tddft_response'")
rt = rt.replace("yn_dc='y'", "yn_dc='n'\n yn_conventional_from_dcdft='y'")
rt = rt.replace("nproc_rgrid=1,1,1", "nproc_rgrid=2,1,1")
rt = rt.replace(" nstate=4", " nstate=2").replace(" temperature_k=300d0\n", "")
rt += """
&tgrid
 nt=4
 dt=0.02d0
/
&emfield
 ae_shape1='impulse'
 e_impulse=0.0001d0
 epdir_re1=1d0,0d0,0d0
/
&analysis
 out_rt_energy_step=1
 nenergy=20
 de=0.01d0
/
"""
coarse = run("rt", rt, rt=True)
fine = run("rt_half_dt", rt.replace("nt=4", "nt=8").replace("dt=0.02d0", "dt=0.01d0"), rt=True)
a = rows(coarse / "H_dc_hse_rt.data")
b = rows(fine / "H_dc_hse_rt.data")
assert len(a) == 4 and len(b) == 8
assert abs(a[-1][0] - b[-1][0]) < 1e-12
current_error = max(abs(a[-1][j] - b[-1][j]) for j in range(13, 16))
assert current_error < 1e-10, current_error
energy = rows(coarse / "H_dc_hse_rt_energy.data")
energy_width = max(row[1] for row in energy[1:]) - min(row[1] for row in energy[1:])
assert energy_width < 1e-9, energy_width
run("reject_restart", rt.replace("&control", "&control\n yn_restart='y'"), rt=True, reject=True)
print(json.dumps(dict(work=str(root), steps=4, half_dt_current_error=current_error,
                      post_impulse_energy_width=energy_width, restart_guard="passed"), indent=2))

mlwf_input = rt.replace("xc='hse06'", "xc='hse06'\n hse_mlwf_maxiter=200")
mlwf = run("rt_mlwf_full", mlwf_input, rt=True, extra_env={"SALMON_LCFO_RT_MLWF": "1"})
assert "LCFO MLWF initial" in (mlwf / "run.log").read_text(), "MLWF path was not used"
assert "LCFO MLWF reuse" in (mlwf / "run.log").read_text(), "U was not reused"
mlwf_rows = rows(mlwf / "H_dc_hse_rt.data")
full_error = max(abs(x-y) for ra, rb in zip(a,mlwf_rows) for x,y in zip(ra,rb))
assert full_error < 1e-10, full_error
large = run("rt_mlwf_large_radius", mlwf_input, rt=True,
            extra_env={"SALMON_LCFO_RT_MLWF":"1", "SALMON_LCFO_RT_RADIUS":"8"})
large_rows = rows(large / "H_dc_hse_rt.data")
assert max(abs(x-y) for ra,rb in zip(mlwf_rows,large_rows) for x,y in zip(ra,rb)) < 1e-12
small = run("rt_mlwf_small_radius", mlwf_input, rt=True,
            extra_env={"SALMON_LCFO_RT_MLWF":"1", "SALMON_LCFO_RT_RADIUS":"3"})
small_rows = rows(small / "H_dc_hse_rt.data")
assert len(small_rows)==4 and "source-mask approximation" in (small / "run.log").read_text()
assert max(abs(ra[13]-rb[13]) for ra,rb in zip(small_rows,mlwf_rows)) > 1e-18
print(json.dumps(dict(mlwf_full_current_parity=full_error, large_radius_identity="passed",
                      small_radius="finite evolving response"),indent=2))

assert (mlwf/'lcfo_mlwf_initial.bin').read_bytes()==(large/'lcfo_mlwf_initial.bin').read_bytes()
assert (mlwf/'lcfo_mlwf_initial.bin').read_bytes()==(small/'lcfo_mlwf_initial.bin').read_bytes()
print("Identical initial U, centers and occupied coefficients across support cases")

def continuity(folder):
    return [[float(x) for x in line.split(":",1)[1].split()[1:]]
            for line in (folder/'run.log').read_text().splitlines()
            if line.startswith('LCFO EXX continuity build/')]
full_c=continuity(mlwf);small_c=continuity(small)
assert full_c and small_c
assert max(row[1] for row in full_c) < 1e-10, full_c
assert max(abs(row[0]) for row in small_c) < 1e-10, small_c
assert max(row[1] for row in small_c) > 1e-18, small_c
print(json.dumps(dict(full_continuity_L1_max=max(row[1] for row in full_c),
                      masked_continuity_L1_max=max(row[1] for row in small_c)),indent=2))
