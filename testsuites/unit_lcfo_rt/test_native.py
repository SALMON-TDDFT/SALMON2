#!/usr/bin/env python3
"""Sequential fragment x orbital MPI native integration test. Needs a USE_HSE binary and H_rps.dat.
Usage: python3 test_native.py --binary /abs/salmon --pseudo /abs/H_rps.dat
Each run owns a fresh temporary directory; no production data is modified.
"""
import argparse
import json
import math
import re
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
parser.add_argument("--omp", type=int, default=1)
args = parser.parse_args()
assert args.omp > 0
root = Path(tempfile.mkdtemp(prefix="lcfo-native-test-", dir=args.work_root))
repo = Path(__file__).resolve().parents[2]
env = dict(os.environ, OMP_NUM_THREADS=str(args.omp), OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
env.pop("SALMON_LCFO_RT", None)
command = [args.mpirun, "-np", "2", str(args.binary.resolve())]

def run(name, text, rt=False, reject=False, extra_env=None, orbital_groups=1):
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
    local_command=command.copy();local_command[2]=str(2*orbital_groups)
    with (folder / "inputfile").open("rb") as inp, (folder / "run.log").open("wb") as log:
        status = subprocess.run(local_command, cwd=folder, env=local_env, stdin=inp,
                                stdout=log, stderr=subprocess.STDOUT, timeout=120)
    log = (folder / "run.log").read_text()
    if reject:
        assert (reject if isinstance(reject,str) else "LCFO RT: restart/checkpoint output is not supported yet") in log, folder
        if reject is True:
            assert "Native LCFO RT active" not in log, folder
    else:
        assert status.returncode == 0 and "end SALMON" in log, folder
        if rt:
            assert "Native LCFO RT active" in log and "LCFO HSE ACE build" in log, folder
            storage=re.findall(r'LCFO distributed storage rank/local/global/halo rows:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',log)
            assert len(storage)==2, (folder,storage)
            assert all(0<int(local)<int(total) and 0<int(halo)<=int(total) for rank,local,total,halo in storage)
            assert 'LCFO timing pack/source/exchange/ACE build' in log, folder
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
 yn_out_dns_rt='y'
 out_dns_rt_step=4
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
# Half box diagonal is sqrt(8**2+4**2+4**2) = 9.798 bohr.
large = run("rt_mlwf_large_radius", mlwf_input, rt=True,
            extra_env={"SALMON_LCFO_RT_MLWF":"1", "SALMON_LCFO_RT_RADIUS":"10"})
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
for width in (4,):
    batched=run(f'rt_fft_batch{width}',mlwf_input,rt=True,
                extra_env={'SALMON_LCFO_RT_MLWF':'1','SALMON_LCFO_RT_RADIUS':'3',
                           'SALMON_LCFO_RT_FFT_BATCH':str(width)})
    batch_rows=rows(batched/'H_dc_hse_rt.data')
    batch_error=max(abs(x-y) for a,b in zip(batch_rows,small_rows) for x,y in zip(a,b))
    assert batch_error<1e-12,batch_error
    stats=re.findall(r'LCFO exchange FFT batches rank/width/calls:\s+(\d+)\s+(\d+)\s+(\d+)',
                     (batched/'run.log').read_text())
    assert stats and all(int(x[1])==width for x in stats),stats
    print('Native batch parity',width,batch_error)
run('reject_fft_batch',mlwf_input,rt=True,reject='LCFO HSE: FFT batch must be an integer from 1 to 32',
    extra_env={'SALMON_LCFO_RT_FFT_BATCH':'0'})


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

# A physical-step interval must reduce expensive exchange builds, not merely
# relabel a cache hit. Interval1 must reproduce the default trajectory.
def builds(folder):
    return (folder/'run.log').read_text().count('LCFO HSE ACE build')
one=run('rt_ace_interval1',mlwf_input,rt=True,
        extra_env={'SALMON_LCFO_RT_MLWF':'1','SALMON_LCFO_RT_ACE_INTERVAL':'1'})
assert rows(one/'H_dc_hse_rt.data')==mlwf_rows
assert builds(one)==builds(mlwf)
for interval in (2,4):
    reused=run(f'rt_ace_interval{interval}',mlwf_input,rt=True,
        extra_env={'SALMON_LCFO_RT_MLWF':'1','SALMON_LCFO_RT_ACE_INTERVAL':str(interval)})
    log=(reused/'run.log').read_text()
    assert 'LCFO HSE ACE retained at step' in log, 'Missing physical-step ACE reuse'
    assert 'LCFO HSE impulse ACE rebuilt before first predictor' in log
    steps=[int(line.split()[-1]) for line in log.splitlines() if line.startswith('LCFO HSE exchange rebuilt at step')]
    assert 1 in steps and all(step==0 or (step-1)%interval==0 for step in steps), steps
    assert builds(reused)<builds(mlwf), (builds(reused),builds(mlwf))
    assert len(rows(reused/'H_dc_hse_rt.data'))==4
    assert len(rows(reused/'H_dc_hse_rt_energy.data'))==5
    assert (reused/'lcfo_mlwf_initial.bin').read_bytes()==(mlwf/'lcfo_mlwf_initial.bin').read_bytes()
    print(json.dumps({'ACE_interval':interval,'builds':builds(reused),'default_builds':builds(mlwf)}))

run('reject_ace_interval',mlwf_input,rt=True,reject='LCFO HSE: ACE interval must be a positive integer',
    extra_env={'SALMON_LCFO_RT_ACE_INTERVAL':'0'})

laser_input=mlwf_input.replace("ae_shape1='impulse'","ae_shape1='Acos2'\n I_wcm2_1=1d8\n omega1=0.5d0\n tw1=100d0").replace('nenergy=20','nenergy=0')
laser=run('rt_ace_laser',laser_input,rt=True,
          extra_env={'SALMON_LCFO_RT_MLWF':'1','SALMON_LCFO_RT_ACE_INTERVAL':'4'})
laser_log=(laser/'run.log').read_text()
assert 'LCFO HSE ACE retained at step        1' in laser_log
laser_steps=[int(line.split()[-1]) for line in laser_log.splitlines() if line.startswith('LCFO HSE exchange rebuilt at step')]
assert 1 not in laser_steps and 4 in laser_steps,laser_steps
print('Impulse first step refresh and smooth laser initial ACE reuse passed')

# Same GS and physics; compare complete trajectories, not truncated zip pairs.
def compare_tables(a, b, columns, tolerance):
    assert len(a) == len(b) and a, (len(a), len(b))
    assert all(len(x) == len(y) for x,y in zip(a,b))
    assert all(abs(x[0]-y[0]) < 1e-12 for x,y in zip(a,b))
    error=max(abs(x[j]-y[j]) for x,y in zip(a,b) for j in columns)
    assert error < tolerance, error
    return error

def density(folder):
    lines=(folder/'H_dc_hse_dns_000004.cube').read_text().splitlines()
    natom=abs(int(lines[2].split()[0]))
    data=[float(x) for line in lines[6+natom:] for x in line.split()]
    assert data and all(math.isfinite(x) for x in data)
    return data

def compare_split(name, baseline, text, extra_env):
    split=run(name,text.replace('nproc_ob=1','nproc_ob=2'),rt=True,
              orbital_groups=2,extra_env=extra_env)
    error=compare_tables(rows(baseline/'H_dc_hse_rt.data'),
                         rows(split/'H_dc_hse_rt.data'),range(13,16),1e-11)
    energy_error=compare_tables(rows(baseline/'H_dc_hse_rt_energy.data'),
                                rows(split/'H_dc_hse_rt_energy.data'),[1],1e-9)
    da,db=density(baseline),density(split)
    assert len(da)==len(db)
    density_error=max(abs(x-y) for x,y in zip(da,db))
    assert density_error < 1e-11,density_error
    assert builds(split)==builds(baseline)
    def schedule(folder):
        return [line for line in (folder/'run.log').read_text().splitlines()
                if line.startswith('LCFO HSE exchange rebuilt at step')]
    assert schedule(split)==schedule(baseline)
    print(json.dumps(dict(case=name,max_current_difference=error,
                          max_energy_difference=energy_error,
                          max_density_difference=density_error)))
    return split

for interval,baseline in [(1,mlwf),(4,root/'rt_ace_interval4')]:
    compare_split(f'rt_orbital2_ace{interval}',baseline,mlwf_input,
                  {'SALMON_LCFO_RT_MLWF':'1','SALMON_LCFO_RT_ACE_INTERVAL':str(interval)})
compare_split('rt_orbital2_small_radius',small,
              mlwf_input.replace('&parallel',"&parallel\n process_allocation='orbital_sequential'"),
              {'SALMON_LCFO_RT_MLWF':'1','SALMON_LCFO_RT_RADIUS':'3'})
compare_split('rt_orbital2_laser',laser,laser_input,
              {'SALMON_LCFO_RT_MLWF':'1','SALMON_LCFO_RT_ACE_INTERVAL':'4'})
# Three occupied target states split 2+1. Reoccupy the same saved basis in
# both runs; this is an algebraic distribution test, not a GS accuracy test.
unequal_input=rt.replace(' nstate=2',' nstate=3').replace(' nelec=4',' nelec=6')
unequal=run('rt_unequal_reference',unequal_input,rt=True)
compare_split('rt_orbital2_unequal',unequal,unequal_input,{})
