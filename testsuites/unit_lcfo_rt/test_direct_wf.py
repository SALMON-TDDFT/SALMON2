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

pt_input=mlwf_input
pt_env={'SALMON_LCFO_RT_MLWF':'1','SALMON_LCFO_RT_DIRECT_WF':'1'}
pt=run('direct_wf',pt_input,rt=True,extra_env=pt_env)
log=(pt/'run.log').read_text()
assert 'LCFO direct WF coefficient Taylor4' in log
assert log.count('LCFO direct WF accepted Gram error:')==5, 'initial and four accepted frames required'
assert log.count('LCFO HSE ACE build')==(mlwf/'run.log').read_text().count('LCFO HSE ACE build'), 'rebase rebuilt unchanged ACE'
ptrows=rows(pt/'H_dc_hse_rt.data');assert len(ptrows)==4
error=max(abs(x[j]-y[j]) for x,y in zip(ptrows,mlwf_rows) for j in (13,14,15))
assert error<1e-11,error
split=run('direct_wf_orbital2',pt_input.replace('nproc_ob=1','nproc_ob=2'),rt=True,extra_env=pt_env,orbital_groups=2)
splitrows=rows(split/'H_dc_hse_rt.data');assert len(splitrows)==4
spliterror=max(abs(x[j]-y[j]) for x,y in zip(ptrows,splitrows) for j in (13,14,15))
assert spliterror<1e-11,spliterror
print(json.dumps({'direct_current_vs_taylor':error,'direct_orbital_split_difference':spliterror,'work':str(root)}))

def compare_density_energy(a,b,step):
    def cube(folder):
        lines=(folder/f'H_dc_hse_dns_{step:06d}.cube').read_text().splitlines()
        nat=abs(int(lines[2].split()[0]));return [float(x) for line in lines[6+nat:] for x in line.split()]
    da,db=cube(a),cube(b);assert len(da)==len(db)
    density_error=max(abs(x-y) for x,y in zip(da,db));assert density_error<1e-10,density_error
    ea,eb=rows(a/'H_dc_hse_rt_energy.data'),rows(b/'H_dc_hse_rt_energy.data')
    assert len(ea)==len(eb)
    energy_error=max(abs(x[1]-y[1]) for x,y in zip(ea,eb));assert energy_error<1e-9,energy_error
    print(json.dumps({'direct_case':a.name,'max_density_difference':density_error,'max_energy_difference':energy_error}))
compare_density_energy(pt,mlwf,4)
compare_density_energy(split,pt,4)
for label,extra in [('r3',{'SALMON_LCFO_RT_RADIUS':'3'}),
                    ('r3_ace4_u2',{'SALMON_LCFO_RT_RADIUS':'3','SALMON_LCFO_RT_ACE_INTERVAL':'4','SALMON_LCFO_RT_U_INTERVAL':'2'})]:
    baseline_env={'SALMON_LCFO_RT_MLWF':'1',**extra}
    base=run('baseline_'+label,mlwf_input,rt=True,extra_env=baseline_env)
    direct=run('direct_'+label,pt_input,rt=True,extra_env={**baseline_env,'SALMON_LCFO_RT_DIRECT_WF':'1'})
    aa,bb=rows(base/'H_dc_hse_rt.data'),rows(direct/'H_dc_hse_rt.data');assert len(aa)==len(bb)==4
    err=max(abs(x[j]-y[j]) for x,y in zip(aa,bb) for j in (13,14,15));assert err<1e-11,err
    compare_density_energy(direct,base,4)
    assert (direct/'run.log').read_text().count('LCFO HSE ACE build')==(base/'run.log').read_text().count('LCFO HSE ACE build')
    print(json.dumps({'direct_case':label,'current_difference':err}))
ptfine=run('direct_half_dt',pt_input.replace('nt=4','nt=8').replace('dt=0.02d0','dt=0.01d0'),rt=True,extra_env=pt_env)
fr=rows(ptfine/'H_dc_hse_rt.data');assert len(fr)==8
assert max(abs(fr[-1][j]-ptrows[-1][j]) for j in (13,14,15))<1e-10
run('reject_direct_without_mlwf',pt_input,rt=True,reject='LCFO MLWF: invalid radius or U interval',
    extra_env={'SALMON_LCFO_RT_DIRECT_WF':'1','SALMON_LCFO_RT_MLWF':'0'})
run('reject_direct_flag',pt_input,rt=True,reject='LCFO direct WF: flag must be 0 or 1',
    extra_env={'SALMON_LCFO_RT_DIRECT_WF':'bad'})
for groups in (1,2):
    measured=run('direct_measured_'+str(groups),pt_input.replace('nproc_ob=1',f'nproc_ob={groups}'),
                 rt=True,extra_env={**pt_env,'SALMON_LCFO_RT_FFT_MEASURE':'1'},orbital_groups=groups)
    measured_rows=rows(measured/'H_dc_hse_rt.data')
    err=max(abs(x[j]-y[j]) for x,y in zip(ptrows,measured_rows) for j in (13,14,15))
    assert err<1e-11,err
    assert 'LCFO HSE FFT measured planning:  1' in (measured/'run.log').read_text()
    compare_density_energy(measured,pt,4)
for flag in ('2','bad','1 trailing','12345678901234567'):
    run('reject_measure_'+flag.replace(' ','_'),pt_input,rt=True,
        reject='LCFO HSE: FFT measure must be 0 or 1',
        extra_env={**pt_env,'SALMON_LCFO_RT_FFT_MEASURE':flag})
print('Direct coefficient Taylor4 regression passed')
