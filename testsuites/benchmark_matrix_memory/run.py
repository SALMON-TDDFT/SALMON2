"""Actual saved-DC -> native impulse RT memory/time comparison, one job at a time."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--work', type=Path, required=True)
p.add_argument('--variants', type=Path, required=True)
p.add_argument('--reference-root', type=Path, required=True)
p.add_argument('--states', type=int, nargs='+', default=[64,128])
p.add_argument('--repeats', type=int, default=1)
p.add_argument('--modes', nargs='+', choices=['replicated','ace','all'], default=['replicated','ace','all'])
p.add_argument('--mpi', type=int, default=8)
p.add_argument('--omp', type=int, default=1)
a = p.parse_args()
root = Path(__file__).resolve().parents[2]
work = a.work.resolve(); work.mkdir(parents=True, exist_ok=True)
variants = json.loads((a.variants.resolve() / 'manifest.json').read_text())
parser_dir = root / 'testsuites/benchmark_h2_pbeh_rt'
sys.path.insert(0, str(parser_dir))
spec = importlib.util.spec_from_file_location('rt_parser', parser_dir / 'run_rt.py')
rt = importlib.util.module_from_spec(spec); spec.loader.exec_module(rt)
probe = root / 'testsuites/benchmark_h2_pbeh/rank_probe.py'
assert (a.mpi,a.omp) in [(8,1),(4,2)]
layout = '1,4,2' if a.mpi == 8 else '1,2,2'
data = dict(complete=False, variants=variants, runs=[], mpi=a.mpi, omp=a.omp,
            blas_threads=1, steps=16, dt=0.02, workload='saved DC -> Gamma PBEh(40) native RT; ACE/support controls in inputfile',
            memory_note='Rank child lifetime peak RSS includes LCFO reconstruction/initialization. Sum of RSS is not unique physical memory.',
            sampling_seconds=0.5)
result_file = work / 'results.json'
if result_file.exists(): raise RuntimeError('Choose a fresh result directory')
def save(): result_file.write_text(json.dumps(data, indent=2)+'\n')
env = dict(os.environ, OMP_NUM_THREADS=str(a.omp), OMP_DYNAMIC='FALSE', OPENBLAS_NUM_THREADS='1',
           MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
references = {}
save()
for states in a.states:
    template = a.reference_root.resolve() / f'adaptive-{states//8}x1x1-mpi{states//8}-r1'
    inp = (template / 'inputfile').read_text()
    inp = re.sub(r'nproc_rgrid=\s*\d+,\d+,\d+', 'nproc_rgrid='+layout, inp)
    inp = re.sub(r'nproc_rgrid_tot=\s*\d+,\d+,\d+', 'nproc_rgrid_tot='+layout, inp)
    seed = (template / 'data_dcdft').resolve()
    assert seed.is_dir()
    for repeat in range(1, a.repeats+1):
        # Reverse order on the second repeat to reduce warm-cache/order bias.
        modes = a.modes if repeat % 2 else list(reversed(a.modes))
        for mode in modes:
            folder = work / f'n{states}-{mode}-mpi{a.mpi}-omp{a.omp}-r{repeat}'
            folder.mkdir(); (folder / 'inputfile').write_text(inp)
            shutil.copy2(root / 'testsuites/pseudo/H_rps.dat', folder / 'H_rps.dat')
            (folder / 'data_dcdft').symlink_to(seed, target_is_directory=True)
            exe = variants['variants'][mode]['executable']
            command = ['mpiexec','--bind-to','none','-n',str(a.mpi),sys.executable,str(probe),exe]
            started = time.monotonic(); samples=[]
            row = dict(states=states, mode=mode, repeat=repeat, directory=str(folder),
                       input_sha256=hashlib.sha256(inp.encode()).hexdigest(), seed=str(seed), command=command)
            print('START',folder.name,flush=True)
            with (folder/'inputfile').open('rb') as fi, (folder/'output').open('wb') as fo:
                process = subprocess.Popen(command,cwd=folder,env=env,stdin=fi,stdout=fo,stderr=subprocess.STDOUT)
                (work/'active.json').write_text(json.dumps(dict(row, launcher_pid=process.pid),indent=2)+'\n')
                while process.poll() is None:
                    sample = subprocess.check_output(['ps','-axo','pid=,ppid=,rss=,command='],text=True)
                    rss = {}
                    for line in sample.splitlines():
                        parts = line.strip().split(None,3)
                        if len(parts)==4 and parts[3]==exe: rss[parts[0]]=int(parts[2])*1024
                    if rss: samples.append(dict(seconds=time.monotonic()-started,rss=rss))
                    time.sleep(.5)
            row['wall_seconds'] = time.monotonic()-started
            row['returncode'] = process.returncode
            (folder/'rss-samples.json').write_text(json.dumps(samples,indent=2)+'\n')
            if process.returncode:
                data['failed'] = row;save();raise RuntimeError('failed run: '+str(folder))
            measured = rt.parse_rt(folder,a.mpi,max_energy_width=None)
            row.update(measured)
            row['sampled_peak_rank_bytes'] = max((max(s['rss'].values()) for s in samples),default=0)
            row['sampled_peak_rank_sum_bytes'] = max((sum(s['rss'].values()) for s in samples),default=0)
            if states not in references and mode=='replicated': references[states]=row
            if states in references:
                reference = references[states]
                row['max_current_difference'] = float(np.max(np.abs(np.array(row['observables'])[:,13:16]-np.array(reference['observables'])[:,13:16])))
                row['max_energy_difference_ha'] = float(np.max(np.abs(np.array(row['energies'])[:,1]-np.array(reference['energies'])[:,1])))
                row['observables_match'] = row['max_current_difference']<2e-8 and row['max_energy_difference_ha']<2e-8
            data['runs'].append(row);save()
            print('END',folder.name,'seconds/step',row['rt_seconds_per_step'],
                  'peak MiB',row['peak_rank_bytes']/2**20,'match',row.get('observables_match'),flush=True)
data['complete']=True;save()
(work/'active.json').write_text(json.dumps(dict(complete=True),indent=2)+'\n')
