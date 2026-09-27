"""Bounded fixed-cell water checks; run in a fresh work directory.

Coarse grids test consistency only, not converged liquid-water physics.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[2]

def run(exe,work,name,inp):
    path=work/name;path.mkdir(parents=True,exist_ok=False)
    (path/'inputfile').write_text(inp)
    for atom in ('H','O'):shutil.copy(ROOT/f'testsuites/pseudo/{atom}_rps.dat',path)
    with (path/'inputfile').open('rb') as f,(path/'output').open('wb') as out:
        p=subprocess.run([str(exe)],stdin=f,stdout=out,stderr=subprocess.STDOUT,cwd=path,
             env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMPI_MCA_btl='self,vader'),timeout=900)
    log=(path/'output').read_text()
    if p.returncode or 'end SALMON' not in log:raise RuntimeError(f'Failed: {path}/output')
    if '#GS does not converged' in log:raise RuntimeError(f'SCF did not converge: {name}')
    if "theory='dft_md'" not in inp and '#GS converged at' not in log:raise RuntimeError(f'Missing convergence: {name}')
    result=dict(name=name)
    if "theory='dft_md'" in inp:
        data=np.loadtxt(path/'H2O_dft_md.data')
        if data.shape[0]!=int(re.search(r'nt=(\d+)',inp)[1])+1:raise AssertionError('Missing MD steps')
        result['md']=data.tolist()
        result['max_energy_drift_eV']=float(np.max(abs(data[:,7])))
    else:
        info=(path/'H2O_info.data').read_text()
        result['energy_eV']=float(re.search(r'Total energy \(eV\)\s*=\s*(\S+)',info)[1])
        forces=info.split('Force [eV/A]')[1].strip().splitlines()[:3]
        result['forces_eV_A']=[[float(x) for x in line.split()[1:]] for line in forces]
    return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('executable',type=Path);parser.add_argument('work',type=Path)
    parser.add_argument('--atom',choices=['H','O'],default='H');parser.add_argument('--xc',default='pbeh40_rvv10');parser.add_argument('--grid',type=int,default=24);parser.add_argument('--md',action='store_true')
    args=parser.parse_args();exe=args.executable.resolve();work=args.work.resolve();work.mkdir(parents=True,exist_ok=True)
    base=(Path(__file__).with_name('water.inp')).read_text().replace('16,16,16',','.join([str(args.grid)]*3))
    base=base.replace("xc ='pbeh40_rvv10'",f"xc ='{args.xc}'")
    base=base.replace('threshold=1d-9','threshold=1d-10').replace('nscf = 300','nscf = 500')
    results=[]
    if args.md:
        for dt in (.1,.05):
            inp=base.replace("theory='dft'","theory='dft_md'")
            inp+=f"\n&tgrid\n dt={dt}\n nt={round(.4/dt)}\n/\n&md\n ensemble='NVE'\n yn_set_ini_velocity='y'\n temperature0_ion_k=300\n step_update_ps=1\n/\n"
            results.append(run(exe,work,f'md_{dt}',inp))
            (work/'results.json').write_text(json.dumps(results,indent=2)+'\n')
            print(json.dumps(results[-1]),flush=True)
        assert all(x['max_energy_drift_eV']<1e-4 for x in results), 'NVE energy drift too large'
        assert results[1]['max_energy_drift_eV']<results[0]['max_energy_drift_eV'], 'Halving dt did not reduce drift'
    else:
        delta=.001
        for name,shift in [('central',0),('plus',delta),('minus',-delta)]:
            coordinate='3.55908403' if args.atom=='H' else '2.61768087'
            inp=base.replace(coordinate,f'{float(coordinate)+shift:.8f}')
            results.append(run(exe,work,name,inp))
            print(json.dumps(results[-1]),flush=True)
        fd=-(results[1]['energy_eV']-results[2]['energy_eV'])/(2*delta)
        err=abs(fd-results[0]['forces_eV_A'][0 if args.atom=='H' else 2][0])
        report=dict(atom=args.atom,grid=args.grid,displacement_A=delta,force_fd_eV_A=fd,error_eV_A=err,runs=results)
        (work/'results.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report),flush=True)
        assert err<.001, f'Force discrepancy {err} eV/A'

if __name__=='__main__':main()
