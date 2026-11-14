"""Diagnostic force/energy consistency along one real-space grid interval.

No solver build or restart/wavefunction/density archive. Inputs and scalar native
reports are retained. Run only with an unchanged built-and-audited runtime record.
A no-NLCC control changes the Hamiltonian and is not a production potential.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
import numpy as np
from salmonpy.dft import DFT, Pseudopotential, SCFSettings, _read_status, _read_forces

SCALAR_FILES = ('gs.inp', 'stdout.log', 'gs_status.dat', 'gs_forces.dat')


def ordered_samples(rows, intervals):
    indices = [r['index'] for r in rows]
    if (any(type(i) is not int or not 0 <= i <= intervals for i in indices)
            or len(set(indices)) != len(indices)):
        raise ValueError('invalid or duplicate sample indices')
    return sorted(rows, key=lambda r: r['index'])


def write_json(path, value):
    tmp = path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(value, indent=2)+'\n')
    tmp.replace(path)


def run_scalar_only(exe, text, pseudo_bytes, directory, env, launcher=(), pseudo_name="pseudo1.upf"):
    """Run in disposable scratch and export ONLY the requested scalar reports.

    SALMON unconditionally writes eigen/k/info side outputs. They never enter
    the results directory and are destroyed with scratch, including on failure.
    """
    directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='salmonpy-force-scalar-') as tmp:
        work = Path(tmp)
        (work/'gs.inp').write_text(text)
        if Path(pseudo_name).name != pseudo_name: raise ValueError('potential filename must be a basename')
        (work/pseudo_name).write_bytes(pseudo_bytes)
        try:
            with (work/'stdout.log').open('w') as log:
                if launcher:
                    from salmonpy.fdtd import _run_fdtd_mpi
                    returncode = _run_fdtd_mpi([*launcher,str(exe)],text,cwd=work,log=log,env=env,timeout=600)
                else:
                    proc = subprocess.run([str(exe)], cwd=work, input=text, text=True,
                        stdout=log, stderr=subprocess.STDOUT, env=env, timeout=600)
                    returncode = proc.returncode
            return returncode
        finally:
            for name in SCALAR_FILES:
                if (work/name).is_file(): shutil.copy2(work/name, directory/name)


def simpson(values, length):
    y=np.asarray(values,dtype=float)
    n=len(y)-1
    if n<2 or n%2 or not np.isfinite(y).all(): raise ValueError('finite samples at an even number of intervals required')
    return float(length/(3*n)*(y[0]+y[-1]+4*y[1:-1:2].sum()+2*y[2:-1:2].sum()))


def summarize(rows, length, natoms=8):
    e=np.array([r['energy_hartree'] for r in rows])
    f=np.array([r['path_force_hartree_per_bohr'] for r in rows])
    de=float(e[-1]-e[0]);integrals={str(stride):simpson(f[::stride],length) for stride in (1,2,4)}
    residual=de+integrals['1']
    return dict(path_length_bohr=length,delta_energy_hartree=de,integrated_force_hartree=integrals,
                consistency_residual_hartree=residual,
                consistency_residual_mHa_per_atom=1000*residual/natoms,
                quadrature_change_hartree=abs(integrals['1']-integrals['2']),
                mean_force_consistency_defect_hartree_per_bohr=-residual/length,
                # This quantity includes all discretization/force implementation errors.
                force_residual_attribution='not uniquely NLCC without controls and convergence')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build_record',type=Path);parser.add_argument('pseudo',type=Path)
    parser.add_argument('output',type=Path);parser.add_argument('--grids',type=int,nargs='+',default=[16,24])
    parser.add_argument('--intervals',type=int,default=16)
    parser.add_argument('--mode',choices=['nlcc','no-nlcc-control'],default='nlcc')
    parser.add_argument('--xc',choices=['pz','libxc_pz','libxc_pbe'],default='libxc_pbe',
                        help='non-PBE with this PBE UPF is a diagnostic implementation comparison only')
    parser.add_argument('--rigid',action='store_true')
    parser.add_argument('--point-indices',type=int,nargs='+',help='selected indices for SCF-tolerance spot checks')
    parser.add_argument('--threads',type=int,default=2)
    parser.add_argument('--tolerance',type=float,default=1e-11)
    parser.add_argument('--pseudo-sha256',required=True)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if args.intervals<8 or args.intervals%8: parser.error('intervals must be a positive multiple of eight')
    if args.point_indices is not None and any(i<0 or i>args.intervals for i in args.point_indices):
        parser.error('point index outside path')
    build=json.loads(args.build_record.read_text())
    if build.get('status')!='built-and-audited': raise ValueError('audited build record required')
    exe=Path(build['executable'])
    if hashlib.sha256(exe.read_bytes()).hexdigest()!=build['sha256']: raise ValueError('recorded binary hash mismatch')
    raw=args.pseudo.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=args.pseudo_sha256: raise ValueError('pseudopotential hash mismatch')
    header=ET.fromstring(raw).find('PP_HEADER')
    if header.get('functional','').strip().upper()!='PBE': raise ValueError('PBE UPF required')
    if ET.fromstring(raw).find('PP_NLCC') is None: raise ValueError('NLCC data required for target potential')
    if args.mode=='no-nlcc-control':
        text=raw.decode()
        text=re.sub(r'<PP_NLCC\b[^>]*>.*?</PP_NLCC>', '', text, flags=re.S)
        text=re.sub(r'core_correction\s*=\s*"[^"]*"', 'core_correction="F"', text)
        raw=text.encode()
    pseudo=Pseudopotential(14,2,'.upf',raw)
    meta=dict(executable=str(exe),executable_sha256=build['sha256'],source_build_record=str(args.build_record),
              original_pseudo_sha256=args.pseudo_sha256,effective_pseudo_sha256=pseudo.sha256,
              xc=args.xc,pseudo_generating_xc='PBE',physical_xc_pairing=args.xc=='libxc_pbe',mode=args.mode,rigid=args.rigid,threads=args.threads,
              grid_list=args.grids,intervals=args.intervals,tolerance=args.tolerance,
              lattice_bohr=10.26,states=32,kmesh=[2,2,2],restart_output=False,
              point_indices=args.point_indices,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if args.resume:
        old=json.loads((args.output/'provenance.json').read_text())
        for key in ('executable_sha256','original_pseudo_sha256','effective_pseudo_sha256','xc',
                    'mode','rigid','threads','grid_list','intervals','tolerance',
                    'lattice_bohr','states','kmesh'):
            if old.get(key)!=meta[key]: raise ValueError('resume settings differ: '+key)
        record=args.output/f'resume-{time.time_ns()}.json'
    else:
        args.output.mkdir(parents=True,exist_ok=False)
        record=args.output/'provenance.json'
        (args.output/'pseudo1.upf').write_bytes(raw)
    write_json(record,meta)
    # Conventional diamond cell, displaced first atom to avoid a symmetric path.
    base=np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0],
                   [.25,.25,.25],[.25,.75,.75],[.75,.25,.75],[.75,.75,.25]],dtype=float)
    base[0,0]+=.12/10.26
    summary={}
    for grid in args.grids:
        length=10.26/grid;rows=[]
        samples=args.output/f'grid{grid}'/'samples.json'
        if args.resume and samples.exists(): rows=ordered_samples(json.loads(samples.read_text()),args.intervals)
        accepted={r['index'] for r in rows}
        for i,shift in enumerate(np.linspace(0,length,args.intervals+1)):
            if args.point_indices is not None and i not in args.point_indices: continue
            if i in accepted: continue
            positions=base.copy()
            if args.rigid: positions[:,0]+=shift/10.26
            else: positions[0,0]+=shift/10.26
            model=DFT((14,)*8,np.eye(3)*10.26,positions,(pseudo,),32,32,(grid,)*3,(2,)*3,args.xc)
            text=model.input_text(SCFSettings(max_iterations=600,tolerance=args.tolerance))
            text=text.replace("write_gs_restart_data='all'", "write_gs_restart_data='no'")
            directory=args.output/f'grid{grid}'/f'point{i:02d}'
            if directory.exists():
                backup=directory/f'interrupted-{time.time_ns()}'
                backup.mkdir()
                for name in SCALAR_FILES:
                    if (directory/name).exists(): shutil.move(directory/name,backup/name)
            env={k:v for k,v in os.environ.items() if k!='SALMONPY_TDDFT_WORKER'}
            env.update(OMP_NUM_THREADS=str(args.threads),OPENBLAS_NUM_THREADS='1')
            start=time.monotonic()
            returncode=run_scalar_only(exe,text,raw,directory,env)
            if returncode or 'end SALMON' not in (directory/'stdout.log').read_text():
                raise RuntimeError(f'failed SCF: {directory}')
            status=_read_status(directory/'gs_status.dat');forces=_read_forces(directory/'gs_forces.dat',8)
            if any(p.name in ('wfn.bin','rho.bin','rho_inout.bin','occupation.bin','info.bin') for p in directory.rglob('*')):
                raise RuntimeError('unexpected GS restart data output')
            row=dict(index=i,shift_bohr=float(shift),**status,seconds=time.monotonic()-start,
                     forces_hartree_per_bohr=forces.tolist(),
                     path_force_hartree_per_bohr=float(forces[:,0].sum() if args.rigid else forces[0,0]))
            rows.append(row)
            rows=ordered_samples(rows,args.intervals)
            write_json(directory.parent/'samples.json',rows)
            print(json.dumps(dict(grid=grid,index=i,energy=row['energy_hartree'],force=row['path_force_hartree_per_bohr'],
                                  iterations=row['iterations'],seconds=row['seconds'])),flush=True)
        summary[str(grid)]=summarize(rows,length) if len(rows)==args.intervals+1 else {'partial_path':True,'points':len(rows)}
        write_json(args.output/'summary.json',summary)
        print(json.dumps(dict(grid=grid,summary=summary[str(grid)])),flush=True)

if __name__=='__main__': main()
