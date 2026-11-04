"""Run only in a fresh scratch directory. The validated SHA comes from the audited build."""
from pathlib import Path
import argparse, hashlib, json, os, re, subprocess, sys, time
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("salmon",type=Path)
    ap.add_argument("sha256")
    ap.add_argument("scratch",type=Path)
    ap.add_argument("--ranks",type=int,choices=(1,2),default=1)
    args=ap.parse_args()
    exe=args.salmon.resolve()
    if hashlib.sha256(exe.read_bytes()).hexdigest()!=args.sha256:
        raise ValueError("validated executable SHA mismatch")
    suite=Path(__file__).resolve().parent
    root=args.scratch.resolve()
    root.mkdir(exist_ok=False,parents=True)
    result={"mpi":args.ranks,"omp":1,"sha256":args.sha256,"runs":{}}
    def run(name,inp):
        d=root/name
        d.mkdir()
        (d/'Si_rps.dat').write_bytes((suite.parent/'pseudo/Si_rps.dat').read_bytes())
        inp=inp.replace('nproc_k=1',f'nproc_k={args.ranks}')
        (d/'input').write_text(inp)
        cmd=[str(exe)] if args.ranks==1 else ['mpiexec','-n',str(args.ranks),str(exe)]
        start=time.monotonic()
        with (d/'run.log').open('w') as log:
            subprocess.run(cmd,cwd=d,input=inp,text=True,stdout=log,stderr=subprocess.STDOUT,
                           env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),check=True)
        result['runs'][name]={'seconds':time.monotonic()-start}
        log=(d/'run.log').read_text()
        if 'end SALMON' not in log:
            raise ValueError(name+': incomplete SALMON execution')
        return d,log
    d,log=run('off',(suite/'si-off.inp').read_text())
    # Upstream regression 131 also compares an empty conduction level. This
    # occupied-only smoke checks the same total energy and valence splitting.
    eoff=float(re.search(r'Total energy \(eV\)\s*=\s*(\S+)',(d/'Si_info.data').read_text())[1])
    bands=np.loadtxt([line for line in (d/'Si_eigen.data').read_text().splitlines()
                      if line.strip() and line.strip()[0].isdigit()])
    if abs(eoff+214.113848465987)>0.001 or abs(bands[1,1]-bands[0,1]-.271266)>4e-5:
        raise ValueError("occupied-only LSDA reference")
    if bands.shape!=(32,3) or not np.all(bands[:,2]==2): raise ValueError("GS occupied-band count")
    d,log=run('warm',(suite/'si.inp').read_text())
    m=re.search(r'#GS converged at\s+(\d+)\s+:\s+(\S+)',log)
    if not m or not (d/'data_for_restart/kli-v1.bin').is_file():
        raise ValueError('GS did not converge with a KLI checkpoint')
    result['gs']={'iteration':int(m[1]),'scf_residual':float(m[2]),
                  'energy_eV':float(re.search(r'Total energy \(eV\)\s*=\s*(\S+)',(d/'Si_info.data').read_text())[1]),
                  'audit':[float(x) for x in re.findall(r'KLI orth/rho/Slater/row/compat/solve/box: (.*)',log)[-1].split()],
                  'localization':[float(x) for x in re.findall(r'KLI spin/status/iter/FFT E/loc/boundary: (.*)',log)[-1].split()]}
    for mode in ['frozen','instantaneous']:
        inp=(suite/('si-'+mode+'.inp')).read_text()
        currents=[]
        for case,nt,dt,intensity in [('zero',32,.04,'0d0'),('pulse',32,.04,'1d8'),('half',64,.02,'1d8')]:
            name=mode+'-'+case
            d,log=run(name,inp.replace('nt = 32',f'nt = {nt}').replace('dt = 0.04',f'dt = {dt}').replace('I_wcm2_1 = 0d0','I_wcm2_1 = '+intensity))
            if 'real-space occupied orbitals per k/spin: 4' not in log: raise ValueError('RT propagated state count')
            rt=np.loadtxt(d/'Si_rt.data');currents.append(rt)
            diag=np.array([[float(x) for x in line.split()] for line in re.findall(r'KLI RT density/orthogonality/potential: (.*)',log)])
            if diag.shape!=(nt+1,3) or len(rt)!=nt: raise ValueError(name+': missing endpoint')
            stats=result['runs'][name]
            stats.update(max_density=float(max(diag[:,0])),max_orth=float(max(diag[:,1])),max_v=float(max(diag[:,2])),max_current=float(np.max(abs(rt[:,13:16]))))
            if mode=='instantaneous':
                loc=np.array([[float(x) for x in line.split()] for line in re.findall(r'KLI spin/status/iter/FFT E/loc/boundary: (.*)',log)])
                if len(loc)!=2*nt: raise ValueError('missing predictor or endpoint')
                stats.update(max_iter=int(max(loc[:,2])),ffts=int(sum(loc[:,3])))
                if case=='zero' and stats['max_iter']!=0: raise ValueError('unexpected localization iterations')
            if not np.all(np.isfinite(diag)) or not np.all(np.isfinite(rt)) or stats['max_orth']>1e-8: raise ValueError(name)
            if case=='zero' and max(stats['max_density'],stats['max_current'],stats['max_v'])>1e-8: raise ValueError(name+': drift')
        coarse,fine=currents[1],currents[2][1::2]
        np.testing.assert_allclose(coarse[:,0],fine[:,0],atol=1e-12,rtol=0)
        delta=float(np.max(abs(coarse[:,13:16]-fine[:,13:16])))
        signal=float(np.max(abs(fine[:,13:16])))
        if signal<=1e-8 or delta>0.01*signal: raise ValueError('unresolved pulse response')
        result[mode+'_current_dt_error']=delta
    # Drive the occupied subspace far enough to exercise both paths with the
    # same GS checkpoint and production tolerance used by the zero-field case.
    name = 'instantaneous-relocalize'
    inp = (suite/'si-instantaneous.inp').read_text().replace('I_wcm2_1 = 0d0', 'I_wcm2_1 = 1d12')
    directory, log = run(name, inp)
    rows = [line.split()[3:] for line in log.splitlines() if line.startswith('KLI RT localization ')]
    stages = ['predictor', 'corrector'] * 32
    if [row[0] for row in rows] != stages:
        raise ValueError('missing or misordered dynamic KLI evaluations')
    for stage in ('predictor', 'corrector'):
        values = np.array([[float(x) for x in row[1:]] for row in rows if row[0] == stage])
        iterations, initial, final, sigma = values.T
        rotated = iterations > 0
        if not np.all(np.isfinite(values)) or not np.all(sigma > 1e-8):
            raise ValueError('invalid transported gauge')
        if not (np.any(rotated) and np.any(~rotated)):
            raise ValueError(stage + ': expected both reuse and additional localization')
        if not np.array_equal(rotated, initial > 1e-5) or np.any(final > 1e-5):
            raise ValueError(stage + ': PZ residual/iteration mismatch')
        result['runs'][name][stage] = dict(rotated=int(sum(rotated)), reused=int(sum(~rotated)),
            max_iter=int(max(iterations)), max_initial=float(max(initial)), max_final=float(max(final)))
    endpoint = np.array([[float(x) for x in line.split()] for line in
                         re.findall(r'KLI RT density/orthogonality/potential: (.*)', log)])
    rt = np.loadtxt(directory/'Si_rt.data')
    if (endpoint.shape != (33, 3) or rt.shape[0] != 32 or not np.isfinite(endpoint).all()
            or not np.isfinite(rt).all() or endpoint[:, 1].max() > 1e-8 or endpoint[:, 2].max() <= 1e-6):
        raise ValueError('invalid driven real-space trajectory')
    result['runs'][name].update(max_orth=float(endpoint[:, 1].max()), max_v=float(endpoint[:, 2].max()))
    (root/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
