"""Isolated-process comparison to a preserved literal legacy Python SBE adapter.

PYTHONPATH=. python tools/benchmark_sbe_handles.py LIBRARY LEGACY_PY NEW_OUTPUT
No build/GS data. Peak RSS sampled before dense diagnostic exports.
"""
import argparse,hashlib,importlib.util,json,os,subprocess,sys,time
from pathlib import Path
import numpy as np


def worker(args):
    from salmonpy._material_handles import SBEAPI as API
    import salmonpy.sbe as current
    api=API(args.library);current._handles_api=lambda:api;current._native=lambda:api.legacy
    if args.mode=='legacy':
        name='salmonpy._sbe_legacy_validation';spec=importlib.util.spec_from_file_location(name,args.legacy)
        module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
        module._native=lambda:api.legacy
    else:module=current
    rng=np.random.default_rng(38219);nb,nk=args.bands,args.kpoints
    p=rng.normal(size=(nb,nb,3,nk))+1j*rng.normal(size=(nb,nb,3,nk))
    p=(p+p.swapaxes(0,1).conj())*.01/np.sqrt(nb)
    eigen=np.linspace(-.5,.7,nb)[:,None]+np.arange(nk)[None,:]*.001
    weights=rng.uniform(.1,1.,nk);values=rng.normal(size=(args.steps+6,args.states,3))*.001
    start=time.perf_counter();model=module.SBEInitialState(eigen,p,p*.13,weights,np.diag([4.,5.,6.]),2)
    states=[model.new_state(dt=.002,initial_A=values[0,i]) for i in range(args.states)]
    setup=time.perf_counter()-start
    for n in range(1,6):
        for i,state in enumerate(states):state.step(values[n,i])
    history=[];start=time.perf_counter()
    for n in range(6,args.steps+6):
        history.append([state.step(values[n,i]).copy() for i,state in enumerate(states)])
    elapsed=time.perf_counter()-start
    try:
        import resource
        raw=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        peak=int(raw if sys.platform=='darwin' else raw*1024)
    except (ImportError,OSError,ValueError):peak=None
    rho_bytes=16*nb*nb*nk;model_bytes=8*nb*nk+8*nk+96*nb*nb*nk
    record={'mode':args.mode,'bands':nb,'kpoints':nk,'states':args.states,'steps':args.steps,
        'setup_seconds':setup,'step_seconds':elapsed,'seconds_per_state_step':elapsed/(args.steps*args.states),
        'peak_RSS_bytes_before_density_export':peak,'RSS_scope':'whole worker process high-water mark; not a hard limit',
        'modeled_explicit_copy_bytes_per_state_step':model_bytes+3*rho_bytes if args.mode=='legacy' else rho_bytes,
        'copy_scope':'explicit Python rho copy + native gs/rho assignments; excludes arithmetic/workspace and vector copies',
        'matrix_model_copies_per_step':1 if args.mode=='legacy' else 0,
        'boundary_dense_arrays_per_step':5 if args.mode=='legacy' else 0,
        'boundary_note':'pointer arguments are not network/device transfer bytes',
        'library_sha256':api.sha256,'legacy_python_sha256':hashlib.sha256(args.legacy.read_bytes()).hexdigest()}
    np.savez(args.output.with_suffix('.npz'),current=np.asarray(history),density=np.stack([s.density for s in states]))
    args.output.write_text(json.dumps(record,indent=2)+'\n')
    for state in states:
        close=getattr(state,'close',None)
        if close:close()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('library',type=Path);p.add_argument('legacy',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--mode',choices=['legacy','persistent']);p.add_argument('--bands',type=int,default=16)
    p.add_argument('--kpoints',type=int,default=32);p.add_argument('--states',type=int,default=1)
    p.add_argument('--steps',type=int,default=60)
    args=p.parse_args()
    if args.mode:return worker(args)
    args.output.mkdir(parents=True,exist_ok=False)
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'}
    records=[];comparisons=[]
    for nb,nk,states in ((8,128,1),(32,8,1),(32,8,2)):
        for repeat in range(3):
            pair={}
            for mode in ('legacy','persistent'):
                path=args.output/f'nb{nb}-nk{nk}-s{states}-r{repeat}-{mode}.json'
                with path.with_suffix('.log').open('w') as log:
                    subprocess.run([sys.executable,str(Path(__file__).resolve()),str(args.library),str(args.legacy),str(path),
                        '--mode',mode,'--bands',str(nb),'--kpoints',str(nk),'--states',str(states),'--steps',str(args.steps)],
                        env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
                row=json.loads(path.read_text());records.append(row)
                with np.load(path.with_suffix('.npz'),allow_pickle=False) as data:pair[mode]={k:data[k].copy() for k in data.files}
            comparisons.append({'bands':nb,'kpoints':nk,'states':states,'repeat':repeat,
                'maximum_current_difference':float(np.max(abs(pair['legacy']['current']-pair['persistent']['current']))),
                'maximum_density_difference':float(np.max(abs(pair['legacy']['density']-pair['persistent']['density']))),
                'bitwise_equal':all(pair['legacy'][k].tobytes()==pair['persistent'][k].tobytes() for k in pair['legacy'])})
    summary=[]
    for nb,nk,states in ((8,128,1),(32,8,1),(32,8,2)):
        rows={mode:[r for r in records if (r['bands'],r['kpoints'],r['states'],r['mode'])==(nb,nk,states,mode)]
              for mode in ('legacy','persistent')}
        old=float(np.median([r['step_seconds'] for r in rows['legacy']]))
        new=float(np.median([r['step_seconds'] for r in rows['persistent']]))
        summary.append({'bands':nb,'kpoints':nk,'states':states,'legacy_median_step_seconds':old,
            'persistent_median_step_seconds':new,'speedup':old/new,
            'legacy_peak_RSS_bytes':max((r['peak_RSS_bytes_before_density_export'] for r in rows['legacy'] if r['peak_RSS_bytes_before_density_export'] is not None),default=None),
            'persistent_peak_RSS_bytes':max((r['peak_RSS_bytes_before_density_export'] for r in rows['persistent'] if r['peak_RSS_bytes_before_density_export'] is not None),default=None),
            'legacy_modeled_copy_bytes_per_state_step':rows['legacy'][0]['modeled_explicit_copy_bytes_per_state_step'],
            'persistent_modeled_copy_bytes_per_state_step':rows['persistent'][0]['modeled_explicit_copy_bytes_per_state_step']})
    result={'cases':summary,'numerical_comparisons':comparisons,'scope':'synthetic isolated-process measurements; modeled copy counts are not a hardware profiler'}
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
