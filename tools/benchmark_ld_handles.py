"""Compare literal prechange LD adapter with persistent handles, isolated workers.

PYTHONPATH=. python tools/benchmark_ld_handles.py LIBRARY LEGACY_PY NEW_OUTPUT
Synthetic E fields only. Copy counts are source-derived, not network traffic.
"""
import argparse,hashlib,importlib.util,json,os,subprocess,sys,time
from pathlib import Path
import numpy as np


def worker(a):
    from salmonpy._material_handles import LDAPI as API
    import salmonpy.ld as current
    api=API(a.library);current._handles_api=lambda:api;current._native=lambda:api.legacy
    if a.mode=='legacy':
        name='salmonpy._ld_legacy_validation';spec=importlib.util.spec_from_file_location(name,a.legacy)
        module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
        module._native=lambda:api.legacy
    else:module=current
    fields=np.random.default_rng(188).normal(size=(a.steps+1,3))*.002
    t=time.perf_counter()
    model=module.LorentzDrude(plasma_frequency=1.,strength=np.ones(a.poles)/a.poles,
                            resonance=np.linspace(0,2,a.poles),damping=np.full(a.poles,.1))
    state=model.new_state(dt=.01,initial_field=fields[0]);setup=time.perf_counter()-t
    # Setup is excluded; the native library was explicitly loaded for both modes.
    p=np.empty_like(fields);j=np.empty_like(fields);p[0]=state.polarization;j[0]=state.current
    t=time.perf_counter()
    if a.route=='step':
        for n,e in enumerate(fields[1:],1):p[n],j[n]=state.step(e)
    else:p,j=state._advance(fields)
    elapsed=time.perf_counter()-t
    try:
        import resource
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        peak=int(rss if sys.platform=='darwin' else rss*1024)
    except (ImportError,OSError,ValueError):peak=None
    # Snapshot export is outside both timing and peak-RSS measurement.
    data=state._snapshot_arrays();np.savez(a.output.with_suffix('.npz'),polarization=p,current=j,
                                        poles_p=data['_p'],poles_j=data['_j'])
    record={'mode':a.mode,'route':a.route,'poles':a.poles,'steps':a.steps,'setup_seconds':setup,
            'advance_seconds':elapsed,'seconds_per_step':elapsed/a.steps,'peak_RSS_bytes':peak,
            'modeled_explicit_pole_copy_bytes_per_call':168*a.poles if a.mode=='legacy' else 48*a.poles,
            'copy_scope':'Python trial P/J + native coefficient inputs/import/export vs native rollback P/J; excludes shared-kernel old-current copy, vector I/O, allocation initialization and coefficient arithmetic',
            'library_sha256':api.sha256,'legacy_python_sha256':hashlib.sha256(a.legacy.read_bytes()).hexdigest()}
    a.output.write_text(json.dumps(record,indent=2)+'\n')
    if hasattr(state,'close'):state.close()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('library',type=Path);p.add_argument('legacy',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--mode',choices=['legacy','persistent']);p.add_argument('--route',choices=['step','batch'],default='step')
    p.add_argument('--poles',type=int,default=8);p.add_argument('--steps',type=int,default=10000)
    a=p.parse_args()
    if a.mode:return worker(a)
    a.output.mkdir(exist_ok=False,parents=True)
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'}
    rows=[];checks=[]
    for poles in (1,8,64):
        for route in ('step','batch'):
            for rep in range(3):
                pair={}
                for mode in ('legacy','persistent'):
                    path=a.output/f'n{poles}-{route}-r{rep}-{mode}.json'
                    with path.with_suffix('.log').open('w') as log:
                        subprocess.run([sys.executable,str(Path(__file__).resolve()),str(a.library),str(a.legacy),str(path),
                            '--mode',mode,'--route',route,'--poles',str(poles),'--steps',str(a.steps)],
                            env=env,stdout=log,stderr=subprocess.STDOUT,timeout=60,check=True)
                    rows.append(json.loads(path.read_text()))
                    with np.load(path.with_suffix('.npz'),allow_pickle=False) as data:pair[mode]={k:data[k].copy() for k in data.files}
                checks.append({'poles':poles,'route':route,'repeat':rep,
                    'bitwise_equal':all(pair['legacy'][k].tobytes()==pair['persistent'][k].tobytes() for k in pair['legacy']),
                    'max_difference':max(float(np.max(abs(pair['legacy'][k]-pair['persistent'][k]))) for k in pair['legacy'])})
    cases=[]
    for poles in (1,8,64):
        for route in ('step','batch'):
            selected={m:[r for r in rows if (r['poles'],r['route'],r['mode'])==(poles,route,m)] for m in ('legacy','persistent')}
            t={m:float(np.median([r['advance_seconds'] for r in selected[m]])) for m in selected}
            cases.append({'poles':poles,'route':route,'legacy_seconds':t['legacy'],'persistent_seconds':t['persistent'],
                          'speedup':t['legacy']/t['persistent'],
                          'peak_RSS_bytes':{m:max(r['peak_RSS_bytes'] for r in selected[m]) for m in selected}})
    summary={'steps':a.steps,'repeats':3,'cases':cases,'comparisons':checks,
             'scope':'macOS isolated synthetic single-thread workers; setup separate, no GS or remote communication'}
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
    if not all(c['bitwise_equal'] for c in checks):raise SystemExit('numerical mismatch')

if __name__=='__main__':main()
