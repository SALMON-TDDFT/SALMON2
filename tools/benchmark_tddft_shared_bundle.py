"""Isolated copy/lifetime benchmark with opaque synthetic payloads, no physical GS.

python tools/benchmark_tddft_shared_bundle.py LEGACY_PY NEW_OUTPUT
The deterministic Python protocol worker tests replay/control, not SALMON physics.
"""
import argparse,hashlib,importlib.util,json,os,subprocess,sys,tempfile,time
from pathlib import Path
import numpy as np

WORKER='''import sys
j=[0.,0.,0.];n=0
for line in sys.stdin:
 p=line.split()
 if not p:continue
 if p[0]=='INIT':dt=float(p[1]);j=list(map(float,p[3:6]))
 elif p[0]=='STEP':
  n+=1;a=list(map(float,p[1:]));j=[.83*x+dt*y+.001*n for x,y in zip(j,a)]
 elif p[0]=='CLOSE':break
 else:continue
 print('SALMONPY_TDDFT_V1',n,*j,flush=True)
'''


def worker(a):
    import salmonpy.tddft as current
    from salmonpy.fdtd import SalmonExecutable
    if a.mode=='legacy':
        name='salmonpy._tddft_copy_baseline';spec=importlib.util.spec_from_file_location(name,a.legacy)
        module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    else:module=current
    with tempfile.TemporaryDirectory(prefix='salmonpy-copy-benchmark-') as temporary:
        root=Path(temporary);bundle=root/'input';restart=bundle/'data_for_restart';restart.mkdir(parents=True)
        for name in ('info.bin','occupation.bin'):(restart/name).write_bytes(b'opaque protocol fixture, not SALMON GS')
        payload=np.random.default_rng(701).bytes(a.payload_mib*1024*1024)
        (restart/'wfn.bin').write_bytes(payload);del payload
        (bundle/'rt.inp').write_text('&calculation\n/\n')
        exe=root/'worker.py';exe.write_text(WORKER);exe=exe.resolve()
        # Explicit interpreter, portable to systems without executable shebang scripts.
        original=module.subprocess.Popen
        def launch(args,**kw):
            if args==[str(exe)]:args=[sys.executable,str(exe)]
            return original(args,**kw)
        module.subprocess.Popen=launch
        engine=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest())
        initial=None;states=[]
        try:
            start=time.perf_counter();initial=module.TDDFTInitialState.from_salmon(bundle,'rt.inp',executable=engine,max_steps=20)
            loader=time.perf_counter()-start;start=time.perf_counter()
            for _ in range(a.states):states.append(initial.new_state(dt=.02))
            setup=time.perf_counter()-start
            roots={initial._snapshot}
            roots.update(s._bundle for s in states);roots.update(Path(s._work.name)/'run' for s in states)
            files=[r/'data_for_restart/wfn.bin' for r in roots]
            logical=sum(p.stat().st_size for p in files)
            allocated=sum(p.stat().st_blocks*512 for p in files) if hasattr(files[0].stat(),'st_blocks') else None
            digest=hashlib.sha256(files[0].read_bytes()).hexdigest()
            assert all(hashlib.sha256(p.read_bytes()).hexdigest()==digest for p in files)
            initial.close()
            values=np.random.default_rng(71).normal(size=(12,3))*.01
            histories=[];tokens=[]
            for s in states:
                histories.append(np.array([s.step(v) for v in values[:4]]));tokens.append(s.snapshot())
            for s in states:
                for v in values[4:8]:s.step(v)
            start=time.perf_counter()
            for s,token in zip(states,tokens):s.restore(token)
            restore=time.perf_counter()-start
            replay=np.stack([[s.step(v) for v in values[4:]] for s in states])
            try:
                import resource
                raw=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                rss=int(raw if sys.platform=='darwin' else raw*1024)
            except (ImportError,OSError,ValueError):rss=None
            np.savez(a.output.with_suffix('.npz'),initial=np.asarray(histories),restored=replay)
            result={'mode':a.mode,'states':a.states,'synthetic_payload_bytes':a.payload_mib*1024*1024,
                    'managed_payload_copies':len(files),'managed_logical_payload_bytes':logical,
                    'summed_st_blocks_bytes':allocated,'allocation_scope':'sum of st_blocks; filesystem shared extents are not deduplicated',
                    'loader_seconds':loader,'state_creation_seconds':setup,'restore_all_states_seconds':restore,
                    'worker_processes':len(states),'parent_peak_RSS_bytes':rss,
                    'payload_sha256':digest,'legacy_adapter_sha256':hashlib.sha256(a.legacy.read_bytes()).hexdigest(),
                    'scope':'opaque synthetic files + deterministic protocol workers; original input excluded from managed copies; no physical GS/TDDFT calculation'}
            a.output.write_text(json.dumps(result,indent=2)+'\n')
        finally:
            for s in states:s.close()
            if initial is not None:initial.close()
            module.subprocess.Popen=original


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('legacy',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--mode',choices=['legacy','shared']);p.add_argument('--states',type=int,default=4)
    p.add_argument('--payload-mib',type=int,default=16)
    a=p.parse_args()
    if a.mode:return worker(a)
    a.output.mkdir(exist_ok=False,parents=True);records=[];comparisons=[]
    for n in (1,4,8):
        for repeat in range(3):
            pair={}
            for mode in ('legacy','shared'):
                path=a.output/f'n{n}-r{repeat}-{mode}.json'
                with path.with_suffix('.log').open('w') as log:
                    subprocess.run([sys.executable,str(Path(__file__).resolve()),str(a.legacy),str(path),
                                    '--mode',mode,'--states',str(n),'--payload-mib',str(a.payload_mib)],
                                   stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60)
                records.append(json.loads(path.read_text()))
                with np.load(path.with_suffix('.npz'),allow_pickle=False) as data:pair[mode]={k:data[k].copy() for k in data.files}
            comparisons.append({'states':n,'repeat':repeat,'bitwise_equal':all(pair['legacy'][k].tobytes()==pair['shared'][k].tobytes() for k in pair['legacy']),
                                'maximum_difference':max(float(np.max(abs(pair['legacy'][k]-pair['shared'][k]))) for k in pair['legacy'])})
    cases=[]
    for n in (1,4,8):
        rows={mode:[r for r in records if (r['states'],r['mode'])==(n,mode)] for mode in ('legacy','shared')}
        c={'states':n}
        for mode in rows:
            c[mode]={key:float(np.median([r[key] for r in rows[mode]])) for key in ('loader_seconds','state_creation_seconds','restore_all_states_seconds')}
            c[mode].update({key:rows[mode][0][key] for key in ('managed_payload_copies','managed_logical_payload_bytes')})
            c[mode]['parent_peak_RSS_bytes']=max(r['parent_peak_RSS_bytes'] for r in rows[mode])
        cases.append(c)
    result={'cases':cases,'numerical_comparisons':comparisons,'scope':'three isolated repeats, synthetic16MiB payload by default; walltime includes protocol worker startup and copy/hash'}
    (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if not all(c['bitwise_equal'] for c in comparisons):raise SystemExit('numerical mismatch')

if __name__=='__main__':main()
