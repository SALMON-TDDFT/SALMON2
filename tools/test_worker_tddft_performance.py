"""Periodic worker performance vs checked and previous validated executable.

python tools/test_worker_tddft_performance.py WORKER_RECEIPT BASELINE_HSE_BUILD FIXTURE NEW_OUT XC
Only fixture storage is substituted: a test copier links existing GS instead of
collecting copies. Production copier is unchanged. Actual worker initialization,
updates, handshake, explicit replay and failure retirement use production code.
XC selects a synthetic update regression, not a converged material response.
"""
import argparse,hashlib,json,re,shutil,tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
from salmonpy.tddft import TDDFTInitialState
from salmonpy.fdtd import SalmonExecutable
from salmonpy.errors import RuntimeFailure
p=argparse.ArgumentParser(description=__doc__)
for name in ('worker','baseline','fixture','output','xc'):p.add_argument(name)
a=p.parse_args()
if a.xc not in ('pz','tbmbj','libxc_pz','libxc_pbe'):raise ValueError('unsupported test XC')
root=Path(a.output);root.mkdir(parents=True,exist_ok=False);fixture=Path(a.fixture).resolve()
new=json.loads(Path(a.worker).read_text());old=json.loads(Path(a.baseline).read_text())
new_exe=Path(new['executable']);old_exe=Path(a.baseline).parent/'build/salmon'
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
assert new['status']=='built-and-object-audited' and digest(new_exe)==new['sha256']
assert digest(old_exe)==old['origin_sha256']
gs=fixture/'data_for_restart';before={p.name:digest(p) for p in gs.iterdir() if p.is_file()}

def links_only(source,dest):
    dest.mkdir();records=[]
    for name in ('rt.inp','pseudo1_rps.dat','sym.dat'):
        if (source/name).is_file():
            shutil.copy2(source/name,dest/name);records.append((name,digest(dest/name)))
    target=(source/'data_for_restart').resolve();(dest/'data_for_restart').symlink_to(target,target_is_directory=True)
    for path in sorted(target.rglob('*')):
        if path.is_file():records.append(('data_for_restart/'+str(path.relative_to(target)),digest(path)))
    return tuple(sorted(records))

with tempfile.TemporaryDirectory(prefix='worker-perf-fixture-') as d:
    bundle=Path(d)/'input';bundle.mkdir()
    (bundle/'rt.inp').write_text(re.sub(r"xc\s*=\s*'[^']+'","xc='"+a.xc+"'",(fixture/'rt.inp').read_text()))
    for name in ('pseudo1_rps.dat','sym.dat'):
        if (fixture/name).is_file():shutil.copy2(fixture/name,bundle/name)
    (bundle/'data_for_restart').symlink_to(gs,target_is_directory=True)
    with patch('salmonpy.tddft._copy_bundle',side_effect=links_only):
        old_model=TDDFTInitialState.from_salmon(bundle,'rt.inp',executable=SalmonExecutable(old_exe,old['origin_sha256']),max_steps=20,threads=1,timeout=30)
        model=TDDFTInitialState.from_salmon(bundle,'rt.inp',executable=SalmonExecutable(new_exe,new['sha256']),max_steps=20,threads=1,timeout=30)
        states=[]
        try:
            baseline=old_model.new_state(dt=.001);checked=model.new_state(dt=.001);fast=model.new_state(dt=.001,execution_mode='performance');states=[baseline,checked,fast]
            traces=[[] for s in states];norms=[[] for s in states]
            def drive(i):return np.array([.001*np.sin(.2*i),.0007*(np.cos(.13*i)-1),.0009*np.sin(.17*i)])
            for i in range(1,17):
                for index,state in enumerate(states):
                    traces[index].append(state.step(drive(i)));norms[index].append(state.electron_number)
                if i==8:checked_saved=checked.snapshot();fast_saved=fast.snapshot()
            for index in (1,2):
                np.testing.assert_array_equal(traces[index],traces[0]);np.testing.assert_array_equal(norms[index],norms[0])
            # Physical snapshot excludes the mode: each receiving state retains its mode.
            checked.restore(fast_saved);fast.restore(checked_saved)
            for state in (checked,fast):
                for i in range(9,17):
                    np.testing.assert_array_equal(state.step(drive(i)),traces[0][i-1])
                    assert state.electron_number==norms[0][i-1]
            assert checked.execution_mode=='checked' and fast.execution_mode=='performance'
            for label,state in zip(('baseline','checked','performance'),states):
                shutil.copy2(Path(state._work.name)/'run/worker.log',root/(label+'.log'))
            fast._proc.kill();fast._proc.wait()
            with patch.object(fast,'_spawn',side_effect=AssertionError('automatic replay forbidden')):
                try:fast.step(drive(17))
                except RuntimeFailure:pass
                else:raise AssertionError('killed worker accepted update')
            assert fast._closed
            # An old worker must reject the requested mode, never silently use checked.
            try:old_model.new_state(dt=.001,execution_mode='performance')
            except RuntimeFailure as error:
                assert 'INIT' in str(error);(root/'old-worker-rejection.txt').write_text(str(error))
            else:raise AssertionError('old worker accepted performance mode')
            np.savez(root/'traces.npz',current=np.array(traces[0]),electron_number=np.array(norms[0]))
        finally:
            for state in states:state.close()
            model.close();old_model.close()
assert before=={p.name:digest(p) for p in gs.iterdir() if p.is_file()}
report={'xc':a.xc,'steps':16,'current_and_electron_number_vs_checked_and_previous':'array_equal',
 'cross_mode_explicit_replay_restore':'array_equal, receiving mode preserved','killed_worker_no_automatic_respawn':True,
 'old_worker_rejects_mode':True,'GS_bytes_unchanged':True,'storage':'test-only GS links; production copier unchanged; no GS collected',
 'scope':'synthetic periodic worker regression; not molecular, physical convergence or speed benchmark'}
(root/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
