"""Validate public TDDFT replay with an audited runtime and existing input bundle.

Usage: PYTHONPATH=. python3 tools/test_tddft_snapshots.py RUNTIME BUNDLE INPUT OUTPUT
The bundle may contain synthetic normalized orbitals (not a converged GS).
No GS files are collected in OUTPUT. No native build is performed.
"""
import json
from pathlib import Path
import sys
import subprocess
import numpy as np
from salmonpy.runtime import load_runtime
from salmonpy.tddft import TDDFTInitialState
from salmonpy.errors import SnapshotError

if len(sys.argv)>1 and sys.argv[1]=='--restore-child':
    runtime,bundle,inp,checkpoint,values_file,output=sys.argv[2:]
    exe=load_runtime(runtime,required=['tddft'])
    with TDDFTInitialState.from_salmon(bundle,inp,executable=exe,
                                      max_steps=25,threads=1,timeout=60) as model:
        with model.new_state(dt=.001) as state:
            state.restore_checkpoint(checkpoint)
            with np.load(values_file,allow_pickle=False) as saved:
                values=saved['A'].copy();expected=saved['original'].copy()
            actual=np.array([state.step(v) for v in values[5:]])
            assert actual.tobytes()==expected.tobytes()
            np.save(output,actual)
    raise SystemExit(0)

runtime, bundle, inp, output = sys.argv[1:]
output = Path(output); output.mkdir(parents=True, exist_ok=False)
exe = load_runtime(runtime, required=['tddft'])  # verifies recorded checksums
values = np.random.default_rng(20261114).normal(size=(25, 3))*1e-3
dt = .001
with TDDFTInitialState.from_salmon(bundle, inp, executable=exe,
                                  max_steps=25, threads=1, timeout=60) as model:
    with model.new_state(dt=dt) as a, model.new_state(dt=dt) as b:
        initial = a.snapshot()
        for v in values[:5]: a.step(v)
        saved = a.snapshot()
        a.save_checkpoint(output/'checkpoint')
        original = np.array([a.step(v) for v in values[5:]])
        # Diverge a live branch before restoring and test an independent branch.
        a.restore(saved); b.restore(saved)
        assert a.time == b.time == 5*dt
        restored = np.array([a.step(v) for v in values[5:]])
        branched = np.array([b.step(v) for v in values[5:]])
        np.testing.assert_array_equal(restored, original)
        np.testing.assert_array_equal(branched, original)
        assert restored.tobytes() == branched.tobytes() == original.tobytes()
        a.restore(initial)
        start_replayed = np.array([a.step(v) for v in values])
        np.testing.assert_array_equal(start_replayed[5:], original)
        assert start_replayed[5:].tobytes() == original.tobytes()
        with model.new_state(dt=dt*2) as wrong:
            with np.testing.assert_raises(SnapshotError): wrong.restore(saved)
        metadata = {'dt': dt, 'saved_steps': 5, 'continued_steps': 20,
                    'threads': 1, 'executable_sha256': exe.sha256,
                    'max_abs_current_difference': float(np.max(np.abs(original-restored))),
                    'bitwise_equal': True,
                    'fixture': 'caller-supplied bundle; no claim of GS convergence'}
        np.savez(output/'currents.npz', A=values, original=original,
                 restored=restored, independent=branched)
        (output/'summary.json').write_text(json.dumps(metadata, indent=2)+'\n')
        print(json.dumps(metadata, indent=2))

# All original states/model are closed before another Python process restores.
subprocess.run([sys.executable,str(Path(__file__).resolve()),'--restore-child',
                runtime,bundle,inp,str(output/'checkpoint'),str(output/'currents.npz'),
                str(output/'disk-restored.npy')],check=True,timeout=300)
disk=np.load(output/'disk-restored.npy',allow_pickle=False)
assert disk.tobytes()==original.tobytes()
metadata.update(disk_fresh_process_bitwise_equal=True,
                disk_max_abs_current_difference=float(np.max(np.abs(disk-original))))
(output/'summary.json').write_text(json.dumps(metadata,indent=2)+'\n')
print(json.dumps(metadata,indent=2))
