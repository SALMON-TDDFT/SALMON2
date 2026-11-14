"""Validate public TDDFT replay with an audited runtime and existing input bundle.

Usage: PYTHONPATH=. python3 tools/test_tddft_snapshots.py RUNTIME BUNDLE INPUT OUTPUT
The bundle may contain synthetic normalized orbitals (not a converged GS).
No GS files are collected in OUTPUT. No native build is performed.
"""
import json
from pathlib import Path
import sys
import numpy as np
from salmonpy.runtime import load_runtime
from salmonpy.tddft import TDDFTInitialState
from salmonpy.errors import SnapshotError

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
