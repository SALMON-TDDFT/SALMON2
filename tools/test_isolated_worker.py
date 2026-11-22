"""Compare explicit experimental E-worker against an isolated batch reference.

Requires output from test_isolated_reference.py and a new output directory.
All protocol values are atomic units. This is not a public current-density API.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np


def run(executable, reference, directory):
    executable = Path(executable).resolve(strict=True)
    reference = Path(reference).resolve(strict=True)
    root = Path(directory).resolve()
    root.mkdir(parents=True, exist_ok=False)
    data = np.loadtxt(reference / 'first/C2H2_rt.data')
    dt = 0.001 / 0.02418884326505
    fields = data[:, 10:13] / (27.21138505 / 0.52917721067)
    contents = (reference / 'first/input.inp').read_text()
    env = dict(os.environ)
    env.update(SALMONPY_TDDFT_WORKER='E2', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    trajectories = []
    for label in ('first', 'repeat'):
        folder = root / label
        folder.mkdir()
        (folder / 'input.inp').write_text(contents)
        for name in ('C_rps.dat', 'H_rps.dat'):
            shutil.copy2(reference / 'gs' / name, folder / name)
        shutil.copytree(reference / 'gs/data_for_restart', folder / 'data_for_restart')
        # Existing native input delimiter separates namelist and worker commands.
        commands = [f'INIT_E {dt:.17e} {len(fields)} 0 0 0']
        commands += ['STEP_E '+ ' '.join(f'{v:.17e}' for v in E) for E in fields]
        commands += ['CLOSE']
        result = subprocess.run([str(executable)], cwd=folder,
                                input=contents+'\n!SALMONPY_END_INPUT\n'+'\n'.join(commands)+'\n', text=True, capture_output=True,
                                env=env, timeout=60)
        (folder / 'worker.log').write_text(result.stdout+result.stderr)
        lines = [line.split() for line in result.stdout.splitlines()
                 if line.startswith('SALMONPY_TDDFT_E_V2 ')]
        if result.returncode or len(lines) != len(fields)+1:
            raise AssertionError(f'E-worker handshake/steps failed: {folder / "worker.log"}')
        assert [int(line[1]) for line in lines] == list(range(len(fields)+1))
        values = np.asarray([[float(v) for v in line[2:]] for line in lines])
        assert values.shape == (len(fields)+1, 7)
        assert np.isfinite(values).all()
        trajectories.append(values)
    np.testing.assert_array_equal(*trajectories)
    # Worker returns electronic dipole change and total charge dipole, in e*bohr.
    expected = data[:, 13:19] / 0.52917721067
    np.testing.assert_allclose(trajectories[0][1:,:6], expected, rtol=1.e-6, atol=1.e-10)
    from isolated_worker_driver import IsolatedWorker
    from salmonpy.field_input import FieldAdapter, FieldContract, FieldPacket, FieldValue
    from salmonpy.timing import TimeSpec, InputPort
    from salmonpy.errors import InputError
    contract = FieldContract(TimeSpec(dt, (InputPort('E'),), ()),
                             units={'E':'E:au'}, position=(0,0,0), gauge='length')
    def packet(E, n):
        return FieldPacket(E=FieldValue(E, time=n*dt, position=(0,0,0), unit='E:au'))
    def prepared(label):
        folder = root / label
        folder.mkdir()
        for name in ('C_rps.dat', 'H_rps.dat'):
            shutil.copy2(reference / 'gs' / name, folder / name)
        shutil.copytree(reference / 'gs/data_for_restart', folder / 'data_for_restart')
        return folder
    # Two live native states, opposite drives, alternating calls. Advancing A
    # cannot change B's clock/observable. Compare both against standalone runs.
    with IsolatedWorker(executable, prepared('interleaved-a'), contents, dt=dt, capacity=20) as a, \
         IsolatedWorker(executable, prepared('interleaved-b'), contents, dt=dt, capacity=20) as b:
        aa, bb = FieldAdapter(a, contract), FieldAdapter(b, contract)
        old = a.dipoles.copy()
        try:
            aa.step_fields(FieldPacket())
        except InputError:
            pass
        else:
            raise AssertionError('missing E accepted')
        np.testing.assert_array_equal(a.dipoles, old)
        assert a.steps == 0
        av, bv = [], []
        for n, E in enumerate(fields, 1):
            old = b.dipoles.copy()
            av.append(aa.step_fields(packet(E, n)))
            np.testing.assert_array_equal(b.dipoles, old)
            assert b.steps == n-1
            bv.append(bb.step_fields(packet(-E, n)))
        np.testing.assert_array_equal(av, trajectories[0][1:])
        if np.array_equal(av, bv):
            raise AssertionError('opposite drive did not change the response')
    with IsolatedWorker(executable, prepared('negative-alone'), contents, dt=dt, capacity=20) as a:
        aa = FieldAdapter(a, contract)
        negative = [aa.step_fields(packet(-E, n)) for n, E in enumerate(fields, 1)]
        np.testing.assert_array_equal(negative, bv)
    # Reconstruction from GS plus accepted E history; not a public snapshot API.
    with IsolatedWorker(executable, prepared('replay'), contents, dt=dt, capacity=20) as a:
        aa = FieldAdapter(a, contract)
        for n, E in enumerate(fields[:10], 1):
            aa.step_fields(packet(E, n))
        continued = [aa.step_fields(packet(E, n)) for n, E in enumerate(fields[10:], 11)]
        np.testing.assert_array_equal(continued, trajectories[0][11:])
    report = dict(schema='salmonpy.isolated-worker-validation/2',
                  capabilities=['isolated','pz','fixed-ions','middlepoint','single-rank','E-to-dipole','measured-electron-number'],
                  interleaved_independence=True, field_packet_missing_E_rejected=True,
                  replay_continuation_max_abs_error=0., executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                  protocol='SALMONPY_TDDFT_E_V2', steps=len(fields),
                  repeat_max_abs_error=0.,
                  reference_max_abs_error=float(np.max(np.abs(trajectories[0][1:,:6]-expected))))
    np.savez(root / 'dipoles.npz', fields=fields, dt=dt, values=trajectories[0])
    (root / 'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('reference')
    parser.add_argument('directory')
    args = parser.parse_args()
    run(args.executable, args.reference, args.directory)
