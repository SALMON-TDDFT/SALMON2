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
    env.update(SALMONPY_TDDFT_WORKER='E3', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
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
                 if line.startswith('SALMONPY_TDDFT_E_V3 ')]
        if result.returncode or len(lines) != len(fields)+1:
            raise AssertionError(f'E-worker handshake/steps failed: {folder / "worker.log"}')
        assert [int(line[1]) for line in lines] == list(range(len(fields)+1))
        values = np.asarray([[float(v) for v in line[2:]] for line in lines])
        assert values.shape == (len(fields)+1, 9)
        assert np.isfinite(values).all()
        trajectories.append(values)
    np.testing.assert_array_equal(*trajectories)
    # Worker returns electronic dipole change and total charge dipole, in e*bohr.
    expected = data[:, 13:19] / 0.52917721067
    np.testing.assert_allclose(trajectories[0][1:,:6], expected, rtol=1.e-6, atol=1.e-10)
    energy_reference=np.loadtxt(reference/'first/C2H2_rt_energy.data')
    energy_indices=np.rint(energy_reference[:,0]/.001).astype(int)
    assert np.all(trajectories[0][energy_indices,8]==1.)
    expected_energy=energy_reference[:,1]/27.21138505
    np.testing.assert_allclose(trajectories[0][energy_indices,7],expected_energy,rtol=0,atol=1.e-10)
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
    run_savepoint(executable,reference,root/'savepoint-proof',saved_step=7)
    proof=json.loads((root/'savepoint-proof/report.json').read_text())
    disk=run_direct_disk(executable,reference,root/'direct-disk-proof')
    report = dict(schema='salmonpy.isolated-worker-validation/5',
                  direct_disk_process_exit_exact=all(v['continuation_exact'] and v['output_streams_exact'] for v in disk['cases']),
                  direct_disk_extra_evolutions=0,
                  savepoint_same_process=proof['same_process'],
                  savepoint_extra_evolutions=proof['restore_additional_evolution_calls'],
                  savepoint_stale_rejected=proof['stale_savepoint_rejected'],
                  savepoint_continuation_max_abs_error=0.,
                  capabilities=['isolated','pz','fixed-ions','middlepoint','single-rank','E-to-dipole','measured-electron-number','internal-energy-freshness','live-savepoint','direct-disk-state'],
                  interleaved_independence=True, field_packet_missing_E_rejected=True,
                  replay_continuation_max_abs_error=0., executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                  protocol='SALMONPY_TDDFT_E_V3', steps=len(fields),
                  repeat_max_abs_error=0.,
                  reference_energy_max_abs_error=float(np.max(abs(trajectories[0][energy_indices,7]-expected_energy))),
                  reference_max_abs_error=float(np.max(np.abs(trajectories[0][1:,:6]-expected))))
    np.savez(root / 'dipoles.npz', fields=fields, dt=dt, values=trajectories[0])
    (root / 'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


def run_direct_disk(executable, reference, directory, *, worker_mode="E3"):
    """Raw-worker process-exit proof; requires no pre-existing runtime receipt."""
    from isolated_worker_driver import IsolatedWorker
    from salmonpy.field_input import FieldPacket, FieldValue
    root=Path(directory);root.mkdir(parents=True,exist_ok=False)
    reference=Path(reference);contents=(reference/'first/input.inp').read_text()
    records=[]
    for saved_step in (0,7,8,20):
        case=root/str(saved_step);case.mkdir()
        def prepare(name):
            folder=case/name;folder.mkdir()
            for asset in ('C_rps.dat','H_rps.dat'):shutil.copy2(reference/'gs'/asset,folder/asset)
            shutil.copytree(reference/'gs/data_for_restart',folder/'data_for_restart')
            if worker_mode=='E3M':shutil.copy2(reference/'first/salmonpy-near-eh.dat',folder/'salmonpy-near-eh.dat')
            if worker_mode=='E3Y':shutil.copy2(reference/'first/salmonpy-ylm.dat',folder/'salmonpy-ylm.dat')
            return folder
        def step(worker,n):
            vector=((.01,.003,.02),(.02,-.004,-.01),(-.01,.005,.006))[(n-1)%3]
            return worker.step_fields(FieldPacket(E=FieldValue(vector,time=n*.003,position=(0,0,0),unit='E:au')))
        first=prepare('first')
        with IsolatedWorker(executable,first,contents,dt=.003,capacity=20,worker_mode=worker_mode) as a:
            for n in range(1,saved_step+1):step(a,n)
            pid=a._proc.pid
            assert a.control('SAVE_E')==(6,0,saved_step,saved_step)
            saved=a.dipoles.copy()
            streams={f.name:f.read_bytes() for f in first.glob('*.data')}
            expected=[step(a,n) for n in range(saved_step+1,21)]
        expected_streams={f.name:f.read_bytes() for f in first.glob('*.data')}
        second=prepare('second');shutil.copyfile(first/'isolated-state.bin',second/'isolated-state.bin')
        with IsolatedWorker(executable,second,contents,dt=.003,capacity=20,worker_mode=worker_mode) as b:
            assert b._proc.pid!=pid
            assert b.control('LOAD_E')==(7,0,saved_step,0)
            for name,raw in streams.items():(second/name).write_bytes(raw)
            assert b.control('RELEASE_E')==(3,0,saved_step,0)
            b.steps=saved_step;b.time=saved_step*.003;b.dipoles=saved
            actual=[step(b,n) for n in range(saved_step+1,21)]
            np.testing.assert_array_equal(actual,expected)
            assert b.control('PROBE_E')==(5,0,20,20-saved_step)
        assert set(expected_streams)=={f.name for f in second.glob('*.data')}
        for name,raw in expected_streams.items():assert raw==(second/name).read_bytes()
        records.append(dict(saved_step=saved_step,new_process=True,restore_extra_evolutions=0,continuation_exact=True,output_streams_exact=True))
    report=dict(schema='salmonpy.isolated-direct-validation/1',cases=records)
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def run_savepoint(executable,reference,directory,*,saved_step=7,worker_mode="E3"):
    from isolated_worker_driver import IsolatedWorker
    if type(saved_step) is not int or not 0<=saved_step<=20:raise ValueError("saved step must be in 0..20")
    from salmonpy.field_input import FieldPacket,FieldValue
    executable=Path(executable).resolve(strict=True);reference=Path(reference).resolve(strict=True)
    root=Path(directory).resolve();root.mkdir(parents=True,exist_ok=False)
    contents=(reference/'first/input.inp').read_text()
    def prepare(name):
        folder=root/name;folder.mkdir()
        for name in ('C_rps.dat','H_rps.dat'):shutil.copy2(reference/'gs'/name,folder/name)
        shutil.copytree(reference/'gs/data_for_restart',folder/'data_for_restart')
        if worker_mode=='E3M':shutil.copy2(reference/'first/salmonpy-near-eh.dat',folder/'salmonpy-near-eh.dat')
        if worker_mode=='E3Y':shutil.copy2(reference/'first/salmonpy-ylm.dat',folder/'salmonpy-ylm.dat')
        return folder
    def step(worker,n):
        samples=((.01,.003,.02),(.02,-.004,-.01),(-.01,.005,.006))
        return worker.step_fields(FieldPacket(E=FieldValue(samples[(n-1)%3],time=n*.003,position=(0,0,0),unit='E:au')))
    with IsolatedWorker(executable,prepare('direct'),contents,dt=.003,capacity=20,worker_mode=worker_mode) as a, \
         IsolatedWorker(executable,prepare('control'),contents,dt=.003,capacity=20,worker_mode=worker_mode) as b:
        original_pid=a._proc.pid
        for n in range(1,saved_step+1):np.testing.assert_array_equal(step(a,n),step(b,n))
        before=b.dipoles.copy();saved=a.capture()
        assert saved['reply'][2:]==(saved_step,saved_step)
        assert a.control('CAPTURE_E')[1]==14
        expected=[step(a,n) for n in range(saved_step+1,21)]
        reply,resumed=a.restore(saved)
        assert reply[2:]==(saved_step,20) and resumed[2:]==(saved_step,20)
        assert a._proc.pid==original_pid
        np.testing.assert_array_equal(b.dipoles,before)
        np.testing.assert_array_equal([step(a,n) for n in range(saved_step+1,21)],expected)
        assert a.control('PROBE_E')[2:]==(20,40-saved_step)
        a.restore(saved)
        np.testing.assert_array_equal([step(a,n) for n in range(saved_step+1,21)],expected)
        assert a.release()[2:]==(20,60-2*saved_step)
        assert a.control('RESTORE_E')[1]==14
        np.testing.assert_array_equal([step(b,n) for n in range(saved_step+1,21)],expected)
        for name in ('C2H2_rt.data','C2H2_rt_energy.data'):
            # Capture flushes both streams before comparing the completed files.
            a.capture();b.capture()
            try:a.restore(saved)
            except ValueError:pass
            else:raise AssertionError('released savepoint was accepted after reuse')
            assert (root/'direct'/name).read_bytes()==(root/'control'/name).read_bytes(),name
            a.release();b.release()
    report=dict(schema='salmonpy.isolated-native-savepoint-validation/1',
        executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
        same_process=True,restore_additional_evolution_calls=0,continuation_bit_exact=True,stale_savepoint_rejected=True,
        independent_state_unchanged=True,output_streams_bit_exact=True,
        saved_step=saved_step,capacity=20,scope='single pending live-worker PZ memory savepoint; not public Python/disk checkpoint')
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('reference')
    parser.add_argument('directory')
    parser.add_argument('--savepoint',action='store_true')
    parser.add_argument('--savepoint-step',type=int,default=7)
    args = parser.parse_args()
    if args.savepoint:run_savepoint(args.executable,args.reference,args.directory,saved_step=args.savepoint_step)
    else:run(args.executable,args.reference,args.directory)
