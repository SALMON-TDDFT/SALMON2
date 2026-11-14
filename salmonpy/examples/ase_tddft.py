"""ASE Atoms -> neutral GS -> A-driven current, using the small Si starter.

python -m salmonpy.examples.ase_tddft RUNTIME_DIRECTORY NEW_OUTPUT_DIRECTORY
This is a connection test, not converged optical physics. One worker/thread.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from salmonpy.ase import BOHR_ANGSTROM, tddft_from_atoms
from salmonpy.crystal import respond
from salmonpy.dft import Pseudopotential, SCFSettings, GroundStateResult
from salmonpy.fields import AcosPulse, FieldSamples
from salmonpy.runtime import load_runtime
from salmonpy.examples.si_intro import _assets, _thread_environment, compare_result


def run(runtime, directory, *, verify_independence=False):
    from ase import Atoms
    assets = _assets()  # Verify the packaged starter input/reference hashes.
    configuration = json.loads(assets['input.json'])
    expected = json.loads(assets['expected.json'])
    engine = load_runtime(runtime, required=('dft', 'tddft', 'pz', 'single-rank'))
    settings = configuration['dft']
    atoms = Atoms(numbers=settings['numbers'],
                  cell=np.array(settings['cell'])*BOHR_ANGSTROM,
                  scaled_positions=settings['fractional_positions'], pbc=True)
    pseudo = Pseudopotential(14, 2, '_rps.dat', assets['Si_rps.dat'])
    rt = configuration['rt']
    drive = AcosPulse(**rt['pulse']).sample(np.arange(rt['steps']+1)*rt['dt'])
    output = Path(directory).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = {'schema': 'salmonpy.ase-si/1', 'setup': {'status': 'passed'},
              'calculation': {'status': 'running'}, 'validation': {'status': 'not_run'},
              'resources': {'threads': 1, 'simultaneous_workers': 1,
                            'peak_memory': {'status': 'not_measured', 'bytes': None},
                            'memory_enforced': False},
              'scope': 'fixed short Si connection regression; not converged optics',
              'executable_sha256': engine.sha256}
    started = time.monotonic()
    try:
        with _thread_environment(1):
            with tddft_from_atoms(atoms, output/'gs', executable=engine,
                                  pseudopotentials=[pseudo], states=settings['states'],
                                  grid=settings['grid'], kmesh=settings['kmesh'],
                                  scf=SCFSettings(**configuration['scf']), max_steps=rt['steps'],
                                  threads=1) as initial:
                # electrons is intentionally omitted: 8 Si atoms * 4 valence electrons.
                response = respond(initial, drive, output_directory=output/'rt')
                if verify_independence:
                    # A different history must not contaminate a subsequent replay.
                    other = respond(initial, FieldSamples(drive.time, -drive.values, kind='A'))
                    replay = respond(initial, drive)
                    np.testing.assert_array_equal(replay.current, response.current)
                    if np.array_equal(other.current, response.current):
                        raise AssertionError('opposite drives produced identical currents')
                    report['independent_replay_exact'] = True
        gs = GroundStateResult.open(output/'gs')  # Also verifies GS unchanged by RT.
        np.savez(output/'current.npz', time=response.time, A=drive.values, current=response.current)
        comparison = compare_result(expected, gs.energy_hartree, response.time, response.current)
        report['calculation'] = {'status': 'passed', 'energy_hartree': gs.energy_hartree,
                                 'iterations': gs.iterations, 'samples': len(response.time),
                                 'electrons': sum(pseudo.valence_electrons for _ in atoms),
                                 'current_convention': response.field_convention}
        report['validation'] = {'status': 'passed' if comparison['passed'] else 'failed', **comparison}
        if not comparison['passed']:
            raise AssertionError('Si starter reference comparison failed; inspect report.json')
        return report
    except BaseException as exc:
        if report['calculation']['status'] == 'running':
            report['calculation'] = {'status': 'failed', 'error': str(exc)}
        raise
    finally:
        report['resources']['wall_seconds'] = time.monotonic()-started
        (output/'report.json').write_text(json.dumps(report, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runtime')
    parser.add_argument('directory')
    parser.add_argument('--verify-independence', action='store_true')
    args = parser.parse_args()
    print(json.dumps(run(args.runtime, args.directory, verify_independence=args.verify_independence), indent=2))


if __name__ == '__main__':
    main()
