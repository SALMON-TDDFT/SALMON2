"""Material ID -> offline Si inputs -> existing DFT -> independent A-driven TDDFT.

Uses the packaged short Si starter, not converged optical physics. Requires an
explicit preinstalled PZ single-rank DFT/TDDFT runtime and optional ASE/f90nml.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np

from salmonpy.ase import BOHR_ANGSTROM
from salmonpy.builtin_materials import builtin_material_registry, material_from_catalog
from salmonpy.config import load_materials
from salmonpy.dft import SCFSettings
from salmonpy.fields import AcosPulse
from salmonpy.preparation import prepare_inputs
from salmonpy.runtime import load_runtime
from salmonpy.examples.si_intro import _assets, _thread_environment, compare_result


def run(runtime, directory):
    assets = _assets()
    configuration = json.loads(assets['input.json'])
    expected = json.loads(assets['expected.json'])
    engine = load_runtime(runtime, required=('dft', 'tddft', 'pz', 'single-rank'))
    output = Path(directory).resolve()
    output.mkdir(parents=True, exist_ok=False)
    source = output/'source'
    source.mkdir()
    settings = configuration['dft']
    cell = np.asarray(settings['cell'])*BOHR_ANGSTROM
    positions = np.asarray(settings['fractional_positions'])@cell
    # Explicit full-precision extxyz, avoiding an additional decimal truncation.
    text = str(len(positions))+'\nLattice="'+' '.join(f'{x:.17g}' for x in cell.ravel())
    text += '" Properties=species:S:1:pos:R:3 pbc="T T T"\n'
    text += ''.join('Si '+' '.join(f'{x:.17g}' for x in position)+'\n' for position in positions)
    (source/'si.xyz').write_text(text)
    (source/'Si_rps.dat').write_bytes(assets['Si_rps.dat'])
    (source/'LICENSE.txt').write_bytes(assets['LICENSE.txt'])
    recipe = {'schema': 'salmonpy.materials/1', 'materials': {'Si-starter': {
        'backend': 'salmonpy.dft', 'model': 'dft', 'parameters': {
            'structure': {'path': 'si.xyz', 'format': 'extxyz'},
            'pseudopotentials': [{'path': 'Si_rps.dat', 'atomic_number': 14, 'local_channel': 2}],
            'xc': settings['xc'], 'states': settings['states'], 'kmesh': settings['kmesh'],
            'grid': {'spacing': {'value': (np.diag(cell)/settings['grid']).tolist(), 'unit': 'angstrom'}}}}}}
    registry = builtin_material_registry()
    bundle = prepare_inputs(recipe, output/'inputs', registry=registry, base_directory=source,
                            files={'si.xyz': 'si.xyz', 'Si_rps.dat': 'Si_rps.dat', 'LICENSE.txt': 'LICENSE.txt'},
                            provenance={'si.xyz': {'source': 'packaged Si starter geometry'},
                                        'Si_rps.dat': json.loads(assets['provenance.json']),
                                        'LICENSE.txt': {'source': 'packaged upstream SALMON license'}})
    catalog = load_materials(bundle/'materials.json', registry=registry)
    material = material_from_catalog(catalog, 'Si-starter', base_directory=bundle)
    rt = configuration['rt']
    drive = AcosPulse(**rt['pulse']).sample(np.arange(rt['steps']+1)*rt['dt'])
    report = {'scope': 'short Si connection regression, not converged optical physics',
              **material.record(), 'status': 'running', 'executable_sha256': engine.sha256,
              'resources': {'threads_per_worker': 1, 'max_simultaneous_workers': 2,
                            'memory_enforced': False, 'peak_memory': None}}
    started = time.monotonic()
    try:
        with _thread_environment(1):
            gs = material.run(engine, output/'gs', scf=SCFSettings(**configuration['scf']), threads=1)
            with gs.tddft_initial_state(executable=engine, max_steps=rt['steps'], threads=1) as initial:
                with initial.new_state(dt=rt['dt']) as a, initial.new_state(dt=rt['dt']) as b:
                    currents = [a.current.copy()]
                    for field in drive.values[1:]:
                        before = b.current.copy()
                        a.step(field)
                        np.testing.assert_array_equal(before, b.current)
                        b.step(-field)
                        currents.append(a.current.copy())
                    assert a.material_id == b.material_id == initial.material_id == gs.material_id == material.material_id
                    if np.array_equal(a.current, b.current):
                        raise AssertionError('opposite drives produced identical currents')
                with initial.new_state(dt=rt['dt']) as replay:
                    values = [replay.current.copy()]
                    for field in drive.values[1:]:
                        values.append(replay.step(field).copy())
                    np.testing.assert_array_equal(values, currents)
        comparison = compare_result(expected, gs.energy_hartree, drive.time, np.asarray(currents))
        report.update(status='passed' if comparison['passed'] else 'failed',
                      comparison=comparison, energy_hartree=gs.energy_hartree,
                      independent_replay_exact=True, grid=list(material.grid))
        np.savez(output/'current.npz', time=drive.time, A=drive.values, current=np.asarray(currents))
        if not comparison['passed']:
            raise AssertionError('Si starter reference comparison failed')
    except BaseException as exc:
        report.update(status='failed', error=str(exc))
        raise
    finally:
        report['resources']['wall_seconds'] = time.monotonic()-started
        (output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runtime')
    parser.add_argument('directory')
    args = parser.parse_args()
    print(json.dumps(run(args.runtime, args.directory), indent=2))


if __name__ == '__main__':
    main()
