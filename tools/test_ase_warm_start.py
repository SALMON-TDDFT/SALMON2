"""Compare opt-in ASE warm start with fresh SCF and native forces for Si.

python tools/test_ase_warm_start.py EXE SI_KY SI_PBE_UPF NEW_OUTPUT [--xc pz ...]
Requires an audited executable with Libxc for the default three-XC matrix.
One thread; 16^3 grid is a connection test, not force/grid convergence.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from ase.build import bulk
from salmonpy.ase import BOHR_ANGSTROM, dft_from_atoms
from salmonpy.ase_calculator import SalmonCalculator
from salmonpy.bands import HARTREE_EV
from salmonpy.dft import Pseudopotential, SCFSettings, GroundStateResult
from salmonpy.fdtd import SalmonExecutable


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable')
    parser.add_argument('ky')
    parser.add_argument('upf')
    parser.add_argument('output')
    parser.add_argument('--xc', nargs='+', default=['pz', 'libxc_pz', 'libxc_pbe'])
    args = parser.parse_args()
    root = Path(args.output).resolve(); root.mkdir(parents=True, exist_ok=False)
    exe = Path(args.executable).resolve()
    engine = SalmonExecutable(exe, hashlib.sha256(exe.read_bytes()).hexdigest())
    os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    summary = {'executable_sha256': engine.sha256, 'threads': 1, 'cases': {},
               'scope': 'fresh/warm numerical regression, not force convergence'}
    for xc in args.xc:
        started = time.monotonic()
        atoms = bulk('Si', 'diamond', a=10.26*BOHR_ANGSTROM, cubic=True)
        pseudo = Pseudopotential.from_file(args.upf if xc=='libxc_pbe' else args.ky,
                                           atomic_number=14, local_channel=2)
        settings = dict(pseudopotentials=[pseudo], states=32, grid=(16,)*3, kmesh=(2,)*3, xc=xc)
        scf = SCFSettings(max_iterations=500, tolerance=1e-11)
        calc = SalmonCalculator(executable=engine, directory=root/xc/'warm',
                                 scf=scf, threads=1, warm_start=True, **settings)
        atoms.calc = calc
        atoms.get_potential_energy()
        first = calc.last_ground_state
        atoms.positions[0,0] += .05*BOHR_ANGSTROM
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        warm = calc.last_ground_state
        cold = dft_from_atoms(atoms, **settings).run(engine, root/xc/'cold', scf=scf, threads=1)
        energy_error = abs(warm.energy_hartree-cold.energy_hartree)
        force_error = float(np.max(np.abs(warm.forces_hartree_per_bohr-cold.forces_hartree_per_bohr)))
        np.testing.assert_allclose(warm.energy_hartree,cold.energy_hartree,rtol=0,atol=1e-7)
        np.testing.assert_allclose(warm.forces_hartree_per_bohr,cold.forces_hartree_per_bohr,rtol=0,atol=2e-6)
        np.testing.assert_array_equal(forces,warm.forces_hartree_per_bohr*HARTREE_EV/BOHR_ANGSTROM)
        assert energy == warm.energy_hartree*HARTREE_EV
        record=json.loads((warm.directory/'ground-state.json').read_text())
        np.testing.assert_allclose(record['fractional_positions'],atoms.get_scaled_positions(wrap=False),rtol=0,atol=1e-14)
        assert record['initial_guess']['mode']=='wavefunction'
        # Native checkpoint atom coordinates must be the new geometry, not the old GS.
        coordinates=warm.directory/'native/data_for_restart/atomic_coor.txt'
        actual=np.loadtxt(coordinates, usecols=(1,2,3))*BOHR_ANGSTROM
        delta=(actual-atoms.positions) @ np.linalg.inv(atoms.cell.array)
        delta-=np.rint(delta)
        np.testing.assert_allclose(delta,0,rtol=0,atol=1e-10)
        checkpoint_files=[p.name for p in (warm.directory/'native/data_for_restart').iterdir()]
        GroundStateResult.open(first.directory)
        summary['cases'][xc] = {'fresh_iterations':cold.iterations,'warm_iterations':warm.iterations,
            'energy_difference_hartree':energy_error,'force_max_difference_au':force_error,
            'native_force_conversion_exact':True,'source_unchanged':True,
            'checkpoint_files':checkpoint_files,'pseudo_sha256':pseudo.sha256,
            'seconds':time.monotonic()-started}
        (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
