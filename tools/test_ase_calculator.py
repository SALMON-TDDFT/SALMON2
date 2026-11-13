"""Real fixed-cell ASE optimizer smoke, not a converged geometry optimization.

Usage: python tools/test_ase_calculator.py EXE SI_RPS_DAT NEW_DIRECTORY
"""
from pathlib import Path
import sys,hashlib,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ase.build import bulk
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from salmonpy.ase_calculator import SalmonCalculator
from salmonpy.ase import BOHR_ANGSTROM
from salmonpy.bands import HARTREE_EV
from salmonpy.dft import Pseudopotential,SCFSettings
from salmonpy.fdtd import SalmonExecutable
exe=Path(sys.argv[1]).resolve();out=Path(sys.argv[3]).resolve();out.mkdir(parents=True,exist_ok=False)
atoms=bulk('Si','diamond',a=10.26*BOHR_ANGSTROM,cubic=True);atoms.positions[0,0]+=.12*BOHR_ANGSTROM
atoms.set_constraint(FixAtoms(indices=list(range(1,8))))
calc=SalmonCalculator(executable=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest()),directory=out/'calculations',
     pseudopotentials=[Pseudopotential.from_file(sys.argv[2],atomic_number=14,local_channel=2)],
     electrons=32,states=32,grid=(24,)*3,kmesh=(2,)*3,scf=SCFSettings(max_iterations=500,tolerance=1e-12),threads=2)
atoms.calc=calc
before=atoms.get_potential_energy();first=calc.last_ground_state
forces=atoms.get_forces(apply_constraint=False)
np.testing.assert_array_equal(forces,first.forces_hartree_per_bohr*HARTREE_EV/BOHR_ANGSTROM)
assert atoms.get_potential_energy()==before and calc.last_ground_state is first
fixed=atoms.positions[1:].copy()
optimizer=BFGS(atoms,logfile=str(out/'optimizer.log'),maxstep=.01)
converged=optimizer.run(fmax=.05,steps=1)
after=atoms.get_potential_energy()
assert after<before,'one-step energy did not decrease'
np.testing.assert_array_equal(atoms.positions[1:],fixed)
assert calc.last_ground_state.directory!=first.directory
summary={'energy_before_ev':before,'energy_after_ev':after,'optimizer_steps':optimizer.nsteps,
         'optimizer_converged':bool(converged),'cache_reused':True,'native_force_conversion_exact':True,
         'fixed_atoms_unchanged':True,'fresh_directory_for_changed_geometry':True,
         'binary_sha256':calc._engine.sha256,'grid':[24]*3}
(out/'validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
