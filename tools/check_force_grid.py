"""Displaced-Si force and egg-box diagnostic; no finite difference is assumed exact.

Usage: python tools/check_force_grid.py EXE SI_RPS_DAT NEW_DIRECTORY --grid 24
The full suite runs independent GS calculations and can take several minutes.
"""
from pathlib import Path
import argparse,hashlib,json,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ase.build import bulk
from salmonpy.ase import dft_from_atoms,BOHR_ANGSTROM
from salmonpy.dft import Pseudopotential,SCFSettings
from salmonpy.fdtd import SalmonExecutable

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('executable');p.add_argument('pseudo');p.add_argument('directory')
p.add_argument('--grid',type=int,required=True)
p.add_argument('--steps',type=float,nargs='+',default=[.02,.01])
p.add_argument('--translations',type=float,nargs='*',default=None,help='rigid x translations in bohr; default half/full grid')
a=p.parse_args();out=Path(a.directory).resolve();out.mkdir(parents=True,exist_ok=False)
if any(not np.isfinite(h) or h<=0 for h in a.steps):p.error('steps must be positive and finite')
exe=Path(a.executable).resolve();engine=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest())
pseudo=Pseudopotential.from_file(a.pseudo,atomic_number=14,local_channel=2)
atoms=bulk('Si','diamond',a=10.26*BOHR_ANGSTROM,cubic=True);atoms.positions[0,0]+=.12*BOHR_ANGSTROM

def calculate(atoms,name):
    model=dft_from_atoms(atoms,pseudopotentials=[pseudo],electrons=32,states=32,grid=(a.grid,)*3,kmesh=(2,)*3)
    return model.run(engine,out/name,scf=SCFSettings(max_iterations=500,tolerance=1e-12),threads=2,timeout=900.)

base=calculate(atoms,'base');forces=[];translations=[]
record={'grid':a.grid,'spacing_bohr':10.26/a.grid,'binary_sha256':engine.sha256,'base_energy_hartree':base.energy_hartree,
        'base_force_ha_per_bohr':base.forces_hartree_per_bohr.tolist(),'finite_differences':forces,'rigid_translations':translations}

def save(): (out/'diagnostic.json').write_text(json.dumps(record,indent=2))
save()
for i,h in enumerate(a.steps):
    energies=[]
    for sign in (1,-1):
        moved=atoms.copy();moved.positions[0,0]+=sign*h*BOHR_ANGSTROM
        energies.append(calculate(moved,f'difference-{i}-{sign}').energy_hartree)
    fd=-(energies[0]-energies[1])/(2*h)
    forces.append({'step_bohr':h,'native':float(base.forces_hartree_per_bohr[0,0]),'energy_derivative':fd})
    save();print(forces[-1],flush=True)
shifts=[10.26/a.grid/2,10.26/a.grid] if a.translations is None else a.translations
for i,shift in enumerate(shifts):
    if not np.isfinite(shift):p.error('translations must be finite')
    moved=atoms.copy();moved.positions[:,0]+=shift*BOHR_ANGSTROM
    result=calculate(moved,f'translation-{i}')
    translations.append({'shift_bohr':shift,'energy_hartree':result.energy_hartree,
                         'total_force_ha_per_bohr':result.forces_hartree_per_bohr.sum(axis=0).tolist()})
    save();print(translations[-1],flush=True)
