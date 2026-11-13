"""Optional ASE Calculator over the independent DFT runner.

This module explicitly requires ASE. It is not imported by the core package or
by salmonpy.ase. Model settings are fixed at construction; geometry is ASE-owned.
"""
from pathlib import Path
from types import MappingProxyType
from uuid import uuid4
import numpy as np
try:
    from ase.calculators.calculator import Calculator, PropertyNotImplementedError, all_changes
except ImportError as exc:
    raise ImportError('SalmonCalculator requires the optional ASE dependency') from exc
from .dft import Pseudopotential, SCFSettings
from .ase import dft_from_atoms, BOHR_ANGSTROM
from .bands import HARTREE_EV
from .fdtd import SalmonExecutable
from .sbe import _integer
from .ld import _positive


class SalmonCalculator(Calculator):
    """Energy (eV) and native forces (eV/Å), fixed ions per SCF, no stress.

    Each changed geometry gets a fresh, retained GS subdirectory. The initial
    supported force scope is PZ with scalar KY pseudopotentials. Native forces
    are not replaced by derivatives of a possibly egg-box-contaminated energy.
    Users must establish grid and SCF convergence for their optimization target.
    """
    implemented_properties=['energy','forces']

    def __init__(self,*,executable,directory,pseudopotentials,electrons,states,grid,kmesh,
                 xc='pz',scf=SCFSettings(tolerance=1e-10),threads=1,timeout=600.,atoms=None):
        if not isinstance(executable,SalmonExecutable): raise TypeError('expected SalmonExecutable')
        pseudos=tuple(pseudopotentials)
        if not pseudos or any(not isinstance(p,Pseudopotential) for p in pseudos): raise TypeError('expected pseudopotentials')
        if xc.lower()!='pz' or any(p.suffix!='_rps.dat' for p in pseudos):
            raise ValueError('initial Calculator force scope is PZ with scalar KY potentials')
        if not isinstance(scf,SCFSettings): raise TypeError('expected SCFSettings')
        grid=tuple(_integer(v,'grid',8,2**20) for v in grid)
        kmesh=tuple(_integer(v,'kmesh',1,2**20) for v in kmesh)
        if len(grid)!=3 or len(kmesh)!=3: raise ValueError('three grid/kmesh dimensions required')
        electrons=_integer(electrons,'electrons',2,2**31-1);states=_integer(states,'states',1,2**31-1)
        if electrons%2 or states<electrons//2: raise ValueError('closed shell and enough states required')
        self._settings=MappingProxyType(dict(pseudopotentials=pseudos,electrons=electrons,states=states,
                                             grid=grid,kmesh=kmesh,xc='pz'))
        self._engine=executable;self._scf=scf
        self._threads=_integer(threads,'threads',1,2**31-1);self._timeout=_positive(timeout,'timeout')
        self._last_ground_state=None
        super().__init__(directory=str(Path(directory).expanduser().resolve()),atoms=atoms)

    @property
    def settings(self): return self._settings

    @property
    def last_ground_state(self): return self._last_ground_state

    def set(self,**kwargs):
        if kwargs: raise ValueError('Calculator settings are immutable; construct a new calculator')
        return {}

    def calculate(self,atoms=None,properties=('energy',),system_changes=all_changes):
        if any(p not in self.implemented_properties for p in properties):
            raise PropertyNotImplementedError('only energy and forces are implemented')
        self.results={};self._last_ground_state=None
        super().calculate(atoms,properties,system_changes)
        model=dft_from_atoms(self.atoms,**self._settings)
        run_directory=Path(self.directory)/('gs-'+uuid4().hex)
        result=model.run(self._engine,run_directory,scf=self._scf,threads=self._threads,timeout=self._timeout)
        energy=float(result.energy_hartree*HARTREE_EV)
        forces=np.array(result.forces_hartree_per_bohr,copy=True)*HARTREE_EV/BOHR_ANGSTROM
        if forces.shape!=(len(self.atoms),3) or not np.isfinite(forces).all() or not np.isfinite(energy):
            raise RuntimeError('invalid native energy/force result')
        self.results={'energy':energy,'forces':forces}
        self._last_ground_state=result
