"""Optional ASE Calculator over the independent DFT runner.

This module explicitly requires ASE. It is not imported by the core package or
by salmonpy.ase. Model settings are fixed at construction; geometry is ASE-owned.
"""
from .errors import InputError, InputTypeError, RuntimeFailure
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


def _has_nlcc(pseudo):
    """Conservatively recognize NLCC in validated NC-UPF bytes.

    SALMON reads PP_NLCC data independently of the header flag. Unknown flags
    and even zero-valued PP_NLCC blocks remain guarded until force validation.
    """
    if pseudo.suffix != '.upf': return False
    import xml.etree.ElementTree as ET
    root = ET.fromstring(pseudo.data)
    flag = root.find('PP_HEADER').get('core_correction', '').strip().lower()
    return flag not in ('', 'f', 'false', '.false.', '0') or root.find('PP_NLCC') is not None


class SalmonCalculator(Calculator):
    """Energy (eV) and native forces (eV/Å), fixed ions per SCF, no stress.

    Each changed geometry gets a fresh, retained GS subdirectory. Supported XC is
    PZ/Libxc-PZ/Libxc-PBE with scalar KY/NC-UPF pseudopotentials. Opt-in
    warm_start reuses compatible wavefunctions, not mixing history. Native forces
    are not replaced by derivatives of a possibly egg-box-contaminated energy.
    NLCC UPF forces require matching completed runtime implementation validation.
    Users must establish grid and SCF convergence for their optimization target.
    """
    implemented_properties=['energy','forces']

    def __init__(self,*,executable,directory,pseudopotentials,states,grid,kmesh,electrons=None,
                 xc='pz',scf=SCFSettings(tolerance=1e-10),threads=1,timeout=600.,atoms=None,
                 warm_start=False):
        if not isinstance(executable,SalmonExecutable): raise InputTypeError('expected SalmonExecutable')
        pseudos=tuple(pseudopotentials)
        if not pseudos or any(not isinstance(p,Pseudopotential) for p in pseudos): raise InputTypeError('expected pseudopotentials')
        if xc.lower() not in ('pz', 'libxc_pz', 'libxc_pbe'):
            raise InputError('Calculator supports PZ, Libxc-PZ and Libxc-PBE with scalar KY/NC-UPF potentials')
        if not isinstance(warm_start, bool):
            raise InputTypeError('warm_start must be a bool')
        self._nlcc = any(_has_nlcc(p) for p in pseudos)
        self._nlcc_forces = xc.lower() in executable.nlcc_force_xc
        self.implemented_properties = ['energy'] if self._nlcc and not self._nlcc_forces else ['energy', 'forces']
        self._warm_start = warm_start
        if not isinstance(scf,SCFSettings): raise InputTypeError('expected SCFSettings')
        grid=tuple(_integer(v,'grid',8,2**20) for v in grid)
        kmesh=tuple(_integer(v,'kmesh',1,2**20) for v in kmesh)
        if len(grid)!=3 or len(kmesh)!=3: raise InputError('three grid/kmesh dimensions required')
        states=_integer(states,'states',1,2**31-1)
        if electrons is not None:
            electrons=_integer(electrons,'electrons',2,2**31-1)
            if electrons%2 or states<electrons//2: raise InputError('closed shell and enough states required')
        self._settings=MappingProxyType(dict(pseudopotentials=pseudos,electrons=electrons,states=states,
                                             grid=grid,kmesh=kmesh,xc=xc.lower()))
        self._engine=executable;self._scf=scf
        self._threads=_integer(threads,'threads',1,2**31-1);self._timeout=_positive(timeout,'timeout')
        self._last_ground_state=None
        self._last_model=None
        super().__init__(directory=str(Path(directory).expanduser().resolve()),atoms=atoms)

    @property
    def settings(self): return self._settings

    @property
    def last_ground_state(self): return self._last_ground_state

    @property
    def warm_start(self): return self._warm_start

    def set(self,**kwargs):
        if kwargs: raise InputError('Calculator settings are immutable; construct a new calculator')
        return {}

    def _check_nlcc_forces(self, properties):
        if self._nlcc and not self._nlcc_forces and 'forces' in properties:
            raise PropertyNotImplementedError(
                'NLCC UPF forces require matching completed runtime validation; '
                'use matching PZ + KY potentials for structural optimization')

    def get_property(self, name, atoms=None, allow_calculation=True):
        self._check_nlcc_forces((name,))
        return super().get_property(name, atoms, allow_calculation)

    def calculate(self,atoms=None,properties=('energy',),system_changes=all_changes):
        self._check_nlcc_forces(properties)
        if any(p not in self.implemented_properties for p in properties):
            raise PropertyNotImplementedError('only energy and forces are implemented')
        previous = self._last_ground_state
        previous_model = self._last_model
        self.results={};self._last_ground_state=None;self._last_model=None
        super().calculate(atoms,properties,system_changes)
        model=dft_from_atoms(self.atoms,**self._settings)
        run_directory=Path(self.directory)/('gs-'+uuid4().hex)
        options = {}
        if (self._warm_start and previous is not None and previous_model is not None
                and model._restart_signature() == previous_model._restart_signature()):
            options['initial_guess'] = previous
        result=model.run(self._engine,run_directory,scf=self._scf,threads=self._threads,timeout=self._timeout,**options)
        energy=float(result.energy_hartree*HARTREE_EV)
        if not np.isfinite(energy): raise RuntimeFailure('invalid native energy result')
        values={'energy':energy}
        if not self._nlcc or self._nlcc_forces:
            forces=np.array(result.forces_hartree_per_bohr,copy=True)*HARTREE_EV/BOHR_ANGSTROM
            if forces.shape!=(len(self.atoms),3) or not np.isfinite(forces).all():
                raise RuntimeFailure('invalid native force result')
            values['forces']=forces
        self.results=values
        self._last_ground_state=result
        self._last_model=model
