"""Common A-driven orchestration for user-provided crystalline backends.

The built-in velocity-gauge SBE adapter is salmonpy.sbe.SBEInitialState.
The TDDFT adapter is not yet implemented. A is SALMON Ac=A/c;
J is macroscopic current density in common laboratory coordinates, atomic units.
"""
from dataclasses import dataclass
from typing import Protocol
import numpy as np
from .fields import FieldSamples, _readonly
from .ld import _array


class AResponseState(Protocol):
    """Backend-owned mutable state; dt/model fixed, step failures must be atomic."""
    @property
    def time(self) -> float: ...
    @property
    def current(self) -> np.ndarray: ...
    def step(self, A: np.ndarray) -> np.ndarray:
        """Advance using the new endpoint Ac; return physical J at the new time."""
        ...


class AInitialState(Protocol):
    """Reusable read-only initial condition and model, e.g. TDDFT or SBE.

    new_state must allocate independent evolution history. No run may mutate the
    initial condition or another state. Loading and validating engine-specific
    initial files is the backend's responsibility; these are not checkpoints.
    """
    @property
    def backend(self) -> str: ...
    def new_state(self, *, dt: float, time0: float, initial_A: np.ndarray) -> AResponseState: ...


@dataclass(frozen=True,slots=True)
class CrystalResponse:
    time: np.ndarray
    current: np.ndarray
    backend: str
    dt: float
    field_convention: str='Ac=A/c; E=-dAc/dt; atomic units; laboratory frame'


def respond(initial: AInitialState, drive: FieldSamples) -> CrystalResponse:
    """Apply an entire A waveform using the backend's same single-step method.

    Each call creates an independent state. E-driven crystalline inputs are
    rejected. The first output is the backend's initial current, not forced zero.
    """
    if not isinstance(drive,FieldSamples) or drive.kind!='A':
        raise ValueError('crystalline response requires FieldSamples(kind="A")')
    dt=np.diff(drive.time)
    if not np.allclose(dt,dt[0],rtol=1e-8,atol=1e-12*abs(dt[0])):
        raise ValueError('uniform A samples are required')
    backend=initial.backend
    if not isinstance(backend,str) or not backend: raise ValueError('backend identifier must be nonempty')
    state=initial.new_state(dt=float(dt[0]),time0=float(drive.time[0]),initial_A=drive.values[0].copy())
    current=np.empty_like(drive.values)
    for i,t in enumerate(drive.time):
        if i==0: current[i]=_array(state.current,'initial current',(3,))
        else: current[i]=_array(state.step(drive.values[i].copy()),'current',(3,))
        # Permit roundoff from repeated time += dt, bounded far below one step.
        scale=max(abs(float(t)),abs(float(drive.time[0])),abs(float(dt[0])))
        roundoff=8*np.finfo(float).eps*(i+1)*scale
        tolerance=min(abs(float(dt[0]))*1e-3,max(abs(float(dt[0]))*1e-8,roundoff))
        if not np.isfinite(state.time) or not np.isclose(state.time,t,rtol=0,atol=tolerance):
            raise ValueError('backend time does not match the common A grid')
    return CrystalResponse(drive.time,_readonly(current),backend,float(dt[0]))
