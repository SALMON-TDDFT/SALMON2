"""SALMON velocity-gauge SBE. Atomic units; driving A denotes A/c.

States return the electric (charge) current density j. The native kernel's
SALMON electronic matter current J is kept internally (``matter_current``); j = -J.

Initial arrays use Fortran index order: (band, k), (band, band, xyz, k).
Closed-shell occupations match SALMON: two electrons in each of ne/2 bands.
"""
from .errors import InputError, NumericalError, RuntimeFailure, SnapshotError
from .snapshots import _SnapshotMixin, _fingerprint
import ctypes as ct
from dataclasses import dataclass, field
from functools import lru_cache
import hashlib
from pathlib import Path
import sys
import re
import numpy as np
from .ld import _array, _positive


def _frozen(value, dtype=float):
    a = np.array(value, dtype=dtype, order='F')
    if not np.isfinite(a).all():
        raise InputError('initial data must be finite')
    return np.frombuffer(a.tobytes(order='F'), dtype=a.dtype).reshape(a.shape, order='F')


def _integer(value, name, low, high):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or not low <= value <= high:
        raise InputError(f'{name} must be an integer in [{low}, {high}]')
    return int(value)


class _WeakInitial:
    __slots__=('__weakref__',)


@dataclass(frozen=True, slots=True, eq=False)
class SBEInitialState(_WeakInitial):
    eigen: np.ndarray
    momentum: np.ndarray
    nonlocal_velocity: np.ndarray
    weights: np.ndarray
    lattice: np.ndarray
    electrons: int
    correction_order: int = 0
    nonlocal_correction: bool = False
    backend: str = field(default='SALMON-SBE-velocity', init=False)
    volume: float = field(init=False)

    def __post_init__(self):
        for name in ('eigen','weights','lattice'):
            object.__setattr__(self,name,_frozen(_array(getattr(self,name),name)))
        for name in ('momentum','nonlocal_velocity'):
            object.__setattr__(self,name,_frozen(getattr(self,name),complex))
        if self.eigen.ndim != 2 or min(self.eigen.shape)<1:
            raise InputError('eigen must have shape (bands, kpoints)')
        nb,nk=self.eigen.shape
        if self.momentum.shape != (nb,nb,3,nk) or self.nonlocal_velocity.shape != (nb,nb,3,nk):
            raise InputError('matrix shape must be (bands, bands, 3, kpoints)')
        if self.weights.shape != (nk,) or np.any(self.weights<0) or not np.isfinite(self.weights.sum()) or self.weights.sum()<=0:
            raise InputError('weights must be nonnegative with positive finite sum')
        if self.lattice.shape != (3,3):
            raise InputError('lattice must contain three column vectors')
        volume=float(np.dot(np.cross(self.lattice[:,0],self.lattice[:,1]),self.lattice[:,2]))
        object.__setattr__(self,'volume',_positive(volume,'right-handed cell volume'))
        ne=_integer(self.electrons,'electrons',0,2*nb)
        if ne%2: raise InputError('SALMON SBE requires an even closed-shell electron count')
        object.__setattr__(self,'electrons',ne)
        object.__setattr__(self,'correction_order',_integer(self.correction_order,'correction_order',0,3))
        if not isinstance(self.nonlocal_correction,(bool,np.bool_)):
            raise InputError('nonlocal_correction must be boolean')

    def __reduce__(self):
        # Process-local handles are never serialized; construction refreezes data.
        return (type(self),(self.eigen,self.momentum,self.nonlocal_velocity,self.weights,
                            self.lattice,self.electrons,self.correction_order,self.nonlocal_correction))

    def new_state(self, *, dt, time0=0., initial_A=(0.,0.,0.)):
        return SBEState(self,dt,time0,initial_A)

    @classmethod
    def from_salmon(cls, directory, sysname, *, bands, kpoints, electrons, lattice,
                    correction_order=0, nonlocal_correction=False, active_bands=None):
        """Read single-spin atomic-unit *_k/eigen/tm.data exports.

        Occupation columns are ignored, exactly as in SALMON SBE. Binary restart,
        spin-polarized exports and eV eigenvalues are not supported.
        """
        nb=_integer(bands,'bands',1,2**31-1); nk=_integer(kpoints,'kpoints',1,2**31-1)
        if not isinstance(sysname,str) or not sysname or Path(sysname).name!=sysname:
            raise InputError('sysname must be a filename stem')
        root=Path(directory)
        def lines(suffix): return iter((root/(sysname+suffix)).read_text().splitlines())
        def skip(it,n):
            for _ in range(n): next(it)
        def row(it,indices,count):
            fields=next(it).replace('D','E').replace('d','e').split()
            if len(fields)!=len(indices)+count or [int(x) for x in fields[:len(indices)]]!=list(indices):
                raise InputError('SALMON export indices or column count mismatch')
            values=np.array(fields[len(indices):],float)
            if not np.isfinite(values).all(): raise InputError('nonfinite export value')
            return values
        def end(it):
            if any(s.strip() for s in it): raise InputError('unexpected trailing data (multiple spins are unsupported)')
        try:
            it=lines('_k.data'); first=next(it)
            if 'nonorthogonal' in first:
                raise InputError('nonorthogonal k exports are not supported by the legacy SBE text loader')
            skip(it,4)
            weights=np.array([row(it,(ik+1,),4)[3] for ik in range(nk)])
            tail=[line for line in it if line.strip()]
            if tail:
                if len(tail)!=2 or tail[0].strip()!='# coefficients (2*pi/a [a.u.]) in kx, ky, kz':
                    raise InputError('unexpected k export footer')
                _array([float(x.replace('D','E').replace('d','e')) for x in tail[1].split()], 'reciprocal coefficients',(3,))
            it=lines('_eigen.data'); header=[next(it) for _ in range(3)]
            if 'esp[a.u.]' not in header[2]: raise InputError('eigen export must explicitly declare esp[a.u.]')
            eigen=np.empty((nb,nk))
            for ik in range(nk):
                block=next(it)
                match=re.fullmatch(r'\s*k=\s*(\d+)\s*,\s*spin=\s*(\d+)\s*',block)
                if not match or tuple(map(int,match.groups()))!=(ik+1,1):
                    raise InputError('expected sequential k blocks with spin=1')
                for ib in range(nb): eigen[ib,ik]=row(it,(ib+1,),2)[0]
            end(it)
            it=lines('_tm.data'); skip(it,3); matrices=[]
            for part in range(2):
                if part: next(it)
                p=np.empty((nb,nb,3,nk),complex)
                for ik in range(nk):
                    for ib in range(nb):
                        for jb in range(nb):
                            v=row(it,(ik+1,ib+1,jb+1),6)
                            p[ib,jb,:,ik]=v[::2]+1j*v[1::2]
                matrices.append(p)
            end(it)
        except (StopIteration,IndexError) as exc:
            raise InputError('truncated SALMON SBE export') from exc
        active=nb if active_bands is None else _integer(active_bands,'active_bands',1,nb)
        return cls(eigen[:active],matrices[0][:active,:active],matrices[1][:active,:active],
                   weights,lattice,electrons,correction_order,nonlocal_correction)


@lru_cache(maxsize=1)
def _native():
    suffix={'darwin':'.dylib','win32':'.dll'}.get(sys.platform,'.so')
    path=Path(__file__).parent/'_native'/('libsalmonpy_ld'+suffix)
    if not path.exists(): raise RuntimeFailure('Build native kernels: python -m salmonpy.build_native')
    digest=path.with_suffix(path.suffix+'.sha256')
    if not digest.exists() or digest.read_text().strip()!=hashlib.sha256(path.read_bytes()).hexdigest():
        raise RuntimeFailure('Native library checksum mismatch; rebuild')
    lib=ct.CDLL(str(path)); fn=lib.salmonpy_sbe
    fn.argtypes=[ct.c_int]*5+[ct.c_double]*2+[ct.c_void_p]*7+[ct.c_int,ct.c_void_p]
    fn.restype=None
    return fn


def _handles_api():
    from ._sbe_handles import default_api
    return default_api()


def _model_handle(initial):
    from ._sbe_handles import model_for
    return model_for(initial,_handles_api())


class SBEState(_SnapshotMixin):
    """Independent density matrix; fixed dt/model; step accepts only the next A/c."""
    __slots__=('_initial','_dt','_time0','_steps','_A','_current','_native_state','_closed')
    _snapshot_kind = 'sbe'
    _snapshot_counter = '_steps'

    def _snapshot_owner(self): return self._initial
    def _snapshot_settings(self): return _fingerprint((self.dt, self._time0))
    def _snapshot_arrays(self):
        return {'_A':self._A,'_rho':self._native_state.density(),'_current':self._current}

    def __init__(self,initial,dt,time0,initial_A):
        self._initial=initial; self._dt=_positive(dt,'dt')
        self._time0=float(_array(time0,'time0',()))
        self._steps=0; self._A=_array(initial_A,'initial_A',(3,))
        self._closed=False
        self._native_state=_model_handle(initial).new_state(self.dt,self._A)
        self._current=self._native_state.current.copy()
    @property
    def dt(self): return self._dt
    @property
    def time(self): return self._time0+self._steps*self._dt
    @property
    def timing(self):
        """salmonpy.time/1 declaration (crystal contract; Taylor step at midpoint A)."""
        from .timing import TimeSpec,InputPort,OutputPort
        return TimeSpec(self._dt,(InputPort('A'),),(OutputPort('j'),),t0=self._time0,
            integrator='4th-order Taylor with midpoint A; J evaluated at the new endpoint A')
    @property
    def density(self): return _frozen(self._native_state.density(),complex)
    @property
    def current(self):
        """Electric (charge) current density j = -J_matter at the current time."""
        return _frozen(-self._current)
    @property
    def matter_current(self):
        """SALMON electronic matter current J (internal representation)."""
        return _frozen(self._current)
    def step(self,Anew):
        if self._closed:raise RuntimeFailure('SBE state is closed')
        A=_array(Anew,'Anew',(3,))
        if not np.isfinite(self._time0+(self._steps+1)*self.dt):
            raise InputError('time overflow')
        current=self._native_state.step(A)
        self._A=A;self._current=current;self._steps+=1
        return self.current

    def snapshot(self):
        if self._closed:raise SnapshotError('SBE state is closed')
        return super().snapshot()

    def restore(self,snapshot):
        if self._closed:raise SnapshotError('SBE state is closed')
        from types import SimpleNamespace
        nb,nk=self._native_state.model.shape
        target=SimpleNamespace(dt=self.dt,_time0=self._time0,_snapshot_kind=self._snapshot_kind,
            _snapshot_counter=self._snapshot_counter,_snapshot_owner=self._snapshot_owner,
            _snapshot_settings=self._snapshot_settings,
            _snapshot_arrays=lambda:{'_A':np.empty(3),'_rho':np.empty((nb,nb,nk),complex,order='F'),
                                     '_current':np.empty(3)})
        _SnapshotMixin.restore(target,snapshot)
        self._native_state.replace(target._rho,target._A,target._current)
        self._A=target._A;self._current=target._current;self._steps=target._steps

    def close(self):
        if not self._closed:
            self._native_state.close();self._closed=True

    def __enter__(self):return self
    def __exit__(self,*args):self.close()
