"""Explicit GS initial estimates; no solver/MPI initialization or file conversion.

Arrays use the GS/2 conventions. Provider readers return u[k,band,z,y,x]
for one-based inclusive native bounds (x,y,z starts/ends, band, k).
"""
from dataclasses import dataclass
import ctypes as c
import numpy as np
from .errors import InputError,InputTypeError

__all__=['GSInitialGuess']

def _owned(value,dtype,shape):
    array=np.asarray(value,dtype=dtype)
    if array.shape!=shape:raise InputError('initial-guess array shape differs: '+str(shape))
    return np.frombuffer(array.tobytes(),dtype=dtype).reshape(shape)

@dataclass(frozen=True,init=False)
class GSInitialGuess:
    model:object
    occupations:object
    kpoints:object
    weights:object
    _reader:object

    @classmethod
    def from_provider(cls,model,*,read_wavefunctions,occupations,kpoints,weights):
        from .dft import DFT
        if not isinstance(model,DFT):raise InputTypeError('explicit immutable source DFT model required')
        if not callable(read_wavefunctions):raise InputTypeError('read_wavefunctions must be callable')
        nk=np.asarray(kpoints).shape[0] if np.asarray(kpoints).ndim==2 else 0
        if nk<1:raise InputError('nonempty explicit k sampling required')
        obj=object.__new__(cls)
        for name,value in dict(model=model,occupations=_owned(occupations,np.float64,(nk,model.states)),
            kpoints=_owned(kpoints,np.float64,(nk,3)),weights=_owned(weights,np.float64,(nk,)),
            _reader=read_wavefunctions).items():object.__setattr__(obj,name,value)
        return obj

    @classmethod
    def from_arrays(cls,model,*,wavefunctions,occupations,kpoints,weights):
        nk=np.asarray(kpoints).shape[0] if np.asarray(kpoints).ndim==2 else 0
        waves=_owned(wavefunctions,np.complex128,(nk,model.states,*reversed(model.grid)))
        def read(bounds):
            x,y,z,xe,ye,ze,o,oe,k,ke=bounds
            return waves[k-1:ke,o-1:oe,z-1:ze,y-1:ye,x-1:xe]
        return cls.from_provider(model,read_wavefunctions=read,occupations=occupations,kpoints=kpoints,weights=weights)

    @classmethod
    def from_archive(cls,archive,*,model):
        from .gs_archive import GroundStateArchive
        a=GroundStateArchive.open(archive.directory if isinstance(archive,GroundStateArchive) else archive)
        m=a.metadata
        if (m['xc']!=model.xc or m['electrons']!=model.electrons or m['n_bands']!=model.states
                or tuple(reversed(m['grid_zyx']))!=model.grid
                or not np.array_equal(a.array('cell_bohr'),model.cell)
                or not np.array_equal(a.array('atomic_numbers'),model.numbers)
                or not np.array_equal(a.array('fractional_positions'),model.fractional_positions)
                or {p['atomic_number']:p['sha256'] for p in m['pseudopotentials']}!={p.atomic_number:p.sha256 for p in model.pseudopotentials}):
            raise InputError('archive and declared source model differ')
        if model.xc=='tbmbj' and m.get('xc_parameters',{}).get('cval')!=model.cval:
            raise InputError('archive and source XC parameters differ')
        def read(bounds):
            x,y,z,xe,ye,ze,o,oe,k,ke=bounds
            return np.stack([a.wavefunctions(i,o-1,oe)[:,z-1:ze,y-1:ye,x-1:xe] for i in range(k-1,ke)])
        return cls.from_provider(model,read_wavefunctions=read,occupations=a.array('occupations'),
                                 kpoints=a.array('kpoints'),weights=a.array('weights'))

    def read_wavefunctions(self,bounds):
        b=tuple(int(v) for v in bounds)
        if len(b)!=10:raise InputError('initial-guess bounds require10 integers')
        x,y,z,xe,ye,ze,o,oe,k,ke=b
        if any(lo<1 or hi<lo or hi>n for lo,hi,n in zip((x,y,z,o,k),(xe,ye,ze,oe,ke),(*self.model.grid,self.model.states,len(self.weights)))):
            raise InputError('initial-guess tile outside source')
        shape=(ke-k+1,oe-o+1,ze-z+1,ye-y+1,xe-x+1)
        return _owned(self._reader(b),np.complex128,shape)

    def _callback(self):
        # Bound only for the synchronous native call, never kept by Fortran.
        prototype=c.CFUNCTYPE(c.c_int,c.c_void_p,*([c.c_void_p]*6))
        failures=[]
        def fill(cookie,bounds,shape,waves,occupations,points,weights):
            try:
                b=tuple(np.ctypeslib.as_array(c.cast(bounds,c.POINTER(c.c_int64)),shape=(10,)))
                dims=tuple(np.ctypeslib.as_array(c.cast(shape,c.POINTER(c.c_int64)),shape=(5,)))
                if dims[3:]!=(self.model.states,len(self.weights)):raise InputError('source band/k counts differ')
                values=self.read_wavefunctions(b)
                def put(pointer,data):
                    flat=np.ascontiguousarray(data).view(np.float64).ravel()
                    np.ctypeslib.as_array(c.cast(pointer,c.POINTER(c.c_double)),shape=(flat.size,))[:]=flat
                put(waves,values);put(occupations,self.occupations);put(points,self.kpoints);put(weights,self.weights)
                return 0
            except Exception as exc:
                failures.append(str(exc));return 11
        return prototype(fill),failures
