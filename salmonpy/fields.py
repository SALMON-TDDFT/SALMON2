"""Pulse tools in atomic units; A denotes SALMON Ac=A/c, not unscaled A."""
from dataclasses import dataclass
import numpy as np
from .ld import _array, _positive


def _readonly(a):
    a=np.ascontiguousarray(a,dtype=float)
    return np.frombuffer(a.tobytes(),dtype=float).reshape(a.shape)


def _times(time):
    time=_array(time,'time')
    if time.ndim != 1 or len(time)<2 or np.any(np.diff(time)<=0):
        raise ValueError('time must be a strictly increasing vector with at least two samples')
    return time


@dataclass(frozen=True,slots=True)
class FieldSamples:
    time: np.ndarray
    values: np.ndarray
    kind: str

    def __post_init__(self):
        t=_times(self.time);v=_array(self.values,'field',(len(t),3))
        if self.kind not in ('A','E'): raise ValueError('kind must be A or E')
        object.__setattr__(self,'time',_readonly(t));object.__setattr__(self,'values',_readonly(v))

    def electric_field(self):
        """Backward difference at time[1:], matching SALMON's discrete E=-dAc/dt."""
        if self.kind!='A': raise ValueError('electric_field requires A samples')
        if len(self.time)<3: raise ValueError('at least three A samples are required')
        return FieldSamples(self.time[1:],-np.diff(self.values,axis=0)/np.diff(self.time)[:,None],'E')

    def write_salmon(self,path):
        """Write four columns for SALMON ae_shape='read', unit_system='au'."""
        if self.kind!='A': raise ValueError('SALMON input file requires A samples')
        # Existing files are preserved; explicit choice of a new path is required.
        with open(path,'x',encoding='utf-8') as f:
            np.savetxt(f,np.column_stack((self.time,self.values)),
                       header='1:Time[a.u.] 2:Ac_ext_x[a.u.] 3:Ac_ext_y[a.u.] 4:Ac_ext_z[a.u.]')


@dataclass(frozen=True,slots=True)
class AcosPulse:
    """SALMON Acos2 with linear polarization; CEP is in cycles, not radians.

    amplitude is SALMON E_amplitude, an envelope scale (not exact peak E).
    """
    amplitude: float
    omega: float
    duration: float
    start: float=0.
    polarization: tuple=(0.,0.,1.)
    cep: float=0.

    def __post_init__(self):
        for key in ('amplitude','omega','duration'):
            object.__setattr__(self,key,_positive(getattr(self,key),key,zero=(key=='amplitude')))
        for key in ('start','cep'):
            raw=np.asarray(getattr(self,key))
            if raw.shape!=() or raw.dtype.kind not in 'iuf': raise ValueError(f'{key} must be a real numeric scalar')
            a=_array(getattr(self,key),key,())
            object.__setattr__(self,key,float(a))
        p=_array(self.polarization,'polarization',(3,))
        if not np.isclose(np.linalg.norm(p),1.,rtol=0,atol=1e-12):
            raise ValueError('polarization must be a unit vector')
        object.__setattr__(self,'polarization',tuple(map(float,p)))

    def sample(self,time):
        t=_times(time);tt=t-self.start-.5*self.duration
        a=np.zeros_like(t);inside=np.abs(tt)<.5*self.duration
        a[inside]=-self.amplitude/self.omega*np.cos(np.pi*tt[inside]/self.duration)**2*np.sin(
            self.omega*tt[inside]+2*np.pi*self.cep)
        return FieldSamples(t,a[:,None]*np.array(self.polarization),'A')


def sample_pulses(pulses,time):
    """Compose independently delayed/polarized Acos2 pulses on a common time grid."""
    t=_times(time);values=np.zeros((len(t),3))
    for pulse in pulses:
        if not isinstance(pulse,AcosPulse): raise TypeError('expected AcosPulse')
        values+=pulse.sample(t).values
    return FieldSamples(t,values,'A')
