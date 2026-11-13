"""SALMON table reading and explicitly normalized Fourier transforms."""
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
import re
import numpy as np
from .ld import _array
from .fields import _readonly, _times


@dataclass(frozen=True,slots=True)
class Table:
    columns: tuple
    units: tuple
    values: np.ndarray
    headers: tuple
    source: str

    def column(self,name):
        return self.values[:,self.columns.index(name)]


def read_table(path):
    """Read numbered SALMON headers, Fortran D exponents and comment lines.

    Units are retained, never inferred or converted. Unlabeled columns are named
    column_1, column_2,... with unknown units represented by None.
    """
    path=Path(path);headers=[];rows=[];labels={}
    for line in path.read_text().splitlines():
        line=line.strip()
        if not line: continue
        if line.startswith(('#','!')):
            headers.append(line)
            for n,label,unit in re.findall(r'(\d+):([^\s\[]+)\[([^\]]+)\]',line):
                labels[int(n)-1]=(label,unit)
            continue
        clean=re.split(r'[#!]',line,maxsplit=1)[0]
        rows.append([float(x.replace('D','E').replace('d','e')) for x in clean.split()])
    if not rows: raise ValueError('table contains no samples')
    data=_array(rows,'table')
    if data.ndim!=2 or data.shape[1]==0: raise ValueError('table must be rectangular')
    names=tuple(labels.get(i,(f'column_{i+1}',None))[0] for i in range(data.shape[1]))
    if len(set(names))!=len(names): raise ValueError('column names must be unique')
    units=tuple(labels.get(i,('',None))[1] for i in range(data.shape[1]))
    return Table(names,units,_readonly(data),tuple(headers),str(path.resolve()))


@dataclass(frozen=True,slots=True)
class Spectrum:
    omega: np.ndarray
    values: np.ndarray
    config: object


def spectrum(time,values,*,window='hann',time_unit='a.u.',value_unit='a.u.'):
    """dt * sum(w[n]*x[n]*exp(+i*omega*t[n])), positive frequencies only.

    No detrending, zero-padding, one-sided doubling or window-gain correction.
    This is a signal spectrum, not an automatically normalized optical observable.
    """
    t=_times(time);x=_array(values,'values');dt=np.diff(t)
    if x.ndim not in (1,2) or len(x)!=len(t): raise ValueError('values must have shape (N,) or (N,M)')
    if not np.allclose(dt,dt[0],rtol=1e-8,atol=1e-12*abs(dt[0])):
        raise ValueError('uniform time samples are required')
    if window=='none': w=np.ones(len(t))
    elif window=='hann':
        if len(t)<3: raise ValueError('hann requires at least three samples')
        w=np.hanning(len(t))
    else: raise ValueError('window must be none or hann')
    xw=x*w if x.ndim==1 else x*w[:,None]
    omega=2*np.pi*np.fft.rfftfreq(len(t),dt[0])
    values=np.conj(np.fft.rfft(xw,axis=0))*dt[0]
    phase=np.exp(1j*omega*t[0]);values*=phase if x.ndim==1 else phase[:,None]
    cfg=MappingProxyType(dict(window=window,sign='exp(+i omega t)',normalization='dt * sum',
                              dt=float(dt[0]),start=float(t[0]),samples=len(t),
                              time_unit=time_unit,value_unit=value_unit,
                              frequency_unit=f'radian/({time_unit})',spectrum_unit=f'({value_unit})*({time_unit})'))
    omega=_readonly(omega)
    values=np.frombuffer(values.tobytes(),dtype=values.dtype).reshape(values.shape)
    return Spectrum(omega,values,cfg)
