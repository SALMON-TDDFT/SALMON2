"""Explicit convex spatial and bracketed temporal sampling, without extrapolation."""
from dataclasses import dataclass
from collections.abc import Mapping
import numpy as np
from .errors import InputError, NumericalError
from .field_input import FieldValue, FieldPacket, FieldContract, _real, _label

__all__ = ['SampledField', 'sample_packet']


@dataclass(frozen=True, slots=True, eq=False)
class SampledField:
    """Immutable point samples values[time,site,xyz]; times are atomic units.

    Positions are in position_unit. No conversion between A/c, E and H, or
    interval averages and point values, is inferred.
    """
    quantity: str
    times: object
    positions: object
    values: object
    unit: str
    position_unit: str = 'bohr'

    def __post_init__(self):
        if self.quantity not in ('A','E','H'): raise InputError('A, E or H quantity required')
        _label(self.unit,'unit');_label(self.position_unit,'position_unit')
        t=np.asarray(self.times);p=np.asarray(self.positions)
        if t.ndim!=1 or len(t)==0 or p.ndim!=2 or p.shape[1]!=3 or len(p)==0:
            raise InputError('nonempty times[T] and positions[P,3] required')
        t=_real(t,t.shape,'times');p=_real(p,p.shape,'positions')
        if np.any(t[1:]<=t[:-1]): raise InputError('strictly increasing sample times required')
        v=_real(self.values,(len(t),len(p),3),'sampled field values')
        for key,value in (('times',t),('positions',p),('values',v)):object.__setattr__(self,key,value)

    def sample(self,*,time,position,weights,temporal):
        """Point value with explicit weights and exact/linear time policy.

        Weights must reproduce position as a convex combination of source sites.
        The caller selects the physical stencil. No extrapolation or automatic
        nearest-neighbor choice is made.
        """
        time=float(_real(time,(),'sample time'))
        position=_real(position,(3,),'sample position')
        w=_real(weights,(len(self.positions),),'spatial weights')
        if temporal not in ('exact','linear'): raise InputError('select exact or linear time sampling')
        if np.any(w<0) or not np.isclose(w.sum(),1.,rtol=0,atol=1.e-12):
            raise InputError('nonnegative spatial weights summing to one required')
        scale=np.maximum(1.,np.maximum(np.max(np.abs(self.positions),axis=0),np.abs(position)))
        if not np.allclose(w@self.positions,position,rtol=0,atol=32*np.finfo(float).eps*scale):
            raise InputError('spatial weights do not reproduce requested position')
        # No tolerance snapping: the caller can pass an actual stored timestamp.
        k=int(np.searchsorted(self.times,time))
        if k<len(self.times) and self.times[k]==time:
            values=self.values[k]
        elif temporal=='exact':raise InputError('exact field timestamp is unavailable')
        elif k==0 or k==len(self.times):raise InputError('requested time is outside available history')
        else:
            with np.errstate(over='raise',invalid='raise',divide='raise'):
                try:
                    fraction=(time-self.times[k-1])/(self.times[k]-self.times[k-1])
                    values=(1-fraction)*self.values[k-1]+fraction*self.values[k]
                except FloatingPointError as exc:raise NumericalError('temporal interpolation overflow') from exc
        with np.errstate(over='raise',invalid='raise'):
            try: result=w@values
            except FloatingPointError as exc:raise NumericalError('field interpolation overflow') from exc
        return FieldValue(result,time=time,position=position,unit=self.unit,
                          position_unit=self.position_unit,sampling='point')


def sample_packet(contract,*,step_index,histories,weights,temporal):
    """Sample each required port at its declared time, with per-port stencils.

    histories and weights must cover exactly the required quantities. Only point
    ports are supported. temporal explicitly selects exact or linear for all ports.
    Missing history is rejected; H is never synthesized from other fields.
    """
    if not isinstance(contract,FieldContract):raise InputError('FieldContract required')
    if type(step_index) is not int or step_index<1:raise InputError('positive integer step index required')
    names=set(contract.units)
    if (not isinstance(histories,Mapping) or set(histories)!=names or
            not isinstance(weights,Mapping) or set(weights)!=names):
        raise InputError('histories and weights must cover exactly required ports')
    values={}
    for port in contract.timing.inputs:
        history=histories[port.quantity]
        if (not isinstance(history,SampledField) or history.quantity!=port.quantity or
                history.unit!=contract.units[port.quantity] or history.position_unit!=contract.position_unit
                or port.sampling!='point'):
            raise InputError('field history quantity/units/point sampling must match contract')
        time=contract.timing.t0+(step_index+float(contract.timing.clock_offset+port.offset)-port.lag)*contract.timing.dt
        values[port.quantity]=history.sample(time=time,position=contract.position,
                            weights=weights[port.quantity],temporal=temporal)
    packet=FieldPacket(**values)
    contract.validate(packet,step_index=step_index)
    return packet
