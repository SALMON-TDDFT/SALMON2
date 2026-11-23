"""Explicit time/space interpolation into A/E packets; synthetic affine fields."""
import numpy as np
from salmonpy.field_sampling import SampledField, sample_packet
from salmonpy.field_input import FieldContract
from salmonpy.timing import TimeSpec, InputPort

positions=np.array([[0.,0.,0.],[2.,0.,0.]])
def history(quantity,unit,times):
    values=np.array([positions+np.array([t,2*t,3*t]) for t in times])
    return SampledField(quantity,times,positions,values,unit)
contract=FieldContract(TimeSpec(.5,(InputPort('A',offset='1/2'),InputPort('E')),()),
                      {'A':'A/c:au','E':'E:au'},(.5,0,0),'custom')
packet=sample_packet(contract,step_index=1,
    histories={'A':history('A','A/c:au',[.25,1.25]),'E':history('E','E:au',[0.,1.])},
    weights={'A':[.75,.25],'E':[.75,.25]},temporal='linear')
np.testing.assert_allclose(packet.A.value,[1.25,1.5,2.25])
np.testing.assert_allclose(packet.E.value,[1.,1.,1.5])
assert packet.H is None
print('A time:',packet.A.time,'E time:',packet.E.time,'position:',packet.E.position)
print('PASS explicit interpolation; no field solver or automatic feedback')
