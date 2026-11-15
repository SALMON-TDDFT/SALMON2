"""Local A/E/H input contracts: native LD parity and custom E/H field delivery.

Synthetic examples, not a native isolated-system TDDFT implementation.
"""
import numpy as np
from salmonpy.field_input import FieldValue, FieldPacket, FieldContract, FieldAdapter
from salmonpy.timing import TimeSpec, InputPort, OutputPort
from salmonpy.ld import LorentzDrude


def main():
    dt=.02
    position=(0.,0.,0.)
    model=LorentzDrude(plasma_frequency=.5,strength=[1.],resonance=[.3],damping=[.1])
    timing=TimeSpec(dt,(InputPort('E'),),(OutputPort('j'),))
    contract=FieldContract(timing,units={'E':'E:au'},position=position,gauge='phenomenological')
    with model.new_state(dt=dt) as state,model.new_state(dt=dt) as reference:
        adapter=FieldAdapter(state,contract,legacy_quantity='E')
        for n in range(1,21):
            value=np.array([.001*np.sin(n*.1),0,0])
            packet=FieldPacket(E=FieldValue(value,time=n*dt,position=position,unit='E:au'))
            actual=adapter.step_fields(packet)
            np.testing.assert_array_equal(actual,reference.step(value))
    # An explicit multiport callback receives separately declared times/units.
    # It only records the supplied values; it is not a magnetic material law.
    class Receiver:
        def step_fields(self, fields):
            return fields.E.value.copy(),fields.H.value.copy()
    timing=TimeSpec(dt,(InputPort('E'),InputPort('H',offset='1/2')),(OutputPort('j'),))
    contract=FieldContract(timing,units={'E':'E:au','H':'SI:A/m'},position=position,gauge='custom')
    adapter=FieldAdapter(Receiver(),contract)
    fields=FieldPacket(E=FieldValue([1,0,0],time=dt,position=position,unit='E:au'),
                       H=FieldValue([0,1,0],time=1.5*dt,position=position,unit='SI:A/m'))
    E,H=adapter.step_fields(fields)
    np.testing.assert_array_equal(E,[1,0,0])
    np.testing.assert_array_equal(H,[0,1,0])
    print('PASS native LD E parity and explicit staggered E/H delivery; no implicit conversions')


if __name__=='__main__':main()
