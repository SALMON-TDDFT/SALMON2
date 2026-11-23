import unittest
import numpy as np
from salmonpy.field_sampling import SampledField, sample_packet
from salmonpy.field_input import FieldContract
from salmonpy.timing import TimeSpec,InputPort
from salmonpy.errors import InputError


class SamplingTests(unittest.TestCase):
    def history(self,q='E',times=(0.,1.)):
        positions=np.array([[0.,0,0],[2.,0,0]])
        values=np.array([positions+np.array([t,2*t,3*t]) for t in times])
        return SampledField(q,times,positions,values,{'A':'A/c:au','E':'E:au','H':'H:au'}[q])

    def test_affine_and_ownership(self):
        h=self.history()
        actual=h.sample(time=.25,position=(.5,0,0),weights=(.75,.25),temporal='linear')
        np.testing.assert_allclose(actual.value,[.75,.5,.75])
        with self.assertRaises(ValueError):h.values.setflags(write=True)
        for t in (-.1,1.1):
            with self.assertRaises(InputError):h.sample(time=t,position=(.5,0,0),weights=(.75,.25),temporal='linear')
        with self.assertRaises(InputError):h.sample(time=.25,position=(.5,0,0),weights=(.75,.25),temporal='exact')
        with self.assertRaises(InputError):h.sample(time=0,position=(1,0,0),weights=(.75,.25),temporal='linear')
        with self.assertRaises(InputError):h.sample(time=0,position=(.5,0,0),weights=(-.1,1.1),temporal='linear')

    def test_staggered_ports_and_missing_H(self):
        contract=FieldContract(TimeSpec(.5,(InputPort('A',offset='1/2'),InputPort('E')),()),
                    {'A':'A/c:au','E':'E:au'},(.5,0,0),'custom')
        histories={'A':self.history('A',(.25,1.25)),'E':self.history()}
        packet=sample_packet(contract,step_index=1,histories=histories,
                             weights={'A':(.75,.25),'E':(.75,.25)},temporal='linear')
        self.assertEqual(packet.A.time,.75);self.assertEqual(packet.E.time,.5)
        self.assertIsNone(packet.H)
        np.testing.assert_allclose(packet.A.value,[1.25,1.5,2.25])
        with self.assertRaises(InputError):sample_packet(contract,step_index=3,histories=histories,
                             weights={'A':(.75,.25),'E':(.75,.25)},temporal='linear')
        with self.assertRaises(InputError):sample_packet(contract,step_index=1,histories={'E':histories['E']},weights={'E':(.75,.25)},temporal='linear')

    def test_large_translation_does_not_hide_off_plane_position(self):
        h=SampledField('E',[0.],[[1.e16,0,0],[1.e16+2,0,0]],[[[1,0,0],[2,0,0]]],'E:au')
        with self.assertRaises(InputError):
            h.sample(time=0.,position=[1.e16,1.,0.],weights=[.5,.5],temporal='exact')

    def test_bad_histories(self):
        for times in ([0,0],[1,0],[0,np.nan]):
            with self.assertRaises(InputError):self.history(times=times)
        with self.assertRaises(InputError):SampledField('B',[0],[[0,0,0]],[[[0,0,0]]],'B:au')
