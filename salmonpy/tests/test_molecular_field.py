import unittest
from dataclasses import replace
import numpy as np
from salmonpy import dipole_source as ds
from salmonpy.field_input import FieldValue
from salmonpy.errors import InputError, NumericalError


def value(v): return FieldValue(v,.1,[1,2,3],'E:au')
def provenance(): return dict(footprint_id='grid-source/1',sample_map_id='map/1',own_source_history_id='history/1')
def reference(): return dict(provenance(),source_ids=['m1'],medium='vacuum',incident=False)


class MolecularField(unittest.TestCase):
    def select(self,total,own,**kw):
        return ds.select_molecular_field(total,molecule_id='m1',policy='environment_scattered_only',
            provenance=kw.pop('provenance',provenance()),reference=own,
            reference_provenance=kw.pop('reference_provenance',reference()),**kw)

    def test_cancellation_external_and_other_molecule(self):
        self.assertTrue(hasattr(ds,'select_molecular_field'))
        own=value([4,5,6])
        np.testing.assert_array_equal(self.select(own,own).field.value,[0,0,0])
        # total = own + illumination + second molecule; subtract only own.
        total=value([7,9,11])
        np.testing.assert_array_equal(self.select(total,own).field.value,[3,4,5])
        np.testing.assert_array_equal(total.value,[7,9,11])

    def test_invalid_reference_and_sample_refused(self):
        for changes in ({'source_ids':['m2']},{'source_ids':['m1','m2']},
                        {'incident':True},{'incident':0},{'medium':'material'},
                        {'footprint_id':'different'},{'sample_map_id':'different'},
                        {'own_source_history_id':'different'}):
            with self.subTest(changes=changes),self.assertRaises(InputError):
                self.select(value([1,0,0]),value([0,0,0]),reference_provenance=dict(reference(),**changes))
        for change in ({'time':.2},{'position':[0,0,0]},{'unit':'V/m'},
                       {'position_unit':'angstrom'},{'sampling':'interval_average'}):
            with self.subTest(change=change),self.assertRaises(InputError):
                self.select(value([1,0,0]),replace(value([0,0,0]),**change))

    def test_provided_drive_and_immutable_metadata(self):
        info={'meaning':'caller-selected cavity field','nested':[1]}
        selected=ds.select_molecular_field(value([1,2,3]),molecule_id='m1',
            policy='provided_drive',provenance=info)
        info['nested'][0]=99
        self.assertEqual(selected.record()['provenance']['nested'],[1])
        record=selected.record();record['provenance']['nested'][0]=88
        self.assertEqual(selected.record()['provenance']['nested'],[1])
        with self.assertRaises(ValueError): selected.field.value[0]=0
        with self.assertRaises(TypeError): selected.metadata['policy']='other'
        with self.assertRaises(InputError):
            ds.select_molecular_field(value([1,0,0]),molecule_id='m1',policy='provided_drive',provenance={})
        with self.assertRaises(InputError):
            ds.select_molecular_field(value([1,0,0]),molecule_id='m1',policy='provided_drive',
                provenance=info,reference=value([0,0,0]))

    def test_overflow_missing_reference_and_unknown_policy(self):
        with self.assertRaises(NumericalError): self.select(value([1e308,0,0]),value([-1e308,0,0]))
        with self.assertRaises(InputError): self.select(value([1,0,0]),None)
        with self.assertRaises(InputError):
            ds.select_molecular_field(value([1,0,0]),molecule_id='m1',policy='automatic',provenance={})
