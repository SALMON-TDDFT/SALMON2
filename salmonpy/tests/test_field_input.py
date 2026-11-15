import unittest
import numpy as np
from salmonpy.errors import InputError
from salmonpy.timing import TimeSpec, InputPort, OutputPort

class FieldInputTests(unittest.TestCase):
    def setUp(self):
        from salmonpy.field_input import FieldValue, FieldPacket, FieldContract, FieldAdapter
        self.Value,self.Packet,self.Contract,self.Adapter=FieldValue,FieldPacket,FieldContract,FieldAdapter

    def value(self, value=(1.,0,0), **kwargs):
        return self.Value(value,time=kwargs.pop('time',.1),position=kwargs.pop('position',(0,0,0)),unit=kwargs.pop('unit','E:au'),**kwargs)

    def contract(self, ports=None, **kwargs):
        return self.Contract(TimeSpec(.1,ports or (InputPort('E'),),(OutputPort('j'),)),
                             units=kwargs.pop('units',{'E':'E:au'}),position=(0,0,0),
                             gauge=kwargs.pop('gauge','length'),**kwargs)

    def test_values_copy_and_none_is_not_zero(self):
        a=np.ones(3);value=self.value(a);a[:]=3
        np.testing.assert_array_equal(value.value,[1,1,1])
        with self.assertRaises(ValueError):value.value[0]=3
        with self.assertRaises(InputError):self.contract().validate(self.Packet(A=self.value(unit='A/c:au')),step_index=1)
        self.contract().validate(self.Packet(E=self.value((0,0,0))),step_index=1)

    def test_metadata_and_numeric_rejection(self):
        for value in ([1j,0,0],[True,False,False],[np.nan,0,0],['1','0','0']):
            with self.subTest(value=value),self.assertRaises(InputError):self.value(value)
        for kwargs in ({'time':True},{'time':.2},{'position':(1,0,0)},{'unit':'H:SI'},{'sampling':'interval_average'}):
            with self.subTest(kwargs=kwargs),self.assertRaises(InputError):
                self.contract().validate(self.Packet(E=self.value(**kwargs)),step_index=1)

    def test_staggered_h_and_custom_callback(self):
        contract=self.contract((InputPort('E'),InputPort('H',offset='1/2')),units={'E':'E:au','H':'SI:A/m'},gauge='custom')
        packet=self.Packet(E=self.value(),H=self.value(time=.15,unit='SI:A/m'))
        class State:
            def step_fields(self,fields):return fields.E.value+fields.H.value
        adapter=self.Adapter(State(),contract)
        np.testing.assert_array_equal(adapter.step_fields(packet),[2,0,0])
        self.assertEqual(adapter.step_index,1)

    def test_legacy_and_failure_contract(self):
        class State:
            def __init__(self):self.calls=0
            def step(self,value):self.calls+=1;return value.copy()
        state=State();adapter=self.Adapter(state,self.contract(),legacy_quantity='E')
        with self.assertRaises(InputError):adapter.step_fields(self.Packet())
        self.assertEqual(state.calls,0)
        np.testing.assert_array_equal(adapter.step_fields(self.Packet(E=self.value())),[1,0,0])
        class Broken:
            def step(self,value):raise RuntimeError('backend')
        adapter=self.Adapter(Broken(),self.contract(),legacy_quantity='E')
        with self.assertRaises(RuntimeError):adapter.step_fields(self.Packet(E=self.value()))
        with self.assertRaises(InputError):adapter.step_fields(self.Packet(E=self.value()))

    def test_gauge_and_port_declarations(self):
        with self.assertRaises(InputError):self.contract(gauge='velocity')
        with self.assertRaises(InputError):self.contract(units={'E':'E:au','H':'H:au'})
        with self.assertRaises(InputError):self.contract((InputPort('j'),),units={'j':'j:au'},gauge='custom')
        with self.assertRaises(InputError):self.Adapter(object(),self.contract(),legacy_quantity='H')

    def test_legacy_rejects_noncanonical_units_and_staggering(self):
        class State:
            def step(self,value):return value
        for contract in (self.contract(units={'E':'SI:V/m'}),
                         self.contract((InputPort('E',offset='1/2'),))):
            with self.assertRaises(InputError):self.Adapter(State(),contract,legacy_quantity='E')

    def test_undeclared_fields_are_not_forwarded(self):
        class State:
            def step_fields(self,packet):
                assert packet.H is None
                return packet.E.value
        adapter=self.Adapter(State(),self.contract())
        adapter.step_fields(self.Packet(E=self.value(),H=self.value(unit='anything',time=99)))

    def test_native_ld_and_sbe_exact_legacy_parity(self):
        import tempfile
        from pathlib import Path
        from salmonpy.ld import LorentzDrude
        from salmonpy.builtin_materials import builtin_material_registry,material_from_catalog
        from salmonpy.config import resolve_materials
        if __package__:
            from .test_builtin_materials import sbe_files
        else:
            from test_builtin_materials import sbe_files
        ld=LorentzDrude(plasma_frequency=.5,strength=[1.],resonance=[.3],damping=[.1])
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            catalog=resolve_materials(sbe_files(root),registry=builtin_material_registry(),base_directory=root)
            sbe=material_from_catalog(catalog,'silicon',base_directory=root)
            for model,q,unit,gauge in ((ld,'E','E:au','phenomenological'),(sbe,'A','A/c:au','velocity')):
                with model.new_state(dt=.1) as direct,model.new_state(dt=.1) as wrapped:
                    timing=TimeSpec(.1,(InputPort(q),),(OutputPort('j'),))
                    contract=self.Contract(timing,units={q:unit},position=(0,0,0),gauge=gauge)
                    adapter=self.Adapter(wrapped,contract,legacy_quantity=q)
                    for n in range(1,5):
                        value=[n*.0001,0,0]
                        expected=direct.step(value)
                        actual=adapter.step_fields(self.Packet(**{q:self.value(value,time=n*.1,unit=unit)}))
                        np.testing.assert_array_equal(actual,expected)

    def test_custom_explicit_metadata_must_match(self):
        class Wrong:
            timing=TimeSpec(1.,(InputPort('H'),),(OutputPort('j'),))
            dt=1.;time=0.
            def step_fields(self,packet):return packet
        with self.assertRaises(InputError):self.Adapter(Wrong(),self.contract())

    def test_custom_clock_offset_and_poststep_check(self):
        timing=TimeSpec(.1,(InputPort('E'),),(OutputPort('j'),),clock_offset='1/2')
        contract=self.Contract(timing,units={'E':'E:au'},position=(0,0,0),gauge='length')
        class State:
            dt=.1;time=.05
            def step_fields(self,packet):self.time+=.1;return packet.E.value
        state=State();state.timing=timing
        adapter=self.Adapter(state,contract)
        adapter.step_fields(self.Packet(E=self.value(time=.15)))
        state.step_fields=lambda packet: packet.E.value  # Broken backend stops its clock.
        with self.assertRaises(InputError):adapter.step_fields(self.Packet(E=self.value(time=.25)))
        with self.assertRaises(InputError):adapter.step_fields(self.Packet(E=self.value(time=.25)))
