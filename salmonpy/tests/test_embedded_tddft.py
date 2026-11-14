import importlib.util
import unittest


class EmbeddedTDDFTTest(unittest.TestCase):
    def test_collective_backend_is_available(self):
        self.assertIsNotNone(importlib.util.find_spec('salmonpy.embedded_tddft'),
                             'embedded TDDFT backend is missing')
        from salmonpy.embedded_tddft import NativeTDDFTBackend
        self.assertEqual(NativeTDDFTBackend.execution, 'collective')
        self.assertEqual(NativeTDDFTBackend.reduction, 'leader')

    def test_loader_detects_changed_assets(self):
        import tempfile
        from pathlib import Path
        from salmonpy.embedded_tddft import NativeTDDFTMaterialProvider
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'data_for_restart').mkdir()
            for name in ('info.bin','occupation.bin','wfn.bin'):
                (root/'data_for_restart'/name).write_bytes(b'loader-only fixture')
            (root/'rt.inp').write_text('&calculation\n theory="tddft_pulse"\n/\n')
            provider=NativeTDDFTMaterialProvider(root,'rt.inp',max_steps=4)
            self.assertEqual(provider.prepare(None).max_steps,4)
            (root/'data_for_restart/wfn.bin').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'changed'):provider.prepare(None)
            with self.assertRaisesRegex(ValueError,'relative'):
                NativeTDDFTMaterialProvider(root,'../rt.inp')

    def test_import_does_not_initialize_mpi(self):
        import subprocess,sys
        code="import sys; import salmonpy.embedded_tddft; assert 'mpi4py.MPI' not in sys.modules"
        subprocess.run([sys.executable,'-c',code],check=True)

    def test_native_eh_field_interface_is_available(self):
        self.assertIsNotNone(importlib.util.find_spec('salmonpy.coupled'),
                             'native E/H coupling module is missing')
        from salmonpy.coupled import NativeEHField
        self.assertEqual(NativeEHField.scheme,'salmonpy.eh-centered-current/1')

    @staticmethod
    def _fake_coupling(matter_record):
        """Single field rank; leaders.allgather also returns one material group's record."""
        from types import SimpleNamespace
        from salmonpy.coupled import NativeEHField
        one=SimpleNamespace(bcast=lambda x,root=0:x)
        leaders=SimpleNamespace(allgather=lambda x:[x,{1:matter_record}],bcast=lambda x,root=0:x)
        context=SimpleNamespace(group=SimpleNamespace(role='field'),plan=SimpleNamespace(point_ids=(1,)),
                                is_leader=True,leaders=leaders,solver_comm=one,field_group=0)
        def agree(error):
            if error:raise ValueError(error)
        response=SimpleNamespace(context=context,dt=.1,step_index=0,closed=False,failed=False,agree=agree,
                                 timing_records=lambda:{1:matter_record})
        field=object.__new__(NativeEHField)
        field.context=context;field.dt=.1;field.point_ids=(1,);field.step_index=0;field.closed=False;field.failed=False
        return response,field

    def test_coupling_checks_matter_timing_on_the_electronic_clock(self):
        from salmonpy.coupled import CoupledStepper
        from salmonpy.material_response import _response_timing
        from salmonpy.timing import TimeSpec,InputPort,OutputPort
        good=_response_timing(.1,'test').to_dict()
        response,field=self._fake_coupling(good)
        stepper=CoupledStepper(response,field)
        self.assertEqual(stepper.timing['electronic_clock_offset'],'1/2')
        self.assertEqual(stepper.timing['field']['outputs']['A']['offset'],'1/2')
        # A backend already on the half-step clock would be shifted to 1: rejected.
        shifted=TimeSpec(.1,(InputPort('A'),),(OutputPort('j'),),clock_offset='1/2').to_dict()
        # A backend that samples A as an interval average does not match the field's point A.
        averaged=TimeSpec(.1,(InputPort('A',sampling='interval_average'),),(OutputPort('j'),)).to_dict()
        for bad in (shifted,averaged):
            response,field=self._fake_coupling(bad)
            with self.assertRaises(ValueError):CoupledStepper(response,field)

    def test_coupling_rejects_out_of_band_advancement_before_field_update(self):
        from salmonpy.coupled import CoupledStepper
        from salmonpy.material_response import _response_timing
        response,field=self._fake_coupling(_response_timing(.1,'test').to_dict())
        def forbidden_update(current):raise AssertionError('field must not be mutated')
        field.advance=forbidden_update
        stepper=CoupledStepper(response,field)
        response.step_index=1
        with self.assertRaisesRegex(ValueError,'synchron'):
            stepper.step()
        self.assertTrue(response.failed)

    def test_coupled_sites_cannot_overlap_incident_source(self):
        from types import SimpleNamespace
        from salmonpy.coupled import NativeEHField
        from salmonpy.fdtd import FDTD,PlaneWave
        from salmonpy.geometry import Grid,Geometry
        from salmonpy.fields import AcosPulse
        comm=SimpleNamespace(size=1,allgather=lambda x:[x])
        context=SimpleNamespace(solver_comm=comm,group=SimpleNamespace(role='field'),plan=SimpleNamespace(point_ids=(1,)))
        calculation=FDTD(Geometry(Grid((40,40,40),(2.,2.,2.))),(),.003,2,
                         source=PlaneWave(AcosPulse(.001,20.,.15),0,-1.))
        with self.assertRaisesRegex(ValueError,'source'):
            NativeEHField(context,calculation,'unused-record','unused-output',points={1:(19,20,20)},grid=(1,1,1))

    def test_field_cannot_retry_a_failed_native_step(self):
        from types import SimpleNamespace
        import numpy as np
        from salmonpy.coupled import NativeEHField
        field=object.__new__(NativeEHField)
        field.context=SimpleNamespace(solver_comm=SimpleNamespace(allgather=lambda x:[x]))
        field.lib=SimpleNamespace(sp_eh_step=object());field.capacity=4;field.point_ids=(1,)
        field.A=np.zeros((1,3));field.E=field.A.copy();field.dt=.1
        field.step_index=0;field.closed=False;field.failed=False;field.manifest={}
        calls=[]
        def broken_native(fn,n,j,e):
            calls.append(1);e[:]=np.nan;return 0
        field._call=broken_native
        with self.assertRaises(ValueError):field.advance(np.zeros((1,3)))
        self.assertTrue(field.failed)
        with self.assertRaises(ValueError):field.advance(np.zeros((1,3)))
        self.assertEqual(len(calls),1)
