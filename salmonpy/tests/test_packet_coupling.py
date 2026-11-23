import unittest
from types import SimpleNamespace
import numpy as np
from salmonpy.coupled import MultiscaleField,PacketCoupledStepper
from salmonpy.maxwell import Multiscale1D
from salmonpy.field_input import FieldContract
from salmonpy.timing import TimeSpec,InputPort,OutputPort
from salmonpy.errors import InputError


class PacketCouplingTests(unittest.TestCase):
    def setup_pair(self,*,port=InputPort('A'),output=OutputPort('j'),position=(2,0,0)):
        comm=SimpleNamespace(bcast=lambda x,root=0:x)
        context=SimpleNamespace(group=SimpleNamespace(role='field'),solver_comm=SimpleNamespace(size=1,bcast=comm.bcast),
                                is_leader=True,leaders=comm,field_group=0,plan=SimpleNamespace(point_ids=(10,)))
        field=MultiscaleField(context,Multiscale1D(1,2.,vacuum=(2,2)).new_state(dt=.003,steps=2))
        def agree(error):
            if error:raise InputError(error)
        contract=FieldContract(TimeSpec(.003,(port,),(output,)),{port.quantity: {'A':'A/c:au','E':'E:au','H':'H:au'}[port.quantity]},position,'custom')
        r=SimpleNamespace(context=context,dt=.003,closed=False,failed=False,step_index=0,_field_contracts={10:contract},agree=agree)
        def step_fields(packets):
            contract.validate(packets[10],step_index=r.step_index+1)
            r.step_index+=1
            return np.zeros((1,3))
        r.step_fields=step_fields
        return r,field

    def test_unavailable_clocks_ports_positions_reject_before_update(self):
        for opts in ({'port':InputPort('E')},{'port':InputPort('H')},{'output':OutputPort('j',offset=.5)},{'position':(0,0,0)}):
            r,f=self.setup_pair(**opts)
            with self.subTest(opts=opts),self.assertRaises(InputError):PacketCoupledStepper(r,f)
            self.assertEqual(f.step_index,0);self.assertEqual(r.step_index,0)

    def test_lagged_E_and_endpoint_A(self):
        for port in (InputPort('A'),InputPort('E',lag=1)):
            r,f=self.setup_pair(port=port)
            driver=PacketCoupledStepper(r,f)
            result=driver.step()
            self.assertEqual(result['current_time'],.003)
            self.assertEqual(result['samples']['E'].times[0],0.)
            driver.step()
            with self.assertRaises(InputError):driver.step()
            self.assertEqual(f.step_index,2)

    def test_failure_after_field_advance_invalidates_response(self):
        r,f=self.setup_pair()
        def fail(packets):raise RuntimeError('injected callback error')
        r.step_fields=fail
        driver=PacketCoupledStepper(r,f)
        with self.assertRaises(RuntimeError):driver.step()
        self.assertTrue(r.failed);self.assertEqual(f.step_index,1)
        with self.assertRaises(InputError):driver.step()
