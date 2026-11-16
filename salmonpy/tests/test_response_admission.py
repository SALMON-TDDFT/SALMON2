"""Admission must run without MPI imports or provider/native callbacks."""
from types import SimpleNamespace
from unittest.mock import Mock
import unittest

from salmonpy.distributed import ParallelPlan, SolverGroup
from salmonpy.material_response import MaterialBinding, ResponseSession


class LocalComm:
    def allgather(self, value): return [value]
    def bcast(self, value, root=0): return value


class ResponseAdmissionTests(unittest.TestCase):
    def context(self):
        field=SolverGroup('field',(0,),role='field')
        material=SolverGroup('material',(1,2),point_ids=(10,))
        return SimpleNamespace(plan=ParallelPlan((field,material)),
                               group=material,solver_comm=LocalComm(),leaders=LocalComm(),is_leader=True)

    def test_unsupported_and_malformed_sizes_fail_before_callbacks(self):
        for sizes in ((1,),(),(True,),(0,),(1,1),'1,2'):
            provider=SimpleNamespace(identity='initial',prepare=Mock(side_effect=AssertionError('prepare called')))
            backend=SimpleNamespace(identity='custom',execution='collective',reduction='leader',
                                    supported_group_sizes=sizes,create_state=Mock())
            with self.subTest(sizes=sizes):
                with self.assertRaisesRegex(ValueError,'supported_group_sizes'):
                    ResponseSession(self.context(),{10:MaterialBinding('custom',provider,backend)},dt=.01)
                provider.prepare.assert_not_called()
                backend.create_state.assert_not_called()

    def test_packet_wrapper_keeps_live_optional_declaration(self):
        from salmonpy.material_response import packetize_A_binding
        backend=SimpleNamespace(identity='custom',execution='leader_only',reduction='leader',
                                supported_group_sizes=(1,))
        binding=MaterialBinding('custom',SimpleNamespace(identity='initial'),backend)
        wrapped=packetize_A_binding(binding,field_phase='endpoint')
        self.assertEqual(wrapped.backend.supported_group_sizes,(1,))
        backend.supported_group_sizes=(2,)
        self.assertEqual(wrapped.backend.supported_group_sizes,(2,))
