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

    def test_invalid_execution_mode_rejected_before_callbacks(self):
        provider=SimpleNamespace(identity='initial',prepare=Mock())
        backend=SimpleNamespace(identity='custom',execution='collective',reduction='leader',create_state=Mock())
        with self.assertRaisesRegex(ValueError,'execution_mode'):
            ResponseSession(self.context(),{10:MaterialBinding('custom',provider,backend)},dt=.01,execution_mode='fast')
        provider.prepare.assert_not_called()

    def test_packet_wrapper_keeps_live_optional_declaration(self):
        from salmonpy.material_response import packetize_A_binding
        backend=SimpleNamespace(identity='custom',execution='leader_only',reduction='leader',
                                supported_group_sizes=(1,))
        binding=MaterialBinding('custom',SimpleNamespace(identity='initial'),backend)
        wrapped=packetize_A_binding(binding,field_phase='endpoint')
        self.assertEqual(wrapped.backend.supported_group_sizes,(1,))
        backend.supported_group_sizes=(2,)
        self.assertEqual(wrapped.backend.supported_group_sizes,(2,))

    def test_optional_group_factory_is_session_owned_and_custom(self):
        from unittest.mock import patch
        events=[]
        context=self.context()
        group=SolverGroup('material',(1,2),point_ids=(10,20))
        context.plan=ParallelPlan((context.plan.groups[0],group));context.group=group
        context._dependents=0
        factory=SimpleNamespace(new_state=lambda p:SimpleNamespace(close=lambda:events.append(('close',p))),
                                close=lambda:events.append(('factory',None)))
        provider=SimpleNamespace(identity='initial',prepare=Mock(return_value='prepared'))
        backend=SimpleNamespace(identity='custom',execution='collective',reduction='leader',
                                open_group=Mock(return_value=factory),create_state=Mock(side_effect=AssertionError('legacy route')))
        binding=MaterialBinding('material',provider,backend)
        with patch.object(ResponseSession,'_initial_current_records',return_value={}):
            session=ResponseSession(context,{10:binding,20:binding},dt=.01)
        provider.prepare.assert_called_once();backend.open_group.assert_called_once()
        self.assertIsNot(session.states[10],session.states[20])
        session.close();session.close()
        self.assertEqual(events,[('close',10),('close',20),('factory',None)])
        self.assertEqual(context._dependents,0)
