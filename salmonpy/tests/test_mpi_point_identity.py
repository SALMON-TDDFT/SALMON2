import copy
import numpy as np
import unittest
from salmonpy.distributed import ParallelPlan, SolverGroup
from salmonpy.errors import InputError


class MPIIdentityTests(unittest.TestCase):
    def setUp(self):
        from salmonpy.material_response import _point_metadata
        self.check = _point_metadata
        self.plan = ParallelPlan((SolverGroup('field', (0,), role='field'),
                                  SolverGroup('left', (1,), point_ids=(10,)),
                                  SolverGroup('right', (2,), point_ids=(90,))))
        self.record = {'schema': 'salmonpy.point-run/1', 'run_id': 'run/α', 'points': [
            {'point_id': name, 'material_id': 'film', 'definition_digest': 'a'*64}
            for name in ('right', 'left')]}
        self.mapping = {10: 'right', 90: 'left'}

    def test_mapping_preserves_explicit_order_and_inputs(self):
        result = self.check(self.plan, self.record, self.mapping)
        self.assertEqual(result['point_run']['points'][0]['point_id'], 'right')
        self.assertEqual(result['mpi_point_map'][1]['transport_point_id'], 90)
        self.assertEqual(result['mpi_point_map'][1]['group_name'], 'right')
        result['point_run']['points'][0]['material_id'] = 'changed'
        self.assertEqual(self.record['points'][0]['material_id'], 'film')

    def test_legacy_path_and_incomplete_input(self):
        self.assertIsNone(self.check(self.plan, None, None))
        for record, mapping in ((None, self.mapping), (self.record, None)):
            with self.assertRaises(InputError): self.check(self.plan, record, mapping)

    def test_missing_duplicated_or_reordered_ids_rejected(self):
        for mapping in ({10:'right'}, {10:'left',90:'right'}, {10:'right',90:'right'},
                        {'10':'right',90:'left'}, {10:np.array('right'),90:'left'},
                        {10:np.str_('right'),90:'left'}):
            with self.subTest(mapping=mapping), self.assertRaises(InputError):
                self.check(self.plan, self.record, mapping)
        bad = copy.deepcopy(self.record)
        bad['points'][1]['definition_digest'] = 'b'*64
        with self.assertRaises(InputError): self.check(self.plan, bad, self.mapping)

    def test_saved_transport_mapping_tamper_rejected(self):
        import tempfile
        from pathlib import Path
        from salmonpy.point_runs import PointRun, open_point_arrays
        from salmonpy.outputs import save_arrays
        metadata = self.check(self.plan, self.record, self.mapping)
        metadata['point_axes'] = {'current': 0}
        metadata['mpi_point_map'][1]['transport_point_id'] = 10
        with tempfile.TemporaryDirectory() as tmp:
            save_arrays(Path(tmp)/'bad', {'current': [[0,0,0],[0,0,0]]}, metadata=metadata)
            with self.assertRaises(InputError): open_point_arrays(Path(tmp)/'bad')

    def test_empty_message_writer_failure_is_not_success(self):
        from types import SimpleNamespace
        from unittest.mock import patch
        from salmonpy.material_response import ResponseSession
        class Comm:
            def allgather(self, value): return [value]
            def bcast(self, value, root): return value
        session = object.__new__(ResponseSession)
        session.closed = session.failed = False
        session._point_metadata = self.check(self.plan, self.record, self.mapping)
        session.context = SimpleNamespace(group=SimpleNamespace(role='field'), is_leader=True,
                                          plan=self.plan, solver_comm=Comm(), leaders=Comm())
        with patch('salmonpy.outputs.save_arrays', side_effect=RuntimeError()):
            with self.assertRaisesRegex(InputError, 'RuntimeError'):
                session.save_arrays('unused', {'j': [[0,0,0],[0,0,0]]}, point_axes={'j':0})
