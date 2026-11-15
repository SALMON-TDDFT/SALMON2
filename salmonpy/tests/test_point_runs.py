import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from salmonpy.errors import InputError


class State:
    def __init__(self):
        self.current = np.zeros(3)
        self.closed = False
    def step(self, value):
        self.current = self.current + value
        return self.current.copy()
    def close(self):
        self.closed = True


class Model:
    def record(self):
        return {'material_id': 'film/α', 'definition_digest': 'a'*64}
    def new_state(self, **kwargs):
        return State()


class PointRunTests(unittest.TestCase):
    def setUp(self):
        from salmonpy.point_runs import PointRun, open_point_arrays
        self.Run = PointRun
        self.open = open_point_arrays

    def test_independent_points_and_lifetime(self):
        model = Model()
        with self.Run({'right': model, 'left': model}, run_id='run/日本') as run:
            a = run.new_state('right', dt=0.1)
            b = run.new_state('left', dt=0.1)
            a.step([1, 2, 3])
            np.testing.assert_array_equal(b.current, [0, 0, 0])
            self.assertEqual(a.record()['point_id'], 'right')
            self.assertEqual(a.run_id, 'run/日本')
            self.assertEqual(a.material_id, 'film/α')
            self.assertEqual(run.point_ids, ('right', 'left'))
            with self.assertRaises(InputError): run.new_state('right')
            a.close()
            with self.assertRaises(InputError): run.new_state('right')
        self.assertTrue(b._state.closed)
        with self.assertRaises(InputError): run.new_state('new')

    def test_records_detached_and_initial_identity(self):
        class Initial(Model):
            def record(self):
                return {**super().record(), 'initial_state_id': 'gs', 'initial_state_digest': 'b'*64}
        with self.Run({'p': Initial()}) as run:
            record = run.record()
            record['points'][0]['material_id'] = 'wrong'
            state = run.new_state('p')
            self.assertEqual(state.record()['material_id'], 'film/α')
            self.assertEqual(state.initial_state_id, 'gs')
            self.assertEqual(state.initial_state_digest, 'b'*64)

    def test_invalid_ids_and_conflicting_definitions(self):
        for points, rid in (({}, 'r'), ({' p': Model()}, 'r'), ({'p': Model()}, ''), ({'p': Model()}, True)):
            with self.subTest(points=points, rid=rid), self.assertRaises(InputError):
                self.Run(points, run_id=rid)
        class Other(Model):
            def record(self): return {**super().record(), 'definition_digest': 'c'*64}
        with self.assertRaises(InputError): self.Run({'p': Model(), 'q': Other()})

    def test_failed_creation_retry_and_shared_state_rejected(self):
        class Fails(Model):
            def new_state(self, **kwargs): raise RuntimeError('failed creation')
        with self.Run({'p': Fails()}) as run:
            for _ in range(2):
                with self.assertRaisesRegex(RuntimeError, 'failed creation'): run.new_state('p')
        raw = State()
        class Shared(Model):
            def new_state(self, **kwargs): return raw
        with self.Run({'p': Shared(), 'q': Shared()}) as run:
            run.new_state('p')
            with self.assertRaises(InputError): run.new_state('q')
            self.assertFalse(raw.closed)

    def test_cleanup_attempts_all_states_and_borrows_models(self):
        calls = []
        class BrokenState(State):
            def close(self):
                calls.append('state')
                raise RuntimeError('cleanup failure')
        class Broken(Model):
            def new_state(self, **kwargs): return BrokenState()
            def close(self): raise AssertionError('model is borrowed')
        run = self.Run({'p': Broken(), 'q': Broken()})
        run.new_state('p')
        run.new_state('q')
        with self.assertRaisesRegex(RuntimeError, 'cleanup failure'): run.close()
        self.assertEqual(calls, ['state', 'state'])
        run.close()

    def test_changed_identity_and_reserved_metadata_rejected(self):
        model = Model()
        with self.Run({'p': model}) as run:
            model.record = lambda: {'material_id': 'other', 'definition_digest': 'a'*64}
            with self.assertRaises(InputError): run.new_state('p')
            with tempfile.TemporaryDirectory() as tmp:
                for metadata in ({'point_run': {}}, {'point_axes': {}}, []):
                    with self.assertRaises(InputError):
                        run.save_arrays(Path(tmp)/'out', {'j': [0]}, point_axes={'j': 0}, metadata=metadata)
                    self.assertFalse((Path(tmp)/'out').exists())

    def test_output_roundtrip_and_axes_rejection(self):
        with tempfile.TemporaryDirectory() as tmp, self.Run({'b': Model(), 'a': Model()}, run_id='r') as run:
            path = Path(tmp)/'result'
            values = np.arange(18).reshape(3, 2, 3)
            archive = run.save_arrays(path, {'j': values, 'time': np.arange(3)}, point_axes={'j': 1}, metadata={'unit': 'test'})
            loaded = self.open(path)
            self.assertEqual(loaded.metadata['point_run']['points'][0]['point_id'], 'b')
            np.testing.assert_array_equal(archive.array('j'), loaded.array('j'))
            for axes in ({}, {'missing': 0}, {'j': 0}, {'j': True}, {'j': -1}):
                with self.subTest(axes=axes), self.assertRaises(InputError):
                    run.save_arrays(Path(tmp)/'bad', {'j': values}, point_axes=axes)
                self.assertFalse((Path(tmp)/'bad').exists())
            manifest = path/'manifest.json'
            data = json.loads(manifest.read_text())
            data['metadata']['point_run']['points'][1]['point_id'] = 'b'
            manifest.write_text(json.dumps(data))
            with self.assertRaises(InputError): self.open(path)
