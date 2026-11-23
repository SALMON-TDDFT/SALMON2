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

    def test_snapshot_rejects_cross_point_before_backend_callback(self):
        class Saved(State):
            def snapshot(self): return tuple(self.current)
            def restore(self, token):
                self.calls=getattr(self,'calls',0)+1
                self.current=np.asarray(token).copy()
        class SnapshotModel(Model):
            def new_state(self, **kwargs): return Saved()
        with self.Run({'p': SnapshotModel(),'q': SnapshotModel()},run_id='run') as run:
            p=run.new_state('p');q=run.new_state('q')
            p.step([1,2,3]);token=p.snapshot()
            with self.assertRaisesRegex(InputError,'identity'):q.restore(token)
            self.assertEqual(getattr(q._state,'calls',0),0)
            np.testing.assert_array_equal(q.current,np.zeros(3))
            p.step([3,2,1]);p.restore(token)
            np.testing.assert_array_equal(p.current,[1,2,3])
            with self.assertRaises(InputError):p.restore((1,2,3))
            p.close()
            with self.assertRaises(InputError):p.snapshot()
            with self.assertRaises(InputError):p.restore(token)

    def test_snapshot_all_identity_fields_version_and_unsupported(self):
        from dataclasses import replace, FrozenInstanceError
        class Saved(State):
            def snapshot(self):return tuple(self.current)
            def restore(self, token):raise AssertionError('must reject before callback')
        class M(Model):
            def new_state(self, **kwargs):return Saved()
            def record(self):
                return {**super().record(),'initial_state_id':'gs','initial_state_digest':'b'*64}
        with self.Run({'p':M()},run_id='r') as run:
            state=run.new_state('p');saved=state.snapshot()
            for key in state.record():
                wrong=dict(saved._identity);wrong[key]='wrong'
                with self.subTest(key=key),self.assertRaises(InputError):
                    state.restore(replace(saved,_identity=tuple(sorted(wrong.items()))))
            for option in ({'schema':'unknown'},{'version':'unknown'}):
                with self.assertRaises(InputError):state.restore(replace(saved,**option))
            with self.assertRaises(FrozenInstanceError):saved.schema='changed'
            detached=saved.record();detached['point_id']='wrong'
            self.assertEqual(saved.record()['point_id'],'p')
        with self.Run({'p':Model()}) as run:
            with self.assertRaisesRegex(InputError,'support'):run.new_state('p').snapshot()

    def test_ld_sbe_rewind_and_cross_point_independence(self):
        from salmonpy.ld import LorentzDrude
        from salmonpy.sbe import SBEInitialState
        def observe(state):
            return {key:np.array(getattr(state,key),copy=True)
                    for key in ('time','current','polarization','density') if hasattr(state,key)}
        ld=LorentzDrude(plasma_frequency=.7,strength=[.3,.8],resonance=[0.,.9],damping=[.02,.06])
        with np.load(Path(__file__).with_name('sbe_reference.npz')) as data:
            sbe=SBEInitialState(**{key:data[key] for key in ('eigen','momentum','nonlocal_velocity','weights','lattice')},electrons=2)
        class Recorded(Model):
            def __init__(self,raw):self.raw=raw
            def new_state(self,**kw):return self.raw.new_state(**kw)
        for raw in (ld,sbe):
            with self.subTest(backend=type(raw).__name__),self.Run({'p':Recorded(raw),'q':Recorded(raw)},run_id='r') as run:
                p=run.new_state('p',dt=.02);q=run.new_state('q',dt=.02)
                p.step([.01,.02,0]);saved=p.snapshot();untouched=observe(q)
                with self.Run({'p':Recorded(raw)},run_id='r') as other:
                    different_dt=other.new_state('p',dt=.03)
                    before=observe(different_dt)
                    with self.assertRaises(InputError):different_dt.restore(saved)
                    for key,value in before.items():np.testing.assert_array_equal(getattr(different_dt,key),value)
                values=[[.02,.01,0],[.03,0,.01],[0,.01,.02]]
                expected=[]
                for value in values:p.step(value);expected.append(observe(p))
                with self.assertRaises(InputError):q.restore(saved)
                for key,value in untouched.items():np.testing.assert_array_equal(getattr(q,key),value)
                p.restore(saved)
                for value,want in zip(values,expected):
                    p.step(value)
                    for key,actual in observe(p).items():np.testing.assert_array_equal(actual,want[key])
