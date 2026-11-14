"""Exercise public replay, branch isolation and rejected-restore atomicity."""
from dataclasses import FrozenInstanceError, replace
from pathlib import Path
import unittest
import numpy as np
from salmonpy.ld import LorentzDrude
from salmonpy.sbe import SBEInitialState
from salmonpy.maxwell import Transverse, Longitudinal, ThinFilm, Multiscale1D


def ld():
    return LorentzDrude(plasma_frequency=.7, strength=[.3,.8],
                       resonance=[0.,.9], damping=[.02,.06])


def sbe(**options):
    with np.load(Path(__file__).with_name('sbe_reference.npz')) as data:
        return SBEInitialState(**{k:data[k] for k in ('eigen','momentum','nonlocal_velocity','weights','lattice')},
                               electrons=2, **options)


def drive(t):
    t=np.asarray(t)
    return np.column_stack((0*t, .001*np.cos(t), .002*np.sin(t)))


def observe(state):
    values = {'time':state.time}
    for name in ('current','polarization','density','A','grid_A','grid_E','induced','reflected'):
        if hasattr(state,name): values[name]=np.array(getattr(state,name),copy=True)
    return values


class Snapshots(unittest.TestCase):
    def same(self,a,b):
        self.assertEqual(set(a),set(b))
        for name in a: np.testing.assert_array_equal(a[name],b[name],err_msg=name)

    def exercise(self, factory, values):
        a,b=factory(),factory()
        initial=a.snapshot()
        a.step(values[0]); a.step(values[1])
        saved=a.snapshot(); at_save=observe(a)
        b.restore(saved)
        self.same(observe(b),at_save)
        first=[]
        for value in values[2:]:
            a.step(value); first.append(observe(a))
        self.same(observe(b),at_save)
        a.restore(saved)
        self.same(observe(a),at_save)
        for value,expected in zip(values[2:],first):
            a.step(value); self.same(observe(a),expected)
            b.step(value); self.same(observe(b),expected)
        # Saved token does not change when either restored branch advances.
        a.restore(saved); self.same(observe(a),at_save)
        b.restore(initial); self.same(observe(b),observe(factory()))
        # Every metadata rejection must preserve the accepted state.
        before=observe(a)
        for invalid in (None,replace(saved,version='future'),replace(saved,schema='future'),
                        replace(saved,step_index=-1),replace(saved,time=float('nan')),
                        replace(saved,schema=np.array(['bad','schema']))):
            from salmonpy.errors import SnapshotError
            with self.assertRaises(SnapshotError): a.restore(invalid)
            self.same(observe(a),before)
        with self.assertRaises((FrozenInstanceError,AttributeError)):
            saved.time=0
        return saved

    def test_ld_replay_and_independent_branches(self):
        model=ld()
        self.exercise(lambda:model.new_state(dt=.02,initial_field=[.01,.02,0]),
                      np.random.default_rng(11).normal(size=(9,3))*.01)

    def test_sbe_replay_all_corrections(self):
        for order in range(4):
            for nonlocal_ in (False,True):
                with self.subTest(order=order,nonlocal_=nonlocal_):
                    model=sbe(correction_order=order,nonlocal_correction=nonlocal_)
                    self.exercise(lambda:model.new_state(dt=.02,time0=1.,initial_A=[.001,0,0]),
                                  np.random.default_rng(12).normal(size=(9,3))*.01)

    def test_all_fields_replay_and_independence(self):
        models=[Transverse(drive),Longitudinal(drive),ThinFilm(2.,drive,epsilon=(1.,2.))]
        models += [Multiscale1D(cells=2,spacing=3.,vacuum=(3,3),incident=incident,boundaries=boundaries)
                   for incident,boundaries in ((drive,('abc','abc')),(None,('pec','pec')),
                                               (None,('periodic','periodic')))]
        for model in models:
            with self.subTest(model=model):
                factory=lambda:model.new_state(dt=.01,steps=9)
                shape=(9,factory().points,3)
                self.exercise(factory,np.random.default_rng(13).normal(size=shape)*.001)

    def test_incompatible_model_dt_time_capacity_and_drive(self):
        from salmonpy.errors import SnapshotError
        l=ld(); s=sbe(); f=Longitudinal(drive)
        for source,targets in (
            (l.new_state(dt=.02),[ld().new_state(dt=.02),l.new_state(dt=.03)]),
            (s.new_state(dt=.02),[s.new_state(dt=.03),s.new_state(dt=.02,time0=1),sbe().new_state(dt=.02)]),
            (f.new_state(dt=.01,steps=9),[f.new_state(dt=.02,steps=9),f.new_state(dt=.01,steps=8),Transverse(drive).new_state(dt=.01,steps=9)])):
            snap=source.snapshot()
            for target in targets:
                before=observe(target)
                with self.assertRaises(SnapshotError): target.restore(snap)
                self.same(observe(target),before)
        offset=[0.]
        model=Transverse(lambda t:drive(t)+offset[0])
        one=model.new_state(dt=.01,steps=9); snap=one.snapshot()
        offset[0]=.2
        two=model.new_state(dt=.01,steps=9)
        with self.assertRaises(SnapshotError): two.restore(snap)

    def test_exhausted_field_can_rewind(self):
        model=ThinFilm(2.,drive)
        a=model.new_state(dt=.01,steps=2); start=a.snapshot()
        for _ in range(2): a.step(np.ones((1,3)))
        end=a.snapshot(); expected=observe(a)
        with self.assertRaises(ValueError):a.step(np.ones((1,3)))
        a.restore(start)
        for _ in range(2): a.step(np.ones((1,3)))
        self.same(observe(a),expected)
        a.restore(end); self.same(observe(a),expected)

    def test_coupled_sbe_field_replay(self):
        matter=sbe().new_state(dt=.01)
        field=Longitudinal(drive).new_state(dt=.01,steps=10)
        def step(): matter.step(field.step(matter.current[None])[0])
        for _ in range(3):step()
        sm,sf=matter.snapshot(),field.snapshot()
        for _ in range(7):step()
        expected=(observe(matter),observe(field))
        matter.restore(sm);field.restore(sf)
        for _ in range(7):step()
        self.same(observe(matter),expected[0]);self.same(observe(field),expected[1])

    def test_corruption_and_failed_step_preserve_state(self):
        from salmonpy.errors import InputError, SnapshotError
        model=ld(); state=model.new_state(dt=.02)
        state.step([.1,.2,.3]); saved=state.snapshot(); before=observe(state)
        with self.assertRaises(InputError): state.step([float('nan'),0,0])
        self.same(observe(state),before)
        # Immutable storage cannot be made writable through a NumPy view.
        records=saved._data
        for _,_,dtype,_,payload in records:
            view=np.frombuffer(payload,dtype=dtype)
            with self.assertRaises(ValueError): view.flags.writeable=True
        bad=list(records)
        name,shape,dtype,order,payload=bad[-1]
        bad[-1]=(name,shape,dtype,order,bytes([payload[0]^1])+payload[1:])
        for invalid in (replace(saved,_data=tuple(bad)),replace(saved,_data=records[:-1]),
                        replace(saved,_settings='different'),replace(saved,kind='other')):
            with self.assertRaises(SnapshotError):state.restore(invalid)
            self.same(observe(state),before)
        state.step([.4,.5,.6]); state.restore(saved)
        self.same(observe(state),before)
