"""Replay and branch LD, SBE and Python fields; no full solver required.

Run: python -m salmonpy.examples.snapshot_restore
"""
import json
import numpy as np
from salmonpy.ld import LorentzDrude
from salmonpy.sbe import SBEInitialState
from salmonpy.maxwell import Transverse, Longitudinal, ThinFilm, Multiscale1D


def _check(model, settings, inputs):
    a = model.new_state(**settings)
    b = model.new_state(**settings)
    for value in inputs[:2]: a.step(value)
    saved = a.snapshot()
    b.restore(saved)
    before = b.snapshot()
    expected = [np.asarray(a.step(value)) for value in inputs[2:]]
    assert b.time == before.time
    a.restore(saved)
    for value, result in zip(inputs[2:], expected):
        np.testing.assert_array_equal(a.step(value), result)
        np.testing.assert_array_equal(b.step(value), result)
    b.restore(saved)
    assert b.time == saved.time and a.time > b.time
    return {'replay': 'exact', 'independent_branches': True, 'saved_time': saved.time}


def run():
    values = np.random.default_rng(14).normal(size=(8,3))*.001
    material = LorentzDrude(plasma_frequency=.5,strength=[1.],resonance=[.4],damping=[.02])
    results = {'ld': _check(material, {'dt':.02}, values)}
    momentum = np.zeros((2,2,3,1),complex)
    momentum[0,1,2,0] = momentum[1,0,2,0] = .3
    initial = SBEInitialState(eigen=[[-.5],[.4]], momentum=momentum,
                              nonlocal_velocity=momentum*.1, weights=[1.],
                              lattice=np.eye(3)*4, electrons=2)
    results['sbe'] = _check(initial, {'dt':.02}, values)
    for name, field in (('transverse',Transverse()),('longitudinal',Longitudinal()),
                        ('thin_film',ThinFilm(2.)),
                        ('multiscale_1d',Multiscale1D(cells=2,spacing=3.,vacuum=(3,3)))):
        count = 2 if name == 'multiscale_1d' else 1
        inputs = np.repeat(values[:,None,:],count,axis=1)
        results[name] = _check(field, {'dt':.01,'steps':8},inputs)
    return results


if __name__ == '__main__':
    print(json.dumps(run(),indent=2))
