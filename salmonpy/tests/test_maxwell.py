import tempfile
import unittest
from pathlib import Path
import numpy as np
from salmonpy.fields import AcosPulse, FieldSamples
from salmonpy.maxwell import (SPEED_OF_LIGHT, Transverse, Longitudinal, ThinFilm,
                              Multiscale1D, couple)
from salmonpy.maxwell import _evaluate
from salmonpy.outputs import OutputArchive

C = SPEED_OF_LIGHT


class FreeElectrons:
    """Test double: diamagnetic-only electrons. SALMON matter current J = n A; returns j = -n A."""
    backend = 'free-electron-test'

    def __init__(self, density):
        self.density = density
        self.closed = []

    def new_state(self, *, dt, time0, initial_A):
        state = _FreeState(self, dt, time0, initial_A)
        return state


class _FreeState:
    def __init__(self, owner, dt, time0, A):
        from salmonpy.timing import crystal_default
        self.owner = owner; self.dt = dt; self.t0 = time0; self.k = 0
        self.timing = crystal_default(dt, time0)
        self.J = -owner.density*np.asarray(A, float)

    @property
    def time(self): return self.t0 + self.k*self.dt

    @property
    def current(self): return self.J.copy()

    def step(self, A):
        self.k += 1; self.J = -self.owner.density*np.asarray(A, float)
        return self.J.copy()

    def close(self): self.owner.closed.append(self)


def impulse(t):
    return np.where(np.asarray(t)[:, None] >= 0, 1.0, 0.0)*np.array([0., 0., 1e-3])


class UniformFieldTests(unittest.TestCase):
    def test_exact_endpoint_samples_cover_all_steps(self):
        dt = .1
        samples = FieldSamples(np.arange(11)*dt, np.zeros((11, 3)), 'A')
        for model in (Transverse(samples), Longitudinal(samples), ThinFilm(1., samples)):
            with self.subTest(model=type(model).__name__):
                state = model.new_state(dt=dt, steps=10)
                for _ in range(10):
                    np.testing.assert_array_equal(state.step(np.zeros((1, 3))), np.zeros((1, 3)))
                self.assertAlmostEqual(state.time, 1.)
                with self.assertRaises(ValueError):
                    state.step(np.zeros((1, 3)))
                with self.assertRaises(ValueError):
                    model.new_state(dt=dt, steps=11)

    def test_incident_requires_left_absorbing_boundary(self):
        pulse = AcosPulse(.001, .2, 10.)
        for boundaries in (('pec', 'abc'), ('periodic', 'periodic')):
            with self.subTest(boundaries=boundaries):
                with self.assertRaisesRegex(ValueError, 'incident.*abc'):
                    Multiscale1D(1, 1., incident=pulse, boundaries=boundaries)
                Multiscale1D(1, 1., boundaries=boundaries).new_state(dt=.001, steps=2)
        Multiscale1D(1, 1., incident=pulse, boundaries=('abc', 'pec')).new_state(dt=.001, steps=2)

    def test_transverse_returns_external_and_ignores_current(self):
        pulse = AcosPulse(amplitude=1e-3, omega=.2, duration=20., polarization=(1, 0, 0))
        state = Transverse(pulse).new_state(dt=.1, steps=300)
        expected = pulse.sample(np.arange(302)*.1).values
        np.testing.assert_array_equal(state.A[0], expected[0])
        for n in range(300):
            A = state.step(np.full((1, 3), 5.))
            np.testing.assert_array_equal(A[0], expected[n+1])
        self.assertEqual(state.config.model, 'transverse')
        with self.assertRaises(ValueError): state.step(np.zeros((1, 3)))  # capacity

    def test_longitudinal_plasma_oscillation(self):
        # J = n A gives d^2 A_ind/dt^2 = -4 pi n A_tot: screening oscillation at omega_p.
        n, dt = .01, .05
        result = couple(FreeElectrons(n), Longitudinal(impulse), dt=dt, steps=4000)
        a = result.A[:, 0, 2]
        self.assertLess(np.abs(a).max(), 1.001e-3)  # bounded: correct sign
        crossings = np.where(np.diff(np.sign(a[10:])) != 0)[0]
        period = 2*np.mean(np.diff(crossings))*dt
        wp2 = 4*np.pi*n
        discrete = np.arccos(1 - wp2*dt**2/2)/dt  # leapfrog dispersion
        self.assertAlmostEqual(period, 2*np.pi/discrete, delta=5e-3*period)

    def test_longitudinal_recurrence_and_initial_current(self):
        state = Longitudinal().new_state(dt=.1, steps=3)
        # SALMON: J(0) is not fed back, Ac_ind(1) = 0.
        np.testing.assert_array_equal(state.step([[1., 2., 3.]]), [[0., 0., 0.]])
        A2 = state.step([[1., 0., 0.]])
        np.testing.assert_allclose(A2, [[4*np.pi*.01, 0, 0]], rtol=1e-15)  # j: charge current
        A3 = state.step([[0., 0., 0.]])
        np.testing.assert_allclose(A3, 2*A2, rtol=1e-15)

    def test_thin_film_fresnel_without_current(self):
        pulse = AcosPulse(amplitude=1e-3, omega=.2, duration=30., start=1.)
        state = ThinFilm(thickness=10., external=pulse, epsilon=(1., 4.)).new_state(dt=.1, steps=400)
        ext = pulse.sample(np.arange(402)*.1).values
        r = (1-2)/(1+2)
        out = [state.step(np.zeros((1, 3)))[0] for _ in range(400)]
        np.testing.assert_allclose(out[1:], (1+r)*ext[2:401] - r*ext[1], rtol=0, atol=1e-14)
        np.testing.assert_allclose(state.reflected, r*(ext[400]-ext[1]), rtol=0, atol=1e-14)

    def test_thin_film_sheet_transmission(self):
        # j = -n A_t: dA_ref/dt = -(2 pi d / c) n A_t for n1=n2=1, t(w) = -i w / (-i w + gamma).
        d, n, w0, dt = 5., .05, .1, .05
        pulse = AcosPulse(amplitude=1e-3, omega=w0, duration=1500., polarization=(0, 0, 1))
        result = couple(FreeElectrons(n), ThinFilm(thickness=d, external=pulse), dt=dt,
                        steps=int(1700/dt))
        t = result.time
        ft = lambda x: dt*np.sum(x*np.exp(1j*w0*t))
        numeric = ft(result.A[:, 0, 2])/ft(pulse.sample(t).values[:, 2])
        gamma = 2*np.pi*d*n/C
        self.assertLess(abs(numeric - (-1j*w0)/(-1j*w0+gamma)), 3e-3)

    def test_field_samples_external_must_match_grid(self):
        samples = FieldSamples(np.arange(12)*.1, np.ones((12, 3)), 'A')
        state = Transverse(samples).new_state(dt=.1, steps=10)
        np.testing.assert_array_equal(state.step(np.zeros((1, 3))), [[1, 1, 1]])
        with self.assertRaises(ValueError): Transverse(samples).new_state(dt=.1, steps=12)
        with self.assertRaises(ValueError): Transverse(samples).new_state(dt=.07, steps=5)
        with self.assertRaises(ValueError): Transverse(FieldSamples(np.arange(3.), np.ones((3, 3)), 'E'))


class Multiscale1DTests(unittest.TestCase):
    def test_vacuum_propagation_at_magic_time_step(self):
        dt = .05; h = C*dt
        pulse = AcosPulse(amplitude=1e-3, omega=.1, duration=200., polarization=(0, 0, 1))
        model = Multiscale1D(cells=5, spacing=h, vacuum=(20, 20), incident=pulse)
        state = model.new_state(dt=dt, steps=6000)
        A = [state.A]
        for _ in range(6000): A.append(state.step(np.zeros((5, 3))))
        A = np.array(A); t = np.arange(6001)*dt
        expected = np.stack([_evaluate(pulse, t - i*h/C) for i in range(1, 6)], axis=1)
        scale = np.abs(expected).max()
        self.assertLess(np.abs(A - expected).max(), 1e-5*scale)
        self.assertLess(np.abs(state.grid_A).max(), 1e-6*scale)  # transmitted, not reflected
        np.testing.assert_array_equal(model.x, np.arange(-19, 26)*h)

    def test_one_cell_matches_thin_film(self):
        dt = .05; h = C*dt; n = .05
        pulse = AcosPulse(amplitude=1e-3, omega=.1, duration=200., polarization=(0, 0, 1))
        r1 = couple(FreeElectrons(n), Multiscale1D(cells=1, spacing=h, vacuum=(20, 20), incident=pulse),
                    dt=dt, steps=6000)
        r2 = couple(FreeElectrons(n), ThinFilm(thickness=h, external=lambda t: _evaluate(pulse, t-h/C)),
                    dt=dt, steps=6000)
        a1, a2 = r1.A[:, 0, 2], r2.A[:, 0, 2]
        self.assertLess(np.abs(a1-a2).max(), 2e-3*np.abs(a2).max())

    def test_validation(self):
        with self.assertRaises(ValueError): Multiscale1D(cells=2, spacing=1.).new_state(dt=1., steps=2)
        with self.assertRaises(ValueError): Multiscale1D(cells=2, spacing=1., boundaries=('periodic', 'abc'))
        with self.assertRaises(ValueError): Multiscale1D(cells=0, spacing=1.)
        with self.assertRaises(TypeError): Multiscale1D(cells=1, spacing=1., incident=FieldSamples(
            np.arange(3.), np.zeros((3, 3)), 'A'))
        with self.assertRaises(ValueError):  # x is the propagation axis
            Multiscale1D(cells=1, spacing=C*.05, incident=AcosPulse(1e-3, .1, 50., polarization=(1, 0, 0))
                         ).new_state(dt=.05, steps=2)
        state = Multiscale1D(cells=2, spacing=1., boundaries=('periodic', 'periodic')).new_state(
            dt=.5/C, steps=5)
        np.testing.assert_array_equal(state.step(np.zeros((2, 3))), np.zeros((2, 3)))

    def test_failed_step_leaves_state_unchanged(self):
        model = Multiscale1D(cells=2, spacing=C*.05, vacuum=(3, 3))
        state = model.new_state(dt=.05, steps=4)
        state.step(np.zeros((2, 3)))
        before = (state.step_index, state.A.copy(), state.grid_A.copy())
        for bad in (np.zeros((3, 3)), np.full((2, 3), np.nan), np.zeros(3)):
            with self.assertRaises(ValueError): state.step(bad)
        self.assertEqual(state.step_index, before[0])
        np.testing.assert_array_equal(state.A, before[1])
        np.testing.assert_array_equal(state.grid_A, before[2])
        with self.assertRaises(ValueError): state.A[0, 0] = 1.  # read-only snapshot


class CoupleTests(unittest.TestCase):
    def test_one_independent_state_per_point_and_closed(self):
        matter = FreeElectrons(.01)
        result = couple(matter, Multiscale1D(cells=3, spacing=C*.05, vacuum=(2, 2)), dt=.05, steps=4)
        self.assertEqual(result.A.shape, (5, 3, 3)); self.assertEqual(result.current.shape, (5, 3, 3))
        self.assertEqual(len(matter.closed), 3); self.assertEqual(len({id(s) for s in matter.closed}), 3)
        self.assertEqual(result.backends, ('free-electron-test',)*3)
        with self.assertRaises(ValueError):
            couple([matter, matter], Multiscale1D(cells=3, spacing=C*.05), dt=.05, steps=2)

    def test_states_closed_on_failure(self):
        class Bad(FreeElectrons):
            def new_state(self, **kw):
                s = super().new_state(**kw)
                s.step = lambda A: (_ for _ in ()).throw(RuntimeError('boom'))
                return s
        matter = Bad(.01)
        with self.assertRaises(RuntimeError): couple(matter, Longitudinal(impulse), dt=.1, steps=3)
        self.assertEqual(len(matter.closed), 1)

    def test_output_archive(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)/'run'
            result = couple(FreeElectrons(.01), Longitudinal(impulse), dt=.1, steps=10, output_directory=out)
            archive = OutputArchive.open(out)
            np.testing.assert_array_equal(archive.array('A'), result.A)
            self.assertEqual(archive.metadata['field']['model'], 'longitudinal')
            self.assertEqual(archive.metadata['timing']['field']['schema'], 'salmonpy.time/1')
            with self.assertRaises(FileExistsError):
                couple(FreeElectrons(.01), Longitudinal(), dt=.1, steps=2, output_directory=out)


class NativeSBECouplingTests(unittest.TestCase):
    """Driver wiring with a real matter backend (bundled SBE kernel)."""
    def setUp(self):
        from salmonpy.sbe import SBEInitialState
        try:
            from salmonpy.sbe import _native
            _native()
        except RuntimeError as exc:
            self.skipTest(str(exc))
        with np.load(Path(__file__).with_name('sbe_reference.npz')) as data:
            self.dt = float(data['dt']); self.A = data['A'].copy()
            self.initial = SBEInitialState(**{k: data[k] for k in ('eigen', 'momentum', 'nonlocal_velocity',
                                           'weights', 'lattice')}, electrons=2, correction_order=1)

    def test_transverse_coupling_equals_crystal_respond(self):
        from salmonpy import crystal
        steps = len(self.A)-2
        drive = FieldSamples(np.arange(len(self.A))*self.dt, self.A, 'A')
        coupled = couple(self.initial, Transverse(drive), dt=self.dt, steps=steps)
        reference = crystal.respond(self.initial, FieldSamples(drive.time[:steps+1], self.A[:steps+1], 'A'))
        np.testing.assert_array_equal(coupled.current[:, 0], reference.current)
        np.testing.assert_array_equal(coupled.A[:, 0], self.A[:steps+1])
        self.assertIn('Taylor', coupled.timing['matter'][0]['integrator'])

    def test_longitudinal_coupling_follows_salmon_recurrence(self):
        steps = len(self.A)-2
        drive = FieldSamples(np.arange(len(self.A))*self.dt, self.A, 'A')
        result = couple(self.initial, Longitudinal(drive), dt=self.dt, steps=steps)
        ind = np.zeros((steps+1, 3))
        for n in range(1, steps):
            ind[n+1] = 2.0*ind[n] - ind[n-1] + 4.0*np.pi*result.current[n, 0]*self.dt**2
        np.testing.assert_array_equal(result.A[:, 0], self.A[:steps+1] + ind)


if __name__ == '__main__':
    unittest.main()
