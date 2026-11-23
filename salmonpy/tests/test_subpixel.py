"""Reference-frequency material matching and native interface geometry."""
import unittest
import numpy as np

class SubpixelTests(unittest.TestCase):
    def test_complex_interface_and_passive_lift(self):
        from salmonpy.subpixel import effective_material
        from salmonpy.ld import LorentzDrude
        left=LorentzDrude(plasma_frequency=1.,strength=[1.],resonance=[0.],damping=[.08])
        right=LorentzDrude(plasma_frequency=0.,strength=[0.],resonance=[0.],damping=[0.],epsilon_infinity=2.3)
        omega=.4;fraction=.35;n=np.array([1.,2.,3.]);n/=np.linalg.norm(n)
        model,report=effective_material(left,right,fraction=fraction,normal=n,reference_omega=omega,units='atomic_gaussian')
        e1=left.permittivity([omega])[0];e2=right.permittivity([omega])[0]
        tangent=fraction*e1+(1-fraction)*e2
        perpendicular=1/(fraction/e1+(1-fraction)/e2)
        target=tangent*np.eye(3)+(perpendicular-tangent)*np.outer(n,n)
        np.testing.assert_allclose(model.permittivity_tensor([omega])[0],target,rtol=3e-13,atol=3e-13)
        self.assertLess(report['reference_error'],1e-12)
        self.assertGreater(np.linalg.eigvalsh(model.config.epsilon_infinity).min(),0.)
        for f in model.config.strength:self.assertGreaterEqual(np.linalg.eigvalsh(f).min(),-1e-13)
        for f,original in [(0.,right),(1.,left)]:
            endpoint,_=effective_material(left,right,fraction=f,normal=n,reference_omega=omega,units='atomic_gaussian')
            np.testing.assert_allclose(endpoint.permittivity_tensor([.2,.6]),original.permittivity_tensor([.2,.6]),atol=1e-13)

    def test_native_plane_and_shape_cut(self):
        from salmonpy.subpixel import _plane_fraction,_shape_cut
        from salmonpy.geometry import Shape
        for n in ([1,0,0],[1,1,1],[1,1e-12,0]):
            self.assertAlmostEqual(_plane_fraction(n,0.,[1,1,1]),.5,places=13)
        self.assertAlmostEqual(_plane_fraction([1,0,0],.25,[1,1,1]),.75,places=13)
        cut=_shape_cut([Shape('ellipsoid',(4.,4.,4.))],[2.,0.,0.],[.1,.1,.1])
        self.assertEqual(cut['flag'],1)
        self.assertAlmostEqual(cut['fraction'],.5,places=13)
        np.testing.assert_array_equal(cut['normal'],[1.,0.,0.])
        self.assertEqual(tuple(cut['materials']),(1,0))
        adjacent=_shape_cut([Shape('box',(2.,2.,1.),center=(0.,0.,-.5),material=1),
            Shape('box',(2.,2.,1.),center=(0.,0.,.5),material=2)],[0.,0.,0.],[1.,1.,1.])
        self.assertEqual(adjacent['flag'],1)
        self.assertEqual(tuple(adjacent['materials']),(2,1))
        self.assertAlmostEqual(adjacent['fraction'],.5,places=13)


    def test_fdtd_carrier_reference_asset(self):
        from salmonpy.subpixel import SubpixelSmoothing
        from salmonpy.fdtd import FDTD,PlaneWave
        from salmonpy.geometry import Grid,Geometry,Shape
        from salmonpy.fields import AcosPulse
        from salmonpy.ld import LorentzDrude
        material=LorentzDrude(plasma_frequency=.3,strength=[1.],resonance=[.4],damping=[.02])
        geometry=Geometry(Grid((40,40,40),(2.,)*3),(Shape('ellipsoid',(8.,)*3),))
        calculation=FDTD(geometry,(material,),.001,10,
            source=PlaneWave(AcosPulse(.001,.5,2.),0,-16.),smoothing=SubpixelSmoothing())
        assets=calculation.input_assets()
        self.assertIn('salmonpy-subpixel.dat',assets)
        self.assertIn('SALMONPY_SUBPIXEL_REFERENCE 2',assets['salmonpy-subpixel.dat'])
        self.assertIn('ld-subpixel',calculation.input_text())

        from dataclasses import replace
        from salmonpy.fdtd import Dielectric
        with self.assertWarnsRegex(RuntimeWarning,'only.*LD'):
            pure=replace(calculation,materials=(Dielectric(2.),))
        self.assertNotIn('salmonpy-subpixel.dat',pure.input_assets())

        from salmonpy.fdtd import Impulse
        with self.assertWarnsRegex(RuntimeWarning,'only.*LD'):
            with self.assertRaisesRegex(ValueError,'impulse excitation requires LD'):
                replace(calculation,materials=(Dielectric(2.),),source=Impulse())
