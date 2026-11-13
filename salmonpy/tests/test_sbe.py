import unittest
import numpy as np
from salmonpy import crystal

class SBE(unittest.TestCase):
    def test_native_sbe_independent(self):
        from salmonpy.sbe import SBEInitialState
        p = np.zeros((2,2,3,1), complex)
        p[0,1,0,0] = p[1,0,0,0] = 0.3
        initial = SBEInitialState(eigen=[[-0.5],[0.4]], momentum=p,
            nonlocal_velocity=p*0.1, weights=[1.], lattice=np.eye(3)*4, electrons=2)
        a = initial.new_state(dt=0.02, time0=0., initial_A=[0,0,0])
        b = initial.new_state(dt=0.02, time0=0., initial_A=[0,0,0])
        rho = b.density.copy()
        j = a.step([0.1,0,0])
        np.testing.assert_array_equal(b.density, rho)
        np.testing.assert_array_equal(b.step([0.1,0,0]), j)
        np.testing.assert_array_equal(a.density,b.density)
        before = a.density.copy()
        with self.assertRaises(ValueError): a.step([np.nan,0,0])
        np.testing.assert_array_equal(a.density,before)
        self.assertEqual(a.time,0.02)
        with self.assertRaises(ValueError): initial.eigen.flags.writeable=True
        with self.assertRaises(AttributeError): a.dt=1
        with self.assertRaises(ValueError): initial.new_state(dt=-1,time0=0,initial_A=[0,0,0])

    def test_legacy_fortran_reference(self):
        from pathlib import Path
        from salmonpy.sbe import SBEInitialState
        from salmonpy.fields import FieldSamples
        with np.load(Path(__file__).with_name('sbe_reference.npz')) as data:
            dt=float(data['dt']); drive=FieldSamples(np.arange(len(data['A']))*dt,data['A'],'A')
            for order in range(4):
                for vnl in (False,True):
                    with self.subTest(order=order,vnl=vnl):
                        initial=SBEInitialState(**{k:data[k] for k in ('eigen','momentum','nonlocal_velocity','weights','lattice')},
                            electrons=2,correction_order=order,nonlocal_correction=vnl)
                        index=2*order+int(vnl)
                        result=crystal.respond(initial,drive)
                        np.testing.assert_allclose(result.current,data['current'][index],rtol=2e-13,atol=2e-15)
                        state=initial.new_state(dt=dt,initial_A=data['A'][0])
                        for it,A in enumerate(data['A'][1:],1):
                            state.step(A)
                            np.testing.assert_allclose(state.density,data['density'][index,it],rtol=2e-13,atol=2e-15)

    def test_salmon_text_loader(self):
        from pathlib import Path
        import tempfile
        from salmonpy.sbe import SBEInitialState
        with np.load(Path(__file__).with_name('sbe_reference.npz')) as d, tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); nb,nk=d['eigen'].shape
            (root/'test_k.data').write_text('# header\n'*5+''.join(f'{k+1} 0 0 0 {d["weights"][k]}\n' for k in range(nk))+'# coefficients (2*pi/a [a.u.]) in kx, ky, kz\n1 1 1\n')
            text='# eigen\n# occupation\n# 1:io, 2:esp[a.u.], 3:occ\n'
            for k in range(nk):
                text+=f'k= {k+1}, spin= 1\n'
                for b in range(nb): text+=f'{b+1} {d["eigen"][b,k]} 0.123\n'
            ep=root/'test_eigen.data'; ep.write_text(text)
            tm='# header\n'*3
            for part,name in enumerate(('momentum','nonlocal_velocity')):
                if part: tm+='# nonlocal\n'
                for k in range(nk):
                    for b in range(nb):
                        for c in range(nb):
                            values=d[name][b,c,:,k]
                            tm+=f'{k+1} {b+1} {c+1} '+ ' '.join(f'{z.real:.17e} {z.imag:.17e}' for z in values)+'\n'
            (root/'test_tm.data').write_text(tm.replace('e-','D-'))
            kwargs=dict(bands=nb,kpoints=nk,electrons=2,lattice=d['lattice'])
            initial=SBEInitialState.from_salmon(root,'test',**kwargs)
            np.testing.assert_array_equal(initial.momentum,d['momentum'])
            np.testing.assert_array_equal(initial.eigen,d['eigen'])
            reduced=SBEInitialState.from_salmon(root,'test',active_bands=2,**kwargs)
            self.assertEqual(reduced.new_state(dt=.1).density.shape,(2,2,nk))
            self.assertEqual(initial.new_state(dt=.1).density[0,0,0],2)
            ep.write_text(text.replace('esp[a.u.]','esp[eV]'))
            with self.assertRaisesRegex(ValueError,'a.u.'): SBEInitialState.from_salmon(root,'test',**kwargs)
            ep.write_text(text+'1 0.2 1.0\n')
            with self.assertRaisesRegex(ValueError,'trailing'): SBEInitialState.from_salmon(root,'test',**kwargs)
