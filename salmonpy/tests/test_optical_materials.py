"""Offline, licensed source integrity and reusable material records."""
import unittest
import numpy as np

class OpticalMaterialsTests(unittest.TestCase):
    def test_sources_catalog_and_models(self):
        from salmonpy.optical_materials import list_materials, load_material
        from unittest.mock import patch
        with patch('urllib.request.urlopen',side_effect=AssertionError('network on listing')):
            items=list_materials()
        self.assertEqual({r['symbol'] for r in items},{'Au','Ag','Cu','Al','Si','SiO2','Al2O3','TiO2','Ta2O5','HfO2'})
        for row in items:
            material=load_material(row['id'])
            self.assertEqual(material.metadata['license'],'CC0-1.0')
            self.assertEqual(len(material.metadata['source_sha256']),64)
            self.assertTrue(material.metadata['doi'].startswith('https://doi.org/'))
            self.assertTrue(np.all(np.diff(material.omega)>0))
            self.assertTrue(np.all(material.nk>=0))
            with self.assertRaises(ValueError):material.nk.setflags(write=True)
            model,report=material.model()
            self.assertGreater(model.config.epsilon_infinity,0)
            self.assertTrue(all(x>=0 for x in (*model.config.strength,*model.config.resonance,*model.config.damping)))
            self.assertEqual(report['source_sha256'],material.metadata['source_sha256'])
            self.assertEqual(model.config.valid_frequency_range,tuple(report['valid_frequency_range']))
        with self.assertRaises(KeyError):load_material('unknown-material')
        # Corruption is checked on a private copy, never by editing original data.
        import tempfile,shutil
        from pathlib import Path
        from salmonpy import optical_materials as api
        with tempfile.TemporaryDirectory() as tmp:
            copied=Path(tmp)/'data';shutil.copytree(api._root(),copied)
            entry=items[0];raw=copied/entry['raw_file'];raw.write_bytes(raw.read_bytes()+b'corrupt')
            with patch.object(api,'_root',return_value=copied):
                with self.assertRaisesRegex(ValueError,'checksum'):load_material(entry['id'])

    def test_wavelength_units_and_offline_import(self):
        from salmonpy.optical_materials import wavelength_to_frequency
        omega,energy=wavelength_to_frequency([1.0,.5])
        np.testing.assert_allclose(energy,[1.239841984,2.479683968],rtol=1e-8)
        np.testing.assert_allclose(omega,[.0455633525,.091126705],rtol=1e-8)
        with self.assertRaises(ValueError):wavelength_to_frequency([0.])


    def test_isotropic_tensor_layout_and_independence(self):
        from salmonpy.ld import LorentzDrude
        from salmonpy.optical_materials import load_material
        model=LorentzDrude(plasma_frequency=1.,strength=[.1],resonance=[.4],damping=[.03],epsilon_infinity=2.)
        omega=np.array([[.1,.2],[.25,.3]])
        tensor=model.permittivity_tensor(omega)
        self.assertEqual(tensor.shape,(2,2,3,3))
        scalar=model.permittivity(omega)
        for i in range(3):
            for j in range(3):
                np.testing.assert_array_equal(tensor[...,i,j],scalar if i==j else np.zeros_like(scalar))
        with self.assertRaises(ValueError):tensor.setflags(write=True)
        with self.assertRaises(ValueError):model.permittivity_tensor([float('nan')])
        with self.assertRaises(ValueError):model.permittivity_tensor([float('inf')])
        no_loss=LorentzDrude(plasma_frequency=1.,strength=[.1],resonance=[.4],damping=[0.])
        with self.assertRaises(ValueError):no_loss.permittivity_tensor([.4])
        source=load_material('sio2-rodriguez-2016')
        measured=source.epsilon_tensor
        self.assertEqual(measured.shape,(len(source.omega),3,3))
        self.assertEqual(source.metadata['tensor_representation'],'isotropic_diagonal_assumption')
        np.testing.assert_array_equal(measured[:,0,0],source.epsilon)
        with self.assertRaises(ValueError):measured.setflags(write=True)


class WaveplateDiagnosticTests(unittest.TestCase):
    def test_native_phasor_and_stokes(self):
        from salmonpy.analysis import waveplate_diagnostic
        t=np.linspace(0.,20*np.pi,4001)
        field=np.column_stack((np.cos(t),np.sin(t),np.zeros_like(t)))
        result=waveplate_diagnostic(t,field,omega=1.,indices=(1.5,2.),thickness=430.,units='atomic_gaussian')
        np.testing.assert_allclose(result['stokes_normalized'],[0.,0.,-1.],atol=1e-12)
        same=np.column_stack((np.cos(t),np.cos(t),np.zeros_like(t)))
        result=waveplate_diagnostic(t,same,omega=1.,indices=(1.5,1.5),thickness=430.,units='atomic_gaussian')
        np.testing.assert_allclose(result['stokes_normalized'],result['analytic_stokes_normalized'],atol=1e-12)
        bad=same.copy();bad[4,1]=float('nan')
        with self.assertRaises(ValueError):waveplate_diagnostic(t,bad,omega=1.,indices=(1.5,2.),thickness=430.,units='atomic_gaussian')
