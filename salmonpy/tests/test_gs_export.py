"""GS/1 export contracts with synthetic converter fixtures, not physical GS."""
from pathlib import Path
import hashlib,json,subprocess,sys,unittest
from unittest.mock import patch
import numpy as np
import test_gs_native as fixture
from salmonpy.dft import GroundStateResult,DFT,SCFSettings
from salmonpy.gs_native import GammaRestartConverter
from salmonpy.errors import InputError

EXPORTER='''#!/usr/bin/env python3
from pathlib import Path
p=Path('.')
if (p/'native_input/info.bin').read_bytes()!=b'complex':raise SystemExit(1)
w=(p/'native_input/wfn.bin').read_bytes()
(p/'portable_real.bin').write_bytes(b''.join(w[i:i+8] for i in range(0,len(w),16)))
(p/'portable_imag.bin').write_bytes(b''.join(w[i+8:i+16] for i in range(0,len(w),16)))
(p/'portable_occupations.bin').write_bytes((p/'native_input/occupation.bin').read_bytes())
(p/'native_kind.txt').write_text('complex')
'''

class GSExportTests(unittest.TestCase):
    def setUp(self):
        fixture.GSNativeTests.setUp(self)
        root=Path(self.tmp.name);self.exporter=(root/'exporter').resolve();self.exporter.write_text(EXPORTER);self.exporter.chmod(0o755)
        record=json.loads(self.converter.build_record.read_text());record.update(exporter=str(self.exporter),
            exporter_sha256=hashlib.sha256(self.exporter.read_bytes()).hexdigest(),native_export_parity=True)
        self.converter.build_record.write_text(json.dumps(record));self.converter=GammaRestartConverter.from_build(self.converter.build_record)
        original=subprocess.run
        def run(args,**kw):
            if args==[str(self.exporter)]:args=[sys.executable,str(self.exporter)]
            return original(args,**kw)
        wrapper=patch('salmonpy._gs_archive.subprocess.run',side_effect=run);wrapper.start();self.addCleanup(wrapper.stop)
        self.source=root/'source-gs';native=self.source/'native';restart=native/'data_for_restart';restart.mkdir(parents=True)
        self.dft=DFT([2],np.eye(3)*4,[[0,0,0]],[self.potential],electrons=2,states=2,grid=(8,8,8),kmesh=(1,1,1))
        (native/'gs.inp').write_text(self.dft.input_text())
        (native/'pseudo1_rps.dat').write_bytes(self.potential.data)
        (native/'gs_status.dat').write_text('SALMONPY_GS_V1\n1 1\n1e-10 1e-9 -1.0\nrho_dne\n')
        (native/'gs_eigen.data').write_text('#esp: single-particle energies (eigen energies)\n#occ: occupation numbers, io: orbital index\n# 1:io, 2:esp[a.u.], 3:occ\nk= 1,  spin= 1\n1 -.5 2\n2 .2 0\n')
        (restart/'info.bin').write_bytes(b'complex')
        (restart/'occupation.bin').write_bytes(np.array([2.,0.],dtype='<f8').tobytes())
        self.waves=np.asarray(self.values['wavefunctions'],dtype='<c16');(restart/'wfn.bin').write_bytes(self.waves.tobytes())
        settings=SCFSettings()
        self.record=dict(self.dft._restart_signature(),schema='salmonpy.gs/1',status='complete',restart_saved=True,
            fractional_positions=self.dft.fractional_positions.tolist(),executable_sha256=self.model.executable.sha256,
            energy_hartree=-1.,iterations=1,tolerance=settings.tolerance,max_iterations=settings.max_iterations,
            input_sha256=hashlib.sha256((native/'gs.inp').read_bytes()).hexdigest())
        self.rehash()
    def rehash(self):
        self.record['files']={p.relative_to(self.source).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in (self.source/'native').rglob('*') if p.is_file()}
        (self.source/'ground-state.json').write_text(json.dumps(self.record))
        self.result=GroundStateResult.open(self.source)
    def export(self,name='export'):
        return self.result.export_archive(Path(self.tmp.name)/name,converter=self.converter,orbital_block=1)

    def test_bits_density_provenance_and_source_unchanged(self):
        original=(self.source/'ground-state.json').read_bytes();a=self.export()
        self.assertEqual(a.wavefunctions(0).tobytes(),self.waves[0].tobytes())
        self.assertEqual(a.array('occupations').tobytes(),np.array([[2.,0.]]).tobytes())
        np.testing.assert_array_equal(a.array('eigenvalues'),[[-.5,.2]])
        np.testing.assert_array_equal(a.array('density'),np.full((8,8,8),2/64))
        self.assertEqual(a.metadata['provenance']['density_source'],'reconstructed_from_native_orbitals')
        self.assertEqual(a.metadata['provenance']['local_channels'][0]['local_channel'],0)
        self.assertEqual((self.source/'ground-state.json').read_bytes(),original)
        GroundStateResult.open(self.source)
        from dataclasses import replace
        self.dft=replace(self.dft,xc='tbmbj',cval=1.2)
        native=self.source/'native'
        (native/'gs.inp').write_text(self.dft.input_text())
        self.record.update(xc='tbmbj',cval=1.2,input_sha256=hashlib.sha256((native/'gs.inp').read_bytes()).hexdigest())
        self.rehash();tb=self.export('tb-export')
        self.assertEqual(tb.metadata['xc_parameters']['cval'],1.2)
        converted=self.converter.write_bundle(tb,Path(self.tmp.name)/'tb-roundtrip',
            executable=self.model.executable,pseudopotentials=[self.potential])
        self.assertIn('cval=1.19999999999999996e+00',(converted/'rt.inp').read_text())

    def test_missing_eigen_units_and_occ_mismatch_no_output(self):
        p=self.source/'native/gs_eigen.data';original=p.read_text()
        for text in (original.replace('a.u.','eV'),original.replace('1 -.5 2','1 -.5 1')):
            p.write_text(text);self.rehash()
            with self.assertRaises(InputError):self.export('bad')
            self.assertFalse((Path(self.tmp.name)/'bad').exists())
        p.unlink();self.rehash()
        with self.assertRaises(InputError):self.export('missing')

    def test_model_mismatch_and_target_rejected(self):
        self.record['grid']=[9,8,8];self.rehash()
        with self.assertRaises(InputError):self.export('wrong-grid')
        self.record['grid']=[8,8,8];self.record['executable_sha256']='0'*64;self.rehash()
        with self.assertRaises(InputError):self.export('wrong-target')

    def test_failure_cleanup_source_protection_and_overwrite(self):
        with self.assertRaises(InputError):self.result.export_archive(self.source/'archive',converter=self.converter)
        (self.source/'native/data_for_restart/info.bin').write_bytes(b'bad');self.rehash()
        with self.assertRaises(RuntimeError):self.export('failed')
        self.assertFalse((Path(self.tmp.name)/'failed').exists())
        dest=Path(self.tmp.name)/'existing';dest.mkdir()
        with self.assertRaises(FileExistsError):self.result.export_archive(dest,converter=self.converter)
        self.assertTrue(dest.exists())

    def test_export_never_reads_entire_native_wavefunction(self):
        original=Path.read_bytes
        def bounded(path):
            if path.name=='wfn.bin':raise AssertionError('full native wavefunction allocation')
            return original(path)
        with patch.object(Path,'read_bytes',bounded):a=self.export('streamed')
        self.assertEqual(a.wavefunctions(0).tobytes(),self.waves[0].tobytes())

    def test_old_import_only_converter_and_bad_provenance_rejected(self):
        from dataclasses import replace
        from salmonpy.gs_archive import save_ground_state
        with self.assertRaises(RuntimeError):self.result.export_archive(Path(self.tmp.name)/'old',converter=replace(self.converter,exporter=None,exporter_sha256=None))
        values=dict(self.values,provenance={'schema':'bad'})
        with self.assertRaises(InputError):save_ground_state(Path(self.tmp.name)/'bad-provenance',**values)
        self.assertFalse((Path(self.tmp.name)/'bad-provenance').exists())

    def test_local_channels_enforced_on_reimport(self):
        from dataclasses import replace
        a=self.export()
        with self.assertRaises(InputError):self.converter.write_bundle(a,Path(self.tmp.name)/'wrong-channel',executable=self.model.executable,pseudopotentials=[replace(self.potential,local_channel=1)])

if __name__=='__main__':unittest.main()
