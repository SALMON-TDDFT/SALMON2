"""Bridge contracts with synthetic converter/protocol fixtures, no physical GS."""
from pathlib import Path
import hashlib,importlib.util,json,os,sys,subprocess,unittest
from unittest.mock import patch
import numpy as np
import test_tddft_snapshots as fixture
from salmonpy.errors import InputError,SnapshotError

CONVERTER='''#!/usr/bin/env python3
from pathlib import Path
import struct
p=Path('.')
r=(p/'portable_real.bin').read_bytes();i=(p/'portable_imag.bin').read_bytes()
out=b''.join(r[j:j+8]+i[j:j+8] for j in range(0,len(r),8))
(p/'data_for_restart/wfn.bin').write_bytes(out)
(p/'data_for_restart/info.bin').write_bytes(b'synthetic converter info')
(p/'data_for_restart/occupation.bin').write_bytes((p/'portable_occupations.bin').read_bytes())
'''

class GSNativeTests(unittest.TestCase):
    def setUp(self):
        fixture.TDDFTSnapshots.setUp(self)
        from salmonpy.gs_native import GammaRestartConverter
        from salmonpy.dft import Pseudopotential
        from salmonpy.gs_archive import save_ground_state
        self.potential=Pseudopotential(2,0,'_rps.dat',b'synthetic potential fixture, not physics')
        root=Path(self.tmp.name);converter=(root/'converter').resolve();converter.write_text(CONVERTER);converter.chmod(0o755)
        record=dict(schema='salmonpy.gs-native-build/1',status='built-audited-and-probed',native_byte_parity=True,
                    converter=str(converter),sha256=hashlib.sha256(converter.read_bytes()).hexdigest(),
                    target_executable=str(self.model.executable.path),target_sha256=self.model.executable.sha256)
        rp=root/'converter-build.json';rp.write_text(json.dumps(record));self.converter=GammaRestartConverter.from_build(rp)
        original=subprocess.run
        def run(args,**kw):
            if args==[str(converter)]:args=[sys.executable,str(converter)]
            return original(args,**kw)
        mock=patch('salmonpy.gs_native.subprocess.run',side_effect=run);mock.start();self.addCleanup(mock.stop)
        self.values=dict(wavefunctions=np.ones((1,2,8,8,8),complex)/8,density=np.full((8,8,8),2/64),
                         eigenvalues=[[-.5,.2]],occupations=[[2.,0.]],kpoints=[[0.,0.,0.]],weights=[1.],
                         cell_bohr=np.eye(3)*4,fractional_positions=[[0.,0.,0.]],atomic_numbers=[2],
                         electrons=2.,energy_hartree=-1.,xc='pz',pseudopotentials=[dict(atomic_number=2,filename='He.dat',
                         sha256=self.potential.sha256,source='synthetic fixture')])
        self.archive=save_ground_state(root/'archive',orbital_block=1,**self.values)

    def initial(self):
        a=self.archive.tddft_initial_state(converter=self.converter,executable=self.model.executable,
                                           pseudopotentials=[self.potential],max_steps=20)
        self.addCleanup(a.close);return a

    def test_retained_archive_loader_converts_once_and_keeps_unknown_provenance(self):
        from types import SimpleNamespace
        from unittest.mock import Mock
        from salmonpy.native_gs import NativeGSSession,NativeGSState
        comm=SimpleNamespace(size=1,rank=0,allgather=lambda v:[v],bcast=lambda v,root:v)
        backend=SimpleNamespace(sha256='b'*64,_io_profile={'source_sha256':None,'compiler_sha256':None,'record_flags':()})
        session=NativeGSSession(backend,SimpleNamespace(_live=lambda:None,context=SimpleNamespace(solver_comm=comm)),None)
        state=NativeGSState(session,1,Path(self.tmp.name)/'state','film-a',None)
        session.load=Mock(return_value=state)
        destination=Path(self.tmp.name)/'retained-import'
        with self.assertWarnsRegex(RuntimeWarning,'compiler'):
            loaded=session.load_archive(self.archive,destination,converter=self.converter,
                executable=self.model.executable,pseudopotentials=[self.potential],material_id='film-a')
        self.assertIs(loaded,state)
        session.load.assert_called_once_with(destination.resolve(),'rt.inp',material_id='film-a')
        self.assertIsNone(loaded.calculation)
        self.assertTrue(loaded.archive_origin['explicit_gamma_complex_table'])
        self.assertIn("file_kw='kpoints.dat'",(destination/'rt.inp').read_text())
        self.assertEqual((destination/'kpoints.dat').read_text(),'1\n1 0 0 0 1\n')
        self.assertEqual(loaded.archive_origin['archive_manifest_sha256'],
                         hashlib.sha256((self.archive.directory/'manifest.json').read_bytes()).hexdigest())
        with self.assertRaises(TypeError):loaded.archive_origin['material_id']='other'
        with self.assertRaises(FileExistsError):
            session.load_archive(self.archive,destination,converter=self.converter,
                executable=self.model.executable,pseudopotentials=[self.potential])
        self.assertEqual(session.load.call_count,1)
        backend._io_profile['record_flags']=('-frecord-marker=8',)
        with self.assertRaisesRegex(InputError,'record format'):
            session.load_archive(self.archive,Path(self.tmp.name)/'bad-profile',converter=self.converter,
                executable=self.model.executable,pseudopotentials=[self.potential])
        self.assertEqual(session.load.call_count,1)
        self.assertFalse((Path(self.tmp.name)/'bad-profile').exists())
        from salmonpy.gs_archive import save_ground_state
        provenance=dict(schema='salmonpy.retained-gs-export/1',material_id='archive-film',
            state_identity='synthetic fixture only',target_sha256='a'*64,converged=False,iterations=0,
            density_source='reconstructed_from_native_orbitals',local_channels=[dict(atomic_number=2,local_channel=0)])
        inherited=save_ground_state(Path(self.tmp.name)/'with-origin',provenance=provenance,**self.values)
        second=NativeGSState(session,2,Path(self.tmp.name)/'second','archive-film',None)
        session.load.return_value=second;backend._io_profile['record_flags']=()
        with self.assertWarns(RuntimeWarning):
            other=session.load_archive(inherited,Path(self.tmp.name)/'inherited-import',converter=self.converter,
                executable=self.model.executable,pseudopotentials=[self.potential])
        self.assertEqual(other.material_id,'archive-film')
        self.assertEqual(other.archive_origin['archive_material_id'],'archive-film')
        self.assertFalse(other.archive_origin['archive_converged'])
        self.assertIsNone(other.archive_origin['scf_converged'])
        self.assertEqual(state.archive_origin['material_id'],'film-a')

    def test_imported_archive_reexport_labels_carried_values_and_reads_native_occupations(self):
        from types import SimpleNamespace
        from unittest.mock import Mock,PropertyMock
        from salmonpy.native_gs import NativeGSSession,NativeGSState,NativeGSInfo,NativeGSArray
        comm=SimpleNamespace(size=1,rank=0,allgather=lambda v:[v],bcast=lambda v,root:v)
        backend=SimpleNamespace(sha256='b'*64,decomposition=(1,1,1,1,1),
            _io_profile={'source_sha256':None,'compiler_sha256':None,'record_flags':()})
        session=NativeGSSession(backend,SimpleNamespace(_live=lambda:None,context=SimpleNamespace(solver_comm=comm)),None)
        state=NativeGSState(session,1,Path(self.tmp.name)/'state','film-a',None)
        state._identity='synthetic imported state';session.load=Mock(return_value=state)
        from salmonpy.gs_archive import save_ground_state
        provenance=dict(schema='salmonpy.retained-gs-export/1',material_id='original-id',
            state_identity='synthetic metadata only',target_sha256='a'*64,converged=False,iterations=0,
            density_source='reconstructed_from_native_orbitals',local_channels=[dict(atomic_number=2,local_channel=0)])
        self.archive=save_ground_state(Path(self.tmp.name)/'source-with-claim',provenance=provenance,**self.values)
        with self.assertWarns(RuntimeWarning):
            session.load_archive(self.archive,Path(self.tmp.name)/'import',converter=self.converter,
                executable=self.model.executable,pseudopotentials=[self.potential],material_id='film-a')
        waves=self.values['wavefunctions']
        raw=waves.transpose(4,3,2,1,0)[:,:,:,None,:,:,None]
        # Native occupations deliberately differ from the retained source values.
        occ=np.array([1.5,.5]).reshape(2,1,1)
        values={'density':NativeGSArray(np.zeros((8,8,8)),(1,1,1),'density','bohr^-3'),
                'orbitals':NativeGSArray(raw,(1,1,1,1,1,1,1),'orbitals','native'),
                'occupations':NativeGSArray(occ,(1,1,1),'occupations','1')}
        state.array=Mock(side_effect=lambda name:values[name])
        import shutil
        shutil.rmtree(self.archive.directory)  # retained small metadata must suffice
        info=NativeGSInfo(987.,0.,0.,0,None,1,2)
        with patch.object(NativeGSState,'info',new_callable=PropertyMock,return_value=info):
            with self.assertWarnsRegex(RuntimeWarning,'unknown'):
                archive=state.export_archive(Path(self.tmp.name)/'reexport')
        np.testing.assert_array_equal(archive.wavefunctions(0),waves[0])
        np.testing.assert_array_equal(archive.array('occupations'),[[1.5,.5]])
        np.testing.assert_array_equal(archive.array('eigenvalues'),self.values['eigenvalues'])
        self.assertEqual(archive.metadata['energy_hartree'],-1.)
        provenance=archive.metadata['provenance']
        self.assertEqual(provenance['schema'],'salmonpy.retained-gs-reexport/1')
        self.assertEqual(provenance['energy_source'],'source_archive')
        self.assertEqual(provenance['eigenvalue_source'],'source_archive')
        self.assertEqual(provenance['sampling_source'],'source_archive')
        self.assertIsNone(provenance['converged'])
        self.assertFalse(provenance['inherited_convergence'])
        self.assertIsNone(state.calculation)
        self.assertNotIn('eigenvalues',[call.args[0] for call in state.array.call_args_list])
        self.assertEqual(state.archive_origin['material_id'],'film-a')

    def test_bridge_entry_and_owned_initial_replay(self):
        self.assertIsNotNone(importlib.util.find_spec('salmonpy.gs_native'))
        initial=self.initial()
        self.assertNotIn('portable_real.bin',dict(initial.manifest))
        with initial.new_state(dt=.02) as a,initial.new_state(dt=.02) as b:
            a.step([.1,0,0]);token=a.snapshot();expected=a.step([.2,0,0]);b.restore(token)
            self.assertEqual(b.step([.2,0,0]).tobytes(),expected.tobytes())
            initial.close();a.restore(token);self.assertEqual(a.step([.2,0,0]).tobytes(),expected.tobytes())

    def test_bad_potential_target_and_gamma_refused_before_conversion(self):
        from dataclasses import replace
        from salmonpy.gs_archive import save_ground_state
        from salmonpy.dft import Pseudopotential
        bad=Pseudopotential(2,0,'_rps.dat',b'wrong potential')
        dest=Path(self.tmp.name)/'bad'
        with patch('salmonpy.gs_native.subprocess.run',side_effect=AssertionError('must not run')):
            with self.assertRaises(InputError):self.converter.write_bundle(self.archive,dest,executable=self.model.executable,pseudopotentials=[bad])
            with self.assertRaises(InputError):replace(self.converter,target_sha256='0'*64).write_bundle(self.archive,dest,executable=self.model.executable,pseudopotentials=[self.potential])
            values=dict(self.values,kpoints=[[.1,0,0]])
            other=save_ground_state(Path(self.tmp.name)/'other',**values)
            with self.assertRaises(InputError):self.converter.write_bundle(other,dest,executable=self.model.executable,pseudopotentials=[self.potential])
        self.assertFalse(dest.exists())

    def test_failed_converter_cleanup_and_overwrite_refusal(self):
        dest=Path(self.tmp.name)/'failed'
        with patch('salmonpy.gs_native.subprocess.run',return_value=subprocess.CompletedProcess([],1)):
            with self.assertRaises(RuntimeError):self.converter.write_bundle(self.archive,dest,executable=self.model.executable,pseudopotentials=[self.potential])
        self.assertFalse(dest.exists())
        dest.mkdir()
        with self.assertRaises(FileExistsError):self.converter.write_bundle(self.archive,dest,executable=self.model.executable,pseudopotentials=[self.potential])
        self.assertTrue(dest.exists())

    def test_environment_change_refuses_spawn_and_preserves_restore(self):
        initial=self.initial()
        with initial.new_state(dt=.02) as s:
            token=s.snapshot();pid=s._proc.pid
            value='swap' if os.environ.get('GFORTRAN_CONVERT_UNIT','')!='swap' else 'native'
            with patch.dict(os.environ,{'GFORTRAN_CONVERT_UNIT':value}):
                with self.assertRaises(InputError):initial.new_state(dt=.02)
                with self.assertRaises(SnapshotError):s.restore(token)
            self.assertEqual(s._proc.pid,pid);s.step([.1,0,0])

    def test_unsupported_record_environment_creates_no_destination(self):
        dest=Path(self.tmp.name)/'bad-environment'
        with patch.dict(os.environ,{'GFORTRAN_CONVERT_UNIT':'swap;native:96'}):
            with self.assertRaises(InputError):self.converter.write_bundle(self.archive,dest,executable=self.model.executable,pseudopotentials=[self.potential])
        self.assertFalse(dest.exists())

    def test_spawn_passes_the_validated_environment_snapshot(self):
        import salmonpy.tddft as tddft
        initial=self.initial();expected=dict(initial._record_environment)
        copy=tddft._copy_bundle;popen=tddft.subprocess.Popen
        def changing(source,dest):
            value=copy(source,dest);os.environ['GFORTRAN_CONVERT_UNIT']='swap';return value
        def checked(args,**kwargs):
            self.assertEqual(kwargs['env'].get('GFORTRAN_CONVERT_UNIT',''),expected['GFORTRAN_CONVERT_UNIT'])
            return popen(args,**kwargs)
        with patch.dict(os.environ,dict(os.environ)),patch.object(tddft,'_copy_bundle',side_effect=changing),patch.object(tddft.subprocess,'Popen',side_effect=checked):
            with initial.new_state(dt=.02) as s:s.step([.1,0,0])

    def test_native_wavefunction_hashing_streams_without_full_read(self):
        original=Path.read_bytes
        def bounded(path):
            if path.name=='wfn.bin':raise AssertionError('full native wavefunction read')
            return original(path)
        with patch.object(Path,'read_bytes',bounded):
            initial=self.initial()
            with initial.new_state(dt=.02) as s:s.step([.1,0,0])

    def test_converter_passes_the_captured_environment(self):
        import salmonpy.gs_native as native
        run=native.subprocess.run;expected=os.environ.get('GFORTRAN_CONVERT_UNIT','')
        def changing(args,**kwargs):
            os.environ['GFORTRAN_CONVERT_UNIT']='swap'
            self.assertEqual(kwargs['env'].get('GFORTRAN_CONVERT_UNIT',''),expected)
            return run(args,**kwargs)
        with patch.dict(os.environ,dict(os.environ)),patch.object(native.subprocess,'run',side_effect=changing):
            self.converter.write_bundle(self.archive,Path(self.tmp.name)/'environment-snapshot',executable=self.model.executable,pseudopotentials=[self.potential])

    def test_checkpoint_cannot_modify_portable_source_archive(self):
        initial=self.initial()
        with initial.new_state(dt=.02) as state:
            dest=self.archive.directory/'checkpoint'
            with self.assertRaises(SnapshotError):state.save_checkpoint(dest)
            self.assertFalse(dest.exists())

    def test_receipt_tamper_and_converter_hash_rejected(self):
        from salmonpy.tddft import TDDFTInitialState
        dest=Path(self.tmp.name)/'native-bundle'
        self.converter.write_bundle(self.archive,dest,executable=self.model.executable,pseudopotentials=[self.potential])
        (dest/'data_for_restart/wfn.bin').write_bytes(b'changed')
        with self.assertRaises(InputError):TDDFTInitialState.from_salmon(dest,'rt.inp',executable=self.model.executable)
        self.converter.path.write_text('changed')
        with self.assertRaises(InputError):self.converter.write_bundle(self.archive,Path(self.tmp.name)/'changed',executable=self.model.executable,pseudopotentials=[self.potential])

if __name__=='__main__':unittest.main()
