"""Portable GS arrays; synthetic values, no physical ground-state fixtures."""
from pathlib import Path
import json,tempfile,unittest
from unittest.mock import patch
import numpy as np
from salmonpy.errors import InputError


def fixture():
    rng=np.random.default_rng(18)
    # nk=2, nb=5, grid(z,y,x)=(2,3,4); last orbital chunk is partial.
    waves=rng.normal(size=(2,5,2,3,4))+1j*rng.normal(size=(2,5,2,3,4))
    return dict(wavefunctions=waves,density=rng.random((2,3,4)),eigenvalues=rng.normal(size=(2,5)),
                occupations=np.tile([2.,2.,0.,0.,0.],(2,1)),kpoints=[[0,0,0],[.5,0,0]],weights=[.5,.5],
                cell_bohr=np.eye(3)*8,fractional_positions=[[0,0,0],[.25,.25,.25]],atomic_numbers=[14,14],
                electrons=4.,energy_hartree=-3.,xc='pz',pseudopotentials=[dict(atomic_number=14,filename='Si.dat',sha256='a'*64,source='caller-supplied')])


class GSArchiveTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
    def save(self,name='gs',**overrides):
        from salmonpy.gs_archive import save_ground_state
        values=fixture();values.update(overrides)
        return save_ground_state(self.root/name,orbital_block=3,**values)

    def test_component_exists(self):
        import importlib.util
        self.assertIsNotNone(importlib.util.find_spec('salmonpy.gs_archive'),'portable GS archive component is missing')

    def test_roundtrip_blocks_and_padding(self):
        from salmonpy.gs_archive import GroundStateArchive
        values=fixture();archive=self.save()
        self.assertEqual(archive.schema,'salmonpy.gs/2');self.assertEqual(archive.wavefunction_shape,(2,5,2,3,4))
        for k in range(2):
            for lo,hi in ((0,5),(1,4),(3,5)):
                actual=archive.wavefunctions(k,lo,hi)
                self.assertEqual(actual.tobytes(),values['wavefunctions'][k,lo:hi].tobytes())
                self.assertFalse(actual.flags.writeable)
        for name in ('density','eigenvalues','occupations','kpoints','weights','cell_bohr','fractional_positions','atomic_numbers'):
            np.testing.assert_array_equal(archive.array(name),values[name])
        with self.assertRaises(TypeError):archive.metadata['xc']='changed'
        self.assertFalse(any(p.suffix in ('.dat','.upf') for p in archive.directory.rglob('*')))
        opened=GroundStateArchive.open(archive.directory);self.assertEqual(opened.metadata['electrons'],4.)

    def test_bit_patterns_and_lazy_selected_chunks(self):
        from salmonpy.gs_archive import GroundStateArchive
        waves=np.empty((2,5,2,3,4),dtype=complex)
        bits=np.resize(np.array([0.,-0.,np.nextafter(0.,1.),-np.nextafter(0.,1.),1e200,-1e200]),waves.shape)
        waves.real=bits;waves.imag=-bits
        a=self.save(wavefunctions=waves)
        original=np.empty
        def bounded(shape,*args,**kwargs):
            if isinstance(shape,tuple) and len(shape)>3:raise AssertionError('full waves allocated on open')
            return original(shape,*args,**kwargs)
        with patch('salmonpy.gs_archive.np.empty',side_effect=bounded):b=GroundStateArchive.open(a.directory)
        self.assertEqual(b.wavefunctions(0).tobytes(),waves[0].tobytes())
        # A stale unrequested k chunk does not force an eager read of every k.
        (a.directory/'wavefunction_real/c/1/0/0/0/0').write_bytes(b'corrupt')
        self.assertEqual(b.wavefunctions(0,1,4).tobytes(),waves[0,1:4].tobytes())
        with self.assertRaises(InputError):b.wavefunctions(1,0,1)

    def test_nonfinite_and_inconsistent_inputs_create_no_archive(self):
        for kwargs in ({'electrons':5.},{'weights':[0.,0.]},{'occupations':np.full((2,5),3.)},
                       {'cell_bohr':np.zeros((3,3))},{'density':np.full((2,3,4),np.nan)},
                       {'wavefunctions':np.full((2,5,2,3,4),np.inf+0j)},
                       {'atomic_numbers':[True,14]},{'xc':''},
                       {'pseudopotentials':[dict(atomic_number=14,filename='../Si.dat',sha256='a'*64,source='caller')]}):
            with self.subTest(kwargs=list(kwargs)),self.assertRaises(ValueError):self.save('bad',**kwargs)
            self.assertFalse((self.root/'bad').exists())

    def test_corruption_and_extra_entries_rejected(self):
        from salmonpy.gs_archive import GroundStateArchive
        a=self.save();chunk=a.directory/'wavefunction_real/c/0/0/0/0/0';raw=chunk.read_bytes()
        chunk.write_bytes(b'x'+raw[1:])
        with self.assertRaises(InputError):GroundStateArchive.open(a.directory)
        with self.assertRaises(InputError):a.wavefunctions(0,0,1)
        chunk.write_bytes(raw);extra=a.directory/'extra';extra.mkdir()
        with self.assertRaises(InputError):GroundStateArchive.open(a.directory)
        extra.rmdir();chunk.unlink();chunk.symlink_to(self.root/'outside')
        with self.assertRaises(InputError):GroundStateArchive.open(a.directory)

    def test_budget_overwrite_and_failed_write(self):
        from salmonpy.gs_archive import GroundStateArchive
        a=self.save()
        with self.assertRaises(FileExistsError):self.save()
        with self.assertRaises(InputError):GroundStateArchive.open(a.directory,max_bytes=1)
        with self.assertRaises(InputError):a.wavefunctions(0,0,5,max_bytes=1)
        with patch('salmonpy.gs_archive._write',side_effect=OSError('full disk')):
            with self.assertRaises(OSError):self.save('failed')
        self.assertFalse((self.root/'failed').exists())

    def test_extreme_dimensions_and_changed_chunk_size_fail_before_allocation(self):
        from salmonpy.gs_archive import GroundStateArchive
        a=self.save();m=a.directory/'manifest.json';raw=m.read_bytes();data=json.loads(raw)
        data['metadata']['n_kpoints']=10**30
        m.write_text(json.dumps(data))
        with patch('salmonpy.gs_archive._chunks',side_effect=AssertionError('must not enumerate impossible file count')):
            with self.assertRaises(InputError):GroundStateArchive.open(a.directory)
        m.write_bytes(raw);chunk=a.directory/'wavefunction_real/c/0/0/0/0/0'
        with chunk.open('ab') as f:f.write(b'extra')
        with patch.object(Path,'read_bytes',side_effect=AssertionError('must reject size before reading')):
            with self.assertRaises(InputError):a.wavefunctions(0,0,1)

    def test_unknown_directory_is_rejected_before_descending(self):
        from salmonpy.gs_archive import GroundStateArchive
        a=self.save()
        def entries(root,pattern):
            yield root/'unexpected'
            raise AssertionError('reader descended into an unknown tree')
        with patch.object(Path,'rglob',entries):
            with self.assertRaises(InputError):GroundStateArchive.open(a.directory)

    def test_malformed_manifest_and_root_symlink(self):
        from salmonpy.gs_archive import GroundStateArchive
        a=self.save();m=a.directory/'manifest.json';raw=m.read_bytes();data=json.loads(raw)
        data['files']['../outside']='a'*64;m.write_text(json.dumps(data))
        with self.assertRaises(InputError):GroundStateArchive.open(a.directory)
        m.write_bytes(raw);link=self.root/'link';link.symlink_to(a.directory,target_is_directory=True)
        with self.assertRaises(InputError):GroundStateArchive.open(link)

if __name__=='__main__':unittest.main()
