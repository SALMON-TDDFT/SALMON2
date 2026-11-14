"""Public C material loader boundary; compiled-plugin tests are explicit tools."""
import importlib
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


class CMaterialBoundary(unittest.TestCase):
    def setUp(self):self.api=importlib.import_module('salmonpy.c_material')

    def test_hash_mismatch_rejected_before_load(self):
        from salmonpy.errors import InputError
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'plugin';path.write_bytes(b'not a shared library')
            with patch('ctypes.CDLL',side_effect=AssertionError('must not load')):
                with self.assertRaises(InputError):self.api.CMaterial(path,'0'*64,b'')

    def test_missing_library_and_nonbytes_config_rejected(self):
        from salmonpy.errors import InputError
        with self.assertRaises((InputError,FileNotFoundError)):
            self.api.CMaterial('/missing/material-plugin','0'*64,b'')
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'plugin';path.write_bytes(b'not a shared library')
            sha=hashlib.sha256(path.read_bytes()).hexdigest()
            with patch('ctypes.CDLL',side_effect=AssertionError('must not load')):
                with self.assertRaises(InputError):self.api.CMaterial(path,sha,'not bytes')

    def test_query_rejects_unsupported_incomplete_or_invalid_metadata(self):
        import ctypes as c
        from types import SimpleNamespace
        from salmonpy.errors import InputError
        for mode in ('version','size','status','callbacks','backend','lag','initial-usage'):
            def query(version,size,out):
                table=c.cast(out,c.POINTER(self.api._API)).contents
                table.abi_version=2 if mode=='version' else 1
                table.struct_size=1 if mode=='size' else c.sizeof(self.api._API)
                table.backend=None if mode=='backend' else b'test/1'
                table.input_lag=-1 if mode=='lag' else 0
                table.uses_initial_A=2 if mode=='initial-usage' else 1
                return 6 if mode=='status' else 0
            with self.subTest(mode=mode),patch('ctypes.CDLL',return_value=SimpleNamespace(salmonpy_material_query=query)):
                with self.assertRaises(InputError):self.api._Library(Path('unused'),None)

    def test_cached_library_path_is_pinned_to_original_hash(self):
        from types import SimpleNamespace
        from salmonpy.errors import InputError
        with tempfile.TemporaryDirectory() as tmp:
            path=(Path(tmp)/'library').resolve();path.write_bytes(b'fixture')
            sha=hashlib.sha256(path.read_bytes()).hexdigest()
            self.api._libraries[path]=SimpleNamespace(sha256='different-build')
            try:
                with patch('ctypes.CDLL',side_effect=AssertionError('must not reload')):
                    with self.assertRaises(InputError):self.api.CMaterial(path,sha,b'')
            finally:self.api._libraries.pop(path)

    def test_failed_abi_load_pins_the_mapped_path(self):
        from types import SimpleNamespace
        from salmonpy.errors import InputError
        with tempfile.TemporaryDirectory() as tmp:
            path=(Path(tmp)/'plugin').resolve();path.write_bytes(b'first fixture')
            sha=hashlib.sha256(path.read_bytes()).hexdigest();identity=(path.stat().st_dev,path.stat().st_ino)
            def query(*args):return 6
            raw=SimpleNamespace(_handle=123456789,salmonpy_material_query=query)
            try:
                with patch('ctypes.CDLL',return_value=raw):
                    with self.assertRaises(InputError):self.api.CMaterial(path,sha,b'')
                path.write_bytes(b'changed fixture');newsha=hashlib.sha256(path.read_bytes()).hexdigest()
                with patch('ctypes.CDLL',side_effect=AssertionError('mapped path must remain pinned')):
                    with self.assertRaises(InputError):self.api.CMaterial(path,newsha,b'')
            finally:
                for name,key in (('_libraries',path),('_pinned_paths',path),('_pinned_files',identity),('_pinned_images',123456789)):
                    getattr(self.api,name,{}).pop(key,None)

    def test_all_libraries_share_the_native_call_lock(self):
        import ctypes as c
        from types import SimpleNamespace
        callbacks=[]
        for name,kind in self.api._API._fields_[5:]:callbacks.append((name,kind(lambda *args:0)))
        def query(version,size,out):
            table=c.cast(out,c.POINTER(self.api._API)).contents
            table.abi_version=1;table.struct_size=c.sizeof(self.api._API);table.backend=b'mock/1'
            table.input_lag=0;table.uses_initial_A=1
            for name,callback in callbacks:setattr(table,name,callback)
            return 0
        raw=SimpleNamespace(salmonpy_material_query=query,refs=callbacks)
        with patch('ctypes.CDLL',return_value=raw):
            one=self.api._Library(Path('first'),None);two=self.api._Library(Path('alias'),None)
        self.assertIs(one.lock,two.lock)

    def test_reused_mapped_image_cannot_claim_a_new_hash(self):
        from types import SimpleNamespace
        from salmonpy.errors import InputError
        with tempfile.TemporaryDirectory() as tmp:
            paths=[(Path(tmp)/name).resolve() for name in ('first','second')]
            for index,path in enumerate(paths):path.write_bytes(bytes([index]))
            raw=SimpleNamespace(_handle=123456790,salmonpy_material_query=lambda *args:6)
            try:
                with patch('ctypes.CDLL',return_value=raw):
                    with self.assertRaises(InputError):self.api.CMaterial(paths[0],hashlib.sha256(paths[0].read_bytes()).hexdigest(),b'')
                    with self.assertRaisesRegex(InputError,'OS loader reused'):
                        self.api.CMaterial(paths[1],hashlib.sha256(paths[1].read_bytes()).hexdigest(),b'')
            finally:
                for path in paths:
                    self.api._libraries.pop(path,None);self.api._pinned_paths.pop(path,None)
                    self.api._pinned_files.pop((path.stat().st_dev,path.stat().st_ino),None)
                self.api._pinned_images.pop(123456790,None)

    def test_invalid_success_outputs_preserve_cache_and_invalidate(self):
        import numpy as np
        from salmonpy.errors import RuntimeFailure
        for clock,current in ((float('nan'),np.zeros(3)),(.1,np.full(3,np.nan)),(.2,np.zeros(3))):
            state=object.__new__(self.api._State)
            state._dt=.1;state._time0=0.;state._time=0.;state._current=np.zeros(3);state._steps=0;state._failed=False
            with self.assertRaises(RuntimeFailure):state._accept(clock,current,1)
            self.assertTrue(state._failed);self.assertEqual(state.time,0.)
            np.testing.assert_array_equal(state.current,np.zeros(3))

    def test_context_cleanup_keeps_primary_exception(self):
        class Owner:
            def close(self):raise RuntimeError('cleanup failure')
        primary=ValueError('original error')
        self.api._finish(Owner(),primary)
        with self.assertRaises(RuntimeError):self.api._finish(Owner(),None)

    def test_import_is_quiet_and_does_not_initialize_mpi(self):
        import subprocess,sys
        result=subprocess.run([sys.executable,'-c',"import sys;import salmonpy.c_material;assert 'mpi4py.MPI' not in sys.modules"],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr);self.assertEqual(result.stdout,'')

if __name__=='__main__':unittest.main()
