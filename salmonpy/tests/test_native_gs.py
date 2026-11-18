"""Public retained GS contracts; real MPI numerical tests live in tools."""
import importlib
import unittest


class NativeGSTests(unittest.TestCase):
    def test_public_backend_and_state_entry(self):
        module=importlib.import_module('salmonpy.native_gs')
        self.assertTrue(callable(module.NativeGSBackend.open))
        self.assertTrue(callable(module.NativeGSSession.run))
        self.assertTrue(callable(module.NativeGSState.array))
        self.assertTrue(callable(module.NativeGSState.close))

    def test_closed_rank_cannot_skip_collective_close(self):
        from types import SimpleNamespace
        from salmonpy.native_gs import NativeGSState
        from salmonpy.errors import InputError
        from unittest.mock import Mock
        comm=SimpleNamespace(allgather=lambda x:[x,(123,False)])
        native=SimpleNamespace(closed=False,_live=lambda:None)
        session=SimpleNamespace(closed=False,comm=comm,_native=native,lib=Mock())
        state=NativeGSState(session,123,None,'A',None);state.closed=True
        with self.assertRaises(InputError):state.close()
        session.lib.sp_gs_close.assert_not_called()

    def test_work_directory_failure_agrees_before_native_call(self):
        from types import SimpleNamespace
        from salmonpy.native_gs import NativeGSSession
        from salmonpy.errors import InputError
        from unittest.mock import Mock
        from pathlib import Path
        session=NativeGSSession(None,SimpleNamespace(_live=lambda:None,context=SimpleNamespace(
            solver_comm=SimpleNamespace(allgather=lambda x:[x,'remote directory missing']))),None)
        fn=Mock()
        previous=Path.cwd()
        with self.assertRaises(InputError):session._call(previous,fn)
        fn.assert_not_called()
        self.assertEqual(Path.cwd(),previous)

    def test_material_identity_cannot_be_relabelled(self):
        from salmonpy.native_gs import NativeGSState
        calculation=object();state=NativeGSState(None,12,None,'film-a',calculation)
        with self.assertRaises(AttributeError):state.material_id='film-b'
        with self.assertRaises(AttributeError):state.calculation=object()
        self.assertEqual(state.material_id,'film-a');self.assertIs(state.calculation,calculation)

    def test_backend_configuration_fixed_without_loading_native_code(self):
        import tempfile,json,hashlib
        from pathlib import Path
        from unittest.mock import patch
        from salmonpy.native_gs import NativeGSBackend
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);library=root/'fixture.bin';library.write_bytes(b'not executable')
            receipt=root/'build.json';receipt.write_text(json.dumps({'schema':'salmonpy.embedded-tddft/1',
                'library':str(library),'sha256':hashlib.sha256(library.read_bytes()).hexdigest()}))
            with patch('salmonpy.native_gs.c.CDLL') as load:
                backend=NativeGSBackend(receipt,root)
                load.assert_not_called()
            for name,value in [('decomposition',(1,2,1,1,1)),('library',Path('/tmp/other')),
                               ('work_root',Path('/tmp/other')),('sha256','different')]:
                with self.assertRaises(AttributeError):setattr(backend,name,value)
