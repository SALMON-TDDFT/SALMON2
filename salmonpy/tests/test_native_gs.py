"""Public retained GS contracts; real MPI numerical tests live in tools."""
import importlib.util
from pathlib import Path
from unittest.mock import patch
import unittest


class NativeGSTests(unittest.TestCase):

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


    def test_imported_force_is_explicitly_unavailable(self):
        from salmonpy.native_gs import NativeGSState
        from salmonpy.errors import InputError
        from types import SimpleNamespace
        from unittest.mock import Mock
        session=SimpleNamespace(_live=lambda:None,comm=SimpleNamespace(allgather=lambda v:[v]),lib=Mock())
        state=NativeGSState(session,1,None,'imported',None)
        with self.assertRaisesRegex(InputError,'forces were not computed'):state.array('forces')
        session.lib.sp_gs_array_info.assert_not_called()


class NativeGSLinkTests(unittest.TestCase):
    def tool(self):
        root=Path(__file__).resolve().parents[2]
        import salmonpy
        path=Path(salmonpy.__file__).parent/'_runtime_source/tools/build_embedded_tddft.py'
        if not path.exists():path=root/'tools/build_embedded_tddft.py'
        spec=importlib.util.spec_from_file_location('gs_link_test',path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        return module

    def test_darwin_mpi_wrapper_does_not_interpose_blas_or_duplicate_sentinels(self):
        tool=self.tool()
        def show(command,**kwargs):
            return 'gfortran\n' if command[-1]=='--showme:command' else '-L/mpi -Wl,-flat_namespace -lmpi_ignore -lmpi\n'
        with patch.object(tool.sys,'platform','darwin'),patch.object(tool.subprocess,'check_output',side_effect=show):
            library=tool.mpi_link_command(['mpifort','object.o','-lopenblas'],shared=True)
            host=tool.mpi_link_command(['mpifort','probe.o','-lnative'])
        self.assertEqual(library[0],'gfortran')
        self.assertNotIn('-Wl,-flat_namespace',library)
        self.assertIn('-Wl,-commons,use_dylibs',library)
        self.assertLess(library.index('-lmpi'),library.index('-lmpi_ignore'))
        self.assertIn('-dynamiclib',library);self.assertNotIn('-dynamiclib',host)

    def test_linux_wrapper_and_flags_stay_intact(self):
        tool=self.tool()
        with patch.object(tool.sys,'platform','linux'),patch.object(tool.subprocess,'check_output') as probe:
            original=['mpifort','object.o','-lopenblas']
            self.assertEqual(tool.mpi_link_command(original,shared=True),['mpifort','-shared','object.o','-lopenblas'])
            self.assertEqual(tool.mpi_link_command(original),original)
            probe.assert_not_called()
