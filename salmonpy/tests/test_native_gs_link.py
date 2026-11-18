"""Darwin host/library namespace regressions; no compiler or MPI startup."""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch


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
