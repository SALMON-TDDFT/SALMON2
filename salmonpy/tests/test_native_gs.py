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

    def test_retained_gamma_archive_roundtrip_and_source_independence(self):
        import tempfile
        import numpy as np
        from types import SimpleNamespace
        from unittest.mock import Mock
        from salmonpy.native_gs import NativeGSState, NativeGSArray, NativeGSInfo
        from salmonpy.dft import DFT, Pseudopotential
        from salmonpy.gs_archive import GroundStateArchive
        model=DFT((2,),np.diag([8.,9.,10.]),[[0.,0.,0.]],
                  (Pseudopotential(2,0,'_rps.dat',b'synthetic only'),),2,2,(8,10,12),(1,1,1))
        comm=SimpleNamespace(size=1,rank=0,allgather=lambda v:[v])
        session=SimpleNamespace(_live=lambda:None,comm=comm,backend=SimpleNamespace(
            decomposition=(1,1,1,1,1),sha256='a'*64))
        state=NativeGSState(session,1,None,'film-a',model);state._identity='synthetic-owner'
        shape=(8,10,12);real=np.arange(8*10*12*2,dtype=float).reshape((*shape,2))/1000
        waves=real+1j*(real+.125)
        padded=np.full((10,14,18,1,2,1,1),complex(999,999))
        padded[1:9,2:12,3:15,0,:,0,0]=waves
        arrays={'density':NativeGSArray(np.ones(shape),(1,1,1),'density','bohr^-3'),
                'orbitals':NativeGSArray(padded,(0,-1,-2,1,1,1,1),'orbitals','native'),
                'occupations':NativeGSArray(np.array([2.,0.]).reshape(2,1,1),(1,1,1),'occupations','1'),
                'eigenvalues':NativeGSArray(np.array([-1.,.5]).reshape(2,1,1),(1,1,1),'eigenvalues','hartree')}
        state.array=Mock(side_effect=lambda q:arrays[q])
        info=NativeGSInfo(-2.,1e-5,1e-9,7,False,1,2)
        with tempfile.TemporaryDirectory() as tmp, patch.object(NativeGSState,'info',new_callable=unittest.mock.PropertyMock,return_value=info):
            target=Path(tmp)/'archive'
            with self.assertWarnsRegex(RuntimeWarning,'not converged'):
                archive=state.export_archive(target,orbital_block=1)
            expected=waves.transpose(3,2,1,0)[None].copy()
            np.testing.assert_array_equal(archive.wavefunctions(0),expected[0])
            self.assertEqual(archive.metadata['provenance']['material_id'],'film-a')
            self.assertFalse(archive.metadata['provenance']['converged'])
            np.testing.assert_allclose(archive.array('density'),2*abs(expected[0,0])**2,rtol=1e-15,atol=0)
            self.assertEqual(padded[0,0,0,0,0,0,0],999+999j)
            with self.assertRaises(FileExistsError):state.export_archive(target)
            old=arrays['orbitals']
            arrays['orbitals']=NativeGSArray(padded,(5,5,5,1,1,1,1),'orbitals','native')
            with self.assertRaisesRegex(ValueError,'halos'):
                state.export_archive(Path(tmp)/'invalid-halo')
            self.assertFalse((Path(tmp)/'invalid-halo').exists())
            arrays['orbitals']=old
            # Real native orbitals are returned by the public reader as complex+0j.
            padded.imag[:]=0
            with self.assertWarns(RuntimeWarning):
                real_archive=state.export_archive(Path(tmp)/'real-archive')
            np.testing.assert_array_equal(real_archive.wavefunctions(0),expected[0].real.astype(complex))
            padded[:]=0;state.closed=True
            reopened=GroundStateArchive.open(target)
            np.testing.assert_array_equal(reopened.wavefunctions(0),expected[0])
            with self.assertRaises(ValueError):reopened.wavefunctions(0).flat[0]=0

    def test_retained_regular_multik_archive_preserves_order_weights_and_independence(self):
        import tempfile
        import numpy as np
        from types import SimpleNamespace
        from unittest.mock import Mock
        from salmonpy.native_gs import NativeGSState,NativeGSArray,NativeGSInfo
        from salmonpy.dft import DFT,Pseudopotential
        from salmonpy.gs_archive import GroundStateArchive
        comm=SimpleNamespace(size=1,rank=0,allgather=lambda v:[v])
        session=SimpleNamespace(_live=lambda:None,comm=comm,backend=SimpleNamespace(
            decomposition=(1,1,1,1,1),sha256='b'*64))
        # Distinct cell/grid/k/XC model values; no solver or real GS data involved.
        with tempfile.TemporaryDirectory() as tmp:
            for mesh in ((2,1,1),(3,1,1),(2,3,1),(2,2,2)):
                with self.subTest(mesh=mesh):
                    nk=int(np.prod(mesh));grid=(8,10,12)
                    model=DFT((2,),np.diag([8.,10.,12.]),[[.1,.2,.3]],
                              (Pseudopotential(2,0,'_rps.dat',b'synthetic'),),2,2,grid,mesh,xc='libxc_pbe')
                    state=NativeGSState(session,2,None,'multik',model);state._identity='synthetic-multik'
                    rng=np.random.default_rng(nk)
                    waves=rng.normal(size=(nk,2,12,10,8))+1j*rng.normal(size=(nk,2,12,10,8))
                    raw=waves.transpose(4,3,2,1,0)[:,:,:,None,:,:,None].copy()
                    occ=np.column_stack((np.linspace(.5,1.5,nk),np.linspace(1.5,.5,nk)))
                    eigen=np.arange(nk*2,dtype=float).reshape(nk,2)/7
                    values={'density':NativeGSArray(np.zeros(grid),(1,1,1),'density','bohr^-3'),
                            'orbitals':NativeGSArray(raw,(1,1,1,1,1,1,1),'orbitals','native'),
                            'occupations':NativeGSArray(occ.T[:,:,None],(1,1,1),'occupations','1'),
                            'eigenvalues':NativeGSArray(eigen.T[:,:,None],(1,1,1),'eigenvalues','hartree')}
                    state.array=Mock(side_effect=lambda name:values[name])
                    info=NativeGSInfo(-2.,1e-12,1e-9,10,True,1,2)
                    target=Path(tmp)/('mesh-'+'-'.join(map(str,mesh)))
                    with patch.object(NativeGSState,'info',new_callable=unittest.mock.PropertyMock,return_value=info):
                        archive=state.export_archive(target,orbital_block=1)
                    for k in range(nk):np.testing.assert_array_equal(archive.wavefunctions(k),waves[k])
                    np.testing.assert_array_equal(archive.array('occupations'),occ)
                    np.testing.assert_array_equal(archive.array('eigenvalues'),eigen)
                    np.testing.assert_array_equal(archive.array('weights'),np.full(nk,1/nk))
                    coordinates=archive.array('kpoints')
                    # Native order has x fastest, then y, then z.
                    for k,coordinate in enumerate(coordinates):
                        ix=k%mesh[0];iy=(k//mesh[0])%mesh[1];iz=k//(mesh[0]*mesh[1])
                        np.testing.assert_allclose(coordinate,[(ix+.5)/mesh[0]-.5,
                            (iy+.5)/mesh[1]-.5,(iz+.5)/mesh[2]-.5],rtol=0,atol=0)
                    expected=sum((occ[k,b]/nk)*abs(waves[k,b])**2 for k in range(nk) for b in range(2))
                    np.testing.assert_allclose(archive.array('density'),expected,rtol=1e-15,atol=0)
                    original=waves.copy();raw[:]=0;state.closed=True
                    reopened=GroundStateArchive.open(target)
                    for k in range(nk):np.testing.assert_array_equal(reopened.wavefunctions(k),original[k])

    def test_retained_export_rejects_unknown_layout_before_writing(self):
        import tempfile
        import numpy as np
        from types import SimpleNamespace
        from unittest.mock import Mock
        from salmonpy.native_gs import NativeGSState
        from salmonpy.errors import InputError
        from salmonpy.dft import DFT,Pseudopotential
        model=DFT((2,),np.eye(3)*8,[[0,0,0]],
                  (Pseudopotential(2,0,'_rps.dat',b'synthetic'),),2,1,(8,8,8),(1,1,1))
        comm=SimpleNamespace(size=1,allgather=lambda v:[v])
        session=SimpleNamespace(_live=lambda:None,comm=comm,backend=SimpleNamespace(decomposition=(1,1,1,1,1)))
        state=NativeGSState(session,1,None,'a',None);state.array=Mock()
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(tmp)/'archive'
            with self.assertRaisesRegex(InputError,'DFT'):state.export_archive(target)
            state._calculation=model;comm.size=2
            with self.assertRaisesRegex(InputError,'decomposition'):state.export_archive(target)
            self.assertFalse(target.exists());state.array.assert_not_called()



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
