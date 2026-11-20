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
            receipt.write_text(json.dumps({'schema':'salmonpy.fast-hse-mpi/1',
                'library':str(library),'sha256':hashlib.sha256(library.read_bytes()).hexdigest(),
                'gs_abi':{'protocol':1,'modes':['import'],'arrays':['density','occupations','orbitals']}}))
            with patch('salmonpy.native_gs.c.CDLL') as load:
                backend=NativeGSBackend(receipt,root)
                load.assert_not_called()
            self.assertFalse(backend._supports_compute)
            computed=json.loads(receipt.read_text());computed['gs_abi']['modes'].append('run')
            receipt.write_text(json.dumps(computed))
            with patch('salmonpy.native_gs.c.CDLL') as load:
                computed_backend=NativeGSBackend(receipt,root)
                load.assert_not_called()
            self.assertTrue(computed_backend._supports_compute)
            self.assertTrue(computed_backend._fast_hse)
            from types import SimpleNamespace
            from unittest.mock import Mock
            import numpy as np
            from salmonpy.native_gs import NativeGSSession
            from salmonpy.dft import DFT,Pseudopotential
            hybrid=DFT((14,),np.eye(3)*8,[[0,0,0]],
                       (Pseudopotential(14,2,'_rps.dat',b'header'),),4,2,(8,8,8),(2,1,1),xc='hse06')
            comm=SimpleNamespace(size=1,allgather=lambda v:[v])
            native=SimpleNamespace(_live=lambda:None,context=SimpleNamespace(solver_comm=comm))
            lib=Mock()
            for backend_options,message in [(dict(_fast_hse=False),'fast-PR'),
                    (dict(_fast_hse=True,decomposition=(1,2,1,1,1)),'k-only')]:
                session=NativeGSSession(SimpleNamespace(_supports_compute=True,**backend_options),native,lib)
                with self.assertRaisesRegex(ValueError,message):session.run(hybrid)
            lib.sp_gs_create.assert_not_called()
            # Old/disabled orbital protocols reject before touching native GS.
            from dataclasses import replace
            gamma=replace(hybrid,kmesh=(1,1,1))
            native2=SimpleNamespace(_live=lambda:None,context=SimpleNamespace(
                solver_comm=SimpleNamespace(size=2,allgather=lambda v:[v,v])))
            orbital_lib=SimpleNamespace(sp_gs_hse_orbital_protocol=Mock(return_value=0),sp_gs_create=Mock())
            session=NativeGSSession(SimpleNamespace(_supports_compute=True,_fast_hse=True,
                decomposition=(1,2,1,1,1)),native2,orbital_lib)
            with self.assertRaisesRegex(ValueError,'matching runtime'):session.run(gamma)
            orbital_lib.sp_gs_create.assert_not_called()
            kspace_lib=SimpleNamespace(sp_gs_hse_kspace_protocol=Mock(return_value=0),sp_gs_create=Mock())
            session=NativeGSSession(SimpleNamespace(_supports_compute=True,_fast_hse=True,
                hse_kspace=True,decomposition=(1,1,2,1,1)),native2,kspace_lib)
            with self.assertRaisesRegex(ValueError,'full mesh'):session.run(hybrid)
            kspace_lib.sp_gs_hse_kspace_protocol.assert_called_once()
            kspace_lib.sp_gs_create.assert_not_called()
            from salmonpy.dft import SymmetryOperations
            operations=np.zeros((1,3,4));operations[0,:,:3]=np.eye(3)
            kspace_lib.sp_gs_hse_kspace_protocol.return_value=1
            with self.assertRaisesRegex(ValueError,'full mesh'):
                session.run(replace(hybrid,symmetry=SymmetryOperations(operations,directions='yyn')))
            kspace_lib.sp_gs_create.assert_not_called()



            kspace_lib.sp_gs_hse_kspace_protocol.return_value=2
            session=NativeGSSession(SimpleNamespace(_supports_compute=True,_fast_hse=True,
                hse_kspace=True,decomposition=(1,2,1,1,1)),native2,kspace_lib)
            with self.assertRaisesRegex(ValueError,'full mesh'):session.run(hybrid)
            kspace_lib.sp_gs_create.assert_not_called()

            malformed=json.loads(receipt.read_text());malformed.pop('gs_abi');receipt.write_text(json.dumps(malformed))
            from salmonpy.errors import InputError
            with self.assertRaises(InputError):NativeGSBackend(receipt,root)
            for name,value in [('decomposition',(1,2,1,1,1)),('library',Path('/tmp/other')),
                               ('work_root',Path('/tmp/other')),('hse_orbital',True),('hse_spatial',True),('hse_kspace',True),('sha256','different')]:
                with self.assertRaises(AttributeError):setattr(backend,name,value)

            # A k-only repartition must reject an empty target k owner before
            # any native GS creation or scientific file decode.
            from salmonpy.native_gs import _scf_disk_manifest
            checkpoint=root/'k-checkpoint';checkpoint.mkdir()
            payload=b'fixture';(checkpoint/'rank-00000000.bin').write_bytes(payload)
            text="&parallel\n nproc_k=1,nproc_ob=1,nproc_rgrid=1,1,1\n/\n&functional\n xc='pz'\n/\n"
            assets=(('Si_rps.dat',b'pp'),)
            manifest={'schema':'salmonpy.native-scf-disk/1','runtime_sha256':'fixture',
                'decomposition':[1,1,1,1,1],'ranks':1,'input_text':text,'material_id':'Si',
                'assets':{'Si_rps.dat':hashlib.sha256(b'pp').hexdigest()},
                'files':[{'name':'rank-00000000.bin','bytes':len(payload),
                          'sha256':hashlib.sha256(payload).hexdigest(),'k_bounds':[1,2]}]}
            (checkpoint/'manifest.json').write_text(json.dumps(manifest))
            probe=Mock(return_value=1)
            local=SimpleNamespace(comm=SimpleNamespace(size=3,rank=0),
                backend=SimpleNamespace(sha256='fixture',decomposition=(3,1,1,1,1)),
                lib=SimpleNamespace(sp_gs_scf_k_repartition_protocol=probe,sp_gs_scf_k_repartition=Mock()))
            with self.assertRaisesRegex(InputError,'more k ranks'):
                _scf_disk_manifest(checkpoint,local,text.replace('nproc_k=1','nproc_k=3'),assets,'Si',repartition=True)
            local.lib.sp_gs_scf_k_repartition.assert_not_called()
            from salmonpy.native_gs import _scf_k_repartition_api,_scf_k_input
            with self.assertRaisesRegex(InputError,'protocol2'):_scf_k_repartition_api(local.lib,hse=True)
            probe.return_value=2
            self.assertIs(_scf_k_repartition_api(local.lib,hse=True),local.lib.sp_gs_scf_k_repartition)
            hse_text=text.replace("xc='pz'","xc='hse06'")
            self.assertIn('functional',_scf_k_input(hse_text,(1,1,1,1,1)))
            # The same manifest admission handles an orbital-only owner, but
            # never treats a fast-HSE orbital layout as a supported profile.
            manifest['files'][0]['orbital_bounds']=[1,2]
            (checkpoint/'manifest.json').write_text(json.dumps(manifest))
            local.backend.decomposition=(1,3,1,1,1)
            local.lib.sp_gs_scf_o_repartition_protocol=Mock(return_value=1)
            local.lib.sp_gs_scf_o_repartition=Mock()
            with self.assertRaisesRegex(InputError,'more orbital ranks'):
                _scf_disk_manifest(checkpoint,local,text.replace('nproc_ob=1','nproc_ob=3'),assets,'Si',repartition=True)
            local.lib.sp_gs_scf_o_repartition.assert_not_called()
            manifest['input_text']=hse_text
            (checkpoint/'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(InputError,'protocol2'):
                _scf_disk_manifest(checkpoint,local,hse_text.replace('nproc_ob=1','nproc_ob=3'),assets,'Si',repartition=True)
            from salmonpy.native_gs import _scf_o_repartition_api
            local.lib.sp_gs_scf_o_repartition_protocol.return_value=2
            self.assertIs(_scf_o_repartition_api(local.lib,hse=True),local.lib.sp_gs_scf_o_repartition)
            manifest['input_text']=text
            (checkpoint/'manifest.json').write_text(json.dumps(manifest))
            local.comm.size=2;local.backend.decomposition=(1,1,2,1,1)
            local.lib.sp_gs_scf_r_repartition_protocol=Mock(return_value=1)
            local.lib.sp_gs_scf_r_repartition=Mock()
            spatial_text=text.replace('nproc_rgrid=1,1,1','nproc_rgrid=2,1,1')
            with self.assertRaisesRegex(InputError,'grid_bounds'):
                _scf_disk_manifest(checkpoint,local,spatial_text,assets,'Si',repartition=True)
            manifest['files'][0]['grid_bounds']=[[1,8],[1,8],[1,8]]
            (checkpoint/'manifest.json').write_text(json.dumps(manifest))
            _,saved=_scf_disk_manifest(checkpoint,local,spatial_text,assets,'Si',repartition=True)
            self.assertEqual(saved['files'][0]['grid_bounds'],[[1,8],[1,8],[1,8]])
            local.lib.sp_gs_scf_r_repartition.assert_not_called()
            from salmonpy.native_gs import _scf_r_repartition_api
            with self.assertRaisesRegex(InputError,'protocol2'):_scf_r_repartition_api(local.lib,hse=True)
            local.lib.sp_gs_scf_r_repartition_protocol.return_value=2
            self.assertIs(_scf_r_repartition_api(local.lib,hse=True),local.lib.sp_gs_scf_r_repartition)
            local.comm.size=8;local.backend.decomposition=(2,2,2,1,1)
            local.lib.sp_gs_scf_m_repartition_protocol=Mock(return_value=1)
            local.lib.sp_gs_scf_m_repartition=Mock()
            mixed_text=text.replace('nproc_k=1','nproc_k=2').replace('nproc_ob=1','nproc_ob=2').replace('nproc_rgrid=1,1,1','nproc_rgrid=2,1,1')
            _,saved=_scf_disk_manifest(checkpoint,local,mixed_text,assets,'Si',repartition=True)
            self.assertEqual(saved['files'][0]['orbital_bounds'],[1,2])
            local.lib.sp_gs_scf_m_repartition.assert_not_called()
            manifest['input_text']=hse_text
            (checkpoint/'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(InputError,'Gamma'):
                _scf_disk_manifest(checkpoint,local,mixed_text.replace("xc='pz'","xc='hse06'"),assets,'Si',repartition=True)
            from salmonpy.native_gs import _scf_m_repartition_api
            with self.assertRaisesRegex(InputError,"protocol"):_scf_m_repartition_api(local.lib,hse=True)
            local.lib.sp_gs_scf_m_repartition_protocol.return_value=2
            self.assertIs(_scf_m_repartition_api(local.lib,hse=True),local.lib.sp_gs_scf_m_repartition)
            with self.assertRaisesRegex(InputError,'protocol'):_scf_m_repartition_api(local.lib,hse=True,kspace=True)
            local.lib.sp_gs_scf_m_repartition_protocol.return_value=3
            self.assertIs(_scf_m_repartition_api(local.lib,hse=True,kspace=True),local.lib.sp_gs_scf_m_repartition)







    def test_imported_force_and_native_complex_array_contract(self):
        from salmonpy.native_gs import NativeGSState
        from salmonpy.errors import InputError
        from types import SimpleNamespace
        from unittest.mock import Mock
        session=SimpleNamespace(_live=lambda:None,comm=SimpleNamespace(allgather=lambda v:[v]),lib=Mock())
        state=NativeGSState(session,1,None,'imported',None)
        with self.assertRaisesRegex(InputError,'forces were not computed'):state.array('forces')
        session.lib.sp_gs_array_info.assert_not_called()
        # All ranks must select the same supported quantity before C/MPI calls.
        for invalid in ([], {}, None, 1):
            with self.subTest(quantity=invalid),self.assertRaises(InputError):state.array(invalid)
        session.lib.sp_gs_array_info.assert_not_called()
        def gathered(value):
            return [value,'occupations'] if value=='density' else [value,value]
        session.comm=SimpleNamespace(allgather=gathered)
        with self.assertRaisesRegex(InputError,'quantities differ'):state.array('density')
        session.lib.sp_gs_array_info.assert_not_called()
        session.comm=SimpleNamespace(allgather=lambda value:[value])
        import ctypes as c
        import numpy as np
        from salmonpy.errors import NumericalError
        expected=np.empty(2,dtype=np.complex128)
        expected.real=[1.,-0.];expected.imag=[-0.,0.]
        def layout(token,quantity,ndim,shape,lower):
            ndim._obj.value=7
            np.ctypeslib.as_array((c.c_int64*7).from_address(shape))[:]=[2,1,1,1,1,1,1]
            np.ctypeslib.as_array((c.c_int64*7).from_address(lower))[:]=1
            return 0
        def read(token,quantity,re,im,length):
            np.ctypeslib.as_array((c.c_double*length).from_address(re))[:]=expected.real
            np.ctypeslib.as_array((c.c_double*length).from_address(im))[:]=expected.imag
            return 0
        session.lib=SimpleNamespace(sp_gs_array_info=layout,sp_gs_read=read)
        session._check=lambda status:None
        result=state.array('orbitals').values
        self.assertEqual(result.tobytes(),expected.tobytes())
        self.assertFalse(result.flags.writeable)
        for value in (np.inf,np.nan,-np.inf):
            expected.imag[0]=value
            with np.errstate(all='raise'),self.assertRaises(NumericalError):state.array('orbitals')

        session.backend=SimpleNamespace(_fast=True)
        with self.assertRaisesRegex(InputError,'eigenvalues were not computed'):state.array('eigenvalues')
        def metadata(token,stats,values):
            np.ctypeslib.as_array((c.c_int*4).from_address(stats))[:]=[0,-1,2,4]
            np.ctypeslib.as_array((c.c_double*3).from_address(values))[:]=0
            return 0
        session.lib.sp_gs_info=metadata
        self.assertIsNone(state.info.energy_hartree)
        self.assertIsNone(state.info.residual)
        self.assertIsNone(state.info.converged)

    @unittest.skipUnless(importlib.util.find_spec('f90nml'),'optional f90nml not installed')
    def test_declared_gs_metadata_normalization_and_independence(self):
        import numpy as np
        from salmonpy.gs_input import describe_gs_input
        from salmonpy.ase import BOHR_ANGSTROM
        text="""&units
 unit_system='A_eV_fs'
/
&system
 yn_periodic='y', spin='unpolarized', nelec=4, nstate=3, natom=1,
 al=4.23341768536,4.23341768536,4.23341768536
/
&rgrid
 num_rgrid=8,10,12
/
&kgrid
 file_kw='tables/k!.dat'
/
&pseudo
 izatom(2)=14, lloc_ps(2)=2, file_pseudo(2)='pot/Si_rps.dat'
/
&functional
 xc='libxc_pbe'
/
&atomic_coor
 'Si' 2.11670884268 0 0 2
/
"""
        assets={'tables/k!.dat':b'2\n1 -.25 0 0 .5\n2 .25 0 0 .5\n','pot/Si_rps.dat':b'synthetic'}
        metadata=describe_gs_input(text,assets=assets,material_id='film-a')
        data=metadata.as_dict();settings=data['settings']
        np.testing.assert_allclose(settings['cell_bohr'],np.eye(3)*8,rtol=0,atol=1e-12)
        np.testing.assert_allclose(settings['fractional_positions'],[[.5,0,0]],rtol=0,atol=1e-15)
        self.assertEqual(settings['atomic_numbers'],[14])
        self.assertEqual(settings['grid_counts'],[8,10,12])
        self.assertEqual(settings['pseudopotentials'][0]['species_index'],2)
        self.assertEqual(settings['pseudopotentials'][0]['local_channel'],2)
        self.assertEqual(settings['k_sampling']['input_weights'],[.5,.5])
        self.assertIsNone(data['evidence']['scf_converged'])
        self.assertEqual(data['evidence']['checkpoint_provenance'],'unknown')
        assets.clear();data['settings']['cell_bohr'][0][0]=999
        self.assertAlmostEqual(metadata.settings['cell_bohr'][0][0],8)
        with self.assertRaises(TypeError):metadata.settings['electrons']=8
        skew=describe_gs_input("&units\n unit_system='au'\n/\n&system\n natom=1,nelec=2,nstate=2,al_vec1=8,0,0,al_vec2=1,9,0,al_vec3=0,2,10,file_atom_coor='coords/a.dat'\n/\n&pseudo\n izatom=2,lloc_ps=0,file_pseudo='p.dat'\n/\n",assets={'coords/a.dat':b"'He' 2.5 6 7.5 1\n",'p.dat':b'synthetic'})
        np.testing.assert_array_equal(skew.settings['cell_bohr'],[[8,0,0],[1,9,0],[0,2,10]])
        np.testing.assert_allclose(skew.settings['fractional_positions'],[[.25,.5,.75]],rtol=0,atol=1e-15)
        self.assertEqual(metadata.source_id,describe_gs_input(text,assets={
            'tables/k!.dat':b'2\n1 -.25 0 0 .5\n2 .25 0 0 .5\n','pot/Si_rps.dat':b'synthetic'},material_id='other-id').source_id)

    @unittest.skipUnless(importlib.util.find_spec('f90nml'),'optional f90nml not installed')
    def test_gs_metadata_unknowns_and_optional_dependency_do_not_claim_defaults(self):
        from salmonpy.gs_input import describe_gs_input
        from salmonpy.native_gs import NativeGSState
        text="&system\n nelec=4,nstate=2,al=8,8,8\n/\n&functional\n xc='pz'\n/\n"
        metadata=describe_gs_input(text,assets={},material_id='unknown')
        self.assertIsNone(metadata.settings['cell_bohr'])
        self.assertIsNone(metadata.settings['fractional_positions'])
        self.assertIsNone(metadata.settings['k_sampling'])
        self.assertTrue(metadata.notes)
        self.assertFalse(describe_gs_input(text,assets={},origin='computed',converged=False).as_dict()['evidence']['scf_converged'])
        state=NativeGSState(None,1,None,'unknown',None)
        state._input_text=text
        first=state.input_metadata()
        state.closed=True
        self.assertEqual(first.as_dict(),state.input_metadata().as_dict())
        with patch.dict('sys.modules',{'f90nml':None}):
            with self.assertRaisesRegex(ImportError,'metadata'):describe_gs_input(text,assets={})

    def test_existing_gs_import_preserves_external_k_table_and_rejects_missing_asset(self):
        import tempfile
        from types import SimpleNamespace
        from unittest.mock import Mock
        from salmonpy.native_gs import NativeGSSession
        from salmonpy.errors import InputError
        comm=SimpleNamespace(size=1,rank=0,allgather=lambda v:[v],bcast=lambda v,root:v)
        native=SimpleNamespace(_live=lambda:None,context=SimpleNamespace(solver_comm=comm))
        def imported(filename,pointer):
            pointer._obj.value=123
            return 0
        fn=Mock(side_effect=imported)
        lib=SimpleNamespace(sp_gs_import=fn)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'source';source.mkdir();work=root/'work';work.mkdir()
            (source/'data_for_restart').mkdir()
            for name in ('info.bin','occupation.bin','wfn.bin'):
                (source/'data_for_restart'/name).write_bytes(b'synthetic placeholder, never passed to Fortran')
            (source/'pseudo_rps.dat').write_bytes(b'synthetic potential')
            (source/'tables').mkdir();table=source/'tables/k!.dat'
            table.write_text('2\n1 -.25 0 0 .5\n2 .25 0 0 .5\n')
            (source/'rt.inp').write_text("&system\n nelec=4\n/\n&pseudo\n file_pseudo(1)='pseudo_rps.dat'\n/\n&kgrid\n file_kw='tables/k!.dat'\n/\n")
            session=NativeGSSession(SimpleNamespace(work_root=work,decomposition=(1,1,1,1,1),sha256='a'*64),native,lib)
            state=session.load(source,material_id='silicon')
            self.assertEqual((state.directory/'tables/k!.dat').read_bytes(),table.read_bytes())
            self.assertIn(('tables/k!.dat',table.read_bytes()),state._assets)
            self.assertEqual(state.material_id,'silicon')
            self.assertEqual(fn.call_count,1)
            session.backend._ordinary_rt=True
            bad=source/'unsupported.inp'
            bad.write_text((source/'rt.inp').read_text()+"&functional\n xc='hse06'\n/\n")
            with self.assertRaisesRegex(InputError,'HSE-disabled'):
                session.load(source,'unsupported.inp')
            self.assertEqual(fn.call_count,1)  # unsupported hybrid never reaches native initialization
            self.assertEqual(len(session._states),1)
            if importlib.util.find_spec('f90nml'):
                description=state.input_metadata()
                self.assertEqual(description.settings['k_sampling']['mode'],'file')
                self.assertEqual(description.as_dict()['evidence']['origin'],'imported')
            table.unlink()
            if importlib.util.find_spec('f90nml'):
                self.assertEqual(description.as_dict(),state.input_metadata().as_dict())
            with self.assertRaisesRegex(InputError,'k!.dat'):session.load(source)
            self.assertEqual(fn.call_count,1)
            self.assertEqual(len(list(work.glob('gs-import-*'))),1)

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

    def test_scf_snapshot_owner_and_step_admission_precede_native_calls(self):
        from types import SimpleNamespace
        from unittest.mock import Mock
        from salmonpy.native_gs import NativeGSState,NativeSCFSnapshot
        from salmonpy.errors import InputError
        comm=SimpleNamespace(allgather=lambda value:[value])
        session=SimpleNamespace(comm=comm,_live=lambda:None,lib=Mock())
        first=NativeGSState(session,1,None,'A',None)
        second=NativeGSState(session,2,None,'B',None)
        saved=NativeSCFSnapshot(first,7)
        with self.assertRaises(InputError):second.restore(saved)
        with self.assertRaises(InputError):first.continue_scf(steps=0)
        session.lib.sp_gs_scf_restore.assert_not_called()
        session.lib.sp_gs_continue.assert_not_called()

    def test_common_scf_target_binds_local_tokens_and_rejects_wrong_owner_or_atomic_fallback(self):
        from types import SimpleNamespace
        from unittest.mock import Mock
        from salmonpy.native_gs import NativeGSState,NativeSCFSnapshot
        from salmonpy.checkpointing import execute_checkpoint,restore_checkpoint,CheckpointState
        from salmonpy.errors import InputError
        comm=SimpleNamespace(allgather=lambda value:[value])
        session=SimpleNamespace(comm=comm,_live=lambda:None,
            backend=SimpleNamespace(decomposition=(1,1,1,1,1)),
            lib=SimpleNamespace(sp_gs_scf_info=lambda *args:0),_check=lambda code:None)
        state=NativeGSState(session,1,None,'film',None)
        state._identity='model/runtime';state._scf_ready=Mock()
        other=NativeGSState(session,2,None,'film',None)
        other._identity=state._identity;other._scf_ready=Mock()
        saved=NativeSCFSnapshot(state,7)
        state.snapshot=Mock(return_value=saved);state.restore=Mock()
        groups=[{'name':'gs','ranks':[0],'point_ids':['film']}]
        target=state.checkpoint_target()
        self.assertIsNone(target.restore);state.snapshot.assert_not_called()
        result=execute_checkpoint(None,groups,states={'film':target})
        self.assertIs(result['states']['film'],saved)
        self.assertNotEqual(state.checkpoint_capability.identity,other.checkpoint_capability.identity)
        with self.assertRaises(InputError):other.checkpoint_target(saved)
        with self.assertRaises(InputError):
            execute_checkpoint(None,groups,states={'film':other.checkpoint_target()},
                               operation='restore',saved=result['compatibility'])
        restore=state.checkpoint_target(saved)
        with self.assertRaises(InputError):
            restore_checkpoint(None,groups,states={'film':restore},
                live=CheckpointState(states={'film':state}),saved=result['compatibility'])
        state.restore.assert_not_called()
        execute_checkpoint(None,groups,states={'film':restore},operation='restore',saved=result['compatibility'])
        state.restore.assert_called_once_with(saved)
        saved.closed=True
        with self.assertRaises(InputError):state.checkpoint_target(saved)

    def test_scf_candidate_reservation_stage_discard_and_publication(self):
        from types import SimpleNamespace
        from unittest.mock import Mock
        from salmonpy.native_gs import _SCFRestoreCandidate
        from salmonpy.errors import InputError
        calls=Mock(return_value=0)
        session=SimpleNamespace(comm=SimpleNamespace(allgather=lambda v:[v]),
            _check=lambda rc:None,_gs_call=lambda directory,fn,*args:fn(*args))
        state=SimpleNamespace(session=session,_scf_ready=lambda:None,directory=None,_token=1,
            _restore_candidate=object(),info=SimpleNamespace(converged=False))
        snapshot=SimpleNamespace(closed=False,_token=2)
        candidate=_SCFRestoreCandidate(state,snapshot,calls)
        with self.assertRaises(InputError):candidate.stage()
        calls.assert_not_called()
        state._restore_candidate=None
        self.assertIs(candidate.stage(),candidate)
        candidate.discard();self.assertIsNone(state._restore_candidate)
        self.assertEqual([call.args[-1] for call in calls.call_args_list],[2,4])
        with self.assertRaises(InputError):candidate.publish()
        candidate=_SCFRestoreCandidate(state,snapshot,calls)
        candidate.stage();self.assertIs(candidate.publish(),state)
        self.assertIsNone(state._restore_candidate)
        self.assertFalse(state._metadata_converged)
        with self.assertRaises(InputError):candidate.stage()
        def reject(code):
            if code:raise InputError('native admission')
        candidate=_SCFRestoreCandidate(state,snapshot,calls)
        with self.assertRaisesRegex(RuntimeError,'host failure'):
            with candidate.stage() as staged:
                self.assertIs(staged,candidate)
                raise RuntimeError('host failure')
        self.assertIsNone(state._restore_candidate)
        candidate.close()  # idempotent after the exceptional context exit
        session._check=reject;failed=Mock(return_value=11)
        candidate=_SCFRestoreCandidate(state,snapshot,failed)
        with self.assertRaises(InputError):candidate.stage()
        candidate.discard()
        self.assertIsNone(state._restore_candidate)
        self.assertEqual([call.args[-1] for call in failed.call_args_list],[2])
