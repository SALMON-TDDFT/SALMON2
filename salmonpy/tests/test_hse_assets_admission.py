import struct
import tempfile
from pathlib import Path
import unittest
import numpy as np
from salmonpy.errors import InputError


def fixture(root):
    restart=root/'data_for_restart';restart.mkdir()
    record=lambda x:struct.pack('=i',len(x))+x+struct.pack('=i',len(x))
    (restart/'info.bin').write_bytes(b''.join(record(struct.pack('=i',x)) for x in [8,2,0,1,0]))
    (restart/'occupation.bin').write_bytes(record(np.full(16,2.,dtype='=f8').tobytes()))
    (restart/'wfn.bin').write_bytes(np.ones(8*8*8*2*8,dtype='=c16').tobytes())
    (root/'Si_rps.dat').write_text('1 .1 0 4\n0 1\n0 0 -1\n.1 0 -1\n0 0\n.1 0\n')
    return {'control':"directory_read_data='./data_for_restart/'",'system':"yn_periodic='y',al=3*8,natom=1,nelem=1,nstate=2,nelec=4",
            'kgrid':'num_kgrid=2,2,2','pseudo':"file_pseudo(1)='Si_rps.dat',izatom(1)=14,lloc_ps(1)=0"}


class HSEAssetsAdmissionTests(unittest.TestCase):
    def check(self,g,root,operations=None):
        from salmonpy._hse_assets import preflight
        return preflight(root,g,operations,np.array([8,8,8]),2)

    def test_valid_single_global_records_and_pp(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.check(fixture(root),root)

    def test_kgrid_and_shift_failures(self):
        for grid in ['num_kgrid=0,2,2','num_kgrid=3,2,2',"file_kw='points.dat'",'num_kgrid=2,2,2,yn_gamma_centered=\'y\'']:
            with self.subTest(grid=grid),tempfile.TemporaryDirectory() as d:
                root=Path(d);g=fixture(root);g['kgrid']=grid
                with self.assertRaises(InputError):self.check(g,root)

    def test_incompatible_star(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);g=fixture(root);g['kgrid']+=' ,dk_shift=.1,0,0'
            ops=np.zeros((2,3,4));ops[:,:,:3]=np.eye(3);ops[1,:,:3]=np.diag([-1,-1,1])
            with self.assertRaisesRegex(InputError,'k'):self.check(g,root,ops)

    def test_missing_or_truncated_potential(self):
        for content in [None,'','1 .1 0 4\n0 1\n']:
            with self.subTest(content=content),tempfile.TemporaryDirectory() as d:
                root=Path(d);g=fixture(root)
                if content is None:(root/'Si_rps.dat').unlink()
                else:(root/'Si_rps.dat').write_text(content)
                with self.assertRaisesRegex(InputError,'pseudo'):self.check(g,root)

    def test_native_ky_list_directed_continuation_and_trailing_fields(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);g=fixture(root)
            (root/'Si_rps.dat').write_text('1 1.0-1\n0 4 header annotation\n0\n1 extra\n0 0\n-1 trailing\n.1 0 -1 ignored\n0 0\n.1 0\n')
            self.check(g,root)

    def test_unfolding_uses_native_namelist_group(self):
        from salmonpy._hse_assets import _mesh
        with self.assertRaisesRegex(InputError,'unfolding'):
            _mesh({'unfolding':"dm_unfold_option='primitive'"},None)

    def test_corrupt_or_mismatched_gs_records(self):
        for name,transform in [('info.bin',lambda b:b[:-1]),('wfn.bin',lambda b:b[:-16]),
                                ('occupation.bin',lambda b:b[:4]+np.full(16,np.nan,dtype='=f8').tobytes()+b[-4:]),
                                ('occupation.bin',lambda b:b[:-1])]:
            with self.subTest(name=name),tempfile.TemporaryDirectory() as d:
                root=Path(d);g=fixture(root);p=root/'data_for_restart'/name;p.write_bytes(transform(p.read_bytes()))
                with self.assertRaises(InputError):self.check(g,root)

    def test_deep_wavefunction_check_is_explicit(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);g=fixture(root);p=root/'data_for_restart/wfn.bin'
            with p.open('r+b') as f:f.write(struct.pack('=d',float('nan')))
            self.check(g,root)  # Default admission does not rescan the full orbital stream.
            from salmonpy._hse_assets import _wavefunction
            with self.assertRaisesRegex(InputError,'wavefunction'):_wavefunction(p,p.stat().st_size)

    def test_occupation_quality_warns_without_changing_data(self):
        for value in [1.999999, -0.001, 2.001]:
            with self.subTest(value=value), tempfile.TemporaryDirectory() as d:
                root=Path(d);g=fixture(root);p=root/'data_for_restart/occupation.bin'
                before=struct.pack('=i',128)+np.full(16,value,dtype='=f8').tobytes()+struct.pack('=i',128)
                p.write_bytes(before)
                with self.assertWarns(RuntimeWarning):self.check(g,root)
                self.assertEqual(p.read_bytes(),before)

    def test_generated_nonuniform_star_weights_and_odd_gamma_shift(self):
        from salmonpy._hse_assets import _mesh,_weights
        ops=np.zeros((2,3,4));ops[:,:,:3]=np.eye(3);ops[1,:,:3]=np.diag([-1,-1,1])
        for text in ['num_kgrid=3,3,1',"num_kgrid=3,3,1,yn_gamma_centered='y',dk_shift=.1,.2,.3"]:
            mesh,matrices,base=_mesh({'kgrid':text},ops)
            np.testing.assert_allclose(_weights(mesh,matrices,base,5,2),[2/9]*4+[1/9],rtol=0,atol=1e-15)

    def test_species_filename_arrays_and_unused_conduction_states(self):
        from salmonpy._hse_assets import _string_array
        self.assertEqual(_string_array("file_pseudo(1:2)='Case_rps.dat','other_rps.dat'",'file_pseudo',2),
                         ['Case_rps.dat','other_rps.dat'])
        self.assertEqual(_string_array("file_pseudo=2*'Same_rps.dat'",'file_pseudo',2),['Same_rps.dat']*2)
        with self.assertRaises(InputError):_string_array("file_pseudo=,'bad_rps.dat'",'file_pseudo',1)

    def test_k_mesh_does_not_impose_polarization(self):
        from salmonpy._hse_assets import _mesh
        ops=np.zeros((2,3,4));ops[:,:,:3]=np.eye(3);ops[1,2,2]=-1
        _mesh({'kgrid':'num_kgrid=2,2,2'},ops)

    def test_general_cell_k_mesh_is_independent_of_polarization(self):
        from salmonpy._hse_assets import _mesh
        g={'kgrid':'num_kgrid=2,2,2','system':'al_vec1=0,0,8,al_vec2=8,0,0,al_vec3=0,8,0'}
        ops=np.zeros((2,3,4));ops[:,:,:3]=np.eye(3);ops[1,:,:3]=np.diag([1,-1,-1])
        _mesh(g,ops)  # Fractional C2x is Cartesian C2z.
        ops[1,:,:3]=np.diag([-1,-1,1])
        _mesh(g,ops)  # Actual field compatibility is checked at each step.

    def test_extra_stored_states_and_real_wavefunction(self):
        for real in [0,1]:
            with self.subTest(real=real),tempfile.TemporaryDirectory() as d:
                root=Path(d);g=fixture(root);restart=root/'data_for_restart'
                rec=lambda b:struct.pack('=i',len(b))+b+struct.pack('=i',len(b))
                (restart/'info.bin').write_bytes(b''.join(rec(struct.pack('=i',x)) for x in [8,3,0,7,real]))
                (restart/'occupation.bin').write_bytes(rec(np.tile([2.,2.,0.],8).astype('=f8').tobytes()))
                dtype='=f8' if real else '=c16'
                (restart/'wfn.bin').write_bytes(np.ones(512*3*8,dtype=dtype).tobytes())
                self.check(g,root)

    def test_assets_reject_before_native_lifetime_or_work_copy(self):
        import hashlib,json
        from types import SimpleNamespace
        from unittest.mock import patch
        from salmonpy.fast_hse import FastHSEMaterialProvider
        from salmonpy.fast_hse_mpi import NativeFastHSEBackend
        class Comm:
            size=1;rank=0
            def allgather(self,value):return [value]
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);bundle=root/'initial';bundle.mkdir();g=fixture(bundle)
            g.update(calculation="theory='tddft_response'",functional="xc='hse06'",
                     propagation="propagator='hse_taylor4'",rgrid='num_rgrid=3*8',
                     atomic_red_coor="\n 'Si' 0 0 0 1\n")
            (bundle/'rt.inp').write_text('\n'.join('&'+k+' '+v+'/' for k,v in g.items()))
            path=bundle/'data_for_restart/info.bin';path.write_bytes(path.read_bytes()[:-1])
            library=root/'library';library.write_bytes(b'not loaded')
            receipt=root/'build.json';receipt.write_text(json.dumps(dict(schema='salmonpy.fast-hse-mpi/1',
                status='built-and-audited',library=str(library),sha256=hashlib.sha256(library.read_bytes()).hexdigest(),origin_sha256='a'*64)))
            backend=NativeFastHSEBackend(receipt,root);provider=FastHSEMaterialProvider(bundle,max_steps=3)
            context=SimpleNamespace(solver_comm=Comm(),point_ids=(0,));before=set(root.iterdir())
            with patch('salmonpy.fast_hse_mpi.c.CDLL',side_effect=AssertionError('no native load')):
                with self.assertRaisesRegex(InputError,'GS'):
                    backend.create_state(provider.prepare(context),0,.001,context)
            self.assertEqual(set(root.iterdir()),before)
