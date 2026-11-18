import tempfile
import unittest
from pathlib import Path
import numpy as np

from salmonpy.errors import InputError


def groups():
    return {'system':"yn_periodic='y',al=8,8,8,natom=1,nelem=1",'rgrid':'num_rgrid=8,8,8',
            'pseudo':'izatom(1)=14','atomic_red_coor':"'Si' 0 0 0 1\n"}


class HSEGeometryAdmissionTests(unittest.TestCase):
    def check(self, config, operations=None, files=None):
        from salmonpy._hse_geometry import preflight
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name,text in (files or {}).items():
                (root/name).write_text(text)
            before=dict(config)
            contents={p.name:p.read_bytes() for p in root.iterdir()}
            try:
                return preflight(root,config,operations)
            finally:
                self.assertEqual(config,before)
                self.assertEqual({p.name:p.read_bytes() for p in root.iterdir()},contents)

    def test_explicit_and_spacing_cells(self):
        for body,grid in [('al=3*8','num_rgrid=3*8'),
                          ('al(1:2)=2*8,al(3)=8','dl=3*1'),
                          ('','num_rgrid=8,8,8,dl=1,1,1'),
                          ('al_vec1=8,0,0,al_vec2=0,8,0,al_vec3=0,0,8','num_rgrid=3*8')]:
            g=groups();g['system']="yn_periodic='y',natom=1,nelem=1,"+body;g['rgrid']=grid
            self.check(g)

    def test_general_cells_follow_native_column_lengths(self):
        for vectors in ['al_vec1=0,8,0,al_vec2=-8,0,0,al_vec3=0,0,8',
                        'al_vec1=0,8,0,al_vec2=8,0,0,al_vec3=0,0,8',
                        'al_vec1=8,0,0,al_vec2=1,8,0,al_vec3=0,0,8']:
            g=groups();g['system']+=','+vectors
            np.testing.assert_array_equal(self.check(g),[8,8,8])

    def test_cell_failures(self):
        for body in ['al=0,8,8','al=-8,8,8','al=nan,8,8','al=,8,8,8','al=8,8,8,,','al=3*1d200',
                     'al_vec1=8,1,0,al_vec2=16,2,0,al_vec3=0,0,8']:
            g=groups();g['system']="yn_periodic='y',natom=1,nelem=1,"+body
            with self.subTest(body=body),self.assertRaises(InputError):self.check(g)

    def test_grid_failures(self):
        for body in ['num_rgrid=0,8,8','num_rgrid=8.5,8,8','dl=-1,1,1',
                     'dl=3*1.1','num_rgrid=2147483647,8,8','num_rgrid=999999,999999,999999']:
            g=groups();g['rgrid']=body
            with self.subTest(body=body),self.assertRaises(InputError):self.check(g)

    def test_atom_failures(self):
        for body in ["'Si' nan 0 0 1", "'Si' 1d308 0 0 1", "'Si' '0' 0 0 1", "'Si' 0 0 0 '1'","'Si' 0 0 0 2","'Si' 0 0 0 1\n'Si' .5 0 0 1",'']:
            g=groups();g['atomic_red_coor']=body
            with self.subTest(body=body),self.assertRaises(InputError):self.check(g)

    def test_cell_grid_and_atom_symmetry(self):
        ops=np.zeros((2,3,4));ops[:,:,:3]=np.eye(3);ops[1,:3,:3]=np.diag([-1,-1,1])
        self.check(groups(),ops)
        g=groups();g['atomic_red_coor']="'Si',0,0,0,1 ! comma input";self.check(g,ops)
        g=groups();g['atomic_red_coor']="'Si' .1 0 0 1"
        with self.assertRaisesRegex(InputError,'atom'):self.check(g,ops)
        swap=ops.copy();swap[1,:,:3]=[[0,1,0],[1,0,0],[0,0,1]]
        g=groups();g['system']=g['system'].replace('8,8,8','8,9,8')
        with self.assertRaisesRegex(InputError,'cell'):self.check(g,swap)
        g=groups();g['rgrid']='num_rgrid=8,10,8'
        with self.assertRaisesRegex(InputError,'grid'):self.check(g,swap)
        translated=np.array([np.column_stack((np.eye(3),[1/3,0,0]))])
        with self.assertRaisesRegex(InputError,'grid'):self.check(groups(),translated)

    def test_case_sensitive_file_units_and_left_handed_cell(self):
        g=groups();del g['atomic_red_coor'];g['system']+=",file_atom_coor='Atoms.dat'"
        self.check(g,files={'Atoms.dat':"'Si' 0 0 0 1\n"})
        g['units']="unit_system='A_eV_fs'"
        with self.assertRaisesRegex(InputError,'unit_system'):self.check(g,files={'Atoms.dat':"'Si' 0 0 0 1\n"})
        g=groups();g['system']="yn_periodic='y',natom=1,nelem=1,al_vec1=0,8,0,al_vec2=8,0,0,al_vec3=0,0,8"
        np.testing.assert_array_equal(self.check(g),[8,8,8])

    def test_grid_accumulated_roundoff_is_checked(self):
        g=groups();g['rgrid']='num_rgrid=1000000,8,8'
        op=np.column_stack((np.eye(3),np.zeros(3)))[None,:,:]
        op[0,0,0]+=1e-11
        with self.assertRaisesRegex(InputError,'grid'):self.check(g,op)

    def test_native_admission_rejects_geometry_before_cdll_or_work_copy(self):
        import hashlib,json
        from types import SimpleNamespace
        from unittest.mock import patch
        from salmonpy.fast_hse import FastHSEMaterialProvider
        from salmonpy.fast_hse_mpi import NativeFastHSEBackend
        class FakeComm:
            size=1
            rank=0
            def allgather(self,value):return [value]
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);bundle=root/'initial';bundle.mkdir()
            restart=bundle/'data_for_restart';restart.mkdir()
            for name in ['info.bin','occupation.bin','wfn.bin']:(restart/name).write_bytes(b'synthetic')
            (bundle/'rt.inp').write_text("&calculation theory='tddft_response' /\n"
                "&functional xc='hse06' /\n&propagation propagator='hse_taylor4' /\n"
                "&system yn_periodic='y',al=0,8,8,nelem=1,natom=1 /\n"
                "&rgrid num_rgrid=8,8,8 /\n&atomic_red_coor\n 'Si' 0 0 0 1\n/\n")
            library=root/'library';library.write_bytes(b'not loaded')
            receipt=root/'build.json';receipt.write_text(json.dumps(dict(schema='salmonpy.fast-hse-mpi/1',
                status='built-and-audited',library=str(library),sha256=hashlib.sha256(library.read_bytes()).hexdigest(),origin_sha256='a'*64)))
            provider=FastHSEMaterialProvider(bundle,max_steps=3)
            backend=NativeFastHSEBackend(receipt,root)
            context=SimpleNamespace(solver_comm=FakeComm(),point_ids=(0,))
            before=set(root.iterdir())
            with patch('salmonpy.fast_hse_mpi.c.CDLL',side_effect=AssertionError('no native load')):
                with self.assertRaisesRegex(InputError,'HSE geometry'):
                    backend.create_state(provider.prepare(context),0,.001,context)
            self.assertEqual(set(root.iterdir()),before)

    def test_grid_transformed_index_range(self):
        ops=np.zeros((2,3,4));ops[:,:,:3]=np.eye(3)
        ops[1,:,:3]=np.diag([-1,-1,1]);ops[1,0,3]=1e20
        with self.assertRaisesRegex(InputError,'grid'):self.check(groups(),ops)

class GeneralCellCapabilityTests(unittest.TestCase):
    def test_old_core_is_only_refused_for_general_symmetry(self):
        from salmonpy.fast_hse_mpi import _require_cell_protocol
        from types import SimpleNamespace
        g=groups();g['system']+=',al_vec1=0,8,0,al_vec2=-8,0,0,al_vec3=0,0,8'
        _require_cell_protocol(SimpleNamespace(),g,False)
        _require_cell_protocol(SimpleNamespace(),groups(),True)
        with self.assertRaisesRegex(InputError,'general.cell'):
            _require_cell_protocol(SimpleNamespace(),g,True)
        class Query:
            def __call__(self):return 1
        _require_cell_protocol(SimpleNamespace(sp_fast_hse_general_cell_protocol=Query()),g,True)

    def test_arbitrary_field_core_capability_uses_cartesian_group(self):
        from salmonpy.fast_hse_mpi import _require_cell_protocol
        from types import SimpleNamespace
        class Query:
            def __call__(self):return 1
        ops=np.zeros((2,3,4));ops[:,:,:3]=np.eye(3);ops[1,:,:3]=np.diag([1,-1,-1])
        with self.assertRaisesRegex(InputError,'arbitrary-field'):
            _require_cell_protocol(SimpleNamespace(),groups(),True,ops)
        _require_cell_protocol(SimpleNamespace(sp_fast_hse_symmetry_field_protocol=Query()),groups(),True,ops)
        g=groups();g['system']+=',al_vec1=0,0,8,al_vec2=8,0,0,al_vec3=0,8,0'
        _require_cell_protocol(SimpleNamespace(sp_fast_hse_general_cell_protocol=Query()),g,True,ops)
