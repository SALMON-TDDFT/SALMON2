import sys,unittest,tempfile,struct,importlib.util
from pathlib import Path
import numpy as np
MODULE=Path(__file__).resolve().parents[2]/'samples/dc_hse/lcfo_rt_reference.py'
sys.path.insert(0,str(MODULE.parent))

def fixture(root,nf=2):
    blocks=[np.diag([1.+2*f,2.+2*f]).astype(complex) for f in range(nf)]
    link=np.array([[.02j,.03],[.04,-.01j]])
    h=np.zeros((2*nf,2*nf),complex)
    for f in range(nf):h[2*f:2*f+2,2*f:2*f+2]=blocks[f]
    if nf==2:
        h[2:4,:2]=link;h[:2,2:4]=link.conj().T
    else:
        for f in range(nf):
            prev=(f-1)%nf
            h[2*prev:2*prev+2,2*f:2*f+2]=link
            h[2*f:2*f+2,2*prev:2*prev+2]=link.conj().T
    e,c=np.linalg.eigh(h)
    def header(kind,f):
        meta=[2*nf,1,1,4,1,1,2,1,1,nf,1,1,1+2*f,1,1,f+1,1,1,2,2*nf]
        geom=[float(2*nf),0,0,0,1.,0,0,0,1.,1.]
        return (b'SLCFO_COMPLEX_V1'.ljust(16)+struct.pack('<6iq',1,0x01020304,32,64,kind,0,336)
              + b'test-run'.ljust(96)+struct.pack('<20i10d',*meta,*geom)+struct.pack('<4d',0,0,0,1))
    def z(a):return np.asarray(a,dtype='<c16').tobytes(order='F')
    def save(kind,f,name,payload):
        data=header(kind,f)+struct.pack('<iq',1,len(payload))+payload
        data+=b'SLCFO_DONE_V1'.ljust(16)+b'test-run'.ljust(96)+struct.pack('<iq',1,len(data)+124)
        path=root/f'{f+1:06d}';path.mkdir(exist_ok=True)
        (path/name).write_bytes(data)
    for f in range(nf):
        save(1,f,'basis_functions.bin',struct.pack('<2i',1,2)+z(np.eye(2)))
        save(2,f,'wavefunctions.bin',struct.pack('<'+'i'*(5+nf),1,2,2*nf,*([2]*nf),1+2*f,2+2*f)+z(c[2*f:2*f+2]))
        halo=link if f==0 else link.conj().T
        payload=struct.pack('<'+'i'*(3+nf),1,2,2*nf,*([2]*nf))+z(blocks[f])
        if nf==2:
            payload+=struct.pack('<6i',1,2-f,1 if f==0 else -1,0,0,2)+z(halo)
        else:
            payload+=struct.pack('<i',2)
            for src,direction,block in [((f-1)%nf,1,link),((f+1)%nf,-1,link.conj().T)]:
                payload+=struct.pack('<5i',src+1,direction,0,0,2)+z(block)
        save(3,f,'hamiltonian_local.bin',payload)
    return h,e,c

class ReferenceTest(unittest.TestCase):
    def setUp(self):
        self.assertTrue(MODULE.exists(),'missing strict complex LCFO reader')
        from lcfo_rt_reference import load_lcfo
        self.load=load_lcfo
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.h,self.e,self.c=fixture(self.root)
    def tearDown(self):
        if hasattr(self,'tmp'):self.tmp.cleanup()
    def test_actual_assembly_and_eigenstates(self):
        data=self.load(self.root)
        np.testing.assert_allclose(data['hamiltonian'],self.h,atol=1e-15)
        np.testing.assert_allclose(data['coefficients'],self.c,atol=1e-15)
        self.assertLess(data['basis_orthogonality_error'],1e-14)
        self.assertLess(data['eigen_residual_Ha'],1e-14)
    def test_native_final_symmetrization(self):
        p=self.root/'000001/hamiltonian_local.bin'
        raw=bytearray(p.read_bytes());struct.pack_into('<d',raw,376,5e-12);p.write_bytes(raw)
        data=self.load(self.root)
        self.assertGreater(data['hermiticity_relative'],0.)
        np.testing.assert_allclose(data['hamiltonian'],self.h,rtol=0,atol=1e-15)
    def test_halo_on_unsplit_axis_rejected(self):
        p=self.root/'000001/hamiltonian_local.bin';raw=bytearray(p.read_bytes())
        struct.pack_into('<3i',raw,440,0,1,0);p.write_bytes(raw)
        with self.assertRaises(ValueError):self.load(self.root)
    def test_wrong_periodic_neighbor_rejected(self):
        root=self.root/'three';root.mkdir();h,e,c=fixture(root,3)
        np.testing.assert_allclose(self.load(root)['hamiltonian'],h,atol=1e-15)
        p=root/'000001/hamiltonian_local.bin';raw=bytearray(p.read_bytes())
        struct.pack_into('<i',raw,444,-1);p.write_bytes(raw)
        with self.assertRaises(ValueError):self.load(root)
    def test_truncation_rejected(self):
        p=self.root/'000001/basis_functions.bin';p.write_bytes(p.read_bytes()[:-1])
        with self.assertRaises(ValueError):self.load(self.root)
    def test_stale_run_rejected(self):
        p=self.root/'000002/wavefunctions.bin';p.write_bytes(p.read_bytes().replace(b'test-run',b'old--run'))
        with self.assertRaises(ValueError):self.load(self.root)
    def test_payload_overflow_rejected(self):
        p=self.root/'000001/basis_functions.bin';raw=bytearray(p.read_bytes());struct.pack_into('<q',raw,340,2**62);p.write_bytes(raw)
        with self.assertRaises(ValueError):self.load(self.root)
    def test_bad_basis_rejected(self):
        p=self.root/'000001/basis_functions.bin';raw=bytearray(p.read_bytes());struct.pack_into('<d',raw,356,2.);p.write_bytes(raw)
        with self.assertRaises(ValueError):self.load(self.root)
    def test_bad_coefficients_rejected(self):
        p=self.root/'000001/wavefunctions.bin';raw=bytearray(p.read_bytes());struct.pack_into('<d',raw,376,7.);p.write_bytes(raw)
        with self.assertRaises(ValueError):self.load(self.root)
if __name__=='__main__':unittest.main()
