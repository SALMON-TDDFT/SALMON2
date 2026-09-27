"""Static DC frozen-orbital force diagnostics; not a certification of truncated DC-MD."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE'), 'SALMON_TEST_EXE is required')
class DCForceTest(unittest.TestCase):
    def run_case(self, ranks=1, buffer=4, delta=0., conventional=False, xc='pbeh40_rvv10'):
        inp=Path(__file__).with_name('dc_hydrogen.inp').read_text()
        inp=inp.replace('nproc_k=2',f'nproc_k={2 if ranks==4 else 1}')
        inp=inp.replace('nproc_rgrid_tot=4,1,1',f'nproc_rgrid_tot={ranks},1,1')
        inp=inp.replace('num_rgrid_buffer=4,0,0',f'num_rgrid_buffer={buffer},0,0')
        inp=inp.replace('threshold=1d-8','threshold=1d-10').replace('nscf=200','nscf=500')
        inp=inp.replace('3.3d0',f'{3.3+delta:.12f}').replace("xc='pbeh40_rvv10'",f"xc='{xc}'")
        if ranks==1:
            inp=inp.replace('num_fragment=2,1,1','num_fragment=1,1,1')
            inp=inp.replace(f'num_rgrid_buffer={buffer},0,0','num_rgrid_buffer=0,0,0')
        if conventional:inp=inp.replace("yn_dc='y'","yn_dc='n'")
        else:inp=inp.replace('&dc\n',"&dc\n yn_dc_force_diagnostic='y'\n")
        return self.run_input(inp,ranks,conventional)

    def run_input(self,inp,ranks=1,conventional=False):
        natom=int(re.search(r'natom\s*=\s*(\d+)',inp)[1])
        command=[os.environ['SALMON_TEST_EXE']]
        if ranks>1:command=[os.environ['SALMON_TEST_MPIEXEC'],'-n',str(ranks)]+command
        with tempfile.TemporaryDirectory(prefix='dc-force-') as tmp:
            for atom in ('H','O'):shutil.copy(ROOT/f'testsuites/pseudo/{atom}_rps.dat',tmp)
            run=subprocess.run(command,input=inp,cwd=tmp,text=True,capture_output=True,timeout=180,
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
            if os.environ.get('SALMON_TEST_SAVE_DIR'):
                record=Path(os.environ['SALMON_TEST_SAVE_DIR'])/Path(tmp).name
                record.mkdir(parents=True)
                (record/'inputfile').write_text(inp)
                (record/'output').write_text(run.stdout+run.stderr)
            self.assertEqual(run.returncode,0,run.stdout[-2000:]+run.stderr)
            self.assertIn('end SALMON',run.stdout)
            info=next(Path(tmp).rglob('*_info.data')).read_text()
            energy=float(re.search(r'Total energy \(eV\) =\s*(\S+)',info)[1])/27.211386245988
            if conventional:
                self.assertTrue('#GS converged' in run.stdout,'\n'.join(line for line in run.stdout.splitlines() if 'iter=' in line or '|rho_i' in line or '#GS' in line)[-3000:])
                lines=run.stdout.split('===== force =====')[1].strip().splitlines()[:natom]
            else:
                self.assertIn('DC frozen-orbital force diagnostic (Ha/bohr)',run.stdout)
                self.last_ts=float(re.search(r'DC occupation TS diagnostic \(Ha\):\s*(\S+)',run.stdout)[1])
                self.last_free_energy=float(re.search(r'DC E-minus-TS diagnostic \(Ha\):\s*(\S+)',run.stdout)[1])
                residual=float(re.search(r'DC occupation electron residual:\s*(\S+)',run.stdout)[1])
                self.assertLess(abs(residual),1e-8)
                differences=re.findall(r'DC #SCF.*diff =\s*(\S+)',run.stdout)
                self.assertLess(float(differences[-1]),1e-10)
                lines=run.stdout.split('DC frozen-orbital force diagnostic (Ha/bohr)')[1].strip().splitlines()[:natom]
            forces=[[float(x) for x in line.split()[1:4]] for line in lines]
            if conventional and "unit_system='A_eV_fs'" in inp:
                forces=[[x*.529177210903/27.211386245988 for x in row] for row in forces]
            return energy,forces

    def test_one_fragment_limit(self):
        _,dc=self.run_case()
        _,ref=self.run_case(conventional=True)
        self.assertLess(max(abs(a-b) for ra,rb in zip(dc,ref) for a,b in zip(ra,rb)),2e-7)

    @unittest.skipUnless(os.environ.get('SALMON_TEST_MPIEXEC'),'MPI launcher required')
    def test_full_buffer_energy_derivative(self):
        for xc in ('pbeh40','pbeh40_rvv10'):
            _,forces=self.run_case(ranks=2,xc=xc)
            plus,_=self.run_case(ranks=2,delta=.001,xc=xc)
            minus,_=self.run_case(ranks=2,delta=-.001,xc=xc)
            self.assertLess(abs(forces[0][0]+(plus-minus)/.002),2e-5)

    @unittest.skipUnless(os.environ.get('SALMON_TEST_MPIEXEC'),'MPI launcher required')
    def test_rank_parity(self):
        _,two=self.run_case(ranks=2)
        _,four=self.run_case(ranks=4)
        self.assertLess(max(abs(a-b) for ra,rb in zip(two,four) for a,b in zip(ra,rb)),2e-7)

    @unittest.skipUnless(os.environ.get('SALMON_TEST_MPIEXEC'),'MPI launcher required')
    def test_periodic_atom_identity(self):
        _,base=self.run_case(ranks=2)
        _,wrapped=self.run_case(ranks=2,delta=16.)
        self.assertLess(max(abs(a-b) for ra,rb in zip(base,wrapped) for a,b in zip(ra,rb)),2e-7)

    def test_water_projector_one_fragment(self):
        base=(ROOT/'samples/pbeh40_rvv10/water.inp').read_text()
        base=base.replace('16,16,16','24,24,24')
        base=base.replace('threshold=1d-9','threshold=1d-10').replace('nscf = 300','nscf = 500')
        _,ref=self.run_input(base,conventional=True)
        inp=base.replace("theory='dft'","theory='dft'\n yn_dc='y'")
        inp=inp.replace('nstate = 4','nstate = 4\n temperature_k=300d0')
        inp+="\n&dc\n num_fragment=1,1,1\n num_rgrid_buffer=0,0,0\n nproc_rgrid_tot=1,1,1\n nstate_frag=4\n yn_dc_lcfo='n'\n yn_dc_force_diagnostic='y'\n/\n"
        _,forces=self.run_input(inp)
        self.assertLess(max(abs(a-b) for ra,rb in zip(forces,ref) for a,b in zip(ra,rb)),2e-6)

if __name__=='__main__':unittest.main()
