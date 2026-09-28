"""Adaptive support SCF: full-support limit and exact compact FFT equivalence."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE') and os.environ.get('SALMON_TEST_MPIEXEC'), 'SALMON_TEST_EXE/MPIEXEC required')
class AdaptiveSCF(unittest.TestCase):
    def run_case(self, fraction, ranks=1, fft='auto', functional='pbeh40', dc=False):
        inp=Path(__file__).with_name('dc_hydrogen.inp').read_text()
        inp=inp.replace("yn_dc='y'", "yn_dc='y'" if dc else "yn_dc='n'")
        inp=inp.replace('nproc_k=2','nproc_k=1').replace('num_kgrid=1,2,1','num_kgrid=1,1,1')
        inp=inp.replace('temperature_k=300d0','temperature_k=0d0' if dc else '')
        inp=inp.replace('nstate=4','nstate=2')
        inp=inp.replace('al=16d0,8d0,8d0','al=48d0,24d0,24d0')
        inp=inp.replace('num_rgrid=16,8,8','num_rgrid=48,24,24')
        inp=inp.replace("xc='pbeh40_rvv10'",f"xc='{functional}'")
        inp=inp.replace('hse_mlwf_maxiter=20',f"yn_hse_wannier='y'\n exx_mlwf_maxiter=100\n exx_mlwf_interval=5\n exx_mlwf_tolerance=1d-7\n exx_mlwf_norm_fraction={fraction}\n exx_local_fft='{fft}'\n pbeh_coulomb_radius=4")
        inp=inp.replace('nscf=200','nscf=1000\n alpha_mb=0.1d0').replace('threshold=1d-8','threshold=1d-10')
        if dc:
            inp=inp.replace('nstate_frag=4','nstate_frag=2').replace('nscf=1000',"nscf=1000\n method_init_wf='gauss10'")
            inp=inp[:inp.index('&atomic_coor')]+"&atomic_coor\n 'H' 11.3d0 12d0 12d0 1\n 'H' 12.7d0 12d0 12d0 1\n 'H' 35.3d0 12d0 12d0 1\n 'H' 36.7d0 12d0 12d0 1\n/\n"
            inp=inp.replace('nproc_rgrid_tot=4,1,1',f'nproc_rgrid_tot={ranks},1,1')
            inp=inp.replace('nproc_rgrid=1,1,1',f'nproc_rgrid=1,{ranks//2},1')
        else:
            inp=inp.replace('nproc_rgrid=1,1,1',f'nproc_rgrid=1,{ranks},1')
        command=[os.environ['SALMON_TEST_MPIEXEC'],'-n',str(ranks),os.environ['SALMON_TEST_EXE']]
        with tempfile.TemporaryDirectory(prefix='adaptive-scf-') as tmp:
            shutil.copy(ROOT/'testsuites/pseudo/H_rps.dat',tmp)
            run=subprocess.run(command,input=inp,cwd=tmp,text=True,capture_output=True,timeout=180,
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
            if os.environ.get('SALMON_TEST_SAVE_DIR'):
                out=Path(os.environ['SALMON_TEST_SAVE_DIR'])/f'{functional}-{dc}-{fraction}-{ranks}-{fft}'
                out.mkdir(parents=True,exist_ok=True)
                (out/'inputfile').write_text(inp);(out/'output').write_text(run.stdout+run.stderr)
            self.assertEqual(run.returncode,0,run.stdout[-3000:]+run.stderr)
            self.assertIn('end SALMON',run.stdout)
            if dc:
                self.assertLess(float(re.findall(r'DC #SCF.*diff =\s*(\S+)',run.stdout)[-1]),1e-8)
            else:self.assertIn('#GS converged',run.stdout)
            energy=float(re.search(r'Total energy \(eV\) =\s*(\S+)',next(Path(tmp).rglob('*_info.data')).read_text())[1])
            if fraction:
                rows=re.findall(r'EXX_ADAPTIVE fraction/max radius/max norm loss:\s*(\S+)\s*(\S+)\s*(\S+)',run.stdout)
                self.assertTrue(rows)
                self.assertLessEqual(max(float(r[2]) for r in rows),1-fraction+1e-10)
                if fraction<1 and fft=='auto':
                    work=re.findall(r'EXX_ADAPTIVE local/global pairs/local FFT points \(orbital group 0\):\s*(\d+)\s*(\d+)\s*(\d+)',run.stdout)
                    self.assertTrue(any(int(r[0])>0 for r in work))
            print('adaptive SCF functional/DC/fraction/ranks/FFT/energy:',functional,dc,fraction,ranks,fft,energy)
            return energy

    def test_full_support_limit(self):
        a=self.run_case(0)
        b=self.run_case(1)
        self.assertLess(abs(a-b),2e-6)

    def test_adaptive_compact_and_rank_parity(self):
        for functional in ('pbeh40','hse06'):
            compact=self.run_case(.999,functional=functional)
            full=self.run_case(.999,fft='off',functional=functional)
            parallel=self.run_case(.999,ranks=2,functional=functional)
            self.assertLess(abs(compact-full),2e-6)
            self.assertLess(abs(compact-parallel),2e-6)

    def test_dc_fragment_rank_parity(self):
        for functional in ('pbeh40','pbeh40_rvv10'):
            a=self.run_case(.999,ranks=2,dc=True,functional=functional)
            b=self.run_case(.999,ranks=4,dc=True,functional=functional)
            self.assertLess(abs(a-b),2e-6)

    def test_tighter_support(self):
        ref=self.run_case(1)
        loose=self.run_case(.999)
        tight=self.run_case(.9999)
        print('support energy errors eV:',loose-ref,tight-ref)
        self.assertLess(abs(tight-ref),abs(loose-ref))

if __name__=='__main__':unittest.main()
