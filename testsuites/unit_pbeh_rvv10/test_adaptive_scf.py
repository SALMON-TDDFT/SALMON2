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
    def run_case(self, fraction, ranks=1, fft='auto', functional='pbeh40', dc=False, pre_scf=0, temperature_k=0, small_cell=False, dc_mlwf='y', orbital_ranks=1, kpoints=1, radius=0, mlwf_maxiter=100, mlwf_tolerance=1e-7):
        inp=Path(__file__).with_name('dc_hydrogen.inp').read_text()
        inp=inp.replace("yn_dc='y'", "yn_dc='y'" if dc else "yn_dc='n'")
        inp=inp.replace('nproc_k=2','nproc_k=1').replace('num_kgrid=1,2,1','num_kgrid=1,1,1')
        inp=inp.replace('temperature_k=300d0',f'temperature_k={temperature_k}' if dc else '')
        inp=inp.replace('nstate=4','nstate=2')
        inp=inp.replace('al=16d0,8d0,8d0','al=48d0,24d0,24d0')
        inp=inp.replace('num_rgrid=16,8,8','num_rgrid=48,24,24')
        inp=inp.replace("xc='pbeh40_rvv10'",f"xc='{functional}'")
        inp=inp.replace('hse_mlwf_maxiter=20',f"yn_hse_wannier='y'\n exx_mlwf_maxiter={mlwf_maxiter}\n exx_mlwf_interval=5\n exx_mlwf_tolerance={mlwf_tolerance}\n exx_mlwf_norm_fraction={fraction}\n exx_local_fft='{fft}'\n pbeh_coulomb_radius=4")
        inp=inp.replace('nscf=200','nscf=1000\n alpha_mb=0.1d0').replace('threshold=1d-8','threshold=1d-10')
        if dc:
            inp=inp.replace('&functional',f"&functional\n yn_exx_dc_mlwf='{dc_mlwf}'")
            inp=inp.replace('nstate_frag=4',f'nstate_frag={4 if temperature_k>0 else 2}').replace('nscf=1000',"nscf=1000\n method_init_wf='gauss10'")
            inp=inp[:inp.index('&atomic_coor')]+"&atomic_coor\n 'H' 11.3d0 12d0 12d0 1\n 'H' 12.7d0 12d0 12d0 1\n 'H' 35.3d0 12d0 12d0 1\n 'H' 36.7d0 12d0 12d0 1\n/\n"
            inp=inp.replace('nproc_rgrid_tot=4,1,1',f'nproc_rgrid_tot={ranks},1,1')
            inp=inp.replace('nproc_rgrid=1,1,1',f'nproc_rgrid=1,{ranks//2//orbital_ranks//kpoints},1')
            inp=inp.replace('nproc_ob=1',f'nproc_ob={orbital_ranks}').replace('nproc_k=1',f'nproc_k={kpoints}')
            inp=inp.replace('num_kgrid=1,1,1',f'num_kgrid=1,{kpoints},1')
        else:
            inp=inp.replace('nproc_rgrid=1,1,1',f'nproc_rgrid=1,{ranks},1')
        if small_cell:
            inp=inp.replace("method_init_wf='gauss10'","method_init_wf='gauss'")
            inp=inp.replace('al=48d0,24d0,24d0','al=16d0,8d0,8d0').replace('num_rgrid=48,24,24','num_rgrid=16,8,8')
            inp=inp[:inp.index('&atomic_coor')]+"&atomic_coor\n 'H' 3.3d0 4d0 4d0 1\n 'H' 4.7d0 4d0 4d0 1\n 'H' 11.3d0 4d0 4d0 1\n 'H' 12.7d0 4d0 4d0 1\n/\n"
        if radius:
            inp=inp.replace('&functional',f'&functional\n exx_mlwf_radius={radius}')
        if pre_scf:
            inp=inp.replace('&functional',f'&functional\n exx_pre_scf_threshold={pre_scf}\n exx_pre_scf_steps=3')
        command=[os.environ['SALMON_TEST_MPIEXEC'],'-n',str(ranks),os.environ['SALMON_TEST_EXE']]
        with tempfile.TemporaryDirectory(prefix='adaptive-scf-') as tmp:
            shutil.copy(ROOT/'testsuites/pseudo/H_rps.dat',tmp)
            run=subprocess.run(command,input=inp,cwd=tmp,text=True,capture_output=True,timeout=180,
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
            if os.environ.get('SALMON_TEST_SAVE_DIR'):
                out=Path(os.environ['SALMON_TEST_SAVE_DIR'])/f'{functional}-{dc}-{fraction}-{ranks}-{fft}-pre{pre_scf}-T{temperature_k}-small{small_cell}-mlwf{dc_mlwf}-ob{orbital_ranks}-k{kpoints}-R{radius}'
                out.mkdir(parents=True,exist_ok=True)
                (out/'inputfile').write_text(inp);(out/'output').write_text(run.stdout+run.stderr)
            self.assertEqual(run.returncode,0,run.stdout[-3000:]+run.stderr)
            self.assertIn('end SALMON',run.stdout)
            if dc and dc_mlwf=='n':
                self.assertIn('EXX_DC canonical full-fragment source',run.stdout)
                self.assertNotIn('EXX_SPATIAL refresh/iterations',run.stdout)
                self.assertNotIn('HSE_WANNIER refresh/iterations',run.stdout)
            if pre_scf:
                before,after=run.stdout.split('EXX_PRE_SCF switch to target hybrid',1)
                self.assertNotIn('EXX_ADAPTIVE',before)
                self.assertNotIn('rVV10 backend:',before)
                self.assertIn('EXX_DC canonical full-fragment source' if dc and dc_mlwf=='n' else
                              ('EXX_FIXED' if radius else ('EXX_ADAPTIVE' if fraction else 'HSE_WANNIER')),after)

            if dc and temperature_k>0:
                charges=re.findall(r'integral\(rho_tot\)=\s*(\S+)',run.stdout)
                self.assertTrue(charges)
                self.assertLess(max(abs(float(q)-4) for q in charges),1e-9)
            if dc:
                self.assertLess(float(re.findall(r'DC #SCF.*diff =\s*(\S+)',run.stdout)[-1]),1e-8)
            else:self.assertIn('#GS converged',run.stdout)
            energy=float(re.search(r'Total energy \(eV\) =\s*(\S+)',next(Path(tmp).rglob('*_info.data')).read_text())[1])
            if radius:
                rows=re.findall(r'EXX_FIXED radius/min retained norm/max loss:\s*(\S+)\s*(\S+)\s*(\S+)',run.stdout)
                self.assertTrue(rows)
                self.assertTrue(all(abs(float(r[0])-radius)<1e-10 for r in rows))
            if fraction and not radius and not (dc and dc_mlwf=='n'):
                rows=re.findall(r'EXX_ADAPTIVE fraction/max radius/max norm loss:\s*(\S+)\s*(\S+)\s*(\S+)',run.stdout)
                self.assertTrue(rows)
                self.assertLessEqual(max(float(r[2]) for r in rows),1-fraction+1e-10)
                if fraction<1 and fft=='auto' and not small_cell:
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
