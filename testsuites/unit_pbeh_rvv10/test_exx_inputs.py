"""Canonical EXX names, legacy aliases and finite-support input contract."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE'),'SALMON_TEST_EXE required')
class ExxInputs(unittest.TestCase):
    def run_case(self, functional, transform=lambda s:s, error=None):
        inp=(ROOT/'testsuites/unit_pbeh_rvv10/dc_hydrogen.inp').read_text()
        inp=inp.replace("yn_dc='y'","yn_dc='n'").replace('nproc_k=2','nproc_k=1')
        inp=inp.replace('hse_mlwf_maxiter=20',functional)
        inp=transform(inp)
        with tempfile.TemporaryDirectory() as tmp:
            shutil.copy(ROOT/'testsuites/pseudo/H_rps.dat',tmp)
            run=subprocess.run([os.environ['SALMON_TEST_EXE']],input=inp,cwd=tmp,text=True,
                capture_output=True,timeout=120,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMPI_MCA_btl='self,vader'))
            if error:
                self.assertNotEqual(run.returncode,0)
                self.assertIn(error,run.stdout+run.stderr)
                return
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            self.assertIn('#GS converged',run.stdout)
            info=next(Path(tmp).rglob('*_info.data')).read_text()
            energy=float(re.search(r'Total energy \(eV\) =\s*([\d.E+-]+)',info)[1])
            return energy,(Path(tmp)/'variables.log').read_text(),run.stdout

    def test_alias_and_matching_dual_values(self):
        old='hse_mlwf_interval=5\n hse_mlwf_maxiter=100\n hse_mlwf_tolerance=1d-7'
        new=old.replace('hse_','exx_')
        a=self.run_case(old);b=self.run_case(new);c=self.run_case(old+'\n'+new)
        self.assertAlmostEqual(a[0],b[0],places=11)
        self.assertAlmostEqual(a[0],c[0],places=11)
        for key,value in [('interval',5),('maxiter',100),('tolerance',1e-7)]:
            self.assertEqual(float(re.search(r'# exx_mlwf_'+key+r'=\s*([\d.E+-]+)',b[1])[1]),value)

    def test_conflicting_aliases(self):
        for key,a,b in [('interval','5','6'),('maxiter','100','101'),('tolerance','1d-7','2d-7')]:
            self.run_case(f'exx_mlwf_{key}={a}\n hse_mlwf_{key}={b}',error='conflicting EXX/legacy MLWF '+key)

    def test_radius_invalid(self):
        self.run_case('exx_mlwf_radius=-1',error='exx_mlwf_radius must be finite and nonnegative')

    def test_radius_md_guard(self):
        self.run_case('exx_mlwf_radius=2',lambda s:s.replace("theory='dft'","theory='dft_md'"),
            error='finite EXX MLWF radius supports static DFT only')

    def test_invalid_mlwf_controls(self):
        for control in ['exx_mlwf_interval=0','exx_mlwf_maxiter=0','exx_mlwf_tolerance=0']:
            self.run_case(control,error='invalid localization controls')

    def test_radius_restart_and_snapshot_guards(self):
        self.run_case('exx_mlwf_radius=2',lambda s:s.replace("sysname='H_dc_hse'", "sysname='H_dc_hse'\n yn_restart='y'"),
            error='restart/snapshot metadata unsupported')
        self.run_case("exx_mlwf_radius=2\n yn_hse_wannier_snapshot='y'",error='restart/snapshot metadata unsupported')

    def test_radius_length_units(self):
        def angstrom_input(s):
            s=s.replace("unit_system='a.u.'", "unit_system='A_eV_fs'")
            # Keep the same physical cell/positions after changing input units.
            factor=.529177210903
            s=re.sub(r'al=16d0,8d0,8d0', 'al='+','.join(str(x*factor) for x in [16,8,8]), s)
            s=re.sub(r"'H' ([0-9.]+)d0 4d0 4d0 1",lambda m:f"'H' {float(m[1])*factor} {4*factor} {4*factor} 1",s)
            return s
        result=self.run_case('exx_mlwf_maxiter=20\n exx_mlwf_radius=52.9177210903',angstrom_input)
        radius=float(re.search(r'# exx_mlwf_radius .*?=\s*([\d.E+-]+)',result[1])[1])
        self.assertAlmostEqual(radius,100.,places=6)

    def test_finite_radius_native_scf(self):
        controls="yn_hse_wannier='y'\n exx_mlwf_interval=5\n exx_mlwf_maxiter=100\n exx_mlwf_tolerance=1d-7"
        for xc in ['pbeh40_rvv10','hse06']:
            transform=lambda s:s.replace("xc='pbeh40_rvv10'",f"xc='{xc}'")
            full=self.run_case(controls,transform)
            cut=self.run_case(controls+'\n exx_mlwf_radius=4',transform)
            self.assertGreater(abs(full[0]-cut[0]),1e-4)
            rows=re.findall(r'EXX_MLWF radius.*?:\s*([\d.E+ -]+)',cut[2])
            self.assertTrue(rows)
            radius,protected,total_loss,max_loss=map(float,rows[-1].split())
            self.assertEqual(radius,4.)
            self.assertGreater(total_loss,0.)
            self.assertLess(total_loss,.02)
            self.assertGreaterEqual(max_loss,total_loss-1e-9)

    def test_zero_and_large_radius(self):
        base=self.run_case('exx_mlwf_maxiter=20\n exx_mlwf_radius=0')
        large=self.run_case('exx_mlwf_maxiter=20\n exx_mlwf_radius=100')
        self.assertAlmostEqual(base[0],large[0],places=11)
        self.assertIn('EXX_MLWF radius',large[2])

if __name__=='__main__':unittest.main()
