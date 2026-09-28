"""Same-seed native HSE RT: pair diagnosis, optional omission and ACE fallback."""
import os
from pathlib import Path
import re
import unittest
import numpy as np
import test_ehrenfest as helpers
import test_exx_inputs as input_helpers

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE') and os.environ.get('SALMON_TEST_MPIEXEC'), 'MPI executable required')
class PairScreenRT(unittest.TestCase):
    functional='hse06'
    setUpClass=classmethod(helpers.RealspaceEhrenfest.setUpClass.__func__)
    execute=classmethod(helpers.RealspaceEhrenfest.execute.__func__)
    rt_input=helpers.RealspaceEhrenfest.rt_input

    def test_initial_localization_guard(self):
        inp=self.rt_input(nt=1,dt=.02,impulse=1e-4,moving=False)
        inp=inp.replace(f"xc='{self.functional}'",f"xc='{self.functional}'\n exx_mlwf_norm_fraction=.999\n exx_mlwf_tolerance=1d-30")
        inp=inp.replace('exx_mlwf_maxiter=100','exx_mlwf_maxiter=1')
        folder,run=self.execute('pair_initial_localization_guard',inp,rt=True)
        self.assertNotEqual(run.returncode,0)
        self.assertIn('Adaptive RT requires an accepted transported MLWF gauge',run.stdout+run.stderr)

    def test_modes_and_fallback(self):
        results={}
        for fraction,mode,tol,ranks in [(1,'off',0,1),(1,'diagnose',1,1),(1,'on',1,1),
                (1,'on',1,2),(.999,'off',0,1),(.999,'diagnose',1,1),(.999,'on',1,1)]:
            inp=self.rt_input(nt=16,dt=.02,impulse=1e-4,moving=False)
            inp=inp.replace(f"xc='{self.functional}'",f"xc='{self.functional}'\n exx_mlwf_norm_fraction={fraction}\n exx_pair_screening='{mode}'\n exx_pair_tolerance={tol}")
            inp=inp.replace('nproc_rgrid=1,1,1',f'nproc_rgrid=1,{ranks},1')
            folder,run=self.execute(f'pair_{fraction}_{mode}_{ranks}',inp,ranks=ranks,rt=True)
            self.assertEqual(run.returncode,0,run.stdout[-4000:]+run.stderr)
            self.assertIn('end SALMON',run.stdout)
            c=np.loadtxt(next(folder.glob('*_rt.data')));e=np.loadtxt(next(folder.glob('*_rt_energy.data')))
            self.assertEqual(len(c),16);self.assertEqual(len(e),17)
            self.assertTrue(np.isfinite(c).all() and np.isfinite(e).all())
            lines=re.findall(r'EXX_PAIR mode/candidates/skipped/action bound/max rank CPU seconds:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)',run.stdout)
            if mode!='off':
                self.assertEqual(len(lines),33)
                self.assertTrue(all(0<=float(v[3])<=tol*(1+1e-10) for v in lines))
                if mode=='diagnose':self.assertTrue(all(int(v[2])==0 for v in lines))
                else:
                    self.assertGreater(sum(int(v[2]) for v in lines),0,'test must exercise omission')
                    accepted=re.findall(r'EXX_PAIR unscreened ACE fallback: ([TF]) accepted action bound:\s+(\S+)',run.stdout)
                    self.assertEqual(len(accepted),33)
                    self.assertTrue(all(0<=float(b)<=tol*(1+1e-10) for flag,b in accepted))
                    self.assertTrue(any(flag=='F' and float(b)>0 for flag,b in accepted),'must accept finite screening')
            results[fraction,mode,ranks]=(c,e)
            print('PAIR_RT',fraction,mode,ranks,'energy width',np.ptp(e[1:,1]),
                'attempted omissions',sum(int(v[2]) for v in lines),
                'fallbacks',run.stdout.count('EXX_PAIR unscreened ACE fallback: T'))
        for fraction in (1,.999):
            for a,b in zip(results[fraction,'diagnose',1],results[fraction,'off',1]):
                np.testing.assert_allclose(a,b,atol=2e-8,rtol=0)
        # The deliberately loose stress tolerance is not a production default.
        # Still verify this fixture's observed short-time energy continuity.
        self.assertLess(np.ptp(results[1,'on',1][1][1:,1]),1e-7)
        self.assertLess(np.ptp(results[.999,'on',1][1][1:,1]),np.ptp(results[.999,'off',1][1][1:,1])+1e-7)
        for a,b in zip(results[1,'on',2],results[1,'on',1]):
            np.testing.assert_allclose(a,b,atol=2e-8,rtol=0)


class PBEhPairScreenRT(PairScreenRT):
    functional='pbeh40'

    def test_orbital_split(self):
        values=[]
        for orbitals in (1,2):
            inp=self.rt_input(nt=4,dt=.02,impulse=1e-4,moving=False)
            inp=inp.replace("xc='pbeh40'", "xc='pbeh40'\n exx_mlwf_norm_fraction=.999\n exx_pair_screening='on'\n exx_pair_tolerance=1d-6")
            inp=inp.replace('nproc_ob=1',f'nproc_ob={orbitals}')
            inp=inp.replace('nproc_rgrid=1,1,1',f'nproc_rgrid=1,{2//orbitals},1')
            folder,run=self.execute(f'pbeh_pair_orbitals_{orbitals}',inp,ranks=2,rt=True)
            self.assertEqual(run.returncode,0,run.stdout[-4000:]+run.stderr)
            self.assertIn('EXX_PAIR generated grid products/catalogue entries:',run.stdout)
            values.append((np.loadtxt(next(folder.glob('*_rt.data'))),np.loadtxt(next(folder.glob('*_rt_energy.data')))))
        for a,b in zip(*values):np.testing.assert_allclose(a,b,atol=2e-8,rtol=0)

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE'),'SALMON_TEST_EXE required')
class PairScreenSCF(unittest.TestCase):
    run_case=input_helpers.ExxInputs.run_case

    def test_screened_scf(self):
        def setup(s):
            return s.replace("xc='pbeh40_rvv10'","xc='hse06'").replace('num_kgrid=1,2,1','num_kgrid=1,1,1').replace('nstate=4','nstate=2').replace('temperature_k=300d0','').replace('nscf=200','nscf=1000\n alpha_mb=.1d0').replace('threshold=1d-8','threshold=1d-10')
        controls="exx_mlwf_norm_fraction=1\n exx_mlwf_maxiter=100\n exx_mlwf_tolerance=1d-7"
        results={}
        for mode,tol in [('off',0),('diagnose',1),('on',1),('on',1e-6)]:
            result=self.run_case(controls+f"\n exx_pair_screening='{mode}'\n exx_pair_tolerance={tol}",setup)
            results[mode,tol]=result
            if mode!='off':self.assertIn('EXX_PAIR',result[2])
            if os.environ.get('SALMON_TEST_SAVE_DIR'):
                folder=Path(os.environ['SALMON_TEST_SAVE_DIR'])/f'pair_scf_{mode}_{tol}'
                folder.mkdir(parents=True,exist_ok=True)
                inp=(input_helpers.ROOT/'testsuites/unit_pbeh_rvv10/dc_hydrogen.inp').read_text()
                inp=inp.replace("yn_dc='y'","yn_dc='n'").replace('nproc_k=2','nproc_k=1')
                inp=setup(inp.replace('hse_mlwf_maxiter=20',controls+f"\n exx_pair_screening='{mode}'\n exx_pair_tolerance={tol}"))
                (folder/'inputfile').write_text(inp);(folder/'output').write_text(result[2])
                (folder/'variables.log').write_text(result[1])
            print('PAIR_SCF',mode,tol,'energy eV',result[0])
        self.assertLess(abs(results['off',0][0]-results['diagnose',1][0]),2e-6)
        self.assertLess(abs(results['off',0][0]-results['on',1e-6][0]),2e-6)

if __name__=='__main__':unittest.main()
