"""Fixed-ion adaptive-support mesh RT, retaining the native Taylor/ACE route."""
import os
import re
import unittest
import numpy as np
import test_ehrenfest as helpers

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE') and os.environ.get('SALMON_TEST_MPIEXEC'), 'MPI executable required')
class AdaptiveRT(unittest.TestCase):
    functional='pbeh40'
    setUpClass=classmethod(helpers.RealspaceEhrenfest.setUpClass.__func__)
    execute=classmethod(helpers.RealspaceEhrenfest.execute.__func__)
    rt_input=helpers.RealspaceEhrenfest.rt_input

    def test_impulse_16(self):
        results={}
        for fraction,ranks,fft in [(0,1,'auto'),(1,1,'auto'),(.999,1,'auto'),(.999,2,'auto'),(.999,2,'off')]:
            inp=self.rt_input(nt=16,dt=.02,impulse=1e-4,moving=False)
            inp=inp.replace(f"xc='{self.functional}'",f"xc='{self.functional}'\n yn_hse_wannier='y'\n exx_mlwf_norm_fraction={fraction}\n exx_local_fft='{fft}'")
            inp=inp.replace('nproc_rgrid=1,1,1',f'nproc_rgrid=1,{ranks},1')
            folder,run=self.execute(f'rt_{fraction}_{ranks}_{fft}',inp,ranks=ranks,rt=True)
            self.assertEqual(run.returncode,0,run.stdout[-3000:]+run.stderr)
            self.assertIn('end SALMON',run.stdout)
            current=np.loadtxt(next(folder.glob('*_rt.data')))
            energy=np.loadtxt(next(folder.glob('*_rt_energy.data')))
            self.assertEqual(len(current),16);self.assertEqual(len(energy),17)
            self.assertTrue(np.isfinite(current).all() and np.isfinite(energy).all())
            if fraction:
                rows=re.findall(r'EXX_ADAPTIVE fraction/max radius/max norm loss:\s*(\S+)\s*(\S+)\s*(\S+)',run.stdout)
                self.assertTrue(rows)
                if fraction<1:self.assertGreater(float(rows[0][2]),0)
                self.assertLessEqual(max(float(r[2]) for r in rows),1-fraction+1e-10)
            results[fraction,ranks,fft]=(current,energy)
        for key,ref in [((1,1,'auto'),(0,1,'auto')),((.999,2,'auto'),(.999,1,'auto')),((.999,2,'off'),(.999,1,'auto'))]:
            for value,reference in zip(results[key],results[ref]):
                np.testing.assert_allclose(value,reference,atol=2e-8,rtol=0)
        c,e=results[.999,1,'auto'];cr,er=results[0,1,'auto']
        print('adaptive RT16 maximum current difference/full and post-kick energy width:',
              np.max(abs(c[:,13:16]-cr[:,13:16])),np.ptp(e[1:,1]),
              'full-support width/current peak:',np.ptp(er[1:,1]),np.max(abs(cr[:,13:16])))

class AdaptiveHSERT(AdaptiveRT):
    functional='hse06'

class AdaptiveRVV10RT(AdaptiveRT):
    functional='pbeh40_rvv10'

if __name__=='__main__':unittest.main()
