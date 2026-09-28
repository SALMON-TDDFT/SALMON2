"""Native mesh RT with ACE trained on masked MLWF support vectors."""
import os,re,unittest
import numpy as np
import test_ehrenfest as helpers

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE') and os.environ.get('SALMON_TEST_MPIEXEC'),'MPI executable required')
class SourceACE(unittest.TestCase):
    functional='pbeh40'
    setUpClass=classmethod(helpers.RealspaceEhrenfest.setUpClass.__func__)
    execute=classmethod(helpers.RealspaceEhrenfest.execute.__func__)
    rt_input=helpers.RealspaceEhrenfest.rt_input

    def test_mesh_rt_and_layouts(self):
        results={}
        for support,fraction,ranks,orbitals in [('occupied',1,1,1),('source',1,1,1),
                ('occupied',.999,1,1),('source',.999,1,1),('source',.999,2,1),('source',.999,2,2)]:
            inp=self.rt_input(nt=16,dt=.02,impulse=1e-4,moving=False)
            inp=inp.replace(f"xc='{self.functional}'",f"xc='{self.functional}'\n exx_ace_support='{support}'\n exx_mlwf_norm_fraction={fraction}")
            inp=inp.replace('nproc_ob=1',f'nproc_ob={orbitals}').replace('nproc_rgrid=1,1,1',f'nproc_rgrid=1,{ranks//orbitals},1')
            folder,run=self.execute(f'source_ace_{support}_{fraction}_{ranks}_{orbitals}',inp,ranks=ranks,rt=True)
            self.assertEqual(run.returncode,0,run.stdout[-4000:]+run.stderr)
            self.assertIn('end complex DC-LCFO wavefunction reconstruction',run.stdout)
            self.assertNotIn('Native LCFO RT active',run.stdout)
            if support=='source':
                self.assertEqual(run.stdout.count('EXX_SUPPORT_ACE accepted: T'),33)
                self.assertNotIn('EXX_SUPPORT_ACE accepted: F',run.stdout)
            c=np.loadtxt(next(folder.glob('*_rt.data')));e=np.loadtxt(next(folder.glob('*_rt_energy.data')))
            self.assertEqual(len(c),16);self.assertEqual(len(e),17)
            self.assertTrue(np.isfinite(c).all() and np.isfinite(e).all())
            results[support,fraction,ranks,orbitals]=(c,e)
        for a,b in zip(results['source',1,1,1],results['occupied',1,1,1]):np.testing.assert_allclose(a,b,atol=2e-8,rtol=0)
        for ranks,orbitals in [(2,1),(2,2)]:
            for a,b in zip(results['source',.999,ranks,orbitals],results['source',.999,1,1]):np.testing.assert_allclose(a,b,atol=2e-8,rtol=0)
        f,s=results['occupied',.999,1,1],results['source',.999,1,1]
        print('SOURCE_ACE current difference / energy difference / energy width',
            np.max(abs(f[0][:,13:16]-s[0][:,13:16])),np.max(abs(f[1][:,1]-s[1][:,1])),np.ptp(s[1][1:,1]))

    def test_source_mode_guards(self):
        for tag,extra,moving,message in [
                ('invalid',"exx_ace_support='invalid'",False,'exx_ace_support must be occupied or source'),
                ('moving',"exx_ace_support='source'\n exx_mlwf_norm_fraction=.999",True,'source-support ACE requires fixed-ion native RT'),
                ('screen',"exx_ace_support='source'\n exx_mlwf_norm_fraction=.999\n exx_pair_screening='on'",False,'source-support ACE uses exact support pairs'),
                ('missing',"exx_ace_support='source'",False,'source-support ACE requires adaptive support')]:
            inp=self.rt_input(nt=1,moving=moving).replace('&functional','&functional\n '+extra)
            _,run=self.execute('source_ace_guard_'+tag,inp,rt=True)
            self.assertNotEqual(run.returncode,0);self.assertIn(message,run.stdout+run.stderr)

class HSESourceACE(SourceACE):
    functional='hse06'

if __name__=='__main__':unittest.main()
