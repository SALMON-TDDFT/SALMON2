"""Fault-inject support-ACE rejection in a temporary executable, never production."""
import os,re,shlex,subprocess,tempfile,unittest
from pathlib import Path
import numpy as np
import test_ehrenfest as helpers

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE') and os.environ.get('SALMON_TEST_MPIEXEC'),'built MPI executable required')
class SupportFallback(unittest.TestCase):
    functional='pbeh40'
    execute=classmethod(helpers.RealspaceEhrenfest.execute.__func__)
    rt_input=helpers.RealspaceEhrenfest.rt_input

    @classmethod
    def setUpClass(cls):
        helpers.RealspaceEhrenfest.setUpClass.__func__(cls)
        cls.original=os.environ['SALMON_TEST_EXE']
        build=Path(cls.original).resolve().parent
        cls.injected=cls.root/'salmon-reject-support'
        source=(helpers.ROOT/'src/xc/hse_native.f90').read_text()
        start=source.index('    subroutine build_source_support_ace(')
        pos=source.index('      if(adaptive_bad/=0)return',start)
        source=source[:pos]+"      adaptive_bad=1 ! test-only rejection after exchange/ACE work\n"+source[pos:]
        probe=cls.root/'hse_native.f90';probe.write_text(source);obj=cls.root/'hse_native.o'
        subprocess.run(['mpifort','-O3','-fopenmp','-cpp','-ffree-line-length-none','-fallow-argument-mismatch',
            '-I'+str(build),'-J'+str(cls.root),'-c',str(probe),'-o',str(obj)],check=True,capture_output=True)
        cmd=shlex.split((build/'src/CMakeFiles/salmon.dir/link.txt').read_text())
        index=cmd.index('CMakeFiles/salmon.dir/xc/hse_native.f90.o');cmd[index]=str(obj)
        cmd[cmd.index('-o')+1]=str(cls.injected)
        subprocess.run(cmd,cwd=build/'src',check=True,capture_output=True)

    def test_rejected_attempt_work_is_counted(self):
        for ranks,orbitals in [(1,1),(2,1),(2,2)]:
            inp=self.rt_input(nt=1,moving=False).replace('&functional','&functional\n exx_mlwf_norm_fraction=1d0')
            inp=inp.replace('nproc_ob=1',f'nproc_ob={orbitals}').replace('nproc_rgrid=1,1,1',f'nproc_rgrid=1,{ranks//orbitals},1')
            reference,base=self.execute(f'fallback_reference_{ranks}_{orbitals}',inp,ranks=ranks,rt=True)
            self.assertEqual(base.returncode,0,base.stdout[-2000:]+base.stderr)
            try:
                os.environ['SALMON_TEST_EXE']=str(self.injected)
                folder,run=self.execute(f'fallback_forced_{ranks}_{orbitals}',inp.replace('&functional',"&functional\n exx_ace_support='source'"),ranks=ranks,rt=True)
            finally:os.environ['SALMON_TEST_EXE']=self.original
            self.assertEqual(run.returncode,0,run.stdout[-2000:]+run.stderr)
            self.assertEqual(run.stdout.count('EXX_SUPPORT_ACE accepted: F'),3)
            pattern=r'EXX_ADAPTIVE local/global pairs/local FFT points \(orbital group 0\):\s+(\d+)\s+(\d+)\s+(\d+)'
            expected=np.array(re.findall(pattern,base.stdout),dtype=np.int64)*2
            actual=np.array(re.findall(pattern,run.stdout),dtype=np.int64)
            np.testing.assert_array_equal(actual,expected,err_msg='both rejected and fallback work must be counted')
            for glob in ('*_rt.data','*_rt_energy.data'):
                np.testing.assert_allclose(np.loadtxt(next(folder.glob(glob))),np.loadtxt(next(reference.glob(glob))),atol=2e-8,rtol=0)

if __name__=='__main__':unittest.main()
