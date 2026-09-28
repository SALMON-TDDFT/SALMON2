"""Instrument live native RT buffers and compare to a frozen pre-change executable."""
import os,shlex,subprocess,unittest
from pathlib import Path
import numpy as np
import test_ehrenfest as helpers

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE') and os.environ.get('SALMON_MEMORY_REFERENCE'),'memory reference executables required')
class MemoryBuffers(unittest.TestCase):
    functional='pbeh40'
    setUpClass=classmethod(helpers.RealspaceEhrenfest.setUpClass.__func__)
    execute=classmethod(helpers.RealspaceEhrenfest.execute.__func__)
    rt_input=helpers.RealspaceEhrenfest.rt_input

    def test_midpoint_live_storage_and_equivalence(self):
        original=os.environ['SALMON_TEST_EXE'];build=Path(original).resolve().parent
        source=(helpers.ROOT/'src/xc/hse_native.f90').read_text()
        output="allocated(output_work)" if 'output_work(:,:,:)' in source else '.false.'
        source=source.replace('      taylor_midpoint=.true.',"      taylor_midpoint=.true.\n      write(*,'(a,3l2)')'MEMORY_PROBE redundant: ',allocated(midpoint_ace%factors),allocated(cached_action),"+output)
        source=source.replace('      taylor_midpoint=.true.', "      taylor_midpoint=.true.\n      write(*,'(a,3l2)')'MEMORY_PROBE packed: ',ace%packed,allocated(ace%factors),allocated(spatial%source)")
        src=self.root/'memory_probe.f90';src.write_text(source);obj=self.root/'memory_probe.o';exe=self.root/'salmon-memory-probe'
        subprocess.run(['mpifort','-O3','-fopenmp','-cpp','-ffree-line-length-none','-fallow-argument-mismatch','-I'+str(build),'-J'+str(self.root),'-c',str(src),'-o',str(obj)],check=True,capture_output=True)
        cmd=shlex.split((build/'src/CMakeFiles/salmon.dir/link.txt').read_text());cmd[cmd.index('CMakeFiles/salmon.dir/xc/hse_native.f90.o')]=str(obj);cmd[cmd.index('-o')+1]=str(exe)
        subprocess.run(cmd,cwd=build/'src',check=True,capture_output=True)
        for ranks,orbitals in [(1,1),(2,1),(2,2)]:
            inp=self.rt_input(nt=16,dt=.02,impulse=1e-4,moving=False).replace('&functional',"&functional\n exx_ace_support='source'\n exx_mlwf_norm_fraction=.999")
            inp=inp.replace('nproc_ob=1',f'nproc_ob={orbitals}').replace('nproc_rgrid=1,1,1',f'nproc_rgrid=1,{ranks//orbitals},1')
            results=[]
            try:
                for tag,binary in [('reference',os.environ['SALMON_MEMORY_REFERENCE']),('new',str(exe))]:
                    os.environ['SALMON_TEST_EXE']=binary
                    folder,run=self.execute(f'memory_{tag}_{ranks}_{orbitals}',inp,ranks=ranks,rt=True)
                    self.assertEqual(run.returncode,0,run.stdout[-2500:]+run.stderr)
                    if tag=='new':
                        lines=[l for l in run.stdout.splitlines() if l.startswith('MEMORY_PROBE redundant:')]
                        self.assertEqual(len(lines),16*ranks)
                        self.assertTrue(all(l.endswith(' F F F') for l in lines),lines[:2])
                        packed=[l for l in run.stdout.splitlines() if l.startswith('MEMORY_PROBE packed:')]
                        self.assertEqual(len(packed),16*ranks)
                        self.assertTrue(all(l.endswith(' T F F') for l in packed),packed[:2])
                    results.append([np.loadtxt(next(folder.glob(pattern))) for pattern in ('*_rt.data','*_rt_energy.data')])
            finally:os.environ['SALMON_TEST_EXE']=original
            for a,b in zip(*results):np.testing.assert_allclose(a,b,rtol=0,atol=2e-8)

if __name__=='__main__':unittest.main()
