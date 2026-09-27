"""HSE spatial SCF and DC-initialized fixed-ion mesh RT parity."""
import os
import re
import unittest
import numpy as np
import test_ehrenfest as helpers

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE') and os.environ.get('SALMON_TEST_MPIEXEC'),
                     'SALMON_TEST_EXE and SALMON_TEST_MPIEXEC required')
class HSESpatial(unittest.TestCase):
    functional='hse06'
    setUpClass=classmethod(helpers.RealspaceEhrenfest.setUpClass.__func__)
    execute=classmethod(helpers.RealspaceEhrenfest.execute.__func__)
    rt_input=helpers.RealspaceEhrenfest.rt_input
    pulse_input=helpers.RealspaceEhrenfest.pulse_input

    def test_serial_spatial_parity(self):
        for field in ('impulse','pulse'):
            results=[]
            for ranks,layout in [(1,'1,1,1'),(2,'1,2,1'),(4,'1,2,2')]:
                inp=(self.rt_input(nt=20,moving=False) if field=='impulse' else
                     self.pulse_input(dt=.08,nt=120,moving=False))
                if ranks==1:
                    inp=inp.replace("xc='hse06'","xc='hse06'\n yn_hse_wannier='y'")
                inp=inp.replace('nproc_rgrid=1,1,1','nproc_rgrid='+layout)
                folder,run=self.execute(field+str(ranks),inp,ranks=ranks,rt=True)
                self.assertEqual(run.returncode,0,run.stdout[-3000:]+run.stderr)
                self.assertIn('end SALMON',run.stdout)
                self.assertIn('EXX_SPATIAL' if ranks>1 else 'HSE_WANNIER',run.stdout)
                self.assertNotIn('Native LCFO RT active',run.stdout)
                energy=np.loadtxt(next(folder.glob('*_rt_energy.data')))
                data=np.loadtxt(next(folder.glob('*_rt.data')))
                self.assertTrue(np.isfinite(energy).all() and np.isfinite(data).all())
                results.append((energy,data))
                for actual,reference in zip(results[-1],results[0]):
                    self.assertLess(np.max(abs(actual-reference)),2e-8)
                print('HSE parity field/ranks/energy/current:',field,ranks,
                      np.max(abs(energy-results[0][0])),np.max(abs(data-results[0][1])))

    def test_hse_md_still_rejected(self):
        _,run=self.execute('md_rejected',self.rt_input(moving=True),rt=True)
        self.assertNotEqual(run.returncode,0)
        self.assertIn('HSE06 requires fixed-ion',run.stdout+run.stderr)

    def test_legacy_full_propagator_admission(self):
        inp=self.rt_input(nt=2,moving=False)+"\n&propagation\n propagator='hse_taylor4_full'\n/\n"
        _,run=self.execute('legacy_full',inp,rt=True)
        # This noncubic fixture reaches the legacy kernel's existing grid guard.
        # It must not be rejected earlier by the new Wannier propagator admission.
        self.assertIn('end complex DC-LCFO wavefunction reconstruction',run.stdout)
        self.assertIn('HSE06: cubic grid and full cubic k mesh required',run.stdout+run.stderr)
        self.assertNotIn('HSE Wannier RT:',run.stdout+run.stderr)

    def scf_input(self):
        return (self.base.replace("yn_dc='y'", "yn_dc='n'")
                .replace(' nstate=4', ' nstate=2')
                .replace(' temperature_k=300d0\n','')
                .replace("xc='hse06'", "xc='hse06'\n yn_hse_wannier='y'"))

    def test_scf_spatial_parity(self):
        results=[]
        for ranks,layout in [(1,'1,1,1'),(2,'1,2,1'),(4,'1,2,2')]:
            inp=self.scf_input().replace('nproc_rgrid=1,1,1','nproc_rgrid='+layout)
            if ranks>1:inp=inp.replace(" yn_hse_wannier='y'\n",'')
            folder,run=self.execute('scf'+str(ranks),inp,ranks=ranks)
            self.assertEqual(run.returncode,0,run.stdout[-3000:]+run.stderr)
            self.assertIn('end SALMON',run.stdout)
            self.assertIn('EXX_SPATIAL' if ranks>1 else 'HSE_WANNIER',run.stdout)

            match=re.search(r'#GS converged at\s+(\d+)\s+:\s+(\S+)',run.stdout)
            self.assertIsNotNone(match,run.stdout[-3000:])
            self.assertLess(float(match[2]),1e-10)
            energy=float(re.search(r'Total energy \(eV\) =\s*(\S+)',
                                  next(folder.glob('*_info.data')).read_text())[1])
            eigen=np.loadtxt(next(folder.glob('*_eigen.data')),skiprows=4)[:,1]
            self.assertTrue(np.isfinite(energy) and np.isfinite(eigen).all())
            results.append((energy,eigen))
            self.assertLess(abs(energy-results[0][0]),1e-7)
            self.assertLess(np.max(abs(eigen-results[0][1])),1e-7)
            print('HSE SCF ranks/iterations/energy eV/eigen error Ha:',ranks,int(match[1]),
                  energy,np.max(abs(eigen-results[0][1])))

    def test_scf_spatial_guards(self):
        base=self.scf_input().replace('nproc_rgrid=1,1,1','nproc_rgrid=1,2,1')
        cases=[
            ('gs_write',"write_gs_restart_data='no'","write_gs_restart_data='all'",'write_gs_restart_data=no'),
            ('self_checkpoint',"sysname='H_dc_hse'","sysname='H_dc_hse'\n yn_self_checkpoint='y'",'yn_self_checkpoint=n'),
            ('eigen_diagnostic',"xc='hse06'","xc='hse06'\n yn_hse_eigen_diagnostic='y'",'diagnostics unsupported'),
            ('empty_states',' nstate=2',' nstate=3','occupied spin pairs'),
            ('temperature',' nstate=2',' nstate=2\n temperature_k=300d0','occupied spin pairs'),
            ('radius',"xc='hse06'","xc='hse06'\n exx_mlwf_radius=1d0",'full support'),
            ('restart',"sysname='H_dc_hse'","sysname='H_dc_hse'\n yn_restart='y'",'restart/snapshot'),
            ('snapshot',"xc='hse06'","xc='hse06'\n yn_hse_wannier_snapshot='y'",'restart/snapshot'),
            ('kpoints','num_kgrid=1,1,1','num_kgrid=1,2,1','Gamma y/z pencils'),
            ('xsplit','nproc_rgrid=1,2,1','nproc_rgrid=2,1,1','Gamma y/z pencils'),
        ]
        for name,old,new,message in cases:
            with self.subTest(name=name):
                _,run=self.execute('scf_guard_'+name,base.replace(old,new),ranks=2)
                self.assertNotEqual(run.returncode,0)
                self.assertIn(message,run.stdout+run.stderr)
