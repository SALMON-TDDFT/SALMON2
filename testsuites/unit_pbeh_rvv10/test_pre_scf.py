"""PBE warmup must switch operators before accepting a hybrid ground state."""
import os
import re
import unittest
import test_exx_inputs as helpers
import test_adaptive_scf as adaptive

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE'), 'SALMON_TEST_EXE required')
class PreSCF(unittest.TestCase):
    run_case = helpers.ExxInputs.run_case

    @staticmethod
    def setup(s, functional='pbeh40'):
        return s.replace("xc='pbeh40_rvv10'",f"xc='{functional}'").replace('num_kgrid=1,2,1','num_kgrid=1,1,1').replace('nstate=4','nstate=2').replace('temperature_k=300d0','').replace('nscf=200','nscf=1000\n alpha_mb=.1d0').replace('threshold=1d-8','threshold=1d-10')

    def test_switch_and_final_energy(self):
        controls='exx_mlwf_norm_fraction=1\n exx_mlwf_maxiter=100\n exx_mlwf_tolerance=1d-7'
        pbe_histories=[]
        for functional in ('hse06','pbeh40','pbeh40_rvv10'):
            with self.subTest(functional=functional):
                transform=lambda s: self.setup(s,functional)
                direct=self.run_case(controls,transform)
                staged=self.run_case(controls+'\n exx_pre_scf_threshold=1d-4\n exx_pre_scf_steps=3',transform)
                self.assertLess(abs(direct[0]-staged[0]),2e-6)
                before,after=staged[2].split('EXX_PRE_SCF switch to target hybrid',1)
                self.assertIn('EXX_PRE_SCF PBE stage',before)
                self.assertNotIn('EXX_ADAPTIVE',before)
                self.assertNotIn('rVV10 backend:',before)
                self.assertIn('EXX_ADAPTIVE',after)
                if functional=='pbeh40_rvv10': self.assertIn('rVV10 backend:',after)
                rows=re.findall(r'EXX_PRE_SCF residual/count:\s*(\S+)\s*(\d+)',before)
                self.assertGreaterEqual(len(rows),3)
                self.assertEqual(int(rows[-1][1]),3)
                self.assertTrue(all(float(r[0])<1e-4 for r in rows[-3:]))
                pbe_histories.append(re.findall(r'iter=\s*\d+\s+Total Energy=\s*(\S+)',before))
                print('PRE_SCF',functional,'direct/staged eV',direct[0],staged[0])

        self.assertTrue(pbe_histories[0])
        self.assertEqual(pbe_histories[0],pbe_histories[1])
        self.assertEqual(pbe_histories[0],pbe_histories[2])

    def test_input_and_incomplete_guard(self):
        self.run_case('exx_pre_scf_threshold=-1',error='exx_pre_scf_threshold must be finite and nonnegative')
        self.run_case('exx_pre_scf_steps=0',error='exx_pre_scf_steps must be positive')
        self.run_case('exx_pre_scf_threshold=1d-4',lambda s:self.setup(s).replace('nscf=1000','nscf=1'),
                      error='PBE pre-SCF unfinished; no hybrid ground state')
        self.run_case('exx_pre_scf_threshold=1d-4',lambda s:self.setup(s).replace("theory='dft'","theory='dft_md'"),
                      error='PBE pre-SCF requires fresh static hybrid SCF')

    def test_requires_real_pbe_residuals(self):
        result=self.run_case('exx_pre_scf_threshold=1d10\n exx_pre_scf_steps=3',self.setup)
        switch=int(re.search(r'switch to target hybrid before iteration (\d+)',result[2])[1])
        self.assertEqual(switch,4)

    def test_auto_mixing_transition(self):
        result=self.run_case('exx_pre_scf_threshold=1d-8\n exx_pre_scf_steps=3',
            lambda s:self.setup(s).replace('&scf',"&scf\n yn_auto_mixing='y'\n update_mixing_ratio=100d0"))
        _,after=result[2].split('EXX_PRE_SCF switch to target hybrid before iteration ',1)
        switch=int(after.split()[0])
        reductions=re.findall(r'decreased from.*?at iter =\s*(\d+)',after)
        self.assertNotIn(str(switch),reductions)

    def test_last_iteration_convergence(self):
        controls='exx_pre_scf_threshold=1d-4'
        result=self.run_case(controls,self.setup)
        last=int(re.search(r'#GS converged at\s*(\d+)',result[2])[1])-1
        self.run_case(controls,lambda s:self.setup(s).replace('nscf=1000',f'nscf={last}'))

    def test_unsupported_modes(self):
        controls='exx_pre_scf_threshold=1d-4'
        for control in ('checkpoint_interval=1', 'time_shutdown=10d0'):
            self.run_case(controls,lambda s:self.setup(s).replace('&control','&control\n '+control),
                          error='PBE pre-SCF stage is not supported by checkpoints or diagnostic snapshots')
        for scf in ("convergence='norm_pot'", "method_mixing='simple_potential'"):
            self.run_case(controls,lambda s:self.setup(s).replace('&scf','&scf\n '+scf),
                          error='PBE pre-SCF requires density mixing and a density convergence metric')
        self.run_case(controls,lambda s:self.setup(s).replace("xc='pbeh40'","xc='pbe'"),
                      error='PBE pre-SCF requires fresh static hybrid SCF')

    def test_dc_localization_controls(self):
        self.run_case("yn_exx_dc_mlwf='x'",
                      error="Bad input: yn_* option only accepts 'y' or 'n'.",error_exit=False)
        for mask in ('exx_mlwf_norm_fraction=.999','exx_mlwf_radius=3'):
            self.run_case(mask,lambda s:self.setup(s).replace("yn_dc='n'","yn_dc='y'"),
                error='DC canonical exchange requires full support')

    def test_dc_requires_thermal_occupations(self):
        self.run_case('exx_pre_scf_threshold=1d-4',
            lambda s:self.setup(s).replace("yn_dc='n'","yn_dc='y'").replace('&system','&system\n temperature_k=0d0'),
            error='PBE pre-SCF in DC requires positive electronic temperature')

    def test_unfinished_hybrid_guard(self):
        self.run_case('exx_pre_scf_threshold=1d10\n exx_pre_scf_steps=1',
                      lambda s:self.setup(s).replace('nscf=1000','nscf=2'),
                      error='Hybrid SCF after PBE pre-SCF not converged; ground state rejected')

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE') and os.environ.get('SALMON_TEST_MPIEXEC'), 'MPI executable required')
class PreSCFDistributed(unittest.TestCase):
    run_case=adaptive.AdaptiveSCF.run_case

    def test_native_adaptive(self):
        for functional in ('hse06','pbeh40'):
            a=self.run_case(.999,ranks=2,functional=functional)
            b=self.run_case(.999,ranks=2,functional=functional,pre_scf=1e-4)
            self.assertLess(abs(a-b),2e-6)

    def test_dc_thermal_rvv10(self):
        a=self.run_case(1,ranks=2,functional='pbeh40_rvv10',dc=True,temperature_k=300,small_cell=True,dc_mlwf='n')
        b=self.run_case(1,ranks=2,functional='pbeh40_rvv10',dc=True,pre_scf=1e-4,temperature_k=300,small_cell=True,dc_mlwf='n')
        c=self.run_case(1,ranks=4,functional='pbeh40_rvv10',dc=True,pre_scf=1e-4,temperature_k=300,small_cell=True,dc_mlwf='n')
        self.assertLess(abs(a-b),2e-6)
        self.assertLess(abs(b-c),2e-6)
        old=self.run_case(1,ranks=2,functional='pbeh40_rvv10',dc=True,temperature_k=300,small_cell=True,dc_mlwf='y')
        self.assertLess(abs(a-old),2e-6)
        orbital=self.run_case(1,ranks=4,functional='pbeh40_rvv10',dc=True,pre_scf=1e-4,
            temperature_k=300,small_cell=True,dc_mlwf='n',orbital_ranks=2)
        self.assertLess(abs(b-orbital),2e-6)

    def test_dc_kmesh_without_localization(self):
        for kpoints,ranks in ((1,2),(2,4)):
            a=self.run_case(0,ranks=ranks,functional='pbeh40_rvv10',dc=True,pre_scf=1e-4,
                temperature_k=300,small_cell=True,dc_mlwf='n',kpoints=kpoints)
            b=self.run_case(0,ranks=ranks,functional='pbeh40_rvv10',dc=True,pre_scf=1e-4,
                temperature_k=300,small_cell=True,dc_mlwf='y',kpoints=kpoints)
            self.assertLess(abs(a-b),2e-6)

if __name__=='__main__': unittest.main()
