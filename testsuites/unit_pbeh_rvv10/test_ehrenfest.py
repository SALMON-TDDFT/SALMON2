"""DC preparation -> mesh orbitals -> native PBEh Ehrenfest, no LCFO projection."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]

@unittest.skipUnless(os.environ.get('SALMON_TEST_EXE') and os.environ.get('SALMON_TEST_MPIEXEC'),
                     'SALMON_TEST_EXE and SALMON_TEST_MPIEXEC required')
class RealspaceEhrenfest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='realspace-ehrenfest-')
        cls.root=Path(cls.temp.name)
        cls.addClassCleanup(cls.temp.cleanup)
        cls.base=(ROOT/'testsuites/422_H_dcdft_hse/inputfile').read_text()
        cls.base=cls.base.replace("xc='hse06'","xc='pbeh40_rvv10'\n exx_mlwf_interval=5\n exx_mlwf_maxiter=100")
        cls.base=cls.base.replace('nproc_k=2','nproc_k=1').replace('nproc_rgrid_tot=4,1,1','nproc_rgrid_tot=2,1,1')
        cls.base=cls.base.replace('num_kgrid=1,2,1','num_kgrid=1,1,1')
        cls.base=cls.base.replace('lmax_ps(1)=0','lmax_ps(1)=1').replace('lloc_ps(1)=0','lloc_ps(1)=1')
        cls.base=cls.base.replace('threshold=1d-8','threshold=1d-10').replace('nscf=200','nscf=500')
        folder,run=cls.execute('gs',cls.base,ranks=2)
        if run.returncode or 'end SALMON' not in run.stdout:
            raise AssertionError(run.stdout[-4000:]+run.stderr)
        differences=re.findall(r'DC #SCF.*diff =\s*(\S+)',run.stdout)
        if not differences or float(differences[-1])>=1e-10:raise AssertionError('DC GS did not converge')

    @classmethod
    def execute(cls,name,inp,ranks=1,rt=False,velocity=None):
        folder=cls.root/name;folder.mkdir()
        for atom in ('H','O'):shutil.copy(ROOT/f'testsuites/pseudo/{atom}_rps.dat',folder)
        (folder/'inputfile').write_text(inp)
        (folder/'velocity.dat').write_text(velocity or '0 0 0\n'*4)
        if rt:(folder/'data_dcdft').symlink_to(cls.root/(rt if isinstance(rt,str) else 'gs')/'data_dcdft',target_is_directory=True)
        cmd=[os.environ['SALMON_TEST_EXE']]
        if ranks>1:cmd=[os.environ['SALMON_TEST_MPIEXEC'],'-n',str(ranks)]+cmd
        run=subprocess.run(cmd,input=inp,cwd=folder,text=True,capture_output=True,timeout=180,
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
        (folder/'output').write_text(run.stdout+run.stderr)
        if os.environ.get('SALMON_TEST_SAVE_DIR'):
            dest=Path(os.environ['SALMON_TEST_SAVE_DIR'])/cls.root.name/name
            shutil.copytree(folder,dest,symlinks=True,dirs_exist_ok=True)
        return folder,run

    def rt_input(self,dt=.04,nt=4,impulse=.01,moving=True):
        inp=self.base.replace("theory='dft'","theory='tddft_response'\n yn_md='"+('y' if moving else 'n')+"'")
        inp=inp.replace("yn_dc='y'","yn_dc='n'\n yn_conventional_from_dcdft='y'")
        inp=inp.replace(' nstate=4',' nstate=2').replace(' temperature_k=300d0\n','')
        inp+=f"""
&tgrid
 dt={dt}
 nt={nt}
/
&emfield
 ae_shape1='impulse'
 e_impulse={impulse}
 epdir_re1=1,0,0
/
&md
 ensemble='NVE'
 file_ini_velocity='velocity.dat'
 step_update_ps=1
/
&analysis
 out_rt_energy_step=1
 out_rvf_rt_step=1
 nenergy=10
 de=.01
/
"""
        return inp

    def run_rt(self,name,**kwargs):
        folder,run=self.execute(name,self.rt_input(**kwargs),rt=True)
        self.assertEqual(run.returncode,0,run.stdout[-2500:]+run.stderr)
        self.assertIn('end SALMON',run.stdout)
        self.assertIn('end complex DC-LCFO wavefunction reconstruction',run.stdout)
        self.assertNotIn('Native LCFO RT active',run.stdout)
        self.assertIn('HSE_WANNIER',run.stdout)
        energy=np.loadtxt(next(folder.glob('*_rt_energy.data')))
        data=np.loadtxt(next(folder.glob('*_rt.data')))
        self.assertTrue(np.isfinite(energy).all())
        self.assertTrue(np.isfinite(data).all())
        charges=re.findall(r'^\s*\d+\s+[\d.]+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s*$',run.stdout,re.M)
        self.assertTrue(charges)
        self.assertLess(max(abs(float(x)-4.) for x in charges),1e-7)
        return folder,energy,data

    def test_realspace_ionic_step(self):
        folder,energy,_=self.run_rt('moving')
        xyz=next(folder.glob('*_trj.xyz')).read_text().splitlines()
        frames=[]
        for i,line in enumerate(xyz):
            if line.strip()=='4':frames.append(np.array([[float(x) for x in row.split()[1:4]] for row in xyz[i+2:i+6]]))
        self.assertGreater(len(frames),2)
        self.assertGreater(np.max(np.abs(frames[-1][:,:3]-frames[0][:,:3])),1e-9)
        self.assertGreater(np.max(energy[:,3]),0.)

    def test_timestep_refinement(self):
        drifts=[];currents=[];coordinates=[]
        for dt,nt in ((.08,20),(.04,40),(.02,80)):
            folder,energy,data=self.run_rt('refine_'+str(nt),dt=dt,nt=nt)
            total=energy[1:,1]+energy[1:,3]  # after the electronic impulse
            drifts.append(float(np.max(np.abs(total-total[0]))))
            currents.append(data[-1,13:16])
            lines=next(folder.glob('*_trj.xyz')).read_text().splitlines()[-4:]
            coordinates.append(np.array([[float(x) for x in line.split()[1:4]] for line in lines]))
        ratio=np.linalg.norm(currents[0]-currents[1])/np.linalg.norm(currents[1]-currents[2])
        print('realspace Ehrenfest drift Ha:',drifts,'current refinement ratio:',ratio)
        self.assertGreater(ratio,3.4)
        self.assertLess(ratio,4.6)
        self.assertLess(drifts[-1],1e-7)
        self.assertLess(drifts[-1],drifts[0]/3)
        self.assertLess(np.max(abs(coordinates[-1]-coordinates[-2])),1e-6)

    def test_initial_force_uses_excited_state(self):
        forces=[]
        for name,impulse in [('initial_zero',0.),('initial_kick',.2)]:
            folder,_,_=self.run_rt(name,impulse=impulse,nt=1)
            lines=next(folder.glob('*_trj.xyz')).read_text().splitlines()[2:6]
            forces.append(np.array([[float(x) for x in line.split('#f=')[1].split()] for line in lines]))
        self.assertGreater(np.max(abs(forces[1]-forces[0])),1e-6,
            'initial nonlocal force must use the post-impulse vector potential')

    def test_frozen_mesh_force_derivative(self):
        values=[];initial_force=None;step=1e-4
        for label,delta in [('minus',-step),('center',0.),('plus',step)]:
            inp=self.rt_input(impulse=0.,nt=1).replace('3.3d0',f'{3.3+delta:.12f}')
            folder,run=self.execute('frozen_'+label,inp,rt=True)
            self.assertEqual(run.returncode,0,run.stdout[-2000:]+run.stderr)
            values.append(np.loadtxt(next(folder.glob('*_rt_energy.data')))[0,1])
            if delta==0:
                line=next(folder.glob('*_trj.xyz')).read_text().splitlines()[2]
                initial_force=float(line.split('#f=')[1].split()[0])
        difference=abs(initial_force+(values[2]-values[0])/(2*step))
        print('frozen mesh force difference Ha/bohr:',difference)
        self.assertLess(difference,2e-6)

    def test_initial_ionic_kinetic_energy(self):
        folder,run=self.execute('initial_velocity',self.rt_input(nt=1),rt=True,velocity='0.001 0 0\n'*4)
        self.assertEqual(run.returncode,0,run.stdout[-2000:]+run.stderr)
        energy=np.loadtxt(next(folder.glob('*_rt_energy.data')))
        self.assertGreater(energy[0,3],1e-3,'initial energy record must include supplied ionic velocities')

    def test_water_mesh_propagation(self):
        gs=(ROOT/'samples/pbeh40_rvv10/water.inp').read_text()
        gs=gs.replace('16,16,16','24,24,24').replace('nscf = 300','nscf = 500')
        gs=gs.replace('threshold=1d-9','threshold=1d-10')
        gs=gs.replace("theory='dft'","theory='dft'\n yn_dc='y'")
        gs=gs.replace('nstate = 4','nstate = 4\n temperature_k=300d0')
        gs+="\n&dc\n num_fragment=1,1,1\n num_rgrid_buffer=0,0,0\n nstate_frag=4\n yn_dc_lcfo='y'\n/\n"
        _,run=self.execute('water_gs',gs)
        self.assertEqual(run.returncode,0,run.stdout[-3000:]+run.stderr)
        differences=re.findall(r'DC #SCF.*diff =\s*(\S+)',run.stdout)
        self.assertLess(float(differences[-1]),1e-10)
        rt=gs.replace("theory='dft'","theory='tddft_response'\n yn_md='y'")
        rt=rt.replace("yn_dc='y'","yn_dc='n'\n yn_conventional_from_dcdft='y'")
        rt=rt.replace(' temperature_k=300d0','')
        # Keep length input in Angstrom; dt is fs in this sample's unit system.
        rt+="\n&tgrid\n dt=.0005\n nt=4\n/\n&emfield\n ae_shape1='impulse'\n e_impulse=.01\n epdir_re1=1,0,0\n/\n"
        rt+="\n&md\n ensemble='NVE'\n file_ini_velocity='velocity.dat'\n step_update_ps=1\n/\n&analysis\n out_rt_energy_step=1\n out_rvf_rt_step=1\n nenergy=10\n de=.01\n/\n"
        folder,run=self.execute('water_rt',rt,rt='water_gs')
        self.assertEqual(run.returncode,0,run.stdout[-3000:]+run.stderr)
        self.assertIn('end SALMON',run.stdout)
        self.assertNotIn('Native LCFO RT active',run.stdout)
        energy=np.loadtxt(next(folder.glob('*_rt_energy.data')))
        self.assertTrue(np.isfinite(energy).all())
        self.assertGreater(energy[-1,3],0.)
        print('water mesh final ionic kinetic energy eV:',energy[-1,3])

    def test_forbidden_modes(self):
        cases=[('pulse',self.rt_input().replace("ae_shape1='impulse'","ae_shape1='Acos2'"),'only impulse'),
               ('skip_ps',self.rt_input().replace('step_update_ps=1','step_update_ps=2'),'every step'),
               ('finite_support',self.rt_input().replace('exx_mlwf_interval=5','exx_mlwf_radius=2'),'static DFT only'),
               ('thermal',self.rt_input().replace(' nstate=2',' nstate=2\n temperature_k=300'),'fixed occupations')]
        for name,inp,expected in cases:
            with self.subTest(name=name):
                _,run=self.execute(name,inp,rt=True)
                self.assertNotEqual(run.returncode,0)
                self.assertIn(expected,run.stdout+run.stderr)

if __name__=='__main__':unittest.main()
