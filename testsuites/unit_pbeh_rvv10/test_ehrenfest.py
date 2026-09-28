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
    functional="pbeh40_rvv10"
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='realspace-ehrenfest-')
        cls.root=Path(cls.temp.name)
        cls.addClassCleanup(cls.temp.cleanup)
        cls.base=(ROOT/'testsuites/422_H_dcdft_hse/inputfile').read_text()
        cls.base=cls.base.replace("xc='hse06'",f"xc='{cls.functional}'\n exx_mlwf_interval=5\n exx_mlwf_maxiter=100")
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
        self.assertIn('EXX_WANNIER',run.stdout)
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
        folder2,run2=self.execute('water_rt_spatial',rt.replace('nproc_rgrid=1,1,1','nproc_rgrid=1,2,1'),
                                  ranks=2,rt='water_gs')
        self.assertEqual(run2.returncode,0,run2.stdout[-2000:]+run2.stderr)
        self.assertIn('EXX_SPATIAL',run2.stdout)
        energy2=np.loadtxt(next(folder2.glob('*_rt_energy.data')))
        np.testing.assert_allclose(energy2,energy,atol=1e-7,rtol=1e-7)
        print('water spatial energy parity eV:',np.max(abs(energy2-energy)))
        # Four occupied water orbitals split unevenly over three orbital groups.
        for ranks,layout in [(3,'1,1,1'),(6,'1,2,1')]:
            inp=rt.replace('nproc_rgrid=1,1,1','nproc_rgrid='+layout).replace('nproc_ob=1','nproc_ob=3')
            folder3,run3=self.execute('water_rt_orbital'+str(ranks),inp,ranks=ranks,rt='water_gs')
            self.assertEqual(run3.returncode,0,run3.stdout[-3000:]+run3.stderr)
            self.assertIn('EXX_ORBITALS',run3.stdout)
            energy3=np.loadtxt(next(folder3.glob('*_rt_energy.data')))
            np.testing.assert_allclose(energy3,energy,atol=1e-7,rtol=1e-7)
            xyz=next(folder3.glob('*_trj.xyz')).read_text().splitlines()[-3:]
            reference_xyz=next(folder.glob('*_trj.xyz')).read_text().splitlines()[-3:]
            vf=lambda lines:np.array([[float(x) for x in line.split('#v=')[1].replace('#f=','').split()] for line in lines])
            np.testing.assert_allclose(vf(xyz),vf(reference_xyz),atol=1e-7,rtol=1e-7)
            print('water orbital ranks/energy parity eV:',ranks,np.max(abs(energy3-energy)))


    def pulse_input(self,dt=.08,nt=120,amplitude=.03,moving=True):
        inp=self.rt_input(dt=dt,nt=nt,moving=moving)
        inp=inp.replace("theory='tddft_response'","theory='tddft_pulse'")
        inp=inp.replace("ae_shape1='impulse'",f"ae_shape1='Acos2'\n E_amplitude1={amplitude}\n omega1=1.9634954084936207\n tw1=6.4\n t1_start=0\n phi_CEP1=0")
        return inp

    def test_finite_pulse_work_balance(self):
        errors=[];end_currents=[]
        for dt,nt in ((.08,120),(.04,240),(.02,480)):
            folder,run=self.execute('pulse_'+str(nt),self.pulse_input(dt,nt),rt=True)
            self.assertEqual(run.returncode,0,run.stdout[-3000:]+run.stderr)
            self.assertIn('end SALMON',run.stdout)
            self.assertNotIn('Native LCFO RT active',run.stdout)
            data=np.loadtxt(next(folder.glob('*_rt.data')))
            energy=np.loadtxt(next(folder.glob('*_rt_energy.data')))
            self.assertTrue(np.isfinite(data).all() and np.isfinite(energy).all())
            # SALMON reports electron matter current: electric current has opposite sign.
            power=16*8*8*np.sum((data[:,16:19]-data[:,13:16])*data[:,10:13],axis=1)
            work=np.cumsum(.5*dt*(np.r_[0.,power[:-1]]+power))
            total=energy[:,1]+energy[:,3]
            error=np.max(abs(total[1:]-total[0]-work))
            errors.append(float(error));end_currents.append(data[-1,13:16])
            self.assertGreater(abs(work[-1]),1e-6,'finite field must transfer energy')
            self.assertLess(np.max(abs(data[-10:,1:13])),1e-13,'field-free tail required')
        ratio=np.linalg.norm(end_currents[0]-end_currents[1])/np.linalg.norm(end_currents[1]-end_currents[2])
        print('pulse energy-work errors Ha:',errors,'current refinement ratio:',ratio)
        self.assertGreater(errors[0]/errors[1],3.3)
        self.assertGreater(errors[1]/errors[2],3.3)
        self.assertLess(errors[-1],1e-5)
        self.assertGreater(ratio,3.3)
        self.assertLess(ratio,4.7)

    def test_pulse_endpoint_observables(self):
        dt=.08
        folder,run=self.execute('pulse_endpoints',self.pulse_input(dt,80),rt=True)
        self.assertEqual(run.returncode,0,run.stdout[-2000:]+run.stderr)
        data=np.loadtxt(next(folder.glob('*_rt.data')))
        # Independent A(t), sampled symmetrically around the reported endpoint.
        def potential(t):
            x=t-3.2
            return np.where(abs(x)<3.2,-.03/1.9634954084936207*np.cos(np.pi*x/6.4)**2*np.sin(1.9634954084936207*x),0.)
        expected=-(potential(data[:,0]+dt)-potential(data[:,0]-dt))/(2*dt)
        np.testing.assert_allclose(data[:,10],expected,atol=2e-12,rtol=2e-10)
        lines=next(folder.glob('*_trj.xyz')).read_text().splitlines()
        velocities=[]
        for i,line in enumerate(lines):
            if line.strip()=='4':
                velocities.append(np.array([[float(x) for x in row.split('#v=')[1].split('#f=')[0].split()] for row in lines[i+2:i+6]]))
        expected_current=np.array([v.sum(axis=0)/1024 for v in velocities[1:]])
        np.testing.assert_allclose(data[:,16:19],expected_current,atol=2e-13,rtol=2e-8)

    def test_zero_pulse(self):
        results=[]
        for name,inp in [('zero_pulse',self.pulse_input(nt=4,amplitude=0.)),
                         ('zero_impulse',self.rt_input(dt=.08,nt=4,impulse=0.))]:
            folder,run=self.execute(name,inp,rt=True)
            self.assertEqual(run.returncode,0,run.stdout[-2000:]+run.stderr)
            results.append(np.loadtxt(next(folder.glob('*_rt_energy.data'))))
        np.testing.assert_allclose(results[0],results[1],atol=1e-12,rtol=1e-10)

    def test_invalid_pulse(self):
        for name,old,new in [('frequency','omega1=1.9634954084936207','omega1=0'),
                             ('width','tw1=6.4','tw1=-1'),
                             ('start','t1_start=0','t1_start=-1')]:
            with self.subTest(name=name):
                _,run=self.execute('invalid_'+name,self.pulse_input().replace(old,new),rt=True)
                self.assertNotEqual(run.returncode,0)
                self.assertIn('positive frequency/width and nonnegative pulse start',run.stdout+run.stderr)

    def test_spatial_mesh_pulse(self):
        reference=None
        for ranks,layout in [(1,'1,1,1'),(2,'1,2,1'),(4,'1,2,2')]:
            inp=self.pulse_input(dt=.08,nt=120)
            inp=inp.replace('nproc_rgrid=1,1,1','nproc_rgrid='+layout)
            folder,run=self.execute('spatial_'+str(ranks),inp,ranks=ranks,rt=True)
            self.assertEqual(run.returncode,0,run.stdout[-3000:]+run.stderr)
            self.assertIn('end SALMON',run.stdout)
            self.assertNotIn('Native LCFO RT active',run.stdout)
            scratch=re.search(r'DC_LCFO_TILE scratch_points/global_points:\s*(\d+)\s+(\d+)',run.stdout)
            self.assertIsNotNone(scratch,'bounded reconstruction diagnostic missing')
            self.assertIn('DC_LCFO_STREAM native payload buffers',run.stdout)
            self.assertEqual(int(scratch[2]),1024)
            self.assertEqual(int(scratch[1]),1024//ranks)
            data=np.loadtxt(next(folder.glob('*_rt.data')))
            energy=np.loadtxt(next(folder.glob('*_rt_energy.data')))
            xyz=next(folder.glob('*_trj.xyz')).read_text().splitlines()[-4:]
            vf=np.array([[float(x) for x in line.split('#v=')[1].replace('#f=','').split()] for line in xyz])
            if reference is None:reference=(data,energy,vf)
            else:
                self.assertIn('EXX_SPATIAL',run.stdout)
                np.testing.assert_allclose(data,reference[0],atol=2e-10,rtol=2e-7)
                np.testing.assert_allclose(energy,reference[1],atol=2e-9,rtol=2e-7)
                np.testing.assert_allclose(vf,reference[2],atol=2e-9,rtol=2e-6)
                print('spatial pulse ranks/data/energy/vf differences:',ranks,
                      np.max(abs(data-reference[0])),np.max(abs(energy-reference[1])),np.max(abs(vf-reference[2])))

    def test_orbital_mesh_pulse(self):
        reference=None
        for ranks,layout,orbitals in [(1,'1,1,1',1),(2,'1,1,1',2),(4,'1,2,1',2),(8,'1,2,2',2)]:
            inp=self.pulse_input().replace('nproc_rgrid=1,1,1','nproc_rgrid='+layout)
            inp=inp.replace('nproc_ob=1','nproc_ob='+str(orbitals))
            folder,run=self.execute('orbital_pulse'+str(ranks),inp,ranks=ranks,rt=True)
            self.assertEqual(run.returncode,0,run.stdout[-3000:]+run.stderr)
            self.assertIn('end SALMON',run.stdout)
            self.assertNotIn('Native LCFO RT active',run.stdout)
            data=np.loadtxt(next(folder.glob('*_rt.data')))
            energy=np.loadtxt(next(folder.glob('*_rt_energy.data')))
            xyz=next(folder.glob('*_trj.xyz')).read_text().splitlines()[-4:]
            vf=np.array([[float(x) for x in line.split('#v=')[1].replace('#f=','').split()] for line in xyz])
            if reference is None:reference=(data,energy,vf)
            np.testing.assert_allclose(data,reference[0],atol=2e-10,rtol=2e-7)
            np.testing.assert_allclose(energy,reference[1],atol=2e-9,rtol=2e-7)
            np.testing.assert_allclose(vf,reference[2],atol=2e-9,rtol=2e-6)
            print('orbital Ehrenfest ranks/energy/vf error:',ranks,np.max(abs(energy-reference[1])),np.max(abs(vf-reference[2])))

    def test_spatial_work_refinement(self):
        errors=[]
        for dt,nt in ((.08,120),(.04,240),(.02,480)):
            inp=self.pulse_input(dt=dt,nt=nt).replace('nproc_rgrid=1,1,1','nproc_rgrid=1,2,2')
            folder,run=self.execute('spatial_work_'+str(nt),inp,ranks=4,rt=True)
            self.assertEqual(run.returncode,0,run.stdout[-2000:]+run.stderr)
            data=np.loadtxt(next(folder.glob('*_rt.data')))
            energy=np.loadtxt(next(folder.glob('*_rt_energy.data')))
            power=1024*np.sum((data[:,16:19]-data[:,13:16])*data[:,10:13],axis=1)
            work=np.cumsum(.5*dt*(np.r_[0.,power[:-1]]+power))
            total=energy[:,1]+energy[:,3]
            errors.append(float(np.max(abs(total[1:]-total[0]-work))))
        print('spatial MPI4 work errors Ha:',errors)
        self.assertGreater(errors[0]/errors[1],3.3)
        self.assertGreater(errors[1]/errors[2],3.3)
        self.assertLess(errors[-1],1e-5)

    def test_spatial_layout_guards(self):
        for name,layout in [('x','2,1,1')]:
            inp=self.rt_input().replace('nproc_rgrid=1,1,1','nproc_rgrid='+layout)
            _,run=self.execute('bad_layout_'+name,inp,rt=True)
            self.assertNotEqual(run.returncode,0)
            self.assertIn('Gamma y/z pencils required',run.stdout+run.stderr)

    def test_reconstruction_bad_coverage(self):
        folder=self.root/'bad_coverage_gs'
        shutil.copytree(self.root/'gs'/'data_dcdft',folder/'data_dcdft')
        path=folder/'data_dcdft/fragments/000001/rgrid_index.bin'
        mapping=np.frombuffer(path.read_bytes(),dtype=np.int32).copy()
        mapping[7]=mapping[6]  # duplicate core point; preserve first-coordinate header
        path.write_bytes(mapping.tobytes())
        inp=self.rt_input(nt=1).replace('nproc_rgrid=1,1,1','nproc_rgrid=1,2,1')
        _,run=self.execute('bad_coverage_rt',inp,ranks=2,rt='bad_coverage_gs')
        self.assertIn('invalid rgrid coverage',run.stdout+run.stderr)
        self.assertNotIn('start complex DC-LCFO wavefunction reconstruction',run.stdout)

    def test_reconstruction_nonfinite_payload(self):
        import struct
        folder=self.root/'bad_payload_gs'
        shutil.copytree(self.root/'gs'/'data_dcdft',folder/'data_dcdft')
        path=folder/'data_dcdft/fragments/000001/wavefunctions.bin'
        wire=bytearray(path.read_bytes())
        header=struct.unpack_from('=q',wire,40)[0]
        nb=struct.unpack_from('=i',wire,header+12+4)[0]
        # Header, k-record header, spin header, two fragment counts, row labels.
        payload=header+12+12+2*4+nb*4
        # Last saved orbital is unused by the occupied-only RT reader.
        struct.pack_into('=d',wire,payload+16*(nb*4-1),float('nan'))
        path.write_bytes(wire)
        inp=self.rt_input(nt=1).replace('nproc_rgrid=1,1,1','nproc_rgrid=1,2,1')
        _,run=self.execute('bad_payload_rt',inp,ranks=2,rt='bad_payload_gs')
        self.assertNotIn('start complex DC-LCFO wavefunction reconstruction',run.stdout)
        self.assertIn('invalid',run.stdout+run.stderr)

    def test_forbidden_modes(self):
        cases=[('pulse',self.rt_input().replace("ae_shape1='impulse'","ae_shape1='Ecos2'\n phi_CEP1=.25"),'impulse or Acos2'),
               ('skip_ps',self.rt_input().replace('step_update_ps=1','step_update_ps=2'),'every step'),
               ('finite_support',self.rt_input().replace('exx_mlwf_interval=5','exx_mlwf_radius=2'),'static DFT only'),
               ('thermal',self.rt_input().replace(' nstate=2',' nstate=2\n temperature_k=300'),'fixed occupations')]
        for name,inp,expected in cases:
            with self.subTest(name=name):
                _,run=self.execute(name,inp,rt=True)
                self.assertNotEqual(run.returncode,0)
                self.assertIn(expected,run.stdout+run.stderr)

if __name__=='__main__':unittest.main()
