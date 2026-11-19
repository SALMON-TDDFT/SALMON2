"""SALMON grid/time campaigns: planning, classification and adjacent comparisons."""
import unittest
import numpy as np
from salmonpy.errors import NumericalError,RuntimeFailure,InputError

class ConvergenceTests(unittest.TestCase):
 def test_grid_plan_and_provisional_campaign(self):
  from salmonpy.convergence import plan_grid_scan,run_grid_scan
  radii={'Si':{'radius_angstrom':.9,'definition':'explicit selected channel cutoff','pp_sha256':'a'*64}}
  plan=plan_grid_scan(np.diag([1.,1.,1.]),radii=radii,material_id='silicon',factor=.95)
  counts=[tuple(t['grid']['counts']) for t in plan['trials']]
  self.assertEqual(len(counts),len(set(counts)))
  self.assertEqual(plan['initial_spacing_angstrom'],.3)
  self.assertLessEqual(max(plan['trials'][-1]['grid']['actual_spacing_angstrom']),.1)
  radii['Si']['radius_angstrom']=4
  self.assertEqual(plan['radii']['Si']['radius_angstrom'],.9)
  visited=[]
  def evaluate(trial):
   visited.append(trial['grid']['counts']);return {'value':-2.,'scf_converged':False}
  with self.assertWarns(UserWarning):result=run_grid_scan(plan,evaluate,atol=1e-8,rtol=0.)
  self.assertEqual(len(visited),len(counts));self.assertTrue(result['provisional'])
  finite=[]
  stopped=run_grid_scan(plan,lambda trial:(finite.append(trial) or {'value':-2.,'scf_converged':True}),atol=1e-8,rtol=0.)
  self.assertEqual(len(finite),2)
  self.assertTrue(stopped['stopped_early'])
  self.assertGreater(stopped['selected_spacing_angstrom'],.1)
  whole=run_grid_scan(plan,lambda trial:{'value':-2.,'scf_converged':True},atol=1e-8,rtol=0.,stop_when_converged=False)
  self.assertEqual(len(whole['trials']),len(counts))
  self.assertFalse(whole['stopped_early'])
  self.assertIsNotNone(result['selected_spacing_angstrom'])
  self.assertFalse(result['kpoint_convergence_established'])
  self.assertEqual(result['plan']['radii']['Si']['pp_sha256'],'a'*64)
  with self.assertRaises(InputError):plan_grid_scan(np.eye(3),radii={'Si':{'radius_angstrom':1.}},material_id='si')

 def test_dt_numerical_runtime_failures_and_independent_impulse(self):
  from salmonpy.convergence import plan_dt_scan,run_dt_scan
  plan=plan_dt_scan(duration_fs=.01,levels=4,material_id='silicon',impulse=[.001,.002,0])
  made=[]
  class State:
   def __init__(self,dt):self.dt=dt;self.n=0;self.current=np.zeros(3);self.closed=False
   def step(self,A):
    self.n+=1
    if self.dt==.2:return np.array([np.nan,0,0])
    if self.dt==.1:raise RuntimeFailure('worker unavailable')
    self.current=np.array(A)*(self.n*self.dt)*(1+self.dt**2);return self.current
   def close(self):self.closed=True
  def factory(dt):s=State(dt);made.append(s);return s
  result=run_dt_scan(plan,factory,electron_number=lambda s:8+s.dt**2*s.n*s.dt,atol=.002,rtol=0.)
  self.assertEqual([t['status'] for t in result['trials']],['numerical_failure','runtime_failure','complete','complete'])
  self.assertEqual(made[0].n,1);self.assertTrue(all(s.closed for s in made))
  self.assertEqual(result['nonbreaking_dt_au'],.05)
  self.assertEqual(result['norm_acceptable_dt_au'],.05)
  self.assertEqual(result['trials'][0]['current'],None)
  self.assertEqual(result['impulse_A_over_c'],[.001,.002,0.])
  self.assertEqual(result['plan']['trials'][0]['dt_au'],.2)
  strict=run_dt_scan(plan,factory,electron_number=lambda s:8+s.dt**2*s.n*s.dt,atol=0.,rtol=0.)
  self.assertIsNone(strict['norm_acceptable_dt_au'])
  with self.assertRaises(InputError):run_dt_scan(plan,factory,electron_number=lambda s:8.,atol=-1,rtol=0)

 def test_dt_exact_sampling_close_failure_and_no_false_convergence(self):
  from salmonpy.convergence import plan_dt_scan,run_dt_scan
  plan=plan_dt_scan(duration_fs=.01,levels=3,material_id='m',impulse=[0,0,.01])
  class State:
   current=np.zeros(3)
   def __init__(self,dt):self.dt=dt;self.t=0
   def step(self,A):self.t+=self.dt;return np.asarray(A)*self.t
   def close(self):raise OSError('measurement/cleanup unavailable')
  report=run_dt_scan(plan,State,electron_number=lambda s:8.,atol=1e-15,rtol=0,current_atol=1e-15,current_rtol=0)
  self.assertEqual(report['norm_acceptable_dt_au'],.2)
  default=run_dt_scan(plan,State,electron_number=lambda s:8*(1-2e-5*s.t))
  self.assertEqual(default['rtol'],1e-5)
  self.assertIsNone(default['norm_acceptable_dt_au'])
  self.assertTrue(all(t['status']=='complete' and t['completed_steps']==t['steps'] for t in default['trials']))
  self.assertAlmostEqual(default['trials'][0]['norm_max_loss_fraction'],1.2e-5)
  self.assertEqual(default['trials'][0]['norm_max_gain_fraction'],0.)
  self.assertTrue(all(t['status']=='complete' and t['cleanup_error'] for t in report['trials']))
  self.assertEqual([t['steps'] for t in report['trials']],[3,6,12])
  self.assertEqual(report['duration_au'],.2*3)
  with self.assertRaises(InputError):plan_dt_scan(duration_fs=1,levels=True,material_id='m')
  for bad in (float('nan'),float('inf'),-float('inf'),0.):
   with self.subTest(norm=bad):
    broken=run_dt_scan(plan,State,electron_number=lambda s:bad if s.t>0 else 8.)
    self.assertTrue(all(t['status']=='numerical_failure' and t['failed_step']==1 and t['current'] is None for t in broken['trials']))

 def test_current_grid_uses_norm_qualified_dt_and_five_percent_rms(self):
  from salmonpy.convergence import plan_grid_scan,run_current_grid_scan
  plan=plan_grid_scan(np.diag([4.,4.,4.]),material_id='si',radii={'Si':{'radius_angstrom':1.2,'definition':'explicit cutoff','pp_sha256':'a'*64}})
  visits=[]
  def evaluate(trial):
   visits.append(trial['index']);i=trial['index']
   dt=.1 if i==0 else .05;steps=2 if i==0 else 4
   values=np.zeros((steps+1,3));values[:,0]=np.linspace(0,1,steps+1)*(1.1 if i==0 else 1.)
   return {'scf_converged':True,'dt_report':{'material_id':'si','duration_au':.2,'impulse_A_over_c':[1e-4,0,0],
      'atol':0.,'rtol':1e-5,'norm_acceptable_dt_au':dt,'trials':[{'dt_au':dt,'status':'complete','norm_within_tolerance':True,'current':values.tolist()}]}}
  result=run_current_grid_scan(plan,evaluate)
  self.assertEqual(visits,[0,1,2])
  self.assertAlmostEqual(result['trials'][1]['comparison']['relative_rms_error'],.1)
  self.assertFalse(result['trials'][1]['comparison']['within_tolerance'])
  self.assertTrue(result['trials'][2]['comparison']['within_tolerance'])
  self.assertTrue(result['stopped_early'])
  self.assertEqual(result['rtol'],.05)
  self.assertIsNotNone(result['selected_spacing_angstrom'])
  # A waveform at zero crossings does not impose pointwise relative thresholds.
  from salmonpy.convergence import _waveform_comparison
  result=_waveform_comparison([[0,0,0],[1,0,0]],[[.01,0,0],[1,0,0]],0.,.05)
  self.assertTrue(result['within_tolerance'])
  boundary=_waveform_comparison([[.95,0,0]],[[1,0,0]],0.,.05)
  self.assertTrue(boundary['within_tolerance'])
  self.assertFalse(_waveform_comparison([[1.0501,0,0]],[[1,0,0]],0.,.05)['within_tolerance'])
  zero=_waveform_comparison([[0,0,0]],[[0,0,0]],0.,.05)
  self.assertFalse(zero['within_tolerance']);self.assertIsNone(zero['relative_rms_error'])
