import unittest
import tempfile,json
from pathlib import Path
import numpy as np
try:
 from compare_linear import same_window,post_impulse
except ImportError:
 same_window=post_impulse=None


class LinearComparisonTests(unittest.TestCase):
 def setUp(self):self.assertIsNotNone(same_window,'matched-window comparison missing')
 def test_two_sampling_rates_and_amplitudes(self):
  t=np.arange(1,1001)*.08;u=np.arange(1,251)*.32;e=np.linspace(2,6,101)
  j=lambda x:np.exp(-.03*x)*np.cos(.15*x)
  a,b=same_window(t,1e-4*j(t),1e-4,u,2e-4*j(u),2e-4,80.,e)
  self.assertLess(np.linalg.norm(a.imag-b.imag)/np.linalg.norm(a.imag),.02)
  c,d=same_window(t,2e-4*j(t),2e-4,u,4e-4*j(u),4e-4,80.,e)
  np.testing.assert_allclose(a,c);np.testing.assert_allclose(b,d)
 def test_delayed_impulse_origin(self):
  time=np.arange(1,101)*.08;a=np.where(time>2.4+1e-10,1e-4,0.)
  current=np.where(time>2.4+1e-10,np.cos(time-2.4)*1e-4,0.)
  t,j,kick,origin=post_impulse(time,current,a)
  self.assertAlmostEqual(origin,2.4);self.assertAlmostEqual(t[0],.08)
  self.assertEqual(kick,1e-4);np.testing.assert_allclose(j,np.cos(t)*1e-4)
 def test_reject_short_mismatched_end_and_wrong_field(self):
  t=np.arange(1,101)*.08;j=np.cos(t);e=np.array([2.,3.])
  for end in (9.,7.99):
   with self.assertRaises(ValueError):same_window(t,j,1.,t,j,1.,end,e)
  with self.assertRaises(ValueError):post_impulse(t,j,np.linspace(0,1,len(t)))
  with self.assertRaises(ValueError):same_window(t,j,0.,t,j,1.,8.,e)

 def test_complete_analysis_and_incomplete_gate(self):
  from compare_linear import analyze
  with tempfile.TemporaryDirectory() as directory:
   root=Path(directory);h=root/'hse';d=root/'td';out=root/'report';h.mkdir();d.mkdir()
   t=np.arange(1,251)*.32;j=1e-4*np.exp(-.03*t)*np.cos(.15*t)
   rows=[dict(step=i+1,time_au=float(x),current=[0.,0.,float(y)]) for i,(x,y) in enumerate(zip(t,j))]
   (h/'trajectory.json').write_text(json.dumps(rows))
   status=dict(status='running',accepted_step=250,amplitude=1e-4,dt_au=.32)
   (h/'status.json').write_text(json.dumps(status));(d/'status.json').write_text(json.dumps(dict(completed=True,exit_code=0)))
   (d/'inputfile').write_text('synthetic oscillator')
   s=np.arange(1,1001)*.08;data=np.zeros((len(s),16));data[:,0]=s;data[:,3]=1e-4
   data[:,15]=1e-4*np.exp(-.03*s)*np.cos(.15*s)
   np.savetxt(d/'Si_rt.data',data,header='Time[a.u.] Jm_z[a.u.]\n'+'header\n'*5)
   with self.assertRaisesRegex(ValueError,'complete'):analyze(h,d,out,end_time=80.)
   status['status']='completed';(h/'status.json').write_text(json.dumps(status))
   result=analyze(h,d,out,end_time=80.)
   self.assertEqual(result['status'],'completed');self.assertEqual(len(result['windows']),3)
   self.assertTrue((out/'spectra.png').exists());self.assertTrue((out/'spectra.csv').exists())

if __name__=='__main__':unittest.main()
