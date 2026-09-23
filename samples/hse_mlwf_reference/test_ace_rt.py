import unittest
import numpy as np
from ace_rt import midpoint_step

class MidpointTests(unittest.TestCase):
 def test_linear_cayley_and_norm(self):
  u=np.array([1.+.2j,.3-.1j]);h=np.array([[.5,.1j],[-.1j,.8]])
  def build(x):return lambda y:-.2*y,-.2*x
  out,stats=midpoint_step(u,.15,lambda x:h@x,build,tolerance=1e-12,inner_tolerance=1e-13)
  total=h-.2*np.eye(2);ref=np.linalg.solve(np.eye(2)+.075j*total,(np.eye(2)-.075j*total)@u)
  np.testing.assert_allclose(out,ref,atol=3e-12)
  self.assertLess(abs(np.vdot(out,out)-np.vdot(u,u)),1e-11)
  self.assertLess(stats['full_residual'],1e-12)
 def test_nonlinear_refresh_and_full_residual(self):
  u=np.array([1.+.2j,.3-.1j]);h=np.array([[1.,.5],[.5,2.]])
  def build(x):
   v=-2*abs(x)**2
   return lambda y:v*y,v*x
  out,s=midpoint_step(u,.2,lambda x:h@x,build,tolerance=1e-11,inner_tolerance=1e-12)
  x=(out+u)/2;r=x-u+.1j*(h@x-2*abs(x)**2*x)
  self.assertLess(np.linalg.norm(r)/np.linalg.norm(u),1e-11)
  self.assertGreater(s['builds'],2)
 def test_stale_seed_real_input_and_limits(self):
  u=np.ones(2)
  def build(x):return lambda y:-.2*y,-.2*x
  out,s=midpoint_step(u,.1,lambda x:.5*x,build,initial_exchange=lambda x:-.4*x)
  ref=(1-.015j)/(1+.015j)*u
  np.testing.assert_allclose(out,ref,atol=1e-9)
  self.assertEqual(s['builds'],2)
  for kw in [dict(max_inner=0),dict(max_outer=-1),dict(max_inner=1.5)]:
   with self.assertRaises(ValueError):midpoint_step(u,.1,lambda x:x,build,**kw)
 def test_fail_closed(self):
  u=np.ones(2,dtype=complex)
  with self.assertRaises(RuntimeError):midpoint_step(u,1.,lambda x:100*x,lambda x:(lambda y:-y,-x),max_inner=2)
  with self.assertRaises(ValueError):midpoint_step(u,-1.,lambda x:x,lambda x:(lambda y:-y,-x))

if __name__=='__main__':unittest.main()
