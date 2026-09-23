import unittest
import numpy as np
from probe_tools import central_response,apparent_width
class ProbeTests(unittest.TestCase):
 def test_background_even_term_and_cubic_convergence(self):
  t=np.arange(20.);background=np.sin(t);linear=np.cos(t);even=np.exp(-t/3);cubic=np.sin(t/2)
  errors=[]
  for h in (.1,.05):
   plus=background+h*linear+h*h*even+h**3*cubic
   minus=background-h*linear+h*h*even-h**3*cubic
   got=central_response(t,plus,t,minus,h)
   np.testing.assert_allclose(got,linear+h*h*cubic,atol=3e-15)
   errors.append(np.linalg.norm(got-linear))
  self.assertAlmostEqual(errors[0]/errors[1],4.,places=10)
 def test_mismatched_grids(self):
  with self.assertRaises(ValueError):central_response(np.arange(3.),np.ones(3),np.arange(3.)+.1,np.ones(3),.1)
  with self.assertRaises(ValueError):central_response(np.arange(3.),np.ones(3),np.arange(3.),np.ones(3),0.)
 def test_width_requires_both_crossings(self):
  e=np.linspace(2,4,2001);a=1/(1+((e-3)/.1)**2)
  self.assertAlmostEqual(apparent_width(e,a),.2,places=12)
  self.assertIsNone(apparent_width(e,np.ones_like(e)))
  self.assertIsNone(apparent_width(e,-a))
if __name__=='__main__':unittest.main()
