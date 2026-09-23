import unittest
import numpy as np
from ptcn import pt_residual,ptcn_step

class PTTests(unittest.TestCase):
 def setUp(self):
  self.h=np.array([[.1,.08j,0],[-.08j,.7,.12],[0,.12,1.2]])
 def test_weighted_projection_and_gauge(self):
  q=np.linalg.qr(np.array([[1,2j],[.2j,1],[.4,.3j]]))[0];u=q.T[None]/np.sqrt(.4);hu=u@self.h.T
  f=pt_residual(u,hu,.4)
  np.testing.assert_allclose(u.conj()@f.transpose(0,2,1)*.4,0,atol=1e-14)
  v=np.array([[1,1j],[1j,1]])/np.sqrt(2)
  np.testing.assert_allclose(pt_residual(v@u,v@hu,.4),v@f,atol=1e-14)
 def test_stationary_occupied_subspace(self):
  _,v=np.linalg.eigh(self.h);u=v[:,:2].T[None]
  def build(x):return lambda y:-.2*y,-.2*x
  out,s=ptcn_step(u,.8,lambda x:x@self.h.T,build,1.)
  np.testing.assert_allclose(out,u,atol=1e-11)
  self.assertLess(s['full_residual'],1e-10)
 def test_second_order_density_and_exact_residual(self):
  u=np.array([[[1.,0,0]]],complex);e,v=np.linalg.eigh(self.h)
  exact=(v@np.diag(np.exp(-.4j*e))@v.conj().T)[:,0]
  pref=np.outer(exact,exact.conj());errs=[]
  def build(x):return lambda y:-.2*y,-.2*x
  for dt in [.2,.1,.05]:
   x=u.copy()
   for _ in range(round(.4/dt)):
    before=x.copy();x,s=ptcn_step(x,dt,lambda y:y@self.h.T,build,1.,tolerance=1e-12,inner_tolerance=1e-13)
    r=x-before+.5j*dt*(pt_residual(x,x@self.h.T-.2*x,1.)+pt_residual(before,before@self.h.T-.2*before,1.))
    self.assertLess(np.linalg.norm(r),2e-12)
   errs.append(np.linalg.norm(np.outer(x[0,0],x[0,0].conj())-pref))
  self.assertGreater(errs[0]/errs[1],3.8);self.assertGreater(errs[1]/errs[2],3.8)
 def test_nonlinear_exchange_and_cached_endpoint(self):
  u=np.array([[[1.,0,0]]],complex);cache={};local=lambda x:x@self.h.T
  def build(x):
   v=-.3*abs(x[0,0])**2
   op=lambda y:v*y
   cache['op']=op;cache['full']=op(x)
   return op,op(x)
  y,s=ptcn_step(u,.2,local,build,1.,initial_action=local(u)-.3*u,
    initial_exchange=lambda x:-.7*x,tolerance=1e-12,inner_tolerance=1e-13)
  f=lambda x:pt_residual(x,local(x)-.3*abs(x[0,0])**2*x,1.)
  self.assertLess(np.linalg.norm(y-u+.1j*(f(y)+f(u))),2e-12)
  self.assertGreater(s['builds'],1)
  cached_action=local(y)+cache['full'];cached_op=cache['op']
  a,sa=ptcn_step(y,.2,local,build,1.,initial_action=cached_action,initial_exchange=cached_op)
  b,sb=ptcn_step(y,.2,local,build,1.)
  np.testing.assert_allclose(a,b,atol=1e-11)
  self.assertEqual(sb['builds'],sa['builds']+1)
 def test_small_stagnating_inner_still_requires_full_gate(self):
  u=np.array([[[1.,0.]]]);calls=[0]
  def local(x):
   calls[0]+=1;h=.5*x.copy();h[0,0,1]+=(-1)**calls[0]*2e-10
   return h
  def build(x):return lambda y:-.2*y,-.2*x
  _,info=ptcn_step(u,.1,local,build,1.,initial_action=.3*u)
  self.assertGreater(info['inexact_inner_exits'],0)
  self.assertLess(info['full_residual'],1e-10)
  def bad_build(x):
   full=-.2*x.copy();full[0,0,1]+=1e-5
   return lambda y:-.2*y,full
  with self.assertRaises(RuntimeError):
   ptcn_step(u,.1,local,bad_build,1.,initial_action=.3*u,max_outer=2)
 def test_stale_seed_and_failure_guards(self):
  u=np.array([[[1.,0,0]]]);local=lambda x:x@self.h.T
  def build(x):return lambda y:-.2*y,-.2*x
  a,_=ptcn_step(u,.1,local,build,1.)
  b,_=ptcn_step(u,.1,local,build,1.,initial_action=local(u)-.2*u,initial_exchange=lambda x:-.7*x)
  np.testing.assert_allclose(a,b,atol=1e-9)
  with self.assertRaises(RuntimeError):ptcn_step(u,10.,local,build,1.,max_inner=1)
  for args in [dict(dt=-1),dict(dt=.1,max_outer=0),dict(dt=.1,dv=-1)]:
   kw=dict(dt=.1,dv=1.);kw.update(args)
   with self.assertRaises(ValueError):ptcn_step(u,local_action=local,build_exchange=build,**kw)

if __name__=='__main__':unittest.main()
