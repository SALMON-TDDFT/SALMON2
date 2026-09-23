"""Driver failure orchestration; scientific kernels are tested separately."""
import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import propagate_ptcn as driver
from checkpoint import load_checkpoint

class DriverFailureTests(unittest.TestCase):
 def test_unknown_exchange_method_rejected(self):
  with self.assertRaisesRegex(ValueError,'exchange method'):
   driver.run('missing','missing','missing',1,exchange_method='wrong')
  with self.assertRaisesRegex(ValueError,'blocked'):
   driver.run('missing','missing','missing',1,exchange_backend=object())
 def test_endpoint_failure_keeps_coherent_accepted_checkpoint(self):
  with tempfile.TemporaryDirectory() as folder:
   p=Path(folder);export=p/'export';export.mkdir();(export/'metadata.txt').write_text('fixture');(export/'complete.txt').write_text('fixture')
   gs=p/'gs';gs.mkdir();u=np.zeros((64,16,12,12,12),complex);flat=u.reshape(64,16,-1)
   for n in range(16):flat[:,n,n]=1
   g=np.broadcast_to(np.eye(16),(64,16,16)).copy();np.savez_compressed(gs/'state.npz',u=u,gauge=g)
   (gs/'result.json').write_text(json.dumps(dict(converged=True,final_pair_tolerance=0.)))
   class Model:
    shape=(12,12,12);nk=64;k=np.zeros((64,3))
    dv=1.;nlcc=np.zeros((12,)*3);native_energies=np.zeros(7);tsymbol=np.zeros((64,12,12,12));h=1.;nab=np.zeros((4,3))
    def __init__(self,*a):self.calls=0
    def set_field(self,a):pass
    def density(self,x):
     rho=np.zeros((12,)*3);rho.flat[:16]=2;return rho
    def core(self,x):
     self.calls+=1
     if self.calls==4:raise RuntimeError('injected endpoint-action failure')
     return np.zeros_like(x)
    def expectation(self,x,y):return float(2*np.vdot(x,y).real/64)
    def current(self,x):return np.zeros(3)
   class Localizer:
    def __init__(self,*a):self.previous=None
   class XC:
    def __init__(self,*a):pass
    def __enter__(self):return self
    def __exit__(self,*a):pass
    def close(self):pass
   class Functional:
    def __init__(self,*a,**kw):pass
    def exchange(self,x,*a):return -x,{}
    def close(self):pass
   def step(x,*a,**kw):return x.copy(),dict(full_residual=0.)
   with patch.object(driver,'NativeModel',Model),patch.object(driver,'Localizer',Localizer),patch.object(driver,'Semilocal',XC),patch.object(driver,'HSEFunctional',Functional),patch.object(driver,'ptcn_step',step),patch.object(driver,'hartree',lambda rho,h:(rho*0,0.)),patch.object(driver,'semilocal_potential',lambda rho,*a:(rho*0,0.)):
    with self.assertRaisesRegex(RuntimeError,'injected endpoint-action failure'):
     driver.run(export,gs/'state.npz',p/'out',1)
    with patch.object(driver,'Localizer',side_effect=AssertionError('blocked exchange must not localize')):
     with patch.object(driver,'DistanceExchange',create=True,return_value=object()):
      with self.assertRaisesRegex(RuntimeError,'injected endpoint-action failure'):
       driver.run(export,gs/'state.npz',p/'blocked',1,exchange_method='blocked')
    with patch.object(driver,'DistanceExchange',side_effect=AssertionError('supplied backend must be used')):
     with self.assertRaisesRegex(RuntimeError,'injected endpoint-action failure'):
      driver.run(export,gs/'state.npz',p/'supplied',1,exchange_method='blocked',exchange_backend=object())
    original_save=driver.save_checkpoint;writes=[0]
    def fail_second_save(*args,**kw):
     writes[0]+=1
     if writes[0]==2:raise OSError('injected checkpoint write failure')
     return original_save(*args,**kw)
    with patch.object(driver,'save_checkpoint',fail_second_save):
     with self.assertRaisesRegex(RuntimeError,'injected endpoint-action failure'):
      driver.run(export,gs/'state.npz',p/'write_failure',1)
    failed=json.loads((p/'write_failure/status.json').read_text())
    self.assertEqual(failed['status'],'failed');self.assertIn('checkpoint write also failed',failed['error'])
    with patch.object(driver,'HSEFunctional',side_effect=RuntimeError('injected constructor failure')):
     with self.assertRaisesRegex(RuntimeError,'injected constructor failure'):
      driver.run(export,gs/'state.npz',p/'constructor_failure',1)
    self.assertEqual(json.loads((p/'constructor_failure/status.json').read_text())['status'],'failed')
   _,_,metadata,rows,_=load_checkpoint(p/'out/restart.npz',{})
   self.assertEqual(metadata['step'],0);self.assertEqual(rows,[])
   status=json.loads((p/'out/status.json').read_text());self.assertEqual(status['status'],'failed')
   self.assertIn('injected endpoint-action failure',status['error'])

if __name__=='__main__':unittest.main()
