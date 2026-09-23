"""Compare native accepted checkpoint/current against the Python PT-CN endpoint."""
import argparse,json,re
from pathlib import Path
import numpy as np
from test_native_exchange import packed

def verify(native,reference):
 native=Path(native)
 with np.load(reference) as z:
  u=packed(z['u']);row=json.loads(str(z['trajectory']))[-1]
 checkpoint=native/f"checkpoint_rt_{row['step']:06d}"
 assert (checkpoint/'hse_restart.bin').exists(),'Missing native HSE restart physics metadata'
 x=np.fromfile(checkpoint/'wfn.bin',np.complex128).reshape(u.shape,order='F')
 error=float(np.linalg.norm(x-u)/np.linalg.norm(u))
 data=np.atleast_2d(np.loadtxt(native/'Si_rt.data'))
 current=data[np.argmin(abs(data[:,0]-row['time_au'])),13:16]
 current_error=float(np.max(abs(current-row['current'])))
 result=dict(step=row['step'],wavefunction_relative_error=error,current_max_abs_error=current_error)
 assert error<2e-9,result
 assert current_error<2e-11,result
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('native');p.add_argument('reference');a=p.parse_args()
 print(json.dumps(verify(a.native,a.reference),indent=2))
