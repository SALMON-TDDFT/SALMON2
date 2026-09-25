"""Check short pulse ACE/full action agreement and time-step convergence."""
from pathlib import Path
import argparse,json
import numpy as np
p=argparse.ArgumentParser();p.add_argument('work');a=p.parse_args();b=Path(a.work)
data={}
for name in ['ace','full','half','quarter']:
 status=json.loads((b/name/'status.json').read_text());assert status['completed'],status
 data[name]=np.loadtxt(b/name/'Si_rt.data')
 assert np.isfinite(data[name]).all()
r=data['ace'];f=data['full'];h=data['half'][1::2];q=data['quarter'][3::4]
for t in [f,h,q]:np.testing.assert_allclose(t[:,0],r[:,0],atol=1e-9,rtol=0)
scale=np.linalg.norm(f[:,15]);relative=float(np.linalg.norm(r[:,15]-f[:,15])/scale)
d1=float(np.linalg.norm(r[:,15]-h[:,15]));d2=float(np.linalg.norm(h[:,15]-q[:,15]))
metrics=dict(ace_full_relative_J=relative,dt_coarse_half_J_L2=d1,dt_half_quarter_J_L2=d2,dt_error_ratio=d1/d2,passed=relative<1e-3 and d2<d1)
(b/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n');print(json.dumps(metrics,indent=2));assert metrics['passed']
