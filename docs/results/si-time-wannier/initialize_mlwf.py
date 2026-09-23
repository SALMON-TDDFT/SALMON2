"""Create the initial bond-centered MLWF gauge used by analyze_dense.py."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import json
import numpy as np
from wannier import read_u,geometry,initial_gauge
from mlwf import overlap_mesh,minimize
root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
r,c,k,b,L,h=geometry(root)
u,_=read_u(root/'calculations/si_tdcdft_k4/gs/data_for_restart')
g,_=initial_gauge(u,k,r,b,L,h**3)
raw,nb,bv,wt=overlap_mesh(u,k,r,L,h**3)
g,log=minimize(g,raw,nb,bv,wt,tolerance=1e-6)
assert log['converged'],log
np.savez_compressed(out/'mlwf_initial.npz',gauge=g,raw_overlap=raw,neighbors=nb,b=bv,weights=wt)
(out/'mlwf_initial_status.json').write_text(json.dumps(log,indent=2)+'\n')
print(json.dumps(log,indent=2))
