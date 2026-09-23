"""Endpoint sensitivity to the ambiguous even-grid Nyquist derivative."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import sys,json
import numpy as np
out=Path(__file__).resolve().parent;root=out.parents[2]
sys.path.insert(0,str(out.parent/'si-time-wannier'))
from wannier import geometry,read_u
from tdelf import gradients,fields
r,c,k,b,L,h=geometry(root)
bond=np.min(np.linalg.norm((r[:,None,:]-b[None,:,:]+L/2)%L-L/2,axis=2),axis=1)<1.
result={}
for case in ['initial','none','weak','strong']:
    p=root/'calculations/si_tdcdft_k4/gs/data_for_restart' if case=='initial' else Path('/private/tmp/salmon-si-time-wannier-dense')/case/'checkpoint_rt_001600'
    u,_=read_u(p);result[case]={}
    for flag in [False,True]:
        z=fields(u,gradients(u,k,L,12,zero_nyquist=flag));n=z['n'];elf=z['elf'];w=n/n.sum()
        result[case][str(flag)]=dict(mean=float(w@elf),bond=float(n[bond]@elf[bond]/n[bond].sum()),
           squared_deviation=float(w@(elf-.5)**2),max_current_correction=float(np.max(abs(elf-z['elf_without_current']))))
(out/'derivative_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
