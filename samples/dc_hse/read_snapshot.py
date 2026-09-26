"""Read version-1 native Wannier diagnostics and convert Gamma factors to NPZ."""
import argparse
import json
from pathlib import Path
import numpy as np


def read_snapshot(path):
    with Path(path).open('rb') as f:
        marker=f.read(4)
        if marker==b'\x04\x03\x02\x01': endian='<'
        elif marker==b'\x01\x02\x03\x04': endian='>'
        else: raise ValueError('invalid snapshot endian marker')
        f.seek(0);header=np.fromfile(f,endian+'i4',14)
        if len(header)!=14 or header[1]!=1: raise ValueError('unsupported snapshot header')
        n=header[2:5];mesh=header[5:8];no=int(header[8]);nk=int(np.prod(mesh));ng=int(np.prod(n))
        if min(*n,*mesh,no)<1: raise ValueError('invalid dimensions')
        real=np.fromfile(f,endian+'f8',9)
        if len(real)!=9: raise ValueError('truncated snapshot metadata')
        def array(count,dtype,shape):
            data=np.fromfile(f,endian+dtype,count)
            if len(data)!=count: raise ValueError('truncated snapshot array')
            return data.reshape(shape,order='F')
        occ=array(no*nk,'f8',(no,nk))
        u=array(no*no*nk,'c16',(no,no,nk))
        phi=array(ng*no*nk,'c16',(ng,no,nk))
        q=array(ng*nk*no,'c16',(ng*nk,no))
        if f.read(1): raise ValueError('unexpected trailing snapshot data')
    for a in [real,occ,u,phi,q]:
        if not np.isfinite(a).all(): raise ValueError('nonfinite snapshot')
    return dict(n=n,mesh=mesh,occupation=occ,u=u,phi=phi,q=q,spacing=real[:3],omega=float(real[3]),
                spread=float(real[4]),gradient=float(real[5]),min_singular=float(real[6]),
                native_exchange_Ha=float(real[7]),scf_residual=float(real[8]),
                refreshes=int(header[9]),localization_iterations=int(header[10]),
                localization_status=int(header[11]),scf_iterations=int(header[12]),converged=bool(header[13]))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('snapshot',type=Path)
    p.add_argument('output',type=Path);p.add_argument('--allow-unconverged',action='store_true')
    a=p.parse_args();s=read_snapshot(a.snapshot)
    if not s['converged'] and not a.allow_unconverged: raise ValueError('SCF not converged; refusing reference export')
    if np.prod(s['mesh'])!=1: raise ValueError('pair diagnostic currently requires Gamma')
    shape=tuple(s['n']);no=s['q'].shape[1]
    q=s['q'].reshape(shape+(no,),order='F').transpose(3,0,1,2)
    metadata={k:v for k,v in s.items() if k not in ['occupation','u','phi','q','spacing','n','mesh']}
    metadata.update(snapshot=str(a.snapshot.resolve()),n=shape,mesh=s['mesh'].tolist(),
                    spacing_bohr=s['spacing'].tolist(),highest_source_occupation=float(np.max(s['occupation'][-1])),
                    occupation_sum=float(np.sum(s['occupation'])))
    np.savez(a.output,q=q,spacing=s['spacing'],omega=s['omega'],occupation=s['occupation'],
             metadata=json.dumps(metadata,default=int))
    print(json.dumps(metadata,indent=2,default=int))

if __name__=='__main__':main()
