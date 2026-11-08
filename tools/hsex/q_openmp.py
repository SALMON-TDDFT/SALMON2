"""Optional audited native q-loop backend; MPI remains on the calling thread."""
import ctypes
import hashlib
import json
from pathlib import Path
import numpy as np


class QKernel:
    def __init__(self,path):
        path=Path(path).resolve()
        record=json.loads((path.parent/'audit.json').read_text())
        if hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256']:
            raise ValueError('native q library SHA256 differs from its validated build record')
        self.library=ctypes.CDLL(str(path))
        try:
            abi = self.library.hsex_q_abi
        except AttributeError as exc:
            raise ValueError('native q kernel predates original-H refinement; rebuild it') from exc
        abi.argtypes=[];abi.restype=ctypes.c_int
        if abi()!=2:
            raise ValueError('unsupported native q kernel ABI; rebuild it')
        self.function=self.library.hsex_q_solve
        self.function.argtypes=[ctypes.c_int]*8+[ctypes.c_void_p]*14+[ctypes.c_int,ctypes.c_void_p]
        self.function.restype=None

    def solve(self,adj,d,e,u,psi,phase,counts,levels,occ,eta,threads,pump,original_h=None):
        arrays=[np.asfortranarray(x,dtype=np.complex128) for x in (adj,u,psi,phase)]
        adj,u,psi,phase=arrays
        d,e,levels,occ,eta=[np.asfortranarray(x,dtype=np.float64) for x in (d,e,levels,occ,eta)]
        counts=np.ascontiguousarray(counts,dtype=np.int32)
        n,nb,nq=psi.shape;mg=phase.shape[1];nc=len(eta);nu=u.shape[1]
        if (n<2 or nb<1 or nq<1 or mg<1 or nc<1 or threads<1 or adj.shape!=(n,n) or
                u.shape!=(n,nu) or phase.shape!=(n,mg,nq) or counts.shape!=(nq,) or
                d.shape!=(n,) or e.shape!=(n-1,) or levels.shape!=(nb,nq) or occ.shape!=(nb,nq) or
                np.any(counts<1) or np.any(counts>mg) or
                (pump and (nu!=0 or np.any(eta<=0))) or (not pump and np.any(eta!=0))):
            raise ValueError('invalid native q response inputs')
        if not all(np.isfinite(x).all() for x in (adj,u,psi,phase,d,e,levels,occ,eta)):
            raise ValueError('nonfinite native q response inputs')
        out=np.zeros((mg,mg,nc,nq),complex,order='F');res=np.zeros(nq)
        status=np.zeros(nq,np.int32);team=np.zeros(1,np.int32)
        refine = original_h is not None
        h = np.asfortranarray(original_h if refine else np.zeros((1,1)), dtype=complex)
        if refine and (not pump or h.shape != (n,n) or not np.isfinite(h).all()):
            raise ValueError('invalid original-H refinement matrix')
        buffers=(adj,d,e,u,psi,phase,counts,levels,occ,eta,out,res,status,team)
        self.function(n,nb,nq,mg,nu,nc,int(threads),int(pump),*[a.ctypes.data for a in buffers],int(refine),h.ctypes.data)
        if np.any(status):raise ValueError(f'native q response failure: {status.tolist()}')
        if int(team[0])!=threads:
            raise ValueError(f'native q OpenMP team {int(team[0])} differs from requested {threads}; check runtime limits')
        if not np.isfinite(out).all() or not np.isfinite(res).all() or res.max()>1e-8:
            raise ValueError('native q response nonfinite output or excessive residual')
        return [out[:g,:g,:,q].transpose(2,0,1).copy() for q,g in enumerate(counts)],float(res.max()),int(team[0])


def response_batches(kernel,adj,d,e,u,nq,load,threads,pump,etas,budget=64*1024**2,original_h=None):
    """Bound padded packed input by 64 MiB, except one indivisible q item.

    load(q) returns normalized orbitals, phases, energies and occupations.
    Target matrices are shared by all threads; native scratch is per worker.
    """
    result=[];worst=0.;team=1;pending=None;q=0
    while q<nq:
        batch=[];max_bands=0;max_g=0
        while q<nq and len(batch)<max(1,threads*2):
            item=pending if pending is not None else load(q)
            bands=max(max_bands,len(item[2]));ngvec=max(max_g,item[1].shape[1])
            size=(16*adj.shape[0]*(bands+ngvec)+16*bands)*(len(batch)+1)
            if batch and size>budget:
                pending=item;break
            pending=None;batch.append(item);max_bands=bands;max_g=ngvec;q+=1
        ng=adj.shape[0];nb=max(len(x[2]) for x in batch);mg=max(x[1].shape[1] for x in batch)
        count=len(batch)
        psi=np.zeros((ng,nb,count),complex,order='F');phase=np.zeros((ng,mg,count),complex,order='F')
        levels=np.zeros((nb,count),order='F');occ=np.zeros_like(levels);counts=[]
        for j,(p,ph,en,oc) in enumerate(batch):
            psi[:,:len(en),j]=p;phase[:,:ph.shape[1],j]=ph
            levels[:len(en),j]=en;occ[:len(en),j]=oc;counts.append(ph.shape[1])
        values,res,team=kernel.solve(adj,d,e,u,psi,phase,counts,levels,occ,etas,threads,pump,original_h=original_h)
        result.extend(values);worst=max(worst,res)
    return result,worst,team
