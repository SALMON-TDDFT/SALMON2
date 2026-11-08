"""Direct projected static Sternheimer reference on SALMON's finite real-space grid.

One Hermitian tridiagonal reduction per destination k, reused at all shifts.
The reduction is a change of grid coordinates, not a conduction eigensolve.
"""
from pathlib import Path
import numpy as np
from scipy.linalg import solveh_banded
from scipy.linalg.lapack import zhetrd, zungqr
from box_factors import integer_phases


def reduce_operator(h, occupied, energies):
    """Return a unitary grid rotation and the diagonals of QHQ + alpha P."""
    if not np.isfinite(h).all() or np.max(abs(h-h.conj().T))>1e-9:
        raise ValueError('non-Hermitian or nonfinite KLI Hamiltonian')
    u=occupied
    np.testing.assert_allclose(u.conj().T@u,np.eye(u.shape[1]),atol=1e-8,rtol=0)
    hu=h@u
    defect=np.linalg.norm(hu-u*energies,axis=0).max()
    if defect>1e-4: raise ValueError(f'Hamiltonian and occupied states mismatch: {defect} Ha')
    alpha=energies.max()+2.
    huq=hu-.5*u@(u.conj().T@hu+alpha*np.eye(len(energies)))
    h-=huq@u.conj().T+u@huq.conj().T
    f,d,e,tau,info=zhetrd(h,lower=1,overwrite_a=1,lwork=64*len(h))
    if info: raise ValueError(f'hetrd {info}')
    rotation=np.eye(len(h),dtype=complex)
    rotation[1:,1:],_,info=zungqr(f[1:,:-1],tau,lwork=64*len(h))
    if info: raise ValueError(f'ungqr {info}')
    return rotation,d,e


def resolvent_gram(b,d,e,energy):
    """B†(T-energy)^-1 B with a true (not recurrent) relative residual."""
    ab=np.zeros((2,len(d)));ab[0]=d-energy;ab[1,:-1]=e
    x=solveh_banded(ab,b,lower=True,check_finite=False)
    r=ab[0,:,None]*x-b
    r[1:]+=e[:,None]*x[:-1];r[:-1]+=e[:,None]*x[1:]
    residual=np.max(np.linalg.norm(r,axis=0)/np.maximum(np.linalg.norm(b,axis=0),1e-12))
    if not np.isfinite(residual) or residual>1e-8: raise ValueError(f'Sternheimer residual {residual}')
    return b.conj().T@x,residual


def equilibrium_operator(directory, j, n, k, cell, eref, vref, magic):
    from operator_io import read_operator
    path = Path(directory) / f'kli-hamiltonian-{j+1:06d}.bin'
    prefix = b'TSH' if magic == b'TBM1' else b'KSH'
    with path.open('rb') as stream:
        tag = stream.read(4)
        if tag not in (prefix+b'1', prefix+b'2'):
            raise ValueError('invalid grid Hamiltonian type')
        np.testing.assert_array_equal(np.fromfile(stream,'i4',4), [*n,j+1])
        for expected in (k[j], cell.ravel(), eref, vref):
            np.testing.assert_array_equal(np.fromfile(stream,'f8',expected.size), np.ravel(expected))
        offset = stream.tell()
    return read_operator(path, offset, int(n.prod()), tag.endswith(b'2'))


def solve_target(job, q_kernel=None, q_threads=1):
    """Occupation-lowering response for all q at a fixed destination k."""
    directory,gs,n,nk,nb,nv,cell,k,vref,eref,frac,order,mesh,R,shells,eig,j,magic=job
    ng=int(np.prod(n));dv=abs(np.linalg.det(cell))/ng
    psi=np.memmap(Path(gs)/'data_for_restart/wfn.bin',dtype='=c16',mode='r',shape=(nk,nb,ng))
    h = equilibrium_operator(directory,j,n,k,cell,eref,vref,magic)
    if not h.flags.writeable:
        h = h.copy(order='F')
    u=psi[j,:nv].T*np.sqrt(dv)
    rot,d,e=reduce_operator(h,u,eig[j,:nv])
    adjoint = np.asfortranarray(rot.conj().T)
    del h, rot
    lookup={tuple(x%mesh):i for i,x in enumerate(order)}
    if q_kernel is not None:
        from q_openmp import QKernel,response_batches
        def load(q):
            qi,g=R[q],shells[q]
            source=lookup[tuple((order[j]-qi)%mesh)]
            indices=g+np.rint(frac[source]+qi/mesh-frac[j]).astype(int)
            if np.any(abs(indices)*2>=n):raise ValueError('dielectric cutoff exceeds FFT grid Nyquist limit')
            return (psi[source,:nv].T*np.sqrt(dv),integer_phases(n,indices),eig[source,:nv],np.full(nv,2.))
        values,worst,team=response_batches(QKernel(q_kernel),adjoint,d,e,u,len(R),load,q_threads,False,[0.])
        if j==0:print(f'Native q OpenMP team: {team}; requested {q_threads}',flush=True)
        return [v[0] for v in values],worst
    result=[];worst=0.
    for qi,g in zip(R,shells):
        source=lookup[tuple((order[j]-qi)%mesh)]
        indices=g+np.rint(frac[source]+qi/mesh-frac[j]).astype(int)
        if np.any(abs(indices)*2>=n):raise ValueError('dielectric cutoff exceeds FFT grid Nyquist limit')
        phase=integer_phases(n,indices)
        b=(psi[source,:nv].T[:,:,None]*np.sqrt(dv)*phase[:,None,:]).reshape(ng,-1)
        b-=u@(u.conj().T@b)
        b=(adjoint@b).reshape(ng,nv,len(g))
        chi=np.zeros((len(g),len(g)),complex)
        for v in range(nv):
            term,res=resolvent_gram(b[:,v,:],d,e,eig[source,v])
            chi+=term;worst=max(worst,res)
        result.append(chi)
    return result,worst


def time_reversal_jobs(jobs,rank=0,size=1):
    """Verify zero-field H, occupied projectors and occupied energy operators."""
    directory,gs,n,nk,nb,nv,cell,k,vref,eref,_,_,_,_,_,eig,_,magic=jobs[0]
    ng=int(n.prod());dv=abs(np.linalg.det(cell))/ng
    psi=np.memmap(Path(gs)/'data_for_restart/wfn.bin',dtype='=c16',mode='r',shape=(nk,nb,ng))
    from scipy.spatial import cKDTree
    pairs = cKDTree(k).query(-k)[1]
    if np.max(abs(k+k[pairs]))>1e-11 or not np.array_equal(pairs[pairs],np.arange(nk)):
        raise ValueError('mesh must list k and -k in the same reciprocal convention; omit --time-reversal')
    selected=[]
    representatives = np.flatnonzero(np.arange(nk) <= pairs)
    for j in representatives[rank::size]:
        jm = pairs[j]
        a=psi[j,:nv].conj().T*np.sqrt(dv);b=psi[jm,:nv].T*np.sqrt(dv);over=a.conj().T@b
        if (np.linalg.norm(b-a@over)>2e-7 or
                np.linalg.norm((eig[j,:nv,None]-eig[jm,None,:nv].reshape(1,-1))*over)>1e-7):
            raise ValueError('occupied space breaks time reversal; omit --time-reversal')
        left,right=(equilibrium_operator(directory,index,n,k,cell,eref,vref,magic) for index in (j,jm))
        for start in range(0,ng,64):
            stop=min(ng,start+64)
            if not np.allclose(left[start:stop].conj(),right[start:stop],rtol=0,atol=1e-10):
                raise ValueError('Hamiltonian breaks time reversal; omit --time-reversal')
        selected.append((jobs[j],j!=jm))
    return selected


def polarizability(jobs,R,mesh,shells,workers,time_reversal=False,solver=solve_target,comm=None,channels=1,distributed=False):
    from concurrent.futures import ProcessPoolExecutor
    if workers<1:raise ValueError('positive worker count required')
    if not isinstance(channels, int) or channels < 1:
        raise ValueError('positive integer channel count required')
    shape = () if channels == 1 else (channels,)
    forward=[np.zeros((*shape,len(g),len(g)),complex) for g in shells];res=0.
    lookup={tuple(x%mesh):i for i,x in enumerate(R)}
    minus=[]
    for i,qi in enumerate(R):
        j=lookup[tuple(-qi%mesh)]
        gm={tuple(g):n for n,g in enumerate(shells[j])}
        idx=[gm[tuple(-g-(qi+R[j])//mesh)] for g in shells[i]]
        minus.append((j,np.ix_(idx,idx)))
    rank, size = (0, 1) if comm is None else (comm.rank, comm.size)
    selected = (time_reversal_jobs(jobs,rank,size) if time_reversal else
                [(job,False) for job in jobs[rank::size]])
    count = len(selected) if comm is None else comm.allreduce(len(selected))
    if rank == 0:
        print(f'Sternheimer target k solves: {count} / {len(jobs)} on {size} MPI ranks',flush=True)
    from contextlib import nullcontext
    with (ProcessPoolExecutor(max_workers=workers) if workers>1 else nullcontext()) as pool:
        mapper = map if pool is None else pool.map
        results=mapper(solver,(job for job,_ in selected))
        for (_,paired),(pieces,err) in zip(selected,results):
            for i,(total,piece) in enumerate(zip(forward,pieces)):
                total+=piece
                if paired:
                    j,idx=minus[i]
                    total+=pieces[j][...,idx[0],idx[1]].conj()
            res=max(res,err)
    if comm is not None:
        from mpi4py import MPI
        for i,matrix in enumerate(forward):
            # Keep reverse-q partners together so assembly needs no replication.
            owner = min(i,minus[i][0]) % size if distributed else 0
            comm.Reduce(MPI.IN_PLACE if rank == owner else matrix,
                        matrix if rank == owner else None, root=owner)
            if distributed and rank != owner:
                forward[i] = None
        res = comm.allreduce(res, op=MPI.MAX)
        if not distributed and rank != 0:
            return None
    if rank == 0:
        print(f'Sternheimer maximum relative solve residual: {res:.3e}',flush=True)
    if distributed:
        return {i:f+forward[j][...,idx[0],idx[1]].conj()
                for i,(f,(j,idx)) in enumerate(zip(forward,minus)) if f is not None}
    return [f+forward[j][...,idx[0],idx[1]].conj() for f,(j,idx) in zip(forward,minus)]
