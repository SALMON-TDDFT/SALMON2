"""Static independent-particle RPA screening from a SALMON KLI GS (atomic units).

Integrate the fully microscopic W over pairs of normalized real-space box
densities. Only these integrals go into RT; no empty orbitals are written.
"""
import argparse
from itertools import product
from pathlib import Path
from time import perf_counter
import numpy as np
from q_mesh import signed_indices, interpolate_boxes
from box_factors import box_forms


def transform_boxes_inplace(values, mesh):
    """Normalized q-to-R FFT. Both indices are stored with x varying fastest."""
    nk = len(values)
    pairs = values.reshape(nk, -1)
    # Bound the temporary FFT workspace independently of the box count.
    width = max(1, (4 * 1024**2) // (nk * values.itemsize))
    for start in range(0, pairs.shape[1], width):
        slab = pairs[:, start:start + width]
        grid = slab.reshape(*mesh[::-1], slab.shape[1])
        slab[:] = np.fft.fftn(grid, axes=(0, 1, 2)).reshape(slab.shape) / nk


def prepare(gs, output, mesh, boxes, cutoff, empty=None, hamiltonians=None, workers=1, time_reversal=False, pump=None, eta=None, qmesh=None, comm=None, head_tolerance=.05, eta_tolerance=.02, q_kernel=None, q_threads=1):
    if pump is not None and (time_reversal or eta is None or not np.isfinite(eta) or eta<=0):
        raise ValueError("pump screening requires positive --eta (Ha) and the full k mesh")
    if time_reversal and hamiltonians is None: raise ValueError("time reversal requires Sternheimer")
    if not 0 < head_tolerance < 1 or not 0 < eta_tolerance < 1:
        raise ValueError('head/eta relative tolerances must lie between zero and one')
    if boxes<1 or not np.isfinite(cutoff) or cutoff<=0 or any(m<1 for m in mesh):
        raise ValueError('box count, mesh and dielectric cutoff must be positive')
    if q_threads<1 or (q_kernel is None and q_threads!=1):
        raise ValueError('positive q threads and an explicit native library required')
    if q_kernel is not None:
        if empty is not None or workers!=1:raise ValueError('native q backend requires Sternheimer and workers=1')
        from q_openmp import QKernel
        QKernel(q_kernel)  # Fail before reading large input or entering collectives.
    from functools import partial
    restart=gs/'data_for_restart'
    refs=[p for p in (restart/'kli-v1.bin',restart/'tbmbj-wref.bin') if p.exists()]
    if len(refs)!=1:raise ValueError('expected exactly one static KLI or TB-mBJ screening reference')
    with refs[0].open('rb') as f:
        magic=f.read(4)
        if magic not in (b'KLI1',b'TBM1'):
            raise ValueError('full-mesh KLI checkpoint required')
        nx,ny,nz,spin,nk,nb,nion = np.fromfile(f,'=i4',7)
        cell = np.fromfile(f,'=f8',9).reshape(3,3)  # lattice vectors as rows
        k = np.fromfile(f,'=f8',3*nk).reshape(nk,3)
        occ = np.fromfile(f,'=f8',nb*nk*spin).reshape(spin,nk,nb)
        ions = np.fromfile(f,'=f8',3*nion)
        vref = np.fromfile(f,'=f8',nx*ny*nz)
        eref = np.fromfile(f,'=f8',1)
    n=np.array([nx,ny,nz]); mesh=np.array(mesh); na=boxes**3
    qmesh = mesh if qmesh is None else np.asarray(qmesh, dtype=int)
    if qmesh.shape != (3,) or np.any(qmesh < 1) or np.any(mesh % qmesh) or qmesh.prod() < 2:
        raise ValueError("q mesh must divide the k mesh on each axis and include nonzero q")
    if comm is not None and (empty is not None or workers != 1):
        raise ValueError("MPI requires a Sternheimer backend and one worker per rank")
    root = comm is None or comm.rank == 0
    nv=int((occ[0,0]==2).sum())
    if sum(x is not None for x in (empty,hamiltonians,pump))!=1: raise ValueError('select exactly one screening solver')
    bands=nv+(empty or 0)
    if spin!=1 or nk!=mesh.prod() or np.any(n%boxes) or (empty is not None and not 0<empty<=nb-nv):
        raise ValueError('need closed shell, full uniform mesh and divisible auxiliary grid')
    if not (np.all(occ[:,:,:nv]==2) and np.all(occ[:,:,nv:]==0)):
        raise ValueError('integer occupations required')
    if nv<1: raise ValueError('empty valence manifold')
    vol=abs(np.linalg.det(cell)); dv=vol/n.prod()
    wavefile=gs/'data_for_restart/wfn.bin'
    if wavefile.stat().st_size!=int(nk*nb*n.prod()*16):
        raise ValueError('incompatible wavefunction layout or checkpoint length')
    psi=np.memmap(wavefile,dtype='=c16',mode='r',
                  shape=(nk,nb,nz,ny,nx))[:,:bands]
    rows=[l for l in (gs/'Si_eigen.data').read_text().splitlines() if l.strip()[:1].isdigit()]
    eig=np.loadtxt(rows)[:,1].reshape(nk,nb)[:,:bands]
    if not np.isfinite(eig).all(): raise ValueError('nonfinite eigenvalues')
    if empty is not None and eig[:,nv:].min()<=eig[:,:nv].max(): raise ValueError('insulating gap required')
    rec=2*np.pi*np.linalg.inv(cell).T
    frac=k@cell.T/(2*np.pi)
    lookup={tuple(v):i for i,v in enumerate(np.rint((frac-frac[0])*mesh).astype(int)%mesh)}
    if nk==1: raise ValueError('cubic head quadrature requires a nonzero q star')
    if len(lookup)!=nk: raise ValueError('unexpected k mesh')
    np.testing.assert_allclose((frac-frac[0])*mesh,np.rint((frac-frac[0])*mesh),atol=1e-8)
    order=np.rint((frac-frac[0])*mesh).astype(int)
    # One representative for each supercell translation (x varies fastest).
    R = signed_indices(mesh)
    shifts = signed_indices(qmesh) * (mesh // qmesh)
    gmax=np.ceil(np.sqrt(2*cutoff)*np.linalg.norm(cell,axis=1)/(2*np.pi)+1).astype(int)
    g=np.array(list(product(*(range(-h,h+1) for h in gmax))))
    grid=np.indices(tuple(n[::-1])).reshape(3,-1)[::-1].T/n
    label=np.ravel_multi_index(tuple(np.floor(grid*boxes+1e-10).astype(int).T),(boxes,)*3,order='F')
    heads=[]; head_q=[]; offdiag=[]
    shells=[g[np.sum(((g+qi/mesh)@rec)**2,axis=1)/2<=cutoff] for qi in shifts]
    response=None
    metallic=False
    report={}
    response_start = perf_counter()
    if pump is not None:
        from sternheimer import polarizability
        from pump_screen import snapshot,solve_pump_target
        from metal_head import population_head
        snapshot_eig=[]; snapshot_occ=[]
        validation_error=None
        try:
            for j in (range(nk) if root else ()):
                step,ac,levels,pop,states,h=snapshot(pump,n,nk,k,cell,j,read_h=False,read_states=False)
                if j==0: pump_step,pump_ac=step,ac
                if step!=pump_step or not np.array_equal(ac,pump_ac):raise ValueError('mixed snapshot times')
                if (not np.isfinite(pop).all() or pop.min() < -1e-10 or pop.max()>2+1e-7 or
                        abs(pop.sum()-2*nv)>1e-4):raise ValueError('increase projection bands: incomplete pump population')
                snapshot_eig.append(levels);snapshot_occ.append(pop)
            if root:
                pops=np.asarray(snapshot_occ)
                head_population=population_head(pops,np.asarray(snapshot_eig),nv)
                metallic=not head_population['insulating']
                if metallic and np.any(qmesh < 6):
                    raise ValueError('pump snapshot without a demonstrated insulating gap needs >=6^3 q mesh for finite-q head validation')
                print(f'Pump step {pump_step}; '
                      f'mean excited population {pops[:,nv:].sum(axis=1).mean():.6e} electrons/cell; '
                      f'mean omitted population {(2*nv-pops.sum(axis=1)).mean():.6e}; '
                      f'maximum omitted at k {np.max(2*nv-pops.sum(axis=1)):.6e}',flush=True)
        except Exception as exc:
            validation_error=str(exc)
        if comm is not None:validation_error=comm.bcast(validation_error,root=0)
        if validation_error is not None:raise ValueError(validation_error)
        if comm is not None:
            metallic,head_population,pump_step=comm.bcast(
                (metallic,head_population,pump_step) if root else None,root=0)
        jobs=[(pump,n,k,cell,frac,order,mesh,shifts,shells,j,(eta,eta/2)) for j in range(nk)]
        response=polarizability(jobs,shifts,mesh,shells,workers,solver=partial(solve_pump_target,q_kernel=q_kernel,q_threads=q_threads),comm=comm,channels=2,distributed=True)
    if hamiltonians is not None:
        from sternheimer import polarizability,solve_target
        jobs=[(hamiltonians,gs,n,nk,nb,nv,cell,k,vref,eref,frac,order,mesh,shifts,shells,eig,j,magic) for j in range(nk)]
        response=polarizability(jobs,shifts,mesh,shells,workers,time_reversal,solver=partial(solve_target,q_kernel=q_kernel,q_threads=q_threads),comm=comm,distributed=True)
    # Each q/-q pair belongs to one rank after the k reduction. No root-only
    # dielectric inversion or full box-matrix replication follows that reduction.
    wR = {}
    w_eta = {} if pump is not None else None
    assembly_start = perf_counter()
    for iq in (range(len(shifts)) if response is None else sorted(response)):
        qi=shifts[iq]
        q=qi/mesh;gsel=shells[iq]
        if not len(gsel): raise ValueError('empty dielectric cutoff sphere')
        # q and -q have identical spheres up to reciprocal translation.
        Q=gsel+q; q2=np.sum((Q@rec)**2,axis=1)
        zero=q2<1e-14
        bare=np.divide(4*np.pi,q2,out=np.zeros_like(q2),where=~zero)
        if response is not None:
            chi=response.pop(iq)
        else:
            chi=np.zeros((len(Q),len(Q)),complex)
            for ik in range(nk):
                jk=lookup[tuple((order[ik]+qi)%mesh)]
                shift=np.rint(frac[ik]+q-frac[jk]).astype(int)
                G=gsel+shift
                if np.any(np.abs(G)*2>=n): raise ValueError('dielectric cutoff exceeds FFT grid Nyquist limit')
                ind=tuple((G%n).T[::-1])
                for i1,i2 in [(slice(0,nv),slice(nv,bands)),(slice(nv,bands),slice(0,nv))]:
                    # Including both occupation directions handles complex local fields at finite q.
                    den=np.abs(eig[ik,i1,None]-eig[jk,None,i2].reshape(1,-1))
                    pair=psi[ik,i1][:,None].conj()*psi[jk,i2][None,:]
                    fft=np.fft.fftn(pair,axes=(-3,-2,-1))
                    m=fft[(slice(None),slice(None),*ind)].reshape(-1,len(Q))*dv
                    m=m/np.sqrt(den.ravel())[:,None]
                    chi+=m.T@m.conj()
        b=np.sqrt(bare)
        eps=np.eye(len(Q))+2/(nk*vol)*(b[:,None]*chi*b[None,:])
        if pump is not None:
            if not np.isfinite(eps).all() or np.linalg.eigvalsh(eps).min() <= 1e-10:
                raise ValueError('unstable or nonfinite static pump dielectric matrix')
        inv=np.linalg.inv(eps)
        W=b[:,None]*inv*b[None,:]
        if pump is not None:
            coarse_W, W = W
            eps, inv = eps[1], inv[1]
        offdiag.append(np.linalg.norm(eps-np.diag(np.diag(eps))))
        if not zero.any():
            head=np.argmin(q2)
            heads.append(inv[head,head].real);head_q.append(q2[head])
        forms=box_forms(n,boxes,Q)
        mat=forms.conj()@W@forms.T/vol
        wR[iq]=mat
        if pump is not None:
            w_eta[iq]=forms.conj()@coarse_W@forms.T/vol
    statistics = [(heads,head_q,offdiag)] if comm is None else comm.allgather((heads,head_q,offdiag))
    heads=[x for item in statistics for x in item[0]]
    head_q=[x for item in statistics for x in item[1]]
    offdiag=[x for item in statistics for x in item[2]]
    qmin=min(head_q); einv=np.mean(np.array(heads)[np.isclose(head_q,qmin)])
    origin=None
    if pump is not None:
        from metal_head import check_eta, metallic_cell, metallic_sample_indices
        indices=[i for i in sorted(wR) if i != 0 or not metallic]
        error=0.
        for iq in indices:
            error=max(error,check_eta(w_eta[iq][None],wR[iq][None],eta_tolerance))
        if comm is not None:
            from mpi4py import MPI
            error=comm.allreduce(error,op=MPI.MAX)
        report=dict(eta_input_Ha=float(eta),eta_used_Ha=float(eta/2),eta_relative_error=error,
                    **head_population,head_tolerance=head_tolerance,eta_tolerance=eta_tolerance,
                    qmesh=qmesh.tolist(),mesh=mesh.tolist())
        if metallic:
            samples=metallic_sample_indices(qmesh,rec)
            selected={i:(wR[i],w_eta[i]) for i in samples if i in wR}
            parts=[selected] if comm is None else comm.gather(selected,root=0)
            result=None
            validation_error=None
            try:
                if root:
                    gathered={i:pair for part in parts for i,pair in part.items()}
                    fine=np.array([gathered[i][0] for i in samples])
                    coarse=np.array([gathered[i][1] for i in samples])
                    point,origin,details=metallic_cell(fine,qmesh,rec,mesh,head_tolerance,indices=samples)
                    _,average_eta,_=metallic_cell(coarse,qmesh,rec,mesh,head_tolerance,indices=samples)
                    details['head_eta_relative_error']=check_eta(average_eta[None],origin[None],eta_tolerance)
                    result=(point,origin,details)
                    del gathered,fine,coarse
            except Exception as exc:
                validation_error=str(exc)
            if comm is not None:validation_error=comm.bcast(validation_error,root=0)
            if validation_error is not None:raise ValueError(validation_error)
            if comm is not None:result=comm.bcast(result,root=0)
            point,origin,details=result
            report.update(details)
            if 0 in wR:wR[0]=point
            del selected,parts,result
        else:
            report['method']='common filled gapped manifold insulating head'
        del w_eta
    from distributed_boxes import q_blocks_to_pairs,interpolate_transform,stream_write
    first,last,coarse=q_blocks_to_pairs(wR,qmesh,na,comm)
    del wR
    def long_range_pairs(q,lo,hi):
        if metallic:return np.zeros(hi-lo,complex)
        Q=g+q
        q2=np.sum((Q@rec)**2,axis=1)
        keep=(q2>1e-14)&(q2<=2*cutoff)
        Q=Q[keep]
        result=np.empty(hi-lo,complex)
        # Bound temporary pair-by-G products independently of the pair-owner size.
        width=max(1,(4*1024**2)//max(64*len(Q),1))
        for begin in range(lo,hi,width):
            end=min(begin+width,hi)
            pairs=np.arange(begin,end)
            left,right=pairs//na,pairs%na
            ids,inverse=np.unique(np.r_[left,right],return_inverse=True)
            forms=box_forms(n,boxes,Q,indices=ids)
            lf=forms[inverse[:len(pairs)]];rf=forms[inverse[len(pairs):]]
            result[begin-lo:end-lo]=np.sum(lf.conj()*rf*(4*np.pi*einv/(vol*q2[keep])),axis=1)
        return result
    fft_start=perf_counter()
    values=interpolate_transform(coarse,qmesh,mesh,long_range_pairs,first,last,
                                 origin=None if origin is None else origin.ravel()[first:last])
    del coarse
    if not metallic:
        qc=(6*np.pi**2/(vol*nk))**(1/3)
        values+=12*np.pi*einv/(qc*qc*vol*nk)
    fft_seconds=perf_counter()-fft_start
    header=tail=b''
    if root:
        header=b'HXP1' if pump is not None else b'HXT1' if magic==b'TBM1' else b'HXR1'
        header+=np.asarray([*n,*mesh,boxes,nv,nion],dtype='=i4').tobytes()
        for value in [cell,k,ions,vref,eref,R@cell]:
            header+=np.asarray(value,dtype='=f8').tobytes()
        if pump is not None:
            tail=np.asarray([pump_step,len(snapshot_eig[0])],dtype='=i4').tobytes()
            for value in (snapshot_eig,snapshot_occ):tail+=np.asarray(value,dtype='=f8').tobytes()
    stream_write(output,header,tail,values,first,last,na,comm,mesh=mesh)
    timing=np.array([assembly_start-response_start,fft_start-assembly_start,fft_seconds])
    if comm is not None:
        from mpi4py import MPI
        comm.Allreduce(MPI.IN_PLACE,timing,op=MPI.MAX)
    if root:
        print(f'Response preparation {timing[0]:.3f} s, dielectric/box assembly {timing[1]:.3f} s, '
              f'pair interpolation/FFT {timing[2]:.3f} s')
        print(f'Nk={nk}, empty={empty}, dielectric cutoff={cutoff} Ha, auxiliary boxes={boxes}³, '
              f'finite-q epsilon_M={1/einv:.6f}, max local-field off-diagonal norm={max(offdiag):.6f}')
        if pump is not None:
            import json
            output.with_suffix(output.suffix+'.screening.json').write_text(json.dumps(report,indent=2)+'\n')
            print('Pump static screening checks:',json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('gs',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--mesh',type=int,nargs=3,required=True)
    p.add_argument('--boxes',type=int,required=True);p.add_argument('--cutoff',type=float,required=True)
    backend=p.add_mutually_exclusive_group(required=True)
    backend.add_argument('--empty',type=int)
    backend.add_argument('--hamiltonians',type=Path,help='directory of SALMON KSH1 grid operators (Sternheimer)')
    backend.add_argument('--pump',type=Path,help='directory containing HSP1 instantaneous projection snapshots')
    p.add_argument('--eta',type=float,help='pump static resolvent broadening (Ha)')
    p.add_argument('--head-tolerance',type=float,default=.05,help='maximum relative small-q fit/nested-stencil error')
    p.add_argument('--eta-tolerance',type=float,default=.02,help='maximum relative W(eta) versus W(eta/2) error')
    p.add_argument('--time-reversal' ,action='store_true',help='verify and reuse conjugate k pairs in static screening')
    p.add_argument('--workers',type=int,default=1)
    p.add_argument('--q-mesh',type=int,nargs=3,help='screening grid, each dimension must divide --mesh')
    p.add_argument('--q-kernel',type=Path,help='audited native OpenMP response library')
    p.add_argument('--q-threads',type=int,default=1,help='outer q OpenMP workers (requires --q-kernel)')
    p.add_argument('--mpi',action='store_true',help='partition destination k across MPI ranks')
    a=p.parse_args()
    comm = None
    if a.mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
    try:
        prepare(a.gs,a.output,a.mesh,a.boxes,a.cutoff,a.empty,a.hamiltonians,a.workers,a.time_reversal,a.pump,a.eta,a.q_mesh,comm,a.head_tolerance,a.eta_tolerance,a.q_kernel,a.q_threads)
    except BaseException:
        if comm is not None and comm.size > 1:
            import traceback
            traceback.print_exc()
            comm.Abort(1)
        raise
