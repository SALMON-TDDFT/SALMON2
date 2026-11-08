"""Static limit at finite broadening for a dephased pump snapshot.

Both directions of the Lindhard response are assembled by polarizability().
Working with the real part of the retarded resolvent avoids singular solves
at crossings. The caller must converge eta and the finite-q head estimate.
"""
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded
from scipy.linalg.lapack import zhetrd, zungqr


from box_factors import BoundedCache, integer_phases

_snapshot_cache=BoundedCache()


def clear_snapshot_cache(max_bytes=64*1024**2):
    global _snapshot_cache
    _snapshot_cache=BoundedCache(max_bytes=max_bytes)


def snapshot_cache_info():
    return _snapshot_cache.info()


def snapshot(directory, n, nk, k, cell, j, read_h=True, read_states=True):
    """Read validated snapshot data; metadata-only calls skip orbital payload.

    State and metadata cache entries are distinct: a metadata hit can never
    return None for a later state request. Neither cache retains descriptors.
    """
    path = (Path(directory) / f'hsex-pump-{j+1:06d}.bin').resolve()
    stat = path.stat()
    ng = int(n.prod())
    # File identity/change stamps and complete validation context prevent stale
    # reuse after snapshot replacement or a different geometry in the same path.
    key = (str(path),stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns,
           tuple(n),nk,j,bool(read_states),np.asarray(k[j]).tobytes(),np.asarray(cell).tobytes())
    cached = _snapshot_cache.get(key)
    if cached is None:
        with path.open('rb') as stream:
            tag = stream.read(4)
            if tag not in (b'HSP1', b'HSP2'):
                raise ValueError('invalid pump snapshot')
            header = np.fromfile(stream, '=i4', 7)
            np.testing.assert_array_equal(header[:5], [*n, nk, j+1])
            nb, step = header[5:]
            if not 0 < nb <= ng or step < 1:
                raise ValueError('invalid snapshot dimensions or step')
            np.testing.assert_array_equal(np.fromfile(stream, '=f8', 9), cell.ravel())
            np.testing.assert_array_equal(np.fromfile(stream, '=f8', 3), k[j])
            ac, energies, occ = (np.fromfile(stream, '=f8', size) for size in (3, nb, nb))
            if (len(ac)!=3 or len(energies)!=nb or len(occ)!=nb or
                    not all(np.isfinite(x).all() for x in (ac,energies,occ))):
                raise ValueError('incomplete or nonfinite snapshot metadata')
            offset = stream.tell()+16*nb*ng
            states = np.fromfile(stream, '=c16', nb*ng).reshape(nb,ng) if read_states else None
        count = ng*(ng+1)//2 if tag==b'HSP2' else ng*ng
        if stat.st_size != offset+16*count:
            raise ValueError('incomplete pump grid Hamiltonian')
        after=path.stat()
        if (after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns) != (
                stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns):
            raise ValueError('snapshot changed while being read')
        arrays=(ac,energies,occ,states) if read_states else (ac,energies,occ)
        for array in arrays:array.flags.writeable=False
        cached=(step,ac,energies,occ,states,offset,tag)
        payload=sum(x.nbytes for x in arrays)+len(str(path).encode())+len(key[-1])+len(key[-2])+256
        _snapshot_cache.put(key,cached,payload)
    step,ac,energies,occ,states,offset,tag=cached
    from operator_io import read_operator
    h = read_operator(path, offset, ng, tag == b'HSP2') if read_h else None
    return step, ac, energies, occ, states, h


def broadened_gram(b, d, e, energy, eta):
    ab = np.zeros((3, len(d)), complex)
    ab[1] = d - energy - 1j*eta
    ab[0, 1:] = ab[2, :-1] = e
    x = solve_banded((1, 1), ab, b)
    residual = ab[1, :, None]*x - b
    residual[1:] += e[:, None]*x[:-1]
    residual[:-1] += e[:, None]*x[1:]
    err = np.max(np.linalg.norm(residual, axis=0) / np.maximum(np.linalg.norm(b, axis=0), 1e-14))
    if not np.isfinite(err) or err > 1e-8:
        raise ValueError(f'pump Sternheimer residual {err}')
    gram = b.conj().T @ x
    return (gram + gram.conj().T)*0.5, err


def solve_pump_target(job, q_kernel=None, q_threads=1):
    directory, n, k, cell, frac, order, mesh, R, shells, j, eta = job
    etas = np.atleast_1d(eta)
    if not np.isfinite(etas).all() or np.any(etas <= 0):
        raise ValueError('positive numerical broadenings required')
    nk = len(k)
    _, _, levels, _, states, h = snapshot(directory, n, nk, k, cell, j)
    np.testing.assert_allclose(h, h.conj().T, atol=1e-9, rtol=0)
    u = states.T * np.sqrt(abs(np.linalg.det(cell)) / len(h))
    np.testing.assert_allclose(u.conj().T @ u, np.eye(len(levels)), atol=1e-8, rtol=0)
    defects = np.linalg.norm(h @ u - u*levels, axis=0)
    if not np.isfinite(defects).all() or defects.max() > 1e-4:
        raise ValueError(f'unconverged projection eigenstates at k={j+1}: {defects}')
    # Keep the pre-existing ZHETRD('L') operator, not an averaged/full noisy H.
    h = np.tril(h,-1)+np.tril(h,-1).conj().T+np.diag(h.diagonal().real)
    f, d, e, tau, status = zhetrd(h, lower=1, lwork=64*len(h))
    if status: raise ValueError(f'hetrd {status}')
    rot = np.eye(len(h), dtype=complex)
    rot[1:, 1:], _, status = zungqr(f[1:, :-1], tau, lwork=64*len(h))
    if status: raise ValueError(f'ungqr {status}')
    # The same adjoint rotates every right hand side for this target k.
    adjoint = np.asfortranarray(rot.conj().T)
    del rot, f
    dv = abs(np.linalg.det(cell))/len(h)
    lookup = {tuple(v % mesh): i for i, v in enumerate(order)}
    if q_kernel is not None:
        from q_openmp import QKernel,response_batches
        def load(q):
            qi,g=R[q],shells[q]
            src=lookup[tuple((order[j]-qi)%mesh)]
            _,_,energies,occ,psi,_=snapshot(directory,n,nk,k,cell,src,read_h=False)
            indices=g+np.rint(frac[src]+qi/mesh-frac[j]).astype(int)
            if np.any(abs(indices)*2>=n):raise ValueError('pump cutoff exceeds Nyquist')
            active=np.flatnonzero(occ!=0)
            if not len(active):active=np.array([0])
            return (psi[active].T*np.sqrt(dv),integer_phases(n,indices),energies[active],occ[active])
        values,worst,team=response_batches(QKernel(q_kernel),adjoint,d,e,np.empty((len(h),0),complex),
                                          len(R),load,q_threads,True,etas,original_h=h)
        if j==0:print(f'Native q OpenMP team: {team}; requested {q_threads}',flush=True)
        return [v[0] for v in values] if np.ndim(eta)==0 else values,worst
    from precise_response import original_gram
    result, worst = [], 0.
    for qi, g in zip(R, shells):
        src = lookup[tuple((order[j]-qi) % mesh)]
        _, _, energies, occ, psi, _ = snapshot(directory, n, nk, k, cell, src, read_h=False)
        indices = g + np.rint(frac[src]+qi/mesh-frac[j]).astype(int)
        if np.any(abs(indices)*2 >= n): raise ValueError('pump cutoff exceeds Nyquist')
        phase = integer_phases(n,indices)
        chi = np.zeros((len(etas), len(g), len(g)), complex)
        active = np.flatnonzero(occ != 0)
        # Bound the pack size while reusing the dense rotation across bands.
        bytes_per_band = len(h) * len(g) * adjoint.itemsize
        block_size = max(1, (8 * 1024**2) // bytes_per_band)
        for start in range(0, len(active), block_size):
            bands = active[start:start + block_size]
            rhs = psi[bands].T[:, :, None] * phase[:, None, :]
            rhs *= np.sqrt(dv)
            for local, band in enumerate(bands):
                for channel, broadening in enumerate(etas):
                    gram, err = original_gram(h, rhs[:, local, :], energies[band], broadening)
                    chi[channel] += occ[band] * 0.5 * gram
                    worst = max(err, worst)
        result.append(chi[0] if np.ndim(eta) == 0 else chi)
    return result, worst
