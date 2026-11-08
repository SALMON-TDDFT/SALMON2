"""Pair-owned box kernels: bounded redistribution, interpolation, FFT and I/O.

Pair indices are C-order flattened (left box, right box). Arrays retain all q
or R values but only the owning rank's pairs. Large arrays use anonymous scratch
mappings, never a gathered Nk*Na**2 array. MPI transfers use typed buffers.
"""
from itertools import product
import os
from pathlib import Path
import tempfile

import numpy as np
from q_mesh import signed_indices

DEFAULT_BYTES = 8 * 1024**2


def pair_bounds(na, size, rank):
    """Half-open balanced C-order pair interval; empty owners are valid."""
    return na*na*rank//size, na*na*(rank+1)//size


def _array(shape, blockbytes, dtype=np.complex128):
    count = int(np.prod(shape))
    if count*np.dtype(dtype).itemsize <= blockbytes or count == 0:
        return np.empty(shape,dtype=dtype)
    # The mapping retains the storage after the temporary file handle closes.
    # It is unlinked by TemporaryFile and disappears when the mapping is freed.
    with tempfile.TemporaryFile(prefix='hsex-pairs-') as file:
        return np.memmap(file,mode='w+',shape=shape,dtype=dtype)


def _errors(comm, error):
    errors = [error] if comm is None else comm.allgather(error)
    if any(errors):
        raise ValueError(next(value for value in errors if value))


def _root_error(comm, error):
    error = error if comm is None else comm.bcast(error,root=0)
    if error is not None:
        raise ValueError(error)


def q_blocks_to_pairs(local_values, qmesh, na, comm=None, blockbytes=DEFAULT_BYTES):
    """Redistribute unique q-owned full matrices to balanced pair owners.

    Peak transfer workspace is O(max(blockbytes, Na**2)); coarse output is
    O(Nq*Na**2/Nranks) and mapped to scratch above blockbytes. The input matrices
    remain caller-owned. All ranks must participate, including empty owners.
    """
    rank,size=(0,1) if comm is None else (comm.rank,comm.size)
    qmesh=np.asarray(qmesh,dtype=int)
    error=None
    if qmesh.shape!=(3,) or np.any(qmesh<1) or na<1 or blockbytes<1:
        error='invalid q mesh, box dimension or workspace limit'
    nq=int(np.prod(qmesh))
    ownership=np.zeros(nq,dtype=np.int32)
    owners=np.zeros(nq,dtype=np.int32)
    if error is None:
        for iq,value in local_values.items():
            if (not isinstance(iq,(int,np.integer)) or not 0<=iq<nq or
                    np.shape(value)!=(na,na) or not np.isfinite(value).all()):
                error='invalid q-owned matrix'
                break
            ownership[iq]=1;owners[iq]=rank
    _errors(comm,error)
    if comm is not None:
        from mpi4py import MPI
        comm.Allreduce(MPI.IN_PLACE,ownership,op=MPI.SUM)
        comm.Allreduce(MPI.IN_PLACE,owners,op=MPI.SUM)
    if np.any(ownership!=1):
        raise ValueError('every q must have exactly one matrix owner')
    bounds=[pair_bounds(na,size,r) for r in range(size)]
    first,last=bounds[rank];width=last-first
    coarse=_array((nq,width),blockbytes)
    # Batch complete q matrices so send/receive buffers never grow with Nq.
    batch=max(1,blockbytes//(na*na*16))
    for begin in range(0,nq,batch):
        end=min(begin+batch,nq)
        selected=[q for q in range(begin,end) if owners[q]==rank]
        per_source=[np.flatnonzero(owners[begin:end]==r)+begin for r in range(size)]
        sendcounts=np.array([len(selected)*(b-a) for a,b in bounds],dtype=np.int64)
        recvcounts=np.array([len(q)*width for q in per_source],dtype=np.int64)
        sdispls=np.r_[0,np.cumsum(sendcounts[:-1])]
        rdispls=np.r_[0,np.cumsum(recvcounts[:-1])]
        send=np.empty(int(sendcounts.sum()),dtype=np.complex128)
        for dest,(a,b) in enumerate(bounds):
            chunk=send[sdispls[dest]:sdispls[dest]+sendcounts[dest]].reshape(len(selected),b-a)
            for j,iq in enumerate(selected):
                chunk[j]=np.asarray(local_values[iq]).reshape(-1)[a:b]
        receive=np.empty(int(recvcounts.sum()),dtype=np.complex128)
        if comm is None:
            receive[:]=send
        else:
            comm.Alltoallv([send,sendcounts,sdispls,MPI.DOUBLE_COMPLEX],
                           [receive,recvcounts,rdispls,MPI.DOUBLE_COMPLEX])
        for source,indices in enumerate(per_source):
            chunk=receive[rdispls[source]:rdispls[source]+recvcounts[source]].reshape(len(indices),width)
            coarse[indices]=chunk
    return first,last,coarse


def interpolate_transform(coarse, qmesh, mesh, long_range_pairs, first, last,
                          blockbytes=DEFAULT_BYTES, origin=None):
    """Pair-slab periodic interpolation and normalized q-to-R FFT.

    long_range_pairs(q, first, last) returns exactly that C-order pair interval.
    Equal q/electronic meshes skip interpolation, matching prepare_w. Optional
    origin replaces the final q=0 cell value BEFORE the FFT (metallic cell
    quadrature); it can be a local pair vector or a full box matrix.
    """
    qmesh=np.asarray(qmesh,dtype=int);mesh=np.asarray(mesh,dtype=int)
    nq=int(np.prod(qmesh));nk=int(np.prod(mesh));width=last-first
    if (qmesh.shape!=(3,) or mesh.shape!=(3,) or np.any(qmesh<1) or np.any(mesh<qmesh) or
            np.any(mesh%qmesh) or width<0 or np.shape(coarse)!=(nq,width) or blockbytes<1):
        raise ValueError('invalid pair-owned interpolation layout')
    output=_array((nk,width),blockbytes)
    if width==0:
        return output
    replace=None
    if origin is not None:
        replace=np.asarray(origin).reshape(-1)
        if len(replace)!=width:
            replace=replace[first:last]
        if len(replace)!=width or not np.isfinite(replace).all():
            raise ValueError('invalid origin cell values')
    same=np.array_equal(qmesh,mesh)
    # Build the interpolated matrix q-first: expensive q-dependent Coulomb
    # forms are not recomputed for every tiny FFT pair slab at large Nk.
    rowwidth=max(1,blockbytes//64)
    if same:
        output[:]=coarse
    else:
        source=_array((nq,width),blockbytes)
        for iq,q in enumerate(signed_indices(qmesh)/qmesh):
            for start in range(0,width,rowwidth):
                stop=min(start+rowwidth,width)
                source[iq,start:stop]=coarse[iq,start:stop]-long_range_pairs(q,first+start,first+stop)
        for iq,site in enumerate(signed_indices(mesh)):
            point=site*qmesh/mesh
            lower=np.floor(point).astype(int);fraction=point-lower
            vertices=[]
            for delta in product(range(2),repeat=3):
                weight=np.prod(np.where(delta,fraction,1-fraction))
                if weight:
                    index=np.ravel_multi_index(tuple((lower+delta)%qmesh),qmesh,order='F')
                    vertices.append((index,weight))
            for start in range(0,width,rowwidth):
                stop=min(start+rowwidth,width)
                values=np.array(long_range_pairs(site/mesh,first+start,first+stop),dtype=complex,copy=True)
                for index,weight in vertices:
                    values+=weight*source[index,start:stop]
                output[iq,start:stop]=values
        del source
    if replace is not None:
        output[0]=replace
    # FFT only a bounded pair slab; allow the output itself to live on scratch.
    slabwidth=max(1,blockbytes//(16*4*nk))
    for start in range(0,width,slabwidth):
        stop=min(start+slabwidth,width)
        values=np.array(output[:,start:stop],dtype=complex,copy=True)
        shape=(*mesh[::-1],stop-start)
        output[:,start:stop]=np.fft.fftn(values.reshape(shape),axes=(0,1,2)).reshape(values.shape)/nk

    return output


def stream_write(output, header_bytes, tail_bytes, localpairs, first, last,
                 na, comm=None, *, mesh):
    """Validate and stream pair-owned R matrices in legacy Fortran box order.

    Root holds only two Na*Na blocks to check W(R)=W(-R).T, realness and
    finiteness. A temporary file is published atomically without overwriting
    an existing output only after every block passes. Failed temporary files
    are retained for diagnosis; the requested usable output is never created.
    """
    rank,size=(0,1) if comm is None else (comm.rank,comm.size)
    mesh=np.asarray(mesh,dtype=int);nk=int(np.prod(mesh))
    bounds=[pair_bounds(na,size,r) for r in range(size)]
    error=None
    if (mesh.shape!=(3,) or np.any(mesh<1) or (first,last)!=bounds[rank] or
            np.shape(localpairs)!=(nk,last-first)):
        error='invalid pair-owned output layout'
    _errors(comm,error)
    counts=np.array([b-a for a,b in bounds],dtype=np.int64)
    displs=np.array([a for a,b in bounds],dtype=np.int64)
    positions=signed_indices(mesh)
    minus=np.ravel_multi_index(tuple((-positions%mesh).T),mesh,order='F')
    block=np.empty(na*na,dtype=complex) if rank==0 else None
    partner=np.empty_like(block) if rank==0 else None
    def gather(index,target):
        send=np.ascontiguousarray(localpairs[index],dtype=complex)
        if comm is None:
            target[:]=send
        else:
            from mpi4py import MPI
            comm.Gatherv(send,[target,counts,displs,MPI.DOUBLE_COMPLEX] if rank==0 else None,root=0)
    handle=None;temporary=None
    if rank==0:
        try:
            output=Path(output)
            if output.exists():
                raise FileExistsError(output)
            handle=tempfile.NamedTemporaryFile(mode='wb',prefix=output.name+'.partial-',
                                               dir=output.parent,delete=False)
            temporary=Path(handle.name)
            handle.write(header_bytes)
        except Exception as exc:
            error=f'cannot create screened-interaction output: {exc}'
    try:
        _root_error(comm,error)
        for ir in range(nk):
            gather(ir,block);gather(int(minus[ir]),partner)
            error=None
            if rank==0:
                try:
                    if (not np.isfinite(block).all() or not np.isfinite(partner).all() or
                            np.max(abs(block.imag))>1e-7):
                        raise ValueError('non-real or nonfinite screened interaction')
                    if not np.allclose(block.reshape(na,na),partner.reshape(na,na).T,atol=1e-8,rtol=0):
                        raise ValueError('screened interaction violates R/-R transpose reciprocity')
                    handle.write(np.ascontiguousarray(block.reshape(na,na).real.T).tobytes())
                except Exception as exc:
                    error=f'screened-interaction block {ir}: {exc}'
            _root_error(comm,error)
        error=None
        if rank==0:
            try:
                handle.write(tail_bytes);handle.close();handle=None
                # link is atomic and fails if output appeared during processing.
                os.link(temporary,output)
                temporary.unlink()
            except Exception as exc:
                error=f'cannot publish screened interaction: {exc}'
        _root_error(comm,error)
    finally:
        if handle is not None:
            handle.close()
