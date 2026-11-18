"""Collective retained-GS array export. No native mutation or implicit MPI import."""
import hashlib
import math
from pathlib import Path
import warnings
import numpy as np
from .errors import InputError,NumericalError
from .native_gs import _agree


def _root_error(comm,error):
    # Standard Python/SALMONpy input and I/O exceptions are picklable. Rank0 is
    # the only filesystem writer; communicate its outcome before the next phase.
    if comm.size>1:error=comm.bcast(error,root=0)
    if error is not None:raise error


def _tile(model,density,orbitals,occupations,eigenvalues):
    nk=math.prod(model.kmesh)
    expected=(model.states,nk,1)
    for array in (occupations,eigenvalues):
        if array.values.shape!=expected or array.lower_indices!=(1,1,1):
            raise InputError('unsupported retained orbital/k/spin layout')
    shape=density.values.shape
    if len(shape)!=3 or len(density.lower_indices)!=3:
        raise InputError('unsupported retained periodic density grid layout')
    if any(n<1 or lo<1 or lo+n-1>total for lo,n,total in zip(density.lower_indices,shape,model.grid)):
        raise InputError('retained density tile outside the physical grid')
    values=orbitals.values
    lower=orbitals.lower_indices
    if (values.ndim!=7 or len(lower)!=7 or values.shape[3]!=1 or values.shape[6]!=1
            or lower[3]!=1 or lower[6]!=1):
        raise InputError('unsupported retained spin/orbital/k/image layout')
    starts=tuple(density.lower_indices[i]-lower[i] for i in range(3))
    if any(start<0 or start+n>values.shape[i] for i,(start,n) in enumerate(zip(starts,shape))):
        raise InputError('retained orbital halos do not cover the physical grid')
    region=tuple(slice(start,start+n) for start,n in zip(starts,shape))
    waves=values[region+(0,slice(None),slice(None),0)].transpose(4,3,2,1,0)
    origin=(lower[5]-1,lower[4]-1,*(v-1 for v in reversed(density.lower_indices)))
    occ=occupations.values[:,:,0].T
    eigen=eigenvalues.values[:,:,0].T
    signature=(origin,waves.shape,hashlib.sha256(occ.tobytes()).hexdigest(),
               hashlib.sha256(eigen.tobytes()).hexdigest())
    return waves,occ,eigen,signature


def _validate_tiles(tiles,full_shape):
    reference=tiles[0][2:]
    for origin,shape,*hashes in tiles:
        if tuple(hashes)!=reference:raise InputError('replicated GS occupations/eigenvalues differ across ranks')
        if len(origin)!=5 or len(shape)!=5 or any(n<1 or lo<0 or lo+n>total for lo,n,total in zip(origin,shape,full_shape)):
            raise InputError('retained orbital tile outside global k/band/grid bounds')
    for i,(origin,shape,*_) in enumerate(tiles):
        for other,extent,*_ in tiles[:i]:
            if all(max(a,b)<min(a+n,b+m) for a,b,n,m in zip(origin,other,shape,extent)):
                raise InputError('retained orbital tiles overlap')
    if sum(math.prod(t[1]) for t in tiles)!=math.prod(full_shape):
        raise InputError('retained orbital tiles do not cover the complete GS')


def _assemble(comm,local,tiles,full_shape):
    if comm.size==1:return local
    waves=None;error=None
    if comm.rank==0:
        try:waves=np.empty(full_shape,dtype=np.complex128)
        except Exception as exc:error=exc
    _root_error(comm,error)
    # One sender per collective limits root scratch to one tile, rather than
    # retaining all rank payloads alongside the assembled wavefunction.
    for sender,(origin,shape,*_) in enumerate(tiles):
        payload=comm.gather(local if comm.rank==sender else None,root=0)
        error=None
        if comm.rank==0:
            try:
                block=payload[sender]
                if block.shape!=shape:raise InputError('received GS tile shape changed')
                region=tuple(slice(lo,lo+n) for lo,n in zip(origin,shape))
                waves[region]=block
            except Exception as exc:error=exc
            finally:block=None
        del payload
        _root_error(comm,error)
    return waves


def _write(gs,destination,block,waves,occ,eigen,info):
    from .gs_archive import save_ground_state
    model=gs.calculation
    nk=math.prod(model.kmesh)
    rho=np.zeros(tuple(reversed(model.grid)))
    try:
        with np.errstate(over='raise',invalid='raise'):
            for k in range(nk):
                for band in range(model.states):
                    rho+=(occ[k,band]/nk)*(waves[k,band].real**2+waves[k,band].imag**2)
    except FloatingPointError as exc:
        raise NumericalError('nonfinite reconstructed retained GS density') from exc
    if not info.converged:
        warnings.warn('SCF not converged; exporting the retained last iterate without '
                      'certifying physical convergence',RuntimeWarning,stacklevel=3)
    provenance=dict(schema='salmonpy.retained-gs-export/1',material_id=gs.material_id,
        state_identity=gs.identity,target_sha256=gs.session.backend.sha256,
        converged=info.converged,iterations=info.iterations,
        density_source='reconstructed_from_native_orbitals',
        local_channels=[dict(atomic_number=p.atomic_number,local_channel=p.local_channel)
                        for p in model.pseudopotentials])
    references=[dict(atomic_number=p.atomic_number,filename=model._pseudo_name(i),
                     sha256=p.sha256,source='retained DFT model')
                for i,p in enumerate(model.pseudopotentials)]
    # Same ordinary full-mesh order as native lattice.f90:init_kvector.
    mx,my,mz=model.kmesh
    points=[[(ix+.5)/mx-.5,(iy+.5)/my-.5,(iz+.5)/mz-.5]
            for iz in range(mz) for iy in range(my) for ix in range(mx)]
    return save_ground_state(destination,wavefunctions=waves,density=rho,
        eigenvalues=eigen,occupations=occ,kpoints=points,weights=np.full(nk,1./nk),
        cell_bohr=model.cell,fractional_positions=model.fractional_positions,
        atomic_numbers=list(model.numbers),electrons=model.electrons,energy_hartree=info.energy_hartree,
        xc=model.xc,pseudopotentials=references,orbital_block=block,provenance=provenance)


def export(gs,directory,orbital_block):
    from .dft import DFT
    from .gs_archive import _integer
    gs._live()
    comm=gs.session.comm
    root=getattr(comm,'rank',0)==0
    error=None
    try:
        model=gs.calculation
        if not isinstance(model,DFT):
            raise InputError('export requires the immutable DFT model of a retained computed GS')
        if math.prod(gs.session.backend.decomposition)!=comm.size:
            raise InputError('GS decomposition does not match subgroup size')
        block=_integer(orbital_block,'orbital_block')
        destination=Path(directory).expanduser().absolute()
        descriptor=(gs._token,gs.identity,gs.material_id,str(destination),block,
                    gs.session.backend.sha256,gs.session.backend.decomposition,
                    model.input_text(),tuple(p.sha256 for p in model.pseudopotentials))
    except Exception as exc:error=str(exc)
    _agree(comm,error)
    _agree(comm,'export selections differ across ranks' if any(v!=descriptor for v in comm.allgather(descriptor)) else None)
    error=None
    if root:
        try:
            if destination.exists() or destination.is_symlink():raise FileExistsError(destination)
        except Exception as exc:error=exc
    _root_error(comm,error)
    info=gs.info
    metadata=comm.allgather(vars(info))
    _agree(comm,'retained GS metadata differs across ranks/model' if (
        info.atoms!=len(model.numbers) or info.bands!=model.states or
        any(v!=vars(info) for v in metadata)) else None)
    arrays=[gs.array(name) for name in ('density','orbitals','occupations','eigenvalues')]
    error=None
    try:local,occ,eigen,signature=_tile(model,*arrays)
    except Exception as exc:error=str(exc)
    _agree(comm,error)
    tiles=comm.allgather(signature)
    shape=(math.prod(model.kmesh),model.states,*reversed(model.grid))
    _validate_tiles(tiles,shape)
    waves=_assemble(comm,local,tiles,shape)
    archive=None;error=None
    if root:
        try:archive=_write(gs,destination,block,waves,occ,eigen,info)
        except Exception as exc:error=exc
    _root_error(comm,error)
    return archive
