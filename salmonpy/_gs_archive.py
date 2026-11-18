"""GS archive import/export internals shared by serial and retained-MPI owners.

The public interfaces remain gs_archive, gs_native and native_gs. This module
owns archive conversion/publication, not GS-to-TDDFT lifetimes or solver steps.
Import never initializes MPI, loads native libraries, compiles or launches jobs.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import warnings
import numpy as np
from .errors import InputError,InputTypeError,RuntimeFailure,NumericalError
from .native_gs import _agree
from .gs_native import _digest,_environment


# Native record compatibility and portable-to-native import.

def _record_flags(tokens):
    selected={}
    for token in tokens:
        if token.startswith(('-frecord-marker=','-fconvert=')):
            key,value=token.split('=',1)
            selected[key]=value
        elif token.startswith(('-fdefault-integer-','-fdefault-real-','-fdefault-double-',
                               '-finteger-4-integer-','-freal-4-real-','-freal-8-real-')):
            selected[token.rsplit('-',1)[0]]=token.rsplit('-',1)[1]
    if selected.get('-frecord-marker')=='4':selected.pop('-frecord-marker')
    if selected.get('-fconvert') in ('native',sys.byteorder+'-endian'):selected.pop('-fconvert')
    return tuple(key+'='+value for key,value in sorted(selected.items()))


def record_profile(record):
    tokens=list(record.get('link_command',()));unknown=False
    for option in record.get('configure_command',()):
        if option.startswith('-DCMAKE_Fortran_FLAGS') and '=' in option:
            try:tokens.extend(shlex.split(option.split('=',1)[1]))
            except ValueError:unknown=True
    return {'source_sha256':record.get('source_hashes',{}).get('src/io/checkpoint_restart.f90'),
            'compiler_sha256':record.get('compiler_sha256'),'record_flags':None if unknown else _record_flags(tokens)}


def _check_profile(profile,converter_record):
    uncertain=[]
    known=profile.get('record_flags')
    converted=_record_flags(converter_record.get('flags',()))
    if known is None:uncertain.append('native record flags are not fully recorded')
    elif tuple(known)!=converted:
        native_flags=dict(item.split('=',1) for item in known)
        converter_flags=dict(item.split('=',1) for item in converted)
        def layout(flags):
            return (flags.get('-frecord-marker','4'),
                    flags.get('-fdefault-integer',flags.get('-finteger-4-integer','4')))
        if layout(native_flags)!=layout(converter_flags):
            raise InputError('converter and native GS record format flags differ')
        uncertain.append('other native record-related flags differ')
    if not profile.get('compiler_sha256') or profile['compiler_sha256']!=converter_record.get('compiler_sha256'):
        uncertain.append('compiler binary match is not recorded/proven')
    if not profile.get('source_sha256') or profile['source_sha256']!=converter_record.get('native_source_sha256'):
        uncertain.append('native record-reader source match is not recorded/proven')
    if uncertain:
        warnings.warn('GS archive import: '+ '; '.join(uncertain)+
                      '; proceeding with the explicitly selected native loader',RuntimeWarning,stacklevel=3)


def load_archive(session,archive,directory,*,converter,executable,pseudopotentials,
                 material_id=None,timeout=300.):
    from .native_gs import _agree
    from .gs_archive import GroundStateArchive
    from .gs_native import GammaRestartConverter,_digest,_environment,_receipt_environment
    from .fdtd import SalmonExecutable
    from .dft import Pseudopotential
    from .config import _identifier,_freeze,_thaw
    from .ld import _positive
    session._live();comm=session.comm
    error=None
    try:
        path=archive.directory if isinstance(archive,GroundStateArchive) else Path(archive)
        path=path.expanduser().absolute();destination=Path(directory).expanduser().absolute()
        if material_id is not None:_identifier(material_id,'material_id')
        if not isinstance(converter,GammaRestartConverter):raise InputTypeError('explicit Gamma/SerialRestartConverter required')
        if not isinstance(executable,SalmonExecutable):raise InputTypeError('explicit paired serial SalmonExecutable required')
        potentials=tuple(pseudopotentials)
        if any(not isinstance(p,Pseudopotential) for p in potentials):raise InputTypeError('explicit Pseudopotential objects required')
        timeout=_positive(timeout,'timeout')
        descriptor=(str(path),str(destination),material_id,converter.sha256,converter.target_sha256,bool(converter._multi_k),
                    executable.sha256,executable.mpi_enabled,tuple((p.atomic_number,p.local_channel,p.sha256) for p in potentials),
                    timeout,_environment(),session.backend._io_profile)
    except Exception as exc:error=str(exc)
    _agree(comm,error)
    _agree(comm,'portable GS import selections differ' if any(v!=descriptor for v in comm.allgather(descriptor)) else None)
    origin=None;export_context=None;error=None
    if comm.rank==0:
        try:
            if destination.exists() or destination.is_symlink():raise FileExistsError(destination)
            source=GroundStateArchive.open(path)
            manifest_hash=_digest(source.directory/'manifest.json')
            selected_id=material_id
            if selected_id is None:selected_id=source.metadata.get('provenance',{}).get('material_id','material')
            _identifier(selected_id,'material_id')
            checked=type(converter).from_build(converter.build_record)
            if (checked.path,checked.sha256,checked.target_sha256)!=(converter.path,converter.sha256,converter.target_sha256):
                raise InputError('converter build selection changed')
            profile=json.loads(converter.build_record.read_text())
            _check_profile(session.backend._io_profile,profile)
            checked.write_bundle(source,destination,executable=executable,pseudopotentials=potentials,timeout=timeout)
            receipt=json.loads((destination/'salmonpy-native-gs.json').read_text())
            if not checked._multi_k:
                # GS/2 and its converter preserve complex amplitudes. A DFT
                # num_kgrid=1 input would instead allocate real orbitals. An
                # equivalent explicit Gamma table selects complex allocation.
                from .bands import _input_groups
                groups=_input_groups((destination/'rt.inp').read_text())
                (destination/'kpoints.dat').write_text('1\n1 0 0 0 1\n')
                groups['kgrid']="\n file_kw='kpoints.dat'\n"
                (destination/'rt.inp').write_text(''.join('&'+key+body+'/\n' for key,body in groups.items()))
                receipt.update(schema='salmonpy.gs-serial-bundle/1',kpoints_file='kpoints.dat',
                               kpoints_sha256=_digest(destination/'kpoints.dat'),
                               scope='explicit Gamma table for complex retained GS initialization')
                (destination/'salmonpy-native-gs.json').write_text(json.dumps(receipt,indent=2)+'\n')
            if receipt['portable_manifest_sha256']!=manifest_hash or _digest(source.directory/'manifest.json')!=manifest_hash:
                raise InputError('portable GS manifest changed during conversion')
            if _receipt_environment(destination,executable)!=_environment():
                raise InputError('converted native record environment differs')
            # Retain only small declared arrays/metadata, never source wavefunctions
            # or a live archive object. Re-export must work after source removal.
            export_context={'manifest_sha256':manifest_hash,'metadata':_thaw(source.metadata),
                            'arrays':{name:source.array(name).tolist() for name in
                                      ('cell_bohr','fractional_positions','atomic_numbers',
                                       'eigenvalues','kpoints','weights')}}
            if _digest(source.directory/'manifest.json')!=manifest_hash:
                raise InputError('portable GS manifest changed while capturing metadata')
            origin={'schema':'salmonpy.retained-gs-import/1','material_id':selected_id,
                    'archive_material_id':source.metadata.get('provenance',{}).get('material_id'),
                    'archive_converged':source.metadata.get('provenance',{}).get('converged'),
                    'explicit_gamma_complex_table':not checked._multi_k,
                    'archive_manifest_sha256':manifest_hash,'converter_sha256':checked.sha256,
                    'converter_runtime_sha256':checked.target_sha256,
                    'native_runtime_sha256':session.backend.sha256,
                    'scf_converged':None,'scope':'array initialization, not SCF mixing/native-state restoration'}
        except Exception as exc:error=exc
    _root_error(comm,error)
    if comm.size>1:origin,export_context=comm.bcast((origin,export_context),root=0)
    error=None
    try:frozen_origin=_freeze(origin);frozen_context=_freeze(export_context)
    except Exception as exc:error=str(exc)
    _agree(comm,error)
    state=session.load(destination.resolve(),'rt.inp',material_id=origin['material_id'])
    state._archive_origin=frozen_origin
    state._archive_export_context=frozen_context
    state._archive_export_potentials=potentials
    return state


# Collective retained-owner export.

def _root_error(comm,error):
    # Standard Python/SALMONpy input and I/O exceptions are picklable. Rank0 is
    # the only filesystem writer; communicate its outcome before the next phase.
    if comm.size>1:error=comm.bcast(error,root=0)
    if error is not None:raise error


def _nk(model):
    return model.n_kpoints if hasattr(model,'n_kpoints') else math.prod(model.kmesh)


def _model(gs):
    from .dft import DFT
    if isinstance(gs.calculation,DFT):return gs.calculation
    source=getattr(gs,'_archive_export_context',None)
    if source is None:
        raise InputError('export requires an immutable DFT model or a captured GS/2 import context')
    metadata=source['metadata'];arrays=source['arrays']
    return SimpleNamespace(archive_context=source,n_kpoints=metadata['n_kpoints'],
        states=metadata['n_bands'],grid=tuple(reversed(metadata['grid_zyx'])),
        numbers=arrays['atomic_numbers'],cell=arrays['cell_bohr'],
        fractional_positions=arrays['fractional_positions'],electrons=metadata['electrons'],
        xc=metadata['xc'],pseudopotentials=gs._archive_export_potentials)


def _description(model):
    source=getattr(model,'archive_context',None)
    return 'GS/2:'+source['manifest_sha256'] if source is not None else model.input_text()


def _tile(model,density,orbitals,occupations,eigenvalues):
    nk=_nk(model)
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


def _write(gs,destination,block,waves,occ,eigen,info,model):
    from .gs_archive import save_ground_state
    source=getattr(model,'archive_context',None)
    nk=_nk(model)
    rho=np.zeros(tuple(reversed(model.grid)))
    try:
        with np.errstate(over='raise',invalid='raise'):
            for k in range(nk):
                for band in range(model.states):
                    factor=occ[k,band]/nk if source is None else source['arrays']['weights'][k]*occ[k,band]
                    rho+=factor*(waves[k,band].real**2+waves[k,band].imag**2)
    except FloatingPointError as exc:
        raise NumericalError('nonfinite reconstructed retained GS density') from exc
    if source is not None:
        warnings.warn('SCF convergence unknown for imported GS; exporting source-labelled arrays',
                      RuntimeWarning,stacklevel=3)
    elif not info.converged:
        warnings.warn('SCF not converged; exporting the retained last iterate without '
                      'certifying physical convergence',RuntimeWarning,stacklevel=3)
    provenance=dict(schema='salmonpy.retained-gs-export/1',material_id=gs.material_id,
        state_identity=gs.identity,target_sha256=gs.session.backend.sha256,
        converged=info.converged,iterations=info.iterations,
        density_source='reconstructed_from_native_orbitals',
        local_channels=[dict(atomic_number=p.atomic_number,local_channel=p.local_channel)
                        for p in model.pseudopotentials])
    if source is None:
        references=[dict(atomic_number=p.atomic_number,filename=model._pseudo_name(i),
                         sha256=p.sha256,source='retained DFT model')
                    for i,p in enumerate(model.pseudopotentials)]
        # Same ordinary full-mesh order as native lattice.f90:init_kvector.
        mx,my,mz=model.kmesh
        points=[[(ix+.5)/mx-.5,(iy+.5)/my-.5,(iz+.5)/mz-.5]
                for iz in range(mz) for iy in range(my) for ix in range(mx)]
        weights=np.full(nk,1./nk);energy=info.energy_hartree
    else:
        inherited=source['metadata'].get('provenance',{})
        provenance.update(schema='salmonpy.retained-gs-reexport/1',converged=None,
            source_archive_manifest_sha256=source['manifest_sha256'],
            inherited_convergence=inherited.get('inherited_convergence',inherited.get('converged')),
            eigenvalue_source='source_archive',energy_source='source_archive',sampling_source='source_archive')
        references=source['metadata']['pseudopotentials']
        # save_ground_state accepts a detached JSON list, not frozen mappings.
        from .config import _thaw
        references=_thaw(references)
        points=source['arrays']['kpoints'];weights=source['arrays']['weights']
        energy=source['metadata']['energy_hartree']
    return save_ground_state(destination,wavefunctions=waves,density=rho,
        eigenvalues=eigen,occupations=occ,kpoints=points,weights=weights,
        cell_bohr=model.cell,fractional_positions=model.fractional_positions,
        atomic_numbers=list(model.numbers),electrons=model.electrons,energy_hartree=energy,
        xc=model.xc,pseudopotentials=references,orbital_block=block,provenance=provenance)


def export_retained(gs,directory,orbital_block):
    from .gs_archive import _integer
    gs._live()
    comm=gs.session.comm
    root=getattr(comm,'rank',0)==0
    error=None
    try:
        model=_model(gs)
        if math.prod(gs.session.backend.decomposition)!=comm.size:
            raise InputError('GS decomposition does not match subgroup size')
        block=_integer(orbital_block,'orbital_block')
        destination=Path(directory).expanduser().absolute()
        descriptor=(gs._token,gs.identity,gs.material_id,str(destination),block,
                    gs.session.backend.sha256,gs.session.backend.decomposition,
                    _description(model),tuple(p.sha256 for p in model.pseudopotentials))
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
    arrays=[gs.array(name) for name in ('density','orbitals','occupations')]
    source=getattr(model,'archive_context',None)
    if source is None:arrays.append(gs.array('eigenvalues'))
    else:
        from .native_gs import NativeGSArray
        carried=np.asarray(source['arrays']['eigenvalues'],dtype=np.float64).T[:,:,None]
        arrays.append(NativeGSArray(carried,(1,1,1),'eigenvalues','hartree (source archive)'))
    error=None
    try:local,occ,eigen,signature=_tile(model,*arrays)
    except Exception as exc:error=str(exc)
    _agree(comm,error)
    tiles=comm.allgather(signature)
    shape=(_nk(model),model.states,*reversed(model.grid))
    _validate_tiles(tiles,shape)
    waves=_assemble(comm,local,tiles,shape)
    archive=None;error=None
    if root:
        try:archive=_write(gs,destination,block,waves,occ,eigen,info,model)
        except Exception as exc:error=exc
    _root_error(comm,error)
    return archive


# Serial native-result export.

def _eigenvalues(text,nb,nk=1):
    lines=[line.strip() for line in text.splitlines() if line.strip()]
    if (len(lines)!=3+nk*(nb+1) or lines[:3]!=['#esp: single-particle energies (eigen energies)',
        '#occ: occupation numbers, io: orbital index','# 1:io, 2:esp[a.u.], 3:occ']):
        raise InputError('unpolarized a.u. eigenvalue table required')
    eigen=np.empty((nk,nb));occ=np.empty((nk,nb))
    try:
        for k in range(nk):
            offset=3+k*(nb+1)
            if re.fullmatch(r'k=\s*'+str(k+1)+r',\s*spin=\s*1',lines[offset]) is None:raise InputError('invalid eigenvalue k/spin order')
            for b,line in enumerate(lines[offset+1:offset+1+nb]):
                parts=line.split()
                if len(parts)!=3 or int(parts[0])!=b+1:raise InputError('invalid eigenvalue band index')
                eigen[k,b]=float(parts[1].replace('D','e').replace('d','e'))
                occ[k,b]=float(parts[2].replace('D','e').replace('d','e'))
    except (ValueError,OverflowError) as exc:raise InputError('invalid eigenvalue numbers') from exc
    if not np.isfinite(eigen).all() or not np.isfinite(occ).all():raise InputError('nonfinite eigenvalue table')
    return eigen,occ.ravel()


def _kpoints(text,nk,cell):
    lines=[s.strip() for s in text.splitlines() if s.strip()]
    if (len(lines)!=nk+7 or lines[0]!='# k-point distribution' or 'Reduced coordinate' not in lines[2]
            or '[none]' not in lines[4] or not lines[-2].startswith('# coefficients (2*pi/a [a.u.])')):
        raise InputError('orthogonal native k-point distribution required')
    points=np.empty((nk,3));weights=np.empty(nk)
    try:
        for k,line in enumerate(lines[5:5+nk]):
            parts=line.split()
            if len(parts)!=5 or int(parts[0])!=k+1:raise InputError('invalid k-point index/order')
            points[k]=[float(v) for v in parts[1:4]];weights[k]=float(parts[4])
        reciprocal=np.array([float(v) for v in lines[-1].split()])
    except (ValueError,OverflowError) as exc:raise InputError('invalid k-point distribution values') from exc
    if reciprocal.shape!=(3,) or not np.isfinite(reciprocal).all() or not np.allclose(reciprocal,2*np.pi/np.diag(cell),rtol=1e-14,atol=0):
        raise InputError('k-point reciprocal basis differs from cell')
    if not np.isfinite(points).all() or not np.isfinite(weights).all() or np.any(weights<0) or not np.isclose(weights.sum(),1.,rtol=0,atol=1e-12):
        raise InputError('invalid native k-point coordinates/weights')
    return points,weights


def _model_and_input(source,record,multi_k=False):
    from .dft import DFT,Pseudopotential,SCFSettings
    from .bands import _input_groups
    if not multi_k and record.get('kmesh')!=[1,1,1]:raise InputError('native export supports Gamma mesh only')
    potentials=[];refs=[]
    for i,p in enumerate(record['pseudopotentials']):
        candidates=[source.directory/'native'/f'pseudo{i+1}{suffix}' for suffix in ('_rps.dat','.upf')]
        candidates=[p for p in candidates if p.is_file()]
        if len(candidates)!=1:raise InputError('one DFT-generated potential file per species required')
        path=candidates[0];potential=Pseudopotential.from_file(path,atomic_number=p['atomic_number'],local_channel=p['local_channel'])
        if potential.sha256!=p['sha256']:raise InputError('potential bytes differ from GS metadata')
        potentials.append(potential);refs.append(dict(atomic_number=potential.atomic_number,filename=path.name,
            sha256=potential.sha256,source='verified native GS/1'))
    model=DFT(record['numbers'],record['cell_bohr'],record['fractional_positions'],potentials,
              record['electrons'],record['states'],record['grid'],record['kmesh'],record['xc'])
    if multi_k and (not np.array_equal(model.cell,np.diag(np.diag(model.cell))) or np.any(np.diag(model.cell)<=0)):
        raise InputError('multi-k native export requires aligned orthogonal cell')
    settings=SCFSettings(record['max_iterations'],record['tolerance'])
    expected=_input_groups(model.input_text(settings));actual=_input_groups((source.directory/'native/gs.inp').read_text())
    if set(actual)!=set(expected):raise InputError('only DFT-generated GS input groups supported')
    for name in expected:
        if name!='control' and actual[name].strip()!=expected[name].strip():raise InputError('GS metadata/input mismatch: '+name)
    control=expected['control'].strip()
    warm=control.replace("yn_restart='n'", "yn_restart='y', read_gs_restart_data='wfn', yn_reset_step_restart='y', directory_read_data='./initial_guess/'",1)
    if actual['control'].strip() not in (control,warm):raise InputError('unsupported native GS control settings')
    if _digest(source.directory/'native/gs.inp')!=record['input_sha256']:raise InputError('GS input hash differs from metadata')
    return model,refs


def export_gamma(converter,result,directory,*,orbital_block=16,timeout=300.):
    from .dft import GroundStateResult
    from .gs_archive import save_ground_state,_integer
    from .ld import _positive
    if not isinstance(result,GroundStateResult):raise InputTypeError('GroundStateResult required')
    orbital_block=_integer(orbital_block,'orbital_block');timeout=_positive(timeout,'timeout')
    if converter.exporter is None or converter.exporter_sha256 is None:raise RuntimeFailure('converter has no probed native exporter; explicitly build a matching new converter')
    if _digest(converter.exporter)!=converter.exporter_sha256:raise InputError('native exporter hash changed')
    build=json.loads(converter.build_record.read_text())
    if _digest(build['target_executable'])!=converter.target_sha256:raise InputError('paired runtime hash changed')
    source=GroundStateResult.open(result.directory);mp=source.directory/'ground-state.json'
    raw_manifest=mp.read_bytes();manifest_hash=hashlib.sha256(raw_manifest).hexdigest();record=json.loads(raw_manifest)
    if record.get('iterations')!=source.iterations or record.get('energy_hartree')!=source.energy_hartree:
        raise InputError('GS metadata changed during export validation')
    if record.get('executable_sha256')!=converter.target_sha256:raise InputError('GS and native exporter require the same runtime')
    try:model,refs=_model_and_input(source,record,converter._multi_k)
    except (KeyError,TypeError) as exc:raise InputError('incomplete DFT-generated GS metadata') from exc
    eigenpath=source.directory/'native/gs_eigen.data'
    if 'native/gs_eigen.data' not in record['files'] or _digest(eigenpath)!=record['files']['native/gs_eigen.data']:raise InputError('missing or changed native eigenvalue output')
    nk=int(np.prod(model.kmesh))
    points=np.zeros((1,3));weights=np.ones(1)
    if converter._multi_k and nk>1:
        kp=source.directory/'native/gs_k.data'
        if 'native/gs_k.data' not in record['files'] or _digest(kp)!=record['files']['native/gs_k.data']:raise InputError('missing or changed native k-point output')
        points,weights=_kpoints(kp.read_text(),nk,model.cell)
        mx,my,mz=model.kmesh
        expected=np.array([[(ix+.5)/mx-.5,(iy+.5)/my-.5,(iz+.5)/mz-.5]
                           for iz in range(mz) for iy in range(my) for ix in range(mx)])
        if not np.allclose(points,expected,rtol=0,atol=1e-14) or not np.allclose(weights,1./nk,rtol=1e-14,atol=0):
            raise InputError('k-point output differs from DFT-generated source mesh')
    eigen,eigen_occ=_eigenvalues(eigenpath.read_text(),model.states,nk)
    root=Path(directory).expanduser()
    if root.resolve().is_relative_to(source.directory):raise InputError('portable archive must be outside native GS result')
    if root.exists() or root.is_symlink():raise FileExistsError(root)
    environment=dict(os.environ);profile=_environment(environment)
    mappings=[]
    with tempfile.TemporaryDirectory(prefix='salmonpy-native-gs-export-') as tmp:
        work=Path(tmp);native=work/'native_input';native.mkdir()
        try:
            for name in ('info.bin','occupation.bin','wfn.bin'):
                relative='native/data_for_restart/'+name
                shutil.copy2(source.directory/relative,native/name)
                if _digest(native/name)!=record['files'][relative]:raise InputError('native GS changed while copying')
            nx,ny,nz=model.grid;nb=model.states;grid=(nz,ny,nx);ng=nx*ny*nz
            (work/'request.txt').write_text(f'{nx} {ny} {nz} {nb}\n{model.electrons}\n{source.iterations}\n'+(f'{nk}\n' if converter._multi_k else ''))
            if nk>1:(work/'portable_weights.bin').write_bytes(np.asarray(weights,dtype='<f8').tobytes())
            with (work/'exporter.log').open('w') as log:
                completed=subprocess.run([str(converter.exporter)],cwd=work,env=environment,stdout=log,stderr=subprocess.STDOUT,timeout=timeout,check=False)
            if completed.returncode:raise RuntimeFailure('native GS export failed: '+(work/'exporter.log').read_text()[-4000:])
            for name,count in [('portable_real.bin',ng*nb*nk),('portable_imag.bin',ng*nb*nk),('portable_occupations.bin',nb*nk)]:
                if (work/name).stat().st_size!=8*count:raise RuntimeFailure('invalid exported stream size')
            kind=(work/'native_kind.txt').read_text().strip()
            if kind not in ('real','complex'):raise RuntimeFailure('unknown exported native orbital kind')
            occ=np.fromfile(work/'portable_occupations.bin',dtype='<f8')
            if not np.array_equal(occ,eigen_occ):raise InputError('native occupation and eigenvalue table disagree')
            occ=occ.reshape(nk,nb)
            real=np.memmap(work/'portable_real.bin',dtype='<f8',mode='r',shape=(nk,nb,*grid));mappings.append(real)
            imag=np.memmap(work/'portable_imag.bin',dtype='<f8',mode='r',shape=(nk,nb,*grid));mappings.append(imag)
            waves=np.memmap(work/'complex-map.bin',dtype=np.complex128,mode='w+',shape=(nk,nb,*grid));mappings.append(waves)
            density=np.zeros(grid)
            try:
                with np.errstate(over='raise',invalid='raise'):
                    for k in range(nk):
                        for b in range(nb):
                            waves[k,b].real=real[k,b];waves[k,b].imag=imag[k,b]
                            density+=(weights[k]*occ[k,b])*(real[k,b]*real[k,b]+imag[k,b]*imag[k,b])
            except FloatingPointError as exc:raise InputError('nonfinite reconstructed density') from exc
            waves.flush()
            provenance=dict(schema='salmonpy.gs-native-export/1',source_manifest_sha256=manifest_hash,
                target_sha256=converter.target_sha256,exporter_sha256=converter.exporter_sha256,
                density_source='reconstructed_from_native_orbitals',eigenvalue_source='native/gs_eigen.data',
                native_real_orbital=kind=='real',local_channels=[dict(atomic_number=p.atomic_number,local_channel=p.local_channel) for p in model.pseudopotentials],
                record_environment=dict(profile))
            if nk>1:provenance['kpoint_source']='native/gs_k.data'
            return save_ground_state(root,wavefunctions=waves,density=density,eigenvalues=eigen,
                occupations=occ,kpoints=points,weights=weights,cell_bohr=model.cell,
                fractional_positions=model.fractional_positions,atomic_numbers=list(model.numbers),electrons=model.electrons,
                energy_hartree=source.energy_hartree,xc=model.xc,pseudopotentials=refs,orbital_block=orbital_block,provenance=provenance)
        finally:
            # Explicitly release mapping handles before TemporaryDirectory cleanup.
            for mapping in reversed(mappings):mapping._mmap.close()
