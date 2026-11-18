"""Explicit GS/2 conversion and existing-GS import, with collective outcomes."""
import hashlib
import json
from pathlib import Path
import shlex
import sys
import warnings
from .errors import InputError,InputTypeError


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
    from ._gs_archive_mpi import _root_error
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
