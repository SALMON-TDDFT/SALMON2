"""Immutable descriptions of declared GS inputs; never certify binary GS provenance.

Only explicit describe calls import optional f90nml. No MPI/native/filesystem reads.
"""
from dataclasses import dataclass
import hashlib
import json
import shlex
import numpy as np
from .errors import InputError
from .config import _freeze,_thaw,_identifier,_copy_json
from .bands import _input_groups
from .structure import _BOHR_ANGSTROM

__all__=['GSInputMetadata','describe_gs_input']


@dataclass(frozen=True)
class GSInputMetadata:
    _record: object

    def __post_init__(self):
        object.__setattr__(self,'_record',_freeze(_copy_json(_thaw(self._record))))

    @property
    def material_id(self):return self._record['material_id']
    @property
    def source_id(self):return self._record['source_id']
    @property
    def settings(self):return self._record['settings']
    @property
    def notes(self):return self._record['notes']
    def as_dict(self):return _thaw(self._record)


def _indexed(group,name):
    values=group.get(name,[])
    values=values if isinstance(values,list) else [values]
    start=group.start_index.get(name,[1])[0] if hasattr(group,'start_index') else 1
    return {start+i:value for i,value in enumerate(values) if value is not None}


def _vector(value,shape):
    array=np.asarray(value,dtype=float)
    if array.shape!=shape or not np.isfinite(array).all():raise ValueError('invalid numeric shape/value')
    return array


def _geometry(groups,system,pseudo,assets,scale,notes):
    cell=None;numbers=None;positions=None
    try:
        keys=['al_vec'+str(i) for i in range(1,4)]
        vectors=_vector([system[k] for k in keys],(3,3)) if all(k in system for k in keys) else None
        lengths=_vector(system['al'],(3,)) if 'al' in system else None
        if vectors is not None and lengths is not None and not np.array_equal(vectors,np.diag(lengths)):
            raise ValueError('ambiguous cell declarations')
        declared=vectors if vectors is not None else np.diag(lengths) if lengths is not None else None
        if declared is None or scale is None or abs(np.linalg.det(declared))<1e-12:
            raise ValueError('unknown cell/length units')
        cell=declared*scale
    except (ValueError,TypeError,KeyError,np.linalg.LinAlgError):
        notes.append('cell_bohr unknown: explicit usable cell and length units are required')
    sources=[]
    for group,field,reduced in [('atomic_coor','file_atom_coor',False),('atomic_red_coor','file_atom_red_coor',True)]:
        if group in groups:sources.append((groups[group],reduced))
        filename=system.get(field)
        if filename is not None and filename!='none':
            if filename in assets:
                try:sources.append((assets[filename].decode('utf-8'),reduced))
                except UnicodeError:notes.append('coordinate asset encoding unknown: '+str(filename))
            else:notes.append('coordinate asset unavailable: '+str(filename))
    try:
        if len(sources)!=1:raise ValueError('coordinate source is not unique')
        rows=[];zs=[];mapping=_indexed(pseudo,'izatom')
        for line in sources[0][0].splitlines():
            lexer=shlex.shlex(line,posix=True);lexer.whitespace_split=True;lexer.commenters='!'
            fields=list(lexer)
            if not fields:continue
            if len(fields)!=5:raise ValueError('unsupported atomic record')
            rows.append([float(x.replace('D','e').replace('d','e')) for x in fields[1:4]])
            z=mapping[int(fields[4])]
            if type(z) is not int or not 1<=z<=118:raise ValueError('invalid species')
            zs.append(z)
        values=_vector(rows,(system['natom'],3))
        numbers=zs
        if sources[0][1]:positions=values
        elif scale is not None and cell is not None:positions=np.linalg.solve(cell.T,(values*scale).T).T
        else:notes.append('fractional_positions unknown: Cartesian lengths/cell unresolved')
    except (ValueError,TypeError,KeyError,np.linalg.LinAlgError):
        numbers=None;positions=None
        notes.append('atomic_numbers/fractional_positions unknown: coordinate/species declarations incomplete')
    return (None if cell is None else cell.tolist(),numbers,None if positions is None else positions.tolist())


def _k_sampling(group,assets,notes):
    filename=group.get('file_kw')
    if filename is not None and filename!='none':
        result={'mode':'file','filename':filename,'input_kpoints':None,'input_weights':None}
        try:
            raw=assets[filename]
            result['sha256']=hashlib.sha256(raw).hexdigest()
            rows=[]
            for line in raw.decode('utf-8').splitlines():
                fields=line.split('!',1)[0].split()
                if fields:rows.append(fields)
            count=int(rows[0][0])
            if count<1 or len(rows)!=count+1:raise ValueError('invalid table length')
            points=[];weights=[]
            for i,row in enumerate(rows[1:]):
                if len(row)!=5 or int(row[0])!=i+1:raise ValueError('unsupported table row')
                values=[float(v.replace('D','e').replace('d','e')) for v in row[1:]]
                if not np.isfinite(values).all():raise ValueError('nonfinite k data')
                points.append(values[:3]);weights.append(values[3])
            result['input_kpoints']=points;result['input_weights']=weights
        except (ValueError,TypeError,KeyError,IndexError,UnicodeError):
            notes.append('k table values unknown: unavailable/unsupported explicit table')
        return result
    if 'num_kgrid' not in group:return None
    return {'mode':'regular','mesh':group['num_kgrid'],'declared_shift':group.get('dk_shift'),
            'declared_gamma_centered':group.get('yn_gamma_centered')}


def describe_gs_input(text,*,assets,material_id='material',origin='declared',converged=None):
    """Describe captured input bytes. Unknown/default/effective native values stay unknown.

    assets maps relative input names to immutable bytes. This reads no paths,
    runs no solver and does not bind a DFT model or authorize GS restoration.
    source_id fingerprints input/assets, not physical equivalence or authenticity.
    """
    try:import f90nml
    except ImportError as exc:
        raise ImportError('GS input descriptions require salmonpy[metadata] (f90nml)') from exc
    _identifier(material_id,'material_id')
    if not isinstance(text,str):raise InputError('input text required')
    if origin not in ('declared','computed','imported','unknown') or (converged is not None and type(converged) is not bool):
        raise InputError('invalid GS metadata evidence')
    assets=dict(assets)
    if any(type(k) is not str or not isinstance(v,bytes) for k,v in assets.items()):
        raise InputError('assets must map names to immutable bytes')
    groups=_input_groups(text)
    try:
        nml=f90nml.reads(''.join('&'+name+body+'/\n' for name,body in groups.items()
                              if name not in ('atomic_coor','atomic_red_coor')))
    except (ValueError,TypeError,AssertionError) as exc:raise InputError('cannot describe GS namelist') from exc
    system=nml.get('system',{});pseudo=nml.get('pseudo',{});notes=[]
    units=nml.get('units',{}).get('unit_system')
    scale=1. if units in ('au','a.u.') else 1./_BOHR_ANGSTROM if str(units).lower()=='a_ev_fs' else None
    if scale is None:notes.append('length unit unknown; no implicit native default assumed')
    cell,numbers,positions=_geometry(groups,system,pseudo,assets,scale,notes)
    filenames=_indexed(pseudo,'file_pseudo');channels=_indexed(pseudo,'lloc_ps')
    potentials=[]
    for index,z in sorted(_indexed(pseudo,'izatom').items()):
        filename=filenames.get(index);raw=assets.get(filename)
        potentials.append({'species_index':index,'atomic_number':z,'local_channel':channels.get(index),
                           'filename':filename,'sha256':None if raw is None else hashlib.sha256(raw).hexdigest()})
        if raw is None:notes.append('pseudopotential bytes unavailable for species '+str(index))
    grid=nml.get('rgrid',{})
    spacing=grid.get('dl')
    try:spacing=None if spacing is None or scale is None else (_vector(spacing,(3,))*scale).tolist()
    except (ValueError,TypeError):spacing=None;notes.append('grid spacing unknown')
    functional={str(k):v for k,v in nml.get('functional',{}).items()}
    if 'xc' in functional and isinstance(functional['xc'],str):functional['xc']=functional['xc'].lower()
    if not potentials:notes.append('pseudopotential species mapping not declared')
    if 'xc' not in functional:notes.append('XC default/effective selection unknown')
    settings={'unit_system':units,'periodic':system.get('yn_periodic'),'spin':system.get('spin'),
        'cell_bohr':cell,'atomic_numbers':numbers,'fractional_positions':positions,
        'electrons':system.get('nelec'),'states':system.get('nstate'),
        'grid_counts':grid.get('num_rgrid'),'grid_spacing_bohr':spacing,
        'k_sampling':_k_sampling(nml.get('kgrid',{}),assets,notes),
        'functional':functional,'pseudopotentials':potentials}
    hashes={name:hashlib.sha256(raw).hexdigest() for name,raw in sorted(assets.items())}
    text_hash=hashlib.sha256(text.encode('utf-8')).hexdigest()
    identity=hashlib.sha256(json.dumps([text_hash,hashes],sort_keys=True).encode()).hexdigest()
    record={'schema':'salmonpy.gs-input/1','material_id':material_id,'source_id':identity,
        'input_sha256':text_hash,'assets_sha256':hashes,'settings':settings,'notes':notes,
        'evidence':{'configuration':'declared input, not effective native state','origin':origin,
                    'scf_converged':converged,'checkpoint_provenance':'unknown'}}
    return GSInputMetadata(record)
