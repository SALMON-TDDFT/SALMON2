"""Structural admission of bounded native fast-PR inputs; not GS provenance.

Legacy records: native-endian4-byte integer/logical/record markers, binary64,
one global real64/complex128 wfn stream. Other binary layouts reject explicitly.
"""
import itertools
import math
from pathlib import Path
import re
import struct
import warnings

import numpy as np

from .errors import InputError
from ._hse_geometry import _assignments, _indices, _numeric, _string, _values, _INT_MAX

_INT = struct.Struct('=i')
_CHUNK = 131072  # binary64 values; bounded reads, no whole-wavefunction copy
_QUOTED = re.compile(r"(?:([0-9]+)\s*\*\s*)?('(?:''|[^'])*'|\"(?:\"\"|[^\"])*\")")


def _fail(reason):
    raise InputError('HSE assets: '+reason)


def _path(root,name):
    p=Path(name)
    if p.is_absolute() or '..' in p.parts or '\\' in name or not name:
        _fail('asset references must be bundle-relative')
    p=(root/p).resolve()
    if not p.is_relative_to(root.resolve()):_fail('asset escapes bundle')
    return p


def _string_array(body,name,size):
    # Do not allocate from an untrusted species count without corresponding text.
    if size>len(body)+1:_fail('missing pseudo entries for declared species')
    result=[None]*size
    for key,index,text in _assignments(body):
        if key!=name:continue
        positions=_indices(index,size,name)
        values=[];end=0
        for match in _QUOTED.finditer(text):
            gap=text[end:match.start()]
            if gap.strip(' ,\t\r\n') or re.search(r',\s*,',gap) or (not values and ',' in gap):
                _fail('invalid pseudo filename list')
            count=int(match[1] or 1)
            if not 1<=count<=len(positions)-len(values):_fail('too many pseudo filenames')
            raw=match[2];values.extend([raw[1:-1].replace(raw[0]*2,raw[0])]*count);end=match.end()
        if not values or text[end:].strip():_fail('invalid pseudo filename literal')
        for i,value in zip(positions,values):result[i]=value
    if any(value is None or value=='none' for value in result):_fail('missing pseudo filename for a species')
    return result


def _ky_record(stream,kinds):
    """List-directed record: allow continuation and ignore unused trailing fields."""
    values=[]
    while len(values)<len(kinds):
        row=stream.readline()
        if not row:_fail('pseudo KY numeric records are truncated')
        for token in row.split('!',1)[0].replace(',',' ').split():
            repeat=re.fullmatch(r'([0-9]+)\*(.+)',token)
            count=int(repeat[1]) if repeat else 1
            literal=repeat[2] if repeat else token
            if count<1:_fail('pseudo KY invalid repetition')
            for _ in range(min(count,len(kinds)-len(values))):
                value=literal
                if not kinds[len(values)]:
                    implicit=re.fullmatch(r'([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+))([+-][0-9]+)',value)
                    if implicit:value=implicit[1]+'e'+implicit[2]
                values.append(_values(value,1,kinds[len(values)])[0])
            if len(values)==len(kinds):break
    return values


def _ky(path,lmax_input,local):
    """Plain scalar KY records following native list-directed consumption."""
    try:
        with path.open() as stream:
            mr,step,angular,charge=_ky_record(stream,[True,False,True,False])
            if not 1<=mr<=50000 or not 0<=angular<=8 or step<=0 or charge<=0:
                _fail('pseudo KY header exceeds supported native bounds')
            if lmax_input < -1 or lmax_input>angular or local < -1 or local> (angular if lmax_input<0 else lmax_input):
                _fail('pseudo KY angular/local channel mismatch')
            radii=_ky_record(stream,[True]+[False]*(angular+1))[1:]
            if any(v<0 for v in radii):
                warnings.warn('HSE assets: negative pseudo KY cutoff radius; continuing', RuntimeWarning, stacklevel=2)
            for width in (angular+3,angular+2):
                for _ in range(mr+1):_ky_record(stream,[False]*width)
    except (OSError,UnicodeError,ValueError) as exc:
        raise InputError('HSE assets: pseudo KY numeric records invalid: '+str(exc)) from exc


def _pseudopotentials(root,groups):
    nelem=int(_numeric(groups.get('system',''),'nelem',integer=True)[0])
    body=groups.get('pseudo','');names=_string_array(body,'file_pseudo',nelem)
    numbers=_numeric(body,'izatom',nelem,True,default=-1,array=True)
    # The bundled native mass/symbol table stops above Bi.
    if np.any(numbers<1) or np.any(numbers>83):_fail('pseudo atomic numbers must be supported by native table (1..83)')
    lmax=_numeric(body,'lmax_ps',nelem,True,default=-1,array=True)
    local=_numeric(body,'lloc_ps',nelem,True,default=-1,array=True)
    for i,name in enumerate(names):
        path=_path(root,name)
        try:
            with path.open('rb') as stream:
                if not stream.read(1):_fail('pseudo file is empty')
        except OSError as exc:raise InputError('HSE assets: pseudo file cannot be read: '+name) from exc
        if name.endswith('_rps.dat'):_ky(path,int(lmax[i]),int(local[i]))
        elif not name.endswith(('.pspnc','.psp8','.fhi','.cpi','.vps','.upf','.UPF')):
            _fail('pseudo suffix is not recognized by native loader')
        # Non-KY content parsing is intentionally still the native loader's job.


def _record(stream,length):
    marker=stream.read(4)
    if len(marker)!=4 or _INT.unpack(marker)[0]!=length:
        _fail('GS record marker/layout mismatch (native-endian int32 records required)')
    payload=stream.read(length)
    if len(payload)!=length or stream.read(4)!=marker:_fail('GS record is truncated or inconsistent')
    return payload


def _info(path):
    try:
        with path.open('rb') as stream:
            values=[_INT.unpack(_record(stream,4))[0] for _ in range(5)]
            if stream.read(1):_fail('GS info has unexpected trailing records')
    except OSError as exc:raise InputError('HSE assets: GS info cannot be read') from exc
    nk,no,step,producers,real=values
    if nk<1 or no<1 or step<0 or producers<1 or real not in (0,1):_fail('GS metadata contains invalid counts/flags')
    return nk,no,bool(real)


def _mesh(groups,operations):
    body=groups.get('kgrid','')
    if _string(body,'file_kw','none')!='none':_fail('k mesh requires generated uniform points; file_kw is unsupported by fast-PR native profile')
    if _string(groups.get('unfolding',''),'dm_unfold_option','no')!='no':_fail('k unfolding is unsupported by this profile')
    mesh=_numeric(body,'num_kgrid',3,True,default=1)
    shift=_numeric(body,'dk_shift',3)
    if np.any(mesh<1) or math.prod(int(v) for v in mesh)>_INT_MAX:_fail('k mesh counts exceed native bounds')
    gamma=_string(body,'yn_gamma_centered','n')
    if gamma not in ('n','y'):_fail('invalid k gamma-centering flag')
    if gamma=='y' and np.any(mesh%2==0):_fail('even gamma-centered k endpoints are incompatible with fast-PR full-q mesh')
    if gamma=='y':shift=np.zeros(3)  # Native odd Gamma branch ignores dk_shift.
    if np.max(abs(shift))>_INT_MAX:_fail('k shift exceeds native integer mapping bounds')
    reduction=_numeric(groups.get('functional',''),'hse_q_reduction',3,True,default=1)
    if np.any(reduction!=1):_fail('k/q profile requires full-q mesh')
    matrices=[];offsets=[]
    if operations is not None:
        # Native representative selection merges delta^2 <1e-5; distinct mesh
        # neighbors must stay outside that tolerance before strict star checks.
        if np.any(1/mesh.astype(float)**2<1e-5):_fail('k mesh is finer than native representative separation tolerance')
        origin=(.5+shift)/mesh-.5
        for operation in operations:
            r=operation[:,:3]
            if np.max(abs(r[:,2]-[0.,0.,1.]))>1e-10:
                _fail('retained k symmetry must preserve the native z-field initialization group')
            matrix=r*mesh[:,None]/mesh[None,:]
            offset=mesh*(r@origin-origin)
            if not np.isfinite(matrix).all() or not np.isfinite(offset).all():_fail('nonfinite k-star mapping')
            rounded=np.rint(matrix);base=np.rint(offset)
            for corner in itertools.product(*[(0,int(n)-1) for n in mesh]):
                mapped=matrix@corner+offset
                if (not np.isfinite(mapped).all() or np.max(abs(mapped-(rounded@corner+base)))>1e-8
                        or np.any(np.rint(mapped)<-_INT_MAX-1) or np.any(np.rint(mapped)>_INT_MAX)):
                    _fail('k shift/mesh is incompatible with retained symmetry stars')
            matrices.append(rounded.astype(np.int64));offsets.append(base.astype(np.int64))
    return mesh,matrices,offsets


def _weights(mesh,matrices,offsets,nk,group_size):
    total=math.prod(int(v) for v in mesh)
    if nk<group_size:_fail('k representatives leave an empty solver rank')
    if not matrices:
        if nk!=total:_fail('GS k count differs from generated full mesh')
        return None  # uniform; avoid a duplicate array for the unreduced case
    if nk>total or nk*len(matrices)<total:_fail('GS k representative count cannot cover generated stars')
    visited=bytearray(total);weights=[]
    for flat in range(total):
        if visited[flat]:continue
        point=np.array([flat%mesh[0],(flat//mesh[0])%mesh[1],flat//(mesh[0]*mesh[1])],dtype=np.int64)
        orbit=set()
        for matrix,base in zip(matrices,offsets):
            p=(matrix@point+base)%mesh
            orbit.add(int(p[0]+mesh[0]*(p[1]+mesh[1]*p[2])))
        if any(visited[i] for i in orbit):_fail('k stars overlap inconsistently')
        for i in orbit:visited[i]=1
        weights.append(len(orbit)/total)
        if len(weights)>nk:_fail('GS k count differs from generated representatives')
    if len(weights)!=nk or abs(math.fsum(weights)-1)>1e-10:_fail('k representative weights/counts inconsistent')
    return weights


def _wavefunction(path,length):
    try:
        if path.stat().st_size!=length:_fail('wavefunction byte length differs from GS/grid metadata')
        with path.open('rb') as stream:
            remaining=length
            while remaining:
                data=stream.read(min(remaining,8*_CHUNK))
                if not data or len(data)%8:_fail('wavefunction stream truncated')
                if not np.isfinite(np.frombuffer(data,dtype='=f8')).all():_fail('wavefunction contains nonfinite values')
                remaining-=len(data)
            if stream.read(1):_fail('wavefunction stream changed length')
    except OSError as exc:raise InputError('HSE assets: wavefunction cannot be read') from exc


def _occupations(path,nk,stored_states,states,electrons,weights):
    length=8*nk*stored_states
    if length>_INT_MAX:_fail('GS occupation record continuation layout is unsupported')
    try:
        if path.stat().st_size!=length+8:_fail('GS occupation record length mismatch')
        with path.open('rb') as stream:
            marker=stream.read(4)
            if len(marker)!=4 or _INT.unpack(marker)[0]!=length:_fail('GS occupation marker/layout mismatch')
            count=0.; outside_range=False
            for ik in range(nk):
                column=0.
                for first in range(0,stored_states,_CHUNK):
                    n=min(_CHUNK,stored_states-first);data=stream.read(8*n)
                    if len(data)!=8*n:_fail('GS occupation stream truncated')
                    values=np.frombuffer(data,dtype='=f8')
                    if not np.isfinite(values).all():
                        _fail('GS occupations must be finite')
                    outside_range |= bool(np.any(values < -1e-12) or np.any(values>2+1e-12))
                    column+=float(np.sum(values[:max(0,min(n,states-first))]))
                count+=column*(1/nk if weights is None else weights[ik])
            if stream.read(4)!=marker or stream.read(1):_fail('GS occupation record ending mismatch')
            if outside_range:
                warnings.warn('HSE assets: GS occupations outside [0,2]; continuing without modifying them', RuntimeWarning, stacklevel=2)
            if not math.isclose(count,electrons,rel_tol=2e-12,abs_tol=1e-8):
                warnings.warn(f'HSE assets: weighted occupations {count:.12g} differ from nelec={electrons}; continuing without renormalization', RuntimeWarning, stacklevel=2)
    except OSError as exc:raise InputError('HSE assets: GS occupations cannot be read') from exc


def preflight(root,groups,operations,grid,group_size):
    """Structural checks only; legacy records cannot authenticate physical GS identity."""
    root=Path(root)
    system=groups.get('system','')
    if _string(system,'spin','auto').lower() not in ('auto','unpolarized') or _string(system,'yn_spinorbit','a').lower() not in ('a','n'):
        _fail('GS layout admission requires the scalar unpolarized native profile')
    if _string(groups.get('control',''),'method_wf_distributor','single')!='single':
        _fail('GS admission supports the single global wavefunction layout only')
    mesh,matrices,offsets=_mesh(groups,operations)
    _pseudopotentials(root,groups)
    directory=_string(groups.get('control',''),'directory_read_data','restart/')
    restart=_path(root,directory)
    nk,stored_states,real=_info(restart/'info.bin')
    states=int(_numeric(groups.get('system',''),'nstate',integer=True)[0])
    electrons=int(_numeric(groups.get('system',''),'nelec',integer=True)[0])
    if states<1 or states>stored_states or electrons<1 or electrons>2*states:
        _fail('GS/input state or electron count mismatch')
    length=math.prod(int(v) for v in grid)*stored_states*nk*(8 if real else 16)
    try:
        if (restart/'wfn.bin').stat().st_size!=length:_fail('wavefunction byte length differs from GS/grid metadata')
    except OSError as exc:raise InputError('HSE assets: wavefunction cannot be inspected') from exc
    weights=_weights(mesh,matrices,offsets,nk,group_size)
    _occupations(restart/'occupation.bin',nk,stored_states,states,electrons,weights)
    # File size is checked above. Do not reread potentially huge orbitals at
    # every state creation; _wavefunction is available for explicit diagnostics.
