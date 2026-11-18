"""Read-only, dependency-light geometry admission for native periodic fast-PR HSE.

Not a general namelist parser, restart validator, or physical convergence check.
Formula/priority references: resolve_cell/init_grid_whole_config, symmetry_init,
and symmetry_validate_atoms in the bundled fast-PR solver.
"""
import math
from pathlib import Path
import re
import shlex

import numpy as np

from .errors import InputError

_INT_MAX = 2**31-1
_REAL = re.compile(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eEdD][+-]?[0-9]+)?')
_TOKEN = re.compile(r"'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"|![^\n]*")
_ASSIGN = re.compile(r'\b([a-zA-Z][a-zA-Z0-9_]*)\s*(\([^)]*\))?\s*=')


def _fail(reason):
    raise InputError('HSE geometry: '+reason)


def _assignments(body):
    clean = _TOKEN.sub(lambda m: ' '*len(m[0]) if m[0].startswith('!') else m[0], body)
    masked = _TOKEN.sub(lambda m: ' '*len(m[0]), clean)
    found = list(_ASSIGN.finditer(masked))
    result=[]
    for i,m in enumerate(found):
        text=clean[m.end():found[i+1].start() if i+1<len(found) else len(clean)].strip()
        # Consume one assignment separator, leaving extra null delimiters for
        # rejection rather than silently treating them as ordinary values.
        if text.endswith(','):text=text[:-1].rstrip()
        result.append((m[1].lower(),m[2],text))
    return result


def _values(text, limit, integer=False):
    if re.search(r',\s*,',text) or text.startswith(',') or text.endswith(',') or not text:
        _fail('null/empty numeric values are unsupported')
    values = []
    for token in text.replace(',',' ').split():
        count = 1
        if '*' in token:
            repeat,token = token.split('*',1)
            if re.fullmatch('[0-9]+',repeat) is None:
                _fail('invalid repeat count')
            count = int(repeat)
        if not 1 <= count <= limit-len(values):
            _fail('too many numeric values')
        pattern = r'[+-]?[0-9]+' if integer else _REAL
        if re.fullmatch(pattern,token) is None:
            _fail('invalid '+('integer' if integer else 'real')+' value: '+token)
        number = int(token) if integer else float(token.replace('D','e').replace('d','e'))
        if not math.isfinite(number) or (integer and abs(number)>_INT_MAX):
            _fail('nonfinite/out-of-range numeric value')
        values.extend([number]*count)
    return values


def _indices(index,size,name):
    if index is None:
        positions = list(range(size))
    else:
        selector = index[1:-1].strip()
        if re.fullmatch('[0-9]+',selector):
            start=int(selector)-1
            if not 0<=start<size:_fail(name+' index outside array')
            positions=list(range(start,size))
        else:
            parts=selector.split(':')
            if len(parts) not in (2,3) or any(x.strip() and re.fullmatch('[0-9]+',x.strip()) is None for x in parts):
                _fail(name+' has unsupported array selector')
            start=int(parts[0]) if parts[0].strip() else 1
            stop=int(parts[1]) if parts[1].strip() else size
            step=int(parts[2]) if len(parts)==3 and parts[2].strip() else 1
            if not 1<=start<=size or not 1<=stop<=size or step<1 or start>stop:
                _fail(name+' invalid array section')
            positions=list(range(start-1,stop,step))
    return positions


def _numeric(body,name,size=1,integer=False,default=0,array=False):
    result = np.full(size,default,dtype=np.int64 if integer else float)
    for key,index,text in _assignments(body):
        if key != name:
            continue
        if index is not None and size==1 and not array:_fail(name+' is scalar')
        positions=_indices(index,size,name)
        values=_values(text,len(positions),integer)
        result[positions[:len(values)]]=values
    return result


def _string(body,name,default):
    value=default
    for key,index,text in _assignments(body):
        if key != name:continue
        if index is not None or len(text)<2 or text[0] not in "'\"" or text[-1]!=text[0]:
            _fail(name+' requires a quoted scalar')
        value=text[1:-1].replace(text[0]*2,text[0])
    return value


def _atoms(root,groups,cell,natom,nelem):
    system=groups.get('system','');sources=[]
    for group,field,reduced in [('atomic_coor','file_atom_coor',False),('atomic_red_coor','file_atom_red_coor',True)]:
        if group in groups:sources.append((groups[group],reduced))
        name=_string(system,field,'none')
        if name!='none':
            path=Path(name)
            if path.is_absolute() or '..' in path.parts or '\\' in name:
                _fail('coordinate file must be bundle-relative')
            resolved=(root/path).resolve()
            if not resolved.is_relative_to(root.resolve()):_fail('coordinate file escapes bundle')
            try:text=resolved.read_text()
            except (OSError,UnicodeError) as exc:raise InputError('HSE geometry: cannot read coordinate file') from exc
            sources.append((text,reduced))
    if len(sources)!=1:_fail('require exactly one explicit atom coordinate source')
    text,reduced=sources[0];positions=[];species=[]
    for line in text.splitlines():
        clean=_TOKEN.sub(lambda m: '' if m[0].startswith('!') else m[0],line)
        masked=_TOKEN.sub(lambda m: ' '*len(m[0]),clean)
        if re.search(r',\s*,',masked) or clean.lstrip().startswith(','):
            _fail('null atom fields are unsupported')
        lexer=shlex.shlex(line,posix=False);lexer.whitespace_split=True;lexer.whitespace+=',';lexer.commenters='!'
        try:fields=list(lexer)
        except ValueError as exc:raise InputError('HSE geometry: invalid atom record') from exc
        if not fields:continue
        if len(fields)!=5:_fail('atom records require label, three coordinates and species ID')
        positions.append(_values(' '.join(fields[1:4]),3))
        if len(positions[-1])!=3:_fail('atom record requires three coordinates')
        kind=_values(fields[4],1,True)[0]
        if not 1<=kind<=nelem:_fail('atom species ID outside 1..nelem')
        species.append(kind)
    if len(positions)!=natom:_fail('natom differs from atom coordinate records')
    positions=np.asarray(positions,dtype=float)
    with np.errstate(over='ignore',invalid='ignore'):
        cartesian=positions@cell.T if reduced else positions
        fractional=np.linalg.solve(cell,cartesian.T).T
    if not np.isfinite(cartesian).all() or not np.isfinite(fractional).all():_fail('nonfinite atom coordinates')
    return fractional,np.asarray(species)


def preflight(root,groups,operations=None):
    """Validate the effective cell/grid and atoms without changing any input."""
    root=Path(root);system=groups.get('system','');rgrid=groups.get('rgrid','')
    if _string(system,'yn_periodic','n')!='y':_fail('native fast-PR requires a periodic cell')
    units=_string(groups.get('units',''),'unit_system','au')
    if units!='au':_fail('native fast-PR requires unit_system=au')
    al=_numeric(system,'al',3)
    vectors=np.column_stack([_numeric(system,'al_vec'+str(i),3) for i in range(1,4)])
    dl=_numeric(rgrid,'dl',3);points=_numeric(rgrid,'num_rgrid',3,True)
    if not all(np.isfinite(a).all() for a in (al,vectors,dl)):_fail('nonfinite cell/grid')
    spacing_used=bool(np.any(dl!=0))
    if spacing_used and np.any(dl<=0):_fail('grid spacing must be positive')
    if not spacing_used and np.any(points<=0):_fail('grid point counts must be positive')
    off=vectors-np.diag(np.diag(vectors))
    if not np.any(off):
        lengths=al.copy()
        if np.any(lengths==0):
            if np.all(np.diag(vectors)!=0):lengths=np.diag(vectors).copy()
            elif np.all(points>0) and np.all(dl!=0):lengths=points*dl
            else:_fail('cell lengths are not defined')
        cell=np.diag(lengths)
    else:
        # Full native TDDFT initialization is not the independent exchange-core
        # cell test. Axis-permuted fixtures hit ACE init STOP in the current
        # runtime; decline before lifetime rather than normalizing the basis.
        _fail('native fast-PR initialization currently requires an axis-aligned orthogonal cell')
    if not np.isfinite(cell).all() or not np.isfinite(lengths).all() or np.any(lengths<=0):
        _fail('cell lengths must be finite and positive')
    magnitude=float(np.max(abs(cell)))
    determinant=float(np.linalg.det(cell/magnitude))
    if not np.isfinite(determinant) or determinant==0 or (operations is not None and abs(determinant)<=1e-14):
        _fail('cell is singular or unsupported')
    with np.errstate(over='ignore',under='ignore',invalid='ignore'):
        volume=abs(float(np.linalg.det(cell)))
    if not math.isfinite(volume) or volume<=0:_fail('cell volume outside finite native range')
    try:inverse=np.linalg.inv(cell)
    except np.linalg.LinAlgError as exc:raise InputError('HSE geometry: singular cell') from exc
    if not np.isfinite(inverse).all() or np.max(abs(cell@inverse-np.eye(3)))>1e-10:
        _fail('cell inverse is not numerically consistent')
    spacing=dl if spacing_used else lengths/points
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        ratio=(lengths+1e-10)/spacing
    if not np.isfinite(ratio).all() or np.any(ratio>_INT_MAX-10) or np.any(ratio<1):
        _fail('grid dimensions outside native integer range')
    grid=np.trunc(ratio).astype(np.int64)
    if math.prod(int(v) for v in grid)>_INT_MAX:_fail('grid total exceeds native integer range')
    if np.sum(abs(dl))<=1e-12:
        if np.any(grid!=points):_fail('grid counts do not match native cell resolution')
    elif np.max(abs(lengths/spacing-grid))>1e-4:
        _fail('grid spacing does not tile the cell')
    natom=int(_numeric(system,'natom',integer=True)[0]);nelem=int(_numeric(system,'nelem',integer=True)[0])
    if natom<1 or nelem<1:_fail('positive natom and nelem required')
    fractional,species=_atoms(root,groups,cell,natom,nelem)
    grid.setflags(write=False)
    if operations is None:return grid
    for operation in operations:
        rotation=operation[:,:3];translation=operation[:,3]
        cartesian=cell@rotation@inverse
        if not np.isfinite(cartesian).all() or np.max(abs(cartesian.T@cartesian-np.eye(3)))>1e-10:
            _fail('symmetry does not preserve the cell metric')
        # Affine grid-index map M*i+b. Integer generators and unimodular
        # rotation establish a bijection without allocating the whole grid.
        matrix=rotation.T*grid[:,None]/grid[None,:]
        offset=-(rotation.T@translation)*grid
        if not np.isfinite(matrix).all() or not np.isfinite(offset).all():
            _fail('nonfinite symmetry grid mapping')
        residual=matrix-np.rint(matrix);constant=offset-np.rint(offset)
        # Bound the affine residual over the complete box, not merely one
        # generator: tiny rotation errors can accumulate on a large grid.
        for bits in range(8):
            corner=np.array([grid[i]-1 if bits&(1<<i) else 0 for i in range(3)])
            mapped=matrix@corner+offset
            if not np.isfinite(mapped).all() or np.any(np.rint(mapped)<-_INT_MAX-1) or np.any(np.rint(mapped)>_INT_MAX):
                _fail('grid symmetry indices outside native integer range')
            if np.max(abs(residual@corner+constant))>1e-8:
                _fail('symmetry does not preserve the grid')
        transformed=fractional@rotation.T+translation
        used=set()
        for position,kind in zip(transformed,species):
            delta=position-fractional;delta-=np.rint(delta)
            distance=np.linalg.norm(delta@cell.T,axis=1)
            candidates=np.flatnonzero((species==kind)&(distance<1e-7))
            if len(candidates)!=1 or int(candidates[0]) in used:
                _fail('symmetry does not uniquely preserve atom positions/species')
            used.add(int(candidates[0]))

    return grid
