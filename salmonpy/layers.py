"""Explicit native mixed-layer assembly; no MPI initialization or solver startup.

Lengths are Bohr. Every electronic grid cell becomes one independent point;
LD cells stay in the field group. This is geometry/rank planning, not convergence
or memory estimation. Resolved models and user bindings are supplied explicitly.
"""
from dataclasses import dataclass, asdict
from decimal import Decimal, ROUND_CEILING
from itertools import product
import math
import warnings
from collections.abc import Mapping
from types import MappingProxyType
import numpy as np
from .config import _copy_json, _object, _identifier, _freeze, _thaw
from .errors import InputError, InputTypeError
from .geometry import Grid, Geometry, Shape
from .fdtd import FDTD, FDTDMPI, PlaneWave
from .ld import LorentzDrude
from .material_response import MaterialBinding, _supported_group_sizes
from .distributed import ParallelPlan, SolverGroup

__all__ = ['NativeLayerAssembly', 'assemble_native_layers']


def _int(value, name, minimum=1):
    if type(value) is not int or value < minimum:
        raise InputError(f'{name}: integer >= {minimum} required')
    return value


@dataclass(frozen=True)
class NativeLayerAssembly:
    calculation: FDTD
    parallel_plan: ParallelPlan
    points: Mapping
    bindings: Mapping
    field_grid: tuple
    _record: object

    def record(self):
        """Detached JSON-safe layout/identity/rank report; not a checkpoint."""
        return _thaw(self._record)


def assemble_native_layers(layout, materials, *, grid, dt, steps, source,
                           field_grid, ranks_per_point, rank_budget):
    """Build FDTD, point bindings and MPI plan from explicit ordered layers.

    ``materials`` maps user IDs to LD models (including configured LD models)
    or existing MaterialBinding instances. Initial-state loaders and backend
    execution capabilities remain caller-supplied. No prepare/create is called.
    ``origin_cell`` is the first axial cell and the transverse anchor; each
    transverse interval starts at anchor - count//2. Reported cell bounds are
    authoritative, including the half-cell shift for an even transverse count.
    Null material_id means a vacuum gap. No implicit interface gap is inserted.
    """
    data = _copy_json(layout)
    _object(data, 'layout', ('schema','length_unit','axis','origin_cell','layers'))
    if data['schema'] != 'salmonpy.native-layers/1' or data['length_unit'] != 'bohr':
        raise InputError('native-layers/1 and explicit bohr lengths required')
    axis = _int(data['axis'], 'axis', 0)
    if axis > 2: raise InputError('axis must be 0,1,2')
    if not isinstance(grid, Grid) or not isinstance(source, PlaneWave):
        raise InputTypeError('Grid and PlaneWave required')
    if source.axis != axis: raise InputError('layer axis must match source propagation axis')
    if not isinstance(materials, Mapping) or not isinstance(ranks_per_point, Mapping):
        raise InputTypeError('material and ranks_per_point mappings required')
    budget = _int(rank_budget, 'rank_budget')
    origin = data['origin_cell']
    if not isinstance(origin,list) or len(origin)!=3: raise InputError('three origin_cell indices required')
    for i,v in enumerate(origin):
        _int(v,'origin_cell',0)
        if v >= grid.counts[i]: raise InputError('origin_cell outside grid')
    parallel = FDTDMPI(tuple(field_grid))
    nf = parallel.ranks
    if not isinstance(data['layers'],list) or not data['layers']: raise InputError('ordered layers required')
    transverse = [i for i in range(3) if i != axis]
    boxes=[];records=[];names=set();ld_models={};electronic={};configured_records={};needed=nf;cursor=origin[axis]
    from .builtin_materials import ConfiguredMaterial
    for item in data['layers']:
        _object(item,'layer',('id','material_id','thickness','transverse_cells'))
        name=_identifier(item['id'],'layer.id')
        if name in names: raise InputError('duplicate layer id: '+name)
        names.add(name)
        value=item['thickness']
        if type(value) not in (int,float) or not math.isfinite(value) or value<=0:
            raise InputError('finite positive thickness required')
        cells=int((Decimal(str(value))/Decimal(str(grid.spacing[axis]))).to_integral_value(rounding=ROUND_CEILING))
        widths=item['transverse_cells']
        if not isinstance(widths,list) or len(widths)!=2: raise InputError('two transverse cell counts required')
        for width in widths: _int(width,'transverse_cells')
        lo=list(origin);hi=list(origin)
        lo[axis]=cursor;hi[axis]=cursor+cells
        for i,width in zip(transverse,widths):lo[i]=origin[i]-width//2;hi[i]=lo[i]+width
        if any(l<0 or h>n for l,h,n in zip(lo,hi,grid.counts)): raise InputError('layer outside grid: '+name)
        cursor=hi[axis];mid=item['material_id'];kind='vacuum';rank_count=0
        if mid is not None:
            _identifier(mid,'material_id')
            if mid not in materials: raise InputError('unknown material ID: '+mid)
            model=materials[mid]
            if isinstance(model,ConfiguredMaterial):
                if model.material_id != mid: raise InputError('configured material ID differs from mapping key: '+mid)
                configured_records[mid]=model.record()
                model=model.model
            if isinstance(model,LorentzDrude):
                kind='ld';ld_models.setdefault(mid,model)
            elif isinstance(model,MaterialBinding):
                if model.material_record is not None and model.material_record.get('material_id') != mid:
                    raise InputError('binding material record differs from mapping ID: '+mid)
                if model.key != mid: raise InputError('binding key must equal material ID: '+mid)
                kind='electronic';rank_count=_int(ranks_per_point.get(mid),'ranks_per_point['+mid+']')
                mode=getattr(model.backend,'execution',None)
                if mode not in ('leader_only','collective'): raise InputError('backend execution declaration required')
                if mode=='leader_only' and rank_count!=1:
                    raise InputError('leader-only layer backend requires one rank per point')
                _supported_group_sizes(model.backend,rank_count,mid)
                from .embedded_tddft import NativeTDDFTBackend
                if isinstance(model.backend,NativeTDDFTBackend) and math.prod(model.backend.decomposition)!=rank_count:
                    raise InputError('TDDFT inner decomposition differs from ranks_per_point')
                electronic[mid]=model
                needed+=cells*math.prod(widths)*rank_count
            else: raise InputTypeError('material must be LD or an explicit MaterialBinding: '+mid)
        records.append({**item,'kind':kind,'cells':cells,'actual_thickness':cells*grid.spacing[axis],
                        'cell_start':lo,'cell_stop':hi,'ranks_per_point':rank_count})
        boxes.append((lo,hi,mid,kind,name,rank_count))
    if set(ranks_per_point)!=set(electronic): raise InputError('rank settings required only for selected electronic materials; LD is field-local')
    if needed>budget:
        warnings.warn(f'required ranks {needed} exceed rank budget {budget}; continuing planning without subsampling; launch with {needed} ranks', RuntimeWarning, stacklevel=2)
    if not electronic: raise InputError('mixed assembly requires electronic points; use FDTD for LD-only calculations')
    # Validate geometry without allocating a full material-ID volume.
    coords=grid.coordinates;source_index=int(np.argmin(abs(coords[axis]-source.position)))
    shapes=[];ld_ids={mid:i+1 for i,mid in enumerate(ld_models)}
    for lo,hi,mid,kind,name,ranks in boxes:
        if kind=='ld':
            if any(l<8 or h>n-8 for l,h,n in zip(lo,hi,grid.counts)):
                raise InputError('LD layer touches absorbing boundary: '+name)
            shapes.append(Shape('box',tuple((h-l)*d for l,h,d in zip(lo,hi,grid.spacing)),
                                tuple(float((coords[i][l]+coords[i][h-1])/2) for i,(l,h) in enumerate(zip(lo,hi))),ld_ids[mid]))
    points={};bindings={};point_records=[];groups=[SolverGroup('field',tuple(range(nf)),role='field')];rank=nf
    for lo,hi,mid,kind,name,ranks in boxes:
        if kind!='electronic': continue
        for cell in product(*(range(l,h) for l,h in zip(lo,hi))):
            if any(v<10 or v>=n-10 for v,n in zip(cell,grid.counts)):
                raise InputError('electronic layer must stay ten cells from boundaries: '+name)
            if abs(cell[axis]-source_index)<=2: raise InputError('electronic layer intersects source support: '+name)
            if any(k=='ld' and all(l-1<=v<h+1 for v,l,h in zip(cell,lows,highs)) for lows,highs,m,k,n,r in boxes):
                raise InputError('electronic point requires an LD-free constitutive neighborhood: '+name)
            point=len(points)+1;points[point]=cell;bindings[point]=electronic[mid]
            groups.append(SolverGroup('point-'+str(point),tuple(range(rank,rank+ranks)),point_ids=(point,)))
            relative=[v-l for v,l in zip(cell,lo)]
            point_records.append({'transport_id':point,'point_id':name+'/'+','.join(map(str,relative)),
                                  'layer_id':name,'material_id':mid,'cell':list(cell),'ranks':list(range(rank,rank+ranks)),
                                  'provider_identity':electronic[mid].provider.identity,'backend_identity':electronic[mid].backend.identity})
            rank+=ranks
    plan=ParallelPlan(tuple(groups));plan.validate(needed)
    geometry=Geometry(grid,tuple(shapes))
    observations=tuple(tuple(float(coords[i][v]) for i,v in enumerate(cell)) for cell in points.values())
    calculation=FDTD(geometry,tuple(ld_models.values()),dt,steps,source=source,observations=observations)
    parallel.input_text(calculation)  # existing native decomposition/admission checks; no process launch
    material_records = {mid: {'kind': 'ld', 'parameters': asdict(model.config), 'material_record': configured_records.get(mid)} for mid, model in ld_models.items()}
    for mid, binding in electronic.items():
        material_records[mid] = {'kind': 'electronic', 'provider_identity': binding.provider.identity,
                                 'backend_identity': binding.backend.identity,
                                 'material_record': _thaw(binding.material_record)}
    report={'schema':'salmonpy.native-layer-assembly/1','layout':data,'layers':records,'points':point_records,
            'constitutive_material_ids':ld_ids,'required_ranks':needed,'rank_budget':budget,
            'field_grid':list(parallel.grid),'length_unit':'bohr','memory_bytes':None,
            'materials':material_records,'grid':asdict(grid),'dt':calculation.dt,'steps':calculation.steps,'source':asdict(source),
            'note':'one independent electronic state per grid cell; LD stays in field group; no implicit gaps'}
    return NativeLayerAssembly(calculation,plan,MappingProxyType(points),MappingProxyType(bindings),
                               tuple(parallel.grid),_freeze(report))
