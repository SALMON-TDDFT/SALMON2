"""Identity envelope for the mixed fixture; not a backend/initial-state factory."""
import hashlib
import json
from salmonpy.config import _copy_json, _identifier
from salmonpy.errors import InputError
from salmonpy.point_runs import _identity


def digest_record(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def sample_records(run_id,signature,points,bindings,*,constitutive,cell_ids,assembly=None,constitutive_records=None):
    """Explicit fixture model/input fingerprints, excluding output paths/budgets.

    General hosts supply their own material records through ResponseSession.
    These generated labels are sample defaults, not chemical identity claims.
    """
    _identifier(run_id,'run_id')
    labels=({p['transport_id']:p['point_id'] for p in assembly.record()['points']} if assembly is not None
            else {p:'point-'+str(p) for p in points})
    assets=signature.get('assets',{})
    entries=[]
    for point in points:
        binding=bindings[point]
        if binding.material_record is not None:
            record=_identity(dict(binding.material_record))
        else:
            mid=binding.key
            if mid=='sbe-point':
                definition={'kind':mid,'model':signature['sbe_model'],'contract':'builtin endpoint charge/1'}
                initial={'material':definition,'preparation':'builtin ground occupation/1'}
                initial_id=mid+':ground'
            else:
                point_assets=(signature.get('hse_assets',assets) if mid=='hse-point' else
                    signature.get('standard_material_inputs',{}).get(mid,{}).get('assets',assets))
                point_inputs={n:sha for n,sha in point_assets.items() if not n.startswith('data_for_restart/') and not n.endswith('.json')}
                definition={'kind':mid,'inputs':point_inputs,
                    'engine':signature.get('fast_hse_backend') if mid=='hse-point' else signature['library']}
                initial={'assets':point_assets}
                initial_id=mid+':gs'
            record={'material_id':mid,'definition_digest':digest_record(definition),
                    'initial_state_id':initial_id,'initial_state_digest':digest_record(initial)}
        entries.append({'point_id':labels[point],**record})
    point_record={'schema':'salmonpy.point-run/1','run_id':run_id,'points':entries}
    materials={mid:{'material_id':mid,'definition_digest':digest_record(parameters),
                    'initial_state_id':mid+':zero','initial_state_digest':digest_record({'P':0,'J':0,'E':0})}
               for mid,parameters in constitutive.items()}
    for mid,record in (constitutive_records or {}).items():
        if record is not None:
            checked=_identity(dict(record))
            if checked['material_id']!=mid:
                raise InputError('constitutive material ID differs from configured ID')
            materials[mid]=checked
    return point_record,labels,{'materials':materials,'cell_ids':dict(cell_ids)}


def checkpoint_identity(coupled):
    callback=getattr(coupled.response,'identity_record',None)
    metadata=callback() if callable(callback) else None
    constitutive=getattr(coupled,'constitutive_identity',None)
    if metadata is None and constitutive is None:
        return None
    if metadata is None:
        raise InputError('mixed identity requires identified response points')
    from salmonpy.material_response import _point_metadata
    metadata=_copy_json(metadata)
    if type(metadata) is not dict or set(metadata)!={'point_run','mpi_point_map'}:
        raise InputError('invalid mixed point identity metadata')
    mapping={entry['transport_point_id']:entry['point_id'] for entry in metadata['mpi_point_map']}
    checked=_point_metadata(coupled.response.context.plan,metadata['point_run'],mapping)
    if digest_record(checked)!=digest_record(metadata):
        raise InputError('mixed identity ownership/map differs from current plan')
    constitutive=_copy_json({'materials':{},'cell_ids':{}} if constitutive is None else constitutive)
    if type(constitutive) is not dict or set(constitutive)!={'materials','cell_ids'}:
        raise InputError('invalid constitutive identity')
    if type(constitutive['materials']) is not dict or type(constitutive['cell_ids']) is not dict:
        raise InputError('constitutive identity requires material and cell mappings')
    definitions={p['material_id']:p['definition_digest'] for p in checked['point_run']['points']}
    initials={p['initial_state_id']:(p['material_id'],p['initial_state_digest'])
              for p in checked['point_run']['points'] if 'initial_state_id' in p}
    for mid,record in constitutive['materials'].items():
        _identifier(mid,'constitutive material_id')
        record=_identity(record)
        if record['material_id']!=mid or definitions.setdefault(mid,record['definition_digest'])!=record['definition_digest']:
            raise InputError('conflicting mixed material identity')
        if 'initial_state_id' in record:
            signature=(mid,record['initial_state_digest'])
            if initials.setdefault(record['initial_state_id'],signature)!=signature:
                raise InputError('conflicting mixed initial identity')
    for number,mid in constitutive['cell_ids'].items():
        if not number.isdecimal() or str(int(number))!=number or int(number)<1 or mid not in constitutive['materials']:
            raise InputError('invalid constitutive material-ID mapping')
    result={'schema':'salmonpy.mixed-identity/1','point_metadata':checked,'constitutive':constitutive}
    if len(json.dumps(result).encode())>16*2**20:
        raise InputError('mixed identity metadata exceeds16MiB')
    return result


def check_identity(current,saved,schema):
    if schema not in ('salmonpy.native-direct/4','salmonpy.native-replay/2'):
        raise InputError('unsupported checkpoint schema; use the current writer')
    if current is None and saved is None:return
    if current is None or saved is None or digest_record(current)!=digest_record(saved):
        raise InputError('checkpoint identity mismatch: run/point/material/initial/ownership must match')
