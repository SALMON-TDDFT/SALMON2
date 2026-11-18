"""Offline per-material declarations for the ordinary mixed sample; no native loading."""
from dataclasses import dataclass
import json
import math
from pathlib import Path
from salmonpy.errors import InputError

@dataclass(frozen=True)
class MaterialInput:
    bundle: Path
    input_file: str = 'rt.inp'
    decomposition: tuple = (1,1,1,1,1)
    points_per_group: int = 1

    @property
    def ranks(self):
        return math.prod(self.decomposition)


def _unique(pairs):
    result={}
    for key,value in pairs:
        if key in result:raise InputError('duplicate material declaration: '+key)
        result[key]=value
    return result


def load_materials(path):
    """Resolve local references against the catalog directory; providers verify contents."""
    path=Path(path).expanduser().resolve()
    data=json.loads(path.read_text(),object_pairs_hook=_unique)
    if not isinstance(data,dict) or set(data)!={'schema','materials'} or data['schema']!='salmonpy.native-tddft-materials/1':
        raise InputError('native-tddft-materials/1 catalog required')
    if not isinstance(data['materials'],dict):raise InputError('materials must map IDs to declarations')
    result={}
    for key,entry in data['materials'].items():
        if not isinstance(key,str) or not key.strip() or key in ('sbe-point','hse-point','drude','lorentz'):
            raise InputError('invalid or reserved TDDFT material ID: '+str(key))
        if not isinstance(entry,dict) or 'bundle' not in entry or set(entry)-{'bundle','input_file','decomposition','points_per_group'}:
            raise InputError('invalid material declaration: '+key)
        bundle=entry['bundle'];name=entry.get('input_file','rt.inp')
        if not isinstance(bundle,str) or not bundle.strip() or '\x00' in bundle:
            raise InputError('local bundle path required: '+key)
        if not isinstance(name,str) or not name or any(x in name for x in ('\x00','\n','\r')) or Path(name).is_absolute() or '..' in Path(name).parts:
            raise InputError('relative input_file required: '+key)
        decomposition=entry.get('decomposition',[1,1,1,1,1]);width=entry.get('points_per_group',1)
        if not isinstance(decomposition,list) or len(decomposition)!=5 or any(type(n) is not int or not 1<=n<2**31 for n in decomposition):
            raise InputError('five positive decomposition integers required: '+key)
        if type(width) is not int or width<1:raise InputError('positive points_per_group required: '+key)
        source=Path(bundle).expanduser()
        if not source.is_absolute():source=path.parent/source
        result[key]=MaterialInput(source.resolve(),name,tuple(decomposition),width)
    return result
