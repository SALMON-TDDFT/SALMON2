"""Exact finite-grid box factors and bounded process-local Fourier reuse."""
from collections import OrderedDict
import numpy as np


class BoundedCache:
    """LRU bounded by array/key payload bytes and number of entries.

    No files or mappings are owned. Python object overhead is additionally
    bounded by max_entries; callers retain immutable arrays after eviction.
    """
    def __init__(self, max_bytes=64*1024**2, max_entries=512):
        self.max_bytes=max_bytes;self.max_entries=max_entries
        self.items=OrderedDict();self.bytes=0;self.hits=0;self.misses=0

    def get(self,key):
        entry=self.items.pop(key,None)
        if entry is None:
            self.misses+=1
            return None
        self.items[key]=entry;self.hits+=1
        return entry[0]

    def put(self,key,value,nbytes):
        old=self.items.pop(key,None)
        if old is not None:self.bytes-=old[1]
        if nbytes>self.max_bytes or self.max_entries<1:return
        while self.items and (self.bytes+nbytes>self.max_bytes or len(self.items)>=self.max_entries):
            _,(_,size)=self.items.popitem(last=False);self.bytes-=size
        self.items[key]=(value,nbytes);self.bytes+=nbytes

    def info(self):
        return dict(entries=len(self.items),bytes=self.bytes,max_bytes=self.max_bytes,
                    max_entries=self.max_entries,hits=self.hits,misses=self.misses)


_phase_cache=BoundedCache()


def clear_phase_cache(max_bytes=64*1024**2):
    global _phase_cache
    _phase_cache=BoundedCache(max_bytes=max_bytes)


def phase_cache_info():
    return _phase_cache.info()


def integer_phases(n,indices):
    """exp(+2 pi i grid·G), x-fastest points, cached modulo integer grid size."""
    n=np.asarray(n);indices=np.asarray(indices)
    if (n.shape!=(3,) or np.any(n<1) or not np.array_equal(n,n.astype(int)) or
            indices.ndim!=2 or indices.shape[1]!=3 or not np.isfinite(indices).all() or
            not np.array_equal(indices,indices.astype(int))):
        raise ValueError('positive grid and integer Fourier indices required')
    n=n.astype(int)
    canonical=np.asarray(indices % n,dtype=np.int64)
    key=(tuple(n),canonical.tobytes())
    cached=_phase_cache.get(key)
    if cached is not None:return cached
    grid=np.indices(tuple(n[::-1])).reshape(3,-1)[::-1].T/n
    phase=np.exp(2j*np.pi*grid@canonical.T)
    phase.flags.writeable=False
    _phase_cache.put(key,phase,phase.nbytes+canonical.nbytes)
    return phase


def box_forms(n,boxes,Q,indices=None):
    """Mean exp(-2 pi i r·Q) on each discrete box, x-fastest box IDs.

    The finite-grid sum factorizes exactly into three one-dimensional sums;
    no continuum box/sinc approximation is introduced. Optional indices are
    zero-based box IDs in any order (duplicates are permitted).
    """
    n=np.asarray(n);Q=np.asarray(Q)
    if (n.shape!=(3,) or np.any(n<1) or not np.array_equal(n,n.astype(int)) or
            not isinstance(boxes,(int,np.integer)) or boxes<1 or np.any(n%boxes) or
            Q.ndim!=2 or Q.shape[1]!=3 or not np.isfinite(Q).all()):
        raise ValueError('finite Q and a box count dividing the grid are required')
    ids=np.arange(boxes**3) if indices is None else np.asarray(indices)
    if ids.ndim!=1 or not np.array_equal(ids,ids.astype(int)) or np.any(ids<0) or np.any(ids>=boxes**3):
        raise ValueError('invalid box indices')
    ids=ids.astype(int)
    result=np.ones((len(ids),len(Q)),complex)
    for axis in range(3):
        points=np.arange(n[axis]).reshape(boxes,-1)/n[axis]
        factors=np.exp(-2j*np.pi*points[:,:,None]*Q[None,None,:,axis]).mean(axis=1)
        result*=factors[(ids//boxes**axis)%boxes]
    return result
