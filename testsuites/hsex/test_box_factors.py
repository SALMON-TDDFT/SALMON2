"""Separable finite-grid factors and bounded phase/snapshot reuse."""
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import os,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/hsex'))
from box_factors import box_forms, integer_phases, clear_phase_cache, phase_cache_info
from pump_screen import snapshot, clear_snapshot_cache, snapshot_cache_info


def run():
    rng=np.random.default_rng(11203)
    n=np.array([6,12,6]);boxes=3
    Q=np.vstack([rng.normal(size=(21,3)),np.zeros((1,3)),[[6,0,0]]])
    grid=np.indices(tuple(n[::-1])).reshape(3,-1)[::-1].T/n
    label=np.ravel_multi_index(tuple((np.floor(grid*boxes+1e-10).astype(int)).T),(boxes,)*3,order='F')
    expected=np.array([np.exp(-2j*np.pi*grid[label==b]@Q.T).mean(axis=0) for b in range(boxes**3)])
    np.testing.assert_allclose(box_forms(n,boxes,Q),expected,atol=2e-15,rtol=2e-14)
    np.testing.assert_allclose(box_forms(n,boxes,Q,[4,0,26]),expected[[4,0,26]],atol=2e-15)
    assert box_forms(n,boxes,Q,[]).shape==(0,len(Q))
    indices=rng.integers(-5,6,(20,3))
    clear_phase_cache()
    phase=integer_phases(n,indices)
    np.testing.assert_allclose(phase,np.exp(2j*np.pi*grid@indices.T),atol=1e-14)
    assert integer_phases(n,indices+n) is phase
    assert not phase.flags.writeable and phase_cache_info()['hits']==1
    # Force eviction without allocating a large cache during routine tests.
    clear_phase_cache(max_bytes=2*int(n.prod())*16+1024)
    for count in range(5):integer_phases(n,np.array([[count,0,0]]))
    assert phase_cache_info()['bytes'] <= phase_cache_info()['max_bytes']
    assert phase_cache_info()['entries']==2
    before=phase_cache_info()['misses'];integer_phases(n,np.array([[0,0,0]]))
    assert phase_cache_info()['misses']==before+1
    clear_phase_cache()
    with TemporaryDirectory() as tmp:
        n=np.array([2,2,2]);k=np.zeros((1,3));cell=np.eye(3);path=Path(tmp)/'hsex-pump-000001.bin'
        def write(pop):
            with path.open('wb') as f:
                f.write(b'HSP1');np.array([*n,1,1,1,2],dtype='i4').tofile(f)
                for a in [cell,k[0],np.zeros(3),[.1],[pop]]:np.asarray(a,dtype='f8').tofile(f)
                np.ones(8,dtype='c16').tofile(f);np.eye(8,dtype='c16').tofile(f)
        write(2.)
        clear_snapshot_cache()
        original_read=np.fromfile;reads=[]
        def track_read(stream,dtype=float,count=-1,*args,**kwargs):
            reads.append((np.dtype(dtype),count))
            return original_read(stream,dtype,count,*args,**kwargs)
        with patch('pump_screen.np.fromfile',side_effect=track_read):
            metadata=snapshot(tmp,n,1,k,cell,0,False,read_states=False)
        assert metadata[4] is None
        assert all(dtype.kind!='c' for dtype,count in reads), 'metadata pass read orbital payload'
        assert sum(dtype.itemsize*count for dtype,count in reads)==164
        full=snapshot(tmp,n,1,k,cell,0,False)
        assert full[4].shape==(1,8) and np.all(full[4]==1)
        assert snapshot_cache_info()['misses']==2, 'metadata entry must not satisfy a state request'
        clear_snapshot_cache()
        original=Path.open;opens=[]
        def track(p,*a,**kw):opens.append(p);return original(p,*a,**kw)
        with patch.object(Path,'open',track):
            first=snapshot(tmp,n,1,k,cell,0,False)
            second=snapshot(tmp,n,1,k,cell,0,False)
        assert len(opens)==1 and first[4] is second[4]
        assert not first[4].flags.writeable and snapshot_cache_info()['hits']==1
        assert not isinstance(first[4],np.memmap)
        try:snapshot(tmp,n,1,k,2*cell,0,False)
        except AssertionError:pass
        else:raise AssertionError('geometry must be revalidated')
        stamp=path.stat().st_mtime_ns
        write(1.5);os.utime(path,ns=(stamp+1000000,stamp+1000000))
        third=snapshot(tmp,n,1,k,cell,0,False)
        assert third[3][0]==1.5 and first[3][0]==2.
        entry_bytes=snapshot_cache_info()['bytes']//2
        clear_snapshot_cache(max_bytes=entry_bytes)
        snapshot(tmp,n,1,k,cell,0,False)
        write(1.4);os.utime(path,ns=(stamp+2000000,stamp+2000000))
        snapshot(tmp,n,1,k,cell,0,False)
        assert snapshot_cache_info()['entries']==1
        assert snapshot_cache_info()['bytes'] <= entry_bytes
        clear_snapshot_cache(max_bytes=32)
        snapshot(tmp,n,1,k,cell,0,False)
        assert snapshot_cache_info()['entries']==0 and snapshot_cache_info()['bytes']==0
        clear_snapshot_cache()
    print('separable forms / phase reuse / bounded snapshot cache PASS')

if __name__=='__main__':run()
