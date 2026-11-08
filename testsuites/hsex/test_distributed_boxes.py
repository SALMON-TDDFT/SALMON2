"""Distributed box interpolation, FFT, reciprocal validation and binary layout."""
from pathlib import Path
import sys
import tempfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/hsex'))
from distributed_boxes import q_blocks_to_pairs, interpolate_transform, stream_write
from q_mesh import signed_indices,interpolate_boxes
from prepare_w import transform_boxes_inplace


def run(comm=None):
    rank,size=(0,1) if comm is None else (comm.rank,comm.size)
    for na in (1,3):
        qmesh=np.array([2,2,2]);mesh=np.array([4,4,4])
        matrix=np.eye(na)+.2*np.ones((na,na))
        anti=np.triu(np.ones((na,na)),1);anti-=anti.T
        def tail(q):
            return .1*(1+np.cos(2*np.pi*q).sum())*matrix+.02j*np.sin(2*np.pi*q[0])*anti
        coarse=np.array([(2+.2*np.cos(2*np.pi*q).sum())*matrix+tail(q)
                         for q in signed_indices(qmesh)/qmesh])
        local={i:coarse[i] for i in range(len(coarse)) if i%size==rank}
        first,last,owned=q_blocks_to_pairs(local,qmesh,na,comm,blockbytes=128)
        np.testing.assert_allclose(owned,coarse.reshape(len(coarse),-1)[:,first:last])
        callback_calls=[]
        def callback(q,a,b):
            callback_calls.append((a,b))
            return tail(q).reshape(-1)[a:b]
        origin=np.eye(na)*3.3
        actual=interpolate_transform(owned,qmesh,mesh,callback,first,last,
                                     blockbytes=256,origin=origin.reshape(-1)[first:last])
        assert len(callback_calls)==(int(qmesh.prod())+int(mesh.prod()))*((last-first+3)//4)
        expected=interpolate_boxes(coarse.copy(),qmesh,mesh,tail)
        expected[0]=origin
        transform_boxes_inplace(expected,mesh)
        np.testing.assert_allclose(actual,expected.reshape(len(expected),-1)[:,first:last],atol=1e-13,rtol=0)
        # Equal meshes must skip the long-range interpolation, as in prepare_w.
        equal=interpolate_transform(owned,qmesh,qmesh,callback,first,last,blockbytes=256)
        equal_expected=coarse.copy();transform_boxes_inplace(equal_expected,qmesh)
        np.testing.assert_allclose(equal,equal_expected.reshape(len(coarse),-1)[:,first:last],atol=1e-13)
        root_tmp=tempfile.TemporaryDirectory(prefix='hsex-boxes-') if rank==0 else None
        directory=root_tmp.name if rank==0 else None
        if comm is not None:directory=comm.bcast(directory,root=0)
        output=Path(directory)/'W.bin'
        stream_write(output,b'HEADER',b'TAIL',actual,first,last,na,comm,mesh=mesh)
        if rank==0:
            raw=output.read_bytes();assert raw[:6]==b'HEADER' and raw[-4:]==b'TAIL'
            file_values=np.frombuffer(raw[6:-4],dtype='=f8').reshape(len(expected),na,na)
            np.testing.assert_allclose(file_values,expected.real.transpose(0,2,1),atol=1e-13)
        def rejects(call):
            try:call()
            except (ValueError,FileExistsError):return
            raise AssertionError('expected collective rejection')
        rejects(lambda:stream_write(output,b'',b'',actual,first,last,na,comm,mesh=mesh))
        bad=np.array(actual,copy=True)
        if last>first:bad[0,0]+=1j
        rejects(lambda:stream_write(Path(directory)/'bad.bin',b'',b'',bad,first,last,na,comm,mesh=mesh))
        if rank==0:assert not (Path(directory)/'bad.bin').exists()
        if na>1:
            bad=np.array(actual,copy=True)
            if first<=1<last:bad[0,1-first]+=.1
            rejects(lambda:stream_write(Path(directory)/'asymmetric.bin',b'',b'',bad,
                                        first,last,na,comm,mesh=mesh))
            if rank==0:assert not (Path(directory)/'asymmetric.bin').exists()

        if comm is not None:comm.Barrier()
        if rank==0:root_tmp.cleanup()
    # Missing and duplicate ownership must reject collectively before Alltoallv.
    missing={i:np.ones((1,1)) for i in range(1,8) if i%size==rank}
    try:q_blocks_to_pairs(missing,[2,2,2],1,comm)
    except ValueError:pass
    else:raise AssertionError('missing q accepted')
    if size>1:
        duplicate={i:np.ones((1,1)) for i in range(8) if i%size==rank}
        if rank==1:duplicate[0]=np.ones((1,1))
        try:q_blocks_to_pairs(duplicate,[2,2,2],1,comm)
        except ValueError:pass
        else:raise AssertionError('duplicate q accepted')
    invalid={i:np.ones((1,1)) for i in range(8) if i%size==rank}
    if rank==0:invalid[0][0,0]=np.nan
    try:q_blocks_to_pairs(invalid,[2,2,2],1,comm)
    except ValueError:pass
    else:raise AssertionError('nonfinite input accepted')
    if rank==0:print(f'distributed boxes PASS: {size} MPI ranks')

if __name__=='__main__':
    if '--mpi' in sys.argv:
        from mpi4py import MPI
        run(MPI.COMM_WORLD)
    else:run()
