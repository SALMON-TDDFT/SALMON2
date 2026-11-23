"""Every accepted-A and field-A component rejects NaN/Inf before staged publication."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.embedded_tddft import NativeTDDFTBackend,NativeTDDFTMaterialProvider


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','fixture','output'):p.add_argument(name,type=Path)
    a=p.parse_args();comm=MPI.COMM_WORLD
    root=a.output.resolve()
    if comm.rank==0:root.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    context=SimpleNamespace(solver_comm=comm,point_ids=(1,2),closed=False,_dependents=0)
    initial=NativeTDDFTMaterialProvider(a.fixture,'rt.inp',max_steps=16).prepare(context)
    backend=NativeTDDFTBackend(a.build,root,decomposition=(1,comm.size,1,1,1))
    model=backend.open_model(initial,.003,context);source=model.new_state(1);target=model.new_state(2)
    cases=0
    try:
        for i in range(1,9):source.step([1e-4*i,2e-4*i,-1e-4*i])
        good=root/f'good-{comm.rank}.bin';source.checkpoint(good)
        original=good.read_bytes()
        assert original[-32:].rstrip()==b'SALMONPY_TDDFT_COMPLETE_V2'
        for offset in (len(original)-256,len(original)-64):
            for axis in range(3):
                for value in (np.nan,np.inf):
                    data=bytearray(original)
                    if comm.rank==comm.size-1:
                        data[offset+8*axis:offset+8*(axis+1)]=np.asarray(value,dtype=np.float64).tobytes()
                    bad=root/f'bad-{cases}-{comm.rank}.bin';bad.write_bytes(data)
                    code=target._call(target.lib.sp_tddft_checkpoint,str(bad).encode(),2)
                    assert comm.allgather(int(code))==[5]*comm.size
                    target._refresh();assert target.steps==0 and not target.failed
                    assert source.steps==8
                    cases+=1
        target.checkpoint(good,restore=True)
        np.testing.assert_array_equal(source.step([9e-4,18e-4,-9e-4]),target.step([9e-4,18e-4,-9e-4]))
    finally:
        target.close();source.close();model.close()
    assert context._dependents==0
    if comm.rank==0:(root/'result.json').write_text(json.dumps({'status':'pass','cases':cases,'ranks':comm.size,'all_components_rejected':True,'valid_resume_exact':True})+'\n')


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback
        traceback.print_exc();MPI.COMM_WORLD.Abort(1)
