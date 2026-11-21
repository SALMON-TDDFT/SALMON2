"""Run the ordinary mixed sample and record owner-only native current routing.

Use its usual arguments; the runtime must expose the spatial-current ABI.
This is an explicitly invoked local contract test, not a CI physics campaign.
"""
import ctypes as c
import json
from pathlib import Path
import sys
from unittest.mock import patch
import numpy as np
from mpi4py import MPI

def main():
    from salmonpy.coupled import NativeEHField
    from salmonpy.examples import mpi_mixed_native as sample
    record=json.loads(Path(sys.argv[1]).read_text())
    lib=c.CDLL(record['library'])
    if not hasattr(lib,'sp_eh_current_mask') or not hasattr(lib,'sp_eh_step_owned'):
        raise RuntimeError('runtime lacks spatial-current routing ABI')
    original_init=NativeEHField.__init__
    original_close=NativeEHField.close
    original_advance=NativeEHField._advance
    records=[]
    def initialize(self,*args,**kwargs):
        original_init(self,*args,**kwargs)
        rows=self.current_rows
        assert rows is not None
        records.append({'rank':MPI.COMM_WORLD.rank,'rows':list(map(int,rows)),
                        'points':len(self.point_ids),'all_rows':[list(map(int,r)) for r in self.current_rows_by_rank]})
    def advance(self,J_previous,*,owned):
        if owned:
            if not self.context.is_leader:
                other=np.ones(len(self.point_ids),dtype=bool);other[self.current_rows]=False
                assert not np.any(J_previous[other])
            if self.step_index==1:
                fn=self.lib.sp_eh_step_owned
                floats=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
                fn.argtypes=[c.c_int,c.c_int,floats,floats,floats];fn.restype=c.c_int
                target=next(i for i,rows in enumerate(self.current_rows_by_rank) if len(rows))
                before=(self.A.tobytes(),self.E.tobytes(),self.step_index)
                for bad in ('nan','inf','count'):
                    local=np.ascontiguousarray(J_previous[self.current_rows]);count=len(local)
                    if self.context.solver_comm.rank==target:
                        if bad=='count':count+=1
                        else:local[0,0]=float(bad)
                    E=np.empty_like(self.E)
                    status=self._call(fn,len(self.point_ids),count,local,self.A,E)
                    assert status==10,(bad,status)
                    assert before==(self.A.tobytes(),self.E.tobytes(),self.step_index)
                records[-1]['invalid_inputs_rejected_without_step']=['nan','inf','count']
        return original_advance(self,J_previous,owned=owned)
    def close(self):
        if not self.closed:
            assert self.manifest.get('current_routing')=='spatial-owned'
            (self.directory/f'spatial-current-rank-{self.context.solver_comm.rank}.json').write_text(json.dumps(records,indent=2)+'\n')
        return original_close(self)
    with patch.object(NativeEHField,'__init__',initialize),patch.object(NativeEHField,'close',close),patch.object(NativeEHField,'_advance',advance):
        sample.main()

if __name__=='__main__':main()
