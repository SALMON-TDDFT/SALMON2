"""Serial controller: MPI FDTD with local, multi-material LD states.
python tools/test_fdtd_mpi_entrypoint.py AUDITED_MPI_BUILD_JSON NEW_OUTPUT
"""
import argparse,json,subprocess
from pathlib import Path
import numpy as np
from salmonpy.fdtd import FDTD,FDTDMPI,SalmonExecutable,PlaneWave
from salmonpy.ld import LorentzDrude
from salmonpy.geometry import Grid,Geometry,Shape
from salmonpy.fields import AcosPulse

p=argparse.ArgumentParser(description=__doc__);p.add_argument('build');p.add_argument('output');a=p.parse_args()
root=Path(a.output).resolve();root.mkdir(parents=True,exist_ok=False)
engine=SalmonExecutable.from_mpi_build(a.build)
models=(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.]),
        LorentzDrude(plasma_frequency=15.,strength=[.7,.3],resonance=[25.,35.],damping=[2.,3.]),
        LorentzDrude(plasma_frequency=12.,strength=[.5,.3,.2],resonance=[0.,20.,30.],damping=[1.,2.,3.]))
grid=Grid((48,40,40),(2.,2.,2.))
def calculation(relabel=False):
    shapes=tuple(Shape('box',(8.,16.,16.),(x,0.,0.),4-i if relabel else i) for i,x in enumerate((-8.,0.,8.),1))
    return FDTD(Geometry(grid,shapes),models[::-1] if relabel else models,.003,160,
        source=PlaneWave(AcosPulse(.001,20.,.15),0,-24.),observations=((-16.,0.,0.),(0.,0.,0.),(16.,0.,0.)))
base=engine.run(calculation(),root/'one',threads=1)
report={}
for name,mpi,relabel in [('four_xy',FDTDMPI((2,2,1)),False),('four_yz',FDTDMPI((1,2,2)),False),('relabel',FDTDMPI((2,2,1)),True)]:
    result=engine.run(calculation(relabel),root/name,mpi=mpi,threads=1)
    for x,y in zip(base.observations,result.observations): np.testing.assert_array_equal(x.values,y.values)
    manifest=json.loads((result.directory/'salmonpy-run.json').read_text())
    assert manifest['mpi_ranks']==4 and manifest['mpi']['grid']==list(mpi.grid)
    report[name]={'bitwise_equal':True,'ranks':4,'grid':mpi.grid}
    print(name,'PASS',flush=True)
try:
    engine.run(calculation(),root/'timeout',mpi=FDTDMPI((2,2,1)),timeout=.01,threads=1)
except subprocess.TimeoutExpired: pass
else: raise AssertionError('timeout was not enforced')
assert json.loads((root/'timeout/salmonpy-run.json').read_text())['status']=='failed'
print('timeout cleanup/record PASS',flush=True)
(root/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
