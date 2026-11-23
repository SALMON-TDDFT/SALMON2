"""Use the same compiled serial C ABI plugin through point-MPI distribution.
Usage: mpiexec -n 3 python tools/test_c_material_mpi.py BUILD_JSON NEW_OUTPUT
"""
import json,struct,sys
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy.c_material import CMaterial
from salmonpy.materials import bind_material
from salmonpy.material_response import ResponseSession
from salmonpy.distributed import ParallelContext
from salmonpy.examples.custom_material import ExampleMaterial,waveform,plan_for
from salmonpy.crystal import respond
import salmonpy

record=json.loads(Path(sys.argv[1]).read_text());world=MPI.COMM_WORLD;plan=plan_for(world.size)
with CMaterial(record['library'],record['sha256'],struct.pack('=dd',1.,.8)) as one, \
     CMaterial(record['library'],record['sha256'],struct.pack('=dd',2.,.8)) as two:
    a=bind_material(one,key='one',identity=one.identity);b=bind_material(two,key='two',identity=two.identity)
    bindings={10:a,20:b,30:a,40:b};history=[]
    with ParallelContext.create(world,plan) as context:
        session=ResponseSession(context,bindings,dt=.02)
        try:
            for step in range(1,9):
                A=np.array([waveform(p).values[step] for p in plan.point_ids]) if world.rank==0 else None
                result=session.advance(A)
                if world.rank==0:history.append(result.copy())
        finally:session.close()
if world.rank==0:
    expected=np.stack([respond(ExampleMaterial(gain=1. if p in (10,30) else 2.),waveform(p)).current[1:]
                       for p in plan.point_ids],axis=1)
    actual=np.asarray(history);np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-14)
    summary={'ranks':world.size,'library_sha256':record['sha256'],
        'maximum_reference_difference':float(np.max(abs(actual-expected))),
        'history':actual.tolist(),'package':salmonpy.__file__,
        'scope':'one serial C plugin state per material point; synthetic interface test'}
    with Path(sys.argv[2]).open('x') as f:json.dump(summary,f,indent=2,allow_nan=False)
    print('C material plugin point-MPI parity PASS',world.size)
