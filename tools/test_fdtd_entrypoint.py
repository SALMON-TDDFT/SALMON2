"""Integration check against an already audited SALMON executable.

Usage: python tools/test_fdtd_entrypoint.py EXECUTABLE SHA256 NEW_OUTPUT_DIRECTORY
Run from the source checkout with PYTHONPATH=. . Never runs a GS calculation.
"""
from pathlib import Path
import os
import subprocess
import sys
import numpy as np
from salmonpy.fdtd import FDTD, Dielectric, PlaneWave, SalmonExecutable
from salmonpy.geometry import Grid, Shape, Geometry
from salmonpy.fields import AcosPulse

exe,sha,out=sys.argv[1:]
root=Path(out).resolve();root.mkdir(parents=True,exist_ok=False)
runner=SalmonExecutable(Path(exe),sha)
for kind,counts,size in [('ellipsoid',(32,32,32),(8,8,8)),('box',(33,32,34),(4,6,8))]:
    geometry=Geometry(Grid(counts,(2,2,2)),(Shape(kind,size),))
    calc=FDTD(geometry,(Dielectric(2.),),.001,200,
              source=PlaneWave(AcosPulse(.001,10.,.05),0,-12.))
    result=runner.run(calc,root/(kind+'-wrapped'))
    lines=(result.directory/'shape.cube').read_text().splitlines()
    natom=int(lines[2].split()[0])
    shape=np.array([float(v) for line in lines[6+natom:] for v in line.split()]).reshape(counts)
    np.testing.assert_array_equal(shape,geometry.material_ids())
    # Independent handwritten input for the same physical problem. Defaults for
    # vacuum, boundaries, pulse start and unrotated shape are deliberately used.
    direct=root/(kind+'-direct');direct.mkdir()
    name='ellipsoid' if kind=='ellipsoid' else 'rectangular-cylinder'
    inp=f'''&calculation
 theory='maxwell'
/
&control
 sysname='fdtd'
/
&units
 unit_system='au'
/
&system
 yn_periodic='n'
/
&emfield
 ae_shape1='Acos2'
 I_wcm2_1=-1
 E_amplitude1=0.001
 omega1=10
 tw1=0.05
 phi_cep1=0.25
 epdir_re1=0,0,1
/
&analysis
 nenergy=100
 de=0.001
/
&maxwell
 al_em={','.join(str(2*n) for n in counts)}
 dl_em=2,2,2
 dt_em=0.001
 nt_em=200
 yn_make_shape='y'
 yn_output_shape='y'
 n_s=1
 typ_s(1)='{name}'
 id_s(1)=1
 inf_s(1,1:3)={','.join(map(str,size))}
 media_num=1
 epsilon_em(1)=2
 wave_input='source'
 ek_dir1=1,0,0
 source_loc1=-12,0,0
 obs_num_em=1
 obs_samp_em=1
 obs_loc_em(1,:)=0,0,0
/
'''
    (direct/'inputfile').write_text(inp)
    with (direct/'stdout.log').open('w') as log:
        subprocess.run([str(runner.path)],input=inp,text=True,cwd=direct,stdout=log,stderr=subprocess.STDOUT,
                       timeout=300,check=True,env={**os.environ,'OMP_NUM_THREADS':'1'})
    assert 'end SALMON' in (direct/'stdout.log').read_text()
    for file in ['obs1_at_point_rt.data','shape.cube']:
        assert (direct/file).read_bytes()==(result.directory/file).read_bytes(),file
    print(kind, 'PASS: 200-step observations and geometry match handwritten SALMON input byte-for-byte')
