"""Explicit FDTD transverse-periodic slab regression (no GS data).

AUDITED_MPI_EXECUTABLE_BUILD_JSON NEW_RESULTS --checker CHECKER
FDTD/LD updates remain native; Fourier sums and the analytic reference use C.
Run from fresh serial Python. Fixture parameters test software, not real optics.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
from salmonpy.fdtd import FDTD, Boundaries, Dielectric, PlaneWave, FDTDMPI, SalmonExecutable
from salmonpy.geometry import Grid, Geometry, Shape
from salmonpy.fields import AcosPulse
from salmonpy.ld import LorentzDrude


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--checker',type=Path,default=Path(__file__).with_name('check_forbidden_eigensolvers.py'))
    parser.add_argument('--cc',default='cc')
    parser.add_argument('--epsilon-infinity',type=float,default=1.0)
    args=parser.parse_args();root=Path(__file__).resolve().parents[1]
    engine=SalmonExecutable.from_mpi_build(args.build)
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=root/'tests/native/test_fdtd_fresnel.c'
    helper=out/'helper';helper.mkdir();obj=helper/'reference.o';exe=helper/'reference'
    subprocess.run([sys.executable,str(args.checker),'--source-root',str(source.parent)],check=True)
    subprocess.run([args.cc,'-std=c99','-c',str(source),'-o',str(obj)],check=True)
    subprocess.run([sys.executable,str(args.checker),'--object-tree',str(helper)],check=True)
    subprocess.run([args.cc,str(obj),'-lm','-o',str(exe)],check=True)
    helper_sha=hashlib.sha256(exe.read_bytes()).hexdigest()
    grid=Grid((256,8,8),(.5,.5,.5));dt=.0015;steps=4000
    boundary=Boundaries(('pml','pml','periodic','periodic','periodic','periodic'))
    source_wave=PlaneWave(AcosPulse(.001,10.,2.),0,-40.)
    observations=((30.,-.25,-.25),(30.,.25,.25))
    traces={};data={};reports=[]
    for ranks in (1,2):
        mpi=FDTDMPI((1,ranks,1))
        for name,shapes,materials in [('vacuum',(),()),('slab',(Shape('box',(20.,4.,4.)),),(Dielectric(2.25),)),
                                      ('ld-slab',(Shape('box',(20.,4.,4.)),),(LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[25.],damping=[2.],epsilon_infinity=args.epsilon_infinity),))]:
            calc=FDTD(Geometry(grid,shapes),materials,dt,steps,boundaries=boundary,
                source=source_wave,observations=observations)
            result=engine.run(calc,out/f'{name}-{ranks}',mpi=mpi,threads=1,timeout=180)
            tables=result.observations
            np.testing.assert_allclose(tables[0].values[:,1:],tables[1].values[:,1:],rtol=1e-10,atol=1e-12)
            traces[(name,ranks)]=tables[0].values
            data[(name,ranks)]=result.directory
            if ranks==2:np.testing.assert_allclose(traces[(name,1)],traces[(name,2)],rtol=1e-10,atol=1e-12)
        # Output column names are part of the SALMON table, not assumed offsets.
        table=tables[0];columns=table.columns
        eindex=next(i for i,label in enumerate(columns) if label.lower() in ('ez','e_z'))
        time=traces[('vacuum',ranks)][:,0]
        for label in ('slab','ld-slab'):
            samples=out/f'fourier-input-{label}-{ranks}.txt'
            np.savetxt(samples,np.column_stack((time,traces[('vacuum',ranks)][:,eindex],traces[(label,ranks)][:,eindex])),fmt='%.17g')
            if hashlib.sha256(exe.read_bytes()).hexdigest()!=helper_sha:raise ValueError('validated C helper changed')
            spectrum=out/f'fresnel-{label}-{ranks}.txt'
            run=subprocess.run([str(exe),str(samples),str(spectrum),*(['ld',str(args.epsilon_infinity)] if label=='ld-slab' else [])],capture_output=True,text=True)
            (out/f'fresnel-{label}-{ranks}.log').write_text(run.stdout+run.stderr)
            run.check_returncode()
            reports.append({'ranks':ranks,'material':label,'maximum_relative_complex_transmission_error':float(np.loadtxt(spectrum)[:,-1].max())})
    (out/'result.json').write_text(json.dumps({'status':'passed','cases':reports,
        'transverse_uniformity':True,'MPI1_2_match':True,'analytic_reference':'vacuum/slab/vacuum, d=20au; dielectric n=1.5 and LD wp=20,w0=25,gamma=2',
        'epsilon_infinity':args.epsilon_infinity,'helper_sha256':helper_sha,'executable_sha256':engine.sha256,'GS_data_collected':False},indent=2)+'\n')
    print('transverse-periodic Fresnel slab PASS',reports)


if __name__=='__main__':main()
