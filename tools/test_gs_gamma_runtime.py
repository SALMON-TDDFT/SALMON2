"""Synthetic-orbital full native TDDFT bridge check; NOT a converged GS test.

python tools/test_gs_gamma_runtime.py CONVERTER_BUILD NEW_OUTPUT
Only existing paired runtime is reused by recorded SHA; never builds a solver.
"""
import argparse,hashlib,json,os,shutil
from pathlib import Path
import numpy as np
from salmonpy.gs_archive import save_ground_state
from salmonpy.gs_native import GammaRestartConverter
from salmonpy.fdtd import SalmonExecutable
from salmonpy.dft import Pseudopotential
from salmonpy.tddft import TDDFTInitialState


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('converter_build',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    record=json.loads(a.converter_build.read_text());converter=GammaRestartConverter.from_build(a.converter_build)
    engine=SalmonExecutable(record['target_executable'],record['target_sha256'])
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    package=Path(__import__('salmonpy').__file__).resolve().parent
    potential=Pseudopotential.from_file(package/'examples/si_intro_data/Si_rps.dat',atomic_number=14,local_channel=2)
    # Two synthetic orthonormal real plane waves in a periodic cell, 4electrons.
    nx=ny=nz=8;volume=8.**3;z,y,x=np.indices((nz,ny,nx))
    waves=np.empty((1,2,nz,ny,nx),complex)
    waves[0,0]=1/np.sqrt(volume);waves[0,1]=np.sqrt(2/volume)*np.cos(2*np.pi*x/nx)
    density=2*(abs(waves[0,0])**2+abs(waves[0,1])**2)
    archive=save_ground_state(root/'synthetic-archive',wavefunctions=waves,density=density,eigenvalues=[[-.5,.2]],
        occupations=[[2.,2.]],kpoints=[[0.,0.,0.]],weights=[1.],cell_bohr=np.eye(3)*8,
        fractional_positions=[[.3,.4,.45]],atomic_numbers=[14],electrons=4.,energy_hartree=-1.,xc='pz',
        pseudopotentials=[dict(atomic_number=14,filename='Si_rps.dat',sha256=potential.sha256,source='bundled licensed potential')],orbital_block=1)
    initial=archive.tddft_initial_state(converter=converter,executable=engine,pseudopotentials=[potential],max_steps=12,timeout=60.)
    histories=[]
    try:
        with initial.new_state(dt=.001) as one,initial.new_state(dt=.001) as two:
            values=np.array([[0,0,.0001*n] for n in range(1,7)])
            for A in values[:2]:one.step(A)
            saved=one.snapshot();two.restore(saved)
            for A in values[2:]:
                j1=one.step(A);j2=two.step(A);assert j1.tobytes()==j2.tobytes();histories.append(j1)
            initial.close();one.restore(saved)
            replay=np.array([one.step(A) for A in values[2:]])
            expected=np.array(histories);assert replay.tobytes()==expected.tobytes()
            result={'scope':'synthetic orthonormal orbitals + real paired native TDDFT; no SCF/physical convergence claim',
                    'bitwise_independent_and_restored':True,'max_current_difference':float(np.max(abs(replay-expected))),
                    'target_sha256':engine.sha256,'converter_sha256':converter.sha256,'time':one.time,
                    'currents':replay.tolist(),'norms':(np.sum(abs(waves)**2,axis=(2,3,4))*volume/(nx*ny*nz)).tolist()}
            (root/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
            # Collect logs/protocol metrics only, never physical GS payloads.
            for label,state in [('one',one),('two',two)]:
                shutil.copy2(Path(state._work.name)/'run/worker.log',root/(label+'-worker.log'))
    finally:initial.close()

if __name__=='__main__':main()
