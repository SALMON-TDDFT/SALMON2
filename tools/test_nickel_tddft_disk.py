"""Explicit Ni SOC TDDFT restart across fresh MPI processes; retains GS/RT data.

Three sequential hosts: SCF4 + continuous RT8; reload GS + RT4 + save;
reload GS + restore RT4 + continue RT5..8. No SCF or RT replay on resume.
The deliberately unconverged GS is a software probe, not a physical result.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def verify_observation_segments(root):
    """Read actual native files; never split/fabricate observations or run SCF."""
    import hashlib
    import numpy as np
    from salmonpy.analysis import read_table
    from salmonpy.outputs import read_spin_dynamics
    root=Path(root)
    record={}
    for suffix in ('rt.data','rt_energy.data','rt_spin.data'):
        paths=[]
        for phase in ('save','resume','continuous'):
            found=list((root/phase).rglob('*_'+suffix))
            if len(found)!=1:raise AssertionError(f'exactly one {phase} {suffix} required')
            paths.append(found[0])
        joined=read_table(paths[:2]);control=read_table(paths[2])
        if (joined.columns,joined.units)!=(control.columns,control.units):
            raise AssertionError('observation definitions differ')
        np.testing.assert_array_equal(joined.values,control.values)
        if suffix=='rt_spin.data':
            spin=read_spin_dynamics(paths[:2]);reference=read_spin_dynamics(paths[2])
            for name in ('time','magnetization','pauli_current'):
                np.testing.assert_array_equal(getattr(spin,name),getattr(reference,name))
        record[suffix]={'exact':True,'rows':len(joined.values),
            'times':joined.column('Time').tolist(),
            'files':[{'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths]}
    return record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--potential',type=Path,required=True)
    p.add_argument('--size',type=int,choices=(1,2),default=1)
    p.add_argument('--worker',choices=('continuous','save','resume'))
    p.add_argument('--verify-outputs-only',action='store_true',help='read an existing completed campaign; run no calculation and write no files')
    args=p.parse_args()
    if args.verify_outputs_only:
        print(json.dumps(verify_observation_segments(args.output),indent=2));return
    if args.worker is None:
        args.output.mkdir(parents=True,exist_ok=False)
        for phase in ('continuous','save','resume'):
            cmd=['mpiexec','--bind-to','none','-n',str(args.size),sys.executable,
                 str(Path(__file__).resolve()),str(args.build.resolve()),str(args.output.resolve()),
                 '--potential',str(args.potential.resolve()),'--size',str(args.size),'--worker',phase]
            with (args.output/(phase+'.log')).open('w') as log:
                subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT,
                               env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
        phases={phase:json.loads((args.output/phase/'result.json').read_text())
                for phase in ('continuous','save','resume')}
        result={'status':'passed','MPI_ranks':args.size,'fresh_MPI_processes':3,
                'electrons':18,'states':36,'temperature_k':300,'kmesh':[2,2,2],
                'phases':phases,'currents_and_electron_numbers_exact':True,
                'observation_segments':verify_observation_segments(args.output),
                'scope':'same binary and topology; short unconverged-GS software probe, not physical convergence'}
        (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
        print('Ni SOC fresh-process TDDFT restart PASS',args.output/'result.json')
        return
    import numpy as np
    from mpi4py import MPI
    from types import SimpleNamespace
    from dataclasses import fields
    from salmonpy.dft import DFT,Pseudopotential,SCFSettings
    from salmonpy.native_gs import NativeGSBackend
    from salmonpy.pseudopotentials import catalog_reference
    from test_nickel_soc import nickel_model
    from test_nickel_scf_mpi import InMemoryDFT
    comm=MPI.COMM_WORLD
    if comm.size!=args.size:raise ValueError('MPI size mismatch')
    work=args.output/args.worker
    if comm.rank==0:work.mkdir(exist_ok=False)
    comm.Barrier()
    pp=Pseudopotential.from_file(args.potential,atomic_number=28,local_channel=2,
                               reference=catalog_reference('pseudodojo-pbe-0.4-fr-standard-ni'))
    base=nickel_model(pp);calc=InMemoryDFT(**{f.name:getattr(base,f.name) for f in fields(DFT)})
    scf=SCFSettings(max_iterations=12,tolerance=1e-12,mixing_method='simple_dm',mixing_rate=.2,
                   initial_wavefunction='random',random_seed=0)
    backend=NativeGSBackend(args.build,work,decomposition=(args.size,1,1,1,1))
    context=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    reference=args.output/'continuous'/f'trace-{comm.rank}.npz'
    checkpoint=args.output/'save'/f'rt-{comm.rank}.bin'
    report={'phase':args.worker,'pid':os.getpid(),'SCF_executed_iterations':None}
    def drive(i):return np.array([1e-5*np.sin(.3*i),7e-6*(np.cos(.2*i)-1),1e-5*np.sin(.17*i)])
    fresh=args.worker=='continuous'
    with backend.open(context) as session:
        with session.run(calc,scf=scf,steps=4 if fresh else None,
                         restart=None if fresh else args.output/'gs',material_id='ni-soc-restart-probe') as gs:
            report['SCF_executed_iterations']=gs.scf_status['executed_iterations']
            if not fresh and report['SCF_executed_iterations']!=0:raise AssertionError('SCF replayed')
            if fresh:gs.save_scf(args.output/'gs')
            with gs.tddft_model(dt=.001,max_steps=10) as model:
                state=model.new_state('ni')
                try:
                    if args.worker=='resume':
                        if state.steps!=0:raise AssertionError('not a fresh RT state')
                        state.checkpoint(checkpoint,True)
                        if state.steps!=4:raise AssertionError('wrong restored RT step')
                        np.testing.assert_array_equal(state.A,drive(4))
                    report['RT_loaded_step']=state.steps
                    start=state.steps
                    currents=[];numbers=[]
                    for i in range(start+1,5 if args.worker=='save' else 9):
                        currents.append(state.step(drive(i)))
                        numbers.append(state.electron_number)
                    if args.worker=='save':state.checkpoint(checkpoint)
                    if fresh:np.savez(reference,current=currents,electron_number=numbers)
                    else:
                        with np.load(reference) as expected:
                            end=4 if args.worker=='save' else 8
                            np.testing.assert_array_equal(currents,expected['current'][start:end])
                            np.testing.assert_array_equal(numbers,expected['electron_number'][start:end])
                    report['RT_steps_executed']=state.steps-start
                    report['RT_final_step']=state.steps
                finally:state.close()
    reports=comm.gather(report,root=0)
    if comm.rank==0:(work/'result.json').write_text(json.dumps(reports[0],indent=2)+'\n')


if __name__=='__main__':
    try:main()
    except BaseException:
        if '--worker' in sys.argv:
            import traceback
            traceback.print_exc()
            from mpi4py import MPI
            MPI.COMM_WORLD.Abort(1)
        raise
