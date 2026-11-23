"""Explicit SOC material GS/response validation; no calculation without --run.

BUILD_JSON MATERIAL_JSON NEW_OUTPUT [--run] [--threads 4]
The input JSON supplies geometry and pinned local FR-UPF references. No download,
source build or job submission. Stores GS and RT data only in explicit run mode.
"""
import argparse
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def load_model(config):
    from salmonpy.dft import DFT,Pseudopotential,SCFSettings,_restart_output_text
    from salmonpy.pseudopotentials import catalog_reference
    path=Path(config).resolve();spec=json.loads(path.read_text())
    if spec.get('schema')!='salmonpy.soc-validation/1':raise ValueError('SOC validation schema required')
    ps=[]
    for entry in spec['pseudopotentials']:
        ref=catalog_reference(entry['catalog_id'])
        ps.append(Pseudopotential.from_file(path.parent/entry['file'],atomic_number=ref.atomic_number,
                  local_channel=entry['local_channel'],reference=ref))
    class ObservedDFT(DFT):
        def input_text(self,scf=SCFSettings()):
            return _restart_output_text(super().input_text(scf),False)+'\n&analysis\n out_rt_spin_step=1\n/\n'
    model=ObservedDFT(spec['numbers'],spec['cell_bohr'],spec['fractional_positions'],ps,
                      electrons=spec['electrons'],states=spec['states'],grid=spec['grid'],kmesh=spec['kmesh'],
                      xc=spec['xc'],spin_orbit=True,electronic_temperature_k=spec['temperature_k'])
    from salmonpy.ld import _positive,_array
    from salmonpy.sbe import _integer
    _integer(spec['response']['steps'],'response steps',2,2**31-2)
    _positive(spec['response']['dt'],'response dt')
    _array(spec['response']['A'],'response A',(3,))
    settings=SCFSettings(**spec['scf'])
    return model,settings,spec


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('config',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--run',action='store_true');p.add_argument('--threads',type=int,default=4)
    p.add_argument('--restart',type=Path,help='load existing SCF state including mixing history')
    p.add_argument('--mixing-rate',type=float,help='change effective mixing coefficient after loading, preserving history')
    p.add_argument('--scf-limit',type=int,help='explicitly review and increase the owned iteration budget after loading')
    p.add_argument('--compatible-runtime',choices=('declared',),default=False)
    args=p.parse_args()
    model,scf,spec=load_model(args.config)
    if not args.run:
        print(model.input_text(scf));print('Input check only; no MPI initialization or calculation.');
        if args.scf_limit is not None:print('Requested runtime SCF cap:',args.scf_limit)
        return
    if args.threads<1:raise ValueError('positive threads required')
    os.environ.setdefault('OMP_NUM_THREADS',str(args.threads))
    os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
    import numpy as np
    from mpi4py import MPI
    from types import SimpleNamespace
    from salmonpy.native_gs import NativeGSBackend
    from salmonpy.errors import SCFReviewRequired
    comm=MPI.COMM_WORLD;root=args.output.resolve()
    if comm.rank==0:root.mkdir(parents=True,exist_ok=False)
    comm.Barrier();started=time.perf_counter()
    context=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    backend=NativeGSBackend(args.build,root,decomposition=(comm.size,1,1,1,1))
    def report(name,value):
        if comm.rank==0:
            (root/(name+'.json')).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    report('configuration',{**spec,'MPI_ranks':comm.size,'OMP_NUM_THREADS':os.environ['OMP_NUM_THREADS'],
                            'runtime_sha256':backend.sha256,'scope':'SCF/interface validation; grid/k/dt/geometry convergence unverified'})
    batch=spec.get('scf_batch',25)
    if type(batch) is not int or batch<1:raise ValueError('positive SCF batch required')
    with backend.open(context) as session:
        with session.run(model,scf=scf,steps=None if args.restart else min(batch,scf.max_iterations),
                         restart=args.restart,compatible_runtime=args.compatible_runtime,
                         allow_unconverged=bool(args.restart and (args.scf_limit or args.mixing_rate is not None)),material_id=spec['material_id']) as gs:
            if args.scf_limit is not None:gs.extend_scf_limit(args.scf_limit)
            if args.mixing_rate is not None:gs.set_scf_mixing(args.mixing_rate)
            limit=gs.scf_status['max_iterations']
            report('resume-control',{'initial_iterations':gs.iterations,'executed_iterations_before_run':gs.scf_status['executed_iterations'],
                                     'mixing_rate':args.mixing_rate,'iteration_limit':limit,'restart':None if args.restart is None else str(args.restart)})
            while True:
                info=gs.info
                value={'iterations':info.iterations,'converged':info.converged,'residual':info.residual,
                       'energy_hartree':info.energy_hartree,'elapsed_seconds':time.perf_counter()-started}
                report('progress',value)
                if comm.rank==0:print('SCF',value,flush=True)
                if info.converged or info.iterations>=limit:break
                try:gs.continue_scf(steps=min(batch,limit-info.iterations))
                except SCFReviewRequired:pass  # state remains live; save diagnostics then stop below
            gs.save_scf(root/'gs-state')
            if not info.converged:
                report('result',{'status':'SCF review required','GS_saved':True,'TDDFT_started':False,**value})
                if comm.rank==0:print('SCF not converged: inspect gs-state/diagnostics; TDDFT not started.',flush=True)
                return
            steps=spec['response']['steps'];dt=spec['response']['dt']
            if type(steps) is not int or steps<2:raise ValueError('at least two RT steps required')
            cut=steps//2;drive=np.asarray(spec['response']['A'],dtype=np.float64)
            result={'status':'passed','GS':value,'MPI_ranks':comm.size,'response_steps':steps,'dt':dt}
            with gs.tddft_model(dt=dt,max_steps=steps) as response:
                checked=response.new_state('checked',execution_mode='checked')
                fast=response.new_state('performance',execution_mode='performance')
                currents=[];numbers=[];fast_currents=[];resumed_currents=[]
                try:
                    for i in range(1,steps+1):
                        currents.append(checked.step(drive));numbers.append(checked.electron_number)
                        fast_currents.append(fast.step(drive))
                        np.testing.assert_allclose(fast_currents[-1],currents[-1],rtol=2e-10,atol=1e-12)
                        if i==cut:checked.checkpoint(root/f'rt-{comm.rank}.bin')
                    resumed=response.new_state('direct-resume')
                    try:
                        resumed.checkpoint(root/f'rt-{comm.rank}.bin',True)
                        if resumed.steps!=cut:raise AssertionError('wrong restored RT step')
                        for i in range(cut,steps):
                            resumed_currents.append(resumed.step(drive))
                            np.testing.assert_allclose(resumed_currents[-1],currents[i],rtol=2e-10,atol=1e-12)
                    finally:resumed.close()
                    result['performance_max_abs_current_difference']=float(np.max(np.abs(np.asarray(fast_currents)-currents)))
                    result['resume_max_abs_current_difference']=float(np.max(np.abs(np.asarray(resumed_currents)-currents[cut:])))
                    result['electron_number_min']=float(np.min(numbers));result['electron_number_max']=float(np.max(numbers))
                finally:
                    fast.close();checked.close()
            if comm.rank==0:np.savez(root/'observations.npz',charge_current=currents,electron_number=numbers,
                                     performance_current=fast_currents,resumed_current=resumed_currents)
            result['elapsed_seconds']=time.perf_counter()-started
            result['physical_convergence_verified']=False
            report('result',result)
            if comm.rank==0:print('SOC material GS/response validation PASS',root/'result.json',flush=True)


if __name__=='__main__':
    try:main()
    except BaseException:
        if '--run' in sys.argv:
            import traceback
            traceback.print_exc()
            from mpi4py import MPI
            MPI.COMM_WORLD.Abort(1)
        raise
