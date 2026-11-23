"""Manual Ni SOC A-to-current state-isolation/checkpoint probe, not physical convergence.

Launch 3*SIZE ranks. GS arrays stay in RAM; electronic RT checkpoints are saved.
The deliberately paused four-iteration GS is used only for software comparison.
"""
import argparse
from dataclasses import fields
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
from mpi4py import MPI
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from salmonpy.dft import DFT,Pseudopotential,SCFSettings
from salmonpy.native_gs import NativeGSBackend
from salmonpy.pseudopotentials import catalog_reference
from test_nickel_soc import nickel_model
from test_nickel_scf_mpi import InMemoryDFT


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--potential',type=Path,required=True)
    p.add_argument('--size',type=int,choices=(1,2),default=1)
    p.add_argument('--baseline',type=Path,help='old audited build for independent controls')
    p.add_argument('--spin-stride',type=int,help='explicit observation stride for this software probe')
    p.add_argument('--without-checkpoint',action='store_true',help='only compare live states; do not save electronic data')
    args=p.parse_args();world=MPI.COMM_WORLD
    if world.size!=3*args.size:raise ValueError('launch 3*size ranks')
    group=world.rank//args.size;comm=world.Split(group,world.rank)
    work=args.output/str(group)
    if comm.rank==0:work.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    pp=Pseudopotential.from_file(args.potential,atomic_number=28,local_channel=2,
                               reference=catalog_reference('pseudodojo-pbe-0.4-fr-standard-ni'))
    class ProbeDFT(InMemoryDFT):
        def input_text(self,scf=SCFSettings()):
            text=super().input_text(scf)
            return text if args.spin_stride is None else text+f'\n&analysis\n out_rt_spin_step={args.spin_stride}\n/\n'
    if args.spin_stride is not None and args.spin_stride<1:raise ValueError('positive spin stride required')
    base=nickel_model(pp);calc=ProbeDFT(**{f.name:getattr(base,f.name) for f in fields(DFT)})
    scf=SCFSettings(max_iterations=12,tolerance=1e-12,mixing_method='simple_dm',mixing_rate=.2,
                   initial_wavefunction='random',random_seed=0)
    backend=NativeGSBackend(args.baseline if group and args.baseline else args.build,work,decomposition=(args.size,1,1,1,1))
    context=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    traces={};resume_exact=False
    def drive(i,k):
        return np.array([1e-5*np.sin(.3*i),7e-6*(np.cos(.2*i)-1),1e-5*np.sin(.17*i)])*(1 if k==1 else -.7)
    with backend.open(context) as session:
        with session.run(calc,scf=scf,steps=4,material_id='ni-soc-probe') as gs:
            with gs.tddft_model(dt=.001,max_steps=10) as model:
                keys=(1,2) if group==0 else (group,)
                states={k:model.new_state(k) for k in keys}
                try:
                    for k in keys:traces[k]=[]
                    for i in range(1,9):
                        for k,state in states.items():
                            traces[k].append(state.step(drive(i,k)))
                            if group==0 and i==4 and not args.without_checkpoint:state.checkpoint(work/f'{k}-{comm.rank}.bin')
                    if group==0 and not args.without_checkpoint:
                        for k in keys:
                            resumed=model.new_state('resume-'+str(k))
                            try:
                                resumed.checkpoint(work/f'{k}-{comm.rank}.bin',True)
                                if resumed.steps!=4:raise AssertionError('wrong restored step')
                                for i in range(5,9):
                                    np.testing.assert_array_equal(resumed.step(drive(i,k)),traces[k][i-1])
                            finally:resumed.close()
                        resume_exact=True
                finally:
                    for state in states.values():state.close()
    files={}
    if comm.rank==0:
        for k,state in states.items():
            files[k]=[[line for line in (state.work/name).read_text().splitlines()
                       if line.strip() and not line.startswith('#')] for name in state.output_files]
    gathered=world.gather(traces if comm.rank==0 else None,root=0)
    all_files=world.gather(files if comm.rank==0 else None,root=0)
    if world.rank==0:
        for k in (1,2):np.testing.assert_array_equal(gathered[0][k],gathered[k*args.size][k])
        for k in (1,2):
            if all_files[0][k]!=all_files[k*args.size][k]:raise AssertionError('independent observation files differ')
        result={'status':'passed','internal_MPI_ranks':args.size,'states':36,'electrons':18,
                'temperature_k':300,'kmesh':[2,2,2],'independent_currents_exact':True,
                'direct_disk_resume_exact':resume_exact,'observation_files_exact':True,
                'baseline_compared':bool(args.baseline),'spin_stride':args.spin_stride,'scope':'short paused-GS software probe; not physical convergence'}
        (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
        print('Ni SOC TDDFT PASS',args.output/'result.json')
    comm.Free()


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback
        traceback.print_exc();MPI.COMM_WORLD.Abort(1)
