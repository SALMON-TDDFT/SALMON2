"""In-memory SCF pause/snapshot/continue; no GS dataset export or replay.

mpiexec -n 2 python -m salmonpy.examples.scf_resume BUILD_JSON NEW_OUTPUT
A small periodic one-Si cell is an execution check, not a physical benchmark.
"""
import argparse,json
from pathlib import Path
from types import SimpleNamespace
from dataclasses import asdict
import numpy as np
from salmonpy.dft import DFT,Pseudopotential,SCFSettings,SymmetryOperations
from salmonpy.native_gs import NativeGSBackend
from salmonpy.checkpointing import execute_checkpoint,restore_checkpoint,CheckpointState


def main():
    from mpi4py import MPI
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--pause-at',type=int,default=3)
    p.add_argument('--iterations',type=int,default=8)
    parallel=p.add_mutually_exclusive_group()
    parallel.add_argument('--k-parallel',action='store_true',help='k-only layout with a fixed 2x2x2 mesh')
    parallel.add_argument('--orbital-parallel',action='store_true',help='Gamma orbital layout for ordinary XC or explicit fast HSE')
    parallel.add_argument('--spatial-parallel',action='store_true',help='ordinary-XC x-domain layout with a fixed Gamma mesh')
    parallel.add_argument('--decomposition',type=int,nargs=5,help='explicit k,orbital,x,y,z layout; fast HSE uses Gamma orbital/spatial mixing')
    p.add_argument('--repartition',action='store_true',help='opt into a supported pure-layout restart with a new rank count')
    p.add_argument('--disk',action='store_true',help='exercise disk save/restore; requires --save-scf')
    p.add_argument('--protected',action='store_true',help='requires native SCF candidate protocol1')
    p.add_argument('--symmetry',choices=('none','c2z'),default='none')
    p.add_argument('--xc',choices=('pz','libxc_pz','libxc_pbe','hse06'),default='pz')
    disk=p.add_mutually_exclusive_group()
    disk.add_argument('--save-scf',type=Path,help='save the paused SCF boundary for a fresh process')
    disk.add_argument('--restart',type=Path,help='load a prior --save-scf boundary, then continue')
    args=p.parse_args();comm=MPI.COMM_WORLD
    if args.disk and args.save_scf is None:p.error('--disk requires --save-scf')
    if args.repartition and args.restart is None:p.error('--repartition requires --restart')
    if comm.size not in (1,2) and args.decomposition is None:p.error('use --decomposition for more than two ranks')
    if not 0<args.pause_at<args.iterations:p.error('pause-at must be between zero and total iterations')
    error=None
    if comm.rank==0:
        try:args.output.mkdir(parents=True,exist_ok=False)
        except Exception as exc:error=str(exc)
    error=comm.bcast(error,root=0)
    if error:raise ValueError(error)
    comm.Barrier()
    pp=Pseudopotential.from_file(Path(__file__).with_name('si_intro_data')/'Si_rps.dat',atomic_number=14,local_channel=2)
    class InMemoryDFT(DFT):
        def input_text(self,scf=SCFSettings()):
            return super().input_text(scf).replace("write_gs_restart_data='all'","write_gs_restart_data='no'")
    symmetry=None
    if args.symmetry!='none':
        ops=np.zeros((2,3,4));ops[:,:,:3]=np.eye(3);ops[1,:2,:2]=-np.eye(2)
        symmetry=SymmetryOperations(ops,directions='yyn')
    mesh=(1,1,1) if args.xc=='hse06' and args.decomposition is not None else (2,2,2) if args.decomposition is not None or symmetry is not None or args.k_parallel or (args.repartition and not (args.orbital_parallel or args.spatial_parallel)) else (comm.size,1,1) if args.xc=='hse06' and not (args.orbital_parallel or args.spatial_parallel) else (1,1,1)
    calculation=InMemoryDFT((14,),np.diag([8.]*3),[[0,0,0]] if symmetry is not None else [[.3,.4,.45]],
        (pp,),4,2,(8,8,8),mesh,xc=args.xc,symmetry=symmetry)
    assert "write_gs_restart_data='no'" in calculation.input_text()
    context=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    decomposition=(comm.size,1,1,1,1) if (args.xc=='hse06' and not (args.orbital_parallel or args.spatial_parallel)) or args.k_parallel or (args.repartition and not (args.orbital_parallel or args.spatial_parallel)) else (1,comm.size,1,1,1)
    if args.spatial_parallel:decomposition=(1,1,comm.size,1,1)
    if args.decomposition:decomposition=tuple(args.decomposition)
    backend=NativeGSBackend(args.build,args.output,decomposition=decomposition,hse_orbital=args.xc=='hse06' and args.orbital_parallel,hse_spatial=args.xc=='hse06' and (args.spatial_parallel or args.decomposition is not None))
    with backend.open(context) as session:
        if args.restart is not None:
            with session.run(calculation,scf=SCFSettings(max_iterations=args.iterations,tolerance=1e-30),
                             material_id='si-example',restart=args.restart,repartition=args.repartition) as state:
                loaded=state.scf_status
                assert loaded['executed_iterations']==0
                state.continue_scf(steps=args.iterations)
                result={'status':'pass','loaded_counters':loaded,'final_counters':state.scf_status,
                        'final_info':asdict(state.info),'snapshot_storage':'same-build pure-layout repartitioned disk checkpoint' if args.repartition else 'same-build/topology disk checkpoint',
                        'GS_datasets_collected':True}
                if comm.rank==0:(args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
            return
        with session.run(calculation,scf=SCFSettings(max_iterations=args.iterations,tolerance=1e-30),
                         material_id='si-example',steps=args.pause_at) as state:
            if args.save_scf is not None and not args.disk:state.save_scf(args.save_scf)
            groups=[{'name':'scf','ranks':list(range(comm.size)),'point_ids':['si-example']}]
            target=state.checkpoint_target(directory=args.save_scf) if args.disk else state.checkpoint_target()
            checkpoint=execute_checkpoint(comm,groups,states={'si-example':target})
            from contextlib import nullcontext
            saved_context=nullcontext(checkpoint['states']['si-example']) if args.disk else checkpoint['states']['si-example']
            with saved_context as saved:
                state.continue_scf(steps=args.iterations-args.pause_at)
                linear=state.scf_status;info=state.info
                expected={q:state.array(q).values for q in ('density','occupations','orbitals','eigenvalues','forces')}
                targets={'si-example':state.checkpoint_target(directory=args.save_scf) if args.disk else state.checkpoint_target(saved)}
                if args.protected:
                    prepared=restore_checkpoint(comm,groups,saved=checkpoint['compatibility'],
                        states=targets,live=CheckpointState(states={'si-example':state}))
                    for candidate in prepared['states'].values():candidate.publish()
                else:
                    execute_checkpoint(comm,groups,operation='restore',saved=checkpoint['compatibility'],states=targets)
                restored=state.scf_status
                assert restored['executed_iterations']==linear['executed_iterations']
                assert restored['initializations']==1
                state.continue_scf(steps=args.iterations-args.pause_at)
                assert state.info==info
                for q,array in expected.items():np.testing.assert_array_equal(state.array(q).values,array)
                result={'status':'pass','linear_counters':linear,'restored_counters':restored,
                        'restored_result_exact':True,'final_info':asdict(info),
                        'snapshot_storage':'same-build/topology disk checkpoint' if args.disk else 'Fortran memory, same live owner/topology',
                        'checkpoint_api':'common protected staging/publication' if args.protected else 'common direct execution','GS_datasets_collected':args.save_scf is not None,'symmetry':args.symmetry,
                        'force_full_reduced_agreement':'not assessed by this execution sample'}
            assert not list(args.output.rglob('wfn.bin'))
            if comm.rank==0:
                (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(args.output)

if __name__=='__main__':main()
