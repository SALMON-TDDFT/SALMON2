"""In-memory SCF pause/snapshot/continue; no GS dataset export or replay.

mpiexec -n 2 python -m salmonpy.examples.scf_resume BUILD_JSON NEW_OUTPUT
A small periodic one-Si cell is an execution check, not a physical benchmark.
"""
import argparse,json
from pathlib import Path
from types import SimpleNamespace
from dataclasses import asdict
import numpy as np
from salmonpy.dft import DFT,Pseudopotential,SCFSettings
from salmonpy.native_gs import NativeGSBackend


def main():
    from mpi4py import MPI
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--pause-at',type=int,default=3)
    p.add_argument('--iterations',type=int,default=8)
    p.add_argument('--xc',choices=('pz','libxc_pz','libxc_pbe'),default='pz')
    args=p.parse_args();comm=MPI.COMM_WORLD
    if comm.size not in (1,2):p.error('one/two internal MPI ranks in this example')
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
    calculation=InMemoryDFT((14,),np.diag([8.]*3),[[.3,.4,.45]],(pp,),4,2,(8,8,8),(1,1,1),xc=args.xc)
    assert "write_gs_restart_data='no'" in calculation.input_text()
    context=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    backend=NativeGSBackend(args.build,args.output,decomposition=(1,comm.size,1,1,1))
    with backend.open(context) as session:
        with session.run(calculation,scf=SCFSettings(max_iterations=args.iterations,tolerance=1e-30),
                         material_id='si-example',steps=args.pause_at) as state:
            with state.snapshot() as saved:
                state.continue_scf(steps=args.iterations-args.pause_at)
                linear=state.scf_status;info=state.info
                expected={q:state.array(q).values for q in ('density','occupations','orbitals','eigenvalues','forces')}
                state.restore(saved)
                restored=state.scf_status
                assert restored['executed_iterations']==linear['executed_iterations']
                assert restored['initializations']==1
                state.continue_scf(steps=args.iterations-args.pause_at)
                assert state.info==info
                for q,array in expected.items():np.testing.assert_array_equal(state.array(q).values,array)
                result={'status':'pass','linear_counters':linear,'restored_counters':restored,
                        'restored_result_exact':True,'final_info':asdict(info),
                        'snapshot_storage':'Fortran memory, same live owner/topology','GS_datasets_collected':False}
            assert not list(args.output.rglob('wfn.bin'))
            if comm.rank==0:
                (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(args.output)

if __name__=='__main__':main()
