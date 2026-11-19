"""Native in-memory Si GS -> dt/norm scan; no GS archive/restart collection.

mpiexec -n 2 python -m salmonpy.examples.tddft_convergence BUILD_JSON NEW_OUTPUT
This small low-k grid is a connection diagnostic, not optical convergence proof.
"""
import argparse,json,hashlib
from pathlib import Path
from types import SimpleNamespace
from dataclasses import asdict
import numpy as np
from salmonpy.convergence import plan_dt_scan,run_dt_scan
from salmonpy.dft import DFT,Pseudopotential,SCFSettings
from salmonpy.native_gs import NativeGSBackend


def main():
    from mpi4py import MPI
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--duration-fs',type=float,default=2.)
    p.add_argument('--levels',type=int,default=6)
    p.add_argument('--grid',type=int,default=12)
    p.add_argument('--scf-iterations',type=int,default=150)
    p.add_argument('--norm-rtol',type=float,default=1e-5)
    p.add_argument('--xc',choices=('pz','tbmbj','libxc_pz','libxc_pbe'),default='pz')
    p.add_argument('--impulse',nargs=3,type=float,default=(1e-4,0,0))
    a=p.parse_args();comm=MPI.COMM_WORLD
    if comm.size not in (1,2):p.error('this example uses one or two internal MPI ranks')
    error=None
    if comm.rank==0:
        try:a.output.mkdir(parents=True,exist_ok=False)
        except Exception as exc:error=str(exc)
    error=comm.bcast(error,root=0)
    if error:raise ValueError(error)
    comm.Barrier()
    pseudo=Path(__file__).with_name('si_intro_data')/'Si_rps.dat'
    pp=Pseudopotential.from_file(pseudo,atomic_number=14,local_channel=2)
    positions=[[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0],
               [.25,.25,.25],[.25,.75,.75],[.75,.25,.75],[.75,.75,.25]]
    class InMemoryDFT(DFT):
        def input_text(self,scf=SCFSettings()):
            # Explicitly prohibit GS restart/wavefunction collection in this example.
            return super().input_text(scf).replace("write_gs_restart_data='all'","write_gs_restart_data='no'")
    calculation=InMemoryDFT((14,)*8,np.diag([10.26]*3),positions,(pp,),32,16,
                    (a.grid,)*3,(1,1,1),xc=a.xc)
    assert "write_gs_restart_data='no'" in calculation.input_text()
    backend=NativeGSBackend(a.build,a.output,decomposition=(1,comm.size,1,1,1))
    context=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
    plan=plan_dt_scan(material_id='silicon',duration_fs=a.duration_fs,levels=a.levels,
                      impulse=a.impulse,expected_electrons=32,
                      metadata={'xc':a.xc,'grid':list(calculation.grid),'kmesh':[1,1,1],
                                'pp_sha256':pp.sha256,'internal_ranks':comm.size,
                                'runtime_build_sha256':hashlib.sha256(a.build.read_bytes()).hexdigest()})
    counts={row['dt_au']:row['steps'] for row in plan['trials']}
    class Trial:
        def __init__(self,gs,dt):
            self.model=gs.tddft_model(dt=dt,max_steps=counts[dt])
            try:self.state=self.model.new_state('dt-trial')
            except BaseException:self.model.close();raise
        @property
        def current(self):return self.state.current
        @property
        def electron_number(self):return self.state.electron_number
        def step(self,A):return self.state.step(A)
        def close(self):
            try:self.state.close()
            finally:self.model.close()
    with backend.open(context) as session:
        gs=session.run(calculation,scf=SCFSettings(max_iterations=a.scf_iterations),material_id='silicon')
        try:
            report=run_dt_scan(plan,lambda dt:Trial(gs,dt),rtol=a.norm_rtol)
            report['gs_info']=asdict(gs.info)
            report['provisional']=not bool(gs.info.converged)
            assert not list(a.output.rglob('wfn.bin')), 'unexpected GS wavefunction output'
            report['GS_restart_collected']=False
            if comm.rank==0:
                (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
                print(json.dumps({key:report[key] for key in ('nonbreaking_dt_au','norm_acceptable_dt_au')},indent=2))
        finally:gs.close()

if __name__=='__main__':main()
