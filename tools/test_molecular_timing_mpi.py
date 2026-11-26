"""Manual half-step audit; physical calculations are not required CI tests.

Response modes use one controller rank plus molecular ranks. Field mode uses
FIELD_RANKS plus one label rank. Existing GS only; no SCF or nested launcher.
"""
import argparse
from pathlib import Path
import json
import shutil
import numpy as np


def field_audit(args,world):
    from salmonpy.distributed import ParallelPlan,ParallelContext,SolverGroup
    from salmonpy.coupled import NativeEHField
    from salmonpy.geometry import Grid,Geometry
    from salmonpy.fdtd import FDTD,PlaneWave
    from salmonpy.fields import AcosPulse
    from salmonpy.examples.molecular_fdtd_mpi_response import abort_on_error
    nf=world.size-1
    plan=ParallelPlan((SolverGroup('field',tuple(range(nf)),role='field'),
                       SolverGroup('labels',(nf,),point_ids=(10,11))))
    with ParallelContext.create(world,plan) as context,abort_on_error(world):
        if context.group.role=='field':
            calculation=FDTD(Geometry(Grid((40,40,40),(2.,)*3)),(),.003,args.steps,
                source=PlaneWave(AcosPulse(0.,20.,.15),0,-20.))
            points={10:(20,20,20),11:(19,20,20)}
            query=np.array([np.array(p)-[x,y,z] for p in points.values()
                            for z in (0,1) for y in (0,1) for x in (0,1)],dtype=np.int32)
            ah=[];eh=[];hh=[];centered=[];previous=None
            with NativeEHField(context,calculation,args.build,args.output/'field',points=points,grid=(nf,1,1)) as field:
                for n in range(1,args.steps+1):
                    old=field.A.copy();current=np.zeros((2,3))
                    if n>1:current[1]=[1.e-6*n,-2.e-6*n,3.e-6*n]
                    field.advance(current)
                    np.testing.assert_allclose(field.A-old,-field.dt*field.E,atol=1.e-20,rtol=2.e-13)
                    samples=field.sample_fields();raw=field.sample_magnetic_field().values[0]
                    component=field.component_fields(query)[1].reshape((2,2,2,2,3))
                    avg=np.empty_like(raw)
                    # Independent oracle for the native transverse spatial sampling.
                    for p in range(2):
                        for d in range(3):
                            rows=[]
                            for z in (0,1):
                                for y in (0,1):
                                    for x in (0,1):
                                        if (x,y,z)[d]==0:rows.append(component[p,z,y,x,d])
                            avg[p,d]=sum(rows)/4.
                    if previous is not None:
                        np.testing.assert_allclose(avg,.5*(previous+raw),atol=1.e-20,rtol=2.e-13)
                    np.testing.assert_array_equal(samples['E'].times,[n*field.dt])
                    np.testing.assert_array_equal(samples['A'].times,[(n+.5)*field.dt])
                    ah.append(field.A.copy());eh.append(field.E.copy());hh.append(raw.copy());centered.append(avg)
                    previous=raw.copy()
            if world.rank==0:
                assert np.max(np.abs(hh))>1.e-12,'magnetic test signal is zero'
                np.savez(args.output/'fields.npz',A=ah,E=eh,H_raw=hh,H_centered=centered)
                (args.output/'summary.json').write_text(json.dumps({'E_half_step':'2*n','A_half_step':'2*n+1',
                    'raw_H_half_step':'2*n+1','centered_H_half_step':'2*n','H_bracket_average_verified_from_step':2,
                    'A_difference_matches_E':True,'nonzero_magnetic_signal':True},indent=2)+'\n')


def response_audit(args,world):
    from salmonpy.isolated_tddft import run_molecular_mpi,MolecularMPIController
    from salmonpy.near_field import GridFieldSeries
    from salmonpy.examples.molecular_fdtd_mpi_response import abort_on_error
    control=world.Dup();group=world.Split(0 if world.rank==0 else 1,world.rank)
    work=args.output/'molecule-work';dt=.003
    try:
        if world.rank==0:
            work.mkdir()
            for name in ('input.inp','C_rps.dat','H_rps.dat'):shutil.copy2(args.bundle/name,work/name)
            shutil.copytree(args.bundle/'data_for_restart',work/'data_for_restart')
            zero=np.zeros((args.steps+2,2,2,2,3))
            GridFieldSeries((-40.,)*3,(80.,)*3,np.arange(args.steps+2)*dt,zero,zero).write(work/'salmonpy-near-eh.dat')
        world.Barrier()
        with abort_on_error(world):
            if world.rank:
                run_molecular_mpi(args.build,group,control,controller_rank=0,work_directory=work,
                                  decomposition=(1,group.size,1,1,1))
            else:
                controller=MolecularMPIController(control,solver_root=1,work_directory=work,dt=dt,
                    max_steps=args.steps,decomposition=(1,world.size-1,1,1,1))
                def row(obs):return np.array([*obs.dipole,*obs.electronic_change,obs.electron_number])
                start=0;values=[row(controller.observables)]
                if args.resume:
                    controller.restore_direct_checkpoint(args.resume/'checkpoint')
                    start=controller.step_index;values=list(np.load(args.resume/'response.npz')['observations'])
                    np.testing.assert_allclose(row(controller.observables),values[-1],atol=2.e-12,rtol=2.e-10)
                restored=controller.worker_statistics()
                stop=args.stop_after or args.steps
                for n in range(start+1,stop+1):
                    left=np.array([1.,-.5,.75])*n*dt
                    right=np.array([1.,-.5,.75])*(n+1)*dt
                    extra=np.zeros(3)
                    if args.mode=='uniform':extra=left;left=np.zeros(3);right=np.zeros(3)
                    elif args.mode=='shifted':extra=.5*(left+right);left=np.zeros(3);right=np.zeros(3)
                    elif args.mode=='right' and not args.magnetic:right=right*100.
                    E=np.empty((2,2,2,2,3));E[0]=left;E[1]=right
                    H=np.zeros_like(E)
                    if args.magnetic:
                        H[0]=np.array([.25,.5,-.25])*1.e-3*n*dt
                        H[1]=np.array([.25,.5,-.25])*1.e-3*(n+1)*dt
                        if args.mode=='right':H[1]*=100.
                    window=GridFieldSeries((-40.,)*3,(80.,)*3,[0.,dt],E,H,time_origin=n*dt)
                    values.append(row(controller.step_grid_field(window,step_index=n,E_extra=extra)))
                if args.stop_after:controller.save_direct_checkpoint(args.output/'checkpoint')
                summary={'mode':args.mode,'magnetic':args.magnetic,'steps':stop,'restore_statistics':restored,
                         'final_statistics':controller.worker_statistics(),'clock':controller.clock,
                         'physical_convergence':False}
                controller.close()
                np.savez(args.output/'response.npz',observations=values)
                (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        assert world.allreduce(1)==world.size
    finally:group.Free();control.Free()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('field','grid','uniform','shifted','right'))
    p.add_argument('build',type=Path);p.add_argument('bundle',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--magnetic',action='store_true',help='also audit orbital magnetic coupling; right mode changes only future H')
    p.add_argument('--steps',type=int,default=12);p.add_argument('--stop-after',type=int);p.add_argument('--resume',type=Path)
    args=p.parse_args()
    from mpi4py import MPI
    world=MPI.COMM_WORLD
    if world.size<2:p.error('at least two MPI ranks required')
    from salmonpy.examples.molecular_fdtd_mpi_response import abort_on_error
    with abort_on_error(world):
        if world.rank==0:args.output.mkdir(parents=True,exist_ok=False)
        world.Barrier()
        if args.mode=='field':field_audit(args,world)
        else:response_audit(args,world)
    if world.rank==0:print('Half-step audit finished:',args.output)


if __name__=='__main__':main()
