"""Shared reduced/full-k HSE numerical probe: BUILD FIXTURE NEW_OUTPUT.

Launch 3*--group-size ranks. Group0 shares A/B; groups1/2 use independent legacy
states as controls. Every context sees a symmetry-preserving Cartesian z drive.
Explicit local build/inputs only; synthetic ownership test, not convergence.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend, _history_info


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','fixture','output'):p.add_argument(name,type=Path)
    p.add_argument('--group-size',type=int,default=2)
    p.add_argument('--symmetry',action='store_true')
    p.add_argument('--steps',type=int,default=32)
    p.add_argument('--baseline',type=Path)
    p.add_argument('--no-chdir',action='store_true')
    args=p.parse_args();world=MPI.COMM_WORLD
    assert world.size==3*args.group_size and args.steps>=16
    group=world.rank//args.group_size;comm=world.Split(group,world.rank)
    out=args.output.resolve();work=out/str(group)
    if comm.rank==0:work.mkdir(parents=True,exist_ok=False)
    comm.Barrier()
    context=SimpleNamespace(solver_comm=comm,point_ids=(1,2) if group==0 else (group,),closed=False,_dependents=0)
    anchor=Path.cwd()
    from unittest.mock import patch
    guard=patch('os.chdir',side_effect=AssertionError('shared HSE host chdir forbidden')) if args.no_chdir and group==0 else None
    if guard:guard.start()
    initial=FastHSEMaterialProvider(args.fixture,max_steps=args.steps).prepare(context)
    backend=NativeFastHSEBackend(args.baseline if group and args.baseline else args.build,work,allow_symmetry=args.symmetry)
    states=[];model=None;traces={};histories={};resumes={}
    def drive(index,which):
        return np.array([0.,0.,1e-4*np.sin(.3*index)*(1 if which==1 else -.7)])
    try:
        if group==0:
            model=backend.open_model(initial,.001,context)
            states=[model.new_state(1),model.new_state(2)]
        else:states=[backend.create_state(initial,group,.001,context)]
        if group==0 and args.symmetry:
            from salmonpy.errors import InputError
            before=[(s.current.copy(),s.A.copy(),s.steps) for s in states]
            try:states[1].step([1e-4,0.,0.])
            except InputError as exc:assert 'breaks' in str(exc)
            else:raise AssertionError('symmetry-breaking A was accepted')
            for state,(j,A,step) in zip(states,before):
                np.testing.assert_array_equal(state.current,j);np.testing.assert_array_equal(state.A,A)
                assert state.steps==step and not state.failed
        for state in states:traces[state.point_id if group==0 else group]=[]
        # Independent control create_state has no point_id attribute.
        for index in range(1,args.steps+1):
            for offset,state in enumerate(states):
                which=offset+1 if group==0 else group
                traces[which].append(np.r_[state.step(drive(index,which)),state.energy])
                if group==0 and index==7+which:
                    state.checkpoint(work/f'{which}-{comm.rank}.bin',False)
        if group==0:
            for state in states:histories[state.point_id]=state._call(_history_info,state.lib).tolist()
            saved_b=(states[1].A.copy(),states[1].current.copy(),states[1].energy,states[1].steps)
            for which in (1,2):
                fresh=model.new_state(which);states.append(fresh)
                fresh.checkpoint(work/f'{which}-{comm.rank}.bin',True)
                trace=[]
                for index in range(8+which,args.steps+1):trace.append(np.r_[fresh.step(drive(index,which)),fresh.energy])
                np.testing.assert_array_equal(trace,traces[which][7+which:])
                resumes[which]=trace
            np.testing.assert_array_equal(states[1].A,saved_b[0]);np.testing.assert_array_equal(states[1].current,saved_b[1])
            assert (states[1].energy,states[1].steps)==saved_b[2:]
        if group==0:
            for state in reversed(states):state.close()
            states=[]
            fresh=model.new_state(3);states.append(fresh)
            np.testing.assert_array_equal(np.r_[fresh.step(drive(1,1)),fresh.energy],traces[1][0])
            model.close()  # A surviving lease still works after creation is disabled.
            np.testing.assert_array_equal(np.r_[fresh.step(drive(2,1)),fresh.energy],traces[1][1])
        np.savez(work/f'trace-{comm.rank}.npz',**{f'point{k}':v for k,v in traces.items()})
    finally:
        for state in reversed(states):state.close()
        if model is not None:model.close()
    assert context._dependents==0
    assert Path.cwd()==anchor
    if guard:guard.stop()
    gathered=world.allgather(traces)
    for which in (1,2):
        ref=gathered[which*args.group_size][which]
        for r in range(args.group_size):np.testing.assert_array_equal(gathered[r][which],ref)
    if world.rank==0:
        assert not np.array_equal(traces[1],traces[2])
        (out/'result.json').write_text(json.dumps({'status':'pass','symmetry':args.symmetry,
            'group_size':args.group_size,'steps':args.steps,'replayed_steps':0,
            'independent_controls_exact':True,'new_runtime_no_chdir':args.no_chdir,'previous_runtime_compared':bool(args.baseline),'resume_exact':True,'histories':histories},indent=2)+'\n')
        print('Shared symmetry HSE: controls and direct resume PASS',flush=True)
    comm.Free()


if __name__=='__main__':
    try:main()
    except BaseException:
        import traceback
        traceback.print_exc()
        MPI.COMM_WORLD.Abort(1)
