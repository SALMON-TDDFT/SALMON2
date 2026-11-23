"""Protected Multiscale1D + LD/SBE/custom example; six ranks, explicit callbacks.

mpiexec -n 6 python -m salmonpy.examples.mpi_protected_step NEW_OUTPUT
Synthetic connection test, not a physical composite-material model. Native whole
mixed-step support is not inferred from component numerical recovery protocols.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from salmonpy.coupled import CoupledStepper
from salmonpy.distributed import ParallelPlan,ParallelContext,SolverGroup
from salmonpy.material_response import MaterialBinding
from salmonpy.timing import TimeSpec,InputPort,OutputPort
from salmonpy.maxwell import Multiscale1D
from salmonpy.fields import AcosPulse
from salmonpy.stepping import StepTarget,StepRecoveryError
from salmonpy.examples.mpi_mixed_multiscale import models


class Provider:
    def __init__(self,kind,model=None):self.kind=kind;self.model=model;self.identity='protected-example/'+kind
    def prepare(self,context):return self


class Backend:
    execution='collective';reduction='leader';identity='protected-example-response/1'
    def create_state(self,initial,point_id,dt,context):return State(initial,point_id,dt,context.solver_comm)


class State:
    def __init__(self,initial,point,dt,comm):
        self.kind=initial.kind;self.point=point;self.dt=dt;self.comm=comm
        self.raw=None if initial.model is None else initial.model.new_state(dt=dt)
        self.A=np.zeros(3);self.current=np.zeros(3);self.n=0
        self.fail_next=False;self.fail_capture=False;self.fail_rollback=False;self.fail_release=False
        self.timing=TimeSpec(dt,(InputPort('A'),),(OutputPort('j'),),integrator='explicit protected example '+self.kind)
    def step(self,A):
        if self.kind=='ld':self.current=self.raw.step(-(A-self.A)/self.dt)[1]
        elif self.kind=='sbe':self.current=self.raw.step(A)
        else:self.current=.9*self.current+.01*A+np.array([0.,0.,1e-6*self.point])
        self.A=A.copy();self.n+=1
        if self.fail_next:
            self.fail_next=False
            raise ValueError('injected late custom mutation')
        return self.current.copy()
    def capture(self):
        if self.fail_capture:raise ValueError('injected capture failure')
        return (None if self.raw is None else self.raw.snapshot(),self.A.copy(),self.current.copy(),self.n,self.A,self.current)
    def rollback(self,token):
        raw,A,current,self.n,original_A,original_current=token
        if not np.array_equal(original_A,A):original_A[:]=A
        if not np.array_equal(original_current,current):original_current[:]=current
        self.A=original_A;self.current=original_current
        if self.raw is not None:self.raw.restore(raw)
        if self.fail_rollback:raise ValueError('injected rollback failure')
    def release(self,token):
        # Collective even on the rank whose capture failed and supplied None.
        self.comm.allgather(token is None)
        if self.fail_release:raise ValueError('injected release failure')
    def close(self):
        if self.raw is not None:self.raw.close()


def targets(coupled):
    if coupled.field is not None:
        f=coupled.field
        def capture():return f.state.snapshot(),f.closed,f.failed
        def restore(t):
            f.state.restore(t[0]);f.closed=t[1];f.failed=t[2]
        return [StepTarget('field',f,capture,restore)]
    return [StepTarget('point-'+str(p),s,s.capture,s.rollback,s.release)
            for p,s in coupled.response.states.items()]


def signature(coupled):
    # Validation-only inspection of internal state; not a user persistence API.
    r=coupled.response
    value=(coupled.step_index,id(coupled.J),coupled.J.tobytes(),r.step_index,id(r._A),r._A.tobytes(),r.failed,
           tuple((p,id(a.A),a.A.tobytes()) for p,a in r._adapters.items()))
    if coupled.field is not None:
        f=coupled.field;s=f.state.snapshot()
        return value+(s._data,s._integrity,f.closed,f.failed)
    items=[]
    for p,s in r.states.items():
        token=s.capture();raw=token[0]
        items.append((p,id(s.A),id(s.current),token[1].tobytes(),token[2].tobytes(),token[3],None if raw is None else (raw._data,raw._integrity)))
    return value+(tuple(items),)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    args=parser.parse_args()
    from mpi4py import MPI
    world=MPI.COMM_WORLD
    if world.size!=6:raise ValueError('launch exactly six ranks')
    plan=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('ld',(1,),point_ids=(10,)),
        SolverGroup('sbe',(2,),point_ids=(20,)),SolverGroup('custom',(3,4,5),point_ids=(30,31))))
    ld,sbe=models();backend=Backend()
    bindings={10:MaterialBinding('ld',Provider('ld',ld),backend),20:MaterialBinding('sbe',Provider('sbe',sbe),backend),
              30:MaterialBinding('custom',Provider('custom'),backend)}
    bindings[31]=bindings[30]
    model=Multiscale1D(4,spacing=5.,vacuum=(300,50),incident=AcosPulse(.001,.5,6.,polarization=(0.,0.,1.)))
    steps=80;dt=.02;trajectories=[];cases=[]
    with ParallelContext.create(world,plan) as context:
        def create():
            f=model.new_state(dt=dt,steps=steps) if world.rank==0 else None
            return CoupledStepper.from_multiscale(context,bindings,f)
        def close(c):
            c.response.close()
            if c.field is not None:c.field.close()
        for protected in (False,True):
            c=create();local=targets(c);frames=[]
            try:
                if protected:
                    before=signature(c)
                    try:c.step_protected(local[:-1] if world.rank==5 else local)
                    except ValueError:pass
                    else:raise AssertionError('missing owner not rejected')
                    assert signature(c)==before
                    # Reject incompatible optional collective cleanup before capture.
                    capture_calls=[]
                    def audit(fn):
                        def capture():
                            capture_calls.append(True)
                            return fn()
                        return capture
                    mismatched=[StepTarget(t.key,t.owner,audit(t.capture),t.rollback,t.release) for t in local]
                    if world.rank==5:
                        t=mismatched[-1]
                        mismatched[-1]=StepTarget(t.key,t.owner,t.capture,t.rollback)
                    try:c.step_protected(mismatched)
                    except StepRecoveryError as exc:assert exc.phase=='admission' and exc.original_preserved
                    else:raise AssertionError('release policy mismatch not rejected')
                    assert not capture_calls
                    assert signature(c)==before
                    # Partial capture on one rank: release(None) still joins collectives.
                    if world.rank==5:c.response.states[31].fail_capture=True
                    try:c.step_protected(local)
                    except StepRecoveryError as exc:assert exc.original_preserved and not exc.cleanup_errors
                    else:raise AssertionError('capture error not returned')
                    if world.rank==5:c.response.states[31].fail_capture=False
                    assert signature(c)==before
                    cases.extend(['missing-owner','release-policy','partial-capture'])
                for step in range(1,steps+1):
                    if protected and step in (1,25,40,50,80):
                        before=signature(c)
                        if step==40 and world.rank==1:c.response.states[10].fail_next=True
                        elif step==50 and world.rank==0:
                            original=c.field.advance
                            def fail_field(J):
                                result=original(J)
                                c.field.advance=original
                                raise ValueError('injected after field mutation')
                            c.field.advance=fail_field
                        elif step not in (40,50) and world.rank==5:c.response.states[31].fail_next=True
                        try:c.step_protected(local)
                        except StepRecoveryError as exc:assert exc.original_preserved and not exc.cleanup_errors
                        else:raise AssertionError('late error not returned')
                        assert signature(c)==before
                        cases.append('late-'+str(step))
                    frame=c.step_protected(local) if protected else c.step()
                    if frame is not None:frames.append(frame)
                trajectories.append(frames)
            finally:close(c)
        # A rollback failure must poison every group, even if a callback did some work.
        for failure in ('rollback','release'):
            c=create();local=targets(c)
            try:
                if world.rank==5:
                    if failure=='rollback':c.response.states[31].fail_next=True;c.response.states[31].fail_rollback=True
                    else:c.response.states[31].fail_release=True
                try:c.step_protected(local)
                except StepRecoveryError as exc:
                    assert not exc.original_preserved
                    assert c.response.failed
                else:raise AssertionError('unsafe failure not reported')
                cases.append('failed-'+failure)
            finally:close(c)
    if world.rank==0:
        for before,after in zip(*trajectories):
            for key in before:np.testing.assert_array_equal(before[key],after[key],err_msg=key)
        out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
        np.savez(out/'response.npz',A=np.asarray([f['A_over_c'] for f in trajectories[1]]),j=np.asarray([f['J_charge'] for f in trajectories[1]]))
        (out/'result.json').write_text(json.dumps(dict(status='passed',ranks=6,steps=steps,cases=cases,
            exact_to_unprotected=True,replayed_steps=0,scope='Multiscale1D + real LD/SBE handles + explicit custom model; not native aggregate recovery'),indent=2)+'\n')
        print('protected mixed step MPI PASS',cases)


if __name__=='__main__':main()
