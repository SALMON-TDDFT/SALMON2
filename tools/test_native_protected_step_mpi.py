"""Run ordinary mixed-native --protect-steps with collective late/deep failures.

Pass the ordinary sample arguments. Optional --kernel LIB for source-only tests.
Use a launcher timeout and OMP/BLAS1. No propagation or numerical kernels mocked.
"""
import argparse
import ctypes as c
from pathlib import Path
import json
import sys
import numpy as np
from mpi4py import MPI


def main():
    parser=argparse.ArgumentParser(add_help=False);parser.add_argument('--kernel');parser.add_argument('--multiscale',action='store_true')
    args,sample_args=parser.parse_known_args()
    if args.kernel:
        from salmonpy import ld,sbe
        from salmonpy._material_handles import LDAPI as API
        from salmonpy._material_handles import SBEAPI
        a=API(args.kernel);b=SBEAPI(args.kernel)
        ld._handles_api=lambda:a;ld._native=lambda:a.legacy
        sbe._handles_api=lambda:b;sbe._native=lambda:b.legacy
    from salmonpy.coupled import CoupledStepper
    from salmonpy.stepping import guarded_step,StepRecoveryError
    if args.multiscale:
        from salmonpy.examples import mpi_multiscale_tddft as sample
    else:
        from salmonpy.examples import mpi_mixed_native as sample
    original=CoupledStepper.step_protected
    visited=set();records=[]
    out=Path(sample_args[3]).resolve()

    def freeze(coupled,label):
        r=coupled.response
        result={'step':coupled.step_index,'J':(id(coupled.J),coupled.J.tobytes()),
            'response':(r.step_index,r.failed,id(r._A),r._A.tobytes()),
            'adapters':tuple((p,id(a.A),a.A.tobytes()) for p,a in r._adapters.items())}
        owners=[]
        if coupled.field is not None:
            if hasattr(coupled.field,'lib'):owners.append((coupled.field,'sp_eh_checkpoint'))
            else:
                t=coupled.field.state.snapshot()
                result['light-field']=(t._integrity,t._data,coupled.field.closed,coupled.field.failed)
        for p,s in r.states.items():
            if hasattr(s,'native'):owners.append((s,'sp_tddft_checkpoint'))
            elif getattr(s,'checkpoint_protocol',None)=='salmonpy.direct-state/1':owners.append((s,'sp_fast_hse_checkpoint'))
            else:
                t=s._raw.snapshot();result[str(p)]=(s._steps,s._failed,s._closed,t._integrity,t._data)
        for i,(s,symbol) in enumerate(owners):
            cp=getattr(s.lib,symbol);cp.argtypes=[c.c_char_p,c.c_int];cp.restype=c.c_int
            path=out/f'{label}-{MPI.COMM_WORLD.rank}-{i}.bin'
            assert s._call(cp,str(path).encode(),0)==0
            result[symbol]=path.read_bytes()
            folder=s.directory if hasattr(s,'directory') else s.work
            result[symbol+'-files']={p.name:p.read_bytes() for p in folder.glob('*.data')}
            names=('A','E','step_index','failed') if hasattr(s,'directory') else ('A','steps','failed')
            for name in names:
                v=getattr(s,name)
                result[symbol+name]=(id(v),v.tobytes()) if isinstance(v,np.ndarray) else v
            if hasattr(s,'current'):result[symbol+'current']=(id(s.current),s.current.tobytes(),s.energy)
        return result

    def exercise(coupled,targets):
        step=coupled.step_index+1
        r=coupled.response
        if step in (1,13,32) and step not in visited:
            visited.add(step)
            # Flush initial headers/accepted records without evolving physics.
            guarded_step(lambda:None,targets,context=r.context)
            if step==1:
                freeze(coupled,'admission-before')
                def probes():
                    owners=[]
                    if coupled.field is not None and hasattr(coupled.field,'lib'):owners.append((coupled.field,'sp_eh'))
                    for s in r.states.values():
                        if hasattr(s,'native'):owners.append((s,'sp_tddft'))
                        elif getattr(s,'checkpoint_protocol',None)=='salmonpy.direct-state/1':owners.append((s,'sp_fast_hse'))
                    for i,(owner,prefix) in enumerate(owners):
                        tx=getattr(owner.lib,prefix+'_step_transaction');tx.argtypes=[c.c_int];tx.restype=c.c_int
                        assert owner._call(tx,0)==14, 'nested native capture accepted'
                        if r.context.solver_comm.size>1:
                            mode=0 if r.context.solver_comm.rank==0 else 1
                            assert owner._call(tx,mode)==11, 'mixed native transaction modes accepted'
                        cp=getattr(owner.lib,prefix+'_checkpoint');cp.argtypes=[c.c_char_p,c.c_int];cp.restype=c.c_int
                        path=out/f'admission-before-{MPI.COMM_WORLD.rank}-{i}.bin'
                        assert owner._call(cp,str(path).encode(),2)!=0, 'checkpoint staged during native transaction'
                guarded_step(probes,targets,context=r.context)
                records.append((0,'transaction-checkpoint-admission'))
            local_sbe=any(not hasattr(s,'native') and getattr(s,'checkpoint_protocol',None)!='salmonpy.direct-state/1' for s in r.states.values())
            flags=MPI.COMM_WORLD.allgather((local_sbe,
                coupled.field is not None and hasattr(coupled.field,'lib'),
                any(hasattr(s,'native') or getattr(s,'checkpoint_protocol',None)=='salmonpy.direct-state/1'
                    for s in r.states.values()),
                any(getattr(s,'symmetry_enabled',False) for s in r.states.values())))
            has_sbe,has_eh,has_native,has_symmetry=(any(row[i] for row in flags) for i in range(4))
            cases=['field-late']+(['sbe-late'] if has_sbe else [])
            if has_native:cases+=['electronic-late','electronic-deep']
            if has_symmetry:cases.append('electronic-field')
            if step==13 and has_eh:cases.append('field-deep')
            for case in cases:
                before=freeze(coupled,f'{step}-{case}-before');replacements=[]
                if case.startswith('field') and coupled.field is not None:
                    s=coupled.field
                    name='advance_owned' if getattr(s,'current_rows',None) is not None else 'advance'
                    fn=getattr(s,name)
                    def call(J,s=s,fn=fn,case=case):
                        if case=='field-deep':
                            s.A[:]=np.finfo(float).max
                            return fn(np.full_like(J,np.finfo(float).max))
                        value=fn(J)
                        if r.context.solver_comm.rank==r.context.solver_comm.size-1:raise ValueError('late field failure')
                        return value
                    replacements.append((s,name,fn));setattr(s,name,call)
                for p,s in r.states.items():
                    native=hasattr(s,'native') or getattr(s,'checkpoint_protocol',None)=='salmonpy.direct-state/1'
                    if case.startswith('electronic') and native:
                        fn=s.advance
                        def call(old,new,fn=fn,case=case,s=s):
                            overflow = [0.,0.,1e100] if getattr(s,'symmetry_enabled',False) else [1e100,0.,0.]
                            drive = np.array(overflow) if case=='electronic-deep' else new
                            if case=='electronic-field' and getattr(s,'symmetry_enabled',False):
                                drive = np.array([1e-250,0.,0.])
                            value=fn(old,drive)
                            if case=='electronic-late' and r.context.solver_comm.rank==r.context.solver_comm.size-1:
                                raise ValueError('late electronic failure')
                            return value
                        replacements.append((s,'advance',fn));s.advance=call
                    if case=='sbe-late' and not native:
                        fn=s.step
                        def call(A,fn=fn):
                            value=fn(A);raise ValueError('late SBE failure')
                        replacements.append((s,'step',fn));s.step=call
                try:
                    try:original(coupled,targets)
                    except StepRecoveryError as exc:
                        assert exc.original_preserved and not exc.cleanup_errors,(exc,exc.rollback_errors,exc.cleanup_errors)
                    else:raise AssertionError('injected failure not returned')
                finally:
                    for owner,name,fn in replacements:setattr(owner,name,fn)
                after=freeze(coupled,f'{step}-{case}-after')
                assert after==before,(step,case,MPI.COMM_WORLD.rank,[k for k in before if before[k]!=after[k]])
                records.append((step,case))
        return original(coupled,targets)

    CoupledStepper.step_protected=exercise
    sys.argv=[sys.argv[0]]+sample_args+([] if '--protect-steps' in sample_args else ['--protect-steps'])
    sample.main()
    if MPI.COMM_WORLD.rank==0:
        (out/'protected-step-validation.json').write_text(json.dumps(dict(status='passed',cases=records,
            all_original_bytes_and_files_preserved=True,restore_replayed_steps=0),indent=2)+'\n')
        print('native protected mixed step PASS',records)


if __name__=='__main__':main()
