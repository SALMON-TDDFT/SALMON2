"""Three-rank native whole-step negative: BUILD GS NEW_OUTPUT MODE.

MODE: old-protocol, output-rejection, unsafe. Use OMP/BLAS1 and launcher timeout.
"""
import sys
from pathlib import Path
import json
from mpi4py import MPI
from salmonpy.distributed import ParallelContext,ParallelPlan,SolverGroup
from salmonpy.material_response import MaterialBinding,ResponseSession
from salmonpy.embedded_tddft import NativeTDDFTMaterialProvider,NativeTDDFTBackend
from salmonpy.coupled import NativeEHField,CoupledStepper
from salmonpy.fdtd import FDTD,PlaneWave
from salmonpy.geometry import Grid,Geometry
from salmonpy.fields import AcosPulse
from salmonpy.native_steps import native_step_targets
from salmonpy.stepping import StepRecoveryError
from salmonpy.errors import InputError


def main():
    build,gs,destination,mode=sys.argv[1:5]
    assert mode in ('old-protocol','output-rejection','unsafe')
    hse_build=sys.argv[sys.argv.index('--hse-build')+1] if '--hse-build' in sys.argv else None
    world=MPI.COMM_WORLD;assert world.size==3
    out=Path(destination).resolve()
    if world.rank==0:out.mkdir(parents=True,exist_ok=False)
    world.Barrier()
    plan=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('electronic',(1,2),point_ids=(1,))))
    with ParallelContext.create(world,plan) as context:
        if hse_build:
            from salmonpy.fast_hse import FastHSEMaterialProvider
            from salmonpy.fast_hse_mpi import NativeFastHSEBackend
            provider=FastHSEMaterialProvider(gs,'rt.inp',max_steps=4)
            backend=NativeFastHSEBackend(hse_build,out)
        else:
            provider=NativeTDDFTMaterialProvider(gs,'rt.inp',max_steps=4)
            backend=NativeTDDFTBackend(build,out,decomposition=(1,2,1,1,1))
        r=ResponseSession(context,{1:MaterialBinding('electronic',provider,backend)},dt=.003)
        field=None
        if context.group.role=='field':
            calc=FDTD(Geometry(Grid((40,40,40),(2.,2.,2.))),(),.003,4,
                      source=PlaneWave(AcosPulse(.001,20.,.15),0,-16.))
            field=NativeEHField(context,calc,build,out/'field',points={1:(20,20,20)},grid=(1,1,1))
        coupled=CoupledStepper(r,field)
        try:
            if mode=='old-protocol':
                try:native_step_targets(coupled)
                except InputError:pass
                else:raise AssertionError('old native transaction protocol accepted')
            elif mode=='output-rejection':
                targets=native_step_targets(coupled)
                try:coupled.step_protected(targets)
                except StepRecoveryError as exc:
                    assert exc.phase=='capture' and exc.original_preserved and not exc.cleanup_errors
                else:raise AssertionError('unaudited output folder accepted')
            else:
                targets=native_step_targets(coupled)
                if context.group.role!='field':
                    state=r.states[1];original=state.advance
                    def unknown(old,new):
                        result=original(old,new)
                        if context.solver_comm.rank==1:
                            state.failed=True
                            raise ValueError('synthetic unknown native lifetime failure')
                        return result
                    state.advance=unknown
                try:coupled.step_protected(targets)
                except StepRecoveryError as exc:
                    assert not exc.original_preserved and exc.rollback_errors and exc.cleanup_errors
                else:raise AssertionError('unsafe native lifetime revived')
                assert r.failed
                if field is not None:assert field.failed
            if mode!='unsafe':
                assert coupled.step_index==0 and r.step_index==0 and not r.failed
                if field is not None:assert field.step_index==0 and not field.failed
                for s in r.states.values():assert s.steps==0 and not s.failed
        finally:
            r.close()
            if field is not None:field.close()
    if world.rank==0:
        (out/'admission-validation.json').write_text(json.dumps(dict(status='passed',case=mode),indent=2)+'\n')
        print('native protected admission PASS',mode)


if __name__=='__main__':main()
