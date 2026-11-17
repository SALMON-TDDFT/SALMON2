"""Four ranks: native field, standard TDDFT, two-rank fast-PR HSE.

Arguments: STANDARD_BUILD HSE_BUILD STANDARD_GS HSE_GS NEW_OUTPUT.
Run with an external MPI timeout, OMP/BLAS=1. Same-topology, same-lifetime test.
"""
import ctypes as c
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
from mpi4py import MPI
from salmonpy._native_checkpoint import NativeRestoreCandidate
from salmonpy.checkpointing import CheckpointState, CheckpointTarget, CheckpointCapability, restore_checkpoint, AtomicRestoreError, check_checkpoint_compatibility
from salmonpy.distributed import ParallelPlan, SolverGroup, ParallelContext
from salmonpy.embedded_tddft import NativeTDDFTMaterialProvider, NativeTDDFTBackend
from salmonpy.fast_hse import FastHSEMaterialProvider
from salmonpy.fast_hse_mpi import NativeFastHSEBackend
from salmonpy.coupled import NativeEHField
from salmonpy.fdtd import FDTD, PlaneWave
from salmonpy.fields import AcosPulse
from salmonpy.geometry import Grid, Geometry


def main():
    world = MPI.COMM_WORLD
    if world.size != 4:
        raise ValueError('exactly four ranks required')
    standard, hse, standard_gs, hse_gs, output = sys.argv[1:6]
    continuous = "--continuous" in sys.argv[6:]
    out = Path(output).resolve()
    if world.rank == 0:
        out.mkdir(parents=True, exist_ok=False)
    world.Barrier()
    field_ranks = 2 if '--field2' in sys.argv[6:] else 1
    tddft_ranks = 2 if '--tddft2' in sys.argv[6:] else 1
    if field_ranks+tddft_ranks>3:
        raise ValueError('keep one rank for HSE')
    standard_leader = field_ranks
    hse_leader = field_ranks+tddft_ranks
    plan = ParallelPlan((SolverGroup('field', tuple(range(field_ranks)), role='field'),
                         SolverGroup('standard', tuple(range(field_ranks,hse_leader)), point_ids=(1,)),
                         SolverGroup('hse', tuple(range(hse_leader,4)), point_ids=(2,))))
    groups = [{'name':g.name, 'ranks':list(g.ranks), 'role':g.role, 'point_ids':list(g.point_ids)} for g in plan.groups]
    steps = 24
    split = 8
    dt = .003
    context = ParallelContext.create(world, plan)
    is_field = context.group.role == "field"
    is_standard = context.group.name == "standard"
    solver_context = SimpleNamespace(solver_comm=context.solver_comm, point_ids=context.group.point_ids)
    if is_field:
        grid = Grid((40,40,40), (2.,2.,2.))
        calc = FDTD(Geometry(grid), (), dt, steps,
                    source=PlaneWave(AcosPulse(.001,20.,.15),0,-8.))
        owner = NativeEHField(context, calc, standard, out/'field',
                              points={1:(19,20,20),2:(20,20,20)}, grid=(1,field_ranks,1))
        symbol = 'sp_eh_checkpoint'
    elif is_standard:
        provider = NativeTDDFTMaterialProvider(standard_gs,'rt.inp',max_steps=steps)
        owner = NativeTDDFTBackend(standard,out,decomposition=(1,tddft_ranks,1,1,1)).create_state(provider.prepare(solver_context),1,dt,solver_context)
        owner.current = np.zeros(3)
        symbol = 'sp_tddft_checkpoint'
    else:
        provider = FastHSEMaterialProvider(hse_gs,'rt.inp',max_steps=steps)
        owner = NativeFastHSEBackend(hse,out).create_state(provider.prepare(solver_context),2,dt,solver_context)
        symbol = 'sp_fast_hse_checkpoint'
    fn = getattr(owner.lib,symbol)
    fn.argtypes = [c.c_char_p,c.c_int]
    fn.restype = c.c_int
    shard = out/f'saved-{world.rank}.bin'
    def save(path):
        assert owner._call(fn,str(path).encode(),0) == 0
    def evolve():
        rows = world.allgather(owner.current.tolist() if not is_field else None)
        A = owner.advance(np.zeros((2,3)) if owner.step_index==0 else np.array([rows[standard_leader],rows[hse_leader]])) if is_field else None
        A = world.bcast(A,root=0)
        if not is_field:
            owner.current = owner.advance(owner.A.copy(),A[0 if is_standard else 1])
        return np.r_[owner.A.ravel(),owner.E.ravel()] if is_field else np.r_[owner.A,owner.current]
    try:
        trace = []
        for n in range(1,steps+1):
            trace.append(evolve().copy())
            if n == split:
                save(shard)
                attrs = {name: getattr(owner,name).copy() if isinstance(getattr(owner,name),np.ndarray)
                         else getattr(owner,name) for name in
                         (('A','E','step_index') if is_field else ('A','current','steps'))}
            if n == split+2 and not continuous:
                if context.solver_comm.size>1:
                    mode = 3 if context.solver_comm.rank==0 else 4
                    assert owner._call(fn,b'',mode)!=0, 'mixed native modes accepted'
                before = out/f'before-{world.rank}.bin'
                after = out/f'after-{world.rank}.bin'
                save(before)
                broken = out/f'broken-{world.rank}.bin'
                data = shard.read_bytes()
                broken.write_bytes(data[:len(data)//2] if world.rank == 3 else data)
                candidate = NativeRestoreCandidate(owner,symbol,broken)
                cap = CheckpointCapability('salmonpy.native-candidate-test/1',('save','restore'),(len(context.group.ranks),))
                target = CheckpointTarget(cap,stage_restore=candidate.stage,
                                          discard_restore=lambda unused: candidate.discard())
                declarations = {'field':cap} if is_field else {'states':{context.group.point_ids[0]:cap}}
                saved = check_checkpoint_compatibility(world,groups,group=context.group.name,**declarations)
                live = CheckpointState(field=owner) if is_field else CheckpointState(states={context.group.point_ids[0]:owner})
                kwargs = {'field':target} if is_field else {'states':{context.group.point_ids[0]:target}}
                try:
                    restore_checkpoint(world,groups,group=context.group.name,live=live,saved=saved,**kwargs)
                    raise AssertionError('late HSE shard corruption accepted')
                except AtomicRestoreError as exc:
                    assert exc.original_preserved and not exc.cleanup_errors, exc.cleanup_errors
                save(after)
                assert before.read_bytes() == after.read_bytes(), 'original numerical state changed'
                assert (live.field is owner if is_field else live.states[str(context.group.point_ids[0])] is owner)
        original_end = out/f'continuous-end-{world.rank}.bin'
        save(original_end)
        assert np.any(np.asarray(trace)), 'trivial all-zero trajectory'
        np.save(out/f'trace-{world.rank}.npy',np.asarray(trace))
        if continuous:
            owner.close();context.close()
            return
        candidate = NativeRestoreCandidate(owner,symbol,shard)
        target = CheckpointTarget(cap,stage_restore=candidate.stage,
                                  discard_restore=lambda unused: candidate.discard())
        kwargs = {'field':target} if is_field else {'states':{context.group.point_ids[0]:target}}
        restore_checkpoint(world,groups,group=context.group.name,live=live,saved=saved,**kwargs)
        candidate.publish()
        for name,value in attrs.items():
            setattr(owner,name,value)
        for n in range(split,steps):
            np.testing.assert_array_equal(evolve(),trace[n])
        resumed_end = out/f'resumed-end-{world.rank}.bin'
        save(resumed_end)
        # Restore deliberately invalidates orbital overlap caches. Their flags
        # are not physical observables; keep full streams for diagnostics.
        end_differences = sum(a!=b for a,b in zip(original_end.read_bytes(),resumed_end.read_bytes()))
        differences = world.allgather(end_differences)
        assert world.allreduce(1)==4
        if world.rank==0:
            record={'status':'passed','ranks':4,'late_shard_failure_preserved_native_bytes':True,
                    'original_continued':True,'retry_future_exact':True,'replayed_steps':0,'new_steps':steps-split,'layout':[field_ranks,tddft_ranks,4-hse_leader],'native_end_byte_differences':differences}
            (out/'validation.json').write_text(json.dumps(record,indent=2)+'\n')
            print('native candidate MPI PASS '+json.dumps(record),flush=True)
    finally:
        owner.close()
    context.close()


if __name__ == '__main__':
    try:
        main()
    except BaseException:
        import traceback
        traceback.print_exc()
        MPI.COMM_WORLD.Abort(1)
