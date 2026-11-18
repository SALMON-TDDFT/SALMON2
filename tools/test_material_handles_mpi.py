"""Shared LD/SBE point-MPI parity, used-model broadcast and ownership checks.

mpiexec -n 3 python tools/test_material_handles_mpi.py LIBRARY NEW_JSON --model ld
Existing libraries only: never builds kernels or generates GS data.
"""
import argparse
import json
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("library")
    parser.add_argument("output", type=Path)
    parser.add_argument("--model", choices=("ld", "sbe"), required=True)
    args = parser.parse_args(argv)
    model = args.model
    import numpy as np
    from mpi4py import MPI
    import salmonpy
    from salmonpy.distributed import ParallelContext
    from salmonpy.material_response import ResponseSession
    from salmonpy.examples.custom_material import plan_for, waveform

    if model == "ld":
        import salmonpy.ld as module
        from salmonpy._ld_handles import API
        from salmonpy.material_response import MaterialBinding, LDMaterialProvider, LDResponseBackend
        state_type = module.LorentzDrude
    else:
        import salmonpy.sbe as module
        from salmonpy._sbe_handles import API
        from salmonpy.materials import bind_material
        from salmonpy.crystal import respond
        state_type = module.SBEInitialState
    api = API(args.library)
    module._handles_api = lambda: api
    module._native = lambda: api.legacy
    world = MPI.COMM_WORLD
    if world.rank == 0:
        if model == "ld":
            initial = module.LorentzDrude(plasma_frequency=1., strength=[.7, .3],
                                         resonance=[0., 1.4], damping=[.1, .2])
            warm_field = [.01, .02, .03]
        else:
            p = np.zeros((2, 2, 3, 1), complex)
            p[0, 1, 2, 0] = p[1, 0, 2, 0] = .3
            initial = module.SBEInitialState([[-.5], [.4]], p, p*.1, [1.], np.eye(3)*4, 2, 1, True)
            warm_field = [0, 0, .01]
        with initial.new_state(dt=.02) as warm:
            warm.step(warm_field)
    else:
        initial = None
    initial = world.bcast(initial, root=0)
    if model == "ld":
        binding = MaterialBinding("ld", LDMaterialProvider(initial), LDResponseBackend())
    else:
        assert not initial.momentum.flags.writeable
        binding = bind_material(initial, key="sbe", identity="synthetic-2band-vnl-order1/1")
    plan = plan_for(world.size)
    finalizers = []
    original = state_type.new_state

    def tracked(self, **kwargs):
        state = original(self, **kwargs)
        finalizers.append(state._native_state.finalizer)
        return state

    history = []
    # Keep LD's explicit cleanup assertion; SBE also checks shared model ownership.
    if model == "ld":
        state_type.new_state = tracked
    try:
        with ParallelContext.create(world, plan) as context:
            session = ResponseSession(context, {p: binding for p in plan.point_ids}, dt=.02)
            try:
                if model == "sbe" and context.group.point_ids:
                    assert len({id(state._raw._native_state.model) for state in session.states.values()}) == 1
                for n in range(1, 9):
                    A = np.array([waveform(p).values[n] for p in plan.point_ids]) if world.rank == 0 else None
                    current = session.advance(A)
                    if world.rank == 0:
                        history.append(current.copy())
            finally:
                session.close()
    finally:
        state_type.new_state = original
    assert all(not f.alive for f in finalizers)
    closed = world.gather(len(finalizers), root=0)
    if world.rank != 0:
        return
    reference = []
    for point in plan.point_ids:
        if model == "ld":
            A = waveform(point).values
            E = np.zeros_like(A)
            E[1:] = -(A[1:]-A[:-1])/.02
            reference.append(initial.respond(E, dt=.02).current[1:])
        else:
            reference.append(respond(initial, waveform(point)).current[1:])
    expected = np.stack(reference, axis=1)
    actual = np.asarray(history)
    assert actual.tobytes() == expected.tobytes()
    summary = {"ranks": world.size, "bitwise_serial_parity": True,
               "maximum_current_difference": float(np.max(abs(actual-expected))),
               "history": actual.tolist(), "library_sha256": api.sha256, "package": salmonpy.__file__}
    if model == "ld":
        summary.update(used_material_broadcast=True, closed_native_states_per_rank=closed)
    else:
        summary.update(used_initial_model_broadcast=True, native_model_shared_per_initial=True)
    with args.output.open("x") as stream:
        json.dump(summary, stream, indent=2, allow_nan=False)
    print(model.upper(), "handles and point-MPI parity PASS", world.size)


if __name__ == "__main__":
    main()
