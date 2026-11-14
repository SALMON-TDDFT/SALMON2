"""One field rank + two independent TDDFT groups, compared to native SALMON.

mpiexec -n 3|5 python tools/test_ms1d_mpi.py BUILD_JSON REFERENCE_DIR NEW_OUTPUT
REFERENCE_DIR is test_longitudinal_salmon.py --trans-longi ms1d --cells 2 output.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy.coupled import CoupledStepper
from salmonpy.distributed import ParallelContext, ParallelPlan, SolverGroup
from salmonpy.embedded_tddft import NativeTDDFTBackend, NativeTDDFTMaterialProvider
from salmonpy.material_response import MaterialBinding
from salmonpy.maxwell import Multiscale1D
from salmonpy.fields import AcosPulse
from test_longitudinal_salmon import DT, NT, LAYER, PULSE, compare_ms1d_reference


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('build'); p.add_argument('reference'); p.add_argument('output')
    args = p.parse_args()
    world = MPI.COMM_WORLD
    if world.size not in (3, 5): raise ValueError('3 or 5 ranks required')
    size = (world.size-1)//2
    reference = Path(args.reference).resolve(); out = Path(args.output).resolve()
    if world.rank == 0: out.mkdir(parents=True, exist_ok=False)
    world.Barrier()
    plan = ParallelPlan((SolverGroup('field', (0,), role='field'),
        SolverGroup('left', tuple(range(1, size+1)), point_ids=(19,)),
        SolverGroup('right', tuple(range(size+1, 2*size+1)), point_ids=(4,))))
    # The immutable input is shared; each group creates its own mutable state.
    provider = NativeTDDFTMaterialProvider(reference/'salmonpy', 'rt.inp', max_steps=NT)
    backend = NativeTDDFTBackend(args.build, out, decomposition=(1, size, 1, 1, 1))
    binding = MaterialBinding('same-GS', provider, backend)
    with ParallelContext.create(world, plan) as context:
        field = Multiscale1D(2, **LAYER, incident=AcosPulse(**PULSE)).new_state(dt=DT, steps=NT) if world.rank == 0 else None
        coupled = CoupledStepper.from_multiscale(context, {19: binding, 4: binding}, field)
        try:
            if field is not None:
                grids = [field.grid_A.copy()]; currents = [np.asarray([coupled.response.initial_current[p] for p in plan.point_ids])]; sites = [field.A.copy()]
            for n in range(NT):
                frame = coupled.step()
                if frame is not None:
                    assert frame['electronic_field_clock'] == (n+1)*DT
                    assert coupled.timing['electronic_clock_offset'] == '0'
                    grids.append(field.grid_A.copy()); currents.append(frame['J_charge']); sites.append(frame['A_over_c'])
            report = coupled.response.communication_report()
        finally:
            coupled.response.close()
            if coupled.field is not None: coupled.field.close()
    assert world.allreduce(1) == world.size, 'outer communicator must remain usable'
    if world.rank == 0:
        grids = np.asarray(grids); currents = np.asarray(currents); sites = np.asarray(sites)
        with np.load(reference/'reference.npz') as ref:
            comparison = compare_ms1d_reference(ref, grids, currents, tolerance=1e-12)
        with np.load(reference/'python.npz') as serial:
            np.testing.assert_allclose(grids, serial['grid_A'], rtol=2e-12, atol=2e-16)
            np.testing.assert_allclose(currents[1:], serial['current'][1:], rtol=2e-10, atol=2e-14)
        assert comparison['passed'], comparison
        assert np.max(np.abs(currents[:, 0]-currents[:, 1])) > 1e-7
        np.testing.assert_array_equal(sites, grids[:, 300:302])
        np.savez(out/'history.npz', grid_A=grids, current=currents, A=sites)
        (out/'comparison.json').write_text(json.dumps(comparison, indent=2)+'\n')
        (out/'communication.json').write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps(comparison), flush=True)


if __name__ == '__main__':
    try: main()
    except BaseException:
        import traceback
        traceback.print_exc()
        MPI.COMM_WORLD.Abort(1)
