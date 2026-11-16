"""Mixed LD/SBE 1D feedback, optionally adding a managed TDDFT point.

mpiexec -n 3 python -m salmonpy.examples.mpi_mixed_multiscale NEW_OUTPUT
Use 4 ranks and --runtime RUNTIME --bundle GS_BUNDLE [--input rt.inp] for TDDFT.
Synthetic LD/SBE parameters: a connection test, not a composite-material prediction.
This endpoint-LD demonstration is not native semi-implicit E/H LD coupling.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from salmonpy.ld import LorentzDrude
from salmonpy.sbe import SBEInitialState
from salmonpy.materials import bind_material
from salmonpy.material_response import (MaterialBinding, LDMaterialProvider, LDResponseBackend,
                                       TDDFTMaterialProvider, TDDFTResponseBackend)
from salmonpy.distributed import ParallelPlan, ParallelContext, SolverGroup
from salmonpy.coupled import CoupledStepper
from salmonpy.maxwell import Multiscale1D
from salmonpy.fields import AcosPulse


def models():
    ld = LorentzDrude(plasma_frequency=.5, strength=[1.], resonance=[.3], damping=[.02])
    momentum = np.zeros((2, 2, 3, 1), complex)
    momentum[0, 1, 2, 0] = momentum[1, 0, 2, 0] = .3
    sbe = SBEInitialState(eigen=[[-.5], [.4]], momentum=momentum,
                         nonlocal_velocity=np.zeros_like(momentum), weights=[1.],
                         lattice=np.eye(3)*4., electrons=2)
    return ld, sbe


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--runtime', type=Path)
    parser.add_argument('--bundle', type=Path)
    parser.add_argument('--input', default='rt.inp')
    parser.add_argument('--steps', type=int, default=800)
    args = parser.parse_args()
    if bool(args.runtime) != bool(args.bundle): parser.error('--runtime and --bundle are required together')
    if args.steps < 1: parser.error('positive steps required')
    from mpi4py import MPI
    world = MPI.COMM_WORLD
    count = 3 if args.runtime else 2
    if world.size != count + 1: raise ValueError(f'launch exactly {count+1} MPI ranks')
    ld, sbe = models()
    bindings = {10: MaterialBinding('ld-layer', LDMaterialProvider(ld), LDResponseBackend()),
                20: bind_material(sbe, key='sbe-layer', identity='synthetic-two-band-v1')}
    provider = None
    if args.runtime:
        provider = TDDFTMaterialProvider(args.bundle, args.input, runtime=args.runtime,
                                         max_steps=args.steps, threads=1)
        bindings[30] = MaterialBinding('tddft-layer', provider, TDDFTResponseBackend())
    points = tuple(bindings)
    plan = ParallelPlan((SolverGroup('field', (0,), role='field'),) + tuple(
        SolverGroup(bindings[p].key, (i+1,), point_ids=(p,)) for i, p in enumerate(points)))
    dt = .02
    field_model = Multiscale1D(count, spacing=5., vacuum=(300, 50),
                              incident=AcosPulse(.001, .5, 6., polarization=(0., 0., 1.)))
    frames = []
    try:
        with ParallelContext.create(world, plan) as context:
            field = field_model.new_state(dt=dt, steps=args.steps) if world.rank == 0 else None
            coupled = CoupledStepper.from_multiscale(context, bindings, field)
            try:
                for _ in range(args.steps):
                    frame = coupled.step()
                    if frame is not None: frames.append(frame)
                communication = coupled.response.communication_report()
            finally:
                coupled.response.close()
                if coupled.field is not None: coupled.field.close()
        # All MPI collectives finish before root-only verification and output.
        if world.rank == 0:
            # Independent serial schedule using the original LD/SBE state calls.
            reference_field = field_model.new_state(dt=dt, steps=args.steps)
            left = right = electronic = None
            try:
                left = ld.new_state(dt=dt)
                right = sbe.new_state(dt=dt)
                if provider is not None: electronic = provider.prepare(None).initial.new_state(dt=dt)
                previous_A = np.zeros((count, 3)); current = np.zeros((count, 3))
                for frame in frames:
                    A = reference_field.step(current)
                    current[0] = left.step(-(A[0]-previous_A[0])/dt)[1]
                    current[1] = right.step(A[1])
                    if electronic is not None: current[2] = electronic.step(A[2])
                    np.testing.assert_array_equal(A, frame['A_over_c'])
                    np.testing.assert_array_equal(current, frame['J_charge'])
                    previous_A = A.copy()
            finally:
                if left is not None: left.close()
                if right is not None: right.close()
                if electronic is not None: electronic.close()
            peaks = np.max(np.abs([f['J_charge'] for f in frames]), axis=(0, 2))
            if not np.all(peaks > 0):
                raise RuntimeError('no response at one or more points; check pulse direction/duration and increase steps')
            output = args.output.resolve(); output.mkdir(parents=True, exist_ok=False)
            np.savez(output/'response.npz', point_ids=points,
                     time=np.arange(1, args.steps+1)*dt,
                     A_over_c=np.asarray([f['A_over_c'] for f in frames]),
                     j=np.asarray([f['J_charge'] for f in frames]))
            report = {'status': 'passed', 'serial_mpi_exact': True,
                      'materials': {str(p): bindings[p].key for p in points},
                      'steps': args.steps, 'dt': dt, 'units': 'atomic',
                      'maximum_abs_current': peaks.tolist(),
                      'communication': communication,
                      'scope': 'endpoint LD plus synthetic SBE; managed TDDFT only when explicitly supplied'}
            (output/'result.json').write_text(json.dumps(report, indent=2)+'\n')
            print(output)
    finally:
        if provider is not None: provider.close()


if __name__ == '__main__': main()
