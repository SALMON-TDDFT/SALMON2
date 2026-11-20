"""Offline job admission only; no calculation is launched.

mpiexec -n 2 python -m salmonpy.examples.mpi_input_preflight INPUT_BUNDLE NEW_REPORT.json
Every rank checks its local view; errors are gathered before collective exit.
Call the same pattern before creating solver states in your batch script.
"""
import argparse
import json
from pathlib import Path

from salmonpy.preparation import preflight_inputs


def check_inputs(comm, directory, report, *, manifest_sha256=None, registry=None):
    """Gather input-only diagnostics before any solver state is created."""
    directory, report = Path(directory), Path(report)
    local = preflight_inputs(directory, manifest_sha256=manifest_sha256,
                             registry=registry)
    reports = comm.gather({'rank': comm.rank, **local}, root=0)
    outcome = None
    if comm.rank == 0:
        outcome = {'schema': 'salmonpy.job-input-preflight/1',
                   'status': 'passed' if all(r['verified'] for r in reports) else 'failed',
                   'ranks': reports, 'solver_started': False}
        try:
            if report.resolve().is_relative_to(directory.resolve()):
                raise ValueError('report must be outside the input bundle')
            with report.open('x', encoding='utf-8') as stream:
                json.dump(outcome, stream, ensure_ascii=False, indent=2)
                stream.write('\n')
        except Exception as exc:
            outcome['status'] = 'failed'
            outcome['report_error'] = str(exc)
        print(json.dumps(outcome, ensure_ascii=False))
    status = comm.bcast(outcome['status'] if outcome else None, root=0)
    return status == 'passed'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('report', type=Path)
    parser.add_argument('--manifest-sha256')
    parser.add_argument('--resolve-builtins', action='store_true')
    args = parser.parse_args()
    from mpi4py import MPI
    registry = None
    if args.resolve_builtins:
        from salmonpy.builtin_materials import builtin_material_registry
        registry = builtin_material_registry()
    passed = check_inputs(MPI.COMM_WORLD, args.directory, args.report,
                          manifest_sha256=args.manifest_sha256, registry=registry)
    raise SystemExit(0 if passed else 2)


if __name__ == '__main__':
    main()
