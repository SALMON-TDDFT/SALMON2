# MPI point identity validation (2026-11-15)

M5d/M6b: explicit user identity over the existing integer-point response service.
ResponseSession and CoupledStepper.from_multiscale accept optional point_record
and point_map, agreed across ranks before callbacks. Callback SolverContext carries
immutable group-local records. The numerical step and current reduction paths are
unchanged. No state wrapper is added to native collective objects.

Collective save_arrays validates caller roles and point axes, writes on field
leader only, then agrees errors. Records include material/initial/point/run IDs,
definition/data digests and transport/group ownership. open_point_arrays validates
the optional MPI ownership mapping. Identities are caller declarations, not proof
of provider content or a solver restart implementation.

## Evidence

- Source tests: 372, one skip (case-sensitive filesystem requirement).
- Installed wheel tests: 372, two skips (same condition plus optional plotting).
- Five focused unit tests: mapping/type/order, absent/partial metadata, conflicting
  definitions, saved mapping tamper, empty-message writer error propagation.
- Local actual MPI runner: seven cases pass both from source and installed wheel:
  3-rank success, 4-rank success with subgroup sum, cross-group mismatch,
  within-group mismatch, malformed mapping, invalid output axis, genuine existing-
  directory write failure. Launches bounded to 45 seconds each; no deadlocks.
- Successful currents agree exactly with the legacy no-metadata session. Direct
  analytic synthetic response agrees within 1e-14 absolute tolerance. IDs and
  array output read back correctly. Callback metadata is immutable.
- Existing 3-rank custom-material adapter regression passes.
- Existing 3-rank Multiscale1D contract regression passes, now including explicit
  identity metadata versus normal-run array equality. Initial-A/timing rejection,
  backend failure, mixed callbacks and resource-clock failure remain covered.

Review caught and corrected empty-message exceptions being treated as success and
non-plain-string point_map values being accepted through NumPy scalar equality.
No full native TDDFT physical calculation was repeated: these are identity,
communication and synthetic integration tests, not physical-accuracy evidence.

## Environment and reproduction

Apple Silicon/macOS, Python 3.14, mpi4py 4.1.2 and local Open MPI. mpi4py installed
explicitly in a temporary dependency directory; no import-time installation.
MPI needed permission to open local communication sockets outside the sandbox.
Build: sdist -> wheel with at most two jobs. MPI: at most four ranks, one thread
per rank. Peak memory not measured/hard-limited. No GitHub Actions or native build.

```sh
# Installed wheel + mpi4py, outside source directory; fresh output directory.
python -I /path/to/checkout/tools/test_mpi_point_identity.py /new/test/output
# User sample, already prepared MPI environment:
mpiexec -n 3 python -m salmonpy.examples.mpi_point_identity /new/result
```

Installed results: `/tmp/salmonpy-mpi-identity-installed-runs-20261115`.
Source MPI results: `/tmp/salmonpy-mpi-identity-final-source-20261115`.
Existing contract logs: `/tmp/salmonpy-mpi-identity-existing-20261115`.
Installed import path verified inside the dedicated venv; no checkout on PYTHONPATH.
The bounded runner loads only its sibling launcher helper, including under -I.

Remaining: automatic configured-model/provider capability integration, verified
native-source identity rather than declarations alone, MPI restart/checkpoint
restoration, full native state isolation and platform/physical convergence.
No commit/push in this continuation.
