# Optional backend group-size admission — 2026-11-17

Layer planning and direct ResponseSession construction now share the same rank-admission helper. A backend may declare supported_group_sizes as a nonempty tuple of unique positive Python integers. Missing/None declarations preserve the existing execution/reduction checks and do not certify unreported capabilities. Fast-PR HSE declares (1,2).

ResponseSession checks the owning group's rank count and agrees normalized declarations across every group before provider.prepare/create_state. Tuple order is irrelevant. A packetized A binding forwards the live underlying declaration rather than hiding it or caching an outdated value. This is rank admission, not a generic physical-validity, persistence or rollback capability registry.

## Evidence

- Red test reproduced direct ResponseSession calling provider.prepare with an unsupported/malformed declaration. It now fails before callbacks. Layer-planning red regression was recorded separately.
- Three-rank real MPI (field1/custom2): nine invalid/declaration-mismatch cases rejected collectively without provider/create calls or dependent registration. Includes a single invalid rank, within-group mismatch, cross-group mismatch and missing-versus-declared disagreement.
- Different tuple orderings agree. Both declared and legacy-undeclared backends produce four-step manufactured currents exactly equal to the serial expression.
- Installed-wheel MPI regression repeats the same checks.
- Installed-wheel native field1/SBE1/HSE-k2 with local Drude/Lorentz,256 steps at dt=.001: point IDs,A/c,E,j and serial-reference arrays byte-identical to the preceding fast-PR HSE continuous baseline. No native kernel, numerical formula, units or update order changed. Native library and prepared inputs were reused unchanged.
- Source and installed-wheel suites:473 tests each,1 skip. Source distribution to wheel, guide synchronization and diff checks passed.

## Reproduction

Use a matching mpi4py/MPI environment with a launcher timeout:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec --timeout 90 -n 3 python tools/test_response_admission_mpi.py
```

For installed-wheel admission, use python -I with the same script after installing the built wheel. The manufactured MPI test needs no native runtime. Native HSE/FDTD reproduction uses COUPLED.md / COUPLED.en.md and the fixed artifacts recorded in fast-hse-kmpi-validation-2026-11-16.md.

Local persistent evidence: /Users/otobetoshihito/.codex/outputs/salmonpy-response-admission-20261117: mpi-source-final.log,mpi-installed.log,source-tests.log,installed-tests.log,package.log,dist/,hse-coupled-regression/. Environment: Apple Silicon/macOS,Python3.14,OpenMPI. No remote Actions or new platform qualification. Memory caps, scaling, native recoverable failures, generic checkpoint/decomposition compatibility and orbital/rgrid HSE remain separate work. Molecular work remains paused.
