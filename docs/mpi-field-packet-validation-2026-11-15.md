# MPI field packet validation — 2026-11-15

ResponseSession accepts explicit per-point field_contracts and step_fields. The
legacy step(A) path remains separate. Contracts agree before providers run;
all points are validated before callbacks and receivers reconstruct immutable
packets. FieldAdapter checks state input timing; session checks endpoint charge-j
output declarations. Callback failure invalidates the session, not global rollback.
No interpolation or physical conversion is inferred.

Apple Silicon macOS/Python 3.14/OpenMPI, three ranks, one thread each, 45-second
launcher bounds; no Actions:

- Source and final installed wheel: distinct A and E/H material groups yield
  expected synthetic currents; required-port/time mismatch refuses every update.
- timing_records correctly reports A versus E/H declarations.
- Injected E/H callback error propagates to all ranks without hanging.
- Backend J_matter declaration contradicting charge j rejects collectively at
  creation; no implicit sign reversal.
- Installed legacy custom A-model serial/MPI comparison passes.
- Installed full tests: 411 tests, two skips (optional plot and case-sensitive FS).
- sdist→wheel build and installation use explicit native kernels. -I keeps source
  checkout imports out of installed MPI validation.

Read-only review identified wrong timing-record fallback and missing output
contract validation; both were corrected and tested. This is transport validation,
not native molecular/FDTD physical validation or production convergence.

Artifacts: /tmp/salmonpy-mpi-fields-final-dist-20261115/.
Tests: tools/test_mpi_field_packets.py success|callback-failure|output-mismatch,
with mpiexec -n 3. Logs /tmp/salmonpy-mpi-fields-final-installed-*.log.
Human guides, AI JSON and mpi_field_inputs sample document explicit limitations.

Final integration guard: CoupledStepper refuses packet sessions before field
advance until explicit sampling/clock integration exists. A unit regression
checks unchanged field step index. Full source suite after the guard: 412 tests,
one skip. The earlier installed-wheel count remains 411 before this extra test;
the new guard was validated in source.
