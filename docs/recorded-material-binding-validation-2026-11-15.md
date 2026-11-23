# Recorded material binding validation (2026-11-15)

M5e subset: bind_recorded_material derives stable adapter identity from a model's
record(), and response_identity constructs point-run/1 from recorded bindings in
integer transport order. MaterialBinding optionally stores a validated detached
immutable record. ResponseSession checks the supplied point record against it.
Record drift is rejected both before prepare and before each state creation.
Known ConfiguredMaterial wrappers cannot bypass the managed TDDFT restriction.

Scope is endpoint-A serial models (configured SBE and user models), one rank per
material group. No LD E conversion, native collective TDDFT factory or restart
restoration was added. Prepared models remain borrowed; existing numerical adapter
sign/timing/clock checks and state cleanup remain in use. Declared hashes are not
independent certification of custom hidden storage or physical accuracy.

Review also exposed inherited empty-message callback handling: str(RuntimeError())
was falsy and could escape collective error agreement during state creation.
ResponseSession callback exception diagnostics now include the nonempty exception
class name. Three-rank fault tests verify whole-session rejection for empty-message
create and step failures, with a 30-second launcher bound per case.

Evidence (Apple Silicon/macOS, Python 3.14):

- Six focused tests: copied binding records, automatic ordered metadata, invalid
  mapping, record drift before prepare/create, configured SBE/direct parity and
  wrapped TDDFT rejection.
- Source suite: 378 tests, one case-sensitive-filesystem skip.
- Installed wheel suite: 378 tests, two skips (same condition and optional plot).
- sdist -> wheel build with max two jobs; fresh venv, isolated package imports.
- Actual installed-wheel MPI, three ranks: custom recorded model and configured
  two-band/one-k-point synthetic SBE both exactly match direct serial currents at
  all eight steps; saved material/point/run identities read back unchanged.
- Empty-message create and step failures rejected on all three ranks, without hang.

Installed test initially found a source-layout-dependent fixture import. Changed
the test to import its sibling in both package and unittest-discovery modes, then
reran the installed suite. Package code and SBE numerical behavior were unaffected.

Local MPI artifacts: `/tmp/salmonpy-recorded-installed-mpi-20261115`.
Synthetic SBE inputs: `/tmp/salmonpy-recorded-sbe-inputs-20261115`.
Build/test logs: `/tmp/salmonpy-recorded-*.log`.
Reproduce user paths with examples/mpi_recorded_material.py as documented in
RECORDED_MATERIALS.md/en; fault worker is tools/test_recorded_material_mpi.py.

No full native runtime build, GitHub Actions, commit or push. One thread/rank,
three MPI ranks. Peak memory not measured; no hard budget enforcement. Synthetic
SBE/custom interface parity is not real-material convergence or platform coverage.
Remaining: dedicated configured TDDFT/LD provider factory, collective native
backend capability integration and restart identity/compatibility.
