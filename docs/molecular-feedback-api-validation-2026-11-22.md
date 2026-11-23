# Single-molecule public feedback validation — 2026-11-22

`salmonpy.dipole_source.MolecularFeedback` and `MolecularFrame` expose the verified
single-molecule schedule in the existing module. The caller owns MPI, both field
groups, one isolated state, mapping and output. start() holds current zero for
one field interval; step() predicts, rolls back the predictor, injects the
Fortran-mapped dipole current and validates the accepted midpoint field. Numerical
subtraction, averaging, current mapping and energy work remain Fortran operations.
No implicit SCF/build/download/retry or MPI initialization/finalization.

The normal interval retains two gathers and two broadcasts. Capture-side errors
are agreed within the field group before any advance, including leader-only
append-journal failures. Whole-application restoration agreement happens only
on failure. A failed candidate returns StepRecoveryError with original_preserved;
retry is allowed only when preservation was established. Electronic restoration
still replays accepted history. This does not certify arbitrary OS/MPI failure
recovery or direct no-recomputation electronic restart.

Setup checks current settings, leaders, time steps, IDs, physical cell volumes,
grid/boundary/sample/deposition footprint and complete energy-observer participation.
The reference's vacuum/no-other-source condition remains a caller declaration.
The scope is one molecule/point, laboratory-frame dipoles and explicit
policy=environment_scattered_only. This removes vacuum radiation reaction.

## Evidence

Existing authorized C2H2 PZ GS, no new SCF: fixed ions, 64 field steps, dt=.003 au.
Vacuum uses MPI4; isotropic LD uses MPI6 with two field ranks in each context.
The field library and isolated E3 worker are unchanged audited artifacts, reused
after their checksum validation.

- Public API continuous and installed-wheel reject30/retry retain every previous
  saved physical/energy/work array and every six-face record exactly.
- Real MPI leader-only capture-journal failure at interval30: no peers advance
  that failed field group, other predictions are rolled back, and retry matches
  the previous continuous arrays and six-face records exactly.
- Vacuum physical arrays also match the previous standalone fixture exactly.
- Complete source and installed-wheel `examples/molecular_feedback.py` execute
  successfully. Their dipoles, E, currents, drives, electron number and sparse
  internal-energy freshness match the reference data exactly.
- Source and installed-wheel suites each run689 tests with one optional skip and
  no failures. The public orchestration test covers rejection, retry, unchanged
  clocks/state, caller ownership, immutable fields and mismatched-footprint
  rejection before updating.
- Installed API modules, paired guides, AI usage/manual JSON and example files
  are byte-identical to the working source.

The development MPI fixture now calls the public entry instead of maintaining
its own coupling transaction loop. Independent response replay, zero-field drift
and archive comparison remain explicit development validation operations. The
example has no hidden control calculation or background subtraction.

Canonical evidence: calc2026/salmonpy-c2h2-feedback-20261122/
public-feedback-release-summary.json and public-feedback-release-artifacts/.
Results are local OneDrive files; cloud synchronization is not independently
confirmed. No optional physical CI job was launched.

Remaining scopes include direct isolated-state persistence, wider media energy
observation, physical convergence, other sites/platforms and multiple molecules.
