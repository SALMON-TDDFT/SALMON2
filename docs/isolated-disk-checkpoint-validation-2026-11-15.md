# Isolated disk replay validation — 2026-11-15

Added save_checkpoint/restore_checkpoint to IsolatedTDDFTState. The schema
salmonpy.isolated-tddft-checkpoint/1 explicitly stores three-component E and
six-component dipoles (electronic change, total charge dipole). GS is separate.
Existing periodic schema/names/layout remain unchanged. The shared serializer
uses explicit format-specific dimensions, file sizes and decoded-data budgets.

Evidence on Apple Silicon macOS/Python 3.14, OMP/BLAS one thread, at most two
simultaneous native workers; no Actions or physical convergence sweeps:

- New tests: reconstructed initial model, corrupted chunk/budget/position refusal
  before spawning, unchanged target, empty history and 1025-row chunk boundary.
- Source suite: 410 tests, one skip, before the final boundary test addition;
  updated isolated tests: eight passed.
- sdist-to-wheel build and isolated installation succeeded. Installed full suite:
  411 tests, two skips (optional plotting and case-sensitive-filesystem test).
- Existing C2H2 E1 fixture: save at step 10, reconstruct model, restore and continue
  to step 20. Dipoles exactly equal uninterrupted continuation (max error 0).
  Source and installed wheel both pass, along with existing isolation/recovery.
- Installed isolated_response --checkpoint writes both the response archive and
  isolated replay directory successfully.
- Independent read-only review found no required correctness fixes; synthetic
  periodic/isolated roundtrips at 0, 1, 1024, 1025, 2050 rows checked exact arrays
  and format-derived size ceilings.

Runtime SHA256:
abb4fce2106b92c9158e6c77d4111480367bad74cfc926947086d6f649b1c62b

Local logs: /tmp/salmonpy-isolated-disk-{full-tests,targeted,build,installed-tests,
installed-native,sample}.log. Native results:
/tmp/salmonpy-isolated-disk-native-20261115/report.json and
/tmp/salmonpy-isolated-disk-installed-native-20261115/report.json.
Wheel: /tmp/salmonpy-isolated-disk-dist-20261115/.

The archive replays accepted history, not native wavefunctions. Package version,
initial data, executable, capacity, dt, threads and position must match. Replay
cost grows with history. max_bytes is not RSS control. MPI whole-session restart,
E runtime publication and other platform verification remain outside this change.
