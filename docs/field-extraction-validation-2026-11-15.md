# Coupled field A/E extraction — 2026-11-15

NativeEHField and MultiscaleField now expose sample_fields(), a local read of
detached SampledField values. No update or collective is performed. Only usable
fields after their first step are accepted; H is absent, not zero-filled.

Native point order/coordinates follow registered macro cells; A is at (n+1/2)dt,
E at n dt under the existing scheme. Multiscale positions use material-cell x
and y=z=0; A is n dt and centered E is (n-1)dt. No Fortran source/ABI changed.
The existing native TDDFT coupling example additionally writes field_samples.npz
with independent clocks, positions and point IDs. That full TDDFT example was
not rerun in this extraction check; native FDTD extraction was exercised below.

Evidence on Apple Silicon macOS/Python 3.14:
- Unit tests: initial/closed rejection, clocks, positions, retained-copy isolation,
  and no mutation of step index. Focused coupled tests: 11 passed.
- Actual native FDTD impulse regression, one/two field ranks (2/3 total ranks),
  one thread, each under 60 seconds. New samples exactly match native A/E arrays
  and intended coordinates; rank-local copies agree. All saved radiation arrays
  are identical between decompositions.
- Source suite: 418 tests, one skip. sdist→wheel build/install successful.
- Installed suite: 418 tests, two skips (optional plotting and filesystem).
- Independent read-only review found no required fixes in extraction changes.

Logs: /tmp/salmonpy-field-extraction-{2,3}.log,
/tmp/salmonpy-field-extraction-tests.log,
/tmp/salmonpy-field-extraction-installed-tests.log.
Native outputs: /tmp/salmonpy-field-extraction-native-{2,3}/.
Wheel: /tmp/salmonpy-field-extraction-dist-20261115/.

Remaining: native H accessor, explicit current-clock mapping, and packet-based
feedback scheduling. History retention is caller-owned; this API returns only
the current A/E samples. No physical convergence claim or Actions run.
