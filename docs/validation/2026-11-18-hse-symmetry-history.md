# Symmetry histories and PT-CN — 2026-11-18

Extends the fast-PR explicit owner route; no conventional-HSE fallback. Existing
rank-local representative histories are connected to native/Python/MPI ownership
and numeric checkpoints. New symmetry combinations require an explicit history
capability; original factor-off Taylor does not. No fitting/propagation formula
or accepted-step order was changed. The PT-CN electron target now uses active
wtk(:nk), matching reduced rocc rather than its padded weight storage.

Native diagnostics expose eight rank-local integers without mutation or MPI:
Taylor history count/ready, CN history count/ready, Taylor prediction-installed,
CN predicted, Taylor factor step and CN accepted step. Tests confirm actual
prediction use, not only exact fallback. Queries are developer diagnostics, not
a promise of a stable user-facing C ABI.

The initial learned-CN test mistakenly requested interval2; native requires1.
The test was corrected and Python checks that actual native requirement before
initialization. No interval or solver setting is silently substituted. Learned CN
still uses predictor_corrector=n. Strict CN is not a native supported mode.

An existing fresh-owner checkpoint probe also exposed excessive optional-cache
matching. Cache presence is optimization state; saved first/count now pin to
trusted native kernel phase_start/phase columns before indexing/allocation,
independent of fresh optional cache allocation. Static map/decomposition,
dimension, offset, finiteness and truncation checks remain; commit stays staged.
Native stream encoding is unchanged, and fixed binaries were not patched.

Final runtime: /private/tmp/salmonpy-symmetry-history-release-build/build.json
from audited prepared /private/tmp/salmonpy-symmetry-history-source, matching
OpenBLAS/O2 profile. Local OMP/BLAS1 unless explicitly stated:

- Eight full8/reduced4-k comparisons: exact PT-CN, strict Taylor, learned Taylor,
  learned PT-CN, each internal MPI1/2,32 steps. Max current difference
  2.4665599424045226e-18; max energy difference1.0658141036401503e-14 atomic units.
- Two independent MPI2 material groups with distinct screening; within-group
  traces agree and material-group traces differ. Splits3/7/13 before/after
  training resume in fresh jobs with identical future current/energy/history
  arrays and final checkpoint bytes. Zero accepted steps are replayed.
- Incompatible field at step4 is refused without changing checkpoint bytes;
  corrected step succeeds. Learned Taylor/CN after training reject actual numeric
  overflow at step13, retain checkpoint bytes/output files and match independent
  no-fault continuous trajectories after retry.
- Learned Taylor/CN ordinary FDTD+twoLD+SBE+reduced-HSE examples match full-k serial
  fast-PR reference. Existing PBE/non-symmetric Taylor/PT-CN/seven-rank protected
  mixed regressions pass against prior numerical references.
- Learned CN27 full/15 representatives with nonuniform weights and uneven MPI
  ownership passes: max current difference1.734723475976807e-18, energy
  difference2.3248070135650778e-12, within the unchanged comparison tolerances.
  Fresh direct resume also matches. Circular-xy learned CN full/reduced passes
  with current difference1.463672932855431e-18, energy3.1086244689504383e-15.
- Native independent-resource and checkpoint probes pass16 combinations:
  full/ACE, stars/orbits, OMP1/2, two independent MPI processes. Both cache-on→
  fresh uncached owner and cache-off→cached owner pass; malformed/truncated/group/
  phase-range mismatches still reject without applying partial state.
- Python source/isolated installed suites and sdist-to-wheel/hash checks are
  recorded in the release evidence. No GitHub Actions launched.

Evidence: developer outputs salmonpy-hse-symmetry-history-20261118. The initial
mistyped test and old fresh-owner cache refusal are preserved separately from
successful final runs. Source snapshots/licensing remain intact; complete-file
and adapter hashes identify the runtime. Synthetic tests are not converged
material validation. One native lifetime/process, same binary/model/MPI partition
on restart, scalar CPU/full-q scope remain. Broader state ownership/handles,
other environments/scaling and portable archives are subsequent tasks.
