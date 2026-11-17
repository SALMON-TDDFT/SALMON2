# Internal fast-PR symmetry TDDFT — 2026-11-18

Explicit internal Fortran `rt_open(...,allow_symmetry=.true.)` connects the existing
owned symmetry exchange/checkpoint core to full TDDFT state. Ordinary worker,
native-C, Python and FDTD calls omit the opt-in and still refuse symmetry. No
new public retry ABI or arbitrary same-process context guarantee is claimed.

Bounded profile: scalar periodic CPU/full-q, one internal k rank, factor-off
Taylor/ACE; full-action only as a same-fast-source comparison. Retained groups
preserve z. One initialization per process. Density/current still consume common
immutable symmetry tables; step/current/checkpoint/staged commit compare flags,
operations, direct/reciprocal lattice and their inverse matrices against owned
values. Invalid initialization assets can still terminate native execution.

Native field status34 rejects incompatible A before backup/propagation. It is an
internal status, not the current public C/Python recoverability contract. Corrected
A can be retried without accepted-step replay. Staged restore rechecks configuration
at commit and retains the pending candidate on rejection for a corrected retry.

Symmetry archive copies canonicalize overlap-cache dirty flags to true: overlap
caches are not serialized and are invalidated after restore. Live flags, all
scientific arrays and non-reduced V1 serialization remain unchanged.

Local Apple Silicon evidence (no Actions):

- Fail-first missing opt-in signature. Separate real executions reproduce shared
  flag, inverse-only drift and staged-commit drift gaps before the added guards.
- Synthetic8^3/full8-k and reduced4-k C2z/mirror-x, origin Si KY potential and two
  orthonormal orbitals; no SCF or physical-convergence claim. External z A drive.
- ACE/full-action × MPI1/OMP1 or MPI2/OMP2 host processes. Each process attaches
  COMM_SELF; MPI2 uses different screening coefficients and12/13 step counts.
  This is not reduced internal-k MPI2. Full/reduced charge current and energy
  match within1e-12 absolute/1e-10 relative comparison tolerances (maximum observed
  current difference1.16e-15, energy difference1.78e-15 atomic units).
- Odd/even6/7 checkpoints: direct and staged fresh-process resume match API arrays,
  orbital/density/potential/ACE diagnostic traces and complete final checkpoint
  bytes exactly. Only remaining accepted steps are called, no elapsed replay.
- Forbidden tiny x A, shared flags/operations/inverse drift and staged-commit
  drift reject while preserving original data. Restoring inputs permits exact
  retry. Truncated/trailing archives reject with identical before/after bytes.
- Default entry, learned mode, PT-CN, a second lifetime and internal-k2 reject.
- Existing non-reduced MPI2 GS, PBE/HSE Taylor/PT-CN and seven-rank field/LD/SBE/
  TDDFT/HSE protected retry/restore comparisons remain exact. PZHEEV audit passes.
- Source and installed545 Python tests(one skip each), sdist-to-wheel/native hashes,
  installed mixed/1D smoke and packaged native tools/probe validation are retained.

Evidence: `~/.codex/outputs/salmonpy-symmetry-tddft-20261118/`.
Memory budgets remain advisory. Next: explicit public native/Python opt-in and
retry protocol with protected mixed testing, then reduced internal-k MPI layouts.
