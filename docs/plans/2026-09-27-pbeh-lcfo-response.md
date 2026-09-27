# PBEh40 LCFO electronic response

Approved scope: reuse MLWF/ACE exchange for large water and solution optical response, retaining the existing separate branch. Fixed nuclei, Gamma, occupied-only Taylor4 response first; DC-MD forces remain a separate milestone.

## Task 1: connect functional and total-density potential
- Add an end-to-end DC ground state → LCFO response regression for PBEh40 and PBEh40+rVV10; observe current rejection.
- Generalize LCFO exchange fraction/kernel and local FFT selection. Keep conventional PBEh RT blocked.
- Gather density and sigma over each real-space communicator for a root rVV10 evaluation, scatter derivatives back through existing discrete gradient adjoint. This is a correctness reference with global scalar storage, not distributed FFT scaling.
- Validate time-step refinement, finite energy/current, orbital-layout parity and HSE regression.

## Task 2: functional provenance
- Write per-fragment functional metadata bound to the binary run ID, containing functional, exchange parameters, rVV10 parameters and source cutoff.
- Require matching metadata for PBEh reconstruction and verify every fragment. Preserve legacy HSE files lacking metadata.
- Reject model mismatch, changed parameters, stale metadata and finite-source-cutoff GS data in the new PBEh RT route.

## Task 3: verification and documentation
- Build HSE on/off, run functional/input/DC/LCFO regressions and one independent review.
- Record numerical evidence and restrictions, commit changes without pushing.

## Review focus
Communicator duplication of nonlocal energy, global grid indexing, exchange fraction in the continuity diagnostic, Coulomb cutoff resolved on fragment grid, metadata association with each fragment/run, unchanged HSE behavior, unsupported paths remaining guarded.
