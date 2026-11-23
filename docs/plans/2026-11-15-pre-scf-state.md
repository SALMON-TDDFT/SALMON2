# Owned GS/DC pre-SCF phase control

Continue the approved explicit-input/state extraction. Preserve original SCF readiness, collective reduction, diagnostics, reset/rebuild and convergence ordering.

- Each system owns pre-SCF threshold, required consecutive steps, count and first iteration, plus the existing owned XC evaluation phase.
- Capture input once at GS construction; other SCF entry points initialize only if absent. Each SCF invocation resets count/start iteration as before.
- Extract local readiness and post-collective counting into internal helpers. Keep the MPI reduction and potential rebuilding in their original places.
- GS, orbital solves, DC fragment energy, DC total-system phase and hybrid metadata use the same system phase. Do not claim removal of DC communicator/material swapping or native shared exchange resources.
- Test two controllers with different settings against frozen original rules, global poisoning and failed configuration retention. Compare a small native GS warmup run and ordinary RT with the immutable fast-PR executable; local tests and wheel checks only.

Follow-up completed locally: MPI-enabled current/reference builds; two DC fragments with 2 and 4 total ranks; disabled/1/2/3 pre-SCF cases all match. The installed wheel's tool also passes the 2-rank cases. See `docs/pre-scf-validation-2026-11-15.md`; physical convergence and general concurrent solver isolation remain outside this completed regression scope.
