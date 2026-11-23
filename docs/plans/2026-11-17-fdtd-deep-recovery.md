# Native FDTD numerical recovery

Continue the approved defined-error/atomic-step objective after TDDFT/HSE. Reuse
existing numerical field/PML/LD payload ownership rather than disk replay. Keep
legacy batch eh_calc arithmetic/order unchanged with optional guarded status.
Guard the native supported plane-wave path only, reject TTM/LR/ASE/ART/energy
observation extensions before mutation until their full side effects are audited.
Collectively validate E after edge update, LD after update and H before ordinary
outputs; return known status and restore numerical payload directly. A checked
endpoint entry also validates sampled E and integrated A before output commit.
Unknown old/native process errors remain conservative; no whole mixed-step claim.

1. Add fail-first Python protocol/error mapping test and native MPI driver.
2. Add optional pre-output gate to eh_calc and detached in-memory backup to the
   native adapter. Validate current/shape/capacity/pending state collectively.
3. Build new external directory explicitly, production and compile-gated one-rank
   fault versions. MPI1/2, independent groups, first/warm/final faults, real
   overflow, strict payload/output preservation, exact old/new accepted fields.
4. Re-run mixed LD/SBE/TDDFT/HSE regressions and protected restoration as relevant.
5. Update paired guides/AI JSON/checklist, source and wheel tests, manifests,
   evidence and persistent artifacts. No remote Actions; warning-only48GiB budget.
