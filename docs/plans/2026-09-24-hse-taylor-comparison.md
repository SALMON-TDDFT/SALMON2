# Native HSE Taylor4 / PT-CN comparison

Execution: superpowers:executing-plans, inline in the existing TDCDFT branch.

Goal: verify that the new propagator approaches the same Si current and density
as SALMON's original fourth-order Taylor propagation as dt is reduced.
User authorized the comparison after the proposed same-HSE, same-state,
time-step-convergence design; no further design approval is needed.

Architecture: call unchanged src/rt/taylor.f90, n_hamil=4, using SALMON's
predictor/corrector local potential average. Extend that average to the Fock
operator, retaining the initial occupied source while forming the predicted
endpoint source. Neither Hpsi trial powers nor predictor states become the
accepted physical source silently. Use named hse_taylor4 (ACE) and
hse_taylor4_full (full operator on every trial) modes; initial conditions,
occupation, field, physical kernel and no-renormalization policy are identical.
Full-source averaging is concatenation / sqrt(2), equivalent to averaging the
two density matrices. ACE averaging similarly concatenates the two factors.
Restart records already distinguish propagator names; require predictor/corrector
and n_hamil=4 for these modes and reject incompatible combinations.

Tasks:
1. Add failing independent midpoint-ACE/kernel tests; implement averaged factors.
2. Add the source lifecycle and full-action audit path to hse_native; connect
   the existing predictor/corrector and input guards without changing Taylor.
3. Build and run identical first-step ACE/full audits, then choose a dt grid
   based on stability. Compare native PT-CN and Taylor over a common physical
   interval at 4^3 k points, reducing dt and reporting density, current, norm
   and energy. Do not mistake short-time agreement for a resolved spectrum.
4. Review, run native/unit and semilocal regressions, document reproducible
   cases and evidence, commit results. Keep the separate long reference job.

Acceptance: evaluate convergence, not just energy conservation. Compare
projectors/density rather than raw gauge-dependent wavefunctions. Use a fine
full-Taylor audit to bound ACE effects before attributing differences to PT-CN.
If the requested agreement fails, report the discrepancy and investigate;
never tune acceptance thresholds to declare success.

## Completed evidence

All four tasks completed. HSE-enabled and disabled builds succeed; 79 numerical
unit tests and six Si/TDCDFT CTests pass. Independent final source review found
no blocker. Restart is bitwise identical. Results, exact inputs, and figures
are archived in docs/results/si-hse-taylor/. Short dt refinement establishes
current convergence; the 32-au comparison gives current differences of 0.473%
(PT dt=.32) and 0.149% (PT dt=.16) against Taylor dt=.08. Induced-density
errors are larger (4.31% and 2.93%); no equivalent-density-accuracy claim is
made. Full-Fock audit at 8 au bounds the measured ACE current difference at
0.00711%. Coarse Taylor dt=.32 failed the electron-count check and is excluded.
The independent long optical-spectrum reference job remains untouched.
