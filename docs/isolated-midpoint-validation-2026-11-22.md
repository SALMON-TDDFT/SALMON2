# Align molecular driving with selected-field midpoint work

Scope: C2H2 fixed-ion PZ development coupling, using the already authorized GS.
No new SCF, no runtime rebuild or production-job changes, no manual heavy CI.
Standalone IsolatedTDDFTState.step(E) and the Fortran solver remain unchanged.

## Change

Fortran averages selected left/right E after own-vacuum subtraction. This constant
interval field drives the unchanged E3 worker. Default development sampling is
midpoint; explicit --isolated-sampling endpoint is a numerical control, not fallback.
Commit accepted right E only after accepted-midpoint validation. Cancel preserves
quantum state AND the previous selected left field. Report drive_time at the
physical field-clock midpoint, while state.time denotes the evolved endpoint.
No extra MPI communication. Source/objects audited for the new light kernel;
unchanged validated FDTD and E3 runtime reused with recorded SHA256 checks.

## Actual results

Electronic grid24×24×32,0.35Å,4 outer MPI ranks (one field rank plus bookkeeping
rank per field),64 field steps,dt=.003au,explicit pulse amplitude.1. Compare at
last fresh electronic-energy sample step60 (.18au), not stale final step63.

| Quantity | endpoint control, LD | midpoint, LD | midpoint, vacuum |
|---|---:|---:|---:|
| actual molecular work (Hartree) |6.98427541762e-6|8.42986769737e-6|8.50830001608e-6|
| selected-field midpoint work |9.85425956741e-6|8.42986769737e-6|8.50830001608e-6|
| internal-energy change |7.04692036635e-6|8.48846212875e-6|8.56687235995e-6|
| energy vs selected-midpoint work relative residual |28.4886%|0.695081%|0.688414%|
| actual vs selected-midpoint work relative difference |29.1243%|0|0|

Zero-field energy drift5.74708671763e-8Hartree is recorded and NOT subtracted.
Endpoint-control old physical/energy/work histories are bit-identical. Corrected
midpoint coupled trajectories intentionally differ because the applied field
sampling changes; do not claim all coupling results unchanged.

Inject rejection at interval30: midpoint quantum/field/norm/energy/work histories
match continuous execution exactly. Source continuous vs installed-wheel retry
matches every final array, including source-work ledgers. Prescribed actual-drive
replay into an independent quantum state remains exactly equal.

## Own-reference ledger and limits

Selected-field work is NOT total-field work. In LD, midpoint total-field work is
8.38021685319e-6Hartree and own-vacuum work is-4.96508441819e-8Hartree. These are
separately retained, alongside selected-field work. Subtracting the full own-vacuum
reference also removes vacuum radiation reaction; account for this declared model
policy before attempting total Maxwell/LD/flux conservation.

Source/wheel687 tests each, one optional skip each. Midpoint cancel/left-field
history, analytic source-work signs, scalar IEEE NaN/Inf rejection and finite
large-value averaging pass. Public dipole mapping, standalone quantum stepping
and the endpoint numerical control are preserved. Paired guides/AI JSON updated.
Full spatial electromagnetic energy/flux conservation, production physical
convergence and a general reusable two-way API remain separate tasks. Electronic
rollback remains accepted-history replay, not direct wavefunction restoration.
