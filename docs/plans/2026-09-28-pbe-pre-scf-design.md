# PBE preconvergence before hybrid SCF

Approved in conversation: stabilize the density and occupied subspace with PBE,
then construct MLWF and converge the requested HSE06/PBEh40(+rVV10) functional.
No gap criterion for switching. DC uses the same specified positive electronic
temperature in both stages (300 K per user clarification), with a common chemical
potential enforcing the global core-weighted charge. This is GS occupancy only. Native real-space states and buffered periodic DC
fragments remain unchanged; DC Hartree remains global.

Use an opt-in `exx_pre_scf_threshold > 0` (default 0, disabled), in units of the
selected density SCF residual, and `exx_pre_scf_steps=3` consecutive qualifying
iterations. The switch is collective across DC fragments. Total nscf includes
both stages. Do not accept or export a pre-stage or unconverged staged result.
Initial scope: fresh static SCF, density mixing/convergence; exclude restarts,
checkpoints, optimization, shutdown timers and eigen/snapshot diagnostics whose
formats do not encode the temporary functional. RT uses the completed hybrid GS.

During PBE, use full PBE exchange+correlation, no rVV10, EXX, ACE or MLWF. Keep
requested xc identity for final output; a separate runtime stage controls actual
operators. At transition rebuild Vxc/local potential and fresh EXX/ACE before the
next solve; reset mixing age and residual. Adaptive support becomes eligible,
but still requires successful localization. No stale ACE exists because no EXX
was initialized during the fresh PBE stage. Existing MLWF tolerance is unchanged
(the input default is 1e-6; previous measurements explicitly used 1e-7).

Alternatives: current full-hybrid warmup costs EXX from the outset; a separate
PBE restart requires user orchestration and cannot ensure consistent cache reset.
The single-run stage is selected for explicit, auditable transition handling.

## User steering: DC fragment SCF without extra localization

DC buffers already bound the exchange domain. Default `yn_exx_dc_mlwf='n'`
bypasses spread minimization, gauge seeding/transport and orbital masking in DC
SCF, keeping full fragment exchange and ACE. `yn_exx_dc_mlwf='y'` retains the
old path for explicit comparisons. Native whole-system SCF/RT is unchanged.
Reject fractional/fixed-radius masks and pair screening with DC MLWF off;
fraction 0 or 1 both mean full support. Cover spatial/orbital distributions
and the existing full-k fragment backend. No zero-T occupation changes.
