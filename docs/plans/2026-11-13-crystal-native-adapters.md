# Next native crystalline adapters: ownership and exact timing

Status: common Python A-only protocols and batch orchestrator implemented;
SALMON-native adapters and initial-state file loaders not implemented.

SBE first adapter should expose the velocity-gauge path already present in
src/ssbe/realtime_ssbe.f90. Its exact interval convention is:
- dt_evolve_bloch(sbe, gs, (A_old+A_new)/2, dt)
- calc_current_bloch(sbe, gs, A_new, J, icomm)
The reusable model owns s_sbe_gs_info and validated band/k-point/occupation data.
Each evolution state owns s_sbe_bloch_solver (rho and every allocated history),
previous A, dt, time index, correction settings and communicator ownership.
new_state must deep-initialize density-matrix state without rereading/altering
another run. No test should collect external GS data merely to exercise the API;
use controlled small numerical fixtures and existing authorized datasets.

The existing length-gauge branch computes E using A(it+1)-A(it-1). It cannot be
silently mapped to an endpoint-only step(A_new) while claiming identical results.
Do not expose it until its look-ahead requirement has an explicit A-only contract
and equivalence tests. Public crystal E input remains prohibited by design.

TDDFT state extraction must cover main_tddft/time_evolution_step's complete context:
system, rt, wavefunction ping-pong buffers, temporary orbitals, rho/spin density,
local potentials, Vh_stock1/2, Poisson/stencil/PP/XC data, communication buffers,
and any predictor-corrector histories. Geometry/basis/functional/MPI layout are
model/context settings, dt/propagator belong to the evolving state, A endpoints
are per-step input. I/O and observations must not own physical state.

Acceptance for each native adapter:
1. Existing standalone and multiscale temporal ordering preserved.
2. One long call equals the same sequence of step calls.
3. Two states from one initial model remain independent when interleaved.
4. Currents match the legacy path under the same MPI/thread configuration.
5. Defined gauge, laboratory/crystal rotation, volume normalization and units.
6. Backend lifecycle/free paths and failure rollback verified before publishing.
