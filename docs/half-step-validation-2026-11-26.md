# FDTD / electronic half-step audit — 2026-11-26

## Clock convention verified

Use integer k counting dt/2 on the FDTD clock. For native Yee update n:

| Quantity | k |
|---|---:|
| old A, old electronic state, previous crystal current, old raw H | 2n-1 |
| E(n), centered H(n) | 2n |
| new A, new electronic state, new crystal current, new raw H | 2n+1 |

The crystal midpoint Hamiltonian reads the mean of old/new A, centered at2n.
A_new-A_old=-dt*S(E_n). Raw H brackets at2n-1 and2n+1 are averaged for the molecular
component-grid reader; the separate macro H sampler returns raw H at2n+1.
Electronic local step n is field half-step2n+1. Initial GS has local step0 / field
half-step1; this assigns an initial time and does not replay an earlier interval.

The isolated worker explicitly rejects propagators other than middlepoint.
Source inspection and actual response tests establish that molecular update
n-1 -> n reads field frame n only. Packet frame n+1 is prefetched, not a second
half-step propagation and not averaged with frame n. The earlier explanation
that the two frames represent two ETRS half-steps does not apply to this interface.
No numerical phase shift was needed in the current supported connection.

## Contract/output corrections

The uniform molecular FieldPacket contract had incorrectly labeled E as the
new electronic endpoint. It now declares offset1/2, lag1: input E is at local
(n-1/2)*dt. A valid first packet at dt/2 is accepted; an endpoint-labeled packet
is rejected. Scalar step(E) and the Fortran numerical update are unchanged.

MPI controller.clock exposes exact integer half-step numbers. One-way examples
retain local time and add electronic_field_time / electronic_field_half_step.
Prepared-history output adds the recorded source time origin to the field clock.
Use those field-clock times for comparisons with E/H instead of treating the
same array index as the same physical time. No per-step global agreement or
numerical interpolation was added.

## Numerical checks

Manual tool: tools/test_molecular_timing_mpi.py. Source execution uses
`python -m tools.test_molecular_timing_mpi` from the checkout; installed wheels
include the tool under salmonpy/_runtime_source/tools. The tools run real
calculations only when explicitly launched in an allocated calculation context.
They are not part of required material-physics CI.

1. Actual two-rank FDTD with nonzero fields from prescribed current: checked
   A differences against -dt*E(n) and compared component-centered H against the
   independently sampled previous/current raw macro H. H comparison starts at
   n=2 to use two separately read raw samples; supported initial fields are zero.
   Sample timestamps are E=n*dt, A/raw H=(n+1/2)*dt. Spatial transverse H averaging
   is independently reconstructed from component samples across the rank boundary.
2. Existing PZ C2H2 GS, dt=.003 atomic time,12steps, spinless fixed ions. A linear
   uniform grid field E_n=n*dt*(1,-.5,.75) matches independent uniform-E-extra
   driving at the same midpoint to8.48e-20 in the seven dipole/change/number
   observations. No new SCF was performed.
3. A deliberately incorrect half-step-ahead drive differs in dipole by9.61e-6
   e bohr, well above the2e-12 absolute comparison tolerance. This detects a
   half-step mistake rather than comparing only constant fields.
4. Replacing only future E(n+1) by100times its value leaves every observation
   bit-identical. A separate orbital-magnetic test with nonzero H_n also remains
   bit-identical when only future H(n+1) is changed by100times. Thus the current
   supported propagator's use of frame n is measured for both E and H.
5. Direct process-exit restart after6steps with molecular MPI1 ->2 preserves the
   half-step labels and matches continuous response within roundoff; measured
   restoration updates are0. Magnetic restart is checked separately.
6. Actual LD/FDTD plus molecular TDDFT direct restart at48/49, molecular MPI4 ->1,
   continuing to96/97: all previous numerical arrays remain bit-identical to the
   earlier accepted run. New labels report electronic half-step193 versus final
   FDTD E half-step194, and zero restoration updates.

The unchanged audited library was verified by SHA256 before reuse:
`d2da24af566630a23efa09812d2ed40b855f71e2cfda438c9dc7ebc37957ecfa`.
Its source inventory is
`47531303883faa8c36fefa76b70ba55269448730843d6b08554d6e52de4f910e`.
No new solver executable or native source/object audit was needed at job startup.
The physics remains in Fortran; Python constructs test input and independent
comparison oracles, validates clocks and writes metadata.

Results, summaries and logs are collected under canonical calc2026 in
salmonpy-half-step-20261126. The early grid/ attempt failed before calculation
because script execution imported an older installed package; grid-final/ and
subsequent runs explicitly selected the source implementation. That attempt is
preserved. No cloud synchronization or physical convergence is claimed.

## Scope

Verified current native E/H/A chronology, isolated middlepoint E/H driving,
contract admission and one-way restart. Crystal midpoint-A selection was inspected
in Fortran and the real field A/E identity was checked; a new material-specific
crystal response campaign was not run. This audit does not validate historical
molecular radiation/current-return timing, ETRS support, SOC, unequal-dt coupling,
subcycling, arbitrary initial nonzero fields or multi-node execution.

## Distribution and final checks

The actual nonzero-field probe has max A-difference/E residual3.31e-24 and
raw-H centering residual1.06e-22, with magnetic signal5.98e-7. Electric restart
MPI1->2 differs by at most5.33e-15; orbital-magnetic restart differs by1.78e-15.
Both report zero restoration updates and electronic field half-step25 at step12.

Built an sdist and then the wheel from the extracted source. Installed into an
independent environment outside the checkout; all79targeted source and installed
tests pass. Installed public modules/guides/AI JSON and the packaged manual audit
tool match the source bytes. The installed tool's actual FDTD probe and magnetic
process-exit restart are bit-identical to the source counterparts. Guide/JSON
synchronization and diff checks pass. No heavy physics CI was dispatched.

## Follow-up: crystal endpoint projector phase

The earlier clock/source audit did not validate ordinary-XC endpoint observations.
A separate frozen-state numerical audit subsequently reproduced a stale midpoint
nonlocal phase in ordinary-XC endpoint current and midpoint-A energy evaluation.
The maintained solver now installs endpoint A and its projectors together before
observations; the propagation-stage A is unchanged. See
[the endpoint observation validation](crystal-endpoint-validation-2026-11-26.md).
The earlier molecular numerical checks are not evidence that the old crystal
observation path was correct.
