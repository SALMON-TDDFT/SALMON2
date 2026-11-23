# Diagonal-background field energy observation — 2026-11-22

## Implemented numerical scope

Extend only the read-only six-face observer to positive nondispersive electric/
magnetic backgrounds and diagonal LD with zero Ohmic conductivity. LD permeability
retains the existing solver requirement1. Do not alter E/H or LD updates, physical
checkpoint state or communication counts. Off-diagonal/subpixel observations and
Ohmic work remain open; unsupported observation warns/unavailable, calculation
continues. This is a completed bounded stage, not full tensor-energy support.

Use actual staggered coefficients. With electric coefficient c1 and current
coefficient cJ, recover c4=-4*pi*dt/cJ and epsilon=c4*(1+c1)/2. For LD,
c4=epsilon+2*pi*dt*sum(c2pole), c5=epsilon-2*pi*dt*sum(c2pole), so the implicit
pole contribution cancels. Recover magnetic mu from the actual H curl coefficient
and its derivative spacing. No material-ID ownership is reconstructed.

At integer step n the modified stored energy uses epsilon*E(n)^2 and
mu*H(n-1/2)*H(n+1/2), summed over staggered components. Reconstruct prior H with
inverse mu. Keep the separate conventional centered-H energy. Poynting face work,
midpoint LD work and external-current work are unchanged. The legacy unit-background
arithmetic path is retained. Native/shared kernel files are byte-identical.

The new C-bound Fortran kernel sp_fdtd_box_energy_diagonal takes nine read-only
coefficient buffers. Scalar finite/positive-mass checks precede output publication.
The existing17-value field observation interface and one owning-group20-value
reduction are unchanged. No new global or per-point communicator operation.

## Verification

A failing missing-symbol test preceded implementation. The final single additional
unit test covers nonuniform dx=(.8,1.2,1.5), all-axis epsilon/mu storage, a
manufactured Yee step including midpoint LD/external-source work, implicit-LD
coefficient cancellation, partition additivity and 0/NaN/Inf rejection with output
unchanged. The original manufactured vacuum test also passes.

Real field-only tests use40^3 cells, h=2bohr, dt=.003au,64 steps, a finite material
box and an incident plane pulse outside the inner observation box. No GS/SCF.
Maximum absolute modified-energy step residuals (Hartree):

| Profile | Maximum residual |
|---|---:|
| dielectric epsilon2.5, mu1.4 | 3.349062563731686e-14 |
| scalar LD epsilon_inf4 | 2.956586798302796e-14 |
| diagonal LD epsilon_inf(2,3,4) | 3.849492339474958e-14 |
| previous unit-background scalar LD | 4.1283057701513304e-14 |
| diagonal LD with epsilon_x=1 | 2.8873222248762044e-14 |

All64 observations are available in every profile. The x-unit diagonal case
checks that mode selection does not assume epsilon_x represents all components.
Conventional centered-H residuals remain separately stored (maximum about
1.68e-3Hartree in these short runs). Modified conservation does not certify
continuum energy accuracy or grid/time convergence.

With/without observation: point E/A histories and complete per-rank final field
binaries (E/H,LD,PML,current) match exactly for all four main profiles. Prior
library diagonal-LD numerical state is byte-identical; prior standard-profile
observations are identical. Existing authorized C2H2 PZ/LD/MPI6 all64 frames and
energy records are unchanged; no new SCF.

Field MPI1/2 samples are exact; energy aggregates agree at rounding tolerance
(atol2e-13,rtol2e-12). Direct save at32 and a separate MPI job restore match all64
observations/energy records and final state bytes exactly. Installed-wheel tests
repeat all five profiles and direct restart with identical source results.

Source and installed wheel suites each ran693 tests, one optional skip. Packaged
API/paired guide/AI JSON bytes match source. Tests run from the checkout against
the installed wheel; tests themselves are intentionally not installed. The first
artifact-comparison harness incorrectly assumed test files were installed; its
packaging-only check was corrected without rerunning successful calculations.

Tracked-only source export/build, native source and application-object prohibition
audits pass. New source local NLCC XC/core-force tests pass at MPI1/2, without GS.
The first source fixture launch lacked an explicit development LD kernel path;
it was corrected, with no calculation failure attributed to the solver. No manual
physical CI dispatch.

## Reproduction and artifacts

Explicit development tool: tools/test_diagonal_field_energy_mpi.py. It executes
calculations and is not an input-only management-node check. Require field ranks
plus one bookkeeping rank. Profiles include dielectric,scalar-ld,diagonal-ld,
diagonal-unit-x,legacy; --observe enables read-only accounting; --checkpoint NEW
and --stop-after32 prepare direct state; --restore NEW resumes in a new job.

Canonical local OneDrive campaign:
calc2026/salmonpy-diagonal-field-energy-20261122/summary.json and artifacts-verified.
C2H2 legacy regression remains in its existing authorized campaign under
diagonal-energy-molecular-legacy. Local storage does not prove cloud synchronization.

Source inventory SHA256:
b2fac8b54846c220e6660263cbe31a50d44aa96195dc5efaab91e0b17f371d26

Field library SHA256:
c0fc2706ec6bf60ce787e689bc001f3f963f0bd24ee158794e56be2c961bf66e

Wheel SHA256:
548f446bb28eadb23409d7a8331292498c6c21ce58cdf848cea868b1558a9635

Remaining: off-diagonal/subpixel constitutive interpolation and any required
boundary exchange accounting; explicit Ohmic dissipative work; broader physical,
platform and facility validation. Do not reuse this diagonal formula for them.
