# Crystal TDDFT endpoint-A observation correction — 2026-11-26

## Defect and change

For ordinary periodic, non-single-scale RT, propagation installed midpoint A/c
in both the kinetic operator and the nonlocal projector phase. That stage was
consistent. At observation time only the kinetic A was changed to the endpoint;
the projector phase remained at the midpoint. Eigen/kinetic/nonlocal energy was
also evaluated before installing endpoint A. The pre-existing HSE observation
branch already installed endpoint A and its phases together.

`time_evolution_step.f90` now installs endpoint A and, for ionic systems, its
nonlocal phases together before projection/energy/current/force observations.
The redundant HSE-only refresh at current readout is removed. Propagation,
density/Hartree/XC updates, molecular E driving and the spatial-A single-scale
branch are not reordered. No additional diagnostic communication is added.

This intentionally corrects current and energy values. Unchanged propagation is
claimed for fixed ions and a prescribed A history, not for feedback calculations
whose future A depends on the corrected current, nor for moving-ion trajectories.

## Independent numerical audit

Manual developer tool: `tools/test_endpoint_observables_mpi.py`. It does not run
as required material-physics CI. The `--instrument PREPARED_SOURCE NEW_SOURCE`
operation adds a test-only Fortran probe to a NEW copied tree. Source/object
prohibition checks remain required. `--link-probe BASE_BUILD_JSON NEW_SOURCE
NEW_OUTPUT` compiles both RT-step and probe objects and audits every borrowed
application object before linking. Other borrowed source hashes must match.
No probe or extra runtime diagnostic API is installed in the production library.

The Fortran probe copies the accepted wavefunction and evaluates:

- Current with the cached phase and independently refreshed endpoint phase.
- Endpoint kinetic and nonlocal energy using the frozen wavefunction.
- Centered differences of frozen-state E_kin+E_nonlocal with A/c +/-1e-5 in
  each of x/y/z, divided by cell volume. This equals the electron current;
  the public state-handle result is the opposite charge-current sign.

Only explicit-A energy terms vary in this frozen-state derivative; it is not a
functional derivative of the evolving density or a proof of TB-mBJ variationality.
Diagnostic copies of the phase, A and orbital work arrays are restored/discarded;
propagation is not performed by the probe. Full owned wavefunction bytes are
hashed in memory; raw probe wavefunction arrays are not saved. Tests use an
existing Si PZ GS in place (no new SCF or GS copying),16^3 grid,32electrons,
32states, original SALMON2x2x2 k mesh, fixed ions, transverse middlepoint,
dt=.02au and8steps. A/c is nonconstant in all three Cartesian components.
Libxc-PZ/PBE and TB-mBJ start from the same existing orbital seed for controlled
software comparison; these are not converged material predictions for each XC.

## Results

| Check | Result |
|---|---|
| Old PZ returned current vs refreshed endpoint current | mismatch2.8071111376e-6au; expected regression failure |
| Old PZ energy component vs endpoint recomputation | mismatch4.3667267680e-4Ha |
| Old PZ total energy correction | up to4.0936447029e-4Ha |
| Corrected current / kinetic / nonlocal energy vs endpoint recomputation | bit-identical |
| Corrected electron current vs independent frozen-energy A derivative | max8.59e-13au |
| Old vs corrected owned wavefunctions, all8steps | bit-identical |
| New-process direct restore at4, then steps5..8 | wavefunctions and all18probe observables bit-identical to continuous |
| Production library vs instrumented library | currents and complete RT/energy streams bit-identical |
| Production direct-resume written rows5..8 vs continuous | bit-identical at matching times |
| Checked vs performance | bit-identical |
| Public state-handle charge current vs legacy electron current | exact opposite sign; direct resume and performance results exact |

PZ covers serial and separate two-rank k/orbital/spatial decompositions; their
observations match with atol2e-12, rtol2e-10. PZ, Libxc-PZ/PBE and TB-mBJ each
cover k-MPI2, old/new wavefunction equality, direct restart, checked/performance,
and production singleton/state-handle routes. The separately initialized
production state restores disk data before any accepted propagation step; no
steps1..4 are replayed. RT current files begin at step1; energy files retain an
initial row, so resumed output comparisons select matching steps5..8, not equal
raw row indices. Changing the MPI decomposition during this RT restore was not
part of this test. Normal/restarted production output uses the uninstrumented
library, so agreement is not based only on diagnostic recomputation.

Source and installed-wheel targeted tests:119passed each (source inventory,
preparation, manual, timing, shared states and checkpoint contracts). An extracted
sdist builds a wheel with the updated source inventory, guides and manual tool.
Installed-wheel PZ handle continuous/direct-resume output and Libxc-PBE frozen-state
audit also match source runs exactly. The installed guides, JSON, audit tool and
solver source bytes match the checkout. No real-material reproduction was added
to required CI or manually dispatched.

## Artifacts and scope

Canonical local records: calc2026/salmonpy-crystal-endpoint-20261126, including
comparison.json, per-rank observations/wavefunction hashes, RT checkpoints,
actual RT/energy files, build/audit/test logs and validated-runtime/build.json.
Initial failed harness attempts are preserved separately from the successful
regression: one obsolete helper argument and an omitted diagnostic timer context
were corrected before accepting numerical evidence. Expected old-code physics
regression failures are also retained. Cloud synchronization is not confirmed.

Maintained source manifest:
`ca157982260f82779077e1ade8786d2207a4bb381ec3de14eb305bff00e91668`.
Production library SHA256:
`1efb6bd8289a56b319eef15ec01ec217c1e23e16a7ba35e24d5aa3960931bc73`.
All source/object eigensolver checks passed before the new libraries were linked;
subsequent unchanged-library reuse checks the recorded SHA256.

Single Apple Silicon Mac, caller-owned MPI, one thread per rank. This is a
numerical/software consistency audit, not optical convergence, multi-node,
GPU, SOC, moving-ion, HSE or unequal-dt validation.
