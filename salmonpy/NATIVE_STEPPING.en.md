# Keep the previous state when a coupled time step fails

[日本語](NATIVE_STEPPING.md)

When FDTD, LD and electrons advance together, a failed update must not leave
some components one step ahead. This operation keeps the previous fields,
material states and time, and restores them together if the update fails.

Use `salmonpy.native_steps.native_step_targets` to collect supported components
as `StepTarget` objects, then advance with `step_protected`. Build the shared
library against the same MPI implementation as the calling program. Unsupported
libraries are rejected before time evolution starts.

```python
from salmonpy.native_steps import native_step_targets
from salmonpy.stepping import StepRecoveryError

targets = native_step_targets(coupled)
try:
    frame = coupled.step_protected(targets)
except StepRecoveryError as error:
    if not error.original_preserved or error.cleanup_errors:
        raise  # This lifetime must not continue.
    # Accepted state is retained. Inspect/fix the cause before explicitly retrying
    # this SAME step; do not increment a loop clock on a rejected attempt.
```

All ranks participate. The factory supports exact built-in types: NativeEHField
(including PML and every field-local LD), one-rank MultiscaleField, standard solver
periodic/unpolarized/fixed-ion middlepoint TDDFT, HSE Fortran HSE calculation and the ordinary
sample's built-in SBE endpoint wrapper. Independent solver groups retain their
own numerical backups; FFT plans, communicators and owning resource handles remain
live. Wavefunctions, density/potentials, XC payload/cval, ACE/history/caches,
field/PML/LD/current and accepted clocks restore directly without replaying past
propagation. Python A/E/current/clock/flag references also restore. This does not
add solver multi-context reentrancy or a stable public C ABI.

For user calculation implementations, pass explicit `custom_targets=(StepTarget(...), ...)`, with
owner identities matching the live local objects. Unknown wrappers/subclasses
are not certified by method-name inspection. Callbacks must cover their complete
state and side effects and return quiescently on every group rank. Rebuild targets
after a successful checkpoint publication replaces state objects; stale owners
reject before capture. Staged saved solver states and active step backups are
mutually exclusive in both Python and the solver entry points.

## Output scope and failure handling

Supported output is the ordinary current/energy RT streams and plain E/H point
observations, directly in each owner's private working folder. Capture flushes
previously accepted buffered records before recording the baseline. On rollback,
solver RT writers close, appended attempt records are truncated to that baseline,
and the same units reopen at the accepted end. Successful retirement flushes
accepted records. Default `step()` retains its previous path and arithmetic.

Only private regular single-link `.data` streams are audited. Different output
roots/subfolders, auxiliary density/microscopic-current/ELF/electrostatic,
projection, ASE/ART, frequency-observation or restart-output configurations are
rejected before stepping where applicable. Options are not silently disabled.
File replacement/removal or an unexpected new stream fails closed instead of
being silently deleted. Do not use concurrent readers/writers or move live
working directories: tentative appends can be visible until rollback. This is
not durable filesystem publication or process/crash/MPI failure recovery. Disk
failure, allocation exhaustion, initialization/third-party STOP and unreturned
calls to the Fortran/C library remain excluded. Unknown standard TDDFT called directly from Python failures now latch
failed just like solver E/H/HSE; a later valid-looking return cannot revive them.

The common StepRecoveryError reports whether every rollback completed. Retry only
if original_preserved is true AND cleanup_errors is empty, after addressing the
cause. No automatic retry, pulse attenuation or timestep change is performed.
Additional whole-step and per-component backups coexist. Memory advice reports
uncalibrated overhead and observed rank peaks; budgets remain warnings, not caps.

## Ordinary examples

Add `--protect-steps` to `salmonpy.examples.mpi_mixed_native`, or to the recommended
one-dimensional `salmonpy.examples.mpi_multiscale_tddft` entry. The ordinary CLI
still stops on an unhandled error; custom hosts can catch a safely rejected step.

The installed `examples/mixed_tddft_hse_layers.json` adds three electronic points:
SBE, TDDFT and HSE, alongside two field-local LD materials. Locate it using
`importlib.resources.files('salmonpy').joinpath('examples/mixed_tddft_hse_layers.json')`.
Use a separate `--hse-bundle` when the HSE initial state differs from positional
BUNDLE. Its files enter HSE-specific identity and checkpoint compatibility;
changing that bundle does not change the TDDFT material identity. No DFT, download
or source conversion happens implicitly.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 7 python -m salmonpy.examples.mpi_mixed_native \
  STANDARD_BUILD_JSON SERIAL_RUNTIME PBE_BUNDLE NEW_OUTPUT \
  --layers LAYOUT_JSON --field-grid 1 2 1 --tddft-ranks 2 \
  --fast-hse-build SERIAL_HSE_BUILD_JSON --fast-hse-mpi-build HOST_HSE_BUILD_JSON \
  --hse-bundle HSE_BUNDLE --hse-ranks 2 --protect-steps \
  --steps 32 --rank-budget 7 --memory-budget-gib 48
```

This 40³ three-dimensional example is a connection regression, not a converged
material/slab calculation. Practical TDDFT multiscale use still starts in 1D;
SBE/user current models are the recommended electronic models for larger 3D runs.
Select ranks/budgets from the actual allocation, not by copying a larger machine's
values. Local tests cover standard PZ, automatic/fixed TB-mBJ, Libxc-PZ/PBE and
HSE Taylor4/learned PT-CN with internal MPI2, plus the 1D two-TDDFT-group
route. See the validation record for measured scope; no fresh cross-platform CI
or general physical-convergence claim follows from these fixtures.


## Measured larger-grid overhead

On the local 64GB/18-core Apple Silicon calling program, 96³ Maxwell + local LD + SBE +
Libxc-PBE2 + HSE PT-CN2 passed 64-step comparisons with/without protection.
Two repetitions each: field2/total7 took 12.00/14.38 seconds mean (plain/protected),
field4/total9 took 7.08/7.66 seconds. Maximum sums of rank peak RSS were
2.012/2.692 GiB and 2.949/4.417 GiB respectively. All arrays and ordinary streams
matched exactly. These are local measured examples, not memory bounds or a rank
recommendation. Only the Maxwell grid grew; electronic precision is unchanged.
Keep unknown memory allowances unknown and compare on the actual target system.

## Multiple points of the same material

The packaged `examples/mixed_multipoint_layers.json` is a larger connection
fixture: SBE4, Libxc-PBE TDDFT2 and HSE2 independent points, plus 1,152
Drude and 1,152 Lorentz cells owned by FDTD ranks. One electronic cell creates one
state and its own solver group; repeating a material ID reuses its definition,
not its evolving state. Distinct logical point IDs and ownership persist through
checkpoint/restart. This sample uses 64³ Maxwell, not the default 40³ domain.

Use the seven-rank command above with `LAYOUT_JSON` pointing at this file, change
to `mpiexec -n 16`, `--field-grid 2 2 1 --field-shape 64 64 64`, and select
`--steps 64 --rank-budget 16` (or the actual larger allocation). This means four
field ranks, four SBE ranks, two TDDFT groups of two ranks and two HSE groups of
two ranks. LD needs no dedicated groups. The grid/geometry/rank counts are coupled;
do not merely increase layer thickness without recalculating the plan.

For a direct restart, run with `--direct-save-at 31` and then start a fresh MPI job
with `--direct-resume PREVIOUS_OUTPUT/direct-checkpoint`, the same inputs/binaries/
topology and a new output directory. Keep `--steps 64` in both jobs. Only the
remaining 33 solver steps run; accepted propagation is not replayed. Ordinary
output files in the new job describe that segment; full response-array history
is retained by the mixed checkpoint. Do not assume historical RT files are merged.

The local 16-rank fixture passes uninterrupted/protected/direct-resume comparisons
and one-group failure/retry. It is a resource/connection test, not converged
three-dimensional TDDFT physics. Start practical TDDFT in 1D and size from actual
user resources. Custom materials continue to use explicit bindings and contracts.
