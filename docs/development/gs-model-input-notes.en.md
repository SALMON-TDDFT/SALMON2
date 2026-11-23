Implementation notes recorded 2026-11-24. Current supported Python interfaces: [retained GS](../../salmonpy/NATIVE_GS.en.md).

## Developer-only standard GS model inputs

The standard runtime adds internal Fortran type salmonpy_gs_model::s_gs_model.
load(path,status) retains a model-owned parser and post-coordinate input values;
run(status) activates them for sequential legacy SCF, and clear() releases input
and scratch. A rejected unsupported replacement preserves the existing model.

Scope: periodic full-k unpolarized fixed-ion PZ/TB-mBJ/Libxc-PZ/PBE. The calling program first
initializes MPI, comm_attach, math/parallel infrastructure and timers. All subgroup
ranks call the same model in matching order. Runs switch/restore caller input and
unit views, calc_mode, actual symmetry switches, modulo tables and XC scratch.
The model-only call retires persistent communication caches before local buffers
expire and releases scoped child communicators after each call. The default run(status) call retains no GS orbitals. Optional retained storage
and Python handles are now available through [retained GS](NATIVE_GS.en.md).
Simultaneous GS/TDDFT, threads and HSE/SOC/+U/DC GS remain separate work.

Status0 means returned normally,1 means empty model/occupied child scope,2 means
unsupported conditions for this internal entry. SCF nonconvergence remains a
warning and outputs the last iterate. General parser/STOP/allocation failure
recovery is not promised. Relative references use the caller's working directory.
Python DFT.run continues to use its existing separate-process execution path.

Developer probe: tools/test_standard_gs_models.py NEW_BUILD OLD_BUILD NEW_OUTPUT
--reference-sha256 SHA --size 2 --split spatial. Supply explicitly audited matching
runtimes; the tool audits probe sources and all linked application objects before
use. It compares A→B→A with independent old-executable controls for8³/12³ grids,
different cell/atom/electron/XC inputs, checking energies, residuals, iterations,
eigenvalues, k output, restored symmetry and repeated communicator counts.
It sets write_gs_restart_data='no', generating/collecting no GS restart bundles.
Four-iteration diagnostic inputs establish implementation parity, not converged
material properties.

