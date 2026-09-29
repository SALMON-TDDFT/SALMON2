# Common distributed k-point exact exchange

**Goal:** Remove the functional-dependent root-only full-support path for ordinary
multi-k hybrid GS/RT before the Si dielectric comparison. User approved this
cleanup on 2026-09-29. Si inputs remain unlaunched and impulse-only.

**Architecture:** Rename hse_exchange/hse_kernel to exx_k_exchange/exx_k_kernel;
retain density-tile MPI transposes and ACE. Kernel construction selects screened
HSE or spherical truncated Coulomb, including its analytic G=0. Full-support
ordinary k-mesh global hybrids use this distributed path. DC/fractional occupations
and finite/local support retain their existing Wannier/spatial paths. Geometry is
an orthorhombic mesh, not a functional restriction. CPU and cuFFT consume the same
kernel and tile layout; accelerator verification still needs NVIDIA hardware.

**Implementation and checks:**
1. Add a small MPI oracle test against retained Wannier full exchange for screened
   and global kernels, shifted/shuffled rectangular k meshes and uneven ranks.
2. Rename the kernel module/API, generalize rectangular grid dimensions, add the
   Coulomb radius and select common routing. Preserve old input spellings.
3. Propagate rectangular dimensions to the owning GPU interface and test oracle.
4. Run MPI1/2/3/4 numerical comparisons, cuFFT CPU-stub tests and normal GS/RT
   regression cases; build MPI and non-MPI; independently review the final diff.
5. Document functional-independent routing and GPU support, then resume the
   approved Si comparison preparation. Do not claim GPU validation on GNU.

**Acceptance:** Exchange action and half-trace energy agree with the old full
support oracle; existing HSE values, Gamma/local/DC behavior, ordinary hybrid
metadata and impulse/Acos2 regression remain valid. No global orbital gather in
ordinary full-support k exchange. No zero-field spectrum work is introduced.

## Verification ledger

- Rectangular shifted/shuffled MPI1/2/3/4 action versus full-support Wannier:
  maximum absolute difference 9.02e-14 (screened), 9.25e-15 (global Coulomb).
- Rectangular/cubic screened/global backend oracle, MPI1/2/3/4, passed; CPU
  cuFFT stub passed. Actual GPU test explicitly skipped: no NVHPC/NVIDIA locally.
- GNU MPI and non-MPI HSE builds passed. Conventional GS/RT 12 tests, input
  validation 20 tests and CTest422-437 (48 stages) passed before final routing
  edge cleanup; final focused verification follows.
- Independent review: keep thermal/extra-state sources on occupation-aware route;
  preserve explicitly requested HSE Wannier snapshots; reject cuFFT for these
  unsupported combinations rather than silently falling back to CPU. Snapshot
  output regression was reproduced before the fix and passed afterward.
- Review clarification: nstate=0 derives occupied count only during RT; ordinary
  GS uses nstate directly. A proposed GS success test with nstate=0 was invalid
  (old and new paths both fail with no states) and was removed. This change does
  not introduce GS default-state support. Routing checks only explicit extra
  states when nstate>0, consistent with the GPU input guard.

- Final conventional GS/RT rerun: 12/12 passed; independent HSE exchange 4/4,
  symmetry 2/2, and MPI uneven/idle-row/OMP regression passed. Final MPI and
  non-MPI HSE builds passed. git diff --check clean.
- Implementation committed as 1bee937f. Si stage-1 runner now snapshots this
  build and starts fresh GS/probes sequentially; no spectrum production has
  started. The initial PBE alias was rejected by Libxc-enabled SALMON; the
  input was corrected to libxc_pbe and the rejected attempt preserved separately.
