# Different-input fast-PR HSE registry validation

Explicit Python `NativeFastHSEBackend.open_registry` and generation-checked native
register/spawn/close select independent model inputs on one host-MPI attachment.
422 global value fields, unit values, symmetry tables and density-symmetrization
maps are model-owned; numerical/XC/FFT/exchange/MPI context owners are not copied.
Legacy single/same-model entry and native checkpoint format remain unchanged.
Registry/model close stops creation while children continue; last-child close
retires its model, and last model lease retires a closed registry attachment.

Six final native cases compare interleaved models to independent preceding-runtime
controls: full Taylor + reduced PT-CN at MPI1/MPI2/OpenMP2, two different density
symmetry maps at MPI1/MPI2, and different general cells at MPI2. Initial-GS fresh
state recreation and direct save7→24 resume are exact, with no replay. Currents,
energies and local history diagnostics match bitwise. MPI model mismatch rejects
before initialization/spawn. Existing same-model mixed/protected-restore probe
still matches all8 previously stored independent trace controls exactly.

Review identified a missed process-global density symmetry map. Before the fix,
C2z versus identity symmetry inputs reproduce a~1.06e-9 energy discrepancy after
switch/recreation. Owning/swapping that map resolves the original case exactly.
Identity input uses an unchanged existing full-k GS with only input/sym.dat changed;
no SCF or new GS reference was generated. Raw GS bundles were not collected.

The existing Python provider rejects non-atomic-unit inputs before native setup;
that extra attempted case was not counted as a passing numerical case. Atomic-unit
scope remains explicit. Changed MPI topology, arbitrary profile/thread/process
failure handling and physical convergence are not established here.

Scalar IEEE old/new A inquiries cover all6 elements and NaN/+Inf/-Inf corrupted
on one inner-MPI rank: reject collectively, accepted endpoints preserved, later
trace equals the continuous old-runtime control. Other affected finite loops retain
every element. GNU compilation verified; Fujitsu remains untested.

Source620 and final-installed620 tests pass, one skip each. The final installed
MPI2 two-model public sample and MPI2 old-independent/direct-restore comparison
pass. Workspace source/all application-object prohibition audits and wheel native
hash checks pass. No Actions dispatched. JP/EN/AI guides, sample and C header updated.

Reproduce with existing matching GS inputs and explicit local runtimes:

```sh
mpiexec -n 6 python tools/test_hse_registry_mpi.py BUILD OLD_BUILD EXISTING_A EXISTING_B NEW_OUTPUT --size 2
mpiexec -n 2 python -m salmonpy.examples.shared_fast_hse BUILD EXISTING_A NEW_OUTPUT --second-initial EXISTING_B --symmetry
```

Canonical selected logs/JSON, a local development runtime/wheel and explicitly
unapproved intermediate-library evidence are stored in
`calc2026/salmonpy-hse-model-registry-20261118`. This runtime is not a portable
published distribution. Local OneDrive write verified; cloud synchronization
unverified. Next: ordinary heterogeneous material-ID/ResponseSession registry
assembly, broader initialization recovery and platforms/archives remain separate.
