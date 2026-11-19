# SCF continuation work after priority change

User order:4→5→3. Fast-HSE norm diagnostics/current-grid campaigns (1/2) deferred;
molecular work paused. No GS dataset collection or routine CI. Scientific work
stays Fortran/C, budgets warn, SCF nonconvergence warns and retains.

## Implemented prerequisite

Shared fast-source s_mixing owns prior density/local potential, normalization,
convergence/no-mix progress, iteration/residual and ion-update flag. Fresh calls
explicitly reset only this progress; mixing arrays/rates/config remain separate.
The iteration body is identical after reversing ownership substitutions. Native
assignment is a deep independent copy, not a checkpoint of an entire live GS.

## Next group4 steps

1. Split the existing SCF lifecycle into initialization, accepted iteration
   boundaries and finalization. A continuation call must bypass initialization,
   initial no-mix work and convergence-workspace reset. Hold numerical order fixed.
2. Enumerate and own every loop-carried value, including CG/orbitals/occupations,
   mixing history/rates, density/spin/potentials/Poisson and fast-HSE warmup/operator
   phase/cache. Keep file/MPI/FFT handles outside serialized numerical payloads.
3. Expose same-topology in-memory pause/continue first; compare uninterrupted
   N iterations with K+(N-K), independent siblings and interleaved models. No SCF
   replay or wavefunction-only warm start substitutes for this comparison.
4. Add versioned numerical snapshots and protected candidate restoration through
   the common checkpoint API. Initial validation uses memory/synthetic payloads;
   do not collect GS wavefunction/density datasets under current instructions.
5. Add portable array mapping/changed decomposition only after complete
   same-topology continuation is measured. Existing portable GS initialization
   and same-topology RT restart are separate completed mechanisms.

Do not advertise a public active-SCF resume API before these controls pass.
Broader shared initialization/recovery/native ABI is group5, followed by NLCC
certification/legacy eV-Hartree consistency as group3.

## Fast-PR HSE memory continuation

Retain the same-owner GS numerical payload plus the explicit HSE numerical owner
(ACE/action/occupation/history/phase/support flags). Keep its FFT/backend kernel
live; never copy opaque resources into snapshots. Preserve actual mixing age and
PBE pre-SCF first-iteration/count across segments. Compare full-k kMPI1/2, Gamma,
and a cut crossing the pre-SCF→HSE transition to fixed continuous binaries.
Restore must not initialize GS or evaluate orbitals/exchange. Same-process memory
only; no disk GS datasets, process-restart/repartitioning or SOC/DC expansion.
