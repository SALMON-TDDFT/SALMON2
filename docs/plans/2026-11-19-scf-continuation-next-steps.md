# SCF continuation work after priority change

User order:4→5→3. Fast-HSE norm diagnostics/current-grid campaigns (1/2) deferred;
molecular work paused. No routine CI; GS saving/reuse is authorized by the latest user clarification below. Scientific work
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

## Common direct execution binding

Expose the existing live owner as CheckpointCapability native-scf-memory/1 and
CheckpointTarget callbacks; keep numerical snapshots native and rank-local.
Identity pins model/runtime/live-owner token; compatibility compares actual MPI
layout. Bind restore only with an open same-owner snapshot. Independent candidate
staging remains absent, so protected restore must reject before callbacks rather
than silently mutate originals. Test actual ordinary/HSE MPI continuation and
unchanged-init/work counts, invalid compatibility and protected-entry refusal.

## Protected same-owner numerical candidates

Stage detached s_scf_payload in each live GS owner (optional candidate1 ABI),
leaving numerical originals and FFT/XC/MPI resources untouched. Reserve mutation
until publish/discard. The common coordinator votes after all stages; publish
then assigns only numerical state, with no init/SCF/exchange or files. Mirror the
existing native field/TDDFT two-phase pattern, not a new runnable solver clone.
Handled staging exceptions discard all pending buffers; test original continuation
before retry and uninterrupted parity. Fatal publication allocation/process/MPI
failure remains outside recovery. Keep direct-only old runtime compatibility.

## Reduced-k HSE extension

Own fixed symmetry/density maps in GS contexts; swap them around continuation,
not into portable payloads or shared globals. Opt native retained GS into the
existing explicit reduced-k HSE owner, leaving CLI/default RT policy intact.
Expose explicit immutable symmetry operations as input assets in existing DFT
module; default None unchanged. Compare resumed/failed/retried paths to fixed
continuous same-representation controls; separately report full/reduced fixed-grid
SCF energy/density and unresolved legacy force discrepancy without suppressing
results or changing the physical force algorithm.

The full-k control located a pre-existing Cartesian-force return defect in
sym_vector_force_xyz (the final Amat transform was missing). Fix the one-line
Fortran transformation with orthogonal/identity/C2/skew real-module probes;
retain all pre-fix evidence. Force comparisons to old binaries explicitly account
for old basis; no tolerance relaxation or broad physical certification. Declare
new force basis through optional C ABI; old runtimes warn/retain unknown-basis
values rather than mislabelled physical forces.

## Ordinary-XC symmetric continuation and GS persistence direction

Reuse the already model-owned symmetry/density maps for PZ/Libxc-PZ/PBE; remove only the obsolete nonsymmetry admission gate. Validate pause/continue/direct/protected restore with interleaved inputs against preceding continuous runtime in k/spatial MPI. No numerical iteration changes.

Latest user clarification: saving/reusing GS is the normal SALMON workflow. The earlier no-GS-collection restriction must not block development of persistence. Distinguish converged-GS reuse from full in-progress SCF restart. Next disk work must retain numerical mixing/HSE payloads, rebuild only runtime resources, and prove no replay in a fresh process before advertising it. Portable/repartitioned restart remains separate.

## Implemented disk v1 and residual history

save_scf / run(restart=...) bind a binary/model/topology-pinned Fortran stream format. Zero-iteration bootstrap skips initialization2, guessed orbitals/density and SCF; static resources are prepared. Read/vote/publish uses detached numerical payloads, with hashes/input validation before creation. Native history carries accepted iteration/residual/threshold/pre-stage columns across pause and disk restart; saving adds CSV and optional PNG, with plotting errors warning only. Validate final history-enabled runtimes and wheels before completion. Repartitioned/portable SCF remains pending.

## 2026-11-20 Installed-wheel completion

After explicit user permission to save validation GS datasets, installed-wheel fast-PR HSE C2 kMPI2 fresh-process restart passes. Save3 + load with0 executed SCF + continue5 exactly matches continuous8, including all residual-history rows. A missing last-rank file during native staging leaves existing state intact; corrupted hashed files are rejected before context creation. Source and installed suites:637 run,636pass,1skip; source/installed guide/API/tool/source-byte parity checked. Eight-iteration residual CSV/PNG verified. No CI or commit/push. Same-binary/model/topology scope remains; no portable/repartitioned SCF or physical-convergence certificate. Evidence: calc2026/salmonpy-scf-disk-history-installed-hse-20261119 and salmonpy-scf-disk-history-record-20261119.
