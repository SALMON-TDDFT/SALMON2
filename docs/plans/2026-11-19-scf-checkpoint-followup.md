# Retained SCF continuation: inspected boundary and next steps

The completed GS/2 path initializes a new GS owner. It is not a native SCF
checkpoint. Keep the current archive and provenance semantics unchanged.

Inspected standard-source boundary:

- src/gs/gs_storage.f90 already owns mixing, CG state, potentials, density,
  orbitals, energy and iteration/residual/threshold scalars.
- src/common/structures.f90 s_mixing contains density/spin and potential
  histories, history length, mixrate, alpha_mb, beta_p and convergence_value_prev.
- src/gs/gs_model.f90 model.run rejects an already active owner; the current
  gs_solve entry performs initialization and is not a continuation entry.
- src/gs/main_dft.f90 has local iteration/loop controls and retires request
  caches before returning retained storage. MPI communicators/XC handles/request
  caches must be reconstructed, never serialized as addresses or shallow copied.
- src/io/checkpoint_restart.f90 checkpoint_gs/restart_gs write/read density
  mixing histories. Their existence does not certify all adaptive/potential/CG
  controls are restored. Audit these fields individually before reusing the path.

Next implementation slices:

1. Separate one explicit SCF iteration continuation context from initialization.
   Preserve the unchanged original entry and arithmetic order. Include every
   scalar and allocated history that the next iteration reads.
2. Add backend-owned serialization into a new versioned checkpoint, recording
   runtime/model/partition and allocation bounds. Recreate resources in a separate
   candidate and publish only after collective validation; preserve original on
   read/allocation/validation failure. Same partition first, no inferred migration.
3. Test synthetic field round trips and failure preservation without new SCF.
   Actual SCF interrupted-versus-continuous validation remains pending applicable
   authorization to generate GS/SCF data. Do not advertise full resume until then.

The prerequisite collective array admission bug is fixed separately: unhashable
quantity arguments formerly escaped before rank agreement; differing valid
quantities could reach different native branches. No numerical solver changed.

## First implementation slice — 2026-11-19

Added s_scf_iteration_progress under s_mixing in the standard source. Previous
rho/Vlocal arrays, electron normalization, no-mixing convergence counter and
convergence flag are now owner-scoped; iteration/residual/rion_update values are
captured at return. No C ABI or call signature changes. Existing fresh entries
reset only this progress; mixing histories/rates are untouched by that reset.

This does not yet save/load progress or continue an active GS. No claims are made
for band esp_old history, DC/molecular continuation, changed partitions or full
SCF interrupted-vs-continuous execution. The stored no-mixing counter alone is
not a resumable no-mixing phase marker. Standard source type layout changed, so
all dependent native objects must be rebuilt together; do not mix old objects.

Verification runner: tools/test_scf_progress.py BUILD ORIGINAL_SOURCE NEW_OUTPUT.
The baseline source is pinned by hash; extracted original/current routines are
compiled against the complete new build. All five convergence criteria, two
spins/electron-count branches, zero/tiny/finite field changes and nontrivial array
bounds are compared exactly. Whole iteration executable body is also compared
after reversing ownership substitutions and excluding final scalar capture.
Synthetic allocations/copy/reset preserve bounds and independent histories.
No new SCF or GS data are generated.
