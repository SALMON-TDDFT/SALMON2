# Same-model reduced-k HSE sharing

Continue approved shared HSE ownership work. One immutable material, symmetry
operations, cell and k mesh per process/model. State-specific wavefunction, HSE
map/kernel/history, clock, dt/capacity, output, transaction and restore candidates
remain owned by each context. No heterogeneous symmetry/input switching.

1. Add a numerical regression entering current open_model with reduced fixtures;
   show current explicit rejection. Unit-test symmetry initializer selection.
2. Remove the redundant blanket shared-symmetry veto, retain opt-in/geometry/
   assets/field checks. Advertise a new optional shared-symmetry capability so old
   binaries are rejected before consuming native lifetime. Bind model symmetry
   at first successful init and pass it to subsequent rt_open. Retain pinned
   settings and exact context symmetry guards; remove only singleton prohibition.
3. Explicitly build NEW runtime from verified bundled sources outside checkout.
   Compare shared and isolated states with actual symmetry-preserving drives,
   full/reduced-k, MPI1/2, Taylor strict/learned and PT-CN learned. Save/fresh-load
   and protected late candidate failure/retry; no replay. Full-k regressions.
4. Opt-in normal mixed sample grouping and reduced/full reference route, same
   IDs/topology and warning-only budgets. Paired guides/AI JSON/checklist.
5. Source/installed tests, sdist→wheel/hash, installed numerical smoke, local
   evidence. No Actions. No changes to original fixed binaries or upstream.
