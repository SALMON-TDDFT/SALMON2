# Retained GS to standard TDDFT in memory

Goal: create independent A/c→j states directly from retained GS, without a restart
bundle or SCF rerun, preserving the existing RT initialization/update order.

Design: one owned seed per RT model (value-only GS inputs, local orbitals and
occupations/layout), no borrowed GS/MPI/Libxc lifetime. initialization_rt optionally
loads this seed at the old restart-read boundary, then follows the identical
Gram-Schmidt/density/potential/current initialization. Model spawn reuses the
immutable seed. A scalar band's spatial box is redistributed across the existing
orbital communicator when GS empty-band layout differs from RT occupied layout.
Same subgroup/inner decomposition initially; repartitioning is separate work.

Expose NativeGSState.tddft_model(dt,max_steps), using the existing registered
NativeTDDFTModel/state step/checkpoint API on a borrowed GS-session attachment.
GS session cleanup owns all RT states/models before detach; GS close alone never
retires already-created RT models. Before a new GS solve, park any active RT model
so the next RT operation restores its correct input/cache view. Standard PZ/Libxc
scope; no HSE/symmetry/SOC/+U/DC or threaded/concurrent calls.

1. Red public bridge entry test.
2. Seed type/capture/load, owned input conversion and optional RT init boundary.
3. Native GS lookup/parking bridge and model registration/spawn with GS seeds.
4. Python borrowed attachment/model connection and cleanup.
5. Actual MPI1/orbital2/spatial2, real/complex and empty-band GS: multiple states,
   GS/RT interleaving and source-close independence. Existing validated GS fixtures only, no new SCF/reference generation;
   synthetic arrays for empty-band redistribution, no GS bundles collected. Source/object audits, no Actions.
6. Direct restore regression, source/wheel tests, paired guides/sample/AI/checklist.
