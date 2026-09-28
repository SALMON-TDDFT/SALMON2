# Native exchange orbital distribution

Goal: distribute native HSE/PBEh mesh exchange over both spatial pencils and orbital groups, without retaining an LCFO basis for propagation.

## Batch 1: exchange action kernel

- Add an optional orbital communicator to the spatial exchange application.
- Keep source and target columns local to their orbital group. Broadcast one source column at a time between matching grid pencils; each group computes its own target actions using its spatial FFT communicators.
- Support unequal and empty orbital partitions. Propagate validation/FFT failures before entering subsequent orbital broadcasts.
- Compare against serial MLWF exchange at Coulomb and screened kernels, fractional/zero occupations, phase transport and spatial/orbital decompositions through 8 MPI ranks.

This batch does not enable nproc_ob for native SCF/RT. Its caller must supply already distributed localized sources. No production-memory or speedup claim follows from the kernel oracle.

## Remaining integration

1. Build MLWF overlap matrices from streamed distributed columns, transport the gauge, and rotate into distributed localized sources; retain only small band matrices replicated.
2. Construct distributed ACE factors from the global small metric and streamed exchange actions. Apply distributed factors to local targets, including averaged Taylor stages.
3. Connect native refresh/cache/energy/core exchange with local orbital indices and suitable communicators; make Gaussian seeds independent of orbital ownership.
4. Relax input guards only for the completed route and validate SCF, fixed-ion response/pulse, and PBEh Ehrenfest trajectories against existing results.

No arbitrary memory cap is introduced. Legacy routes remain until their replacement is validated.

## Batch 1 validation

The oracle passes 9 spatial/orbital layouts × 3 screening values × 4 occupation/transport stages, including empty source and target groups. An invalid FFT coordinate in one orbital group tests failure propagation after streaming starts. HSE integration passes all 8 existing tests. A first implementation copied the source even without orbital distribution and altered sensitive SCF convergence; retaining direct source expressions and avoiding that extra buffer on the existing route restores baseline iteration counts (109/466/281), without changing tolerances or iteration limits.
