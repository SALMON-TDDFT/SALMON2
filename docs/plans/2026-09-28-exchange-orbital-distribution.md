# Native exchange orbital distribution

Goal: distribute native HSE/PBEh mesh exchange over both spatial pencils and orbital groups, without retaining an LCFO basis for propagation.

## Batch 1: exchange action kernel

- Add an optional orbital communicator to the spatial exchange application.
- Keep source and target columns local to their orbital group. Broadcast one source column at a time between matching grid pencils; each group computes its own target actions using its spatial FFT communicators.
- Support unequal and empty orbital partitions. Propagate validation/FFT failures before entering subsequent orbital broadcasts.
- Compare against serial MLWF exchange at Coulomb and screened kernels, fractional/zero occupations, phase transport and spatial/orbital decompositions through 8 MPI ranks.

Batch 1 alone did not enable nproc_ob for native SCF/RT. Its caller must supply already distributed localized sources. No production-memory or speedup claim follows from the kernel oracle.

## Integration (implemented)

1. Build MLWF overlap matrices from streamed distributed columns, transport the gauge, and rotate into distributed localized sources; retain only small band matrices replicated.
2. Construct distributed ACE factors from the global small metric and streamed exchange actions. Apply distributed factors to local targets, including averaged Taylor stages.
3. Connect native refresh/cache/energy/core exchange with local orbital indices and suitable communicators; make Gaussian seeds independent of orbital ownership.
4. Relax input guards only for the completed route and validate SCF, fixed-ion response/pulse, and PBEh Ehrenfest trajectories against existing results.

No arbitrary memory cap is introduced. Legacy routes remain until their replacement is validated.

## Batch 1 validation

The oracle passes 9 spatial/orbital layouts × 3 screening values × 4 occupation/transport stages, including empty source and target groups. An invalid FFT coordinate in one orbital group tests failure propagation after streaming starts. HSE integration passes all 8 existing tests. A first implementation copied the source even without orbital distribution and altered sensitive SCF convergence; retaining direct source expressions and avoiding that extra buffer on the existing route restores baseline iteration counts (109/466/281), without changing tolerances or iteration limits.

## Integrated native route

Implemented streamed small overlaps, temporal polar gauge transport and weighted source rotations; distributed ACE metric construction, factor rotations/application and local midpoint concatenation; native orbital offsets, collective cache invalidation and exchange-energy reductions. Gaussian initialization uses the same global centers across orbital groups. The original non-orbital arithmetic path remains separate to preserve its SCF behavior.

Admission covers Gamma full-support HSE/PBEh conventional and DC SCF, plus the existing DC-initialized native mesh RT routes (fixed-ion HSE, fixed-occupation PBEh Ehrenfest). Conventional BOMD and legacy finite-support/multi-k/projected routes retain previous restrictions. Although the internal kernels support empty partitions, native input/runtime guards reject orbital groups with no states because SALMON propagation skips such groups; native uneven nonempty partitions are tested.

Mesh arrays are local in grid rows and orbital columns; gauge, raw overlap and ACE metric matrices remain replicated quadratic band arrays. This does not distribute the LCFO global eigensolver and does not establish production speedup or total-process memory scaling.

Validation adds SCF/impulse/pulse comparisons through 8 ranks for HSE/PBEh, fractional DC through 16 ranks with 6 states across 4 orbital groups, PBEh Ehrenfest pulse trajectories, and water with 4 states over 3 orbital groups. Compare final energy/eigenvalues/charge and trajectory outputs at unchanged tolerances. Independent review identified the native empty-rank collective hazard; guards and regression close it.

Final verification: HSE-enabled and HSE-disabled builds passed. The exchange/MLWF/ACE oracle passed 108 cases (maximum exchange difference 5.52e-15); `python3 -m unittest test_hse_spatial test_ehrenfest -v` passed all 30 tests. Native empty-rank rejection and internal collective error propagation were exercised. Tests do not measure large-system runtime or peak resident memory.
