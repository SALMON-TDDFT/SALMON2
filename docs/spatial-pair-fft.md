# Spatial-pair distribution of symmetry exchange FFTs

Set `SALMON_EXX_FFT_ORBITS=pair` to distribute independent real-space pairs
`(r_a, r_prime)` over the k-point MPI communicator in bounded batches.
This CPU route requires crystal symmetry. Each task retains the complete
k-mesh transform, kernel, phase convention, and FFT normalization.
Empty task owners skip local FFTs while participating in collective transposes.
The default route and `SALMON_EXX_FFT_ORBITS=1` remain available.
Learned ACE is selected independently with `SALMON_FACTOR_HISTORY=learned`.

Validation: 54 exchange-action cases cover MPI 1/2/8, shifted/unshifted
meshes, translations, both FFT layouts, partial batches, empty ranks, and
collective rejection of inconsistent layouts. Maximum relative action error:
1.187e-15. Matched native RT checks cover conventional ACE and learned ACE
through the first post-warmup correction (88 steps).
The local native diagnostic build uses a separate compiler workaround for
pre-existing MATMUL failures; no such workaround is included in this change.
A Fugaku 160-step Nk=4^3 run completed successfully with learned ACE.
