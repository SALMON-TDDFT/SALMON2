# Exchange scratch lifetime reduction

Goal: release temporary allocations at their last use and avoid buffers in mutually exclusive branches. Keep numerical expressions, physical route choice and persistent caches unchanged.

- exx_k_orbits: free forward after orbit ownership, lookup after k permutation, inverse after slot map; build slot before large communication/FFT buffers.
- exx_k_exchange: free inverse after slot construction.
- exx_spatial: allocate source_column only after wave-first return; release multiplier in wave branch after copied kernels are built; release source norm scratch after scalar envelope; allocate global FFT density/spectrum only on actual fallback.
- Verify bounded exchange oracles (HSE and Coulomb, shifted/unshifted multi-rank), SR wave/pair and normal fallback inputs, serial/MPI builds and regression tests. No production runs or performance claims.
