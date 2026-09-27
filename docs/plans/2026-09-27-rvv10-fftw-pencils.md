# Cached FFTW pencil transforms

User approved: retain y/z decomposition, replace local FFTs with cached FFTW plan_many, compare warm execution against FFTE before selecting the native default.

1. Add a direct distributed FFT test against FFTE for multiple channels, forward/inverse normalization, xyz modes and mixed rank layouts; verify cached reuse.
2. Implement x/y/z pencil redistribution using existing axis all-to-all communicators. Cache measured FFTW plans and packing maps by grid, process coordinates and batch size; batch at most four channels to bound workspace. Use serial FFTW within each MPI rank initially, and compare timing at OMP=1 fairly.
3. Connect the rVV10 callback with optional backend selection for tests; preserve the current supported-grid/fallback contract. Keep FFTE as reference. Benchmark plan creation separately from warm FFT, communication/packing and complete rVV10 evaluation.
4. Run direct scalar/full-potential, MPI/native regressions and build checks; independent review. Record results and practical limitations; commit without push.
