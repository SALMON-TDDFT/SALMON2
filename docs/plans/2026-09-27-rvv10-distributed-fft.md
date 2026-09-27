# rVV10 distributed FFT implementation

Approved continuation: remove root-only nonlocal convolution bottleneck, following native spatial-grid conventions. Keep existing branch.

1. Add MPI kernel regression comparing all local energy/vrho/vsigma entries to serial FFTW, with x/y/z and combined domain layouts, nonuniform density, unequal global dimensions, q-channel variations and invalid layout fallback. Observe missing implementation failure.
2. Factor the serial evaluator through an optional convolution callback; leave local spline/functional algebra identical. Add a native FFTE convolution adapter using y/z pencil distribution and x-only gathering. Keep local channel arrays O(nq*Nx*ny*nz); x-only decomposition does not reduce channel storage. Validate 2/3/5 factors, FFTE dimension bound and transpose divisibilities before entry. Own FFT tables so native Poisson tables cannot be invalidated.
3. Connect conventional/LCFO rVV10 to the adapter on compatible grids and retain reference fallback; native gradient/divergence halo semantics remain unchanged. Connect DC total-grid correction if the same halo interface is sufficient, preserving global fragment mapping.
4. Verify direct MPI derivatives and energy against serial, finite-difference potential, existing 23 tests, LCFO integration and HSE-off build. Independent review, record limitations and commit, no push.

Review focus: FFT normalization and global reciprocal indices, communicator rank ordering, identical fallback decision, no overwritten Poisson FFT tables, x replication explicitly documented, collective failure safety, DC energy counted once.
