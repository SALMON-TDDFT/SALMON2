# Independent band workflow implementation plan

Goal: execute the user-approved GS snapshot → k path → checked band data → plot workflow.

Architecture: BandInitialState owns a copied, hashed self-contained GS bundle. Each
run copies it into a new output directory and overrides only calculation/control,
band convergence and path settings. BandPath owns fractional vertices, intervals,
and labels. Native SALMON emits a versioned convergence/reciprocal-lattice report;
the adapter requires every block to converge, checks actual Cartesian coordinates,
and removes only verified endpoint padding. Existing solve/update order is retained.

Scope: serial process launch, periodic unpolarized PZ/Libxc-PZ/PBE GS; no MPI launcher
or automatic symmetry path. Files and inputs remain available after failures.
Independent paths cannot write into the original bundle. Plotting remains optional.

1. Add tests for fractional-to-Cartesian transformation in a skew lattice, exact
   batch division and padding, explicit nonconvergence, source snapshot isolation.
2. Add native convergence report without altering numerical iterations; test path
   generation directly using a Fortran harness compiled against minimal type stubs.
3. Implement immutable path, initial snapshot, isolated batch run and checked result.
4. Run actual existing Si GS on two paths and one intentionally unconverged case;
   compare to direct SALMON on the same input, verify original hashes unchanged.
5. Update human topic guide and AI JSON task with implemented API and environment
   requirements. Run unit tests, audits and validate saved example plot.
