# Performance-log zero-timer regression

`python3 testsuites/unit_perflog/run.py` compiles the production `timer`,
`perflog`, `math_constants`, and `salmon_global` modules. Test adapters provide
MPI collectives, the clock, and unused checkpoint interfaces. No physics solver
or production binary is executed. Source and every linked application object
are audited before the test executable runs.

GNU Fortran, OpenMP, and MPI are required. The 1/2/4-rank fixtures set timers
through the public timer API and write both human-readable and CSV output.
Floating-point invalid, zero-division, and overflow exceptions are trapped.
Coverage includes all-zero timers, zero-minimum/nonzero-maximum timers,
ordinary balanced/unbalanced ratios, an unrepresentable relative ratio, and
variance for nearly equal positive timings.

The relative column is `max/min`. When the minimum is zero, or the ratio cannot
be represented, it is **N/A** in both output formats. Other measured columns are
retained. Completely zero timer rows remain omitted. CSV consumers must allow
`N/A` in the final relative column; it is not a physical nonfinite value.
