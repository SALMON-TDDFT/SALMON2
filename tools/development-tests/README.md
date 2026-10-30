# Development-only regression checks

These checks support local development and are not public input-file examples.
They are disabled by default and excluded from `git archive` source distributions
through `.gitattributes`. Release tooling that copies the working tree instead of
using `git archive` must also exclude `tools/development-tests/`.

Enable with `cmake -S . -B <build> -DSALMON_ENABLE_DEVELOPMENT_TESTS=ON`.
Run with `ctest --test-dir <build> -R '462_' --output-on-failure`.
The existing compiled probes remain under `testsuites/462_collinear_workers`.

The Python drivers retain the layout oracles, all-rank GS/RT checks, no-PC direct
Taylor4 and learned-ACE comparisons, and unsupported-input rejection checks.
The native no-PC checks compare two owner ranks with 41 balanced ranks at two
OpenMP threads, without changing numerical acceptance thresholds.

Public integration tests should use fixed SALMON input files and ordinary
GS/restart dependencies; do not add these development drivers to that suite.
