# Non-diagonal LD pole energy observation

Diagonal background epsilon (2,3,4), reciprocal pole strength with xy=.25, xz=.1, yz=.2; 40³ grid, 64 steps, no GS. Complete LD component currents include cross response. Background cross terms and subpixel observation remain unavailable. No field update or communication change.

Validated locally: field MPI1/2; observation/no-observation and prior solver final state bytes identical; step32 direct checkpoint/new-job resume gives exact history and final state. Maximum absolute modified energy residual: 2.1963285357085027e-11 hartree. Python suite: 693 tests, one skipped. The sdist-built wheel was installed into a separate environment and imported with Python -I. Its 64-step history, energy observations and final per-rank field checkpoint bytes exactly match the source run. Direct continuation from the source step32 checkpoint in a separate MPI job also matches exactly. The installed-wheel suite ran 693 tests, one skipped. No remote CI was manually run.

Results: canonical calc2026/salmonpy-diagonal-field-energy-20261122/tensor-pole-summary.json. Local OneDrive storage; cloud sync not checked.

Wheel SHA256: `91d8f5782b3918067aef81b82c4e21861f75ae4ca63685b5f9d3b0c14f272644`. Wheel implementation and paired manuals/AI JSON match the checkout byte-for-byte. Dependency packages were provided through the existing validation environment; this checks wheel installation and execution, not offline dependency installation. The field library is the separately built, audited runtime, not a library embedded in the lightweight wheel. Artifacts and logs: `tensor-pole-wheel-artifacts/` under the campaign directory.
