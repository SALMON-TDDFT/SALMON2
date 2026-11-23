# Explicit packet feedback validation — 2026-11-15

PacketCoupledStepper supports exact registered samples, with explicit per-point
contracts. Native A/H and output j have offset1/2, E offset0; Multiscale A/j
have offset0 and E lag1. All state-local clocks start at zero. Constructor checks
positions/units, availability and current phase. No automatic interpolation,
dipole normalization or magnetic physics. Initial feedback is explicitly zero.
Legacy CoupledStepper and A-only providers remain unchanged.

Apple Silicon/macOS/Python 3.14, one thread per rank, 3–4 ranks, each job <=60s:
- Source MPI: Multiscale A legacy/packet, native A legacy/packet, native E/H
  (one/two field ranks), and wrong-phase refusal for both field types.
- A/E/j histories are exactly equal between legacy and packet A cases, and
  between one/two native field ranks for synthetic E/H current response.
- Native E/H maximum synthetic current: 1.1369560202549155e-5 atomic units.
- Pure tests: wrong clock/position/H availability rejected before field update;
  lagged E and endpoint A admitted; post-advance material error invalidates response.
- Review found sample startup failure could strand material ranks. The sample now
  aborts MPI before context cleanup; missing build-record test exits without timeout.
- Source suite: 423 tests, one skip. Fresh sdist→wheel install succeeds.
- Installed suite: 423 tests, two skips (optional plot and filesystem).
- Installed 3-rank 1D and 4-rank native E/H samples pass; histories exactly equal
  source runs. All tests local, no Actions or physical convergence study.

Logs /tmp/salmonpy-packet-driver-*.log, outputs with corresponding names.
Wheel /tmp/salmonpy-packet-driver-dist-20261115/.
Native H-enabled library: /tmp/salmonpy-h-sample-build-b-20261115/build.json.

This validates scheduling/transport using synthetic current models, not magnetic
material physics or native TDDFT packet adaptation. Existing TDDFT/SBE backends
continue via their established A-driven path. General interpolation/subcycling
and an explicit dipole-current material adapter remain follow-up work.
