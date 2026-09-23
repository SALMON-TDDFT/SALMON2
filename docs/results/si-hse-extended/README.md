# Extended HSE response: restart verification and progress

The accepted dt=0.32 a.u. PT-CN-ACE setting is being advanced beyond the short pilot. First stage: zero field for13 steps (0.1006256fs), including a global-step10 MLWF update. Next stage, conditional on that validation: positive impulse A_z=1e-4 a.u. toward130 steps (1.00626fs). This is a stability/response trajectory, not a resolved exciton spectrum.

Restart is separately checked by splitting the previously tested dt=0.16, two-step impulse run after its first step and comparing against the uninterrupted `ptcn_dt016` reference. Both use the same exported SALMON Hamiltonian and all exchange pairs. Short validation jobs may run concurrently on separate cores; their wall times are not new isolated performance benchmarks.

Checkpointing retains the exact localization reference, original energy and full history, with global-step localization cadence. Atomic archive replacement and an output lock protect restart consistency. Tests include interrupted writing, fingerprint mismatch, explicit localization-reference round-trip, injected endpoint-action failure, secondary checkpoint-write failure and initial Hamiltonian-construction failure. See the sample `README_EXTENDED.md` for commands and status semantics.

The propagation remains an independent Python/Libxc/FFTW reference on the native exported grid; native SALMON HSE is not enabled. No k-point convergence study is added.

## Completed validation

Zero field reached 0.10062559 fs (13 steps). Maximum energy change: 7.105e-15 Ha; electron-number error: 4.619e-14; Gram error: 2.389e-12; current norm: 8.694e-12 a.u. The MLWF update before step11 converged in4 iterations. All full residuals stayed below1e-10. This validates this short interval only.

The split dt0.16 restart matches uninterrupted propagation to orbital relative difference9.62e-15, density difference6.55e-15, and current difference8.04e-12. All50 Python tests pass after the accepted-state transaction fix; independent review found no remaining blocker for this scope.

The positive-impulse target is130 total steps (1.00626fs). Current status is in `calculations/si_hse_reference/linear_response/status.json`; authoritative restart/history are in its `restart.npz`. Do not infer completion or a resolved spectrum from this report.
