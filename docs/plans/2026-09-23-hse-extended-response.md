# Extended PT-CN HSE response plan

User accepted dt0.32 and requested continuation toward longer-time linear response.

1. Add atomic self-contained restart archives (orbitals, MLWF gauge, step/dt/field, original energy, history, Hamiltonian/initial-state fingerprints). Test round trip, validation, interrupted-write preservation. Never checkpoint an unaccepted endpoint.
2. Add an extended constant-A driver using the reviewed PT-CN/full-exchange acceptance, same norm gates, energy/current history at every accepted endpoint, MLWF updates at global steps10,20,..., and checkpoints every5 steps/final/failure. Resume rebuilds ACE and validates all physics settings; it does not reapply an impulse.
3. Verify numerical restart parity with a short split run versus uninterrupted prior state, and run zero field to13 steps (~0.1fs). If gates pass, advance positive impulse toward1fs (130steps); keep running status explicit until completion. No spectrum from a sub-fs trajectory. Full spectral/pump–probe work remains later.
4. Independent review, tests, retain local checkpoints and compact results. No k convergence or native hybrid enablement.

Checkpoint/archive and driver transaction tests pass (50 total). Real split-run parity is~1e-14 in orbitals. Zero field passed13steps/0.1006fs including the global-step10 localization update. The next positive-impulse target is130steps/1.006fs, a continuing calculation rather than an already completed spectrum.
