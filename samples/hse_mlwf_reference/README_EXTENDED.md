# Resumable constant-A PT-CN-ACE trajectories

The user-selected pilot step is dt=0.32 a.u. This driver extends time, not k sampling. It uses the same full-exchange endpoint acceptance (1e-10), inexact-inner stagnation handling, and electron-number/Gram gates as the short PT-CN tests. Energy and current are recorded at every accepted endpoint. No normalization or density rescaling is applied.

Start from the converged all-pair HSE ground state:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 samples/hse_mlwf_reference/propagate_ptcn.py calculations/si_hse_reference/export calculations/si_hse_reference/scf/state.npz calculations/si_hse_reference/linear_response 130 --dt 0.32 --amplitude 0.0001
```

To extend or resume, use the same arguments plus `--resume` and a target step larger than the saved step. `target_steps` means the **total** number since the initial state, not an additional number. Reapplying `set_field` sets the same constant A; it does not reset the wavefunctions or apply a new kick.

`restart.npz` is the authoritative saved state. It embeds wavefunctions, gauge, the exact previous localized basis used to initialize MLWF minimization, initial energy, global step, history, dt/field, and SHA256 fingerprints of the exported Hamiltonian and original ground state. Resume rejects mismatched physics. ACE is rebuilt rather than serialized. MLWF updates use the global step modulo10, preserving cadence across resume.

The archive is written to a temporary file, flushed/fsynced, then atomically replaced every5 steps (configurable), at completion, and on handled failure. Atomic replacement protects against partial archive writes; it is not a guarantee against all power-loss/filesystem failures. An advisory file lock disallows simultaneous writers. Existing trajectories require explicit `--resume`.

Only a candidate that passes solver, norm and observable checks enters the single accepted-state snapshot. Error handling saves that coherent snapshot and preserves the original error even if checkpoint writing also fails. Initial Hamiltonian construction failures are recorded as failed status as well.

`status.json` describes the current run. `result.json` describes the last completed requested segment and may remain from an earlier target during a resumed run. `trajectory.json` is a convenience copy and may be ahead of the last periodic checkpoint after a hard kill. On resume the archive's embedded history is authoritative. Large arrays stay local under the ignored calculation directory; compact verification reports are stored in docs/results.

This driver currently handles constant A after an impulse only. It does not implement a laser waveform, nor does reaching1fs constitute a resolved exciton spectrum. The full impulse spectrum and subsequent pump–probe comparison need a longer validated time record.
