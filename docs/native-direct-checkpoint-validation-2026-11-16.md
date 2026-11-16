# Native direct internal-state checkpoint validation — 2026-11-16

Implemented native `sp_eh_checkpoint(path,restore)` and `sp_tddft_checkpoint(path,restore)` plus an internal same-topology coordinator in the mixed sample. Each rank writes its own stream and JSON state. A single completion manifest is atomically published after all payloads/frames exist. Input/GS/native library/topology hashes and each file checksum are verified before native load. Unsupported backends are refused before checkpoint creation. Native load failures invalidate the job; cross-group rollback is not promised.

## What is restored

- Field: all153 allocated E/H/PML/configuration/observation arrays, dimensions/bounds, iteration, maxima, current injection buffer and step. Static grid/communicator resources are initialized from unchanged input.
- LD: all25 allocated coefficient/history arrays, P,J,old J/E and Q; same material grid/native update schedule.
- TDDFT: both orbital buffers plus scratch, density/local/XC/Hartree potentials and histories, RT arrays and scalar state, Poisson numerical buffers, energy, vector potential, accepted A,current and iteration/parity. Reconstruct nonlocal phase/cache data from loaded A. Native direct profile is periodic unpolarized PZ, transverse middlepoint, no MD/projection/unfold/singlescale.
- SBE: density matrix,A,current and step through validated snapshot payload, with original model/settings retained.
- Coupler: integrated A,E,last charge-current feedback, response input/adapter histories and synchronized clocks.

This is direct dynamic-state loading after static initialization, without propagating from0 to the saved step. It is a native same-binary checkpoint, not a portable interchange archive. Stored observable frames allow complete output comparison; independent serial response replay remains a separate validation calculation, not reconstruction of the resumed trajectory.

## Actual local tests

Apple Silicon/macOS, Python3.14, OpenMPI, OMP/BLAS1, two build workers, at most four MPI ranks. Si starter PZ GS16³/2³k, field40³,dt.003,128steps, two native LD materials and an SBE point. New verified source directory `direct-native-source` and host library `direct-native-build-v4`, using existing OpenBLAS and Release `-O2 -fallow-argument-mismatch`. Failed preliminary builds were retained separately; no existing binaries or prepared source were modified.

| Field ranks | SBE ranks | TDDFT ranks | Saved step | Propagation steps in resumed job | A/c,E,j vs same-topology continuous run |
|---|---|---|---|---|---|
| 1 | 1 | 2 | 64 | 64 | All bit-identical |
| 1 | 1 | 2 | 63 | 65 | All bit-identical |
| 2 along y | 1 | 1 | 64 | 64 | All bit-identical |

Each save job closed and a fresh MPI job/output restored the checkpoint. The odd63 test covers ping-pong buffer parity. The two-field test splits both LD regions. The field1 results also equal the prior pre-checkpoint-library continuous fixture exactly. Different TDDFT decompositions differ at floating-point roundoff and were not used for the strict field2 claim.

Sample's clock/lag checks and independent serial electronic replay pass. TDDFT2 versus serial maximum current difference7.657177843178875e-19; SBE exact. A four-rank64 checkpoint contains8 payload files totaling82,190,142 bytes plus the manifest. This is a measured fixture size, not a general memory/storage estimate or4GiB cap.

Source suite and installed-wheel suite each446 tests,1 skipped,no failures. Four new admission regressions reject unsupported backend, missing completion, changed signature and corrupted payload. Final focused checkpoint/manual/source-manifest tests22 pass. Sdist-to-wheel build and installed MPI execution pass. Japanese/English/AI guides synchronized. Source manifests and base-profile hash references updated, shared kernels unchanged. No GitHub Actions used.

Persistent developer evidence under `~/.codex/outputs/salmonpy-state-inputs-20261116`: `direct-{half64,resumed64,half63,resumed63}`, `direct-field2-{half64,resumed64,continuous}`, matching `.log` files, `direct-{source,installed}-tests.log`, `direct-final-dist`, `direct-native-build-v4/build.json`. User recipes use caller-supplied paths; these evidence paths are not prerequisites.

## Remaining scope

Fresh native text output files are opened in the new work directory; old text-file positions and complete prior native output streams are not resumed. The saved Python observable frames retain the full checked A/E/j history.

General user-backend persistence API, changed rank decomposition, other XC profiles, portable cross-version state conversion, resource budgeting and whole-session atomic rollback remain open. Initial GS/static resource preparation still occurs in a fresh job; elapsed time evolution does not. Physical convergence/production scaling are separate. Molecular work remains paused.
