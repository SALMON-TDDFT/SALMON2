# Mixed-session identity validation — 2026-11-17

Local Apple Silicon validation only; no remote Actions were used. Resource request: 4 GiB advisory aggregate memory budget, at most four MPI ranks, one OpenMP and BLAS thread per rank. The estimate is not a hard allocation limit or a physical-convergence guarantee.

## Change and acceptance

Ordinary mixed native samples now preserve run/point/material/definition/initial identities and MPI ownership, including LD material-number mapping. Direct envelope `salmonpy.native-direct/3` and replay envelope `salmonpy.native-replay/2` leave numerical payloads and kernel update order unchanged. Explicit configured material records take precedence over fixture-generated fingerprints. Output/temp paths do not define model identity.

- Source suite using an explicitly selected built LD/SBE library: **514 tests, passed, one skip**. Running source tests without selecting a built kernel initially failed as expected; no import-time build was added.
- sdist → wheel build succeeded; wheel installation in the dedicated test environment succeeded. An initial installation lacked permission and therefore tested the previous wheel; that run is not counted. Reinstalled wheel suite: **514 tests, passed, one skip**.
- Source ordinary FDTD + two LD materials + SBE + Libxc-PBE TDDFT (two internal ranks): step 64 save, injected late staging failure at step 66, original native bytes and Python state preserved, old state continues, retry and following three frames exactly equal; zero replayed steps.
- Installed wheel ordinary FDTD + two LD materials + SBE + fast-PR HSE PT-CN (two internal k ranks): same acceptance, exact three-frame continuation, zero replayed steps.
- Both fixtures reject run, point, material, definition, initial ID, initial digest, LD mapping, MPI ownership, and JSON boolean-vs-integer rank mismatch before any native staging/load. A cross-rank current-identity mismatch is rejected before checkpoint-directory creation.
- A fresh installed-wheel CLI restores the source TDDFT checkpoint into another output directory: all eight recorded observables across three continuation frames are exactly equal.
- Legacy native-direct/2 checkpoint restores with explicit `--legacy-checkpoint-identity`, protected restore enabled; saves through step 68. No unidentified-to-identified upgrade is claimed.
- Three-rank two-SBE/two-LD ordinary run writes `identity.json`, result metadata and a readable `point-response` outputs/1 archive. Replay/2 save at step 32 and resume to step 64 reproduces every response.npz array exactly.
- Japanese/English guides and AI JSON guide copies synchronized; whitespace check passed.

## Evidence location

`/Users/otobetoshihito/.codex/outputs/salmonpy-mixed-identity-20261117/`

Main records: `source-tddft/protected-restore-validation.json`, `installed-hse-v2/protected-restore-validation.json`, corresponding logs and complete manifests, `installed-fresh/direct-checkpoint`, `installed-legacy/direct-checkpoint`, `installed-normal`, `replay-save`, `replay-resume`. These paths are local validation provenance, not required installation paths. Build and suite logs are copied alongside the records; distribution SHA256 inventory is recorded there.

## Limits

Same fixed binaries, inputs, time step, run identity and MPI topology are required. Default run label `mixed-example` is for fixtures; experiments should choose explicit labels. Identity metadata is a caller declaration, not authentication or proof of physical convergence. No cross-topology restart, arbitrary backend archive migration, or general transactional MPI-failure recovery is claimed. Numerical solver sources were not changed.
