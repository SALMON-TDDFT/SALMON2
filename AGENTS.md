# SALMONpy repository

Current direction, decisions and next steps: docs/plans/2026-11-14-handoff.md
(read it first; it records the project lead's decisions and how to verify).

Read salmonpy/AGENTS.md and salmonpy/AI_USAGE.md before changing public APIs.
This repository is for independent SALMONpy distribution. Do not merge it into
SALMON upstream. The project lead authorized bundled solver sources on 2026-11-15.

Bundled solver sources live in runtime/sources; stage-specific complete files live
in runtime/overrides. Runtime preparation verifies manifests and copies to a NEW
external directory; do not apply compatibility patches or modify fixed binaries.
HSE support is the fast PR profile only; comparison baselines are not supported engines.
Preserve shared-kernel byte equality between native/ and the standard solver.
Retain original commits, file inventories, licensing and provenance.

Build kernels explicitly; never compile on import. Use new output directories
outside the repository. Run Python tests and the source-distribution-to-wheel
check after changing packaging or layout. Historical docs/plans and validation
records may reference the old full-SALMON layout; current instructions are in
README.md and runtime/README.md.
