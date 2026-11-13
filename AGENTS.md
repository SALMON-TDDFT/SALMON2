# SALMONpy repository

Read salmonpy/AGENTS.md and salmonpy/AI_USAGE.md before changing public APIs.
This repository is for independent SALMONpy distribution. Do not merge it into
SALMON upstream or restore the full SALMON application source tree here.

Bundled native code is limited to native/ LD/SBE kernels. Full solver integration
uses runtime/source.json and the pinned compatibility patch; prepare the runtime
outside this checkout with tools/prepare_salmon_runtime.py. Preserve shared-kernel
byte equality between native/ and the prepared patch. Retain licensing/provenance.

Build kernels explicitly; never compile on import. Use new output directories
outside the repository. Run Python tests and the source-distribution-to-wheel
check after changing packaging or layout. Historical docs/plans and validation
records may reference the old full-SALMON layout; current instructions are in
README.md and runtime/README.md.
