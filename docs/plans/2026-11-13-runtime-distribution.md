# Explicit runtime acquisition

User-approved goal: prefer compatible prebuilt runtimes, with an explicit source
build fallback. Import must perform neither network access nor compilation.

Implement `python -m salmonpy.runtime status/install/build` and Python
`inspect_runtime`/`install_runtime`/`load_runtime` APIs. Select a catalog entry by
OS, architecture, minimum OS/libc, interface revision and required capabilities.
Verify archive and every installed file, reject unsafe archives, install into a
new user-selected directory, and return the existing SalmonExecutable adapter.
Never equate a newly computed hash with scientific validation.

Ship pinned preparation inputs in the wheel so explicit builds work without a
SALMONpy Git checkout. Build requires an upstream checkout, a new work directory,
and explicit compiler/CMake options; audit objects before recording provenance.
No automatic dependency installation, MPI selection or numerical changes.

Add bundle packaging and Linux CI that tests the relocated executable before
producing downloadable artifacts. Catalogs distinguish actual published assets
from CI-only build targets. No claim of Windows/macOS prebuilt availability until
assets have been built, tested and hosted.

Tests: unsupported environment does no network/build; feature/platform mismatch;
checksum failures; unsafe ZIP members; revalidation on reuse; explicit build
command validation; installed-wheel resource presence; native relocated TDDFT.
Update bilingual guide and AI JSON task with selection/fallback examples.
