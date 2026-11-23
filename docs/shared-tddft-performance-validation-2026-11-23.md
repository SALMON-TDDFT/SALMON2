# Shared and registered TDDFT performance mode — 2026-11-23

NativeTDDFTBackend execution_mode now supplies model defaults in shared/registry
routes. NativeTDDFTModel.new_state(..., execution_mode=...) fixes the choice per
state; explicit checked/performance siblings can coexist. Registry.open_model
can set its model default, and session registry reuse preserves each binding's
selected default. Defaults do not modify the identity of physical model data.

The new state-token C update uses local generation/bounds checks. Model-selection
and token diagnostic votes are omitted only in this explicit performance route.
Required input/workspace activation, physical communication and evolution remain.
Checked calls retain previous collective validation/rollback. Explicit restore
is checked and retains the receiving mode; automatic failed-update preservation
is not promised in performance mode. No new numerical Python algorithm.

Source and separately installed sdist-built wheel each passed 12 MPI cases:
PZ/TB-mBJ/Libxc-PZ/Libxc-PBE shared states × 1/2 ranks, and heterogeneous registry
pairs PZ↔TB-mBJ and Libxc-PZ↔PBE × 1/2 ranks. Interleaved checked/performance
currents exactly matched independent previous-runtime controls. Siblings retained
their states. Explicit in-memory step8 restoration and steps9–16 continuation
matched; modes were retained. Closed handle updates were rejected. Installed/
source traces were also exactly equal. Tiny pre-existing synthetic GS was read
via symlinks; no GS calculation/copy or disk electronic checkpoint was collected.
This is regression, not a material convergence study or performance benchmark.

Git-only source preparation/new build and source/application-object prohibition
audits passed. Source and installed-wheel suites: 702 tests each, 701 passed,
one optional skip. Wheel code/header, paired guides and AI JSON match source.
No remote CI was manually dispatched. Fast HSE, process workers, molecules and
outer MPI performance mode remain separate work.

Records/artifacts: calc2026/salmonpy-shared-tddft-performance-20261123.
Local OneDrive storage does not confirm cloud synchronization.
Runtime SHA256: `d3a9d279368e5c9a2a95ca12f8aee4186fddb735f35dcdc00bb8966ff5314a31`.
Wheel SHA256: `9be6e6b4decedfaaca1fc00376cc57015a32369e9994ec03d6a01a14ed100f7c`.
