# Heterogeneous layer groups implementation plan

Goal: explicit material-ID groups in normal layer planning and ResponseSession,
using one native registry per shared solver group without changing default layouts.

Architecture: optional model_groups maps user group IDs to material-ID lists. Each
listed material keeps its binding and decomposition; ranks per group must agree.
All its points use the named group; unlisted materials retain points_per_group.
NativeTDDFTBackend(registry_model=True) opts into registry construction through
open_group. A session-owned resource dictionary shares the attachment across
bindings and releases it after all state/model factories. Custom backends keep
existing prepare/open_group/create_state contracts and can use the group context.
LD stays field-local; HSE retains its separate fast-PR path. Fixed coupled dt.

1. Red planning tests: two IDs one group, rank count and identity, invalid group
   membership and incompatible native mode. Add session-resource ownership tests.
2. Implement planning and opt-in session registry factory; preserve defaults.
3. Real MPI normal ResponseSession vs independent previous-runtime controls,
   interleaved material currents, direct restore and sibling independence. Exercise
   mixed field connection with existing FDTD fixtures where available.
4. Update normal sample, paired guides and AI JSON; run full source/wheel tests.
No native source changes/builds required; verify existing runtime SHA. No Actions
or GS calculation/collection. Store selected evidence in canonical calc2026.

## Verification

Normal mixed native sample with local two-material LD, SBE and two PBE/two TB-mBJ
points: separate10 ranks vs combined4 ranks has exactly equal point IDs and A/E/j
for128 steps. The two TDDFT materials use orbital2 and spatial2 respectively.
Direct save63/resume65 agrees with continuous execution, without replay. A late
last-TB-mBJ candidate failure leaves every original state/output unchanged and
retry succeeds. Same-topology/binary/input restart only. Native runtime is reused
by recorded SHA; no new SALMON build.

Review identified two HSE admission holes: two native HSE owners in one group and
one non-shared HSE owner with multiple points. Both now reject during planning,
with red→green tests. Custom backend/subclass factories remain available.

Final source and wheel suites:628 each, one skip. Installed continuous128 steps
and final-wheel direct resume65 steps both equal the source continuous result
exactly. Paired guides, AI JSON and packaged example group JSON verified. Final
review has no blockers. Selected evidence at canonical
calc2026/salmonpy-layer-registry-20261118, with no GS collection or Actions use.
Local writes verified; cloud synchronization unverified. This stage remains
uncommitted; previous registry stage was pushed as05926ba0.
