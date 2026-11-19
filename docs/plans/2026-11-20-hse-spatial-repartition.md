# Fast-PR HSE spatial SCF direct rank-count restart

Add immutable `hse_spatial=True` to NativeGSBackend, exclusive with hse_orbital.
The existing full-support Gamma collective algorithm is selected on both hosts,
including one rank; inject the same existing selector into the input identity.
Require HSE-spatial compute protocol1 and spatial-repartition protocol2; defaults
and original CLI spatial ownership stay unchanged. Pure spatial, scalar fully
occupied fixed-ion CPU orthogonal-cell Gamma profile only.

Extend the existing spatial donor map in Fortran: dense ACE factor rows,
source/action/work rows and spatial source/previous rows use the old/new global
physical grid boxes and x-fast flattened index. Gauge and occupation/phase/SCF
mixing histories remain replicated numerical state. No accepted ACE evaluation,
eigenstate solve, density evaluation or SCF replay. New static resources and
communication handles are prepared/rebound as for ordinary spatial migration.
Packed/distributed ACE, tiled gauges, sparse localization, multi-k/mixed layouts,
changed-runtime portability and HSE spatial GS→RT are separate.

Validate1↔2, pre-PBE-stage1↔2, changed x/y/z axes and multidimensional2↔4 layouts,
exact physical cut interiors and continued source-initialized tolerance1e-9/1e-10;
late missing donor must preserve original state and sources. Check default
ordinary/k-HSE and owned orbital controls, protected independent owners,
source/installed637 tests, wheel/sample, JP/EN/AI/header and checklist. No new
Python module or routine CI.

## Current status

Both runtime builds and prohibition audits passed. Python637 tests (636pass,1skip)
and a compiled tiny artificial-array row merge/slice/independence/rejection probe
passed. This probe performs no GS/SCF initialization or dataset export.
Actual HSE spatial GS/restart export was automatically rejected as outside the
earlier campaign's data permission. New authorization was requested; no actual
HSE job was executed, rerouted or disguised. Real numerical/wheel acceptance
remains pending. Temporary evidence: /private/tmp/hse-spatial-*.log and
/private/tmp/hse-spatial-synthetic-probe/result.json.

## Acceptance progress after authorization

The user authorized the HSE spatial source/wheel dataset export on2026-11-20.
Actual1↔2, pre-PBE1↔2,16³2↔4 and2×2×1→1×1×2 passed, physical cuts exact,
SCF bootstrap0 and continued source-initialized tolerance unchanged. Late donor
read failure preserved originals and source file hashes. Two different materials
and protected restoration matched same-build independent controls exactly.
Ordinary/k/orbital old-runtime parity passed; ordinary/k-HSE16-step RT exact.
The ordinary control's first rejection was resolved by inspecting the same harness
and prior artifacts: GS output is disabled, arrays stay in RAM and only scalar
logs/records are written. The same unchanged command was then approved.
Installed-wheel acceptance remains to run. The earlier pending entry is historical.

## Completed numerical acceptance

Installed HSE2→1 direct continuation and sample1→2 passed. Source/installed637
tests each636pass/1skip. Native/installed physical cuts exact, zero SCF replay,
continuation at unchanged1e-9/1e-10; independent protected restoration exact.
Prior ordinary/k/orbital runtime controls exact; ordinary/k-HSE16-step RT exact.
Evidence: calc2026/salmonpy-hse-spatial-record-20261120. Paired guides/AI/header
now describe the validated scope; earlier pending entries are historical. No CI.
