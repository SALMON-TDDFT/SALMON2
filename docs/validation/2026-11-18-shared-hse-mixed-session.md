# Shared fast-PR HSE in ordinary mixed MPI sessions

## Change and scope

Normal ResponseSession can own an optional per-binding factory. Native fast-PR
HSE opts in with shared_model=True; explicit points_per_group groups independent
full-k macropoints without changing their IDs, number, cells or material inputs.
Default one-point groups retain legacy behavior and output names. Multiple binary
states use per-point filenames and native-direct/4; versions1–3 remain readable.
No native numerical kernel or fixed binary was modified in this change.

Local environment: Apple Silicon macOS26, Python3.14, Open MPI, CPU64GiB/18 cores;
OMP/BLAS1. Existing explicitly built standard and fast-PR runtimes reused. This is
one-platform implementation regression, not a physical-convergence certificate.
No GitHub Actions were dispatched.

## Numerical campaign

| Case | Independent groups | Shared groups | Steps / saved at | Result |
|---|---:|---:|---:|---|
| LD/SBE/Taylor-learned HSE two points |6 ranks|4 ranks|256 /127|exact fields and currents|
| LD/SBE/Libxc-PBE/PT-CN-learned HSE two points |8 ranks|6 ranks|128 /63|exact fields and currents|

Both HSE points have nonzero, different responses. The normal sample also checks
point responses against independent serial HSE/SBE/PBE references. Direct resumed
A/c, E and j equal continuous results; both HSE energy streams agree exactly after
the split. Only remaining steps execute:129 or65; no accepted-step replay.
Each schema4 archive contains four HSE worker binaries (two points ×two ranks).

For each profile, inject a handled error after the last point's restore candidate
has staged on the last internal rank. Compare every original field/LD binary,
electronic binary, SBE snapshot, Python arrays/clocks/flags, object identity and
output data streams before/after rejection. All stay unchanged. Retry the same
live session successfully, then compare with continuous execution. This tests
handled staged restoration, not MPI process crashes or arbitrary native failures.

Additional successful runs: combined case with two field ranks (7 total), exactly
matching the one-field-rank run; two HSE points separated by an SBE layer sharing
one group; four legacy protected-step regressions (PBE, Taylor, PT-CN, combined),
32 steps with two field ranks, against preceding singleton baselines.

## Distribution checks

Source and installed-wheel suites:602 tests each, one external-condition skip.
Source distribution→wheel and native checksum verification passed. Installed wheel
runs the normal shared MPI sample (4 ranks,128 steps), saves at63 and directly
resumes with the same late-candidate failure/retry check, comparing to continuous.
The two layouts and restoration probe are included in the wheel. Paired guides
and AI catalog are synchronized. Local archives/logs are retained outside Git;
raw native checkpoint workspaces are not claimed as canonical portable archives.

## Remaining boundaries

One immutable non-symmetric full-k model/process; HSE internal MPI1/2. Symmetry
multi-context and heterogeneous same-process GS inputs remain separate work.
Same topology/binary/input required for direct restart. Resource estimates remain
heuristic warnings and never modify inputs. Larger-scale and other-platform
verification, canonical archives and broad failure recovery remain open.
