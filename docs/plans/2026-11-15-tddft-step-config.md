# TDDFT explicit step configuration implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** First native ownership extraction: RT carries its timestep/Taylor/field integration configuration instead of those kernels reading active global input.

**Architecture:** An opt-in external overlay introduces s_rt_step_config under s_rt, initialized once from validated legacy input. Taylor order, propagation mode, field dt/mode/film parameters and core-step diagnostic dt read the supplied RT object's config. Default runtime patch stays unchanged while this extraction is validated. Scalar IEEE/status-return config initialization rejects invalid input without changing an accepted config. Existing arithmetic order remains unchanged.

**Tech Stack:** Small external-source diff, audited source preparation, GNU native config/field-kernel probes, fresh source/object-audited solver build and existing synthetic GS fixtures. No new public native multi-context claim.

## Accepted scope

This is a staged implementation of the approved different-material/state-dt ownership design. FDTD coupling still requires common dt; no subcycling. This stage does not finish parser/XC/communicator/output/error ownership. No source import, real GS collection, remote jobs or fixed executable changes.

## Tasks

1. Add a native acceptance harness that requires the config helper/type, validates two different configs plus finite/NaN/Infinity rejection, and invokes the actual calc_emfields kernel with distinct state dt and poisoned legacy globals. Run against old source and observe missing API failure.
2. Implement overlay on a fresh safe source snapshot. Add config type/helper under structures, capture config in initialization_rt, replace native Taylor/global field config/core-step references. Do not change numerical loop order. Audit source and all probe application objects before linking.
3. Build a fresh serial solver, resolving explicit OpenBLAS/Libxc dependencies from the known reference profile before configure; retain CMake prohibition checks and build only salmon. Compare actual synthetic2/3/6/8/weighted currents at dt.001/.002 against unchanged SHA-validated reference worker, at least26steps. Validate all field modes using the extracted native field-kernel probe and replay parity.
4. Package opt-in overlay/tool and explain limits. Full source/installed suites and sdist→wheel checks where packaging changes. Save reports/logs only under canonical calc2026, retain failed evidence, commit locally. Do not report full native multi-context support.
