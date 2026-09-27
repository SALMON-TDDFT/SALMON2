# Stream complex DC input into mesh tiles

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Eliminate resident fragment basis/coefficient payloads from native DC-to-mesh reconstruction.

**Architecture:** During existing complete-file preflight, validate finite values and record per-spin/k basis/coefficient payload offsets. Native reconstruction opens fragment files for each requested tile/orbital and reads one requested coefficient column plus contiguous selected grid runs via64-bit stream offsets. No grid-times-band or band-times-orbital allocation. Retain existing resident loader for projected LCFO RT because it retains a basis for propagation. Rooted destination reductions remain; per user correction, remove the arbitrary65536-point cap and derive allocations from destination size. Repeated seeks/opens are an I/O tradeoff, not a speed claim.

**Tech Stack:** Fortran stream IO, pure tile maps, MPI native regression.

1. Add failing stream unit oracle with nontrivial offsets, complex coefficients, wrapped/sparse mappings, >4096 grid points and >16 bands. Add native diagnostic expectation and NaN payload rejection test.
2. Add lcfo_mesh_stream module with demand-sized readers and tile contraction. Record preflight offsets in fragment descriptors, strengthen finite preflight, and choose streaming only for nonprojected reconstruction. Propagate read failure collectively before reductions.
3. Run independent oracle, native impulse/pulse/water and malformed-data tests, before/after trajectory comparison, projected LCFO/provenance regression and HSE ON/OFF builds. Independent review, evidence and local commit. Document metadata/map arrays and I/O tradeoffs separately from payload memory.

User correction: fixed caps are not justified by the system or workload. Remove the introduced point/band caps, keep no new input control, and report actual domain-dependent storage rather than a total-memory limit.


## Execution record

Missing stream-module oracle and late rejection of a NaN in an unused saved
orbital were observed before implementation. Implemented64-bit validated
payload offsets, complete finite preflight and native selected-column/run
reads. User correction supersedes the initial cap design: remove the arbitrary
point cap and derive arrays from the domain, retained bands and file runs.
Independent review found no blockers and requested explicit accounting for
index and wire buffers; documentation includes them. Two unit tests including
70,000 contiguous points,16 native tests,8 input tests, both builds and projected
LCFO regression passed. Before/after native outputs agree to file precision or
9.67e-14 eV for water MPI2. No RSS/performance limit or speed claim, no push/merge.
