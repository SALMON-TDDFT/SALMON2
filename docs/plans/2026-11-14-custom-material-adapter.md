# Custom material contract and MPI adapter implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Let a user implement one serial A-driven initial/state contract and use it directly or through the existing MPI point-response service, with an explicit conformance checker and paired guides.

**Architecture:** Add independent `salmonpy.materials` with `bind_material(initial, *, key, identity)`, `check_material(initial, drive, *, rtol=1e-12, atol=1e-14)` and immutable `MaterialCheck`. Binding creates existing MaterialBinding/provider/backend plumbing for one-rank material groups; it does not parallelize a single serial state. Stable identity is explicitly provided by the model author and must describe algorithm/parameters/data. User code uses existing new_state(dt,time0,initial_A), state.time/current/timing, step(A)->current and optional close. Explicit timing supports endpoint A and j, or declared J_matter converted once to j. Native TDDFT keeps its managed-runtime adapter. Existing collective MPI backends are unchanged.

**Tech Stack:** Python, NumPy, existing timing/crystal/material_response/distributed contracts, optional mpi4py/OpenMPI for tests. No new native solver code or GS collection.

## Design choices and limits

The user approved unifying serial/MPI user models after prioritizing TDDFT/component functionality over structural optimization. A generic wrapper around serial models is selected rather than asking authors to duplicate provider/backend implementations. Replicating a serial model across ranks would duplicate work; initial adapter requires one rank per material group and parallelizes points only. True internally parallel models continue to implement collective group callbacks.

check_material explicitly runs only a caller-supplied representative FieldSamples(kind='A') waveform, creating two independent states. It checks declarations/clocks/shapes/finite values, returned-vs-current consistency, untouched input vectors, independent branches and repeatability, and closes its states. It neither invents physical test amplitudes nor proves the material law, units/sign semantics, failure atomicity, MPI collectives, model immutability or all future inputs. No implicit production-time conformance runs. Adapter failure invalidates its state; no generic rollback is claimed.

## Task1: Contract tests first

Create tests/test_materials.py with a small deterministic declared A-law fixture. Confirm missing-API failures, then implement minimal module. Cover serial/respond parity, explicit sign conversion, wrong clocks/timing, nonfinite/wrong shapes, aliasing, input mutation, deterministic repeats, lifecycle cleanup and explicit stable identity. Keep MPI imports lazy and preserve existing native/LD adapters.

## Task2: MPI integration and shared example

Add examples/custom_material.py with one reusable immutable example model, a main-guarded serial/checker example and an MPI mode using the same binding. Synthetic demonstration only. Add tools/test_material_adapter_mpi.py for 2/3-rank single-rank material groups, multiple independent points per group, serial parity, descriptor mismatch and collective rejection of a failed material step. Do not initialize MPI on module import.

## Task3: Guides/catalog/packaging

Add paired MATERIALS.md/en with minimal authoring contract, serial and MPI usage, identity responsibilities, units/timing, independent ownership, lifecycle/failure limitations and existing collective backend path. Register public names and manual task; install module/guides in wheel. Correct stale guide wording about process-TDDFT disk persistence where relevant.

## Task4: Verification and completion

Run focused/full Python tests, actual 2/3-rank MPI tests with launcher timeout, and compare all point currents to serial results. No native application builds. Rebuild sdist→wheel with required kernel prohibition checks and run installed-package suite and MPI sample outside checkout. Independent review, fix findings, save records under canonical calc2026 and commit without push or changes to existing jobs/binaries.
