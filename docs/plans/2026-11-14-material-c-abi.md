# Public material C ABI implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Provide one versioned serial endpoint-A material ABI consumable by Python or an external C/C++/Fortran host, with independent model/state handles and a tested reference plugin.

**Architecture:** A public C header declares ABI1 query/table and callbacks for model_create/destroy and state_create/observe/step/destroy. Opaque positive int64 handles are owned by the library; stale/invalid handles fail, no pointers/arrays are retained across calls. Each state owns its full mutable history. Input is Ac=A/c and output is physical charge-current density j, both three doubles in atomic units at endpoints. Fixed dt/t0, initial A, declared input lag/uses_t0. Every failure leaves accepted state unchanged. Status codes are returned rather than crossing language boundaries with exceptions or process termination.

**Tech Stack:** C99 header/reference plugin, C++/Fortran ISO_C_BINDING host smoke tests, Python ctypes+NumPy wrapper, existing bind_material/check_material. No full SALMON source import, GS collection or change to TDDFT globals.

## Scope and decisions

The user asked to continue after the serial/MPI adapter and previously requested public C ABI for plugging custom models in and calling components from external hosts. This implements that common A-material boundary first. A callback table selected through version+struct-size handshake avoids guessing exported function signatures. Model/state tokens separate read-only setup and mutable evolution, with no fixed state-count cap. Configuration is opaque caller-provided bytes interpreted by the plugin, whose documentation supplies the schema. Reference plugin uses a small explicitly documented two-double config. Native TDDFT remains one lifetime/process; multiple native TDDFT contexts, LD/SBE native-handle migrations, collective MPI ABI, disk/native checkpoints and GPU are separate work.

Python CMaterial verifies the caller-supplied library SHA before loading, pins a loaded path to its original SHA, copies configuration, keeps the library/model alive while states exist, and maps the common model/new_state/step contract to callbacks. Library calls are serialized; no MPI/compilation on import. It validates finite inputs before native calls and finite clock/current outputs before committing cached state. A success response violating the ABI invalidates that Python state; no generic native rollback can be invented. Nonzero callbacks obey the ABI atomic-failure contract. Model close prevents new states and defers native destruction until live states close.

## Task1: Header and contract tests

Create include/salmonpy_material.h packaged under salmonpy/include. Add an explicit standalone reference-plugin build/check tool and native source under examples/c_material. Audit source before build and every newly linked application object. Unit tests first cover missing wrapper, hash/version/size/callback rejection, independent handles, lifetime, atomic bad inputs, failed callbacks, stale handles, and serial/checker/bind adapter parity.

## Task2: Reference plugin and Python model

Implement dynamically allocated handle registries with monotonic IDs, allocation checks, scalar finite checks, overflow-safe clock/counter arithmetic, copied configuration, no fixed capacity. Implement CMaterial/CABIError and owned state adapter in c_material.py. Do not change existing Python LD/SBE implementation or historical executables.

## Task3: Cross-language verification

Compile C reference plugin and C++ and Fortran hosts in a new temporary build directory. Compile/link objects in explicit audited stages. Both hosts query ABI1, create independent states, interleave different histories, reject invalid/stale handles and compare to the same manufactured recurrence. Python runs check_material and existing MPI point binding on the plugin. These are interface tests, not physical validation.

## Task4: Docs/distribution

Paired C_ABI.md/en specify exact ownership/layout/status/timing/error contracts and supported scope. Package header, ISO declarations and example sources. Register public API/manual. Run full source tests, sdist→wheel with retained prohibition checks, installed-package/host/plugin smoke tests. Review, record canonical calc2026 results and commit without push or altering existing jobs.
