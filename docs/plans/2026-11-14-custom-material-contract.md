# Custom material integration: review disposition

Status: next design work; no new API is implemented by this record.

The current `ResponseSession` accepts `step(A)` and adapts legacy
`advance(A_old, A_new)`. Construction remains different: serial
`new_state(dt, time0, initial_A)` versus provider/backend
`prepare`/`create_state` plus execution/reduction declarations. The review's
claim that MPI necessarily requires `advance` is therefore outdated, but its
construction and user-experience concern remains valid.

`MultiscaleField` currently requires the built-in `_Multiscale1DState`, one field
rank, fresh state, and zero initial A at matter cells. A generic user field
cannot yet be substituted there. Do not advertise arbitrary MPI field support.

Next work, after the approved ASE warm-start/XC-force milestone:

1. Define a public adapter from the serial material contract to MPI execution,
   retaining solver-group ownership, stable identities, reduction choice,
   zero-initial-A restriction and failure/close semantics. Do not silently
   assume that a serial Python model is collective or internally parallel.
2. Add paired custom-model guides, a runnable minimal current model, and a
   proposed `check_material` validator. Validate shape, finiteness, declared
   timing/units/sign, independent states, initial values and lifetime cleanup.
   Interface validation is not a test of physical correctness.
3. Design a versioned C ABI before implementation: opaque state handles,
   float64 three-vectors, explicit create/step/destroy, errors without abort,
   buffer ownership, thread/reentrancy limits and MPI communicator ownership.
   Fortran `bind(C)` and C++ wrappers should use the same contract.
4. Extra observables may be optional named outputs. E/B/nonlocal input ports
   require a separate capability and timing design, not an implicit extension
   of the crystalline contract. Crystal TDDFT/SBE remain A/c-driven only.

Transverse prescribed-A operation is the basic bulk workflow. Longitudinal
feedback is an explicit field selection, not the default material behavior.
Acceptance must compare serial and MPI use of the same custom model and retain
existing SALMON update-order and two-point TDDFT regressions.
