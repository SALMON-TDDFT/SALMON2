# LD input lifecycle implementation plan

Goal: make model definition, evolution settings and per-call fields unambiguous.
Architecture: frozen, tuple-valued resolved configuration records; read-only public model/state properties; validated copied input arrays; transactional state update. No kernel or units change.

1. Add behavioral tests for configuration inspection, attempted public mutation, caller-array independence and rejected-step rollback.
2. Add LDModelConfig and LDStateConfig; model and state expose config without setters. Slots reject silently ignored public attributes. Fixed units/method and initial conditions are recorded.
3. Retain E as the only step input, reject invalid scalar inputs and preserve state after validation/native failures.
4. Run all LD Python tests and example; document shared input contract and implemented scope. Restart import/export remains unimplemented, with compatibility requirements documented.
