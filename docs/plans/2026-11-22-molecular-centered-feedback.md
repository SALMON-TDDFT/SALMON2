# Molecular centered feedback implementation plan

**Goal:** Resume the approved single fixed-molecule E-driven two-way FDTD work.
**Architecture:** Fortran midpoint oscillator/work kernel; explicit field trial,
rollback and corrected-current advance using existing whole-step backups. The
oscillator is a development substitute, not molecular TDDFT or a physical model.
Keep total-minus-own-vacuum field selection and one declared initial zero-current
hold. Numerical kernels stay Fortran/C; Python owns MPI scheduling/records only.
**Tech stack:** Fortran iso_c_binding, existing lightweight shared kernel, ctypes,
NativeEHField/whole-step transaction, existing four-rank regression.

1. Resume molecular status in the current checklist; preserve historical pause.
2. Add failing midpoint/work, finite/nonfinite and failure-input-preservation tests.
3. Add a pure Fortran midpoint/work kernel to the lightweight library; private
   binding in existing dipole_source module, no new public TDDFT claim. Replace
   the old Python oscillator arithmetic with the same Fortran kernel.
4. Existing MPI fixture gains explicit centered scheduling: capture both fields,
   predict endpoint fields, compute midpoint response, rollback, advance accepted
   current and compare actual endpoint to predictor. Keep backups until acceptance.
   The held schedule remains an explicitly selected regression control, not a
   compatibility fallback. Label any required nonlinear iteration as unresolved.
5. Validate vacuum cancellation against illumination-only reference and an LD
   environment with the same vacuum source reference. Check source-work and
   oscillator balance; do not claim closed-system energy/flux conservation or
   vacuum radiative reaction. Do not attach isolated TDDFT before this passes.
6. Build/audit lightweight sources/objects, run source/wheel suites and small MPI
   tests. No real GS data or electronic calculation needed for this stage.
7. Update paired radiation guides/AI JSON/checklist. Next stages: replace the
   development oscillator with isolated TDDFT, retain its independent state and
   supported rollback/retry contract, then broaden physical validation as needed.
