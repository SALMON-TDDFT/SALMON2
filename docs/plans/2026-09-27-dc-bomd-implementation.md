# DC-BOMD implementation plan

Approved spec: 2026-09-27-dc-bomd-design.md. Execute inline; final independent review.

1. Add an opt-in static force diagnostic and executable regression: a one-fragment
   and full-cell-buffer DC force must match the conventional analytic force and
   converged DC energy differences. Watch missing input/diagnostic fail first.
2. Preserve global atom IDs and periodic images in init_fragment. Add dc_force
   module computing global electrostatic terms and core-weighted projector
   derivatives. Diagnose frozen-orbital contributions explicitly; never expose
   these as certified MD forces before response checks pass. User control
   yn_dc_force_diagnostic defaults n and requires static PBEh, full MLWF support,
   automatic fragment atom lists, complex orbitals, supported parallel layout.
3. Compare full buffers and truncated buffers at two SCF/displacement tolerances,
   including MPI and dispersion parity. If response is significant, quantify it,
   derive the adjoint/response requirements and keep MD guard intact until solved.
   Numerical derivatives are validation references, not a production MD substitute.
4. Once forces pass, connect total-system velocity Verlet, moving atom maps and
   cache invalidation; reject failed SCF/charge convergence before ionic update.
   Add crossing, dt refinement and NVE conservation tests before enabling DC-MD.
5. Independent review, build ON/OFF, regression tests, record actual achieved
   scope and any unresolved force gate, commit locally on the existing branch.

## Force-gate milestone

Implemented atom identity/image maps and opt-in static frozen-orbital force,
charge and TS diagnostics. Missing-input and missing-contraction tests failed
before implementation;32 tests now pass, including one-fragment24^3 water,
H4 full-buffer derivatives and1/2/4-rank consistency. HSEON/OFF builds pass.
Independent review found no formula/reduction defect; its asymmetric-projector
and thermal-attribution concerns led to separate contraction FD tests and E-TS
audits. Native truncated-oxygen and moving-boundary certification remain open.

Truncated buffers fail the force/energy gate by~0.05–0.08eV/angstrom. MD remains
disabled. User choice of finite-temperature/free-energy versus zero-temperature
occupations is pending before implementing the dependent response equations;
see2026-09-27-dc-force-response.md. This milestone is not DC-MD completion.
