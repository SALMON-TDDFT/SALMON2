# Insulating zero-limit verification

Use the vanishing insulating DC conductivity as a constraint, not as a replacement
for the finite dielectric coefficient. Verify the equivalent static density
response: the uniform-density source has zero interband response, the head is
O(q^2), and wings are O(q). The existing optical Sternheimer path already
projects away occupied components. Do not impose this constraint on fractional
occupations, metals, or pumped Drude responses.

Add an independent four-site, gapped spectral reference. Construct finite-q
density vertices from exp(i q.r), and velocity vertices from i[H,r], rather than
constructing finite-q chi from the optical answer. Compare head/q^2 and odd
wing/q against Sternheimer coefficients for several q, both signs and oblique
directions, require quadratic convergence, and test occupied-vertex invariance.
Run the audited native tests and scalar/SOC MPI regressions. No production
correction, new tolerance, or input variable is needed if the existing path
satisfies the constraint. Preserve existing executables and production jobs.

## Verified outcome

The existing production path satisfies the static density-response constraint;
no production correction was necessary. The independent model's contact and
paramagnetic current terms cancel to 2.22e-16. Both oblique directions show
quadratic convergence of normalized finite-q head/wing coefficients. The
occupied-component invariance test and original optical/local-field/IEEE tests
pass. Scalar screening with uneven five-MPI decomposition, spinor screening,
metallic-head and inverse finite/NaN/Infinity tests pass. The tests are registered
through the existing longwave standalone target and hsex_fortran_screening suite.

These tests do not establish actual Si RT DC cancellation, full optical current
normalization, or the Si nonlocal-pseudopotential f-sum. They validate the model
identity and the production optical Sternheimer response used on that model.
No extra acceleration is claimed from the zero condition; the quantitative
q^2 and q coefficients still require the existing response solves.
