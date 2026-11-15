# Molecular self/local-field boundary after one-way tests

Status: proposed next implementation contract, not validated two-way physics.
The one-way refinement shows source-cell fields strongly depend on spatial grid.
Simply connecting sampled total E to the length-gauge worker is not sufficient.

Keep molecular E->p and the conservative p->j map separate. A field-selection
adapter must declare which field is returned to the molecule, with its time,
position, source footprint and subtraction provenance. H/B and higher multipoles
are outside this initial electric-dipole contract.

Candidate first experiment: one fixed molecule in a linear environment, with
an identical vacuum reference grid/source deposition and sample map. Compute
E_drive=E_environment_total - E_vacuum_own_source. The vacuum reference has no
incident source; environmental total E may include prescribed illumination.
This retains environmental scattering feedback while removing the discretized
vacuum self contribution. It also removes the vacuum radiative self term;
therefore call the policy environment_scattered_only and do not claim that
vacuum radiation reaction or spontaneous emission is included. Do not add a
phenomenological damping term without declaring whether it already accounts for
radiative decay.

For more than one molecule, subtract only each molecule's own vacuum source
response. Subtracting the summed vacuum field would incorrectly remove direct
vacuum interactions between molecules. A reusable response kernel assumes a
linear time-invariant reference and must retain its causal time layout.

Alternative: a finite physical polarization footprint plus a documented
local-field/cavity prescription. Refinement must preserve its physical size;
shrinking a one-cell footprint is not the same model. Its relation to the
molecular TDDFT Hartree response must be specified rather than inferred.

Before enabling a public two-way loop:
1. Show vacuum own-source subtraction yields zero feedback (including deposition
   and sampling), while external illumination remains.
2. Show another molecule's vacuum field is not subtracted.
3. Verify causality, time placement, source-work accounting and timestep behavior
   with a linear oscillator/reference environment.
4. Distinguish environment-scattered feedback from full radiative dynamics in
   runtime/response metadata and human/AI guides.
5. Only then connect the experimental molecular TDDFT worker; keep its independent
   state/replay requirements and computational budgets explicit.

No new default feedback policy is selected by the current public APIs.
