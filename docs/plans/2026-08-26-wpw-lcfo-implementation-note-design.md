# WF+PW Divided-SCF/LCFO Implementation Note Design

## Purpose

Create a Japanese, self-contained LaTeX technical note for future SALMON implementers. The note must remain useful without access to the development conversation and must explain every mathematical transition needed to connect the Kohn--Sham problem, the divided WF+PW basis, fragment SCF, and the final distributed LCFO generalized eigensolve.

## Audience

The primary reader is a computational-physics developer who understands electronic-structure theory and Fortran/MPI but has not worked on this branch. The text must also be readable by a physicist who needs to audit the numerical meaning of the implementation.

## Governing Structure

Use a mathematics-and-architecture-first narrative. Present the physical and numerical construction before naming source files. Put development chronology, failure analysis, test commands, and commit mapping in implementation appendices.

The main text will cover:

1. scope and invariants, including the unchanged DC+LCFO/Wannier90 route;
2. Kohn--Sham equations and the WF+PW divided trial space;
3. windows, overlap regions, quadrature weights, and physical-point identities;
4. fragment generalized eigenproblems and density reconstruction;
5. the divided SCF fixed-point map, mixing, and convergence gate;
6. distributed LCFO Hamiltonian and metric assembly;
7. the single post-SCF generalized eigensolve;
8. the state-machine proof that no post-LCFO density SCF or density gate exists;
9. MPI ownership, redistribution, fingerprints, and memory scaling;
10. source-level implementation mapping and callback contracts;
11. TDD and Si64 validation requirements;
12. observed failures, root causes, current status, and future maintenance checks.

## Mathematical Standard

Every symbol is defined before use. Derivations explicitly distinguish continuous integrals from discrete weighted sums, local fragment indices from global physical-point IDs, ordinary eigenproblems from generalized eigenproblems, and density convergence from orbital residual checks. Assumptions and tolerances are stated adjacent to the equations they control.

## Implementation Evidence

The note will cite repository-relative source and test paths, relevant commit identifiers, route markers, and JSON evidence fields. It will label claims as implemented, unit-tested, integration-tested, or still awaiting successful Si64 completion. Failed validation attempts are retained as diagnostic history rather than presented as successful evidence.

## Artifact and Verification

The final artifact will be a standalone `.tex` file under `docs/notes/`. It will avoid generated build products in version control. Verification will include a LaTeX build when a suitable engine is available, log inspection for undefined references and overfull boxes, and textual checks for all design invariants and current-status caveats.
