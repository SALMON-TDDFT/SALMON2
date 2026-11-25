# Magnetic spatial self-field feedback

Extend the existing approved total/own-vacuum scheme. Read raw E and time-centered H on their respective Yee component grids in one owning-field-group query. Fortran averages adjacent H half steps; no extra application-wide field broadcasts. Subtract total/own fields in Fortran, stage the selected E/H pair n,n+1, advance quantum n-1->n and return dP/dt plus c*curl((M[n-1]+M[n])/2) exactly once to each field group.

Implementation sequence: failing API geometry/query tests; owned Fortran E/H read and Python binding; extend existing feedback tool with explicit orbital-magnetic mode and padded component geometry; MPI continuous/rejection/restart at matching field/molecule step numbers; compare both E/H and ordinary streams; preserve electric-only results; update bilingual guides/AI JSON/checklist.

Use fixed-ion spinless PZ C2H2 and existing authorized GS only. Trial zero-current drive and accepted current drive must agree under the declared own-vacuum policy. Save both field groups and quantum direct state; after restoration obtain the previous M from the restored quantum state without propagation. This subtracts the vacuum self-field while retaining environment-scattered self-field. It does not claim inclusion of vacuum radiation reaction, spin response, convergence or an arbitrary failure-atomic production coordinator.

### 2026-11-25: orbital magnetic own-vacuum bidirectional feedback

- Done: one owning-field-group E/H query with independent raw H half-step cache, direct persistence and rollback; --orbital-magnetic connects selected total-minus-own-vacuum E/H and one-time electric+magnetic interval current to both groups. No full-grid application-wide replication added.
- Verified:192 molecular/193 field steps; save96/97, process-exit direct continuation and retry64/128 reproduce all histories, observables and six ordinary files exactly. One/two field ranks per group agree; existing electric-only arrays unchanged. Environment return is nonzero in E and H.
- Git-only audited build, source106targeted tests, sdist->wheel, installed106tests and installed6-rank restore/retry passed. See docs/magnetic-self-field-feedback-validation-2026-11-25.md.
- Remaining: arbitrary-layout/failure-atomic production orchestration, spin magnetization, vacuum radiation reaction policy, full charge-continuity/energy and physical convergence checks, performance/scaling. The electric-and-orbital-magnetic own-vacuum response connection is now implemented and validated under the stated conditions.
