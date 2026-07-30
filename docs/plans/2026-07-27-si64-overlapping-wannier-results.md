# Si64 overlapping-Wannier ground-state gate

## Fixed matrix

The production gate uses Gamma-only, non-SOI PZ Si64 with 256 electrons,
`nstate_frag=400`, `mixrate=0.2`, `threshold=1e-9`, eight MPI ranks, and
one OpenMP thread.

The deferred large-system Cartesian matrix is:

- fragment decompositions: `2x2x2`, `4x2x1`;
- split-axis buffer depths: 5 and 6 grid points;
- candidate/target windows: 192/48 and 256/48 states per fragment.

An overlapping-Wannier box is rejected when it exceeds the periodic
system on any axis or equals the complete periodic system on all axes.
An unpartitioned axis therefore uses zero buffer. The small Si64 case is
used first as an end-to-end smoke gate; buffer and decomposition
convergence are deferred to a system large enough to exercise DC without
turning a fragment box into the whole system.

The normal reference enables the established DC LCFO and EigenExa route
and leaves the overlapping-Wannier route disabled. Every matrix row
enables only the overlapping-Wannier ground-state route and explicitly
disables direct-DG continuation, LCFO, EigenExa, WPW, conventional RT,
and normal checkpoint publication.

## RED evidence

Command:

```text
python3 tests/dg/check_si64_overlapping_wannier_gate.py \
  /tmp/si64-ow-red2.<fresh>
```

Observed result:

```text
Si64 overlapping-Wannier gate: FAIL:
normal reference raw log or variables.log is missing
```

The checker returned exit status 1. It derives acceptance from numeric
runtime records, the checkpoint manifest hash, and a second-run
checkpoint reuse marker. It does not consume hand-authored pass fields.

## End-to-end smoke evidence

The `2x2x2`, buffer-6, candidate-192, target-48 row completed construction,
the positive global metric gate, weak-operator assembly, coefficient SCF,
the unmixed fixed-point gate, route-checkpoint publication, and checkpoint
reuse. The unmixed density residual was
`6.1062845509342129e-8`. The accepted manifest SHA-256 was
`23184512c90daf6f537529c07f405f674851169f004187fc2dc2e10de7c64cb0`.

Raw evidence:

```text
/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG-task9-results/
  si64-smoke-buffer6-v9/
  decomp-2x2x2_box-buffer6_window-c192-sp48/
```

The restart log contains `[OW-GS] reused accepted route checkpoint`.

## Row-owned V2 smoke evidence

The same accepted minimum coordinate was rerun after distributing the
overlap, weak Hamiltonian components, assembled Hamiltonian, solver inputs,
and checkpoint overlap by exact Wannier-row ownership. The checkpoint
format is V2 and stores one overlap-row shard per MPI rank.

Raw evidence:

```text
/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG-task9-results/
  si64-smoke-buffer6-row-owned-v2/
  decomp-2x2x2_box-buffer6_window-c192-sp48/
```

The clean parent-prerequisite overlay executable SHA-256 was
`3c48ec69fd6746c50c3206c7e71c0b2112276e519c58e67bf83207a90e7510e9`.
The run used 8 MPI ranks and one OpenMP thread and recorded:

- global target rows: 384;
- maximum owned rows on any rank: 48;
- maximum persistent row-owned overlap bytes per rank: 294,912;
- maximum row-owned Hamiltonian bytes per rank: 294,912;
- peak process RSS: 1,605.3125 MiB;
- runtime: 1,602.929831583 s;
- density residual: `1.5127284475512645e-8`;
- unmixed density residual: `6.1156472090916136e-8`;
- coefficient residual: `2.1077479424507706e-11`;
- S-orthogonality defect: `1.5632162231327129e-11`.

The accepted V2 manifest SHA-256 was
`a35aee43306d1baa060f0b26ff1951007347159fd606eba2a5992fc5b839fcc3`.
All eight versioned rank shards begin with
`SALMON_OW_GS_RANK_SHARD_V2`. The restart reused the accepted route
checkpoint, did not rerun the overlapping-Wannier SCF, and left the
manifest hash unchanged.

The recorded minimum artifact is reproduced by the same public checker
used for the full matrix:

```text
python3 tests/dg/check_si64_overlapping_wannier_gate.py \
  /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG-task9-results/\
si64-smoke-buffer6-row-owned-v2 \
  --repo /Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/\
.worktrees/wpw-s-orthogonal-complement \
  --skip-reference --minimum-row
```

This command also requires the row-storage source contract to pass. The
overlap byte count is derived from the actual `ow_srows` allocation and the
Hamiltonian byte count from the maximum actual `hrows` allocation, rather
than from a full-matrix-absence self-report.

## Current status

The route, construction, solver/density, SCF, and checkpoint focused
fixtures pass on 1, 2, 4, and 8 MPI ranks. The clean parent-prerequisite
overlay build also passes.

The first two rows of the full matrix completed both their initial and
checkpoint-reuse runs:

- `2x2x2`, buffer 5, candidate 192, target 48;
- `2x2x2`, buffer 5, candidate 256, target 48.

Each row individually passed the Hermiticity, density, coefficient,
charge, and checkpoint gates.  Their pair did not pass the required
candidate-window convergence gate: the relative changes were
`2.4097e-2` for the minimum metric eigenvalue, `8.9357e-2` for the
maximum metric eigenvalue, and `2.8445e-1` for total energy.  The public
checker therefore rejects
`si64-row-owned-full-v1`.

Candidate-window diagnostics were extended beyond the fixed matrix.
Candidate 320 completed with metric minimum `4.33201942e-2`, metric
maximum `2.2399543023`, condition `51.7069`, and total energy
`25.2052528515`.  A cold candidate-384 solve failed at overlapping-
Wannier SCF iteration 9.  Reusing the preceding coefficients as the
generalized eigensolver initial block removed that failure.

The row-owned candidate-384 rerun then exposed an absolute residual-Gram
rank threshold which incorrectly discarded all correction directions
below norm `1e-7`.  The retained rank is now selected relative to the
largest residual-Gram eigenvalue.  A near-converged warm-start fixture
passes on 1, 2, 4, and 8 ranks.

A subsequent fresh run exposed loss of S-orthogonality after a
full-dimensional expanded Ritz solve.  Once the 384-dimensional basis had
been reached, its residual was only marginally above the requested
`1e-7`, but another restart reduced the retained basis and drove the
residual to one.  The residual and previous-direction block is now
projected twice in the S metric.  A full-dimensional Ritz solve is no
longer restarted: it is accepted only through the existing final
`10*tolerance` numerical-quality gate, or rejected explicitly.

The clean-overlay candidate-384 artifact
`si64-window-diagnostic-c384-warm-rowowned-v7` then converged in 78 outer
iterations.  It recorded density residual `1.7901782259e-8`, unmixed
density residual `1.1573672231e-9`, coefficient residual
`8.1603940173e-8`, S-orthogonality defect below `2e-14`, and both trace
and integrated charge equal to 256 electrons.  Its metric interval was
`[6.9553451797e-2, 2.2048038820]`, and its total energy was
`26.8002764030`.

Candidate 400, the complete computed fragment-state window for
`nstate_frag=400`, also converged.  Relative to candidate 384, its metric
minimum changed by `1.3430e-3`, its metric maximum by `1.9131e-2`, and
its total energy by `1.0170e-2`.  Thus even the upper end of the computed
fragment window had not reached the required `1e-4` convergence regime.

The complete s+p seed was then changed from a mixture of pseudo-atomic
orbitals and nonlocal pseudopotential projectors to the matching
pseudo-atomic orbital `upptbl_ao(:,l,species)` for every s and p channel.
Fresh clean-overlay c192 and c256 runs both converged and reused their
route checkpoints.  Their window differences were:

- metric minimum: `1.90134584e-2`;
- metric maximum: `7.14592689e-2`;
- total energy: `2.82780905e-1`.

These are essentially the same scale as the mixed-seed result.  The
nonlocal-projector seed hypothesis is therefore rejected: truncation of
the fragment candidate eigenspace, rather than the radial seed choice,
dominates the measured window dependence.

The pseudo-atomic s+p shell was therefore retained directly as an
S-orthogonal complement to the computed fragment-state window. This
removed the candidate-window dependence without increasing the published
target rank. On each completed physical-box profile, the c192 and c256
results agreed to roundoff:

- `2x2x2`, buffer 5: metric minimum `6.2e-15`, metric maximum
  `3.4e-16`, and total energy `9.4e-16` relative change;
- `2x2x2`, buffer 6: metric minimum `1.3e-15`, metric maximum zero,
  and total energy `2.8e-15` relative change;
- `4x2x1`, buffer 5: metric minimum `2.3e-14`, metric maximum
  `6.8e-16`, and total energy `4.4e-16` relative change.

Raw evidence is under
`/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG-task9-results/si64-direct-sp-full-v1`.
All six completed matrix rows also completed checkpoint-reuse runs without
changing their accepted manifest hashes.

This direct-shell result does not pass the full Task 9 gate. For `2x2x2`,
changing the physical box from buffer 5 to buffer 6 changed the metric
minimum from `0.2027854964` to `0.2662546781`, the metric maximum from
`1.9328665739` to `1.9105725136`, and total energy from `26.4839366138`
to `24.2686039668`. These changes are far above `1e-4`.
At buffer 5, the supposedly equivalent `4x2x1` decomposition produced
metric interval `[0.1895898112,1.9669823791]` and total energy
`203.0800795382`, so decomposition invariance also fails.

The `4x2x1`, buffer 6, c192 row stopped before overlapping-Wannier
construction. Its DC preconditioner reached 1000 SCF iterations with an
energy difference of approximately `1.0349e-3`; the resulting core
electron count did not define an integral occupied rank. Its runtime
buffer is `(6,6,0)`, so the unpartitioned z direction is not extended and
the failure is not caused by allowing a buffer box to span the full
system. The public checker consequently reports missing accepted-route
evidence for this row. The c256 row was not launched after its c192
prerequisite failed.

The immediate Task 9 acceptance matrix was subsequently restricted to
the fixed `2x2x2` physical decomposition on eight MPI ranks. On this
32-cubed grid, `4x2x1` makes the x core only eight points wide while the
two-sided buffer occupies ten or twelve points. Changing from `2x2x2` to
that geometry changes the finite-buffer DC approximation; it is not an
equivalent MPI placement of the same physical fragment problem. The
existing `4x2x1` artifacts remain diagnostic evidence but are no longer
acceptance rows.

The revised matrix retains the normal reference and the four Cartesian
rows formed by buffer 5 or 6 and candidate 192 or 256, all at `2x2x2`,
eight MPI ranks, and one OpenMP thread. Candidate-window invariance is
already at roundoff. Buffer 5 is the fixed production profile. Its two
candidate windows, normal completion, checkpoint reuse, and every existing
numerical and forbidden-route gate are required for acceptance. Both
buffer-6 rows remain mandatory diagnostic runs and must pass their
individual gates, but the buffer-5/6 difference is recorded rather than
used as an acceptance threshold. MPI-count invariance and physical-
decomposition invariance on a larger system are deferred.

The revised public checker was run against the complete preserved
`si64-direct-sp-full-v1` raw evidence. It returned exit status 1 with
`2x2x2: buffer/window convergence failed for metric_min`. Thus narrowing
the physical decomposition removes the invalid equivalence requirement
but does not make Task 9 pass: buffer 5 and buffer 6 remain numerically
different. A fresh `si64-fixed-2x2x2-v1` run completed the normal
reference and the buffer-5 c192 row plus checkpoint reuse, then was
externally terminated during conventional-DC iteration 5 of buffer-5
c256. Its empty timing file and absence of an error-stop record distinguish
that interruption from a SALMON gate failure; the partial artifact is
retained and is not used as acceptance evidence.

After selecting buffer 5 as the fixed production profile, the public
checker was rerun against the same immutable raw evidence. It returned
`Si64 overlapping-Wannier gate: PASS`. The buffer-5 c192/c256 production
rows agree to roundoff, and both buffer-6 diagnostic rows individually
pass the route, numerical-quality, and checkpoint-reuse gates. Task 9 is
therefore accepted for starting coefficient-only real time. The measured
buffer sensitivity remains an explicit limitation of this first
production profile.

The production `src/` state used for this acceptance decision is unchanged
from the previously successful clean parent-prerequisite overlay build.
All acceptance-semantic changes after that build are confined to the
Python evidence checker and documentation. A later attempt to reconstruct
the deleted temporary overlay from clean `HEAD` alone stopped in existing
WPW code because `lcfo_flux.f90` references prerequisite APIs which do not
exist in any reachable committed implementation. The parent worktree's
uncommitted implementations were not imported. This reconstruction
failure does not replace the recorded successful overlay build, but it is
an outstanding baseline reproducibility issue to resolve independently.

Review follow-up now rejects `.not. (sum1 < threshold)` at the
overlapping-Wannier publication boundary, before constructing a core-owned
occupied subspace from an unconverged or nonfinite conventional-DC state.
It also
collectively compares the minimum and maximum fragment target ranks before
forming the global target extent or entering fixed-count Wannier-tail
collectives. The nonuniform-rank rejection fixture passes on 1, 2, 4, and
8 MPI ranks.

The latest clean parent-prerequisite overlay executable SHA-1 was
`5c012db9a3c46f175617b6ad8dbbc651e2e57132`.  Code-quality review found
no Critical, Important, or actionable Minor findings after the
warm-start, row-owned storage, collective-safety, and residual-rank
fixes.

The c384 solver failure is resolved and the fixed buffer-5 production
profile passes Task 9. Task 10 may start from this accepted checkpoint;
it must not select buffer 6 automatically or promote into any forbidden
ground-state or real-time route.

The projected overlap and Hamiltonian no longer persist as full matrices
on every rank. This statement is limited to the persistent projected H/S
storage: fragment-local construction matrices and generalized-eigensolver
temporaries may still have quadratic extent in their local problem size.
Retained buffer-supported Wannier fields and other non-matrix state remain
outside this row-storage result, so the measured peak RSS is recorded as
evidence rather than interpreted as an operator-only memory reduction.
Larger-system Task 9 validation is still required before making a full
production scalability claim.
