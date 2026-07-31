# Obsolete DG Route Removal Results

## Accepted boundary

The branch retains exactly two production surfaces:

1. buffer-supported, symmetry-preserving overlapping-Wannier GS, its
   `SALMON_OW_GS_CHECKPOINT_V3` checkpoint, and generalized-eigenvalue
   exponential coefficient RT with dedicated restart; and
2. normal DC LCFO with EigenExa.

Direct real-space DG GS, DG-Fragment RT, Nodal RT, WPW GS/RT, mixed-Z,
full-H seed, adaptive DG promotion, and experimental ExpDiag routes and
their public inputs, tests, samples, and executable documentation were
removed. The only file below `src/rt/dg/` is
`rt_dg_overlapping_wannier.f90`.

## Source and focused verification

The final source contracts passed:

```text
python3 tests/dg/check_obsolete_dg_routes_removed.py
  obsolete DG route removal contract: PASS
python3 tests/dg/check_dg_overlapping_wannier_route.py
  overlapping-Wannier route contract: PASS
git diff --check
  exit 0
```

The retained metadata, complete-s+p projection, metric, construction, weak
operators, physical matrices, solver/density, SCF, V3 checkpoint, and
coefficient-RT runners each passed on 1, 2, 4, and 8 MPI ranks. Normal
LCFO/EigenExa source contracts and the 34-test SAWF alignment suite also
passed.

## Clean-first parent-prerequisite overlay

A fresh source was created from `git archive HEAD`. The current Task 8 diff
was then applied. The enumerated parent-prerequisite overlay was empty:
all retained prerequisites are committed in the branch ancestry, while the
parent worktree's remaining uncommitted files are forbidden WPW sources and
were deliberately not imported. This prevents a dirty-parent overlay from
silently restoring a removed route.

Configuration enabled MPI, ScaLAPACK, and EigenExa with GNU Fortran 15.2:

```text
cmake -S <overlay-source> -B <overlay-build> \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/mpifort \
  -DCMAKE_Fortran_FLAGS=-fallow-invalid-boz \
  -DUSE_MPI=ON -DUSE_SCALAPACK=ON -DUSE_EIGENEXA=ON
cmake --build <overlay-build> --clean-first -j1
  [100%] Built target salmon
```

The final executable SHA-256 was
`bff463844fc351195268a88281463f523da9f815d006def96647d3d6a71030e8`.

## Production gates

The unchanged normal Si64 case converged at DC SCF iteration 86 with total
energy `-0.99454347E+04`, density difference `0.83895929E-09`, and 256
electrons. Its log contains, in order, `start DC-LCFO-Flux`,
`eigenexa diag, #dim=552`, `eigen_sx: done`, and `end DC-LCFO-Flux`.
Runtime was 613.29 s and maximum resident memory was 920,240,128 bytes.

The fresh 8-rank buffer-6/c192 Si64 overlapping-Wannier run passed the
public minimum-row gate. Selected evidence was:

- metric minimum `2.6625467811430542E-01`, condition `7.1757331256424566`;
- symmetry center defect `2.2204460492503131E-15`;
- density residual `1.9772241369239652E-08` and unmixed residual
  `8.5339113507206519E-10`;
- coefficient residual `1.8098975373646166E-12` and S-orthogonality defect
  `1.9317880628478099E-14`;
- trace/integrated charge `256.00000000000017` / `256.00000000000068`;
- runtime 1415.29 s and peak memory 1857.23 MiB.

The V3 manifest SHA-256 was
`1a38778c27350ae07173ba891e94b8a9fa0b67f5ce223c83ca821b77e8f0269c`.
Restart logged `[OW-GS] reused accepted route checkpoint`, did not rerun OW
SCF, and left this hash unchanged.

Production coefficient RT used 8 MPI ranks and the accepted V3 checkpoint:

| Gate | Runtime | Peak bytes | Restart SHA-256 |
|---|---:|---:|---|
| 2-step field-off | 0.36 s | 131,104,768 | `29e7943cbbe16f019cab5edc2ad88f28c5b3a1b30d202fdf277f45b24bbe81c9` |
| 2-step impulse | 0.46 s | 137,412,608 | `88a8ef5d9c33a6e1dcefd667a9034deb7bcaa81b9d0dc562268450c7e8cf838e` |
| 4-step one-shot | 0.38 s | 131,629,056 | `2402169a0db39fe9b11340d2db26188d429a12c6afe8b07b4a62fef7b393b532` |

The impulse and field-off artifacts differ, proving that the impulse acts
on the coefficients. The 2+2 restart-split artifact is byte-identical to
the 4-step one-shot artifact and has the same SHA-256.

No production log contained direct-DG, WPW, DG-Fragment, Nodal, normal
checkpoint publication, conventional RT initialization, or conventional
time-evolution markers.

## Review disposition

Specification review found and resolved one Critical production issue:
the coefficient driver originally differenced two already-jumped impulse
vector-potential samples, producing a zero impulse. A RED source contract
and production field-off/impulse artifact comparison were added; the fresh
first interval now uses the zero pre-impulse vector potential. The focused
1/2/4/8 coefficient fixture, clean-first build, and production gates passed
after the fix.

Code-quality review found and resolved two Important test-harness issues:
ignored `.dat` fixtures were force-tracked for clean archives, and retained
SAWF tests were detached from the deleted experimental sample tree. It also
made the SAWF compiler-build prerequisite explicit. After these resolutions,
the final branch review found no remaining Critical or Important findings.
