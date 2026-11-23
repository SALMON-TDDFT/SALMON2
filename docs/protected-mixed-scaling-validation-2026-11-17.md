# Protected mixed-step resource campaign — 2026-11-17

This extends the [native protected-step validation](native-protected-step-validation-2026-11-17.md)
to a larger Maxwell grid. Solver sources and scientific update order are unchanged.
The previously validated wheel and relocated byte-identical production libraries
are used; there is no build on import, download or remote Actions run.

## Scope and acceptance

The user host has 64GB RAM and 18 cores. Runs are sequential with OMP/BLAS1,
48GiB advisory memory budget and at most nine ranks. Memory budgets remain
warnings; the 4GiB/four-core beginner profile is unchanged.

Every case uses a 96³ Maxwell grid (884,736 cells), two local LD materials, one
SBE point, one Libxc-PBE TDDFT point with two orbital ranks and a separate fast-PR
HSE learned PT-CN point with two k ranks. Field2 and field4 therefore use seven
and nine ranks in total. There are 64 accepted steps and two repetitions for each
layout and protected/unprotected mode (eight ordinary runs).

Only the Maxwell grid is enlarged. Electronic grids, k meshes and material
point counts are unchanged. The larger box adds vacuum; comparisons use the same
96³ inputs, never the earlier 40³ physical domain. This is a connection/resource
regression, not a grid-convergence or real-material accuracy result. Practical
TDDFT multiscale calculations still start with the 1D route.

Acceptance: all response arrays and the multiset of ordinary `.data` file names
and SHA256 hashes must equal the first uninterrupted case, across both layouts,
repetitions and protection modes. Time and memory measurements cannot decide
scientific pass/fail; measurement failure is reported separately.

## Measurements

All eight ordinary runs passed exact array and ordinary output-byte comparisons.

| Field ranks / total | Plain mean (range) s | Protected mean (range) s | Time ratio | Plain / protected maximum summed peaks GiB |
|---|---:|---:|---:|---:|
| 2 / 7 | 12.00 (11.85–12.15) | 14.38 (14.27–14.50) | 1.198× | 2.012 / 2.692 |
| 4 / 9 | 7.08 (6.97–7.19) | 7.66 (7.63–7.69) | 1.082× | 2.949 / 4.417 |

Loop time includes communication waits and protected transaction work; setup,
final output and serial response replay are excluded. End-to-end and per-rank
communication measurements are retained separately. No affinity/frequency or
background-load control is imposed. Two repetitions show local variability;
they do not establish production scaling or an optimal rank allocation.

Memory is the sum of independently observed rank peak RSS before the later serial
reference computation, not simultaneous aggregate RSS. Peak differences between
separate runs are empirical observations, not isolated allocation accounting or
a proven heap upper bound. Existing memory-advice unknowns remain unknown; these
results do not justify replacing them with a universal multiplier.

## Failure and retry

The seven-rank 96³ warm-step13 campaign passed both late electronic rejection
and real field overflow. Original numerical payload hashes, Python state and
ordinary output bytes were preserved. Explicit valid retry completed all64 steps;
arrays and every ordinary stream match uninterrupted execution exactly. No
elapsed propagation was replayed. Summed rank peak RSS in this diagnostic path
was 3.943GiB; hashing/fault instrumentation is included and this is not a
normal-run memory estimate.

The dedicated large fault harness derives from
`tools/test_native_protected_step_mpi.py` and tests warm step13, electronic-late
and real native field-overflow failure. To avoid retaining many large duplicate
checkpoint files, it compares SHA256 of the complete numerical payloads and
removes only its own diagnostic checkpoint files after hashing. Python state
references/values/clocks and ordinary stream bytes retain the original checks.
The exact harness, recipes and build receipts remain in the local evidence.
This does not extend recovery to I/O/process/MPI failures or uncontrolled STOP.

## Reproduction and evidence

Use the existing installed mixed sample and the documented separate HSE bundle:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
mpiexec -n 7 python -m salmonpy.examples.mpi_mixed_native \
  STANDARD_BUILD_JSON SERIAL_RUNTIME PBE_BUNDLE NEW_OUTPUT \
  --layers LAYOUT_JSON --field-shape 96 96 96 --field-grid 1 2 1 \
  --tddft-ranks 2 --hse-ranks 2 --hse-bundle HSE_BUNDLE \
  --fast-hse-build SERIAL_HSE_BUILD_JSON --fast-hse-mpi-build HOST_HSE_BUILD_JSON \
  --steps 64 --memory-budget-gib 48 --rank-budget 9 --protect-steps
```

Remove `--protect-steps` for the paired baseline. For field4, use nine ranks and
`--field-grid 2 2 1`. Give every run a fresh output directory; add `--compare`
with the first run's directory. Runtime inputs/binaries are the same as in the
preceding validation, not replacements downloaded during execution.

Persistent evidence: `/Users/otobetoshihito/.codex/outputs/salmonpy-protected-scaling-20261117/`.
`campaign.json` contains commands, all measurements and paired summaries;
per-run arrays, streams/hash lists, identities and memory reports remain alongside
logs. Harness/recipes and fault comparison are retained. This local path records
provenance; user instructions do not require it.
