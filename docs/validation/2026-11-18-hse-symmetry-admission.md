# Fast-PR HSE symmetry initialization admission (2026-11-18)

Base: de59e45c (internal MPI2), committed/pushed before this work. Scoped next
item from the initialization/failure checklist. Python-only admission change;
no Fortran source, fixed native binary, exchange formula or field is modified.

## Contract

Before native loading/lifetime/work-copy creation, reject malformed sym.dat
records, nonfinite/invalid Fortran literals, invalid repetition counts and invalid
canonical yn_symmetry flags. Use the same disabled-axis filter and bounded
fast-PR signed-permutation/identity/duplicate/closure rules. Supported records
have three lines per operation and four finite reals per line; D exponents,
commas, unsigned decimal repeat counts and inline ! comments are accepted.
Inputs are never normalized, repaired or rewritten by admission.

All ranks in the solver subgroup agree on errors through the existing collective
entry. No native lifetime is acquired on rejection. Users repair the source and
construct a new provider; previously prepared-provider hash guards remain intact.
Cell/atom/species/grid/k-star/restart binary validation and other native STOP
paths are NOT covered. This does not establish general initialization recovery.

## Local verification

Apple Silicon/macOS, GNU Fortran15, OpenMPI5, Libxc5.2.3; reused the unchanged
explicitly built internal-MPI fast-PR library. No GitHub Actions.

- New parser regressions fail before implementation; final source suite560 tests,
  one skip. Covers empty/truncated/wrong-width/nonfinite records, invalid numeric
  and repeat syntax, canonical flags, valid D/repeat forms, disabled-axis filter,
  identity, duplicate modulo translation and non-closed groups.
- Backend test asserts no CDLL call or new work directory on malformed input.
- MPI4: two independent groups of two ranks, omega=.11/.12. Only the last rank
  in each group receives truncated, nonfinite or nonclosed data. All ranks reject
  with InputError and native lifetime flags remain false. Repaired inputs then
  initialize in the same processes. Continuous current/energy traces match the
  previous normal32-step traces bit-for-bit for all three cases.
- Own save7→32 fresh-process direct resume matches current/energy traces and final
  checkpoint bytes exactly, without elapsed-step replay. Cross-initialization
  checkpoint workspace bytes are not claimed invariant; equality here is against
  the matching continuous run from the same saved state.
- Normal FDTD + two local LD materials + SBE + internal-HSE2 mixed sample passes
  the existing same-physics comparison. Final installed-wheel mixed arrays match
  the source run. No normal input is rewritten.
- Source-distribution→wheel/native-hash checks pass. The final wheel runs560 tests
  (one skip), the normal mixed sample, and the packaged MPI admission/retry tool.
- Independent read-only review found and prompted a repeat-count lexical guard;
  unsigned decimal counts are now checked before Python int conversion.

Run the explicit save tool with --group-size2 (four total MPI ranks) and add
--admission-fault truncated|nonfinite|nonclosed. The tool recreates a provider for
its deliberately altered input, then restores the original file before normal
initialization; it does not bypass immutable input manifests.

Local evidence:
`/Users/otobetoshihito/.codex/outputs/salmonpy-hse-symmetry-admission-20261118/`.
The bundle retains failed development logs, final logs, checkpoint/trace evidence,
wheel/source package, changed sources and hashes. This is numerical/admission
verification, not material convergence or cross-platform validation.
