# Fast-PR k/PP/legacy GS admission — 2026-11-18

Scope is structural pre-native admission, not certification of physical GS identity
or material convergence. Existing native binaries/kernels are unchanged. Admission
runs within existing MPI error agreement before work copy, CDLL load and native
lifetime consumption. No import-time build, network or MPI initialization added.

Checks cover generated shifted uniform meshes, native odd-Gamma shift handling,
affine retained stars/representative counts and weights; readable recognized PP
references and plain scalar KY numeric records; legacy GS metadata, occupation
records and effective grid/state/k orbital byte length. KY accepts continued
records, repetitions, ignored trailing fields and implicit signed exponents.
Other PP contents remain owned by native. GS producer rank count need not equal
the new runtime group. Extra stored conduction states and real orbitals are supported.

Finite occupation range/electron-sum deviations warn without clipping or changing
inputs. Default startup does not rescan the full orbital stream; explicit private
deep diagnostics are separate. This does not detect arbitrary NaN orbitals before
native startup. Resource rank/memory budgets are advisory. SCF threshold failure
also warns; completion and convergence are distinct recorded quantities.

Current native full-q profile still excludes arbitrary file_kw, unfolding, even
Gamma endpoint meshes and non-z-preserving retained initialization groups. These
are execution limitations, not rejection solely because a case is untested.
Legacy records cannot prove cell/XC/PP/k-order physical identity. Native-endian4
record markers/integers/logicals and global binary64 real/complex orbitals are the
supported legacy layout; no conversion or fixed-binary patching is performed.

Local verification (OMP/BLAS1):

- 583 source tests pass, one optional skip.
- Six asset fault kinds: incompatible k shift, k count, empty PP, truncated info,
  nonfinite occupations and truncated orbitals. On each of two independent MPI2
  groups, only the last rank has the fault. Collective rejection occurs before
  lifetime use; repair permits startup, all32-step current/energy traces match
  the prior reference bit-for-bit. Final reviewed parser rerun passes all six.
- Fresh four-process direct resume from the same saved split7 checkpoint reaches
  step32 with exact continuation trace and final checkpoint bytes, no temporal replay.
- Normal FDTD/LD/SBE/reduced-HSE mixed sample passes the existing full-k reference
  and prior-result comparison.
- Ordinary PBE, non-symmetric fast-PR Taylor/PT-CN and combined seven-rank mixed
  protected-step regressions pass against existing references.
- Reviewed KY grammar regressions pass; the unfolding guard now reads the native
  unfolding namelist rather than analysis.
- Source-distribution-to-wheel/native checksum and isolated installed-package
  suite are recorded with the release evidence. No GitHub Actions launched.

Raw logs/scripts/checkpoints and final packaging evidence are retained outside
this repository under the developer outputs directory:
`salmonpy-hse-assets-admission-20261118`.

General-cell native initialization, learned/PT-CN retained symmetry, richer state
ownership and broader environment/archives remain separate next work.
