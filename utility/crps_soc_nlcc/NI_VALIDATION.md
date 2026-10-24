# Ni SOC+NLCC KY candidate (2026-10-24)

## V2 same-shell NLCC extension

The selected input now enables V2 auxiliary-Pauli core primitives. The new file
is in `calc2026/ni-soc-nlcc-20261024/v2-aux-pauli-validated/Ni/` under canonical
OneDrive. Its ionic potentials, waves and all legacy NLCC density columns are
byte-identical to the selected V1 numeric rows. See
[CORE_AUX_PAULI_V1.md](CORE_AUX_PAULI_V1.md) for model and validation details.
All V1 files remain preserved.

## Current larger-radius candidate

Following the requested reduction of internal potential oscillations, the input
now selects s=3.2, p=2.8, d=3.2 bohr, s-local, NLCC radius0.45 Å.
Actual mesh-rounded radii are 1.69/1.48/1.69 Å. The generated file is:

`calc2026/ni-soc-nlcc-20261024/larger-rps/s3.2-p2.8-d3.2-local0/Ni_soc_nlcc_rps.dat`

under the canonical OneDrive root described below. SHA256:
`3754b1a364d81828629fa3ef4eaf8514395b80d66f7ec8f936dd471749683f2e`.

Fifteen larger-radius candidates were generated without relaxing the unique-c2
root requirement. Both initial and selected potentials are monotonic inside
each Rps; the selected candidate has no added stationary points. The RMS radial
curvature drops by 43.2% (s), 45.8% (both p channels), and 71.2–71.3% (d).
The shape metric excludes the forced origin duplicate and uses second radial
differences; no monotonicity or smoothness is inferred from radii alone.
The initial potential is therefore described as steeper, not falsely labeled
as oscillatory.

The larger p=3.2-bohr candidates introduce a small ~0.2–0.3-eV wiggle with two
stationary points. More importantly, combined p=3.0 and p=3.2 candidates yield
additional shallow screened KB p states (e.g. about -0.67 eV for p3.0, compared
with a positive semilocal non-reference state); those candidates are rejected
for selection despite no deep-ghost flag. The independent radial tool now
records additional-bound-state flags separately from the deep-state test.

The selected 3.2/2.8/3.2 candidate passes 60 sampled radial checks with zero
deep flags and matching screened bound-state counts on every grid/box.
At h=0.005 Å, maximum reference-eigenvalue error is 3.47 meV (semilocal),
3.45 meV (KB); KB/semilocal reference difference is at most 2.12e-5 eV.
The generator norm-equation residual is at most 3.07e-15.

Evidence: `larger-rps/larger-radius-validation.json`, each candidate
`radial-screen.json`, selected `radial-converged-screen.json`, and
`larger-rps/ni-selected-potential-curvature-comparison.png`. The plot was
visually inspected. `ni-potential-curvature-comparison.png` is the prior p3.0
plot and remains as rejected-candidate evidence.

No full all-electron scattering or solid transferability certification is
claimed. Initial smaller candidates and all prior records remain unchanged.
The records below describe the earlier smaller-radius validation stage.

## Earlier 2.2/2.0/2.0-bohr stage

`inputs/Ni.nml` is a **candidate**, not a transferable production pseudopotential certification.
It uses neutral spherical Ni, Z=28, [Ar]3d8 4s2; the d8 occupations are split
3.2/4.8 across j=3/2,5/2. The 16 explicit electrons are 3p6 3d8 4s2;
3s2 remains in the core. The selected local channel is s (`lloc=0`); initial d-local and p-local candidates fail the radial ghost screen. Reference radii are
s=2.2, p=2.0, d=2.0 bohr (rounded on the 0.01 Å writer mesh), with
NLCC smoothing radius 0.45 Å. These are exploratory physical inputs.
The generator uses its existing relativistic LDA/PZ atomic reference.

All generated tables, atomic-generation logs, unsuccessful-input evidence, and
records are in the canonical local OneDrive directory:

`/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/ni-soc-nlcc-20261024`

The selected candidate file is `candidate-2.2-2.0-2.0-local0/Ni_soc_nlcc_rps.dat`, SHA256
`d6596e68fb0cea6fc7a3dbd27c404122aaf2c1aef794de01d49f22575da3413d`.
The original d-local file remains preserved as failed evidence.
Cloud synchronization was not confirmed. No SALMON GS data was collected.

## Checks completed

- Source prohibition audit passed before the utility build. Both linked application
  objects passed the object-tree audit. No PZHEEV was executed.
- GNU Fortran 15.2.0, legacy input syntax, O0, runtime bounds checks/backtrace.
  Compiler and source/executable SHA256 are recorded in `manifest.json`.
- Seven candidates: central radii and independent ±0.2 bohr changes of s, p, d.
  All five channels in each candidate have exactly one bracketed norm root
  within the unchanged [-10.5,5] interval. This is uniqueness in the sampled
  search scheme, not a mathematical global uniqueness proof.
- `validate_ni.py` verifies header/units/grid, every numeric entry finite,
  nonnegative core density, occupations, j ordering, nonzero j potential
  differences, zero core derivative at the origin, and all interior NLCC
  first/second derivatives against the writer's five-point formula.
- Central maximum generator norm-equation residual: 2.03e-15. This checks
  the generator's norm constraint, not an independent all-electron reconstruction.
- Atomic reference eigenvalue splittings carried in the header: p 2.21274630 eV,
  d 0.24045601 eV. These are reference-atom splits, not solid-state bands or
  independent pseudo-atom eigenvalue comparisons.
- Narrow scalar IEEE fixes preserve all elements checked in generator inputs
  and radial writer. Regenerated W/Se tables are byte-identical to fixed
  distributed references. Eight NaN/Infinity cases (occupation/radius/c2lo/c2hi)
  are rejected before any output file is created.
- The generator reports IEEE underflow after ordinary successful runs; every
  output entry is finite. It is recorded rather than treated as an accuracy test.

## Independent radial physics screening

`screen_ni_radial.py` reconstructs the generator's spherical valence Hartree and
relativistic LDA/PZ screening from the tabulated waves/occupations/NLCC. It builds
second-order finite-difference radial Hamiltonians in Angstrom/eV using the
original generator constants. Each j-channel is checked with both the semilocal
potential and the reader-equivalent separable operator:

`beta_j=(V_j-Vlocal) u_j; N_j=integral u_j beta_j dr; V_KB=|beta_j><beta_j|/N_j`.

The norm uses the native table's trapezoidal rule, matching `read_ps_ky_soc`.
The interpolation between table points is four-point cubic; this standalone
radial interpolation is not a claim to reproduce SALMON's 3D grid interpolation.
NumPy dense Hermitian diagonalization was used because SciPy was unavailable.

The original d-local candidate has a screened s-channel ghost near -95.16 eV,
while its semilocal ground state is near -5.857 eV. It is rejected. The p-local
alternative also gives deep-extra-state flags and is rejected. Neither is an
approved input for a laser calculation.

The selected s-local candidate has zero deep-extra-state flags in 60 checks:
all five j channels, ionic/screened operators, boxes 8/12 Å, steps
0.04/0.02/0.01/0.005 Å (the complete grid list is recorded in the JSON).
The diagnostic flag is KB ground energy more than 0.5 eV below semilocal ground;
this explicit threshold detects deep ghosts, not arbitrary shallow ghosts.
Screened bound-state counts agree. At h=0.005 Å, the maximum reference-energy
error is 5.72 meV, and the maximum KB/semilocal reference-energy difference is
5.28e-5 eV. Reference-wave large-component norms are 0.99965–0.99984;
they should not be silently renormalized to one as a test of Dirac total norms.
The largest h=0.01 Å reference-energy change between boxes8/12 Å is approximately
2.5 micro-eV. The h=0.01→0.005 Å reference errors shrink approximately as h².

The script additionally compares screened finite-radius logarithmic derivatives
at R=8 Å, h=0.01 Å and E=0,2,5,10,20,50 eV. These differ between KB and
semilocal operators (maximum p-channel difference 2.51 Å^-1, d 0.234 Å^-1).
The ionic p-channel KB operator has one fewer negative finite-box level than
its semilocal counterpart. These differences are retained as transferability
warnings rather than presented as successful scattering equivalence.
All-electron scattering was not available for independent comparison.

Results are in `candidate-2.2-2.0-2.0-local0/radial-converged-screen.json` and
`physics-validation-manifest.json`. The older `manifest.json` describes the
initial structural-only d-local stage; it has not been rewritten.

## Validation still required

Atomic scattering logarithmic derivatives against all-electron data across
energies and alternate occupations; independent eigenvalue/norm comparison;
exhaustive separable-projector ghost-state tests beyond this sampled spherical
screen and actual SALMON 3D interpolation; core-radius and local-channel sensitivity; real-space mesh and solid
SOC convergence. The s-local candidate passes the documented sampled deep-ghost screen; no
candidate has full transferability certification. Do not infer
absence of ghosts from smoothness, successful generation, or nonzero SOC splits.
The ±0.2 bohr scan establishes generator robustness only, not a quality ranking.
Fujitsu compilation was not available. The file has not been accepted for a
production laser-demagnetization comparison.

To repeat the structural checks:

```sh
python3 utility/crps_soc_nlcc/validate_ni.py /absolute/path/Ni_soc_nlcc_rps.dat /absolute/path/stdout.txt
```

Keep every candidate's file and log together. Output names are intentionally
not installed over any existing pseudopotential.
