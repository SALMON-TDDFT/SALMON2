# Corrected noncollinear MGGA kernel verification

Run `python3 testsuites/unit_nc_mgga23/run.py` from the active source tree.
The runner audits the entire source before compilation, compiles the kernel and
unit program with GNU Fortran runtime checks and invalid/zero/overflow traps,
audits every application object, and runs the test. It uses a temporary build
directory and requires no production calculation or ground-state data collection.
`FC` selects the compiler (the supplied flags are GNU Fortran flags).
`SALMON_FORBIDDEN_CHECKER` may specify the canonical workspace checker.

Tests cover:

- All 36 energy partial derivatives against independent centered differences.
- A proper global spin rotation.
- Every scalar input independently replaced by NaN and positive Infinity.
- Explicit vacuum return.
- Exact hydrogenic BR exchange at three radii and vanishing fully polarized correlation.
- All partial derivatives on both sides of, and exactly at, zero BR curvature.
- Independent explicit evaluation of the 2024-corrected CS correlation formula,
  and agreement with unaffected original Eqs. (14) and (15).
- U(1) phase-gradient transformation of current and kinetic density.
- Local SU(2) spin rotation using explicit spinor orbitals, their first derivatives,
  and Laplacians, including all induced gradient, Laplacian, current and tau changes.

The kernel takes Pauli **trace** components `(0,x,y,z)`: the underlying Hermitian
matrix is `(x0 I + x.sigma)/2`. Its kinetic density is the paper's definition
`T_ab = sum_i grad(psi_ia).grad(psi_ib)*`, **without a factor 1/2**.
Returned derivatives are with respect to these trace components. Consequently,
the corresponding action matrix is `d0 I + d.sigma`, without a further half.
The exchange curvature scaling parameter is fixed at the paper's default 0.8.
The BR root is differentiated implicitly; no finite differences enter production
energy derivatives. The correlation contains both corrected `2*bar_tau` and
`(m.lap(m)-n*lap(n))/(4*n)` contributions, with prefactor `b*r_s**5/2`.

Error codes: 1 nonfinite input; 2 negative density or magnetization exceeding
density; 3 unrepresentable BR root/invariant; 4 nonfinite energy/derivatives.
For density at or below `1e-25` atomic units, the explicitly defined vacuum
convention returns zero energy and derivatives. This cutoff requires spatial and
vacuum convergence testing in production applications.

These tests verify a local energy kernel and its derivatives; they do not establish
solid-state accuracy, full SALMON gKS action correctness, temporal convergence,
frozen-core/NLCC consistency, or Fujitsu compiler compatibility.

References: Tancogne-Dejean, Rubio, Ullrich, Phys. Rev. B **107**, 165111 (2023),
DOI: 10.1103/PhysRevB.107.165111, and Erratum, Phys. Rev. B **110**, 199902 (2024),
DOI: 10.1103/PhysRevB.110.199902. The original arXiv v2 is pre-erratum.
