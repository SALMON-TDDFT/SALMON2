# Ni spinor screening and HSE range fit

Use the existing fcc Ni LSDA+SOC+NLCC restart with Ni_CA19H.vps
(SHA256 3d9e0ae32810b55cc7795b1f044f83da7c2327224388bfb02fb0cf25e9893d4b).
Require SCF convergence and choose the real-space grid from the hard-potential
convergence runs. Preserve the Fermi occupations and k weights.

1. Validate the two-component charge vertex with a small reference Hamiltonian.
   Check spin rotation invariance and reduction to two independent spin blocks.
   For equilibrium metallic occupations verify total electron count across k.
2. Connect the SOC GS restart to the existing spinor Hamiltonian action.
   Test the charge vertex and both spinor components against an explicit
   finite-temperature band sum, including the degenerate occupation derivative.
   The existing broadened post-pump solver is useful at finite q. Its
   eta convergence and the metallic q limit need separate verification.
3. Use the dense solver only as a small-system oracle. At Nr=18 one complex
   spinor matrix is about 2.03 GiB. Multiple copies per rank make this unsuitable
   for the existing four ranks per node. Develop the production resolvent
   using the spinor Hamiltonian action without forming dense grid matrices.
4. Retain the microscopic dielectric matrix and its local-field terms.
   Check q, k, eta and dielectric cutoff convergence; compare finite-q
   extrapolation against the intraband static limit at fixed smearing.
5. Build the same Ni 3d pair densities for both kernels. Fit the exchange
   matrix elements of W to 0.25 times erfc(omega*r)/r, with a common
   reciprocal-space and cell convention. Report the fit objective,
   per-pair residuals, and residuals at 0.11 bohr^-1. Inspect sensitivity
   to the 3d subspace, mesh and screened interaction cutoff.

Production computation stays in Fortran. Do not use the closed-shell GS W
entry point on Ni. The first patch prepares the internal response driver to
receive two spin components. It does not connect the Ni GS reader or provide
a production solver or a fitted omega.

The SALMON adapter now flattens both components in spatial-grid-first order,
uses the spatial volume element, and builds Hamiltonian columns with the
corresponding spinor component selected. A bounded SOC GS diagnostic is available through the existing preparation flag. The existing scalar Si GS
restart produced a bitwise identical W after the adapter change.

Run the spin rotation check with the audited pipeline probe:

```sh
bash testsuites/hsex/build_fortran_screen.sh /scratch/ni-screen-probe
python3 testsuites/hsex/test_spinor_screen.py /scratch/ni-screen-probe/pipeline 4 12
```

The comparison uses the same spin-degenerate Hamiltonian and occupations in
scalar and rotated spinor representations. The test checks the full output W
after local-field assembly and metallic head closure. It does not validate the
finite-temperature degenerate derivative or a physical Ni SOC calculation.

For the dense SOC diagnostic, use `yn_hsex_prepare_w='y'` with a
finite-temperature PZ, noncollinear SOC GS, a full uniform k mesh and k-only
MPI. The q dimensions must be at least six and divide the k dimensions. The
SOC orbitals are exchanged in q batches. After SCF convergence the intended
output is `hsex-spinor-w.bin` (HXN1). Its reference-potential fields are
placeholders and the existing RT W reader rejects the tag. This diagnostic
uses the finite-eta response and metallic extrapolation; an independent
finite-temperature intraband reference and a successful SOC GS W run remain
to be verified. It is not the production route for the existing Ni grids.

The small SOC GS route now uses the static Sternheimer callback instead of
assembling a grid Hamiltonian. The explicit target subspace uses Fermi divided
differences, including the derivative for coincident energies. The projected
empty complement is solved by conjugate gradients and checked against the
Hamiltonian action. This assumes the computed bands cover all appreciably
occupied states and that the complementary operator is positive. The four-state
band-sum test covers a degenerate intraband term and a non-diagonal complement.

The SOC GS screening driver preserves the k decomposition and exchanges
source orbitals in bounded q batches. Scalar GS and RT screening use the same exchange. Measure the total workspace and runtime at the production grid.
The fitted omega is pending.
