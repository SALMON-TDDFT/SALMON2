# Formal checks of SALMON formulas (Lean 4)

Each file transcribes formulas exactly as SALMON v2.3.0 implements them
(prefactors, signs, unit conversions, the arrays that enter each term), with
file:line references in the docstrings. Next to each one is the physically
correct formula, followed by one of:

- a proof that the two are equal (`*_ok`), or
- a proof that they differ (`*_bug`, `*_counterexample`, `*_ne`). This is either a
  general relation such as `code = 8 * correct`, or a counterexample with explicit numbers.

Grids, stencil weights, occupations and similar inputs are abstracted to
arbitrary types or reals. Each file says what was abstracted. No proof uses
`sorry`, `axiom` or `native_decide`.

| Directory | Area |
|---|---|
| `XC`   | built-in functionals (PZ, PZM, PW92, r2SCAN, TB-mBJ), libxc interface, E_xc in the total energy |
| `GS`   | FD stencil, Ewald, Hartree, local/nonlocal pseudopotential, total energy, forces, occupations, DOS/PDOS |
| `RT`   | propagators, laser fields, gauges, macroscopic current, response functions, output units |
| `EM`   | classical FDTD, Weyl-gauge multiscale, TTM coupling |
| `SBE`  | semiconductor Bloch equations (velocity and length gauge) |
| `MISC` | spin-orbit / noncollinear, DFT+U, symmetrization, TTM, MD, unit conversion |

Build (downloads the Mathlib cache):

```
cd formal
lake exe cache get
lake build
```
