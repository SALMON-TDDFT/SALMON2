# Use k-point symmetry in HSE calculations

[日本語](FAST_HSE_SYMMETRY.md)

Group symmetry-equivalent k points to reduce computation. Both the crystal and the driving field must preserve the retained symmetry.

## Supported conditions

- Periodic, fixed-ion, unpolarized HSE calculations on CPUs. This does not imply SOC or molecular support.
- One or two ranks per electronic MPI group, using k-point decomposition.
- Taylor and PT-CN propagation: `hse_taylor4` supports factor off/strict/learned; `hse_ptcn` supports off/learned.
- A calculation library supporting symmetry, with a consistent GS, k-point and symmetry input bundle.
- Keep the same MPI decomposition for save/resume. Resuming does not replay accepted time steps.

## Enable it from Python

Symmetry is disabled by default. Enable it explicitly:

```python
from salmonpy.fast_hse_mpi import NativeFastHSEBackend

backend = NativeFastHSEBackend(build_record, work_root, allow_symmetry=True)
```

`build_record` identifies the compatible calculation library; `work_root` is a new working directory.
See [independent HSE states](SHARED_HSE.en.md) for creating multiple calculations from a GS.
There is no automatic switch to a symmetry-disabled calculation or conversion of a different GS.

## Field direction and symmetry

Supply Cartesian three-component A/c. x, y, z and circular polarization are supported, provided every retained symmetry operation leaves A/c unchanged at every step.
For example, C2z can preserve a z-directed field, and mirror-z can preserve an xy-plane field or circular polarization. Use only operations also compatible with the actual atoms and grid.
An incompatible field is rejected before updating. The waveform is not projected or redirected automatically. Full-k calculations are not restricted to z-directed fields.

## Input and save/resume

Check `sym.dat` and the operations selected by `yn_symmetry`. Each operation in `sym.dat` has three rows of four real values. Input values are not filled in or rewritten.
Prepare the cell, atoms, grid and representative k-point weights as one consistent material description. General cells also require compatible grids and retained symmetry operations. Input units are atomic units.

The `mpi_mixed_native` verification sample accepts `--hse-symmetry`. It also needs an equivalent full-k GS for reference comparison, supplied with `--hse-reference-bundle FULL_K_BUNDLE`. Use `--hse-ranks 1` or `--hse-ranks 2`.
Save the field, LD, SBE and electronic states together at the same completed step. See [checkpointing](CHECKPOINTING.en.md).

## Check your results

Small tests compare full and reduced k points, independent states and resumed calculations. They do not establish physical convergence for arbitrary materials. Check k-point, grid-spacing and time-step convergence for your material and waveform.

[Implementation and validation details for developers](../docs/development/hse-symmetry-notes.en.md)
