# Custom FDTD materials and interfaces

[日本語](FDTD_MATERIALS.md) · Topic: `fdtd_materials` · Paired revision: `2026-11-13.4`

Run `python -m salmonpy.examples.fdtd_layers MANAGED_RUNTIME NEW_DIRECTORY`
with a verified runtime v2. This one-process, one-thread example contains two
synthetic finite layers. It is not a converged infinite-film optical calculation.

Define each isotropic dispersive material with `LorentzDrude(plasma_frequency=...,
strength=[...], resonance=[...], damping=[...])`. All quantities use atomic units;
frequency and damping are angular frequencies. Each material may have its own
number of poles. A zero resonance describes a Drude pole. Use validated material
coefficients; SALMONpy does not fit experimental data or guarantee convergence.
For a constant isotropic medium use `Dielectric(epsilon=..., mu=..., sigma=...)`.
See [LD conventions](README.en.md) before translating coefficients from another
code. The example coefficients are deliberately synthetic.

Pass the materials as a tuple to `FDTD`. Shape ID 1 refers to the first material,
ID 2 to the second; ID 0 is vacuum. Boxes and ellipsoids can be combined. If IDs
are renumbered, change both the tuple and shapes consistently. At overlapping
shape cells the last shape wins: this geometrical order is intentional. Material
IDs do not encode physical priority. Inspect the generated material-ID cube.
Keep material clear of the vacuum PML, and test grid/time convergence. Current
limits are 200 materials, 100 poles per LD material, and 1000 shapes.

The interface rule is a deterministic staircase assignment, without averaging
unlike dispersive models:

- An electric-field edge between two nonvacuum cells takes the lower-coordinate
  cell's complete material model. If either cell is vacuum, retain the legacy
  vacuum-edge coefficients.
- A magnetic-field site takes the first nonvacuum cell in fixed spatial order.
  For Hx this is `(i,j,k)`, `(i,j+1,k)`, `(i,j,k+1)`, `(i,j+1,k+1)`; Hy and Hz
  use the corresponding transverse planes. A fully vacuum stencil stays vacuum.
- Native MPI initialization exchanges faces, edges and corners of the material
  grid. The existing E/H and LD time-update order remains unchanged.

This rule preserves the tested single-material/vacuum calculations. Results at
interfaces between unlike materials intentionally change from the old numeric-ID
priority. It does not assert subpixel accuracy, rotational invariance or improved
convergence order. At observation points the legacy printed material parameters
refer to adjacent **cells**, not effective interpolated E-edge coefficients.

Runtime v1 is rejected for multi-material runs. Obtain/build v2 explicitly; imports
never compile. `salmonpy-run.json` records `interface_assignment`. A new `run`
creates independent solver state. Recreate the calculation when changing material
parameters or timestep.

This public batch interface supports LD and constant dielectric materials. It does
not yet accept arbitrary current callbacks, anisotropic tensors, tabulated optical
response, periodic boundaries, oblique incidence or TDDFT/SBE materials in E/H
FDTD. The Python batch launcher can also start the existing native MPI FDTD with
`FDTDMPI`; see below. TDDFT/SBE A-to-J interfaces and the native A-form multiscale path remain
separate capabilities.


## MPI FDTD with many local LD states

FDTD was already MPI-parallel in SALMON. This API selects its existing spatial
partition; it does not add a new Maxwell parallel algorithm. Each FDTD rank owns
its local field grid and LD P/J/Q histories, selecting coefficients by material
ID. Several LD materials and many LD cells can reside on one rank. No dedicated
MPI group is allocated per LD cell. Native implicit E/LD/H update order is kept;
the separate endpoint-A LD response-service adapter is not used here.

```python
from salmonpy.fdtd import FDTDMPI, SalmonExecutable
solver = SalmonExecutable.from_mpi_build("MPI_BUILD_DIRECTORY/build.json")
result = solver.run(calculation, "NEW_OUTPUT", mpi=FDTDMPI((2, 2, 1)), threads=1)
```

Run this from a fresh serial Python controller, not under mpiexec or initialized
mpi4py. FDTDMPI's three integers define the spatial process grid; their product is
the MPI rank count. `launcher=("mpiexec",)` can be replaced with an argument tuple
for the host launcher; there is no shell interpolation. Each rank needs at least
two cells per direction. MPI launching currently targets macOS/Linux.

Build explicitly using `python -m salmonpy.runtime build` with USE_MPI=ON and the
host MPI C/Fortran compilers. See RUNTIME.en.md; choose the build job count from
your memory/CPU budget. The successful build record binds MPI capability to the
binary hash. An older record without that capability requires a fresh build;
editing CMakeCache does not confer MPI capability. Portable serial runtime bundles
are not MPI builds. For externally audited builds outside this builder, the
low-level SalmonExecutable(..., mpi_enabled=True) is the caller's explicit
capability assertion, not automatic binary detection.

`salmonpy-run.json` records the launcher arguments, rank grid/count, threads and
completion/failure/cancellation. Timeout and interruption terminate the local
launcher process group; remote job cleanup still depends on the MPI/scheduler
launcher. Do not share one output directory between jobs. Threads is the OpenMP
setting; account for any separately configured BLAS threads in CPU planning.

```sh
python -m salmonpy.examples.fdtd_mpi_layers MPI_BUILD_JSON NEW_DIRECTORY --grid 2 2 1
```

This example uses three synthetic LD materials with 1, 2 and 3 poles. It prints
material-cell counts and leaves memory explicitly unmeasured. Production sizes
need resource planning and numerical convergence tests.

For electronic feedback into E/H, use the separate [native coupling entry](COUPLED.en.md). Do not put TDDFT objects in the batch API material list.
