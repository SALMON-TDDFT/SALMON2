# Native HSE implementation

This merge candidate extracts the native HSE path from the experimental TDCDFT
branch. It does not contain TDCDFT functionals, ELF feedback, their input variables,
or the accumulated material-specific calculation outputs.

## Functional and input

HSE is enabled by default for CPU builds. CMake finds or builds FFTW, Libxc
(C ABI), BLAS and LAPACK automatically; see [build details](hse-build.md).
The option defaults to OFF. See [input settings](inputs/hse.md) and the
[Si example](../samples/hse_native/README.md).

`xc='hse06'` combines 25% short-range Fock exchange with the matching Libxc
semilocal remainder. `hse_omega` defaults to 0.11 bohr^-1; A_eV_fs input accepts
inverse Angstrom and converts before any exchange evaluation. Both Libxc screening
parameters and the native exchange kernel receive the same internal value.
Custom omega defines a modified HSE functional and requires GS reconvergence.

## Exchange and ACE

The native path uses Bloch orbitals, not MLWF support truncation. The sampled
short-range kernel uses the distance-dependent erfc(omega*r)/r interaction.
K-point sources, targets, phase arrays and ACE factors remain on their MPI owners.
Density-matrix row tiles are exchanged for the FFT work; occupied orbital arrays
are not replicated over all k ranks.

At each exchange refresh, full exchange computes W=K[U]U on occupied orbitals U.
ACE factors are constructed from the occupied metric -U†W with grid weights.
Subsequent trial-vector actions use -B(B†X), two BLAS products. ACE reproduces
exchange on the construction space; general trial directions are approximate.
It neither freezes exchange permanently nor replaces it with a local potential.

GS potential updates refresh exchange from the current occupied orbitals. RT
uses the step-start ACE for prediction and the average of initial/predicted
operators for correction, then refreshes at the accepted state. This is an
average of operators, not of occupied wavefunctions.

## Propagation and fields

Omitting the propagation namelist selects Taylor4 + ACE with predictor/corrector.
The polynomial order is four; the self-consistent predictor/corrector has generally
second-order global time accuracy. Tests must examine time-step convergence,
orbital norms, overlaps and currents. There is no forced orbital rescaling.

Transverse impulse response is supported. Full-k, linearly polarized Acos2
pulses may include a weak impulse probe. Propagation uses midpoint A; current
uses endpoint A with refreshed nonlocal pseudopotential phases. Pulse RT restart
is rejected, including attempts to reopen a pulse checkpoint as an impulse.

Developer-only explicit propagators remain for validation: hse_taylor4_full
(full Fock action on each trial vector) and hse_ptcn (parallel-transport CN).
They are not required in ordinary user inputs. The latter retains endpoint
residual checks; these must not be attributed to Taylor propagation.

## Symmetry and restart

Supported symmetry-reduced impulse calculations validate group closure, atomic
positions/species, grid transforms and the field-preserving z subgroup. Source
states are reconstructed with their Bloch phases before exchange. Density rows
are batched across destination ranks to avoid repeated reconstruction.

RT metadata checks the method, time step, kick, omega, geometry, k mesh,
occupations and pseudopotential bytes, plus retained symmetry. Ordinary GS
wavefunction files do not certify the input functional: the user must supply a
GS converged with the same omega. Non-HSE writes remove stale HSE metadata.
The generic RT buffer alternation is relative to the loaded step so odd-step
restarts use the initialized wavefunction buffer.

## Performance and scope

BLAS calls remain outside application OpenMP regions. The linked library controls
its own threading; runtime settings are platform dependent. FFTW/MPI calls are
outside application OpenMP regions. Packing loops distribute work over k points.

A bounded initialization trial selects the FFT data layout; it does not impose
material-specific input choices. Developer overrides are SALMON_HSE_FFT_LAYOUT
(auto/strided/contiguous) and SALMON_HSE_BLOCK_ROWS. Detailed timing is opt-in via
SALMON_HSE_PROFILE=1; the developer CN solver also uses timing totals internally.

Current scope: fixed ions, periodic, unpolarized, occupied-only fixed occupations,
cubic real-space/uniform k meshes and k-only MPI decomposition. No certified
NLCC, spin-orbit, DFT+U, ionic dynamics or accelerator path. Serial CPU and the
HSE-disabled build are tested separately. Fugaku performance is not established
by Mac measurements.

See [merge verification](merge-preparation.md) for the actual tested configurations.

## Nonorthogonal full-k exchange geometry (2026-10-12)

The full-k CPU exchange initializer accepts the real and reciprocal primitive
vectors already computed by SALMON. It checks `transpose(A)*B=2*pi*I`.
Grid positions are `A*(i/n)`, k-mesh indices are
`transpose(A)*(k-k0)*mesh/(2*pi)`, and reciprocal vectors of the
Born--von Karman cell are `B*(m/mesh)`. Coordinate differences used for indexing
remain integer lattice coordinates. The unscreened spherical-cutoff limit is
half the smallest distance between opposite supercell faces.

The real inverse-transform kernel pairs conjugate discrete frequencies.
On skew even grids this averages the positive multipliers at paired Nyquist
indices. This choice preserves the real, Hermitian discrete operator and the
existing real-kernel backend interface. Orthogonal-grid results are unchanged
to roundoff. No real/reciprocal lattice inverse is recomputed in exchange.

Local validation uses independent NumPy Fourier kernels, Bloch phases and a
Gamma direct pair-density convolution. Errors were 2.8e-16, zero at the tested
points and 8.4e-17, respectively. A Cartesian rotation, the legacy diagonal-cell
API, and rejection of inconsistent dual vectors were checked separately.
Eight-rank Si tests use two atoms, eight states, an unreduced 2^3 k mesh and
12^3 spatial points. GS converges; ordinary and learned K=8 ACE RT complete
192 steps at dt=0.02 a.u., including 98 predicted steps. Matched direct exchange also completes.
Relative current errors are 4.055e-07 (ACE) and 2.689e-05 (learned ACE).
The initial coarse
8^3 test at dt=0.16 became unstable and is retained as a failed diagnostic.
No grid, k-mesh or long-time accuracy convergence is claimed by these checks.

This change enables the integer-occupation full-k route. Fractional-occupation,
Wannier, spatial-neighborhood/SR localization and symmetry-reduced routes keep
their existing geometry restrictions. DC remains orthogonal only.
