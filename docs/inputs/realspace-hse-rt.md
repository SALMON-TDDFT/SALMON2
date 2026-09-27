# Real-space HSE time propagation

DC-SCF/LCFO supplies initial mesh orbitals only. The time-dependent state is the
complex SALMON grid wavefunction. Standard mesh Taylor4 applies the kinetic,
local, nonlocal and HSE operators without any fixed-LCFO projection. Occupied
MLWF rotations only select an exchange-source representation; they do not define
the propagation space. ACE compresses the exchange operator, not the complete
Hamiltonian or the wavefunction space.

For the supplied DC restart inputs, keep `yn_conventional_from_dcdft='y'` and use:

```fortran
&functional
 xc='hse06'
 yn_hse_realspace_rt='y'
 yn_hse_wannier='y'
 hse_rt_wf_radius=0d0
 hse_rt_ace_interval=1
 hse_rt_u_interval=1
 hse_rt_fft_batch=1
 yn_hse_rt_fft_measure='n'
 yn_hse_rt_seed_distributed='y'
/
```

All controls are namelist entries. Radius is always bohr; 0 means full periodic
support. Positive radius masks MLWF sources only, without renormalization; the
propagating orbitals, Hamiltonian action, density and Hartree field are never
masked. Initial geometric norm coverage is saved in `grid_mlwf_radius.dat` and
warns below 99.9%. Unreliable periodic centers are protected and left uncut.
Centers are fixed from initial localization; accuracy needs validation for long
or strongly driven trajectories. Radius is not selected automatically.

ACE interval counts physical steps; impulse step 1 always rebuilds starting ACE.
U interval is independent. Predictor U/anchor changes are rolled back before
corrector propagation. Standard SALMON Hartree/XC predictor-corrector remains.
ACE factors are real-space rows, spatially distributed, and shared across orbital
groups. Application uses a spatial collective of occupied-space overlaps.

Current support: Gamma, fixed fully occupied spin pairs, unpolarized, periodic
orthogonal cells, spatial and orbital MPI, Taylor4. Restart/checkpoint metadata
for the time-dependent gauge/ACE state are not yet supported; wavefunction
output for analysis is allowed. `yn_hse_lcfo_rt='y'` and
`yn_hse_lcfo_direct_wf='y'` now produce an explicit error. Remove these old flags.
Other old `hse_lcfo_*` controls do not configure this new route.

## Exchange and scaling limits

The reference exchange kernel uses the full system periodic screened Coulomb
kernel, with source columns distributed among spatial ranks. Target columns are
communicated one at a time; no rank retains a full-grid copy of every source.
FFTW works on complete periodic grids per pair. Finite source supports reduce
point work but do not shrink this FFT domain. This is a correctness baseline,
not the previous fragment-FFT scaling result.

Large occupied-space matrices (MLWF links, U, ACE metric/transform) and global
column collectives still limit memory and scaling. No new large 3D performance,
Fugaku execution, or converged dielectric spectrum is claimed. Do not use old
fixed-LCFO timing or radius-error tables as results for this implementation.
