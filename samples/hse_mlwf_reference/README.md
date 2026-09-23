# Independent Si HSE06 reference

This prototype evaluates self-consistent HSE06 on the discrete Hamiltonian exported by SALMON. It does not enable native SALMON `xc='hse'`. The initial supported case is fixed-ion, unpolarized Si8, 16 doubly occupied orbitals, a shifted uniform 4³ k mesh, an isotropic 12³ primitive real-space grid, and no nonzero NLCC. Numerical agreement with the exported Hamiltonian is checked before SCF.

HSE06 uses 25% short-range Fock exchange, omega=0.11 bohr^-1, and the matching Libxc HSE06 semilocal remainder. The reciprocal screened kernel includes its finite q=0 limit. Hartree and ionic interactions remain unchanged. Exchange uses all translation images of occupied MLWFs, with the k-mesh twist retained. Reciprocal orbital pairs reuse conjugate-translated potentials. Pair selection uses normalized density overlap; zero threshold retains every pair. Spatial support boxes are experimental operator profiling only and are not used in SCF or reference dynamics.

The final SCF residual and energy use all pairs. RT uses RK4 with the current density and exchange recomputed at every stage. It propagates original Bloch orbitals; the unitary MLWF gauge is an exchange representation and is refreshed every 10 steps using the preceding gauge. Stage orbitals are not artificially orthonormalized. Current includes the kinetic and nonlocal pseudopotential contributions. The RT driver currently supports only a constant vector potential following an impulse, not a laser waveform or pump–probe spectrum.

## Running the local prototype

Dependencies: Python/NumPy, Libxc 7, FFTW3, and the existing Si restart/localization artifacts. The localization initializer currently reads `calculations/si_tdcdft_k4/gs/data_for_restart` and `docs/results/si-time-wannier/mlwf_initial.npz`. This is a case-specific research driver, not a standalone general SALMON interface. Benchmark scripts also contain the local snapshot paths used for this experiment.

Create an export directory and set `SALMON_HSE_REFERENCE_EXPORT` to it when running a supported serial ground-state calculation. Export requires a complete marker written last. Without the environment variable, the native exporter returns immediately.

From the repository root:

```sh
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 samples/hse_mlwf_reference/verify_native_export.py /path/to/export
python3 samples/hse_mlwf_reference/scf.py /path/to/export /path/to/scf-output
python3 samples/hse_mlwf_reference/rt.py /path/to/export /path/to/scf-output/state.npz /path/to/rt-output 1 0.08 0.0001
python3 -m unittest discover -s samples/hse_mlwf_reference -p 'test_*.py'
```

The RT state requires the adjacent SCF `result.json` to certify convergence with no pair truncation. Results must be checked for the requested duration; a short pilot is not a long-time stability or spectral validation. Export binaries are same-machine native-endian streams, with metadata guards. Large restart data are not versioned.

Accuracy, timing, and scope of the current experiment are recorded in `docs/results/si-hse-mlwf/README.md`.

The optional ACE extension and self-consistent implicit-midpoint pilot are described in [README_ACE.md](README_ACE.md). They retain fresh full-exchange residual checks; they do not enable native SALMON hybrid dynamics or implement PT-CN yet.

Parallel-transport Crank–Nicolson with ACE is a separate option documented in [README_PTCN.md](README_PTCN.md). It preserves the full-exchange endpoint residual gate and reports its norm drift without renormalizing the solution.

For longer constant-A response runs with validated restart and live status, see [README_EXTENDED.md](README_EXTENDED.md). The authoritative checkpoint embeds both wavefunctions and history; failed trial endpoints are not saved as accepted states.

The fixed spatial-support experiment is documented in [local-support results](../../docs/results/si-hse-local-support/README.md). It has a consistent fixed-gauge energy gradient, but its orbital-specific cutoffs fail the common-Hermitian ACE gate. It is not integrated into the propagation driver.
