# H₂ PBEh(40) supercell scaling

Run the SALMON executable itself to converged ground states. This is a periodic
H₂ array, not the synthetic exchange-kernel probe. No production code changes
or memory caps are required.

```sh
python3 benchmarks/benchmark_h2_pbeh/run.py \
  --binary /absolute/path/to/salmon \
  --output /absolute/path/to/new-results-directory --repeat 3
```

An MPI/HSE-enabled GNU build is used; the recorded run also enables ScaLAPACK.
No DC/LCFO calculation is performed here, so LCFO's ScaLAPACK solver is not
measured. The launcher defaults to `mpiexec --bind-to none` (OpenMPI); override
`--mpiexec` for other launchers or placement. `rank_probe.py` supports OpenMPI
and PMI rank environment variables on macOS/Linux. Each SALMON child inherits
the MPI environment and stdin. Its `wait4` peak RSS excludes the Python wrapper
and launcher. An output directory containing results is never overwritten.
`--pilot` runs only 4×4×1 at MPI16 once. `--timeout` is a per-launch timeout.

## Fixed physical and numerical settings

- One x-oriented H₂ per 8×8×8 bohr primitive cell; bond 1.4 bohr.
- Grid spacing 0.5 bohr (16³ points per primitive cell), Γ only, occupied states
  only (one doubly occupied orbital per molecule), `xc='pbeh40'`, no rVV10.
- Coulomb cutoff **4 bohr for every shape**, avoiding a changed cutoff in 2×2×2.
  Full source support (`exx_mlwf_radius=0`), MLWF interval/maxiter/tolerance
  5/100/1e-7. This does not benchmark finite-radius sparsity.
- SCF `rho_dne` threshold 1e-10, maximum 2000 iterations, four CG steps,
  `method_init_wf='gauss10'`, Broyden `alpha_mb=0.1` (default 0.75 diverged for 16×1×1),
  no restart output. FFTE Hartree; distributed exchange uses FFTW pencils.
- Repository H pseudopotential (lmax=1, lloc=1); each run retains its input and
  pseudopotential. This grid/cutoff is a performance fixture, not a convergence
  study of periodic full-range PBEh. Supercell Γ sampling also changes the
  primitive-cell equivalent k sampling as shapes grow.
- OpenMP and common BLAS thread controls are all one.

## Matrix of cases

Weak scaling uses one MPI rank per molecule: 1×1×1, 2×1×1, 4×1×1, 8×1×1,
16×1×1, 2×2×1, 4×4×1, 2×2×2. Strong scaling fixes 4×4×1 (16 molecules),
with MPI 1, 2, 4, 8, 16. The common 4×4×1/MPI16 case is measured once per repeat
and used in both series (12 distinct cases, 36 launches for three repeats).

All runs use one orbital group and y/z spatial ranks:

| MPI | `nproc_rgrid` |
|---:|---|
| 1 | 1,1,1 |
| 2 | 1,2,1 |
| 4 | 1,2,2 |
| 8 | 1,4,2 |
| 16 | 1,4,4 |

Thus weak scaling keeps the number of owned grid points at 4096 per rank.
However, all occupied columns and dense band matrices still grow with system
size. Full-support exchange evaluates orbital pairs. Constant physical volume
per rank does **not** imply constant work or constant wavefunction storage.
The x dimension is never divided, including for 16×1×1. This is a property of
the present x-complete pencil implementation, not a general best decomposition.
MPI1 selects the legacy single-rank Wannier path; MPI>1 selects spatial EXX.
The measured speedup therefore includes this implementation switch. None of
these cases independently measures orbital distribution or DC decomposition.

## Timing, memory and validation

Runs are sequential and start fresh processes. `results.json` records the binary
hash, compiler build cache, MPI version, exact benchmark sources, pseudo hash,
settings, raw runs and median/min/max summaries. Each run directory retains the
full stdout, generated input, info/eigen files and per-rank measurements.

The main timing is SALMON's **maximum-rank `scf iterations` timer**. It includes
final forces and any output in that region, not just minimizer arithmetic.
Divide by the actual completed iteration count in `h2_info.data`, not the
printed convergence index (which is one larger). Also record root total
calculation time and launcher wall time (including output parsing); keep these separate. Per-iteration
averages include initialization effects within SCF and periodic MLWF refreshes,
and are not a fixed-state microbenchmark.

For weak efficiency use T(1 molecule,1 rank)/T(N molecules,N ranks). For strong
speedup use T(16 molecules,1 rank)/T(16 molecules,P ranks); divide by P for
parallel efficiency. Report per-iteration variants separately because iteration
counts can differ. Median ratios are summaries, not uncertainty estimates.

Peak RSS is the child's lifetime high-water resident memory, in bytes (macOS
native bytes, Linux `ru_maxrss` converted from KiB). Use the maximum rank peak.
It includes libraries and MPI; it is neither array-only storage nor simultaneous
node RSS. Small jobs can be dominated by process/communication overhead.

A successful case requires normal SALMON completion, convergence below 1e-10,
finite energy/times, and successful rank exits. All 4×4×1 energies must agree
within 1e-5 eV across rank counts. Localization diagnostics are preserved:
nonconverged gauges at full support do not by themselves invalidate exchange,
but the result is not evidence for converged finite-radius MLWF approximations.

## Representative values on a shared workstation

The user reports contention with other applications, especially at MPI16.
For this run, prioritize the minimum observed **SCF seconds per iteration**
over three converged trials; also report minimum time to convergence. Retain
medians, ranges, and all raw trials. Minima of different metrics can come from
different trials: do not divide minimum total time by minimum iteration count.
Use the already normalized per-run time when taking a minimum. These minima
are observed best cases, not estimates of guaranteed uncontended performance.

Recorded results: [Japanese measurement note](../../docs/h2-pbeh-scaling-ja.md).
