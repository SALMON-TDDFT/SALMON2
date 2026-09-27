# Si 3D DC-HSE / LCFO RT inputs

Siの4³/6³/8³/10³用GS・RT入力セットです。RTの`&functional`で
`hse_lcfo_wf_radius=9d0`を指定済み。初期WFの球内ノルムが99.9%未満なら
Warningを出します。半径の自動変更はしません。大規模計算は未実行です。

Si conventional diamond cubic cells, lattice constant **10.26 bohr** (the same
value as the existing Si64 chain input). These are generated weak-scaling inputs,
not converged numerical results. No large jobs have been submitted.

| Cells / DC fragments | Atoms | MPI ranks (GS and RT) | Occupied RT states | Global mesh |
|---|---:|---:|---:|---|
| 4×4×4 | 512 | 64 | 1024 | 64³ |
| 6×6×6 | 1728 | 216 | 3456 | 96³ |
| 8×8×8 | 4096 | 512 | 8192 | 128³ |
| 10×10×10 | 8000 | 1000 | 16000 | 160³ |

Each case contains `atom.dat`, `gs/inputfile`, `rt/inputfile`, and an RT restart
symlink `rt/data_dcdft -> ../gs/data_dcdft` (initially dangling, populated by GS).
The root contains `Si_rps.dat`, generator, validator and manifest.
`gs-env.sh` and `rt-env.sh` are no-op compatibility files; do not source them. Preserve this directory layout when extracting/copying.

Common settings: 8 atoms per core, 16³ core mesh; buffer 8 grid points per side,
so each periodic DC fragment has 32³ mesh and 64 atoms. `nstate_frag=256`,
Gamma, k/orbital MPI=1. Thus rank count is fragment count, not node count.
Use the same MPI/OMP placement policy at all sizes. The supplied Si pseudopotential
comes from `samples/dc_hse/si64-chain/Si_rps.dat`, Z=14, local channel=2.

DC-SCF uses HSE06, simple mixing 0.01, one CG iteration, 1800 maximum SCF iterations
and threshold 1e-7. LCFO eigensolver is CheFSI (degree60, max200, residual1e-7).
GS state count is 32 per core, RT occupied states 16 per core. These SCF settings
are a starting protocol; convergence has not been established for the 3D cases.

RT: Taylor4 direct WF, dt=0.02 a.u., 16 steps, x impulse1e-4, ACE1/U1,
FFT batch1/ESTIMATE. This short run tests correctness/performance, not a converged
dielectric spectrum. Existing GS/RT density differences are retained.

## Radius and norm warning

The RT input now includes:

```fortran
&functional
 xc='hse06'
 yn_hse_lcfo_rt='y'
 yn_hse_wannier='y'
 yn_hse_lcfo_direct_wf='y'
 yn_hse_lcfo_seed_distributed='y'
 yn_hse_lcfo_continuity='n'
 yn_hse_lcfo_fft_measure='n'
 hse_lcfo_ace_interval=1
 hse_lcfo_u_interval=1
 hse_lcfo_fft_batch=1
 hse_lcfo_wf_radius=9d0
/
```

`hse_lcfo_wf_radius` is always in **bohr**, independent of `unit_system`.
0 means full support; positive values specify the fixed global periodic sphere.
The default is 0 (full support); negative values are rejected. Algorithm settings
are read only from the namelist. The former `SALMON_LCFO_*` variables are obsolete.
`yn_hse_wannier='y'` enables MLWF sources for LCFO RT; the same flag also
selects Wannier exchange in the non-LCFO HSE path. DC-SCF enables it internally.

Initial sphere coverage is evaluated separately for each WF using |w|² on disjoint
core grids and 3D minimum-image distances. If any WF has sphere norm / total norm
<0.999, root prints `WARNING LCFO MLWF radius: sphere norm below 99.9%` with the
count; minimum coverage and WF index are also logged. `lcfo_mlwf_radius.dat` saves
all initial total norms, sphere norms, fractions and protected flags. Protected
WFs have unreliable centers and remain uncut, even if their geometric fraction
is below threshold. Radius and normalization are not changed automatically.
This checks the initial WFs only; the existing RT total discarded-norm diagnostic
continues. 99.9% norm retention is not a bound on current or dielectric error.

9 bohr is a starting radius, not a demonstrated 99.9%-coverage radius for 3D Si.
Compare separate RT runs from the same GS with radius0/full and chosen finite radii.
Do not change radius between sizes if measuring fixed-work weak scaling.

## Run order

Requires the HSE/MPI/ScaLAPACK build with the complete LCFO namelist-control update (the earlier radius-only patch is insufficient). On Fugaku,
include the `-Nalloc_assign` toolchain correction before numerical tests.
Use the existing batch allocation/launcher; the commands below describe working
directories, not a new scheduler submission script:

```sh
# For example, in 4x4x4/gs:
# Launch SALMON with 64 MPI ranks and inputfile as standard input.
# Check SCF and LCFO eigensolver convergence and GS output before proceeding.

# Then in 4x4x4/rt:
# Launch SALMON with 64 MPI ranks and inputfile as standard input.
```

Repeat with 216, 512 and 1000 MPI ranks for the other sizes. Distributed seed QR
is enabled in each RT input and requires MPI+ScaLAPACK. GS explicitly selects
`yn_hse_lcfo_rt='n'`; RT selects `'y'`. No algorithm environment setup is needed.

## Remaining memory limits

Distributed seed QR reduces root concentration; root six-link localization storage
is still 0.094/1.068/6.000/22.888 GiB for these four sizes. Dense SVD, U, ACE, basis
and wavefunctions require additional memory. Input generation is not confirmation
that 8³ or 10³ fits memory. Large production runs remain pending peak-memory checks.

Regenerate/validate without submitting any calculations:

```sh
python3 generate.py
python3 validate.py
```

## Fugaku archive-source update

The bundled patches bring the cumulative Gamma-memory source to the current
namelist-only version. Preserve the existing compiler/POSIX fixes and apply
`-Nalloc_assign` separately if not already present. Back up the source first.
Choose the starting step matching the source; do not apply an earlier patch twice.

1. If only the cumulative Gamma-memory patch is installed, check
   `patches/lcfo-radius-before.sha256`, apply `lcfo-radius-namelist.patch`, then
   check `lcfo-radius-after.sha256`.
2. If the radius namelist patch is installed but the loop cleanup is not, apply
   `lcfo-radius-loop-cleanup.patch`.
3. After the loop cleanup (`3c4ff577` / `eb23bfea` equivalent), apply the new
   namelist migration below. Its hashes check all 10 affected source files.

Each earlier patch also requires a successful `--dry-run --fuzz=0` before applying.
For step 3, from the SALMON source root:

```sh
(
  set -e
  si_inputs=/absolute/path/to/si-3d-weak-scaling
  sha256sum -c "$si_inputs/patches/lcfo-namelist-only-before.sha256"
  patch --batch --forward --fuzz=0 --dry-run -p1 < "$si_inputs/patches/lcfo-namelist-only.patch"
  patch --batch --forward --fuzz=0 -p1 < "$si_inputs/patches/lcfo-namelist-only.patch"
  sha256sum -c "$si_inputs/patches/lcfo-namelist-only-after.sha256"
  cmake -S . -B build
  cmake --build build -j 8
)
```

If a hash or dry-run differs, stop; do not force the patch. The complete patch
was applied to clean copies and the resulting files matched the expected hashes.
Local GNU/MPI build and small-system regressions pass. Fugaku compilation and
3D numerical runs of this version remain unverified. Older executables cannot
read these new namelist entries.
