# Conventional hybrid ground state to native real-time propagation

> 2026-10-06：本記録中の`developer_tests`は当時のローカル開発検証です。GitHub配布から除外しました。通常の回帰試験は`testsuites`を使用します。

PBE0, PBEh(40), and PBEh(40)+rVV10 can start fixed-ion native real-space RT
from an ordinary converged GS. This route does not require divide-and-conquer
(DC) preparation or LCFO reconstruction.

Use the same physical system, functional, k mesh, and hybrid parameters in the
GS and RT inputs. The supported functional names are `pbe0`, `pbeh40`, and
`pbeh40_rvv10`.

## Ground-state producer

```fortran
&calculation
 theory='dft'
 yn_dc='n'
/
&control
 write_gs_restart_data='wfn'
 yn_restart='n'
/
&functional
 xc='pbeh40'
 exx_mlwf_radius=0d0
/
```

Use fully occupied, spin-unpolarized periodic orbitals: `nstate=nelec/2`, with
no finite-temperature occupations or empty states. Converge the GS before
starting RT. Gamma spatially distributed GS export uses `write_gs_restart_data='wfn'`.
The producer writes the usual `data_for_restart` wavefunction files together
with `hybrid_gs.bin` metadata.

## Real-time consumer

```fortran
&calculation
 theory='tddft_response'
 yn_dc='n'
 yn_conventional_from_dcdft='n'
/
&control
 directory_read_data='../ordinary_gs/data_for_restart/'
 yn_restart='n'
/
&functional
 xc='pbeh40'
 exx_mlwf_radius=0d0
/
&tgrid
 dt=.02d0
 nt=16
/
&emfield
 trans_longi='tr'
 ae_shape1='impulse'
 e_impulse=1d-4
 epdir_re1=1,0,0
/
```

The time step above is in atomic units when `unit_system='a.u.'`; otherwise use
the selected input units. Omit `propagator` to use the supported default hybrid
Taylor4/ACE propagation. MLWF and ACE factors are rebuilt from the loaded GS.
For a finite pulse use `theory='tddft_pulse'`, `ae_shape1='Acos2'`, and the
corresponding pulse amplitude, frequency, duration, start time and polarization.

Supported parallel layouts are Gamma spatial decomposition and full uniform
multi-k meshes with k-point parallelism only. Symmetry-reduced meshes and
simultaneous multi-k spatial decomposition are not supported by this route.
The numbered examples use full exchange support and two MPI ranks. Gamma RT
also accepts adaptive MLWF support. Switching from full-support GS to truncated
RT changes the approximation and may introduce an initial transient; compare
with full support before using the result for spectra.

## Data validation and scope

New conventional GS data must contain `hybrid_gs.bin`. Older ordinary GS files
without this metadata are rejected for this route. The versioned metadata
validates the functional and exchange fraction, screening/effective Coulomb
radius, rVV10 settings, cell/grid/k mesh, ionic and pseudopotential configuration,
and occupations. The stored `info.bin` and `occupation.bin` are checked before
wavefunctions are read, and loaded occupations are checked again.

Copy the complete GS output directory; do not mix files from different runs.
Changing the functional, physical system, occupations, or validated parameters
requires a new GS. Existing wavefunction binary formats and the legacy HSE/DC
routes are retained.

This is a fresh RT calculation initialized from GS, not RT checkpoint
continuation. Conventional-route moving nuclei, RT continuation, and fixed-radius
exchange truncation are outside the supported scope.

## Small regression cases

- 431 → 432/437: PBE0 full two-point k mesh, impulse and Acos2, two k ranks.
- 433 → 434: PBEh40 Gamma GS and Acos2 RT, two spatial ranks.
- 435 → 436: PBEh40+rVV10 Gamma GS and impulse RT, two spatial ranks.

These CPU/MPI/HSE-enabled CTest cases generate fresh GS data and require its
verification before preparing RT. They check completion, finite currents and
energies, and electron-number conservation. `developer_tests/653_functional/test_conventional_rt.py`
also compares one-rank and two-rank propagation for all three functionals,
checks zero-field stability and Gamma adaptive 99.9% source ACE, and rejects
missing/mismatched saved data.
Set `SALMON_TEST_EXE` and `SALMON_TEST_MPIEXEC` to run that integration suite.
These small tests do not establish k-mesh, time-step, or spectral convergence.

## Common full-support k-point exchange

Ordinary fully occupied multi-k GS and native RT now share `exx_k_exchange`
with HSE06. Density tiles are transposed across k ranks; the full orbital set
is not gathered to the k-root for exchange. Rectangular orthorhombic grids are
accepted. PBE0/PBEh retain the analytic spherical-Coulomb G=0 term and the existing
`pbeh_coulomb_radius` convention (zero selects half the shortest BvK cell length).
This is a kernel choice, not an MLWF support radius.

DC, fractional/extra-state GS, explicitly localized support and requested HSE
Wannier snapshots retain the occupation-aware Wannier/spatial paths. Gamma
spatial exchange is unchanged. Legacy `hse_block_rows`, `hse_fft_layout` and
`yn_hse_profile` input names still control the common k engine for compatibility.
The `EXX_DISTRIBUTED_K` log marker identifies the shared path.
