# LCFO RT development status

The DC input guard remains unchanged: `yn_dc=y` currently accepts `theory=dft`
only. An experimental fixed-LCFO-subspace RT adapter now reuses native
`initialization_rt`, density, Hartree, semilocal XC, pseudopotential, Taylor4
predictor/corrector and current routines. No additional SCF is performed.
Si128/16-fragment self-consistent RT has passed a four-step integration check;
there is still no converged dielectric spectrum.

## Experimental native path

Set environment `SALMON_LCFO_RT=1` and use `theory='tddft_response'`,
`yn_dc='n'`, `yn_conventional_from_dcdft='y'` in `&calculation`.
Use Gamma, unpolarized HSE06, `yn_hse_wannier='n'`, default `hse_taylor4`,
one real-space MPI rank per fragment, and no orbital/k distribution.
Each rank's grid must coincide exactly with its fragment core. The existing
`./data_dcdft/fragments` LCFO records provide both the initial orbitals and fixed
orthonormal core bases. Native `hpsi` is projected after all Hamiltonian terms;
Hartree and density are evaluated on the whole physical system by existing code.

The current implementation requires occupied-only fixed occupations (for Si128,
256 states, 512 electrons, omit temperature). It loads the first256 saved LCFO
states without reconverging the accepted DC-to-LCFO density difference.
Restart input/output, checkpoints and time_shutdown are rejected until LCFO
basis and exchange state metadata are supported.

`hse_lcfo_rt.f90` reconstructs periodic fragment+buffer bases and density factors,
uses the existing screened FFT exchange kernel, core-weights and Hermitianizes
the projected exchange, and sums contributions once across spatial ranks.
The global trace energy is passed to the existing native energy bookkeeping.
Existing ACE is built in coefficient space; a rejected indefinite/singular
metric falls back explicitly to the full projected operator without clipping.
The predictor/corrector averages endpoint operators. By default, full support
uses density eigenfactors (not MLWFs). Its relative density eigenvalue threshold
1e-14 only removes numerical null modes, with discarded trace logged.

### Opt-in transported MLWF sources

Set `SALMON_LCFO_RT_MLWF=1`. Initial localization uses global occupied states in
the distributed fixed LCFO basis, not independently diagonalized fragment density
factors. A Gamma-specific SU(2) Jacobi optimizer minimizes the same six-link MV
functional as the existing gauge code. `hse_mlwf_maxiter` is the maximum sweep
count and `hse_mlwf_tolerance` the gradient tolerance; Si128 uses200 and1e-6.
The original general-k optimizer is unchanged. A coefficient-space pivoted trial
seed avoids gathering the complete real-space wavefunctions on one rank.

After this one localization, polar overlap transport aligns each new occupied
frame with the previous accepted WFs. This removes band dynamical phases;
merely holding the numerical U fixed would spread WFs even at stationary density.
Predictor and corrected endpoint frames share the accepted step-start reference;
trial updates are rolled back before accepting the corrected frame. Cache keys
include both coefficients and occupations. No repeated MV minimization occurs.

`SALMON_LCFO_RT_RADIUS=0` (default) retains full support. A positive value is a
halfwidth in bohr along global periodic x, around fixed INITIAL WF centers.
Distances use the whole Si128 cell, not the shorter fragment period. Centers with
circular reliability below0.1 are protected at full support. No source
renormalization, occupation cutoff or density/Hartree truncation is introduced.
Initial localization must converge before a finite-width test proceeds.

Diagnostics `lcfo_mlwf_initial.bin` and `lcfo_mlwf_links.bin` preserve the initial
occupied coefficients, unitary/centers and optimizer inputs. They are not restart
files. The active fragment-source count and global discarded source norm fraction
are logged. Existing zero-source skipping avoids unnecessary FFTs after masking;
no estimated speedup is inferred from source counts alone.

### Optical and continuity interpretation

Source-only masks are explicitly nonvariational. Hermitian assembly and ACE do
not by themselves restore local charge continuity. Set
`SALMON_LCFO_RT_CONTINUITY=1` to evaluate the exchange density source before LCFO
output projection, including both core-weighted adjoint halves and summing the
overlapping fragment buffers. Logs contain the signed integral and L1 integral
(electrons per atomic time) and maximum density-rate magnitude. Full-support Fock
cancellation should hold to roundoff; a masked case can conserve total charge
while producing a nonzero local defect. The additional diagnostic FFT work must
be held equal when comparing runtimes.

Native kinetic-plus-pseudopotential current remains the recorded observable.
Its masked transformed spectrum is an approximate-model diagnostic, not yet a
certified dielectric function: cutoff acceptance requires assessing the
continuity defect/field-current consistency as well as current and spectral
errors. Do not add an unwrapped-position exchange commutator in this periodic
fragment model as an ad hoc correction.

The environment switch is a development opt-in, not a new production input.
Without it the ordinary reconstruction and native RT paths are unchanged.


`src/rt/lcfo_rt_core.f90` supplies a fixed orthonormal complex basis kernel:

- `lcfo_cayley_step(H,C,dt,next,status)` solves
  `(I+i dt H/2) next=(I-i dt H/2) C`. It does not diagonalize H or normalize C.
- `lcfo_density(C,f,P,status)` forms `P=C diag(f) C†`, retaining complex coherence
  and fractional occupations. f is the physical occupation; no implicit factor2.
- `lcfo_midpoint_step` iterates a Hamiltonian callback on `(P_start+P_end)/2`
  at `time+dt/2` to the specified relative Frobenius density residual. Each trial
  propagates from the unchanged starting C. It is a small reference solver;
  fixed-point convergence is not guaranteed for arbitrary dt or functionals.
- `lcfo_grid_density` reconstructs density on a core's grid from its basis and
  the corresponding coefficient rows of ALL global occupied states.
- `lcfo_project_potential` projects one real grid potential using grid quadrature.

Outputs use separate storage from inputs. Failure retains the input C in `next`
when its shape permits; status1 invalid input,2 non-Hermitian H,3 linear solve
failure,4 fixed point not converged,5 callback failure. A caller must reject the
step on any nonzero status. No silent symmetrization of caller H is performed.

The independent Python reader `samples/dc_hse/lcfo_rt_reference.py` accepts only
Gamma, unpolarized V1 complex LCFO data. It verifies completion/footer lengths,
run identity, geometry, periodic halo mapping, orthonormal core bases and saved
eigenvectors. It reconstructs the native symmetrized matrix, retaining a separate
raw-directed anti-Hermitian diagnostic. Basis arrays use Fortran grid ordering.
The `lcfo_frozen_probe.f90` helper tests real saved LCFO eigenstates against exact
Cayley phases and time reversal; this is deliberately NOT a TDHSE driver.

Still required: longer-time/dt/LCFO-basis convergence, zero-field subtraction,
field-current assessment of the source-mask approximation, and dielectric comparisons.
The accepted initial density mismatch is retained; stationarity is not a gate
for starting RT.
A Hermitian matrix alone establishes neither energy-functional consistency nor
optical-response accuracy. In particular, the static source-cutoff scripts are
not silently promoted into a variational time-dependent functional.

Tests: CTest `lcfo_rt_core`; `OPENBLAS_NUM_THREADS=1 python3
testsuites/unit_lcfo_rt/test_reference.py`. The standalone check.py compiles
against the local Homebrew BLAS for the current development machine; CTest uses
the build's selected BLAS/LAPACK and is the portable verification path.

Native integration regression: `python3 testsuites/unit_lcfo_rt/test_native.py
--binary /absolute/path/to/salmon --pseudo /absolute/path/to/H_rps.dat`.
It runs MPI2 jobs sequentially in a fresh temporary directory: Gamma DC-SCF,
LCFO RT, half-dt RT, and rejection of unsupported restart. It checks normal
completion, finite output, endpoint-current agreement and post-impulse energy
width. This small integration test is not a production convergence criterion.

Additional regressions: `testsuites/unit_hse_wannier/test_gamma.py` (known
anisotropic Gamma solution and noncommuting links), and
`testsuites/unit_lcfo_rt/test_transport.py` (stationary-density phase invariance,
predictor rollback and corrected-state cache acceptance). These standalone tests
currently use the local Homebrew BLAS path. `test_native.py` also checks full
source parity, large-radius identity, initial-U identity across support cases,
and full/masked exchange continuity diagnostics with MPI2 jobs run sequentially.

### ACE reuse across physical time steps (2026-09-26)

`SALMON_LCFO_RT_ACE_INTERVAL=N` accepts a positive integer (default1).
For an impulse, step1 always rebuilds ACE **before the first predictor**, as
well as at predicted/corrected endpoints; later refresh steps are1+N,1+2N,….
For a smooth laser starting at zero field, initial ACE is retained and endpoint
refreshes occur at N,2N,… . The present LCFO path constructs initial ACE from
loaded GS orbitals; it does not deserialize fragment GS ACE factors.

Between refresh steps, only the nonlocal exchange operator is frozen. Density,
Hartree and semilocal XC continue to update through native routines. MLWF frames
still undergo polar transport for each predictor/corrected state. An invalid ACE
or changed occupations bypass reuse. Unchanged endpoint operators retain a
single set of ACE factors instead of doubling their rank by averaging duplicates.
Exact-coefficient cache keys are invalidated when a retained step transports WFs.

The skipped-refresh exchange energy is evaluated as0.5 sum f<C|K_ACE|C> at current
orbitals, explicitly logged as a frozen-operator trace diagnostic. It is neither
the instantaneous self-consistent HSE energy nor a conserved frozen-Hamiltonian
energy; compare its width only as a diagnostic. Local-field updates are not skipped.
The existing continuity diagnostic evaluates rebuilt exchange only and cannot
certify continuity of a retained ACE at evolving orbitals.

MPI2 regression checks default/interval1 identity, fewer builds at2/4, identical
initial MLWF data, finite output, invalid-interval rejection, impulse pre-predictor
rebuild and smooth-laser first-step reuse. Si128 wall-time/current comparisons
are recorded below once measured; a short pilot does not certify dielectric
accuracy or select a production refresh interval.

<!-- ACE_REUSE_20260926 -->
## ACEを複数ステップで再利用した測定（2026-09-26）

実装 `a9e3716`。Si128・16フラグメント、MLWF全範囲、MPI16/OpenMP1/BLAS1、16ステップ、dt=0.16 a.u.（合計0.06192 fs）。同じ初期状態・U・励起1e-4 a.u.、1計算ずつ実行。連続性診断は全条件で無効とし、診断コストを含めずに比較。

impulseでは第1ステップの最初の予測伝播前に必ずACEを再構築し、そのステップの予測・修正状態でも更新する。以後は1+N、1+2N…ステップで更新。再利用中も密度・Hartree・半局所XCとMLWFの極分解による位相追跡は継続する。滑らかなレーザーでは初段更新を省き、GS軌道から初期化時に構築したACEを再利用できる（保存済みGSのACE因子の読み込みではない）。

| 更新間隔 N | ACE構築回数 | 全実測時間 s | 全時間の高速化率 | RT部分 s | RT部分の高速化率 | x電流最大差 / 基準ピーク | 診断エネルギー幅 Ha |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 35 | 336.73 | 1.000倍 | 324.19 | 1.000倍 | 0.00000% | 3.73290e-08 |
| 2 | 19 | 219.16 | 1.536倍 | 203.07 | 1.596倍 | 0.31990% | 1.67213e-04 |
| 4 | 11 | 153.20 | 2.198倍 | 136.51 | 2.375倍 | 0.96453% | 4.29557e-04 |
| 8 | 7 | 119.28 | 2.823倍 | 103.02 | 3.147倍 | 2.26090% | 7.39916e-04 |

旧版と更新間隔1の最初の8ステップの電流差は最大 5.477e-18 a.u.。小規模試験で更新間隔1と既定動作の一致、impulseの予測前更新、レーザー初段再利用、無効間隔の拒否を確認。既存回帰7件とMLWF位相追跡試験も通過。

再利用中の交換エネルギーは、現在の軌道と保持したACEによる半トレースの診断値。瞬時の自己無撞着HSEエネルギーや保存すべき凍結Hamiltonianエネルギーではないため、この幅を物理的なエネルギー保存誤差と解釈しない。

電流差は無励起差引き前の同一初期状態の軌道比較。短時間・各条件1回の測定であり、速度の統計誤差、長時間の誘電関数の誤差、採用可能な更新間隔は未確定。既存の局所連続性診断は交換を再構築した時点を測るものであり、再利用中ACEの連続性を検証したことにはならない。長時間の範囲比較は更新間隔1を基準として維持する。
