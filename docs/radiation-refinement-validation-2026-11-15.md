# One-way radiation: finite-window refinement study

## Reference and fixed quantities

A prescribed physical charge dipole p_z(t)=a sin^4(pi(t-t0)/T) for t0<t<t0+T,
zero otherwise, with a=1e-3 e*bohr, t0=.03 au and T=.12 au. Source position is
(0,0,0); observation position is (6,0,0) bohr in all cases. The first interval is
intrinsically constant, so this study adds no extra sample hold/delay. The native
source/observation macro-cell to Yee maps are unchanged.

At retarded time u=t-r/c, the full transverse point-dipole reference is
E_z(t)=-p(u)/r^3-p'(u)/(c r^2)-p''(u)/(c^2 r), with c=137.035999084.
The near and induction terms are retained: this observer is not assumed far away.
The continuum total radiated energy is W=2/(3c^3) integral p''^2 dt,
which for this pulse equals 8*pi^4*a^2/(3*c^3*T^3)=5.8414474429781124e-8 Ha.
Analytic first/second derivatives were checked against finite differences and
this closed-form integral against independent numerical quadrature.

Reference background: [MIT, electric dipole fields and radiation](https://web.mit.edu/6.013_book/www/chapter12/12.2.html).
The expressions above are written in the Gaussian atomic-unit conventions used
by the native vacuum E/H solver, not SI.

## Results

Relative field error is the sampled L2 norm of E_z-error divided by that of the
reference, over 0<=t<=.30 au. Source work is evaluated independently as
-dt sum_n,points V*j(n+1/2) dot (E(n)+E(n+1))/2. This is **source work**, not a
measurement of Poynting flux through a closed surface, nor an exact proof of the
solver's discrete energy identity.

| dx (bohr) | dt (au) | Grid | Relative E error | Source-work vs continuum-energy difference |
|---|---|---|---|---|
| 2 | .003 | 41³ | 44.8151% | 13.9566% |
| 2 | .0015 | 41³ | 46.4888% | 15.7603% |
| 1 | .0015 | 81³ | 8.00055% | 2.48394% |
| .75 | .0015 | 109³ | 4.28259% | 1.14449% |
| 1 | .00075 | 81³ | 8.28731% | 2.83798% |

At fixed dt=.0015, both field error and energy discrepancy decrease with dx.
At fixed dx, halving dt does not monotonically improve agreement with the
continuum result; do not label this temporal convergence. The sampled-field
change on halving dt is 2.33307% at dx=2 and .373980% at dx=1, each normalized by
the reference field norm. Possible cancellation between spatial and temporal
errors must not be mistaken for a better converged coarse calculation.

The finest observed 4.3% field discrepancy is still material. These tests do not
establish production precision or an emission spectrum of a converged molecule.
The domain's coordinate extrema are approximately +/-40 bohr (40.5 for dx=.75);
PML thickness is eight cells, so its physical width changes with dx. The short
window is chosen ahead of nominal boundary-return times, but PML/domain-size
independence has not been separately demonstrated. Source deposition width also
shrinks with dx, as appropriate to this point-source refinement, not a molecule
of fixed finite physical size.

## Self-field observation

At fixed dt=.0015, the maximum sampled source E_z is 1.53875e-4, 1.09323e-3,
and 2.81724e-3 au for dx=2,1,.75 respectively. This strong source-grid dependence
prevents treating raw source-cell E as a validated molecular driving field.
The study does not determine a radiation-reaction/local-field prescription.
See plans/2026-11-15-molecular-self-field.md for the next design boundary.

## Reproduce and resources

```sh
python tools/test_dipole_refinement_mpi.py --check-reference
mpiexec -n 2 python tools/test_dipole_refinement_mpi.py BUILD_JSON NEW_OUTPUT DX DT
python tools/summarize_dipole_refinement.py NEW_JSON COARSE TIME SPACE FINE FINE_TIME
```

Five runs used an explicitly verified existing embedded library SHA256
1e9abfa895dc216399377c164c0cbefa606c436c7e7eee920f4aa8901fbdae2a, one field rank
plus one label-only rank, OMP/BLAS threads=1. Each was launched sequentially with
run_bounded and a 55-second timeout. No native build, network acquisition or
GitHub Actions were run. The tool bounds dx>=.75 and dt>=.00075 to limit this
study to at most 109³ cells and 400 steps. This is not a hard memory cap.
A separate ps-based RSS inspection was denied by the local sandbox; computation
and analytic checks succeeded, but peak memory was not measured. Do not report
resource-inspection failure as calculation failure or infer measured 4-GiB use.
