# Full TDDFT context validation — 2026-11-15

## Result and scope

An internal Fortran context now owns all arguments of native initialization/time_evolution_step, including orbitals, densities, potentials, pseudopotential work, XC/HSE owners and per-state output handles. Two states are initialized independently from the same synthetic restart and advanced with different endpoint A/c waveforms. No active global-state pointer or initialized-context copy is used.

Scope is serial CPU/FFTE, periodic unpolarized fast HSE, fixed ions, zero initial A, transverse fields, common immutable process settings. Different grid/k/XC/dt per state, concurrency, MPI, general failure atomicity, a stable Python/C ABI and physical convergence remain unverified/unimplemented. Shared XC scratch and optimization tables still exist; basic captured settings are guarded, not fully extracted. Energy sampling is required at every step.

## Local checks

- 8 configurations (k=2/3; exact ACE Taylor, learned Taylor, direct Taylor, learned PT-CN), each at 1 and 2 OMP threads, 12 steps.
- Separate-process A/B versus interleaving and reversed order: current/energy, returned API values and binary orbital/density/potential/HSE-history traces match exactly. A/B traces are confirmed different.
- Close/reopen either state halfway through: the surviving state matches its uninterrupted separate-process reference, and the reopened state matches a fresh process prefix.
- Invalid A and changed dt are rejected before advancement; interleaved traces still match clean isolated references.
- Native unchanged full-driver Acos2 comparison: max current difference 1.0587911840678754e-22; max energy-table difference 0; required rtol=1e-12, atol=1e-14. Returned charge current has the opposite sign to native matter current.
- The native reference rejects arbitrary sampled-field input. That failed attempt was retained; supported Acos2 was used without modifying the reference.
- A deliberate external mutation removing per-system HSE binding fails on interleaved a_rt.data. This confirms the test detects real full-step contamination.
- Source Python suite: 305 pass. Installed suite: 305 run, one optional plotting skip, no failure. sdist-to-wheel passes. Final installed context/probe/runner match repository hashes; external-output guard passes before any write.
- Source and all linked application object audits precede every new probe link. The audited full runtime is unchanged.
- No GitHub Actions launched. Previous HSE owner work was pushed as c094f726 with [skip ci]; REST readback showed no run for that commit. This full-context change remains uncommitted.

## Evidence identities

Context SHA256: `2b648162b98f69641760b01a55d4fe71346872db26741929b3e1636d35b83764`

- 1 OMP result SHA256: `24b4ef8292d0dcdb23ef7ddc7e848f009c453d221a7a522661f27bf2f3bcf5e5`
- 2 OMP result SHA256: `af8fcca21c2ba873a25faf9802e9ebd279c343f57e98cf11810c5962e8ccf158`

Evidence lives under /tmp/salmonpy-full-context-verified1-20261115, /tmp/salmonpy-full-context-final2-20261115, and /tmp/salmonpy-full-context-mutation-20261115. Synthetic traces are diagnostics, not restart files; no real GS payload was collected.
