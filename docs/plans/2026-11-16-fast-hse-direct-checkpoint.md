# Fast-PR HSE direct persistence

User direction: commit prior validated changes, then implement direct numerical-state save/restore for fast-PR HSE, with fresh-process continuation equality and no elapsed propagation replay. Conventional HSE and comparison baselines are not execution engines.

Implementation uses the existing serial internal Fortran tddft_context, not the standard semilocal host-MPI middlepoint entry. Separate hse_checkpoint module saves numeric ACE/history/cache state while retaining freshly initialized native FFT resources. The outer context serializes orbital ping-pong buffers, density/potentials/Hartree/RT/Poisson/energy histories, accepted step and fields. Check same settings/grid/cell/dt/capacity/initial energy, and require caller runtime/adapter/GS/archive hash authentication. Do not copy a context or opaque resource owner. Late restore failure quarantines further stepping.

Acceptance completed:2/3k,4 fast-PR modes,split6/7,OMP1/2;32 positive fresh-process runs and48 negative truncated/settings/trailing-byte runs. Compare API current/energy and listed internal numerical trace; verify6/5 remaining step calls. Prior8 two-context independence/lifetime regressions preserved. Package adapters/tools/probe and synchronize runtime JP/EN/AI guides. Full evidence: docs/fast-hse-direct-checkpoint-validation-2026-11-16.md.

Bounds: synthetic serial CPU full-support scalar/fixed-ion/occupation profile. Not host-MPI/FDTD-coupled HSE, stable Python/C API, GPU/symmetry/distributed/localized/general portable archive. Step-zero persistence not established. Physical convergence remains user-owned.
