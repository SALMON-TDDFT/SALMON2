# Non-PZ direct checkpoint validation — 2026-11-16

Base1851c533. Native checkpoint format V2 extends the direct numerical state route to builtin TB-mBJ (automatic/fixed cval), Libxc-PZ and Libxc-PBE. All prior field/PML/LD/SBE/TDDFT dynamic arrays remain saved. New XC metadata binds XC/x/c names, internal type/spin, derivative flags and input cval; scalar cval is restored only after collective validation. Libxc handles remain those created by fresh initialization and are never serialized. Builtin `pbe` is unavailable; use `libxc_pbe`. HSE is not admitted by this native checkpoint profile.

## Preparation and environment

Local Apple Silicon, Python3.14, OpenMPI, OMP/BLAS1; builds use two workers, tests at most four cores. Libxc4.3.4 official GitLab tag resolved to c45ee4310327089aab3e1379c7fc158a510ec6ec, built in an independent external prefix with GCC/gfortran15. Existing Homebrew Libxc remains unchanged. This dedicated installation supplies SALMON's xc_f90 interface.

New source prepared in `xc-direct-source`, native host-MPI build `xc-direct-build-v2`, separate serial reference build `xc-direct-serial-build`. Source/object forbidden-eigensolver audits pass. The serial reference executable is intentionally MPI-OFF; an MPI-linked executable is unsuitable as a child of the outer MPI sample. Preliminary failed builds/runs remain separate and are not accepted evidence.

`tools/prepare_xc_checkpoint_fixture.py AUDITED_EXECUTABLE NEW_ROOT INITIAL_GUESS_RESTART` prepares each functional's own Si SCF and matching RT input. Same16³ grid,2³k,8-atom/32-electron starter. PZ wavefunctions are an explicit initial guess, followed by SCF using the target functional. Convergence threshold1e-10: automatic TB-mBJ45 iterations, fixed cval1.2 TB-mBJ50, Libxc-PZ2, Libxc-PBE41. The direct-from-simple-guess TB-mBJ attempt produced NaN and was rejected; only the finite, converged warm-start fixtures are used. The preparation tool now requires convergence and finite wavefunction/log data.

These are small numerical continuation tests; SCF convergence is not grid/k-point/time convergence or converged optical response.

## Acceptance

Each case runs native field with two local LD materials, one SBE group and one TDDFT group,128 steps atdt.003. Continuous execution is compared to a separate save job stopping at63, then a fresh job loading internal state and performing65 propagation steps. The sample records restored_step=63/native_steps_this_job=65. Independent serial electronic response replay remains a separate validation; it does not reconstruct the native resumed trajectory.

The root comparison script explicitly asserts bit equality of complete A/c,E,j arrays, in addition to the sample's existing tolerance comparison and clocks/lag1/serial-response checks. All five cases passed with zero maximum differences:

| XC | Saved step | Resumed propagation calls | Max differences A/c,E,j |
|---|---|---|---|
| PZ regression | 63 | 65 | 0,0,0 |
| TB-mBJ automatic cval | 63 | 65 | 0,0,0 |
| TB-mBJ fixed cval1.2 | 63 | 65 | 0,0,0 |
| Libxc-PZ | 63 | 65 | 0,0,0 |
| Libxc-PBE | 63 | 65 | 0,0,0 |

All use three outer ranks (field1,SBE1,TDDFT1). Independent serial SBE and electronic currents also match exactly for each case. This campaign does not newly validate internally parallel non-PZ TDDFT groups. Full machine-readable results: `xc-direct-comparison.json` under the persistent evidence root.

Source and installed-wheel suites each446 tests,1 skipped,no failures. Sdist-to-wheel construction and installed MPI example are exercised. English/Japanese/AI documentation synchronized. Additional review found no important issue in XC metadata restoration or the retained fresh Libxc handles. No GitHub Actions used.

Persistent local evidence under `~/.codex/outputs/salmonpy-state-inputs-20261116`: `xc-direct-fixtures-warm`, `xc-{pz-v2,tbmbj-auto,tbmbj-fixed,libxc-pz,libxc-pbe}-{continuous,save,resume}`, their `.log` files, `xc-direct-{source,installed}-tests.log`, Libxc source/build/install logs, native/serial build directories and `xc-direct-dist`. Paths are developer evidence locations, not user prerequisites.

## Limits

Same binary/environment/topology, periodic unpolarized transverse middlepoint, no MD/projection/unfold/singlescale. Opaque resources rebuild from authenticated inputs. Other XC names, HSE resources, changed decomposition and generic user-backend archives remain unsupported. Native historical text-file positions are not restored. Aggregate RSS/memory caps and production scaling remain unmeasured. Molecular work stays paused.
