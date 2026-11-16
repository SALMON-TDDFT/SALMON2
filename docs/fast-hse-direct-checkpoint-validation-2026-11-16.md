# Fast-PR HSE direct internal-state checkpoint — 2026-11-16

Base:b475288d. Only the pinned fast-PR HSE engine is used, through the previously audited explicit-routing build. No conventional HSE or comparison baseline is introduced as an execution engine. Numerical kernels and existing immutable compiled binaries were not modified. New adapter modules/probes were explicitly compiled and linked into new external outputs from the audited fast objects, after source/object forbidden-eigensolver checks.

## Implemented

- `salmonpy/native_adapters/hse_checkpoint.f90`: numeric-only HSE payload. Four ACE owners (including packed numerical representation), factor/PT-CN histories and coefficients, sources/actions/occupations, flags/counters, exchange energy and pair-FFT count. Validate initial native-kernel dimensions, mesh, phases and static numeric layout. Preserve fresh FFT handles/accelerator ownership/communicators; reject unsupported accelerator/symmetry/distributed ACE/reduced-q profiles.
- `tddft_context.f90`: `rt_checkpoint(ctx,path,restore,status)` and `rt_accepted_step(ctx)`. Both orbital buffers/scratch, densities/local/XC/Hartree histories, RT and Poisson numeric arrays, energy including original E_tot0, system fields, input/run config and Mit/accepted are saved. Rebuild nonlocal phase caches after loading. Rebuild static exchange resources if needed; invoke no time evolution during restoration.
- Require fresh initialized original GS/input for restore, same native format/binary/configuration and caller authentication of all assets/archive. Metadata mismatch refuses restoration; a late failure after HSE loading marks the context unusable for further rt_step (close remains valid).
- `tools/test_hse_direct_checkpoint.py` packages hash receipts binding checkpoint, runtime, compiled probe, input, initial restart files and pseudopotential. New native driver starts a NEW process for each continuous/save/resume, loops from the accepted saved step and records actual remaining calls. Existing context test tool compiles the new dependency too.

## Measured results

Apple Silicon/macOS, gfortran/OpenBLAS/fast-PR dependencies from existing audited build, OMP1 and2 with BLAS1. Source fixtures are the established synthetic normalized orbitals on8³ grid with2/3k points, not converged physical GS. dt.001,12 total steps. Four fast-engine configurations: ACE Taylor, learned-ACE Taylor, full-action Taylor within fast PR, learned PT-CN.

At split6 and7, all16 cases per OpenMP setting (32 cases total) produce bit-identical charge current and total energy after6/5 remaining calls. Concatenated pre-save and resumed `rt_probe_write` traces exactly equal continuous traces. The trace covers accepted step,current,energies,both orbital buffers,density,Hartree histories,local/XC potentials,exchange energy,factor/PT-CN counters,ACE factors,cached source and factor-history samples. It is a selected internal trace, not a claim that every scratch field was compared. Future learned predictions are exercised (e.g. steps9/11 predicted after split7).

Twenty-four rejection cases per OpenMP setting (48 total): half-truncated archive, different dt/settings, extra trailing byte after HSE payload. All refuse with zero successful propagation calls. Late failure additionally yields accepted_step=-1 and rejects an attempted rt_step. Eight prior full-step two-context regressions (different dt/capacity/OpenMP route,solo/interleaved/reverse/close-reopen) still pass.

Checkpoint size426,917–708,649 bytes for these tiny synthetic inputs. This is not a production resource or memory-limit claim. Saving is at accepted steps with initialized exchange resources; step-zero save is not established. Old native text file positions are not resumed.

## Evidence

Persistent developer root: `~/.codex/outputs/salmonpy-state-inputs-20261116` (not an installation prerequisite).

- `hse-direct-verified` and `hse-direct-verified-omp2`: result.json with16 positive and24 rejection records each, compiled sources,build audit logs,hash receipts,traces/API/calls.
- `hse-direct-context-regression`: eight existing independent-state cases.
- `hse-input-build.json`: audited fast engine SHA256 d0c47cd74c4911cb18fd5f8eb668b65bfac83e3293ce9471f8e4e2aef1422728; same binary retained.
- `hse-direct-installed`: installed-wheel adapters/tool/probe repeat16 positive and24 rejection cases successfully, using installed relative paths.
- `fast-hse-checkpoint-final-dist` and the package/source/installed test logs: source and wheel450 tests,1 skipped,no failures; sdist-to-wheel built successfully. Runtime JP/EN guides and AI JSON synchronized.

Read-only reviews found no missing evolving state for the admitted profiles. Review also prompted the tested late-failure quarantine and precise description of selected trace coverage. No Actions used.

## Scope

Internal serial CPU fast-HSE Fortran context, fixed occupation/ions,periodic scalar full-support profile. No new stable Python/C handle or host-MPI/FDTD-HSE connection is claimed. General native same-input restrictions, broader independent GS input ownership, generic archives, other resources/topologies/platforms and physical convergence remain separate. Molecular work stays paused.
