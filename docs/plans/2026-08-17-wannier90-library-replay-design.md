# Wannier90 Library Replay Design

## Purpose

The Si64 diagnosis currently spends roughly twenty minutes constructing the
DC/LCFO state space and the Wannier90 library inputs before the site-symmetry
iteration starts.  Export the already assembled standard Wannier90 inputs so
that convergence experiments can be replayed with `wannier90.x` without
repeating the SALMON preparation.

The normal in-process library calculation remains unchanged.  Export is an
explicit diagnostic side effect and is disabled by default.

## Interface

Set `SALMON_DG_W90_REPLAY_DIRECTORY` to a non-empty directory before running
SALMON.  Immediately before `run_dg_w90_gamma_library`, rank zero writes the
standard files for seed `overlapping_wannier_mlwf`:

- `overlapping_wannier_mlwf.win`, produced by the existing setup path;
- `overlapping_wannier_mlwf.dmn`, produced by the existing symmetry path;
- `overlapping_wannier_mlwf.mmn`, from the assembled neighbor overlaps;
- `overlapping_wannier_mlwf.amn`, from the assembled projections;
- `overlapping_wannier_mlwf.eig`, from the assembled eigenvalues.

The requested directory must already exist.  This avoids hidden directory
creation and shell-dependent behavior inside the MPI program.  If export was
requested, invalid data or any I/O failure is a collective fatal error rather
than a silently incomplete replay bundle.

## Data and Ordering

Reuse `write_sawf_local_eig_amn_mmn` rather than introducing another file
formatter.  Its conventions match the Gamma-only library payload:

- `.eig`: `(band, kpoint=1, eigenvalue)`;
- `.amn`: `(band, projection, kpoint=1, real, imaginary)`;
- `.mmn`: one block per neighbor, followed by `(band_i, band_j)` values;
- neighbor reciprocal lattice vectors come from `w90_nncell`.

The writer receives the root-owned `w90_m_matrix`, `w90_a_matrix`, and
`w90_eigenvalues` directly.  It streams formatted records and does not make a
second copy of any dense matrix.  Non-root ranks only participate in agreeing
the option and export status.

The `.win` and `.dmn` files used for replay must be the same files prepared for
the library call.  The export hook therefore copies those small text files to
the replay directory when their source directory differs.  Copying is done by
Fortran stream I/O with checked status, not by invoking a shell command.

## Replay Driver

Add a small Python driver that accepts a replay directory and optional edits
to `symmetrize_eps` and the patched site-symmetry iteration limit.  It copies
the five standard inputs into a temporary run directory, edits only the copied
`.win`, invokes the selected standalone `wannier90.x`, and retains the `.wout`
and exit status.  The original export is immutable, so several convergence
settings can be compared safely and quickly.

The iteration limit remains a build-time Wannier90 patch for now.  The driver
reports the limit encoded by the executable/build metadata and never claims
that editing `.win` changes it.

## Failure and MPI Semantics

Before any rank branches on export, all ranks agree on whether the environment
variable is present and on its exact value.  Only rank zero performs file I/O.
Its result and diagnostic message are broadcast before any rank proceeds into
the library call.  Thus a missing directory, partial write, or inconsistent
environment cannot leave ranks on different collective paths.

The existing writer validates dimensions and finite values.  The integration
also checks that `w90_nncell` has one vector for every overlap neighbor and
that all generated files are present and non-empty before reporting success.

## Verification

1. Extend the focused writer test with the exact Gamma-only library shapes and
   verify `.eig/.amn/.mmn` headers, neighbor vectors, and complex ordering.
2. Add an export/replay fixture that supplies a `.win/.dmn`, runs the export
   hook, and verifies all five files without entering a long Wannier solve.
3. Add route checks proving export occurs after matrix assembly and before
   `run_dg_w90_gamma_library`, and that the disabled path performs no export.
4. Add a standalone replay smoke test against the bundled `wannier90.x` when
   available; otherwise report a clean skip.
5. Run the focused MPI fixture on 1, 2, 4, and 8 ranks and compare exported
   payload hashes.

## Non-goals

- No change to the numerical library call or production defaults.
- No checkpointing of the earlier DC/LCFO computation.
- No new dense matrix allocation or rank-wide replication.
- No automatic retry or convergence-policy change inside SALMON.
