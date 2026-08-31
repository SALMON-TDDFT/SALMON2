# Hybrid Variational Payload Replay Design

## Purpose

The eight-rank Si64 Hybrid continuation reaches the LCFO symmetry handoff but
fails before the first continuation iteration because the fixed variational
payload has an invalid extent.  Conventional SALMON GS restart cannot shorten
the diagnosis: DC mode rejects `yn_restart='y'`, and the failed run did not
serialize the converged DC density and fragment orbitals.

Capture the already-assembled Hybrid payload at the failure boundary once, then
replay its validation without repeating DC or Wannier construction.

## Scope

The diagnostic bundle contains only the inputs to
`freeze_dg_hybrid_variational_payload`:

- global basis count;
- local row IDs;
- metric, kinetic, nonlocal, and interface matrix rows;
- basis, metric, and interface fingerprints;
- MPI rank count and local rank identity.

It does not attempt to restart DC, Wannier90, Hybrid continuation, or RT.  The
bundle is diagnostic evidence and remains untracked under the existing Si64
verification directory.

## Data flow

An opt-in environment variable names a bundle directory.  Immediately before
the fixed-payload freeze, every rank writes one versioned unformatted stream
file.  Rank zero writes a small manifest only after all rank files have been
published and an MPI barrier succeeds.  Existing production behavior is
unchanged when the variable is unset.

A standalone MPI runner reads the manifest and the file matching each current
rank.  It validates the version and rank count, reconstructs the arrays, prints
the exact observed dimensions and row-ID range, and calls the production freeze
routine.  It must reproduce the original acceptance or rejection collectively.

## Safety and error handling

- Never overwrite an existing bundle directory or rank file.
- Use temporary rank files followed by atomic rename.
- Reject truncated data, version mismatch, rank-count mismatch, invalid array
  allocation sizes, or inconsistent global metadata.
- Preserve the original Si64 logs and all captured files on success or failure.
- Do not stage generated payload bundles.

## Testing

TDD covers a small valid distributed payload and fixtures with a row-ID range
error and a matrix-shape error.  The writer/reader round trip must preserve all
values and fingerprints on 1, 2, and 4 MPI ranks.  A parser-only test confirms
that the replay diagnostic identifies the individual extent condition.

After focused tests and a SALMON build pass, run Si64 once with capture enabled.
All subsequent diagnosis uses the captured eight-rank replay bundle.  The
actual source defect is fixed only after replay supplies a minimal failing
case, followed by the original protected runners and a final fresh Si64 run.
