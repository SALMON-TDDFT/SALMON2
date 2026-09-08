# Si8 MPI8 Memory Diagnostic Design

## Goal

Exercise the complete overlapping-Wannier ground-state route on the eight-atom
silicon conventional cell with eight MPI ranks while measuring whether resident
memory grows after Wannier90.  This is a diagnostic precursor to the Si64/C64
acceptance run and must not use periodic process suspension.

## Physical and parallel model

Use the established silicon conventional cell (`al=10.26` a.u., eight Si atoms,
32 electrons) and the tracked Si pseudopotential.  Use a `2 x 2 x 2` fragment
decomposition on eight MPI ranks, so every rank owns one physical fragment.  A
`16 x 16 x 16` real-space grid keeps every fragment extent integral.  The
overlapping-Wannier target is derived by production code from the occupied
space plus the complete pseudopotential-backed `s+p` complement; no
material-specific target rank is forced.

## Runtime safety

Run with one OpenMP and one BLAS thread per rank.  Use `nice` and background
quality of service only; never send periodic `SIGSTOP`/`SIGCONT`.  Set
`GFORTRAN_UNBUFFERED_ALL=1` so phase diagnostics are durable.  A separate
monitor records timestamp, phase marker, available memory, swap use, and RSS
for every SALMON rank at 30-second intervals.  The monitor terminates the whole
MPI job once if available memory drops below 8 GiB or if an explicit configurable
per-rank RSS ceiling is crossed.  It never pauses and resumes ranks.

## Evidence and acceptance

The run directory is persistent and records hashes of the executable, input,
atomic coordinates, and pseudopotential.  Acceptance requires the expected
8-atom/8-fragment route, successful DC-SCF and LCFO ranks, completed Wannier90,
translation-character and point-cogroup receipts, and V3 publication.  Memory
evidence must show the post-Wannier90 peak and whether RSS is bounded or grows
monotonically.  A failure is still useful when the last flushed phase marker
and memory series identify the failing stage.

## Scope

This diagnostic does not replace the full Si64/C64 acceptance case.  It tests
the same MPI rank count and principal algorithmic route with substantially
smaller state and grid extents.  No response/HHG run is started until this GS
memory diagnostic is understood.
