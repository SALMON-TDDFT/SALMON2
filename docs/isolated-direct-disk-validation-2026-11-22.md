# Molecular direct disk restart — 2026-11-22

## Result and scope

The public isolated PZ fixed-ion state now has save_direct_checkpoint and
restore_direct_checkpoint. Fortran writes and reads both orbital buffers,
temporary orbital, density/potential and RT histories, Poisson work, energy,
XC payload, field vectors and accepted index. Array codecs are shared with the
periodic context; existing codec statements are unchanged except for replacing
the caller-owned limit reference by a caller-owned accessor. No numerical Python
fallback or additional evolution is used.

A new-directory manifest (salmonpy.isolated-direct/1) binds the executable,
input/settings, output streams and numerical state with per-file size/SHA256.
Import is into an untouched candidate process. Complete header, bounded array
reads and trailer are required before installing numerical state; output streams
are reconnected before switching the Python owner. Import failure preserves the
current process and observations. Imported states never silently replay on failure.
Explicit history replay archives are separate and are not a fallback.

## Measured evidence

Authorized existing C2H2 GS reused; no SCF. PZ, fixed ions, dt=.003 atomic units,
capacity20, three supplied x/y/z field vectors. Original process exits before
resume at steps0,7,8,20. A different worker imports with zero actual evolution
calls, then uses precisely20-saved_step calls. Dipole change, total dipole,
electron number, energy/freshness and flushed output streams match uninterrupted
execution exactly. Both the low-level worker proof and final installed-wheel
public API proof passed. A corrupted completion trailer with an updated file
checksum is rejected by Fortran in a candidate; the original PID, step, update
counter and dipole are preserved.

The installed response example ran as separate invocations: save at step7, exit,
resume remaining input samples. Concatenated time, E, dipole and electronic-change
arrays match the continuous example exactly. The comparison harness initially
used the dataclass constructor instead of OutputArchive.open; this analysis call
was corrected and the existing successful calculation outputs were compared.
No calculation rerun was needed for that harness correction.

Source and installed wheel suites each ran691 tests, one optional skip. Packaged
code, paired guides, AI JSON and sample bytes match source. Synthetic checks cover
process exit, explicit read limit, corruption before spawn, model/position mismatch,
non-overwrite and candidate failure preserving the original state.

Tracked-only source export preparation/build and source/application-object
prohibition audits passed. Shared codec byte-expression comparison passed.
NLCC prescribed-density XC and force tests passed at MPI1/2. C/Fortran ABI query
checks passed; these query checks do not claim MPI initialization or state stepping.
The regular E-worker comparison retained its previous maximum errors:
dipole3.85523175341039e-12 e*bohr; energy2.1316282072803006e-14 Hartree.

## Artifacts and limits

Canonical local OneDrive campaign: calc2026/salmonpy-c2h2-feedback-20261122.
Use direct-disk-verified-summary.json, direct-disk-verified-artifacts,
direct-disk-git-worker-validation, direct-disk-verified-wheel-public,
direct-disk-verified-wheel-complete and direct-disk-verified-wheel-example.
Preliminary development runs are not the final executable validation receipt.
The final receipt is revision5 and contains real process-exit proof.

Source inventory SHA256:
c85b972a6ede6baebd240532893b4a57cb3cc434f8e6b4651fcab5577c90e137

Runtime executable SHA256:
473f3f8f91cc44e61995e55457e484b6fb8c1120f695e39ce28c83cd3bb750e8

Wheel SHA256:
315691f1957dbcfd69cef74ab436334f5e926653582d66768592dfeeb4b40459

Same runtime, GS/input, timestep, capacity, position and threads are required.
The original GS is loaded for initialization but elapsed propagation is not
recomputed. A temporary additional process/state is needed during import.
This is molecular-state persistence, not coordinated durable restart of both
FDTD fields and their accepted field/electronic step indices. No isolated inner-MPI, cross-compiler
portability, arbitrary OS/MPI failure or converged-optics claim. Local storage
does not prove cloud synchronization. No extra physical CI was dispatched.

The final wheel removes the development-only fixed manifest size cap; only an
explicit caller max_bytes budget can reject a larger metadata file. The full
source/wheel suites and installed real-state restart/example checks were repeated.

A final raw-worker check also compared the complete set of .data files after
both processes had closed (including end-of-run output). All bytes match at
saved steps0/7/8/20; direct-disk-postclose-proof/report.json records this check.
