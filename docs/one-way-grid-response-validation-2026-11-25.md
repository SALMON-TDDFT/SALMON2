# One-way spatial molecular excitation (2026-11-25)

Scope: saved LD + 3D FDTD component E/H -> spinless fixed-ion PZ C2H2.
There is no current return, own-vacuum subtraction, radiation or rescattering calculation.
Python selects frames and invokes the existing Fortran spatial interpolation/path coupling.
The numerical worker and electromagnetic solver were unchanged.

Validated worker SHA256: `14557305ee88f46a99d3e15352f11ffb88c7c1c80626b93b21b2659a1b035897`.
Input fields came from the existing two-field-rank illumination-only run
`magnetic-feedback-ld-external/history.npz`, array `total` (not a feedback trajectory).
Component origins came from its matching prepared grid. Molecule position:
(-2.54,1.62,1.94)bohr, off the field grid; patch24^3, spacing2bohr, dt0.003au.

The public `molecular_grid_response` example completed96molecular steps using
frames0..97. Checkpoint at48, process exit, then restore and continue48steps:
all continuous output arrays (step numbers, time, dipole, electronic change,
electron number) matched bit-for-bit after concatenating the two segments once.
Restoration updates:0. Peak electronic dipole change:0.0007455745791570146au.
The initial internal energy is unavailable; the example records availability and
omits the energy array when it is not available at every output frame.

Canonical local record: `calc2026/salmonpy-c2h2-nearfield-20261125/one-way-ld-fdtd-c2h2/comparison.json`.
An additional four-step manufactured spatial E/H history also passed process-exit
restart at2. Failed intermediate sample output checks are retained in the v2 directory.
No new GS or FDTD run was performed. This validates data transfer and restart,
not physical convergence, realistic nanostructure spectroscopy or strong-back-action accuracy.

Distribution checks: source sdist -> macOS arm64 CPython3.14 wheel succeeded,
including source and application-object eigensolver prohibition audits. The
build used existing Make because Ninja was absent; no source audit was skipped.
The installed wheel passed32related tests; source passed the same32tests.
From outside the checkout, the installed sample restored the96-step real-field
case at48 with0restore evolutions, executed48new steps and matched all source
resumed arrays bit-for-bit. See `installed-comparison.json` beside the local record.
Guide synchronization and whitespace checks passed. No CI run was requested.
