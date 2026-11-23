# Native fast-PR geometry admission

Extend pre-native host-MPI HSE admission to the effective periodic cell, grid and
explicit atomic source. Read-only, NumPy-only, collective rejection before work
copy/native lifetime. Retain SALMON cell priority (al, diagonal vectors, n*dl;
nonorthogonal-vector detection) and actual int((length+1e-10)/spacing) resolution; never
ceil, normalize, change the native input or certify restart provenance.

A focused numeric namelist reader handles only relevant scalar/vector settings,
indexed/ranged assignments and finite Fortran literals/repeats. Unsupported/null
forms fail clearly rather than guessing; comments/quoted paths stay intact.
Support canonical inline/file Cartesian or reduced atom records with exact natom/species
references. Keep the existing fast-PR au-only input contract.

Symmetry checks: effective Cartesian cell isometry, retained grid mapping via its
affine index generators (no grid-sized allocation), unique same-species atom
permutations modulo cell using native tolerances. Validate finite/positive actual
spacing/counts and native integer bounds. Do not estimate binary restart layouts
or imply k-star, pseudopotential-content or general initialization recovery.

Test first: invalid/singular cell, anisotropic forbidden rotation, invalid grid
or translation, bad atom count/species/coordinates and symmetry mismatch. MPI4
(two groups2) corrupt one rank, reject without lifetime, repair/retry; normal
current/energy and own direct resume remain unchanged. Update JP/EN/AI docs,
checklist; source and installed suites + sdist-to-wheel, no Actions.


Measured scope correction: equivalent axis-permuted left/right-handed fixtures
hit ACE initialization STOP in the existing native library. Decline off-axis
vector cells before native lifetime; do not silently rotate/reorder a basis or
claim the exchange-core geometry tests prove TDDFT initialization support.
A_eV_fs is already rejected by FastHSEMaterialProvider; retain that contract.
