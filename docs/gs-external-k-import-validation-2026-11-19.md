# Existing-GS external k-table import — 2026-11-19

Preceding regular-k/MPI export committed and pushed asd67ce847 to the salmonpy
remote branch; remote SHA verified. [skip ci], no manual Actions dispatch.

NativeGSSession.load omitted actual native input file_kw from its asset scanner
(it recognized file_kgrid instead). The retained work directory and subsequent
memory-GS→TDDFT model could therefore lack the referenced k table; native file
open could fail. Add file_kw while retaining the prior recognized names, copying
nested paths and verifying bytes on each subgroup rank. Comment masking now
preserves ! inside quoted filenames. No Fortran/numerical update was changed.

A mock-native regression reproduces missing copying with the preceding version;
fixed code preserves nested tables/k!.dat and rejects its absence before any
native call or new import directory. The existing shared no-SCF comparison
runner gains --external-k, avoiding another separate native test executable.
It stages working copies of existing GS records in scratch only and converts the
regular mesh to an identical explicit table; original GS inputs remain unchanged.

Existing GS/native tests, internal MPI2 (four total ranks) cover PZ2-k orbital/
spatial, Libxc-PBE8-k orbital, TB-mBJ8-k spatial, plus unchanged regular control.
All16-step currents equal the independent file route exactly. Independent RT
states and checkpoint-at8/resume8 agree exactly with zero replay. A missing k
table during active RT is rejected collectively before native initialization;
GS resource counts/current remain unchanged, then valid loading/RT continue.
The final installed-wheel native external-k/orbital2 comparison also passes.

Source624 and installed624 tests pass with one skip each. Wheel built offline
from sdist, source/all-application-object audits and native wheel hashes pass.
Existing electronic library was reused after recorded SHA256 verification,
not rebuilt; no new SCF/reference GS or real GS collection into results storage.
JP/EN/AI guides, shared runner packaging and checklist synchronized.

This fixes referenced-asset ownership, not canonical physical metadata binding:
imported calculation stays None, convergence stays unknown and forces unavailable.
Imported-owner portable export, full GS/native SCF mixing checkpoints and changed
MPI topology restoration remain separate. These existing-GS parity tests do not
certify physical convergence. Measurement host is Apple Silicon/GNU/OpenMPI.

Selected summaries/logs/hash inventory/wheel only:
/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/salmonpy-gs-k-table-import-20261119
No GS restart/orbital/density arrays collected here. Local writes verified;
cloud synchronization unverified. Continued changes remain uncommitted.
