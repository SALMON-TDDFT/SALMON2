# Fast-PR HSE orbital-parallel SCF direct restart

User-authorized next step after ordinary orbital repartition. First expose the
existing Gamma orbital-distributed fast-PR exchange for retained scalar fixed-ion
GS, without changing its arithmetic. Move spatial exchange numerical workspace
from module SAVE into the owned HSE state; copy/serialize it with ACE, wavefunction,
CG, density and mixing history. Do not serialize CPU procedure pointers/backend
plans; only a CPU full-support Gamma profile is admitted. Preserve prior k-MPI
and CLI policies and numerical results.

First acceptance: same-layout orbital MPI2 fresh-process restore, including a
cut in PBE pre-stage, no SCF bootstrap/replay, continuous result and independent
state equality. Existing ordinary/k-HSE regression and installed wheel required.
Changing HSE orbital rank count needs additional ACE/workspace reindexing and is
separate from this same-layout contract. Update paired guides, AI JSON and queue.

## Verification

MPI2 same-layout HSE/PBE-stage direct restart passed: every serialized numerical
field exact before continuation, zero SCF replay, five new iterations exactly
match continuous8. Two different-input states, independent MPI groups, protected
late-failure/continue/retry passed. Installed direct restart and sample passed.
Source/installed637 tests each636pass/1skip. Ordinary and k-HSE old-binary SCF
and16-step RT parity exact; no old-binary orbital equivalence is claimed.
Initial unsupported-status15 evidence is retained: SALMON auto-selects its Gamma
orbital exchange via yn_hse_wannier; the explicit owner admission now recognizes
that existing full-support route. Numerical formulas were not changed.
Evidence: calc2026/salmonpy-hse-orbital-record-20261120. No CI dispatch.
