# Source tracking repair after external review

The current manifest lists739 files.25 correct local assets were not tracked:
16 cmakefiles and9 version headers/pseudopotential/sample inputs. The root ignore
rules and the bundled nested .gitignore both hide them on this case-insensitive
checkout. Root negations alone cannot override deeper rules. No source manifest
or vendored bytes were changed to repair this omission. Only hash-verified
manifest paths were force-added; build/CMakeFiles remains ignored.

The existing inventory test now checks Git tracking in developer checkouts.
It reproduced the25 missing entries, then passed all8 source tests after staging.
A Git-index-only export (tree recorded in the validation receipt) prepared and
built USE_HSE=OFF from scratch with mandatory source/object/link audits. Its new
2-rank FDTD A/E/H arrays and .data outputs exactly equal the validated R2 baseline.
A fresh sdist/wheel made from that export passed641 tests (640pass/1skip) with
matching newly built LD/SBE/convergence kernels. No Python numerical fallback,
implicit build, new GS collection or new CI run was used. HSE ON was not newly
validated in this repair. Heavy CI avoidance is not a substitute for clean-export
validation after a source-layout change.

This packaging-only repair is independent of the pending native-input scope
changes in the working tree. Evidence at canonical results root:
calc2026/salmonpy-review-source-integrity-20261120. Local storage verified;
OneDrive cloud synchronization unconfirmed.
