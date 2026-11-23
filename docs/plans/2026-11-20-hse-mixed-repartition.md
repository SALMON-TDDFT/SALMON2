# Fast-PR HSE mixed SCF rank-count restart

Extend the measured dense Gamma profile to simultaneous orbital and Cartesian
spatial ownership. Keep k groups one; multi-k spatial is a separate backlog item.
Use existing hse_spatial selector with compute protocol2 and mixed migration2.
Map accepted wave/CG data in five-dimensional ownership boxes. Map HSE cached
source/action, ACE, workspace and previous/source matrices by global grid rows
and global orbital columns. Replicated gauge/occupation and scalar mixing use
one donor. Rebind new communicators only after complete collective staging; no
SCF, density or exchange rebuild on restore. Retain pure-layout arithmetic.

Verify 2↔4 and changed-layout4→4, uneven orbital blocks, pre-PBE phase, source
and installed-wheel continuation vs source-initialized continuous8 at unchanged
rtol1e-9/atol1e-10, exact cut interiors and zero bootstrap. Late missing donor
keeps original state/files. Verify independent states, old pure-profile parity,
source/object audits, existing Python suite, docs/JSON/sample parity. No CI.

Accepted: source2↔4,4→8,changed-axis4→4,uneven6→4,pre-PBE2↔4 and
installed-wheel4→2 passed exact physical wave/density cut interiors, zero SCF
bootstrap and continuation within unchanged1e-9/1e-10. Late donor failure preserves
original arrays and unchanged source hashes. Independent different-input mixed
owners and protected late-failure/continue/retry exact; previous pure orbital/
spatial HSE and ordinary-XC results exact. Public wheel sample4→2 and source/
installed637 tests each636pass/1skip passed. Source/object audits and source/
wheel/docs/native-byte parity checked. User explicitly approved this source/wheel
HSE mixed GS/restart dataset campaign. Evidence:
calc2026/salmonpy-hse-mixed-record-20261120. No CI dispatched.
