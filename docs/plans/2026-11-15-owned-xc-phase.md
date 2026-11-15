# Explicit XC evaluation phase

Continue input/state separation after commit 65cc34ce. Native TDDFT and explicitly owned HSE evaluations must not inherit another calculation's PBE pre-SCF flag.

1. Reproduce optimized native HSE action changing when legacy exx_pre_scf_active is poisoned after owner creation.
2. Add per-system phase ownership plus explicit internal set/query entries. Route native HSE enablement/thermal damping and XC semilocal/rVV10 evaluation through the same phase query.
3. Native RT initialization explicitly owns production phase before potential evaluation. Legacy GS/DC drivers retain their existing global-phase compatibility path until their SCF state machines are extracted; do not silently freeze their phase transitions.
4. Test independently selectable phases, unchanged unrelated state, explicit transition back with refresh, legacy fallback, full RT interleave and original-driver preservation. Update English/Japanese/AI documentation and package locally. No Actions.

This does not remove legacy SCF convergence/phase globals from the GS/DC driver. It establishes explicit phase ownership for RT and internal evaluation clients. Switching from pre-SCF to production requires rebuilding exchange from the current orbitals before using its action. Other algorithm/global-mode restrictions remain.

Review correction: a phase transition clears the explicitly bound owner (including frozen prediction flags) before future refresh; otherwise refresh could silently retain a stale operator. Test transition with changed orbitals and a frozen owner.
