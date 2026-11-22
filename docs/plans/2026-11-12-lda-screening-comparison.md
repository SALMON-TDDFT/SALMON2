# LDA screening at a common TB-mBJ RT reference

Approved scope: compare LDA and TB-mBJ screening on Mac at small k, keeping the TB-mBJ propagation reference. Start with Nr=12^3, Nk=Nq=2^3, boxes=6^3, cutoff=3 Ha; GS has four occupied and four empty orbitals and RT propagates four occupied orbitals only. This is a mechanism test, not k convergence.

1. Observe rejection of scalar LDA internal screening in the unchanged executable.
2. Permit unpolarized PZ LDA GS internal screening. Store its generating potential in a distinct HXL1 file, retaining HXR1 KLI and HXT1 TB formats unchanged.
3. Add explicit `hsex_screening_w` initial kernel override for TB-mBJ HSEX RT. First validate the normal TB reference file against the TB GS, then independently validate the HXL1 source geometry, k mesh, box layout and finite metadata/kernel. Keep the TB reference potential and initial density; do not rewrite either file's header. For this first implementation reject simultaneous pump W updates/replacement.
4. Test same-k baseline, LDA W and zero-field controls at the same dt, duration, impulse and cubic window, plus metadata/nonfinite rejection and restart behavior. Verify default path remains unchanged.
5. Run mandatory source/object audits, GNU MPI and serial builds, existing current tests and the Mac calculation. Record source/binary/pseudo hashes. Save only RT/W/logs/analysis records in canonical calc2026, never GS restart data.

The existing W format has no pseudopotential fingerprint; this controlled fixture verifies identical pseudo bytes in the run manifest. General cross-source portability is not inferred. Keep production jobs and their fixed executables unchanged.
