# Real-space HSE RT implementation plan

Goal: propagate mesh wavefunctions with SALMON Taylor4, without a fixed LCFO projector.
User approved replacement of the fixed-basis RT design. DC/LCFO supplies initial orbitals only.

Architecture: new mesh-row-distributed HSE adapter. Occupied-space MLWF rotation affects exchange sources only; ACE factors have mesh rows. Full kinetic/local/nonlocal action remains on the mesh. Global periodic screened kernel establishes the full-support reference; masked MLWF sources remain an explicit approximation. FFT work distributed by source columns, with streamed mesh target columns. This is a correctness-first implementation, not a claim of restored weak scaling; full-domain FFT/column collectives remain a scaling limitation.

1. Implement typed distributed grid exchange kernel using existing FFTW Wannier kernel, no basis projection. Validate against independent direct action, Hermiticity, arbitrary targets outside occupied subspace, MPI layouts.
2. Implement mesh MLWF initialization, polar U transport, source-only sphere mask, initial 99.9% coverage diagnostics and explicit accepted/predictor frame rollback.
3. Add namelist-only realspace RT controls. Route HSE refresh/action/Taylor stages to mesh adapter. Disable fixed-LCFO RT flags with actionable error. Existing DC import must not configure retained basis.
4. Build ACE from mesh-row C and W; apply via distributed overlap, never project local Hamiltonian. Impulse step1 must rebuild; cadence independent of U transport. Update cached energy consistently.
5. Integration tests: DC import, spatial/orbital MPI, all-grid propagation outside initial LCFO span, full/localized source, ACE/U cadence, reference parity. Sequential small jobs only.
6. Replace active 3D RT inputs and notes; mark all old fixed-basis speed/spectra as historical. Provide portable verified patch, compile/test limits, commit/push.
