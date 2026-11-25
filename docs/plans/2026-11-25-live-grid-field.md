# Live grid-field input and electric feedback

Continue the approved next stage: explicitly stage two E/H grid frames for native indices n,n+1, then advance molecular step n using the unchanged integrator. Validate target step, fixed geometry, dt and all samples before installing either frame. A failed stage must preserve state. Reuse direct savepoints/checkpoints; never reconstruct spatial input by E-only replay. The public API requires executable-bound live-grid evidence.

First verify streamed input against the same prescribed field history, including orbital magnetic input and process-exit restoration. Then use total and own-vacuum FDTD groups for an electric-only density-current feedback prototype: field steps n,n+1 drive molecular step n; start with a declared field hold. Trial zero-current field updates predict the drive; restore and inject the density-derived current; reject/retry if the accepted selected field differs. Keep step numbers explicit, domain ownership and no application-wide broadcasts of full grids. Full magnetic backreaction and physical convergence remain separate work.

### 2026-11-25: installed streamed-field and electric-feedback verification

- Git-only solver build and source-distribution-to-wheel passed; installed targeted129tests passed. Public isolated_response --grid-series matches the accepted-step reference.
- Final executable and installed-wheel6-rank LD feedback: direct restore molecular96/FDTD97, retry128 and192-step continuous control agree; six ordinary output streams are byte-identical, with zero extra restoration evolution.
- See docs/grid-field-feedback-validation-2026-11-25.md. Magnetic-current return, general production coordinator, rolling storage and physical convergence remain open; this test does not close those items.
