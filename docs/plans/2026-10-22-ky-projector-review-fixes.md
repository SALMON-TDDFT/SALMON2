# KY projector review fixes

Human request: repair all three reviewed defects, verify Coding Rules, commit and push the existing branch. Existing calculations and binaries remain untouched.

1. Add numbered CPU regression 453, exercising the actual compiled input_pp/salmon_pp path. Reproduce absent scalar-KY minus-j projectors, unequal masked equal-j radial tables, and double-expanded mask radius before editing production code.
2. For legacy scalar KY under SOC, copy the completed scalar projector tables and normalization into the minus-j branch. Do not infer that an extended j-resolved file is scalar from zero-valued channels.
3. Apply the existing regular-solid-harmonic radial conversion to the masked minus-j tables and derivatives. Set the masked support once per species before filtering either branch; retain the old single-pass scalar support rule.
4. Verify source prohibition before builds and all linked application objects before probe execution. Run bounds-checked GNU serial and MPI2 fixtures plus existing 449–452 regressions. Check extended W/Se reader and mixed-species full-application smoke. No full GS/5fs resumption.
5. Independent code review, narrow diff/Coding Rules audit, then commit and push SALMON-TDDFT/SALMON2 codex/fix-cartesian-soc-symmetry. Record actual verification and optional configurations not tested.

## Verified outcome

- RED: production-object fixture reproduced missing scalar minus-j projector signs, unequal masked equal-j tables and double-expanded support. The final fixture treats local zero-norm channels as valid and permits floating-point rounding at a grid boundary; numerical agreement tolerances are unchanged.
- GNU serial Debug and MPI Release, both with `-fcheck=all -fbacktrace`: source and application-object prohibition guards PASS. Focused CTest 449–453: 13/13 PASS in each build, including five real SALMON reader smoke configurations and MPI2 rank-table identity.
- A first sandboxed MPI invocation could not open local communication sockets (exit213), before SALMON began. The same suite with local sockets permitted passed; no numerical threshold change.
- Independent read-only code review found no critical/important issue. `implicit none`, `use, only`, explicit interfaces, input namelists, MPI broadcasts and free-form line limits retained. No new input/environment switch or mutable shared state.
- New regression453 is GNU/Makefiles CPU-only because it links the real application's objects using its existing link command; all linked application object directories are inspected. Ninja, accelerators, HSE/LibXC and ScaLAPACK/EigenExa configurations were not exercised for this change.
- No GS/5fs physical-accuracy acceptance or production restart. Existing calculation inputs, pseudopotential data and running/fixed production binaries were not changed.
