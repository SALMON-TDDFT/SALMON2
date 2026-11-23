# Direct mixed-session checkpoint implementation inventory

User requests direct field/LD/electronic save and fresh-job restore without replay. Replay implementation is retained independently at09045027. Direct support is now implemented for the same-topology periodic PZ native mixed fixture; see docs/native-direct-checkpoint-validation-2026-11-16.md for actual acceptance results. Other profiles and generic user-backend archives remain open.

## Located persistence boundaries

- `runtime/sources/standard/src/maxwell/embedded_eh.f90`: field_context owns fs,fe,ld,points,current,step,capacity. Constructor currently rejects yn_restart other than n. A/PML histories are not captured by Python A/E.
- `src/maxwell/fdtd_eh.f90`: eh_calc already writes native restart records starting at the flag_o_restart block, including split E/H, LD and observation data; eh_init reads them for yn_restart=y. Extract reusable save/load routines instead of invoking an extra physics step. Retain halo/bounds and native order.
- `src/rt/embedded_tddft.f90`: rt_context includes both wavefunction buffers, Hartree history, RT state, XC and worker current; accepted_A lives separately. Preserve accepted_A and buffer parity explicitly.
- `src/io/checkpoint_restart.f90`: checkpoint_rt/restart_rt exist, but their arguments and format alone do not establish complete embedded context restoration. Audit density/local potentials, iteration/parity, pp k-vector, XC/TB-mBJ/fast-HSE derived resources before advertising support.
- `salmonpy/coupled.py`: NativeEHField A,E,step_index and CoupledStepper J,step_index must be restored together. Restoring only native arrays loses accumulated vector potential and lag1 feedback.
- ResponseSession point states, previous drive and step index must share one checkpoint epoch. SBE uses its existing direct snapshot; native TDDFT needs an explicit collective hook.

## Required delivery order

1. Extract collective native E/H save/read callable without eh_calc; verify E/H/PML/LD-only continuation first, with nonzero LD memory and MPI halos.
2. Add embedded TDDFT direct save/read, restoring current, A, buffer parity and Hartree history. Start with supported PZ fixture; refuse unverified XC paths explicitly.
3. Add backend capability hooks and session coordinator; reject unsupported user backends before writing or changing live states.
4. Write new temporary archive with version, completed step, native library/input hashes, material/point IDs, rank topology and per-file hashes; collectively agree success before publishing completion manifest. Same topology initially.
5. Fresh-job initialization may rebuild static coefficients from original inputs, but must perform zero propagation steps before loading dynamic arrays. Restore native and Python clocks/current in a synchronized epoch.
6. Acceptance: field1+SBE1+TDDFT2 with local two-material LD, checkpoint64/continue128 vs continuous128; instrument resumed propagation count to equal64 (not128). Compare A/E/j and relevant internal arrays; test corrupted/incomplete archives and rank/identity mismatches before mutation. Separate field2 check under four-core limit. No routine Actions.

Original sources/manifests/provenance must be updated through the repository source workflow; never edit previously prepared external trees or fixed binaries. Build a new external output directory, preserve shared-kernel byte equality, and run packaging checks for layout changes. Direct restore remains unchecked in the active checklist until the acceptance tests actually pass.
