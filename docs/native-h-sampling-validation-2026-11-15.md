# Native H accessor validation — 2026-11-15

Added optional sp_eh_sample_h to bundled standard embedded_eh. No eh_calc changes:
raw split H at (n+1/2)dt is averaged over four transverse Yee faces to each macro
point. Each rank contributes owned faces only, followed by global field-group
sum. Python sample_magnetic_field is explicitly collective, not a local property.
Old ABI rejection and native-liveness preflight happen before native collectives.
Source inventory/hash and RT-step overlay base hash were synchronized.

Apple Silicon Mac, GNU15/OpenMPI, O2, OpenMP/BLAS one thread, native build jobs2:
- Build explicitly selected OMPI_CC=gcc-15: default mpicc chose Clang, which did
  not accept -fopenmp. No global toolchain configuration was changed.
- Native tests at one/two field ranks (2/3 total ranks) compare mean of consecutive
  raw H values to existing saved H at n dt. Max error 4.632211430296955e-23.
- Initial comparison test mistakenly averaged the already spatially averaged
  observation output again; fixed the reference, not native sampling arithmetic.
- A/E/H arrays match exactly across one/two field ranks and observation strides
  1,2,64. Existing FDTD wrapper requires observations and stride<=steps; no claim
  is made for an observation-free wrapper configuration. H reads on non-output
  steps are covered. Old library's E history matches new library exactly.
- Repeated reads do not change A/E or step count.
- Installed-wheel three-rank native test includes one-rank liveness failure,
  verifies collective refusal, restores valid context and continues successfully.
- Source suite before final liveness regression: 419 tests, one skip. Final
  focused coupling tests: 13 passed. Installed full suite: 420, two skips.
- sdist→wheel build/install successful. Review's collective-preflight issue fixed
  and exercised with real MPI. No Actions or physical convergence study.

Source/build: /tmp/salmonpy-h-sample-source-20261115,
/tmp/salmonpy-h-sample-build-b-20261115/build.json.
Native outputs: /tmp/salmonpy-native-h-fixed-result-{2,3},
/tmp/salmonpy-h-stride-valid-result-{2,64}, /tmp/salmonpy-h-installed-native-result.
Wheel: /tmp/salmonpy-native-h-reviewed-dist-20261115/.
Tests: tools/test_native_h_samples.py BUILD_JSON NEW_OUTPUT [STRIDE].

The existing coupled TDDFT example adds opt-in --sample-h, called on every field
rank and saved with H_time separately. That full TDDFT example was not rerun;
actual H sampling above used native FDTD. Packet feedback integration, magnetic
material physics and prebuilt-runtime publication remain separate work.
