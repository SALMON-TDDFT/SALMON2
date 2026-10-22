# KY SOC+NLCC integration and Coding Rules follow-up

Human approved the audit follow-up. Keep running calculations and their binaries untouched.

1. Integrate the candidate KY header dispatcher, j-resolved reader and per-species scratch reset into the existing input_pp_sub module. Preserve legacy KY; reject extended-format SOCoff explicitly. No unvalidated scalar averaging.
2. Remove only new generator COMMON blocks crps_search/crps_jdata via explicit inputs/outputs to legacy Pseudo_Potential. Provide an explicit interface for that call. Convert the new writer to a default-private module with one public entry point. Keep unrelated legacy routines unchanged.
3. Add numbered CPU regression450 with inputfile/preparation/verification/CMake registration: actual production reader extraction, legacy and extended data, local-channel zero norm, NLCC units and malformed/NaN/SOCoff refusal; writer polynomial origin regression; generator W/Se reproduction. Add a bounded real SALMON serial/MPI reader smoke after full build.
4. Check source PZHEEV prohibition before either serial/MPI builds and all application objects before running rebuilt executables. Run registered tests and standalone generator reproduction. Independent review before commit. No existing binary/job/source snapshots overwritten.
5. Update instructions/manifests for writer module compilation order and reader now integrated. Commit/push this focused follow-up to the previously approved same repository/branch after checks pass. Physical 5fs mesh acceptance remains separate.

## Verified outcome

- GNU serial Debug with bounds checks and MPI Release builds passed source/object prohibition guards. CTest450/451: 6/6 passed in each build; full smoke includes legacy Si, W/Se, W/Si both species orders, and Se (one SCF iteration only).
- Existing449 covariance tests: 5 passed. Generator: W/Se byte-identical with Homebrew GCC15.2; portable CTest uses finite numeric comparison (relative/absolute1e-8) with exact discrete format metadata.
- Actual bounds-checked integration exposed a pre-existing calc_uv_so indexing inconsistency: set_nlma_so/set_lma_tbl exclude channels with both projector norms zero, but application did not. Matching exclusion fixes the W local-s trigger. Independent review confirmed the predicate and active single-j behavior.
- Human approved stopping the old local MPI2 campaign; all16 owned processes disappeared and evidence/outputs remain. Those results are unaccepted; no automatic restart. Old binary/source snapshots retained.
- Optional GPU/OpenACC, ScaLAPACK, EigenExa, FFTW and LibXC configurations were not exercised. Physical GS/grid/current accuracy and production performance are not established by these smoke tests.
