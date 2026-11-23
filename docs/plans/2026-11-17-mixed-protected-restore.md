# Protected restore in the ordinary mixed-native coordinator

Default direct resume stages native field (including all LD state), native
TDDFT/fast-PR HSE and detached built-in endpoint snapshots before a common MPI
ready vote. Prepare coupling arrays, endpoint adapters, ownership lists and
clocks before publication. Handled staging failures leave the live session usable
and do not set failed or change checkpoint diagnostics. Successful publication
installs native candidates, then prepared Python bindings, then retires old
snapshot handles; cleanup errors are separate diagnostics.

Allow restoration into paused advanced sessions. Preserve exact native-direct/2
format/topology/inventory admission. Reject old candidate-less native libraries
before staging. Legacy restore is an explicit option, never a fallback.

Test the actual ordinary mixed sample: save, advance, inject a late electronic
stage failure after loading, compare full field/LD/native shards plus snapshots
and all coupling/adapter clocks before/after failure, continue originals, retry
and compare the next frames exactly to continuous evolution. Cover internal MPI
TDDFT, field MPI, and native fast-PR HSE; source and installed-wheel verification.
Keep healthy-MPI/paused-host and allocation/process-failure exclusions explicit.
