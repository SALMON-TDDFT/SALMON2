# Automatic MPI factorization recovery

### Recoverable automatic MPI factorization failure

Native E/H returns1117 when the existing automatic MPI distributor cannot assign
its rank count. The current automatic domain routine factors2/3/5; this is a limit
of that algorithm, not a blanket rank-count restriction. A supported explicit
spatial partition remains usable (tested with7 ranks). No Python/outer factor
whitelist is added and the distribution algorithm is unchanged. Correct the
partition and create a new work directory on the same caller communicator.

The internal Fortran helper has an optional status (0success,1unassigned prime
factor,2no allowed non-orthogonal k/orbital split). Without status it retains the
legacy diagnostics/STOP behavior. Matching Fortran modules must be used when
building callers; this is source compatibility, not stable Fortran binary ABI.
The C initialization signature and recovery protocol1 are unchanged. Python
releases a failed lifetime only for a recognized unanimous1117 with protocol1.
The existing cold cleanup/status agreement is reused; no new collective or
per-step work is added. This does not handle allocation, Fujitsu topology, MPI
failure or every specialized initialization stage.

Local GNU tests cover unchanged automatic choices on MPI1/2 and helper failure/
continuation on1/2/7. Native MPI7 repeats raw/public automatic-failure rejection,
keeps parent MPI alive, releases owned children and retries explicit decomposition;
32-step fields, point bytes and direct checkpoint/resume match the prior runtime.
This is local contract/regression evidence, not multi-node scaling certification.
