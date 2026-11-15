# Configured TDDFT/LD response bindings

Continuation of approved M5 adapter work. Public bind_configured_ld_endpoints
makes A-endpoint-to-E conversion explicit in its name and reuses existing LD
provider/backend. Public bind_configured_tddft requires an explicit managed runtime,
borrows an existing configured initial snapshot, checks executable identity and
reuses TDDFTResponseBackend. Both attach material_record for response_identity.

Tests first: wrong models/runtime mismatch/no implicit state creation; LD direct
E-step equality and independent states; prepared TDDFT borrows without extra copy.
Identity includes material record plus numerical provider/source manifest, input,
executable and worker settings, never object address or absolute runtime path.
No worker on binding/prepare; state creation starts existing worker. Model owner
must outlive session; no provider close may close borrowed model. No implicit DFT.

Local native/3-rank installed-wheel tests: synthetic LD and fresh packaged short
Si TDDFT, compare MPI with direct paths and verify IDs. Max three MPI ranks,
one thread per worker, no more than two simultaneous TDDFT workers. Preserve
native update/sign contracts; no full runtime build or Actions. Paired docs/AI
JSON/example, full tests and sdist->wheel. No commit/push.
