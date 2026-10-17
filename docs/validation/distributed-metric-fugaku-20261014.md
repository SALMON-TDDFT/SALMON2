# Distributed ACE metric driver on Fugaku

Base commit: 63ab85f164519f92de1c04ceca211502523a2165.

Use PZHEEVD for the distributed Hermitian metric eigensolve, round queried workspaces upward, and allocate real/integer workspace with conservative local bounds. Distributed matrix ownership and the existing layout are retained.

Evidence observed through authenticated Fugaku OOD:
- Independent Release build with MPI/HSE/ScaLAPACK enabled: exit 0. Binary SHA256 ea46edcde1424afb5830e232bdbde0c71d4f1f5735663373b2c8a6a94710ec91.
- Minimal 128-state, 8-rank PZHEEVD test (52258924): all 8 ranks INFO=0; eigenvalue maximum error 4.440892098500626e-16.
- Real Gamma validation (52258948), 8 nodes, 8 MPI, 48 OMP: spatial and orbital layouts, both conventional and learned K8 ACE, completed 160 steps and passed the completion gate (all-rank exit, memory records, end SALMON, finite current/energy).
- Combined spatial/orbital conventional ACE completed 160 steps and passed the same gate at 07:47 JST. Combined learned ACE was still running at step 33; its final result is not established here.

This evidence does not establish all-layout numerical accuracy, SVD quality, production capacity, or scaling performance. The original heap-corruption origin is not proven; the ScaLAPACK abort location alone is not a root-cause diagnosis. No precision thresholds were relaxed. Mac heavy validation was not run.

Original source revision, remote binaries, failed output and run locks are preserved. This commit records the tested metric-driver change; production acceptance remains subject to the planned numerical/capacity gates.
