# MPI point identity implementation plan

Approved material-backend design M5/M6 continuation; preserve existing integer
transport IDs and numerical communication/update order. User IDs are opaque strings.

1. Add tests of point-run record plus explicit integer->string point_map validation,
   model definition conflict, ordered arrays, detached callback metadata and collective
   rank disagreement rejection before prepare/create. Legacy no-metadata path remains.
2. ResponseSession accepts optional point_record and point_map together. Validate
   against plan; all ranks supply identical records. Bindings sharing a cache key
   must share full material/initial identity. Collective descriptor comparison covers
   metadata before callbacks. Group-only SolverContext.point_records exposes immutable
   per-transport-point identities plus run ID; user callbacks retain integer point_id.
3. Collective save_arrays: field leader alone supplies output arguments; preflight
   agreement before write, then propagate I/O failure to all ranks. Add point_run,
   point_axes and MPI ownership metadata to outputs/1. open_point_arrays remains reader.
   No restart/ABI/capability certification. Coupled factory forwards optional metadata.
4. Run real local 3-rank sample (one field, two material groups), compare direct and
   legacy paths and rejection/deadlock cases; optionally 4 ranks for group reduction.
   Use bounded MPI runner, one thread/rank, no Actions/native rebuild.
5. Update paired guides, AI JSON, sample, full tests, sdist->wheel and installed
   sample. No commit/push. Keep shared native state work and restart work open.
