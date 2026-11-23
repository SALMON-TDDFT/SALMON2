# Point/run bindings implementation plan

Continuation of approved material configuration design sections 3, 6 and 10.
Goal: serial point identity and explicit output mappings without numerical changes.
Architecture: PointRun owns one state per point; caller owns prepared models.
Models implement record() and new_state(**kwargs), including custom Python models.
State proxy delegates numerical operations unchanged. No implicit E/A conversion,
MPI allocation, restart restoration or coupled-step rollback. Tech: Python/NumPy.

- Add failing tests for independent states, exact IDs/order, duplicate creation,
  failed creation retry, borrowed models, duplicate state object, output roundtrip,
  explicit point axes and malformed mappings.
- Implement salmonpy.point_runs.PointRun/PointState/open_point_arrays. IDs use
  existing validator; definition and optional initial digests are lowercase SHA256.
  For a material ID, all occurrences require the same definition digest. For an
  initial ID, all occurrences require the same material and initial digest.
  One run+point means one state lifetime; reopening a point requires a new run.
  close attempts every state, never closes models. No concurrency guarantee.
- Save mappings under metadata.point_run in existing outputs/1; require explicit
  array-name -> nonnegative point axis mapping. Open validates identity and shapes;
  no automatic reorder. User retains responsibility for supplied array ordering.
- Add native LD comparison/independence sample and custom-model tests, paired
  POINT_RUNS guides, public API/install entries and AI manual JSON.
- Run focused/full tests, sdist->wheel and isolated installed sample. Retain local
  verification record. No native rebuild, GitHub Actions, commit or push.
