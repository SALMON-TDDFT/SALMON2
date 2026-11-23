# Common checkpoint admission API

User-approved scope: promote save/restore capability and compatibility inspection to a reusable API, usable by custom backends; retain numerical results and update paired/AI guides.

- Add `salmonpy.checkpointing.CheckpointCapability(format, modes, decomposition, identity='')`, a validated immutable declaration, with detached JSON conversion.
- Add `check_checkpoint_compatibility(comm, groups, *, group=None, field=None, states=None, operation='save', saved=None, required_files=(), available_files=None, local_error=None)`. Groups are ordinary rank/point mappings, not solver objects. `comm=None` supports standalone serial hosts. Capabilities may be declarations or JSON mappings, allowing validation inside the collective.
- No MPI initialization/import, native loading, file I/O, stepping or persistence callbacks. The API validates declarations, group/point coverage, actual decomposition products, operation/saved-record agreement, within-group capabilities, and optional per-rank required-file inventories. Every error is agreed before return.
- An optional identity pins caller-selected implementation/model/input identity. Absent identity preserves existing records; metadata checks cannot infer physical equivalence or authenticate file contents. Same topology only.
- Preserve the existing direct-compatibility/1 record shape (omit an empty identity), so current native-direct/2 archives remain readable. The sample keeps solver-specific lifecycle/ABI checks and uses the common validator for declaration/topology agreement and inventory admission.
- Test serial immutability/identity/layout/invalid declarations, real MPI one-rank errors and differing operations/saved records, plus a working user-material restart example. Verify installed-wheel use and byte equality of existing native continuation paths. No new native kernels or remote Actions.
