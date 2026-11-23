# Common backend-owned checkpoint execution

Approved user scope: call save/restore through a shared entry, keep formats in each backend, support custom code, verify unchanged continuation and synchronize paired/AI guides.

- Add `CheckpointTarget(capability, save=None, restore=None)`: explicit zero-argument callbacks; no backend inheritance or prescribed file format. Validate callbacks inside collective admission, not during binding construction.
- Add `execute_checkpoint` alongside metadata admission. Inputs mirror `check_checkpoint_compatibility`, with targets instead of declarations. Return compatibility JSON plus process-local field/state callback results; never transmit opaque payloads/results over MPI.
- Agree declarations/callback availability before the first callback. Iterate group-local field/point slots, agreeing errors globally after each slot. Different groups may own different point counts. Collective backend callbacks remain responsible for returning/raising consistently within their own MPI group.
- Preflight errors raise InputError with no callbacks. Callback errors become RuntimeFailure on all ranks and stop subsequent slots; already performed mutations are not rolled back. No communicator ownership changes, implicit compiler, native loader or generic archive writer.
- Route actual custom-example I/O and native mixed field/electronic/snapshot operations through the common execution entry. Retain the sample's existing file hashing, lifecycle checks, completion publication and native-direct/2 format. An earlier planning preflight may precede directory setup; execution rechecks admission.
- Verify unit no-callback failures and local opaque results, real MPI unequal-group scheduling/error agreement, custom file restart and old native archives. Use <=4 compute ranks, local packaging/tests, no routine remote Actions.
