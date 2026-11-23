
### Recoverable material-ID CUBE input errors

Native E/H initialization returns1115 for unavailable, malformed, truncated,
nonfinite or grid-count-mismatched SALMON material-ID CUBE input. Reported tiny
reader staging-buffer allocation failure returns1116; this allocation path has
not been fault-injected and does not cover general OOM or solver allocations.
The reader owns its file unit and checks every value used in integer conversion,
including the final partial record. Existing streaming order and material-ID
conversion are preserved; this is not a general electronic-density CUBE reader.

For external shape input only, one field-group initialization agreement precedes
halo exchange, so a read failure on one rank returns on all ranks. No per-step
collective or generated-shape collective is added. Known unanimous statuses and
recovery protocol1 permit cold cleanup and retry with a new work directory on
the same caller communicator; status-absent standalone behavior is retained.

Local GNU MPI1/2 tests cover retry, caller MPI survival, direct resume and exact
agreement with prior generated-shape fields, point output and checkpoint bytes.
A host-owned unit600 survives rejected shape reads. Successful initialization's
ordinary point-output headers still use a fixed unit600; general output-unit
ownership remains a separate open task. No claim covers arbitrary device failure,
external STOP, process failure or all initialization allocations.
