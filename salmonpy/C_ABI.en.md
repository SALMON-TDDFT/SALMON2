# Custom material C ABI (preview 1)

[日本語](C_ABI.md) · Topic: `c_abi` · Paired revision: `2026-11-14.1`

Use one native material contract when plugging a model into Python or driving it
from a C/C++/Fortran host. Public declarations are in
[include/salmonpy_material.h](include/salmonpy_material.h) and
[include/salmonpy_material_iso.f90](include/salmonpy_material_iso.f90).
Python loads plugins through `salmonpy.c_material.CMaterial`.

Scope: fixed-dt, CPU, serial endpoint-A materials. Preview ABI1 does not expose
all existing SALMON/LD/SBE internal C symbols. Native TDDFT's shared globals,
one-lifetime limit and multiple-context implementation are unchanged.

## Data and ownership

- Input is Ac=A/c; output is physical charge-current density j. Laboratory frame,
  atomic units, three contiguous IEEE754 binary64 doubles.
- Separate positive int64 model/state handles are library-owned; 0 is invalid.
  Never transfer a handle to another library's callback table.
- Models copy/own configuration; states own all mutable history. Do not retain
  caller configuration/input/output buffers. Input/output buffers must not overlap.
- Hosts serialize calls. ABI1 supplies no internal MPI collectives or thread-safety guarantee.
- Time is initially time0 and then time0+n·dt. dt is fixed and positive. Initial
  A is explicit and initial current must be physical.
- Declare input lag and use of initial A in the table; Python preserves both in TimeSpec.

Configuration is opaque bytes whose schema the plugin defines and documents.
The shared ABI never guesses model-specific parameter or file formats.

## Query and callbacks

`salmonpy_material_query(version, table_capacity, &table)` negotiates version and
size before writing the table. Hosts check returned status, abi_version,
struct_size, mandatory callbacks, backend and timing metadata before calling.
This is the platform C ABI, not an on-disk/wire format. Table metadata/backend
string remain valid until library unload.

| Callback | Responsibility |
|---|---|
| model_create | Copy configuration and create model; failure sets handle0 |
| model_destroy | Return BUSY with live states, otherwise release model |
| state_create | Independent history from model,dt,time0,initial A; failure sets handle0 |
| state_observe | Accepted clock/current without advancing |
| state_step | Next endpoint A to current/clock at that endpoint |
| state_destroy | Release all resources; stale handles are rejected |

Callbacks return int32 status: 0 success, 1 argument, 2 handle, 3 memory, 4 numerical,
5 busy, 6 ABI, 7 internal. Failures leave accepted state unchanged; hosts ignore
clock/current outputs on error. Failed constructors release partial resources.
Optional caller-sized error text is truncated and NUL-terminated when capacity
permits. Never throw C++ exceptions across the boundary or stop/error stop a
Fortran plugin's host process.

## Python client

```python
from salmonpy.c_material import CMaterial
from salmonpy.materials import check_material, bind_material
from salmonpy.crystal import respond

# The provider specifies configuration bytes, library and SHA256.
with CMaterial(library_path, library_sha256, configuration_bytes) as model:
    report = check_material(model, drive)
    response = respond(model, drive)
    binding = bind_material(model, key='sample', identity=model.identity)
    # Create/close ResponseSession inside this model lifetime.
```

CMaterial verifies SHA before explicit loading and copies configuration. Its
identity binds library/configuration hashes and backend name. Import never loads
plugins, compiles or initializes MPI. Avoid concurrent operations on Python
models/states. Native calls are serialized process-wide, including library
aliases; this does not make Python object lifecycle operations thread-shareable.

Model.close prevents new states and defers native destruction until existing
states close. Existing states remain usable. Use context managers or explicit
close; garbage collection does not release forgotten handles. State.close is
idempotent. Direct C hosts destroy states before models.

Nonfinite Python inputs fail before native invocation. Nonzero callbacks raise
CABIError(RuntimeFailure), retaining the status code and cached accepted clock/
current. Plugins must also provide atomic failure. Success with invalid clock/
current invalidates the Python state; no generic native rollback is invented.

Plugins are trusted executable code. Hashes do not authenticate publishers or
isolate native faults. Do not preload outside CMaterial or replace loaded plugin
files; use a new immutable file for each build. Pinning is retained after a
successful dlopen even when ABI validation fails. Path/file/image aliasing that
mixes a loaded image with a different asserted build hash is rejected.

## C++ and Fortran host examples

Explicit developer build, into a new destination:

```sh
python tools/build_c_material_example.py NEW_BUILD --cc cc --cxx c++ --fc gfortran
```

Installed wheels carry this tool under salmonpy/_runtime_source/tools and include
headers, ISO declarations and example sources. Building the example needs C/C++/
Fortran compilers; ordinary wheel use/import never invokes one. The example build
uses GNU-compatible Fortran and Unix tooling; actual current tests ran on Apple
Silicon macOS.

The example configuration is two native C doubles {gain,decay}, supplied as
`struct.pack('=dd',gain,decay)` in Python or a contiguous pair in Fortran.
Its j_new=decay*j_old+gain*(A_new−A_old) recurrence is a manufactured interface
fixture, not a real-material model. The explicit build audits sources and all
application objects before linking, preserves source snapshots/logs/hashes, and
runs separately compiled C++/Fortran hosts.

Both hosts passed independent state, copied configuration, invalid/stale handle,
BUSY and post-NaN atomicity checks. Python reference difference was 0; the same
C plugin on 2/3rank point-MPI distribution also differed by 0. This is not physical
validation or proof of native TDDFT reentrancy.
