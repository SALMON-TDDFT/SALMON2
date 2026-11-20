# Prepare inputs for offline compute nodes

[日本語](PREPARATION.md) · Topic: `offline_inputs` · Paired revision: `2026-11-15.1`

Run preparation on a workstation or management/login node where the facility
permits it. `salmonpy.preparation` collects explicitly selected **local files**,
validates material parameters through your registered resolver, and records an
integrity manifest. It never downloads, installs, builds or submits a job.

## Python host

```python
from salmonpy.preparation import prepare_inputs, verify_inputs

# recipe and registry follow CONFIG.en.md. A registered resolver must validate
# file references relative to the staged base_directory, not the original cwd.
bundle = prepare_inputs(
    recipe, "new-job-inputs", registry=registry, base_directory="downloaded",
    files={"structures/cell.cif": "my-supercell.cif"},
    provenance={"structures/cell.cif": {
        "source": "user-provided supercell", "license": "user-managed"}}
)
verify_inputs(bundle)
```

`files` maps portable relative destination paths to local source paths. Relative
source paths use `base_directory`; absolute source paths are allowed but are not
saved automatically. Parameters must use the destination paths. Include every
required structure, potential, initial-state and nested reference file explicitly.
Directories/symlinks, path traversal, reserved metadata names, case collisions
and existing destinations are rejected. Missing/nonregular sources are reported
together before callbacks/output creation. There is no recursive dependency scan.

`provenance` is optional JSON metadata per asset; absent metadata is marked unknown.
Record source/provider ID, version and license where available. Supplied metadata
is descriptive, not authenticated. Credential secrets do not belong in it.
The current package does not load CIF/ASE structures or resolve grid spacing
through this helper; use the explicit structure adapter linked below. A custom resolver may do so
explicitly under its own validated contract. Opaque parameters cannot be checked
for hidden file dependencies by this collector.

The result contains `materials.json` (recipe), `resolved.json` (immutable catalog
export), copied assets at the requested paths, and `manifest.json` with SHA256,
byte size and provenance for every other file. There are no file symlinks or
source machine paths inserted by the collector. Copying preserves contents,
not executable permissions or a full virtual environment. You may include wheel
and runtime archives explicitly; this does not resolve dependencies or certify
OS/Python/MPI compatibility. Follow [RUNTIME.en.md](RUNTIME.en.md).

The trusted resolver runs against the temporary staged directory. It must not
retain its temporary absolute path, use original source paths, modify staged
files, launch solvers, or acquire data. Only portable relative references and
resolved values should be returned. Resolver side effects are not sandboxed or
rolled back. The directory is published only after validation; failures clean
staging and never replace an existing destination.

## Transfer and verify

Copy the whole directory, preserving names. Record the SHA256 of `manifest.json`
separately when preparing it. On the compute node:

```sh
python -m salmonpy.preparation transferred-inputs --manifest-sha256 RECORDED_SHA256
```

The Python equivalent is `verify_inputs(directory, manifest_sha256=...)`.
Verification reads only local files and does not import or call registered
backends. It rejects missing/extra files or directories, symlinks, malformed
manifests and hash/size mismatches. It returns a detached manifest dictionary;
the CLI prints a short integrity-only summary. Omit `manifest_sha256` only when
an internal consistency check is enough. Hashes within a modified manifest cannot
authenticate its contents; keep the recorded manifest hash through a trusted path.
OS errors propagate; invalid content raises `InputError`.

Passing verification does not prove all scientific inputs/dependencies are
present, that the material is physically appropriate, or that a solver will run.
Load/validate with the chosen backend and test a small case on the target node.
Use a new directory for outputs: adding results inside the bundle deliberately
fails its exact-inventory check. No network blocking is imposed on arbitrary user
code; the verifier itself makes no network calls.

## Executable example

Create a local `coefficient.json` with `{"coefficient": 2.0}`, then run:

```sh
python -m salmonpy.examples.prepare_inputs coefficient.json new-inputs
python -m salmonpy.preparation new-inputs
```

This example has an explicit custom resolver and records a synthetic scalar
response file. It demonstrates staging and integrity, not a Si/DFT calculation.

For an explicitly registered structure/grid resolver see [STRUCTURE.en.md](STRUCTURE.en.md) and `examples.prepare_structure`. The generic collector does not automatically interpret structures.

## Batch jobs: optional user-run checks and startup failures

The user may explicitly run this optional check on a management/login node only
where facility policy allows it. SALMONpy does not access that node, arrange the
check automatically, require its report for submission, or submit the job.
Skipping this optional step does not prevent an ordinary calculation. To check
local files and registered material parameters manually:

```sh
python -m salmonpy.preparation new-inputs --resolve-builtins --report preflight.json
```

`--resolve-builtins` is for built-in material recipes only. For your own backend,
call `preflight_inputs(bundle, registry=your_explicit_registry)` in Python; it
returns a JSON-ready report with `status`, errors (`phase`, `type`, `message`),
manifest hash and declared limitations. Trusted resolvers must follow the
no-solver/no-download preparation contract above. With no registry, the check
is file integrity only. Runtime/MPI compatibility, native namelist semantics,
GS compatibility, convergence and memory capacity are not certified.

The CLI prints a JSON report and exits 0 on success or 2 on invalid input.
`--report` creates a new file outside the immutable bundle; it never overwrites
an existing report. Report-writing errors also fail the command. Keep this
report with the submitted job records. This changes the CLI reporting format;
`verify_inputs` itself retains its exception and manifest-return contracts.

For a batch script, check after transfer before launching the calculation:

```sh
set -eu
python -m salmonpy.preparation transferred-inputs --manifest-sha256 RECORDED_SHA256 --report job-input-check.json
# Launch your calculation here only if the previous command succeeded.
```

If every rank needs to check its local filesystem view, the runnable example is:

```sh
mpiexec -n 2 python -m salmonpy.examples.mpi_input_preflight transferred-inputs job-ranks-check.json --manifest-sha256 RECORDED_SHA256
```

This example only checks inputs; it starts no solver. It gathers all rank
reports before deciding success/failure and exits every rank with the same
status. A missing file on one rank prevents solver creation for the whole job.
Use the same pattern before state creation in a production MPI script. Invalid
input is recorded, fixed outside the job, then resubmitted. No automatic input
editing or interactive retry is needed. Numerical warnings such as unconverged
SCF or an exceeded advisory memory budget remain warnings under their existing
policies; this check introduces no new numerical stopping thresholds.

### Ordinary mixed calculation

`examples.mpi_mixed_native` accepts optional `--input-check STAGED_INPUT_BUNDLE`
and `--input-manifest-sha256 RECORDED_SHA256`. Every world rank verifies the
bundle before creating solver states, and the new output directory receives
`input-preflight.json`. Failure exits every rank with code2. Output-directory
creation failures are also agreed before any barrier, so peers do not wait for
a rank that could not create it. Existing output directories are preserved.

This startup check is integrity-only. It does not automatically take calculation
arguments from the bundle or prove that separate `--layers`, runtime and GS
arguments match its recipe. Keep staged references and command arguments aligned;
use the material-resolver check before submission. Later runtime/numerical
failures retain the sample's existing fatal/transaction policies.

The ordinary mixed sample also records startup declaration reads even when no
`--input-check` is supplied. Each rank reads its explicitly selected layers,
model-group and TDDFT material-catalog files; `OUTPUT/startup-inputs.json` records
rank, input phase, path, exception type and message. Read/JSON errors are gathered
and every rank exits2 before solver creation. Successful startup records describe
these reads only, not later native validation or overall calculation success.
This requires no management-node action or prior report. Other native/runtime
failures retain their separate policies.
