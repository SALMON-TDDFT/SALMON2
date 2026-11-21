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

## Separate code-and-input checker (optional)

The user may run this on a management/login node where facility policy permits:

```sh
python -m salmonpy.preparation input.json --code calculation.py --report check.json
```

For a staged input directory, use that directory instead of `input.json`, with
`--manifest-sha256 RECORDED_SHA256` if available. Python equivalent:
`check_job(code, inputs, manifest_sha256=None)`.

This command **only checks and exits**, even on success. It never imports or
executes the supplied calculation, follows its imports, invokes material plugins,
starts MPI, builds a runtime, launches a solver or submits a job. The report is a
new file; existing reports are preserved. Run your calculation separately through
your scheduler on allocated compute nodes. No prior report is required to run.
The mixed calculation sample has no `--input-check` option or checker pipeline.
The previous MPI preflight example was removed to avoid confusing a check with a
calculation. Ordinary backend validation and runtime error handling remain.

### Deliberately limited scope

- Python source: syntax/compilation only, with file hash. Compilation does not
  execute the compiled code or check installed imports.
- A JSON file: syntax, duplicate keys and finite JSON values only; no scientific
  parameter schema is inferred.
- A staged bundle: inventory, size and hashes through `verify_inputs`.

The report lists `checks`, `errors` and `not_checked`. Exit0/status passed means
only these checks passed; it is **not permission to compute on a management
node**, a promise of successful execution, or a physical-validity certificate.
It does not prove that code uses the supplied input, discover computed paths,
validate arbitrary custom model combinations, certify native namelist semantics,
verify runtime/MPI/GS compatibility, or establish convergence/memory capacity.
No automatic repair or execution follows a successful check. Invalid code and
input are reported independently; exit2 indicates a detected problem.

Do not attempt to support every user-program combination in this checker.
Future scientific checks should name their exact supported schema/backend and
scope, with unknown combinations reported as unverified rather than forbidden.
Existing `preflight_inputs(..., registry=...)` is a separate explicitly trusted
resolver-based material check; `--code` deliberately does not combine with
`--resolve-builtins`, so arbitrary resolver code cannot run through this entry.


## Preparation-only FDTD on a management node

`prepare_fdtd(calculation, new_directory, mpi=None)` writes native inputfile and
material assets to a new directory and stops. It needs no executable, creates no
subprocess, initializes no MPI, builds nothing, runs no solver and submits no job.
It returns a Path without an execution method; successful checks do not continue
automatically into a calculation.

```python
from salmonpy.preparation import prepare_fdtd
from salmonpy.fdtd import FDTDMPI
from salmonpy.examples.prepare_fdtd import calculation  # small LD connection fixture

folder = prepare_fdtd(calculation(), "fdtd-input", mpi=FDTDMPI((2, 1, 1)))
# Stop here. Transfer, submission and execution are separate actions.
```

```sh
python -m salmonpy.examples.prepare_fdtd NEW_DIRECTORY --mpi-grid 2 1 1
```

`salmonpy-prepared-fdtd.json` records prepared-not-run status, decomposition, rank
count and input hashes/sizes, without site-specific executable/launcher paths.
Existing directories are never overwritten. This is a native SALMON input folder,
not the prepare_inputs material archive; verify_inputs does not accept this
format. It is not a physical convergence or runtime-compatibility certificate.

Transfer the folder if required and explicitly launch a matching audited SALMON
from the site's allocated batch job. Verify its recorded executable SHA256 at
startup rather than repeating source/object audits. Run in the input directory
and pass inputfile through stdin. For example, the launch portion of an allocated
Slurm job might use `srun --ntasks=2 "$SALMON_EXECUTABLE" < inputfile > stdout.log 2>&1`.
This is a site-specific illustration, not a management-node execution instruction
or automatic submission. Set queues, accounts, limits and threads according to
the site's rules in the job script.

## Site-specific launcher arguments

FDTDMPI defaults to rank_arguments=('-n','{ranks}'), preserving existing commands.
command(executable) only returns argv; it neither launches nor checks a binary.

```python
from salmonpy.fdtd import FDTDMPI
mpi = FDTDMPI((2, 1, 1), launcher=("site-launch",),
              rank_arguments=("--tasks={ranks}",))
argv = mpi.command("/site/salmon")  # preview only
# ('site-launch', '--tasks=2', '/site/salmon')
```

Use rank_arguments=() to append no rank-count arguments when the site obtains
counts from its allocation. The actual count must still match the input
decomposition. No arbitrary launcher semantics are inferred or certified. Only
{ranks} is substituted; other strings remain literal, without shell evaluation or
environment expansion. For jsrun, configure resource sets and tasks per set
consistently with MPI ranks according to site rules. Renaming the launcher alone
does not necessarily configure a different site's execution model.

SalmonExecutable.run remains the explicit child-process path from serial Python;
nested-launch protection inside MPI ranks stays enabled. For Python itself as a
caller-owned MPI application, use ParallelContext and native embedded components
collectively as described in COUPLED.en.md. This addition does not unify the
separate DFT/TDDFT batch paths.
