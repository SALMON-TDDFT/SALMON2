# Explicit pseudopotential downloads

[日本語](PSEUDOPOTENTIALS.md)

Download with the user's authorization on a connected preparation machine, then transfer files with the inputs to compute nodes. Calculations and imports never trigger this operation.

Fetch the pinned Apache-2.0 Si file used by the bundled example:

```sh
python -m salmonpy.pseudopotentials NEW_DIRECTORY --silicon-example --allow-download
```

Or use Python; create the parent directory beforehand:

```python
from salmonpy.pseudopotentials import download, silicon_reference
path = download(silicon_reference(), 'NEW_DIRECTORY', allow_download=True)
```

The new directory contains the file and `download.json` recording provenance, license, URLs and SHA256. Existing destinations are rejected. HTTPS is required, including redirects; downloads are limited to 64 MiB and must match the pinned SHA256. A mismatch does not publish a calculation file or update the reference.

For other files, construct `PseudopotentialReference(filename, atomic_number, url, sha256, license, source)`. A local JSON with these fields is accepted by `--reference FILE.json`. Obtain the expected hash from a trusted source; do not automatically trust a hash computed from an unknown download. License metadata records the publisher's declared terms and does not grant redistribution rights for arbitrary URLs. The caller must check permission.

This verifies file identity, not format, XC, valence count, local channel or grid convergence. Supply those conditions through `Pseudopotential.from_file(path, ...)`. GS/2 continues to store references without triggering downloads.

## Selecting PseudoDojo data and SOC

The current catalog pins PBE v0.4: scalar-relativistic H, C, N, O, Al, Si, Ti, Se and W; fully relativistic W and Se for WSe₂, plus Ni with 3s/3p semicore (18 valence electrons) for the [short SOC test](NICKEL_SOC.en.md). Catalog entries record file identity/header metadata, not successful physical calculations for every element.

```sh
python -m salmonpy.pseudopotentials --list --element W --xc PBE --relativity fully-relativistic
python -m salmonpy.pseudopotentials NEW_W_DIRECTORY --catalog-id pseudodojo-pbe-0.4-fr-standard-w --allow-download
```

```python
from salmonpy.pseudopotentials import catalog_reference, download
from salmonpy.dft import Pseudopotential
reference = catalog_reference('pseudodojo-pbe-0.4-fr-standard-w')
path = download(reference, 'NEW_W_DIRECTORY', allow_download=True)
potential = Pseudopotential.from_file(path, atomic_number=74, local_channel=3,
                                     reference=reference)
```

Version, XC, relativity, valence count, NLCC, SOC flag and upstream commit are recorded. Header metadata and hash are checked during download and Python loading. No element-only choice silently substitutes a different version. Grid convergence and local-channel suitability require separate checks.

PseudoDojo data are CC-BY-4.0: credit PseudoDojo and ONCVPSP and cite their work. See the [official project](https://github.com/PseudoDojo/pseudodojo) and [PseudoDojo paper](https://doi.org/10.1016/j.cpc.2018.01.012). The wheel includes references, not a bulk copy of these UPFs.

**End-to-end SOC support is in progress.** FR acquisition/loading and explicit SOC DFT inputs are implemented. [Ni validation](NICKEL_SOC.en.md) covers independent SCF states and fresh-process SCF restart, plus independent TDDFT states and direct restoration within the same process or after process exit with unchanged MPI layout. Changed-layout SOC RT restart remains unvalidated. A SOC request never falls back to scalar-relativistic propagation. Monolayer WSe₂ is the first verification target.
