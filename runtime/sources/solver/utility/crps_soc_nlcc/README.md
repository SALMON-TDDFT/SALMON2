# CRPS SOC/NLCC pseudopotential generator

The retained utility contains the generator, output module
and example generation inputs. Material scans, measured outputs and standalone
tests are archived outside the submitted source. The public regression uses only
small-cell Si; it does not validate these generated pseudopotentials.

- `crps_soc_nlcc.f`: multi-element relativistic generator with j-resolved projectors.
- `crps_output.f90`: extended KY output writer.
- `inputs/{Ni,W,Se}.nml`: generation inputs, not additional SALMON regression systems.

## Generate into a new empty directory

GNU Fortran example, from this directory:

```sh
gfortran -std=legacy -O0 -fcheck=all -fbacktrace -Wno-align-commons \
  crps_output.f90 crps_soc_nlcc.f -o crps_soc_nlcc
mkdir W-run
(cd W-run && ../crps_soc_nlcc ../inputs/W.nml > stdout.txt 2> stderr.txt)
```

The generator refuses to overwrite existing output. It inherits relativistic
LDA exchange and PZ correlation; it is not a PBE generator. Changes to valence
configuration, cutoff radii or core treatment require independent validation.
The radial grid is0.01Angstrom, so requested radii are rounded to that grid.
The c2 search requires a unique bracket and does not silently choose a branch.

Output uses `KY_SOC_NLCC_V1` with j-resolved ionic projectors and the scalar
NLCC density. NC-MGGA23 auxiliary kinetic-energy and spin-current fields are
outside this integration. SALMON can still read the ordinary ionic/NLCC data
from existing V2 tables; the auxiliary trailer is not used.

SALMON detects the extended KY header and resolves omitted noncollinear/SOC
settings; inconsistent explicit settings are rejected rather than silently
scalar-averaged. Rebuild SALMON to use the extended reader. No automatic copying
into shared pseudopotential directories or changes to running binaries occur.

Physical inputs remain explicit: element, occupations, reference orbitals,
cutoff radii, NLCC core radius and any required local channel. Input conveniences
may derive Lmax/NA/lloc and output filenames, but do not replace physical checks.
