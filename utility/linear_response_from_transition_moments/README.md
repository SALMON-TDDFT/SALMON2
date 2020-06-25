# tm2lr.f90
generates the electric conductivity \sigma_{\mu,\nu}(omega) and the dielectric function \varepsilon_{\mu,\nu}(omega)) from the transition-dipole moments calculated by SALMON.
cf. Eq. (13.37) of Ashcroft-Mermin, Solid State Physics.

## Compilation
```
$ ifort -o tm2lr tm2lr.f90 
```

## How to use

Execution should be performed in the directory which includes an input file *input*, the transition-dipole moment data file *<sysname>_tm.data*, and the single-particle energies data file *<sysname>_eigen.data*.
The calculation of SALMON must be in the atomic unit (unit_system = 'a.u.').
Run the following command, and then you get *<sysname>_sigma.data*.
```
$ ./tm2lr < input
```
For more details, see *example/README.md*.
