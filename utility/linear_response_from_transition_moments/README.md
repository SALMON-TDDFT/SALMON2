# tm2sigma.f90
generates the electric conductivity \sigma_{\mu,\nu}(\omega) and the dielectric function \varepsilon_{\mu,\nu}(\omega) from the transition-dipole moments calculated by SALMON [cf. Eq. (13.37) of Ashcroft-Mermin, Solid State Physics].

## Compilation
```
$ ifort -o tm2sigma tm2sigma.f90 
```

## How to use

Execution should be performed in the directory which includes an input file "*input*" (specifying *sysname*), the transition-dipole moment data file "*sysname_tm.data*", and the single-particle energies data file "*sysname_eigen.data*".
The calculation of SALMON must be in the atomic unit (unit_system = 'a.u.').
Run the following command, and then you get the conductivity *sysname_sigma.data* and the dielectric function *sysname_epsilon.data*.
```
$ ./tm2sigma < input
```
For more details, see *example*.
