# Tow-Temperature Model (TTM) calculation

TTM calculation is available during FDTD simulation.

The implementation is based on the following paper:

> "Three-dimensional thermal response of a metal subwavelength tip under femtosecond laser illumination",
> J. Houard, A. Vella, F. Vurillot, and B. Deconihou,
> Phys. Rev. **B** 84, 033405 (2011).

## Usage
In addition to the usual input files for FDTD simulation in SALMON ( *input_fdtd* and *shape.cube* ),
**_ttm.dat_** file is also necessary for FDTD with TTM. An example of the **_ttm.dat_** file is
```
24.5  / g        (J/m^3)
2.44  / Cl       (J/m^3/K)
235   / kappa_e0 (J/m/s/K)
135   / Ce'      (J/m^3/K^2)
46    / zeta_bal (nm)
80    / Tini     (K)
```
If this file exists in the execution directory, SALMON automatically run in the TTM calculation mode.
You can see an example in *src/ttm/example/* directory.

## Output files
Temperature distribution files are generated. Each file corresponds to the distribution at one sampling time.
The control parameter of the sampling interval is common with ordinary FDTD simulation (*obs_samp_em*).
The following files are generated according to the sampling parameter *obs_samp_em*:
```
fort.4001, fort.4002, fort.4003, ...
```
To visualize the data with gnuplot, please see a document in the utility directory (*utility/ttm/README.md*). 
