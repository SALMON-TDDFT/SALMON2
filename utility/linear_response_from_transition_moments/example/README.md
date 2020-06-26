# Example: calculation of Lithium

## 1. SALMON
Using an input file for SALMON (*input_salmon.dat*) and a pseudopotential file for Li (*Li_rps.dat*), calculate the transition-dipole moments (*Li_tm.data*) and the energy eigenvalues (*Li_eigen.data*).

## 2. tm2sigma

Run the following command, and then you get *Li_sigma.data* and *Li_epsilon.data*.
```
$ ./tm2sigma < input_tm2sigma.dat
```

## 3. Comparison with *theory = 'TDDFT_response'*
You can calculate the linear response functions using *input_salmon.dat* with *theory = 'TDDFT_response'* [cf. Bertsch et al., PRB 62, 7998 (2000)].
Compare *Li_response.data* with *Li_sigma.data* or *Li_epsilon.data*.
