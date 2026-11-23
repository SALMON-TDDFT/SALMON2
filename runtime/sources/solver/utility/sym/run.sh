#!/bin/bash

rm ./a.out
gfortran get_atomic_number_module.f90  cif_format_module.f90   main.f90

./a.out  <  Si_sample.cif
#./a.out  <  fort.970.cif
