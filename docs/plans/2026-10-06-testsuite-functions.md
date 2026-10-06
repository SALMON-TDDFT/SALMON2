# Functional test consolidation

Approved design: one entry point per exchange/ACE, SR/FFT, WF localization, DC-LCFO, LCFO RT, and PBEh/rVV10 function. Retain independent numerical oracles and MPI comparisons inside each function. Move platform/GPU checks outside ordinary regression. Validate roots, syntax, CMake registration and representative compiled functions; record unexecuted cases.
