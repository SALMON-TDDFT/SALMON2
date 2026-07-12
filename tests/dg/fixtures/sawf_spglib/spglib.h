#ifndef SPGLIB_H
#define SPGLIB_H

int spg_get_multiplicity(const double lattice[3][3],
                         const double position[][3], const int types[],
                         int num_atom, double symprec);

int spg_get_symmetry(int rotation[][3][3], double translation[][3],
                     int max_size, const double lattice[3][3],
                     const double position[][3], const int types[],
                     int num_atom, double symprec);

#endif
