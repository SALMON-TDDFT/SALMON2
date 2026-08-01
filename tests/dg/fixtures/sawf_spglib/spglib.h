#ifndef SPGLIB_H
#define SPGLIB_H

typedef struct {
  int spacegroup_number;
  int hall_number;
  char pointgroup_symbol[6];
} SpglibDataset;

SpglibDataset *spg_get_dataset(const double lattice[3][3],
                               const double position[][3], const int types[],
                               int num_atom, double symprec);
void spg_free_dataset(SpglibDataset *dataset);

int spg_get_multiplicity(const double lattice[3][3],
                         const double position[][3], const int types[],
                         int num_atom, double symprec);

int spg_get_symmetry(int rotation[][3][3], double translation[][3],
                     int max_size, const double lattice[3][3],
                     const double position[][3], const int types[],
                     int num_atom, double symprec);

#endif
