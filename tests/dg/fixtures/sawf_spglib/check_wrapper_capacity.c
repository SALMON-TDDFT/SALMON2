#include <stdio.h>

int salmon_sawf_spglib_get_symmetry(
    const double *lattice_columns, const double *fractional_positions,
    const int *species, int num_atoms, double tolerance, int output_capacity,
    int *rotations, double *translations, int *num_operations);

int main(void) {
  const double lattice[9] = {1.0, 0.0, 0.0, 0.0, 1.0,
                             0.0, 0.0, 0.0, 1.0};
  const double positions[3] = {0.0, 0.0, 0.0};
  int species[1] = {8};
  int rotations[18];
  double translations[6];
  int count = 0;
  int status = salmon_sawf_spglib_get_symmetry(
      lattice, positions, species, 1, 1.0e-8, 1, rotations, translations,
      &count);

  if (status != 3 || count != 2) {
    printf("capacity status=%d count=%d\n", status, count);
    return 1;
  }
  species[0] = 6;
  status = salmon_sawf_spglib_get_symmetry(
      lattice, positions, species, 1, 1.0e-8, 2, rotations, translations,
      &count);
  if (status != 0 || count != 2 || rotations[12] != -1 ||
      rotations[10] != 1) {
    printf("layout status=%d count=%d r01=%d r10=%d\n", status, count,
           rotations[12], rotations[10]);
    return 1;
  }
  puts("PASS capacity");
  return 0;
}
