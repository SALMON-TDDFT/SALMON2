#include <spglib.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  enum { max_atoms = 4096, max_operations = 2048 };
  double lattice[3][3] = {{0}}, positions[max_atoms][3];
  int types[max_atoms], rotations[max_operations][3][3];
  double translations[max_operations][3];
  char symbol[32];
  double cell, x, y, z;
  int n = 0, operation, row, count;
  FILE *input;

  if (argc != 3 || (cell = strtod(argv[2], NULL)) <= 0.0) {
    fprintf(stderr, "usage: %s atom_file cubic_cell_bohr\n", argv[0]);
    return 2;
  }
  input = fopen(argv[1], "r");
  if (!input) return 2;
  while (n < max_atoms && fscanf(input, " '%31[^']' %lf %lf %lf %*d", symbol, &x, &y, &z) == 4) {
    positions[n][0] = x / cell;
    positions[n][1] = y / cell;
    positions[n][2] = z / cell;
    types[n++] = 1;
  }
  fclose(input);
  lattice[0][0] = lattice[1][1] = lattice[2][2] = cell;
  count = spg_get_symmetry(rotations, translations, max_operations, lattice,
                           positions, types, n, 1.0e-6);
  if (count <= 0) return 1;
  fprintf(stderr, "operations=%d atoms=%d\n", count, n);
  for (operation = 0; operation < count; ++operation)
    for (row = 0; row < 3; ++row)
      printf("%2d %2d %2d %.15f\n", rotations[operation][row][0],
             rotations[operation][row][1], rotations[operation][row][2],
             translations[operation][row]);
  return 0;
}
