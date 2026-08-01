#include <spglib.h>

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

enum {
  SAWF_SPGLIB_SUCCESS = 0,
  SAWF_SPGLIB_INVALID_INPUT = 1,
  SAWF_SPGLIB_COUNT_FAILURE = 2,
  SAWF_SPGLIB_CAPACITY_TOO_SMALL = 3,
  SAWF_SPGLIB_API_FAILURE = 4,
  SAWF_SPGLIB_ALLOCATION_FAILURE = 5
};

int salmon_sawf_spglib_get_dataset_metadata(
    const double *lattice_columns, const double *fractional_positions,
    const int *species, int num_atoms, double tolerance,
    int *spacegroup_number, int *hall_number, char *pointgroup_symbol,
    int pointgroup_capacity) {
  SpglibDataset *dataset = NULL;
  double spg_lattice[3][3];
  double(*spg_positions)[3] = NULL;
  int atom;
  int cartesian;
  int vector;
  int status = SAWF_SPGLIB_INVALID_INPUT;

  if (spacegroup_number == NULL || hall_number == NULL ||
      pointgroup_symbol == NULL || pointgroup_capacity < 2 ||
      lattice_columns == NULL || fractional_positions == NULL ||
      species == NULL || num_atoms <= 0 || !isfinite(tolerance) ||
      tolerance <= 0.0) return SAWF_SPGLIB_INVALID_INPUT;
  *spacegroup_number = 0;
  *hall_number = 0;
  pointgroup_symbol[0] = '\0';
  if ((size_t)num_atoms > SIZE_MAX / sizeof(*spg_positions))
    return SAWF_SPGLIB_ALLOCATION_FAILURE;
  spg_positions = malloc((size_t)num_atoms * sizeof(*spg_positions));
  if (spg_positions == NULL) return SAWF_SPGLIB_ALLOCATION_FAILURE;
  for (vector = 0; vector < 3; ++vector) {
    for (cartesian = 0; cartesian < 3; ++cartesian) {
      double value = lattice_columns[3 * vector + cartesian];
      if (!isfinite(value)) goto cleanup;
      spg_lattice[vector][cartesian] = value;
    }
  }
  for (atom = 0; atom < num_atoms; ++atom) {
    if (species[atom] <= 0) goto cleanup;
    for (cartesian = 0; cartesian < 3; ++cartesian) {
      double value = fractional_positions[3 * atom + cartesian];
      if (!isfinite(value)) goto cleanup;
      spg_positions[atom][cartesian] = value;
    }
  }
  dataset = spg_get_dataset(spg_lattice,
      (const double(*)[3])spg_positions, species, num_atoms, tolerance);
  if (dataset == NULL) {
    status = SAWF_SPGLIB_API_FAILURE;
    goto cleanup;
  }
  *spacegroup_number = dataset->spacegroup_number;
  *hall_number = dataset->hall_number;
  if (snprintf(pointgroup_symbol, (size_t)pointgroup_capacity, "%s",
               dataset->pointgroup_symbol) >= (size_t)pointgroup_capacity) {
    status = SAWF_SPGLIB_CAPACITY_TOO_SMALL;
    goto cleanup;
  }
  status = SAWF_SPGLIB_SUCCESS;

cleanup:
  if (dataset != NULL) spg_free_dataset(dataset);
  free(spg_positions);
  return status;
}

int salmon_sawf_spglib_get_symmetry(
    const double *lattice_columns, const double *fractional_positions,
    const int *species, int num_atoms, double tolerance, int output_capacity,
    int *rotations, double *translations, int *num_operations) {
  double spg_lattice[3][3];
  double(*spg_positions)[3] = NULL;
  int(*spg_rotations)[3][3] = NULL;
  double(*spg_translations)[3] = NULL;
  int required_operations;
  int returned_operations;
  int atom;
  int vector;
  int cartesian;
  int operation;
  int row;
  int column;
  int status = SAWF_SPGLIB_INVALID_INPUT;

  if (num_operations == NULL) return SAWF_SPGLIB_INVALID_INPUT;
  *num_operations = 0;
  if (lattice_columns == NULL || fractional_positions == NULL ||
      species == NULL || num_atoms <= 0 || output_capacity < 0 ||
      !isfinite(tolerance) || tolerance <= 0.0) {
    return SAWF_SPGLIB_INVALID_INPUT;
  }

  /* Fortran stores A(cartesian,vector), so each lattice vector is a column.
   * Spglib stores each lattice vector as a C row. Copying with swapped
   * mathematical indices preserves that convention independently of layout. */
  for (vector = 0; vector < 3; ++vector) {
    for (cartesian = 0; cartesian < 3; ++cartesian) {
      double value = lattice_columns[3 * vector + cartesian];
      if (!isfinite(value)) return SAWF_SPGLIB_INVALID_INPUT;
      spg_lattice[vector][cartesian] = value;
    }
  }

  if ((size_t)num_atoms > SIZE_MAX / sizeof(*spg_positions)) {
    return SAWF_SPGLIB_ALLOCATION_FAILURE;
  }
  spg_positions = malloc((size_t)num_atoms * sizeof(*spg_positions));
  if (spg_positions == NULL) return SAWF_SPGLIB_ALLOCATION_FAILURE;
  for (atom = 0; atom < num_atoms; ++atom) {
    if (species[atom] <= 0) goto cleanup;
    for (cartesian = 0; cartesian < 3; ++cartesian) {
      double value = fractional_positions[3 * atom + cartesian];
      if (!isfinite(value)) goto cleanup;
      spg_positions[atom][cartesian] = value;
    }
  }

  required_operations = spg_get_multiplicity(
      spg_lattice, (const double(*)[3])spg_positions, species, num_atoms,
      tolerance);
  *num_operations = required_operations;
  if (required_operations <= 0) {
    status = SAWF_SPGLIB_COUNT_FAILURE;
    goto cleanup;
  }
  if (output_capacity == 0) {
    status = SAWF_SPGLIB_SUCCESS;
    goto cleanup;
  }
  if (output_capacity < required_operations) {
    status = SAWF_SPGLIB_CAPACITY_TOO_SMALL;
    goto cleanup;
  }
  if (rotations == NULL || translations == NULL) goto cleanup;
  if ((size_t)output_capacity > SIZE_MAX / sizeof(*spg_rotations) ||
      (size_t)output_capacity > SIZE_MAX / sizeof(*spg_translations)) {
    status = SAWF_SPGLIB_ALLOCATION_FAILURE;
    goto cleanup;
  }
  spg_rotations = malloc((size_t)output_capacity * sizeof(*spg_rotations));
  spg_translations =
      malloc((size_t)output_capacity * sizeof(*spg_translations));
  if (spg_rotations == NULL || spg_translations == NULL) {
    status = SAWF_SPGLIB_ALLOCATION_FAILURE;
    goto cleanup;
  }

  returned_operations = spg_get_symmetry(
      spg_rotations, spg_translations, output_capacity, spg_lattice,
      (const double(*)[3])spg_positions, species, num_atoms, tolerance);
  *num_operations = returned_operations;
  if (returned_operations != required_operations) {
    status = SAWF_SPGLIB_API_FAILURE;
    goto cleanup;
  }

  for (operation = 0; operation < returned_operations; ++operation) {
    for (row = 0; row < 3; ++row) {
      translations[3 * operation + row] = spg_translations[operation][row];
      for (column = 0; column < 3; ++column) {
        rotations[9 * operation + 3 * column + row] =
            spg_rotations[operation][row][column];
      }
    }
  }
  status = SAWF_SPGLIB_SUCCESS;

cleanup:
  free(spg_translations);
  free(spg_rotations);
  free(spg_positions);
  return status;
}
